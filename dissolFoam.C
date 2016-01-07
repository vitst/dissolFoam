/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-201X OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    dissolFoam

Description
    Solves for flow (Stokes) and transport (steady-state) and moves
    the mesh according to the reactant flux.
\*---------------------------------------------------------------------------*/

#include <gsl/gsl_sf_hyperg.h>

//// OF includes
// common OFe solver and OF simpleFoam
#include "fvCFD.H"

// OF simpleFoam
#include "singlePhaseTransportModel.H"
#include "fvIOoptionList.H"
//#include "RASModel.H"  // currently no turbulence model

// OFe
#include "pointPatchField.H"
#include "dynamicFvMesh.H"
#include "syncTools.H"

// dissol prj
#include "steadyStateControl.H"
#include "meshRelax.H"
#include "fieldOperations.H"

//// auxiliary includes
// mesh search
#include "interpolation.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "meshSearch.H"

// interpolation and tables
#include "interpolationTable.H"

#include "memInfo.H"

/*
 #######################################################################################
 *    Main program body
 * 
 * Schematic fracture geometry:
 * 
 *           <------ outlet ------>
 *        ||                         |          |
 *        ||                         |          |
 *        ||<----- walls             |          |
 *        ||                         |          |
 * z     /  \                   z    |          |
 * |_ y      <------ inlet -->  |_ x
 * 
 * controlDictionary should have a subdictionary: CONVECTION_DIFFUSION, e.g.:
 * 
 *      CONVECTION_DIFFUSION
 *      {
 *        convergence                 1e-9;
 *        maxIter                     2000;
 *      }
 *
 * @TODO
 * - create all the lists for time 0
 * - check the parameter preservePatches for cyclic BC
 * 
 * 
 #######################################################################################
*/
int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createDynamicFvMesh.H"

  // OF simpleFoam  
  #include "createFvOptions.H" 
  
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  // reading dissolFoam and transportProperties dictionary
  const word dissolDictName("dissolFoamDict");

  IOdictionary dissolProperties
  (
    IOobject
    (
      dissolDictName,
      runTime.system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );
  
  IOdictionary transportProperties
  (
      IOobject
      (
          "transportProperties",
          runTime.constant(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      )
  );

  // if true it switches on the convection term in Navier-Stokes eqn
  // moved to transport properties
  bool NStokesInertia(transportProperties.lookupOrDefault<bool>("NStokesInertia", false));
  bool gradCwrite(dissolProperties.lookupOrDefault<bool>("gradCwrite", false));
  
  #include "createFields.H"
  
  // l_T=D/(k*h_0)
  // moved to transport properties
  scalar l_T( transportProperties.lookupOrDefault<scalar>("l_T", 1.0) );
  
  // Reynolds number
  // moved to transport properties
  scalar Re( transportProperties.lookupOrDefault<scalar>("Re", 1.0/nu.value()) );
  
  // moved to transport properties
  bool constFlux( transportProperties.lookupOrDefault<bool>("constFlux", false) );
  
  Info << "*****************************************************************"<<nl;
  Info << "dissolFoamDict, NStokesInertia:  " << NStokesInertia <<nl;
  Info << "dissolFoamDict, Re:  " << Re <<nl;
  Info << "dissolFoamDict, constFlux:  " << constFlux <<nl;
  Info << "*****************************************************************"<<nl;

  // Get patch ID for boundaries we want to move ("walls" "inlet")
  label wallID  = mesh.boundaryMesh().findPatchID("walls");
  //label inletID = mesh.boundaryMesh().findPatchID("inlet");
  //label outletID = mesh.boundaryMesh().findPatchID("outlet");
  
  Info<< "Setup mesh relaxation class" << endl;
  meshRelax* mesh_rlx = new meshRelax(mesh); // pointer to the mesh relaxation object
  Info<< "Setup field operation class" << endl;
  fieldOperations* fieldO = new fieldOperations(); // pointer to the mesh relaxation object
  
  // calculating initial area of the inlet in order to scale U later
  scalar areaCoef0 = fieldO->getInletArea(args);
  
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  /*
   * A variable run0timestep is used to skip the Navier-Stokes and 
   * the convection-diffusion solvers if current time is not 0.
   * At the end of each cycle it is set to true.
   */
  bool run0timestep = (runTime.value()==0) ? true : false;
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  // MAIN TIME LOOP //
  while (runTime.run())
  {
    Info<< "Begin cycle: Time = " << runTime.timeName() << nl << endl;
    
    // At every restart it skips the Stokes and the convection-diffusion solvers
    // if timestep is not 0 because it was done within the last timestep.
    if( run0timestep )
    {
// *********************************************************
// *    Stokes flow
// *********************************************************
      steadyStateControl simple(mesh);
      while ( simple.loop() )
      {
        if(NStokesInertia)
        {
          #include "UEqn.H"
          #include "pEqn.H"
        }
        else
        {
          #include "UEqnStokes.H"
          #include "pEqn.H"
        }
      }
      
// *********************************************************
// *    Keeping flow rate constant
// *********************************************************
      if(constFlux){
        scalar nU = fieldO->getConstFlowRateFactor(mesh, U, areaCoef0);
        U==U/nU;
      }
      
// *********************************************************
// *    Danckwerts boundary condition loop, Convenction-Diffusion
// *********************************************************
      Info << "Steady-state convection-diffusion"<< endl;
      //int inlet_count = 0;
      //while ( true )
      //while(inlet_count<3)
      //{
        #include "convectionDiffusion.H"
      //  scalar tttol = fieldO->calcDanckwerts(mesh, U, C, inletID, D.value());
      //  inlet_count+=1;
      //  Info<<"Inlet C iter "<<inlet_count<<"  tolerance: "<<tttol<< nl;
      //  if(tttol<convCritCD) break;
      //}
      
      if(gradCwrite) maggradC == mag(-fvc::grad(C));
      
// *********************************************************
// *    Write Output data
// *********************************************************
      (runTime.value()==0) ? runTime.writeNow() : runTime.write();
      Info<<"Write data, after conv-diff"<<nl<<nl;
    }
    else
    {
      Info<<"dissolFoam: Skip the Stokes and the convection-diffusion solvers."<<nl<<nl;
    }
    runTime++;
    
// *********************************************************
// *    Mesh motion & relaxation
// *********************************************************
    // calculate mesh motion
    vectorField pointDispWall = fieldO->getWallPointMotion(mesh, C, l_T, wallID);
    // move and relax the mesh (@TODO separate surface motion and relaxation)
    mesh_rlx->meshUpdate(pointDispWall, runTime);
    
    Info<<"Mesh update: ExecutionTime = "<<runTime.elapsedCpuTime()<<" s"
        <<"  ClockTime = " << runTime.elapsedClockTime()<<" s"<<nl<<nl;
    
    run0timestep = true;
    Info<<"Process: Time = " << runTime.timeName() << nl << nl;
  }

  Info<<"End"<<nl;
  return 0;
}

// **************************** End of the solver ******************************** //
