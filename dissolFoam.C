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

  bool gradCwrite(dissolProperties.lookupOrDefault<bool>("gradCwrite", false));
  
  #include "createFields.H"
  
  // l_T=D/(k*h_0)
  scalar l_T;
  if( !transportProperties.readIfPresent<scalar>("l_T", l_T) ){
    SeriousErrorIn("main")
            <<"There is no l_T parameter in transportProperties dictionary"
            <<exit(FatalError);
  }
  
  bool inertia;
  if( !dissolProperties.readIfPresent<bool>("inertia", inertia) ){
    SeriousErrorIn("main")
            <<"There is no constFlux parameter in transportProperties dictionary"
            <<exit(FatalError);
  }

  bool constFlux;
  if( !dissolProperties.readIfPresent<bool>("constFlux", constFlux) ){
    SeriousErrorIn("main")
            <<"There is no constFlux parameter in transportProperties dictionary"
            <<exit(FatalError);
  }

  Info << "*****************************************************************"<<nl;
  Info << "transportProperties, l_T:  " << l_T <<nl;
  Info << "dissolFoamDict, inertia:   " << inertia <<nl;
  Info << "dissolFoamDict, constFlux: " << constFlux <<nl;
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
  scalar areaCoef = fieldO->getInletArea(args, constFlux);
  
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  bool runTimeIs0 = (runTime.value()==0) ? true : false;
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  // MAIN TIME LOOP //
  while (runTime.run())
  {
    Info<< "Begin cycle: Time = " << runTime.timeName() 
            << "    dt = " << runTime.deltaTValue()
            << nl << endl;

    if( !runTimeIs0 )
    {
// *********************************************************
// *    Mesh motion & relaxation
// *********************************************************
      // calculate mesh motion
      vectorField pointDispWall = fieldO->getWallPointMotion(mesh, C, l_T, wallID);
      // move and relax the mesh (@TODO separate surface motion and relaxation)
      mesh_rlx->meshUpdate(pointDispWall, runTime);
    }  
// *********************************************************
// *    Stokes flow
// *********************************************************
    steadyStateControl simple(mesh);
    while ( simple.loop() )
    {
      if(inertia){
        #include "UEqn.H"
        #include "pEqn.H"
      }
      else{
        #include "UEqnStokes.H"
        #include "pEqn.H"
      }
    }
      
// *********************************************************
// *    Keeping flow rate constant
// *********************************************************
    if(constFlux){
      areaCoef = fieldO->getInletArea(args, constFlux);
    }
    scalar nU = fieldO->getConstFlowRateFactor(mesh, U, areaCoef);
    U   == U   / nU;
    phi == phi / nU;
      
// *********************************************************
// *    Danckwerts boundary condition loop, Convenction-Diffusion
// *********************************************************
    Info << "Steady-state convection-diffusion"<< endl;
    #include "convectionDiffusion.H"

    if(gradCwrite) maggradC == mag(-fvc::grad(C));
      
// *********************************************************
// *    Write Output data
// *********************************************************
    Info<<"Write fields"<<nl<<nl;
    (runTimeIs0) ? runTime.writeNow() : runTime.write();
    runTime++;
    runTimeIs0 = false;
    
    Info<<"Mesh update: ExecutionTime = "<<runTime.elapsedCpuTime()<<" s"
        <<"  ClockTime = " << runTime.elapsedClockTime()<<" s"<<nl<<nl;
    
    Info<<"Process: Time = " << runTime.timeName() << nl << nl;
  }

  Info<<"End"<<nl;
  return 0;
}

// **************************** End of the solver ******************************** //
