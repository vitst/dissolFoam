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
  #include "createFields.H"

  Info << "*****************************************************************"<<nl;
  Info << "transportProperties, l_T:  " << l_T <<nl;
  Info << "dissolFoamDict, inertia:   " << inertia <<nl;
  Info << "dissolFoamDict, constFlux: " << constFlux <<nl;
  Info << "*****************************************************************"<<endl;
  
  // Get patch ID for boundaries we want to move ("walls" "inlet")
  label wallID  = mesh.boundaryMesh().findPatchID("walls");
  //label inletID = mesh.boundaryMesh().findPatchID("inlet");
  //label outletID = mesh.boundaryMesh().findPatchID("outlet");
  
  Info<< "Setup mesh relaxation class" << endl;
  meshRelax* mesh_rlxPtr = new meshRelax(mesh, args);       // pointer to the mesh relaxation object
  Info<< "Setup field operation class" << endl;
  fieldOperations* fieldOpPtr = new fieldOperations(); // pointer to the mesh relaxation object
  
  // calculating initial area of the inlet in order to scale U later
  scalar areaCoef = fieldOpPtr->getInletArea(args, constFlux);
  
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  bool runTimeIs0 = (runTime.value()==0) ? true : false;
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  // MAIN TIME LOOP //
  while (runTime.run())
  {
    if( !runTimeIs0 )
    {
      runTime++;                         // Don't update at t=0
      Info << "Begin cycle: Time = " << runTime.timeName() 
           << "    dt = " << runTime.deltaTValue()
           << nl << endl;


// *********************************************************
// *    Mesh motion & relaxation
// *********************************************************
      // calculate mesh motion
      vectorField pointDispWall = fieldOpPtr->getWallPointMotion(mesh, C, l_T, wallID);
      // move and relax the mesh (@TODO separate surface motion and relaxation)
      
      /*
      for(int i=0;i<10;i++){
        Info<<i<<"  "<<pointDispWall[i]<<nl;
      }
      */
      
      mesh_rlxPtr->meshUpdate(pointDispWall, runTime);
      
      Info << "Mesh update: ExecutionTime = " << runTime.elapsedCpuTime()
           << " s" << "  ClockTime = " << runTime.elapsedClockTime()
           << " s"<< nl<< endl;
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
      
    Info << "Flow solver: "
         << "ExecutionTime = " << runTime.elapsedCpuTime() << " s "
         <<"ClockTime = "<< runTime.elapsedClockTime() << " s"
         << nl << nl << endl;

// *********************************************************
// *    Keeping flow rate constant
// *********************************************************
    if(constFlux){
      areaCoef = fieldOpPtr->getInletArea(args, constFlux);
    }
    scalar nU = fieldOpPtr->getConstFlowRateFactor(mesh, U, areaCoef);
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
    Info << "Write fields: Time = " << runTime.timeName() << nl <<endl;
    (runTimeIs0) ? runTime.writeNow() : runTime.write();
    runTimeIs0 = false;
  }

  Info << "End" << endl;
  return 0;
}

// **************************** End of the solver ******************************** //
