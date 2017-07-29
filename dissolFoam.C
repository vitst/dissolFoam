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
 
    Works with OpenFOAM 2.x.x.
\*---------------------------------------------------------------------------*/

// common and simpleFoam
#include "fvCFD.H"

// OF
#include "pointPatchField.H"
#include "dynamicFvMesh.H"
#include "syncTools.H"

// dissolFOAM project
#include "steadyStateControl.H"
//#include "meshRelax.H"
#include "fieldOperations.H"
#include "dissolMotionPointPatchVectorField.H"
#include "pointPatchField.H"

/*####################################################################
 *    Main program body
 * 
 * Schematic fracture geometry:
 * 
 *           <------ outlet ------>
 *        ||                         |          |
 *        ||                         |          |
 *        ||<------- walls           |          |
 *        ||                         |          |
 * z     /  \                   z    |          |
 * |_ y      <------ inlet -->  |_ x
 * 
 * fvSolution has a subdictionary: convDiff, e.g.:
 * 
 *      convDiff
 *      {
 *        convergence                 1e-9;
 *        maxIter                     2000;
 *      }
 *
 #####################################################################*/

int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createDynamicFvMesh.H"

  #include "readDicts.H"
  #include "createFields.H"

  // Get patch ID for moving boundaries ("walls")
  //const label patchMotion  = mesh.boundaryMesh().findPatchID("walls");
  // Get patch ID for scaling flow rate
  const label patchScaling = mesh.boundaryMesh().findPatchID("outlet");
  
  //meshRelax mesh_rlx(mesh, args);
  //Info << "Setup mesh relaxation object. meshRelax version is "
  //        << mesh_rlx.get_version() << endl;

  fieldOperations fieldOp(args, patchScaling);
  Info << "Setup field operation object" << endl;
  
  Info<<"Patch \""<<mesh.boundaryMesh().names()[patchScaling]
          <<"\" is used for scaling U"<<nl;
  
  
  
  
  
  //Info<<nl<<"Patches  "<<mesh.boundaryMesh().names() <<nl;


  /*
  Info<<"Check access to dissolMotion boundary parameters"<<endl;
  const label patchMotion  = mesh.boundaryMesh().findPatchID("walls");
  pointVectorField& pointMotionU = const_cast<pointVectorField&>(
    mesh.objectRegistry::lookupObject<pointVectorField>( "pointMotionU" )
  );
   
  dissolMotionPointPatchVectorField& dM = 
      dynamic_cast<dissolMotionPointPatchVectorField&>
      (
        const_cast<pointPatchField<vector>&>( pointMotionU.boundaryField()[patchMotion]) 
      ); 
 
  Info<< "surfaceRlx: "<< dM.getSurfaceRlx()<<nl;
  dM.setSurfaceRlx(false);
  Info<< "surfaceRlx: "<< dM.getSurfaceRlx()<<nl;

  Info<< "edgeRlx: "<< dM.getEdgeRlx()<<nl;
  dM.setEdgeRlx(false);
  Info<< "edgeRlx: "<< dM.getEdgeRlx()<<nl;
  */

  
//  return 0;
  
  
// * * * * *   MAIN LOOP   * * * * * * * * * * * * * * * * * * * * * //

  runTime.functionObjects().execute();     // Execute cntlDict functions
  bool runTimeIs0 = (runTime.value()==0) ? true : false;

  while (runTime.run())
  {
    if( !runTimeIs0 )                      // No mesh update at t=0
    {
      runTime++;
      Info << "Begin cycle: Time = " << runTime.timeName() 
           << "    dt = " << runTime.deltaTValue()
           << nl << endl;

// *********************************************************
// *    Mesh motion & relaxation
// *********************************************************

      // calculate mesh motion 
      //vectorField pointDisplacement = 
      //      l_T * fieldOp.getWallPointMotion(mesh, C, patchMotion);
      
      //mesh_rlx.meshUpdate(pointDisplacement, runTime);
      
      mesh.update();
      //runTime.write();
      Info << "Mesh update: ExecutionTime = " << runTime.elapsedCpuTime()
           << " s" << "  ClockTime = " << runTime.elapsedClockTime()
           << " s"<< nl<< endl;
//      runTime++;
//      runTime.write();
    }


//    return 0;
    
/*###############################################
 *   Steady-state flow solver
 *###############################################*/

    steadyStateControl simple(mesh);
    while ( simple.loop() )
    {
      if(inertia)
      {
        #include "UEqn.H"        // Navier Stokes
        #include "pEqn.H"
      }
      else
      {
        #include "UEqnStokes.H"  // Stokes
        #include "pEqn.H"
      }
      
      if(limitFlux)
      {
        scalar Q = fieldOp.getScalingFlowRate(phi);
        scalar compareTo = (constFlux) ? SMALL : limitValue;
        scalar scaleFactor = ( Q > compareTo ) ? limitValue / Q : 1.0;
                
        p == scaleFactor * p;
        U == scaleFactor * U;
      }
    }
    
    Info << "Final flow rate: " << fieldOp.getScalingFlowRate(phi)<<nl<<nl;

    Info << "Flow solver: "
         << "ExecutionTime = " << runTime.elapsedCpuTime() << " s "
         << "ClockTime = "<< runTime.elapsedClockTime() << " s"
         << nl << nl << endl;

/*##########################################
 *   Steady-state convection-diffusion solver
 *##########################################*/

    Info << "Steady-state concentration solver"<< endl;

    int  iter = 0;
    while ( true )
    {
      iter++;

      double residual = solve
      (
        fvm::div(phi, C) - fvm::laplacian(D, C)
      ).initialResidual();
      
      if( residual < convCrit )
      {
        Info << "Convection-diffusion: "
             << "ExecutionTime = " << runTime.elapsedCpuTime() << " s "
             << "ClockTime = " << runTime.elapsedClockTime() << " s" <<nl 
             << "Converged in " << iter << " steps.  Residual="<< residual
             << nl << endl;

        if(iter >= maxIter)
        {
          Info << nl << "dissolFoam Runtime WARNING:"
               << "Convection-diffusion solver did not converge." << nl
               << "Maximum number of iterations"
               << "  iter: "<< iter << endl;
        }
        break;
      }
      else{
        Info << " Step " << iter
             << " residual: "<< residual << " > " << convCrit << endl;
      }
    }
  
// *********************************************************
// *    Write Output data
// *********************************************************

    if(gradCWrite) gradC == -fvc::snGrad(C);
    
    Info << "Write fields: Time = " << runTime.timeName() << nl <<endl;
    (runTimeIs0) ? runTime.writeNow() : runTime.write();
    runTimeIs0 = false;
  }

  Info << "End" << endl;
  return 0;
}

// ********************* End of the solver ************************** //
