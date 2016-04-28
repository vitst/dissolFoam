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

// OF includes

// common and simpleFoam
#include "fvCFD.H"

// simpleFoam
//#include "fvIOoptionList.H"

// OF
#include "pointPatchField.H"
#include "dynamicFvMesh.H"
#include "syncTools.H"

// dissol project
#include "steadyStateControl.H"
#include "meshRelax.H"
#include "fieldOperations.H"

// auxiliary: includes mesh search
#include "interpolation.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "meshSearch.H"

// interpolation and tables
#include "interpolationTable.H"
#include "memInfo.H"

/*######################################################################
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
 * fvSolution has a subdictionary: convDiff, e.g.:
 * 
 *      convDiff
 *      {
 *        convergence                 1e-9;
 *        maxIter                     2000;
 *      }
 *
 * @TODO
 * - create all the lists for time 0
 * - check the parameter preservePatches for cyclic BC
 * - make mesh relax independent on the inletID and outletID
 * 
 #####################################################################*/

int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createDynamicFvMesh.H"

  // OF simpleFoam  
  //#include "createFvOptions.H" 
  
  #include "readDicts.H"
  #include "createFields.H"

  // Get patch ID for moving boundaries ("walls" "inlet")
  const label wallID  = mesh.boundaryMesh().findPatchID("walls");
  //const label inletID = mesh.boundaryMesh().findPatchID("inlet");
  const label outletID = mesh.boundaryMesh().findPatchID("outlet");
  
  meshRelax mesh_rlx(mesh, args);
  Info << "Setup mesh relaxation object. meshRelax version is "
          << mesh_rlx.get_version() << endl;

  fieldOperations fieldOp(args, outletID);
  Info << "Setup field operation object" << endl;
  
  Info<<"Patch \""<<mesh.boundaryMesh().names()[outletID]<<"\" is used for scaling U"<<nl;
  
  //const scalar scalingAreaT0 = fieldOp.getScalingAreaT0();
  Info << "Initial patch area for the flow rate scaling: " 
          << fieldOp.getScalingAreaT0() << nl;
  Info << "U field is scaled so that <U_inlet> = 1 at t=0" << nl;
  Info << "phi field is unscaled" << nl
       << "Used to recover scale factor after restart" << endl;
  
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
      vectorField pointDispWall = 
            fieldOp.getWallPointMotion(mesh, C, l_T, wallID);
      
      mesh_rlx.meshUpdate(pointDispWall, runTime);
      
      Info << "Mesh update: ExecutionTime = " << runTime.elapsedCpuTime()
           << " s" << "  ClockTime = " << runTime.elapsedClockTime()
           << " s"<< nl<< endl;
    }
    
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
    }

    Info << "Flow solver: "
    << "ExecutionTime = " << runTime.elapsedCpuTime() << " s "
    <<"ClockTime = "<< runTime.elapsedClockTime() << " s"
    << nl << nl << endl;

/*###############################################
 *   Scale U and phi (using outlet value)
 *   constFlux == true  -> constant flow rate
 *   constFlux == false -> constant pressure drop
 * 
 *   limitFlux == true  -> limit a flow rate in case of
 *                         constant pressure drop
 *   limitValue = 3.0   -> limit a flow rate to 3.0 * rate at time 0
 * -----------------------------------------------------------------------
 *   
 * 
 * 
 *###############################################*/
    scalar Q  = fieldOp.getScalingFlowRate(phi)  ;
    scalar Q0 = fieldOp.getScalingFlowRateT0(phi);
    scalar A0 = fieldOp.getScalingAreaT0();
    
    scalar nU = 1.0;
    if(rescale && !inertia)
      nU = Q0/A0;
      
    if(constFlux)
    {
      nU *= Q / Q0;
    }
    else
    {
      if(limitFlux && Q > Q0*limitValue)
        nU *= Q / (Q0*limitValue);
    }
      
    //scalar nU = Q / Q0;
    Info << "U and phi scale factor: " << nU << "   Q: "<< Q << nl << endl;

    U   == U   / nU;
    phi == phi / nU;

/*############################################
 *   Steady-state convection-diffusion solver
 *##########################################*/

  Info << "Steady-state concentration solver"<< endl;

  int iter = 0;
  while ( true ){
    iter++;
    
    double residual = solve
    (
      fvm::div(phi, C) - fvm::laplacian(D, C) // == fvOptions(C)
    ).initialResidual();
    
    if( residual < convCrit ){
      Info << "Convection-diffusion: "
           << "ExecutionTime = " << runTime.elapsedCpuTime() << " s "
           << "ClockTime = " << runTime.elapsedClockTime() << " s" <<nl 
           << "Converged in " << iter << " steps" << nl << endl;

      if(iter >= maxIter){
        Info << nl << "dissolFoam Runtime WARNING:"
             << "Convection-diffusion solver did not converge." << nl
             <<"Maximum number of iterations"
             <<"  iter: "<< iter << endl;
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
    
    // Rescale phi for restart and for next cycle
    phi == phi * nU;

    Info << "Write fields: Time = " << runTime.timeName() << nl <<endl;
    (runTimeIs0) ? runTime.writeNow() : runTime.write();
    runTimeIs0 = false;

    // Rescale U for next cycle
    U == U * nU;
  }

  Info << "End" << endl;
  return 0;
}

// ********************* End of the solver ************************** //
