/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

// common OFe solver and OF simpleFoam
#include "fvCFD.H"

// OF simpleFoam
#include "singlePhaseTransportModel.H"
#include "fvIOoptionList.H"

// OFe
#include "pointPatchField.H"
#include "dynamicFvMesh.H"
#include "syncTools.H"

// dissol prj
#include "dissolControl.H"
#include "primitivePatchInterpolationSync.H"
#include "DissolMeshRlx.H"

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
 #######################################################################################
*/
int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createDynamicFvMesh.H"
  #include "createFields.H"
  //#include "initContinuityErrs.H"

  // OF simpleFoam  
  #include "createFvOptions.H" 
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  // reading dissolFoam dictionary
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

  bool fixInletConcentration
  (
    readBool( dissolProperties.lookup("fixInletConcentration") )
  );
  scalar rlxTol( dissolProperties.lookupOrDefault<scalar>("relaxationTolerance", 0.1) );
  
  scalar inigradingZ( dissolProperties.lookupOrDefault<scalar>("inigradingZ", 1.0) );
  scalar timeCoefZ( dissolProperties.lookupOrDefault<scalar>("timeCoefZ", 1.0) );
  int Nz( dissolProperties.lookupOrDefault<int>("numberOfCellsZ", 10) );
  
  scalar inigradingY( dissolProperties.lookupOrDefault<scalar>("inigradingY", 1.0) );
  scalar timeCoefY( dissolProperties.lookupOrDefault<scalar>("timeCoefY", 1.0) );
  int Ny( dissolProperties.lookupOrDefault<int>("numberOfCellsY", 4) );
  
  Info << "*****************************************************************"<<nl;
  Info << "dissolFoamDict, fixInletConcentration:  " << fixInletConcentration <<nl;
  Info << "dissolFoamDict, rlxTol:  " << rlxTol <<nl;
  Info << "dissolFoamDict, inigradingZ:  " << inigradingZ <<nl;
  Info << "dissolFoamDict, timeCoefZ:  " << timeCoefZ <<nl;
  Info << "dissolFoamDict, numberOfCellsZ:  " << Nz <<nl;
  Info << "dissolFoamDict, inigradingY:  " << inigradingY <<nl;
  Info << "dissolFoamDict, timeCoefY:  " << timeCoefY <<nl;
  Info << "dissolFoamDict, numberOfCellsY:  " << Ny <<nl;
  Info << "*****************************************************************"<<nl;
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
  // Get patch ID for boundaries we want to move ("walls" "inlet")
  label wallID  = mesh.boundaryMesh().findPatchID("walls");
  label inletID = mesh.boundaryMesh().findPatchID("inlet");
  label outletID = mesh.boundaryMesh().findPatchID("outlet");
  
  Info<< "Setup mesh relaxation class" << endl;
  // reference to the mesh relaxation object
  DissolMeshRlx* mesh_rlx = new DissolMeshRlx(mesh);
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
  /*
   * run0timestep is used for skipping the Stokes and the convection-diffusion solvers
   * if timestep is not 0. At the end of each cycle it is set to true.
   */
  bool run0timestep = false;
  if( runTime.value() == 0 ){  run0timestep = true; }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  // main time loop
  while (runTime.run())
  {
    Info<< "Begin cycle: Time = " << runTime.timeName() << nl << endl;
    
    // Skip the Stokes and the convection-diffusion solvers if timestep is not 0
    // because it was done within the last timestep
    if( run0timestep )
    {
      /*##########################################
       *   Stokes flow
       *##########################################*/
      dissolControl simple(mesh);
      while ( simple.loop() ){
        // Pressure-velocity SIMPLE corrector
        #include "UEqn.H"
        #include "pEqn.H"
        //turbulence->correct();
      }

      /*############################################
      *   Steady-state convection-diffusion solver
      * ##########################################*/
      Info << "Steady-state convection-diffusion"<< endl;
      // define dictionary for convection diffusion solver
      dictionary conv_diff = mesh.solutionDict().subDict("CONVECTION_DIFFUSION");
      // initialize values for convergence checks
      double convCritCD = 0;   // convergence criteria
      int maxNumIterCD = 0; // maximum number of iterations
      // read values for convergence checks from the dictionary
      conv_diff.readIfPresent("convergence", convCritCD);
      conv_diff.readIfPresent("maxIter", maxNumIterCD);

      // Steady-state convection-diffusion solver main loop
      int counter = 0;
      while ( true ){
        counter++;

        double residual = solve
        (
          fvm::div(phi, C) - fvm::laplacian(D, C) == fvOptions(C)
        ).initialResidual();

        if( residual < convCritCD ){
          Info << "Convection-diffusion: ExecutionTime = " << runTime.elapsedCpuTime() << " s"
               << "  ClockTime = " << runTime.elapsedClockTime() << " s"<< nl 
               << "Converged in "<< counter <<" steps" << nl << endl;
          if(counter >= maxNumIterCD){
            Info<<" dissolFoam Runtime WARNING:"
                    << " Steady-state convection-diffusion solver did not converge."
                    << nl
                    << "It reached maximum number of iterations"
                    << "  counter: "<< counter
                    << "  maxNumIterCD: " << maxNumIterCD
                    << nl;
          }
          break;
        }
        else{
          Info<< nl <<" Step "<< counter
                    <<"  residual: "<< residual
                    <<" > "<< convCritCD <<nl;
        }
      }
      
      /*############################################################################
       *   Write Output data
       * ###########################################################################*/
      if( runTime.value() == 0 ){ 
        runTime.writeNow(); 
      }
      else{
        runTime.write();
      }
      Info<< "Write data, after conv-diff" << nl << nl;
    }
    else{
      Info<< "dissolFoam: Skip the Stokes and the convection-diffusion solvers."<<nl<<nl;
    }
    runTime++;
    
    /*##############################################################################
    *   Mesh relaxation
    * ##############################################################################*/
    
//  Mesh update 1: 
    pointVectorField& pointVelocity = const_cast<pointVectorField&>(
      mesh.objectRegistry::lookupObject<pointVectorField>( "pointMotionU" )
    );

//  Mesh update 2.1: Calculate new displacements (interpolate C field to points)  

    // interpolate the concentration from cells to wall faces
    primitivePatchInterpolationSync patchInterpolator( mesh.boundaryMesh()[wallID], mesh );
    
    scalarField pointCface = -C.boundaryField()[wallID].snGrad();
    vectorField pointNface = mesh.boundaryMesh()[wallID].faceNormals();
    vectorField motionVec = pointCface * pointNface;

    vectorField pointDispWall = patchInterpolator.faceToPointInterpolate(motionVec);
    
    vectorField& pdw = pointDispWall;
    
    if(fixInletConcentration){
      Info << "Fix concentration on the edge between walls and inlet"<< nl;
      mesh_rlx->fixEdgeConcentration(pdw);
    }
    
//  Mesh update 2.2: Inlet displacement
    vectorField pointDispInlet = mesh_rlx->calculateInletDisplacement(pdw);
    //pointDispInlet += mesh_rlx->inletRlx(mesh.boundaryMesh()[inletID], pdw);
            
    vectorField pointDispOutlet = mesh_rlx->calculateOutletDisplacement(pdw);
    //pointDispOutlet += mesh_rlx->outletRlx(mesh.boundaryMesh()[outletID], pdw);
      
//  Mesh update 4: Update boundary and relax interior mesh
//    Info<< nl << "Update boundary and relax interior mesh" <<nl;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    //pointVelocity.boundaryField()[wallID] == pointDispWall;
    //pointVelocity.boundaryField()[inletID] == pointDispInlet;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    //pointVelocity.boundaryField()[outletID] == pointDispOutlet;
    
    //mesh.update();
    
    // ****************************************************************
    //vectorField zeroOutlet( pointDispOutlet.size(), vector::zero );
    //vectorField zeroInlet( pointDispInlet.size(), vector::zero );
    //vectorField zeroWall( pointDispWall.size(), vector::zero );
    //pointVelocity.boundaryField()[inletID] == zeroInlet;
    //pointVelocity.boundaryField()[wallID] == zeroWall;
    //pointVelocity.boundaryField()[outletID] == zeroOutlet;
    // ****************************************************************
    
//  Mesh update 5: boundary mesh relaxation
    
    Info<<nl<<"Calculating new grading...."<< nl;
            
    scalar Gz = inigradingZ / (timeCoefZ * runTime.value() + 1.0);
    scalar lambdaZ = 1/static_cast<double>(Nz) * std::log( Gz );
    
    scalarListList wallWeights = mesh_rlx->calc_weights( mesh.boundaryMesh()[wallID], lambdaZ, 3);
    const scalarListList& wW = wallWeights;
    
    scalar Gy = inigradingY / (timeCoefY * runTime.value() + 1.0);
    scalar lambdaY = 1/static_cast<double>(Ny) * std::log( Gy );
    
    scalarListList inletWeights = mesh_rlx->calc_weights( mesh.boundaryMesh()[inletID], lambdaY, 2);
    const scalarListList& iW = inletWeights;
    
    scalarListList outletWeights = mesh_rlx->calc_weights( mesh.boundaryMesh()[outletID], lambdaY, 2);
    const scalarListList& oW = outletWeights;
    
    Info<<nl<<"Boundary mesh relaxation. Z grading is "<< Gz
            <<"   Y grading is "<< Gy <<nl<<nl;
    
    pointField savedPointsAll = mesh.points();
    
    const vectorField& wR = pointDispWall;
    pointField mpW = mesh_rlx->doWallDisplacement(wR * runTime.deltaTValue() );
    mesh.movePoints( mpW);
    
    const vectorField& iR = pointDispInlet;
    pointField mpI = mesh_rlx->doInletDisplacement(iR * runTime.deltaTValue());
    mesh.movePoints( mpI );
    
    const vectorField& oR = pointDispOutlet;
    pointField mpO = mesh_rlx->doOutletDisplacement(oR * runTime.deltaTValue());
    mesh.movePoints( mpO );
    //runTime.write();
    //runTime++;
    
    Info<<"Relaxing wall... time: "<< runTime.cpuTimeIncrement() <<nl;
    vectorField wallRelax = mesh_rlx->wallRelaxation( mesh.boundaryMesh()[wallID], wW, rlxTol);
    Info<<"Wall relaxation time: " << runTime.cpuTimeIncrement() << " s" << nl;

    Info<<"Relaxing inlet..."<<nl;
    vectorField inlRelax = mesh_rlx->wallRelaxation( mesh.boundaryMesh()[inletID], iW, rlxTol);
    Info<<"Inlet relaxation time: " << runTime.cpuTimeIncrement() << " s" << nl;

    Info<<"Relaxing outlet..."<<nl;
    vectorField outRelax = mesh_rlx->wallRelaxation( mesh.boundaryMesh()[outletID], oW, rlxTol);
    Info<<"Outlet relaxation time: " << runTime.cpuTimeIncrement() << " s" << nl;
    
    mesh.movePoints( savedPointsAll );

    wallRelax /= runTime.deltaTValue();
    pointVelocity.boundaryField()[wallID] == wallRelax + pointDispWall;

    inlRelax /= runTime.deltaTValue();
    pointVelocity.boundaryField()[inletID] == inlRelax + pointDispInlet;

    outRelax /= runTime.deltaTValue();
    pointVelocity.boundaryField()[outletID] == outRelax + pointDispOutlet;
    
    mesh.update();

    //runTime.write();
    //runTime++;
    
    Info<< "Mesh update: ExecutionTime = " << runTime.elapsedCpuTime()
          << " s" << "  ClockTime = " << runTime.elapsedClockTime()
          << " s" << nl << nl;
    
    //std::exit(0);
  
    run0timestep = true;
    Info<< "Process: Time = " << runTime.timeName() << nl << nl;
  }

  Info << "End" << nl;
  return 0;
}

// **************************** End of the solver ******************************** //
