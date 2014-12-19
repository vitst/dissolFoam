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
  #include "initContinuityErrs.H"

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
  
  Info << "dissolFoamDict, fixInletConcentration:  " << fixInletConcentration <<nl;
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
  // Get patch ID for boundaries we want to move ("walls" "inlet")
  label wallID  = mesh.boundaryMesh().findPatchID("walls");
  label inletID = mesh.boundaryMesh().findPatchID("inlet");
  
  Info<< "Setup mesh relaxation class" << endl;
  // reference to the mesh relaxation object
  DissolMeshRlx* mesh_rlx = new DissolMeshRlx(mesh);
  
  scalarListList wallWeights = mesh_rlx->calc_weights( mesh.boundaryMesh()[wallID] );
  scalarListList inlWeights = mesh_rlx->calc_weights( mesh.boundaryMesh()[inletID] );
  
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
      int counter = 0;
      while ( simple.loop() ){
        counter++;
        Info<< "Time = " << runTime.timeName() 
            << "; Iteration: "<< counter 
            << nl << nl;

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
      counter = 0;
      while ( true ){
        counter++;

        double residual = solve
        (
          fvm::div(phi, C) - fvm::laplacian(D, C) == fvOptions(C)
        ).initialResidual();
         
        // @TODO make general function for convergence check
        if( residual < convCritCD ){
          Info << "Convection-diffusion: ExecutionTime = " << runTime.elapsedCpuTime() << " s"
               << "  ClockTime = " << runTime.elapsedClockTime() << " s"<< nl 
               << " Converged in "<< counter<<" steps" << nl << endl;
          if(counter >= maxNumIterCD){
            Info<<" dissolFoam Runtime WARNING:"
                    " Steady-state convection-diffusion solver did not converge.\n"
                << "It reached maximum number of iterations"
                    << "  counter: "<< counter
                    << "  maxNumIterCD: " << maxNumIterCD
                    << nl << endl;
          }
          break;
        }
        else{
          // @TODO probably here we should switch relTol to 0
          Info<< nl <<" Step "<< counter
                    <<"  residual: "<< residual
                    <<" > "<< convCritCD <<nl<<endl;
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

    pointDispInlet += mesh_rlx->inletRlx(mesh.boundaryMesh()[inletID], pdw);
      
//  Mesh update 4: Update boundary and relax interior mesh
    Info<< nl << "Update boundary and relax interior mesh" <<nl;
    pointVelocity.boundaryField()[wallID] == pointDispWall;
    pointVelocity.boundaryField()[inletID] == pointDispInlet;
    mesh.update();
    
    
/*    

    vectorField& pdd = pointDispInlet;
    mesh_rlx->doInletDisplacement(pdd);
    
    //std::exit(0);
    pointField newPoints = mesh.points();
    const labelList& inletToAll = mesh.boundaryMesh()[inletID].meshPoints();
    forAll(pointDispInlet, i){
      newPoints[ inletToAll[i] ] += pointDispInlet[i] * runTime.deltaTValue();
    }
    mesh.movePoints( newPoints );
    
    runTime.write();
    runTime++;
    */
    // ****************************************************************
    vectorField zeroInlet( pointDispInlet.size(), vector::zero );
    vectorField zeroWall( pointDispWall.size(), vector::zero );
    pointVelocity.boundaryField()[inletID] == zeroInlet;
    pointVelocity.boundaryField()[wallID] == zeroWall;
    // ****************************************************************
    
      //runTime.write();
      //runTime++;
//  Mesh update 5: boundary mesh relaxation
    Info<<nl<<"Boundary mesh relaxation"<<nl<<nl;
    Info<<"Wall"<<nl;
    for(int i=0; i<10; i++){
      vectorField wallRelax = mesh_rlx->wallRelaxation( mesh.boundaryMesh()[wallID], wallWeights );
      //pointVelocity.boundaryField()[wallID] == 10.0 * wallRelax;
      
      //vectorField inlRelax = mesh_rlx->wallRelaxation( mesh.boundaryMesh()[inletID], inlWeights );
      //pointVelocity.boundaryField()[inletID] == 10.0 * inlRelax;
      
      pointField newPoints = mesh.points();
      const labelList& wallToAll = mesh.boundaryMesh()[wallID].meshPoints();
      forAll(wallRelax, j){
        newPoints[ wallToAll[j] ] += wallRelax[j];
      }
      mesh.movePoints( newPoints );
      
      
      //mesh.update();
      
      //runTime.write();
      //runTime++;
    }
    mesh.update();
    
    /*
    scalarListList inlWeights11 = mesh_rlx->calc_weights( mesh.boundaryMesh()[inletID] );
    
    Pout << nl << nl;
    forAll(inlWeights, i){
      point curP = mesh.boundaryMesh()[inletID].localPoints()[i];
      if( curP.x()>23.5 && curP.x()<24.5 && curP.y()>.4 ){
        Pout << curP << "  " << inlWeights[i] << "  " << inlWeights11[i]<< nl;
      }
    }
    std::exit(0);
    */
    /*
    Info<<"Inlet"<<nl;
    for(int i=0; i<100; i++){
      vectorField boundaryRelax = mesh_rlx->wallRelaxation();
      pointVelocity.boundaryField()[wallID] == boundaryRelax;
      mesh.update();
    }
    */
    
    // *******************************
/*    
    // Internal points relaxation
    pointVelocity.boundaryField()[inletID] == zeroInlet;
    pointVelocity.boundaryField()[wallID] == zeroWall;
    
    for(int ij=0;ij<5; ij++){
      mesh.update();
    }
 */
    Info<< "Mesh update: ExecutionTime = " << runTime.elapsedCpuTime()
          << " s" << "  ClockTime = " << runTime.elapsedClockTime()
          << " s" << nl << nl;
  
    run0timestep = true;
    Info<< "Process: Time = " << runTime.timeName() << nl << nl;
  }

  Info << "End" << nl;
  return 0;
}

// **************************** End of the solver ******************************** //
