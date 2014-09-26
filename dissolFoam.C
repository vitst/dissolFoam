/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    dissolFoam

Description
    Solves for flow (Stokes) and transport (steady-state) and moves
    the mesh according to the reactant flux.
  
 * 
 *  Surface relaxation added
 * 

\*---------------------------------------------------------------------------*/



// OF includes
#include "fvCFD.H"
#include "pointPatchField.H"
#include "primitivePatchInterpolation.H"
#include "dynamicFvMesh.H"
#include "primitivePatch.H"
#include "partialSlipModPointPatchFields.H"


#include "interpolation.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// General includes
#include <list>
#include <math.h>
#include <stdlib.h>
#include <map>
#include <iostream>
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Local includes
#include "DissolMeshRlx.H"


// #######################################################################################
// ##     Auxiliary functions
// #######################################################################################


/*
 * Description
 *   It calculates in how many substeps one time-step should be divided.
 * TODO
 *   parallel computation
 */
int setSubStep(const edgeList& edges, const pointField& points, const scalarField& concField){
  // define the maximum length of the displacement vector
  scalar max_mesh_step = 0;
  forAll(concField, ii){
    max_mesh_step = max(max_mesh_step, concField[ii]);
  }
  
  // look for a smallest edge
  scalar min_edge = edges[0].mag(points);
  forAll(edges, i){
    min_edge = min(min_edge, edges[i].mag(points) );
  }
  
  return static_cast<int>(max_mesh_step / min_edge + 1);
}


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
 *      CONVECTION_DIFFUSION{
 *        convergence                 1e-9;
 *        maxIter                     2000;
 *      }
 *
 * 
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  // Get patch ID for boundaries we want to move ("walls" "inlet")
  label wallID  = mesh.boundaryMesh().findPatchID("walls");
  label inletID = mesh.boundaryMesh().findPatchID("inlet");

  // mesh rlx
  Foam::Time runTime1
  (
      Foam::Time::controlDictName,
      args.rootPath(),
      args.caseName(),
      "system",
      "constant",
      !args.optionFound("noFunctionObjects")
  );
  
  Foam::instantList timeDirs = Foam::timeSelector::select0(runTime1, args);

  Foam::fvMesh mesh1
  (
    Foam::IOobject
    (
      Foam::fvMesh::defaultRegion,
      runTime1.timeName(),
      runTime1,
      Foam::IOobject::MUST_READ
    )
  );

  // reference to the mesh relaxation object
  DissolMeshRlx* mesh_rlx = new DissolMeshRlx(mesh, mesh1);
  
  //mesh.
  //std::exit(0);

  /*
   * run0timestep is used for skipping the Stokes and the convection-diffusion solvers
   * if timestep is not 0. At the end of each cycle it is set to true.
   */
  bool run0timestep = false;
  if( runTime.value() == 0 ){  run0timestep = true; }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  /*
  const fvBoundaryMesh& patches = mesh.boundary();
  Info<< nl << nl << " !!Patchy: " << patches.size() << nl;
  Info<< " !!Patchy0: " << patches[0].patch() << nl;
  Info<< " !!Patchy1: " << patches[1].patch() << nl;
  Info<< " !!Patchy2: " << patches[2].patch() << nl;
  Info<< " !!Patchy3: " << patches[3].patch() << nl;
  Info<< " !!Patchy4: " << patches[4].patch() << nl;
  std::exit(0);
   */

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
      int counter = 0; // counter for steady state solver iterations
      while ( true ){
        counter++;

        #include "readSIMPLEControls.H"

        // initialize values for convergence checks
        scalar eqnResidual = 1, maxResidual = 0;
        scalar convergenceCriterion = 0;
        int maxNumIter = 0;
        // read convergence criterion from  system/fvSolution section SIMPLE
        simple.readIfPresent("convergence", convergenceCriterion);
        simple.readIfPresent("maxIter", maxNumIter);

        // @need_description
        p.storePrevIter();

        // Pressure-velocity SIMPLE corrector
        #include "UEqn.H"
        #include "pEqn.H"

        // check convergence
        // TODO introduce check if maxNumIter==0
        if( maxResidual < convergenceCriterion || counter >= maxNumIter ){
          Info<< nl <<" Convergence info: maxResidual: "<< maxResidual
                  <<"  convergence criterion: "<< convergenceCriterion <<nl
                  << "Stokes flow: ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                  << " converged in "<< counter<<" steps"
                  << "  ClockTime = " << runTime.elapsedClockTime() << " s"<< nl << endl;
          if(counter >= maxNumIter){
            Info<<" dissolFoam Runtime WARNING:"
                    " steady state Stokes flow solver did not converge.\n"
                << "It reached maximum number of iterations"
                    << "  counter: " << counter
                    << "  maxNumIter: " << maxNumIter
                    << nl << endl;
          }
          break;
        }
        else{
          // @TODO probably here we should switch relTol to 0
          Info<< nl << " Step "<< counter <<"   maxResidual: "<< maxResidual
                    <<" > "<< convergenceCriterion <<nl<<endl;
        }
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
      counter = 0; // set iteration counter to 0

      // Steady-state convection-diffusion solver main loop
      while ( true ){
        counter++;

        double residual = solve(  fvm::div(phi, C) == fvm::laplacian(D, C)  ).initialResidual();
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
      Info<< "Write data, after conv-diff" << nl << endl;
    }
    else{
      Info<< "dissolFoam: Skip the Stokes and the convection-diffusion solvers." << nl << endl;
    }
    runTime++;
    
    /*##############################################################################
    *   Mesh relaxation
    * ##############################################################################*/
    
//  Mesh update 1: 
    pointVectorField& pointDisplacement = const_cast<pointVectorField&>(
      mesh.objectRegistry::lookupObject<pointVectorField>( "pointDisplacement" )
    );

//  Mesh update 2: Calculate new displacements (interpolate C field to points)  

    // interpolate the concentration from cells to wall faces
    primitivePatchInterpolation patchInterpolator(mesh.boundaryMesh()[wallID]);
    scalarField pointC = patchInterpolator.faceToPointInterpolate(C.boundaryField()[wallID]);

    // exponential extrapolation of the concentration on the edge vertices
    scalarField& point_conc = pointC;
    mesh_rlx->fixEdgeConcentration( point_conc );
    
    /*
    // check interpolation
    pointField wwalls = mesh.boundaryMesh()[wallID].localPoints();
    autoPtr<interpolation<scalar> >interpolatorC(interpolation<scalar>::New( "cellPointFace",C ));
    meshSearch searchEngine(mesh, true);
    forAll(wwalls, iii){
      label cellI = searchEngine.findNearestCell( wwalls[iii] );
      label faceI = searchEngine.findNearestFace( wwalls[iii] );
      scalar intpFieldC = interpolatorC->interpolate(wwalls[iii], cellI, faceI);
      
      
      /*
      if( std::abs(point_conc[iii] - intpFieldC) > 0.0001)
      //if( wwalls[iii].z() < 0.0001)
        Info<< "coord[" << iii << "]=" << wwalls[iii]<<" : "
                << oldCField[iii] << "   :   "
                << pointC[iii] <<"  intrpl: "<<intpFieldC<<nl;
      * / 
      point_conc[iii] = intpFieldC;
    }
    */
    
    /*
    forAll(pointC, iii){
      Info<< pointC[iii] << nl;
    }
    exit(0);
     */
    
    
    
//  Mesh update 2.1: Adjust timestep to mesh resolution. UPD: timestep means number of 
    // time sub steps when we modify boundary using the same concentration field 
    // in order not to deform mesh too mush.
    
    //Pout<< nl << nl;
    
    //Pout << "min "<< min( mesh.boundaryMesh()[wallID].localPoints() ) <<nl;
    //Pout << "max "<< max( mesh.boundaryMesh()[wallID].localPoints() ) <<nl;
    /*
    forAll(mesh.boundaryMesh()[wallID].localPoints(), i){
      Pout << "coord["<<i<<"]= "<< mesh.boundaryMesh()[wallID].localPoints()[i]<<nl;
    }
    */
    
    //Pout<< "edges0: "<< mesh.edges()[0]<<nl;
    
    Info<< "Starting calculation of mesh move substeps"<<nl<<endl;
    
    //int nsubsteps_aux = setSubStep( mesh.edges(), mesh.points(), point_conc );
    
    //int nsubsteps = static_cast<int>( nsubsteps_aux / 1.0 );
    int nsubsteps = 1;
    
    //double nsubsteps_inv = 1.0 / double(nsubsteps);
    double nsubsteps_inv = 1.0 / 1.0;
    
    scalarField pointCsub = pointC * nsubsteps_inv;
    
    Info<< nl << "Number of substeps  "<< nsubsteps << "   factor: "<< nsubsteps_inv << nl << endl;
    
    for(int subi = 0; subi < nsubsteps; subi++)
    {
      Info<< nl << "Substep nr. "<< (subi+1) << "  of total "<< nsubsteps << endl;
      
      vectorField pointNormY = mesh.boundaryMesh()[wallID].pointNormals();
      
      /*
      forAll(pointNormY, iii){
        pointNormY[iii].x() = 0.0;
        pointNormY[iii].z() = 0.0;
      }
       */
      
      
      //vectorField pointDispWall = pointCsub * mesh.boundaryMesh()[wallID].pointNormals(); //*runTime.deltaTValue()
      vectorField pointDispWall = pointCsub * pointNormY * runTime.deltaTValue();
      
      vectorField& pointDispWall_ = pointDispWall;
      
      mesh_rlx->fixWallDisplPeriodic( pointDispWall_ );
      

//  Mesh update 2.2: Inlet displacement
      vectorField& pdw = pointDispWall;
      vectorField pointDispInlet = mesh_rlx->calculateInletDisplacement(pdw);
      
//  Mesh update 4: Update boundary and relax interior mesh
      Info<< nl << "Update boundary and relax interior mesh" << endl;
      
      vectorField &wallDisp = refCast<vectorField>(pointDisplacement.boundaryField()[wallID]);
      pointDisplacement.boundaryField()[wallID] == wallDisp + pointDispWall;

      Info<< nl << "Wall displacement" << endl;
      mesh.update();
      
      vectorField &inletDisp = refCast<vectorField>(pointDisplacement.boundaryField()[inletID]);
      pointDisplacement.boundaryField()[inletID] == inletDisp + pointDispInlet;

      Info<< nl << "Inlet displacement" << endl;
      mesh.update();
      
      //runTime++;
      //runTime.write();
      
      Info<< nl << "Inlet rlx" << endl;
      partialSlipModPointPatchVectorField &p_dispPatch = 
                dynamic_cast< partialSlipModPointPatchVectorField& >
                (pointDisplacement.boundaryField()[inletID]);
      p_dispPatch.valueFraction() = false;
      for(int ii=0; ii<5; ii++){  
        mesh.update();
        
      //runTime++;
      //runTime.write();
        
      }
      Info<< nl << "Back to patch fixedValue" << endl;
      p_dispPatch.valueFraction() = true;
      mesh.update();
      
      //runTime++;
      //runTime.write();
      
//  Mesh update 5: boundary mesh relaxation

      Info<< nl << "Boundary mesh relaxation" << nl << endl;
      
      for(int i=0; i<20; i++){
        vectorField boundaryRelax = mesh_rlx->wallRelaxation();
        
        vectorField& boundaryRelax_ = boundaryRelax;
        mesh_rlx->fixWallDisplPeriodic( boundaryRelax_ );

        vectorField &wallDispRelx = refCast<vectorField>(pointDisplacement.boundaryField()[wallID]);
        pointDisplacement.boundaryField()[wallID] == wallDispRelx + boundaryRelax;
        
        //Info<<"Before mesh update"<<nl;
        mesh.update();
      //runTime++;
      //runTime.write();
        //Info<<"After mesh update"<<nl;
        
      //mesh_rlx->boundaryCheck();
      
      }
      
      
      Info<< "Mesh update: ExecutionTime = " << runTime.elapsedCpuTime()
            << " s" << "  ClockTime = " << runTime.elapsedClockTime()
            << " s" << nl << endl;
      
      //runTime++;
      //runTime.write();
    }
  
    run0timestep = true;
    Info<< "Process: Time = " << runTime.timeName() << nl << endl;
  }

  Info << "End" << nl << endl;
  return 0;
}


// ************************************************************************* //
