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

#include <gsl/gsl_sf_hyperg.h>

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
#include "steadyStateControl.H"
//#include "primitivePatchInterpolationSync.H"
#include "coupledPatchInterpolation.H"
#include "DissolMeshRlx.H"


// mesh search
#include "interpolation.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "meshSearch.H"

// interpolation and tables
#include "interpolationTable.H"

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

  bool dissolDebug
  (
    dissolProperties.lookupOrDefault<bool>("dissolDebug", false)
  );
  
  bool fixInletConcentration
  (
    readBool( dissolProperties.lookup("fixInletConcentration") )
  );

  
  word newInletConcentration
  (
    dissolProperties.lookupOrDefault<word>("newInletConcentration", "none")
  );
  
  /*
  bool cInletHypergeometric
  (
    dissolProperties.lookupOrDefault<bool>("hypergeometric", false)
  );
  bool modifyInletC
  (
    dissolProperties.lookupOrDefault<bool>("modifyInletC", false)
  );
  */
  
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
  
    scalar Gy = inigradingY / (timeCoefY * runTime.value() + 1.0);
    scalar lambdaY = 1/static_cast<double>(Ny) * std::log( Gy );
    
    scalarListList inletWeights = mesh_rlx->calc_weights1( mesh.boundaryMesh()[inletID], lambdaY, 2);
    const scalarListList& iW = inletWeights;
    
    /*
    scalar Gz = inigradingZ / (timeCoefZ * runTime.value() + 1.0);
    scalar lambdaZ = 1/static_cast<double>(Nz) * std::log( Gz );
    scalarListList wallWeights = mesh_rlx->calc_weights0( mesh.boundaryMesh()[wallID], lambdaZ, 3);
    //scalarListList wallWeights = mesh_rlx->calc_weights( mesh.boundaryMesh()[wallID], lambdaZ, 3);
    const scalarListList& wW = wallWeights;
    */

	const surfaceScalarField& magSf1 = mesh.magSf();
	scalar areaCoef0 = 0.0;
        forAll(magSf1.boundaryField()[inletID], facei){
          areaCoef0 += magSf1.boundaryField()[inletID][facei];
        }
	reduce(areaCoef0, maxOp<scalar>());
 
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
      steadyStateControl simple(mesh);
      while ( simple.loop() ){
        // Pressure-velocity SIMPLE corrector
        #include "UEqn.H"
        #include "pEqn.H"
        //turbulence->correct();
      }
      
      // ***************************************************************************
      if(newInletConcentration=="parabolic"||newInletConcentration=="hypergeometric" ){
        
        labelHashSet includePatches(1);
        includePatches.insert(wallID);
        triSurface wallTriSurface
        (
          triSurfaceTools::triangulate( mesh.boundaryMesh(), includePatches )
        );
        const triSurfaceSearch querySurf(wallTriSurface);
        const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();
        bool pm = false, mp = false;

        pointField pointFaceCntr = mesh.boundaryMesh()[inletID].faceCentres();
        scalarField Cinlet( pointFaceCntr.size() );
        
        interpolationTable<scalar> sherwoodNumTab ("Sh.dat");

        scalar lp = 0.5;

        forAll(Cinlet, i){
          point p = pointFaceCntr[i];
          point searchStart = p;
          searchStart.y() = 0.0;
          searchStart.z() += 0.01;
          point searchEnd = searchStart;
          searchStart.y() = 500.0;

          point maxY = p;
          pointIndexHit pHit = tree.findLine(searchStart, searchEnd);
          if ( pHit.hit() )
          {
            maxY =  pHit.hitPoint();
            pm = true;
          }
          else{
          }

          point minY = p;
          searchStart.y() = -500.0;
          pHit = tree.findLine(searchStart, searchEnd);
          if ( pHit.hit() )
          {
            minY = pHit.hitPoint();
            mp = true;
          }
          else{
          }

          if(pm && !mp){
            Info<<"WARNING! The ray did not find the surface from + direction"<<nl;
          }
          else if(!pm && mp){
            Info<<"WARNING! The ray did not find the surface from - direction"<<nl;
          }
          else if(!pm && !mp){
            Info<<"WARNING! The ray did not find the surface from both direction"<<nl;
          }

          scalar h = mag(maxY-minY);
          if(h == 0.0){
            Info<<"h=0  "<< p << "   "<< min(mesh.boundaryMesh()[wallID].localPoints().component(vector::Z))<<nl;
          }
        
          if( newInletConcentration == "hypergeometric" ){
            // the concentration profile based on Hypergeometric function
            //scalar Sh = 8.0;
            scalar G = 2*h/lp;
            scalar Sh = sherwoodNumTab( Foam::log10(G) );
            scalar lam = Sh / (G+Sh);
            scalar r  = Foam::sqrt( 3*G*lam/8 );
            scalar yh = (p.y() + h/2.) / h;
            scalar gz = gsl_sf_hyperg_1F1( 0.25*(1-r), 0.5, r*(2*yh-1)*(2*yh-1) ) * Foam::exp(-2 * r *yh*(yh-1));
            yh = 0.0;
            scalar gz0 = gsl_sf_hyperg_1F1( 0.25*(1-r), 0.5, r*(2*yh-1)*(2*yh-1) ) * Foam::exp(-2 * r *yh*(yh-1));
            scalar fA = lam / gz0;

            Cinlet[i] = fA * gz;
          }
          else{
            // parabolic function
            scalar a=8.0/(h*h+4*h*lp);
            Cinlet[i] = 1-a/2.0 * p.y()*p.y();
            scalar G = 2*h/lp; // = 2kh/D
            scalar alpha = 4*G/(8+G);

            Cinlet[i] /= 1-alpha/20.0;
          }

        }

        C.boundaryField()[inletID] == Cinlet;
      }
        
      if(newInletConcentration=="distanceFunction" ){
        vectorField ndist = mesh_rlx->normalsOnTheEdge();
        scalarField ddist = mesh_rlx->distanceToTheEdge(ndist);
        
	const surfaceScalarField& magSf = mesh.magSf();
	scalar areaCoef = 0.0;
        forAll(ddist, facei){
          areaCoef += magSf.boundaryField()[inletID][facei];
        }
	//areaCoef = Foam::sqrt(areaCoef);
	areaCoef /= areaCoef0;
	//areaCoef = 1.0 ;

        //scalar lp = 1.025;
        scalar lp = 0.5;
        scalarField newC(ddist.size(), 0.0);
        forAll(newC, ii){
          //newC[ii] = Foam::exp(ddist[ii]/lp/areaCoef);
          //newC[ii] = 1 + 1-Foam::exp(-ddist[ii]/lp/areaCoef);
          newC[ii] = 1.0/areaCoef + 1-Foam::exp(-ddist[ii]/lp/areaCoef);
	  //Info<<ii<<"  "<<ddist[ii]<<"  "<<newC[ii]<<"  "<<nl;
        }
	//Info<<"areaCoef:   "<<areaCoef<<nl;
	//std::exit(0); 
        scalar totArea = 0.0;
        scalar totCArea = 0.0;
        forAll(newC, facei){
          totArea += phi.boundaryField()[inletID][facei];
          totCArea += newC[facei] * phi.boundaryField()[inletID][facei];
        }
        
        scalar corrF = 0.0;
        if(totCArea!=0.0){
          corrF = totArea/totCArea;
        }
        
        C.boundaryField()[inletID]==corrF*newC;
      }

      // ***************************************************************************
      
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
    coupledPatchInterpolation patchInterpolator( mesh.boundaryMesh()[wallID], mesh );
    
    scalarField pointCface = -C.boundaryField()[wallID].snGrad();
    vectorField pointNface = mesh.boundaryMesh()[wallID].faceNormals();
    
    /*
    vectorField motionVec = pointCface * pointNface;
    vectorField pointDispWall = patchInterpolator.faceToPointInterpolate(motionVec);
    */
    
    
    scalarField motionC = patchInterpolator.faceToPointInterpolate(pointCface);
    vectorField motionN = patchInterpolator.faceToPointInterpolate(pointNface);
    forAll(motionN, ii){
      motionN[ii]/=mag(motionN[ii]);
    }
    vectorField pointDispWall = motionC*motionN;
    
    
    vectorField& pdw = pointDispWall;
    
    //mesh_rlx->fixEdgeConcentration();
    //std::exit(0);
    
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
    //scalarListList wallWeights = mesh_rlx->calc_weights0( mesh.boundaryMesh()[wallID], lambdaZ, 3);
    scalarListList wallWeights = mesh_rlx->calc_weights( mesh.boundaryMesh()[wallID], lambdaZ, 3);
    const scalarListList& wW = wallWeights;
    
    /*
    scalar Gy = inigradingY / (timeCoefY * runTime.value() + 1.0);
    scalar lambdaY = 1/static_cast<double>(Ny) * std::log( Gy );
    
    vectorFieldList inletWeights = mesh_rlx->calc_weights2( mesh.boundaryMesh()[inletID], lambdaY, 2);
    const vectorFieldList& iW = inletWeights;
    */
    
    /*
    forAll(inletWeights, ii){
      Info<< ii << "  " << inletWeights[ii]<<nl;
    }
    std::exit(0);
     */
    
    //vectorFieldList outletWeights = mesh_rlx->calc_weights2( mesh.boundaryMesh()[outletID], lambdaY, 2);
    //const vectorFieldList& oW = outletWeights;
    scalarListList outletWeights = mesh_rlx->calc_weights1( mesh.boundaryMesh()[outletID], lambdaY, 2);
    const scalarListList& oW = outletWeights;
    
    Info<<nl<<"Boundary mesh relaxation. Z grading is "<< Gz
            <<"   Y grading is "<< Gy <<nl<<nl;
    
    pointField savedPointsAll = mesh.points();
    
    const vectorField& wR = pointDispWall;
    pointField mpW = mesh_rlx->doWallDisplacement( wR * runTime.deltaTValue() );
    mesh.movePoints( mpW);
    
    const vectorField& iR = pointDispInlet;
    pointField mpI = mesh_rlx->doInletDisplacement( iR * runTime.deltaTValue() );
    mesh.movePoints( mpI );
    
    const vectorField& oR = pointDispOutlet;
    pointField mpO = mesh_rlx->doOutletDisplacement( oR * runTime.deltaTValue() );
    mesh.movePoints( mpO );
    
    if( dissolDebug ){
      runTime.write();
      runTime++;
    }
    
    
    Info<<"Relaxing inlet-wall edge..."<<nl;
    vectorField wiEdgeRlx = mesh_rlx->edgeRelaxation( mesh.boundaryMesh()[inletID], rlxTol);
    const vectorField& werI = wiEdgeRlx;
    pointField mpWIE = mesh_rlx->doWallDisplacement( werI );
    mesh.movePoints( mpWIE );
    
    Info<<"Relaxing outlet-wall edge..."<<nl;
    vectorField woEdgeRlx = mesh_rlx->edgeRelaxation( mesh.boundaryMesh()[outletID], rlxTol);
    const vectorField& werO = woEdgeRlx;
    pointField mpWOE = mesh_rlx->doWallDisplacement( werO );
    mesh.movePoints( mpWOE );
    

    Info<<"Relaxing wall... time: "<< runTime.cpuTimeIncrement() <<nl;
    //vectorField wallRelax = mesh_rlx->wallRelaxation12( mesh.boundaryMesh()[wallID], wW, rlxTol);
    vectorField wallRelax = mesh_rlx->wallRelaxation11( mesh.boundaryMesh()[wallID], wW, rlxTol);
    Info<<"Wall relaxation time: " << runTime.cpuTimeIncrement() << " s" << nl;
    
    mesh.movePoints( savedPointsAll );
    const vectorField wR1 = wallRelax/runTime.deltaTValue() 
            + wiEdgeRlx/runTime.deltaTValue()
            + woEdgeRlx/runTime.deltaTValue()
            + pointDispWall;
    
    vectorField vvff1 = wR1 * runTime.deltaTValue();
    const vectorField &vvff = vvff1;
    
    
    Info<<"Relaxing inlet..."<<nl;
    vectorField inlRelax = mesh_rlx->wallRelaxation42( mesh.boundaryMesh()[inletID], iW, rlxTol, vvff);
    Info<<"Inlet relaxation time: " << runTime.cpuTimeIncrement() << " s" << nl;

    Info<<"Relaxing outlet..."<<nl;
    vectorField outRelax = mesh_rlx->wallRelaxation42( mesh.boundaryMesh()[outletID], oW, rlxTol, vvff);
    Info<<"Outlet relaxation time: " << runTime.cpuTimeIncrement() << " s" << nl;
    
    /*
    pointField mpW1 = mesh_rlx->doWallDisplacement(wR1 * runTime.deltaTValue() );
    mesh.movePoints( mpW1 );
    runTime.write();
    runTime++;
    
    const vectorField iR1 = inlRelax/runTime.deltaTValue(); // + pointDispInlet;
    pointField mpI1 = mesh_rlx->doInletDisplacement(iR1 * runTime.deltaTValue());
    mesh.movePoints( mpI1 );
    runTime.write();
    runTime++;
     */
    /*
    const vectorField oR1 = outRelax/runTime.deltaTValue() + pointDispOutlet;
    pointField mpO1 = mesh_rlx->doOutletDisplacement(oR1 * runTime.deltaTValue());
    mesh.movePoints( mpO1 );
    runTime.write();
    runTime++;
     */
    //std::exit(0);
    
    
    mesh.movePoints( savedPointsAll );

    wiEdgeRlx /= runTime.deltaTValue();
    woEdgeRlx /= runTime.deltaTValue();
    wallRelax /= runTime.deltaTValue();
    pointVelocity.boundaryField()[wallID] == wallRelax + pointDispWall + wiEdgeRlx + woEdgeRlx;

    inlRelax /= runTime.deltaTValue();
    //pointVelocity.boundaryField()[inletID] == inlRelax; // + pointDispInlet;
    pointVelocity.boundaryField()[inletID] == inlRelax + pointDispInlet;
    
    outRelax /= runTime.deltaTValue();
    pointVelocity.boundaryField()[outletID] == outRelax + pointDispOutlet;
    
    mesh.update();

    Info<< "Mesh update: ExecutionTime = " << runTime.elapsedCpuTime()
          << " s" << "  ClockTime = " << runTime.elapsedClockTime()
          << " s" << nl << nl;
    
    if( dissolDebug ){
      runTime.write();
      runTime++;
      std::exit(0);
    }
  
    run0timestep = true;
    Info<< "Process: Time = " << runTime.timeName() << nl << nl;
  }

  Info << "End" << nl;
  return 0;
}

// **************************** End of the solver ******************************** //
