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
#include "RASModel.H"

// OFe
#include "pointPatchField.H"
#include "dynamicFvMesh.H"
#include "syncTools.H"

// dissol prj
#include "steadyStateControl.H"
#include "coupledPatchInterpolation.H"
#include "DissolMeshRlx.H"

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
  #include "createFields.H"

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

  // if true the solver will be stopped after a single iteration, writing to the disk
  // three times: overwriting 0, mesh with moved surface, relaxed mesh.
  bool dissolDebug
  (
    dissolProperties.lookupOrDefault<bool>("dissolDebug", false)
  );
  
  // if true the concentration for the points at the edge between wall and inlet
  // will be modified
  bool fixInletConcentration
  (
    readBool( dissolProperties.lookup("fixInletConcentration") )
  );

  // defines the inlet boundary for the C field
  word newInletConcentration
  (
    dissolProperties.lookupOrDefault<word>("newInletConcentration", "none")
  );
  
  // if true it switches on the convection term in Navier-Stokes eqn
  bool NavierStokesConvection
  (
    dissolProperties.lookupOrDefault<bool>("NavierStokesConvection", false)
  );
  
  // l_T=D/(k*h_0)
  scalar l_T( dissolProperties.lookupOrDefault<scalar>("lT", 1.0) );
  
  // Reynolds number
  scalar Re( dissolProperties.lookupOrDefault<scalar>("Re", 1.0/nu.value()) );
  
  // a tolerance for the relaxation cycles
  scalar rlxTol( dissolProperties.lookupOrDefault<scalar>("relaxationTolerance", 0.1) );
  
  // if true the grading in Z direction will change with time in accordance to the
  // formula G = inigradingZ/(timeCoef*t+1)
  bool varG
  (
    dissolProperties.lookupOrDefault<bool>("varG", false)
  );
  scalar inigradingZ( dissolProperties.lookupOrDefault<scalar>("inigradingZ", 1.0) );
  scalar timeCoefZ( dissolProperties.lookupOrDefault<scalar>("timeCoefZ", 1.0) );
  int Nz( dissolProperties.lookupOrDefault<int>("numberOfCellsZ", 10) );
  
  
  Info << "*****************************************************************"<<nl;
  Info << "dissolFoamDict, fixInletConcentration:  " << fixInletConcentration <<nl;
  Info << "dissolFoamDict, inlet concentration based on:  " << newInletConcentration <<nl;
  Info << "dissolFoamDict, rlxTol:  " << rlxTol <<nl;
  Info << "dissolFoamDict, inigradingZ:  " << inigradingZ <<nl;
  Info << "dissolFoamDict, timeCoefZ:  " << timeCoefZ <<nl;
  Info << "dissolFoamDict, numberOfCellsZ:  " << Nz <<nl;
  
  Info << "dissolFoamDict, NavierStokesConvection:  " << NavierStokesConvection <<nl;
  Info << "dissolFoamDict, Re:  " << Re <<nl;
  Info << "*****************************************************************"<<nl;
  
  // Get patch ID for boundaries we want to move ("walls" "inlet")
  label wallID  = mesh.boundaryMesh().findPatchID("walls");
  label inletID = mesh.boundaryMesh().findPatchID("inlet");
  label outletID = mesh.boundaryMesh().findPatchID("outlet");
  
  Info<< "Setup mesh relaxation class" << endl;
  // pointer to the mesh relaxation object
  DissolMeshRlx* mesh_rlx = new DissolMeshRlx(mesh);
  
  // calculating initial area of the inlet in order to scale U later
  scalar areaCoef0 = 0.0;
  const surfaceScalarField& magSf2 = mesh.magSf();
  forAll(magSf2.boundaryField()[inletID], facei){
    areaCoef0 += magSf2.boundaryField()[inletID][facei];
  }
  reduce(areaCoef0, sumOp<scalar>());
  
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
  /*
   * run0timestep is used to skip the Navier-Stokes and the convection-diffusion solvers
   * if timestep is not 0. At the end of each cycle it is set to true.
   */
  bool run0timestep = false;
  if( runTime.value() == 0 ){  run0timestep = true; }
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  //const cellList & cc = mesh.cells();
  //Info<<cc<<nl; 
  //std::exit(0);
  
  vectorFieldList wallWeightsV;
  if( !varG ){
    wallWeightsV = mesh_rlx->calc_weights2( mesh.boundaryMesh()[wallID]);
  }
  const vectorFieldList& wWv = wallWeightsV;
  
  /*
  pointField ppp = mesh.boundaryMesh()[wallID].localPoints();
  forAll(ppp, iii){
    if(ppp[iii].z()<2.0)
      Info<<wWv[iii]<<nl; 
  }
  std::exit(0);
  */
  
  
  vectorFieldList inletWeightsV = mesh_rlx->calc_edge_weights( mesh.boundaryMesh()[inletID]);
  vectorFieldList outletWeightsV = mesh_rlx->calc_edge_weights( mesh.boundaryMesh()[outletID]);
  vectorFieldList& iwv = inletWeightsV;
  vectorFieldList& owv = outletWeightsV;

  // MAIN TIME LOOP //
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
        if(NavierStokesConvection){
          #include "UEqnNavierStokesConv.H"
          #include "pEqn.H"
          turbulence->correct();
        }
        else{
          #include "UEqnStokes.H"
          #include "pEqn.H"
        }
        //#include "pEqn.H"
        //turbulence->correct();
      }
      
      /*##########################################
       *   Keeping flow rate constant
       *##########################################*/
      scalar magU = mag( fvc::domainIntegrate( U ).value() );
      
      vectorField po = mesh.points();
      scalar maxZ = max( po.component(vector::Z) );
      reduce(maxZ, maxOp<scalar>());
      scalar minZ = min( po.component(vector::Z) );
      reduce(minZ, minOp<scalar>());
      
      scalar nU = magU / (maxZ-minZ) / areaCoef0;
      
      U==U/nU;
      
      // ***************************************************************************
      // Calculating different boundary condition for the C field on the inlet
      if(newInletConcentration=="parabolic"||newInletConcentration=="hypergeometric" ){
        
        // if the "parabolic" or "hypergeometric" option, we need to know the aperture
        // of the fracture at the inlet
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
        coupledPatchInterpolation patchInterpolator(mesh.boundaryMesh()[inletID], mesh);
        scalarField phph = patchInterpolator.faceToPointInterpolate(phi.boundaryField()[inletID]);
        scalarField & phphr = phph;
        scalarField newConc = mesh_rlx->newC(phphr, 0.000000001);
        scalarField Cf = patchInterpolator.pointToFaceInterpolate(newConc);
        C.boundaryField()[inletID]==Cf;
      }
      
      // ***************************************************************************

      int inlet_count = 0;
      while ( true ){
      
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
      
        if(newInletConcentration != "Dankwerts"){
          break;
        }
        
        scalarField& oldC = C.boundaryField()[inletID];
        
        scalarField newC(oldC.size(), 0.0);
        
        const labelList& fc = mesh.boundaryMesh()[inletID].faceCells();
        const vectorField& fcr = mesh.boundaryMesh()[inletID].faceCentres();
        const vectorField& ccr = mesh.cellCentres();
        forAll(fc, ii){
          //newC[ii] = C[fc[ii]];
          vector vdel = fcr[ii]-ccr[fc[ii]];
          scalar del = std::abs(vdel.z());
          //scalar del = mag(vdel);
          
          scalar aa = D.value() / ( std::abs(U[fc[ii]].z()) * del);
          newC[ii] = (1+aa*C[fc[ii]])/(1+aa);
        }
        
        scalarField diff = newC - oldC;
        
        C.boundaryField()[inletID]==newC;
        
        scalar ttt = std::abs( gSum( diff ) );
        reduce(ttt, sumOp<scalar>());
        scalar tttn =std::abs( gSum( C.boundaryField()[inletID] ) );
        reduce(tttn, sumOp<scalar>());
        
        scalar tttol = 0.0;
        if( tttn!=0.0 ){
          tttol = ttt/tttn;
        }

        Info<<"Inlet C iter "<<inlet_count<<"  tolerance: "<<tttol<< nl;
        inlet_count+=1;
        
        if(tttol<convCritCD){
          break;
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
    *   Mesh motion & relaxation
    * ##############################################################################*/
    
    pointVectorField& pointVelocity = const_cast<pointVectorField&>(
      mesh.objectRegistry::lookupObject<pointVectorField>( "pointMotionU" )
    );

    // Mesh update: Calculate new displacements (interpolate C field to points)  

    // interpolate the concentration from cells to wall faces
    coupledPatchInterpolation patchInterpolator( mesh.boundaryMesh()[wallID], mesh );
    
    scalarField pointCface = -C.boundaryField()[wallID].snGrad();
    vectorField pointNface = mesh.boundaryMesh()[wallID].faceNormals();
    
    //vectorField motionVec = pointCface * pointNface;
    //vectorField pointDispWall = patchInterpolator.faceToPointInterpolate(motionVec);
    
    scalarField motionC = patchInterpolator.faceToPointInterpolate(pointCface);
    vectorField motionN = patchInterpolator.faceToPointInterpolate(pointNface);
    
    // normalize vertex normals to 1
    forAll(motionN, ii){
      motionN[ii]/=mag(motionN[ii]);
    }
    
    vectorField pointDispWall = l_T * motionC*motionN;
    vectorField& pdw = pointDispWall;
    
    
    if(fixInletConcentration){
      Info << "Fix concentration on the edge between walls and inlet"<< nl;
      mesh_rlx->fixEdgeConcentration(pdw);
      
      /*
      vectorField& motionNp = motionN;
      scalarField& inlC = C.boundaryField()[inletID];
      scalarField& walC = C.boundaryField()[wallID];
      mesh_rlx->fixEdgeConcentration1(pdw, inlC, motionNp, walC);
       */
    }
    
//  Mesh update 2.2: Inlet and outlet displacement
    //vectorField pointDispInlet = mesh_rlx->calculateInletDisplacement(pdw);
    //vectorField pointDispOutlet = mesh_rlx->calculateOutletDisplacement(pdw);
      
//  Mesh update 5: boundary mesh relaxation
    
    // *********************************************************************************
    // @TODO make it nice
    
    //scalarListList wallWeightsS;
    if( varG ){
      Info<<nl<<"Calculating new Z grading...."<<nl;

      scalar Gz = inigradingZ / (timeCoefZ * runTime.value() + 1.0);
      scalar lambdaZ = 1/static_cast<double>(Nz-1) * std::log( Gz );
      wallWeightsV = mesh_rlx->calc_weights( mesh.boundaryMesh()[wallID], lambdaZ);
      //scalarListList wallWeights = mesh_rlx->calc_weights_faces( mesh.boundaryMesh()[wallID]);
    }
    //const scalarListList& wWs = wallWeightsS;
    
    //Info<<nl<<"Boundary mesh relaxation. Z grading is "<< Gz <<nl<<nl;
    
    // *********************************************************************************
    
    pointField savedPointsAll = mesh.points();
    
    /*
    if( dissolDebug ){
      const vectorField& wR = pointDispWall;
      pointField mpW = mesh_rlx->doWallDisplacement( wR * runTime.deltaTValue() );
      mesh.movePoints( mpW);

      runTime.write();
      runTime++;
      
      mesh.movePoints( savedPointsAll );
    }
    */
    
    vectorField pointDispInlet = mesh_rlx->calculateInletDisplacement(pdw);
    vectorField pointDispOutlet = mesh_rlx->calculateOutletDisplacement(pdw);
    
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
    
    // *********************************************************************************
    // Relaxing edges. 1D
    Info<<"Relaxing the inlet-wall edge..."<<nl;
    vectorField wiEdgeRlx = mesh_rlx->edgeRelaxation( mesh.boundaryMesh()[inletID], rlxTol, iwv);
    const vectorField& werI = wiEdgeRlx;
    pointField mpWIE = mesh_rlx->doWallDisplacement( werI );
    mesh.movePoints( mpWIE );
    
    Info<<"Relaxing the outlet-wall edge..."<<nl;
    vectorField woEdgeRlx = mesh_rlx->edgeRelaxation( mesh.boundaryMesh()[outletID], rlxTol, owv);
    const vectorField& werO = woEdgeRlx;
    pointField mpWOE = mesh_rlx->doWallDisplacement( werO );
    mesh.movePoints( mpWOE );
    // *********************************************************************************
    

    // *********************************************************************************
    // Relaxing surfaces. 2D
    Info<<"Relaxing the wall... time: "<< runTime.cpuTimeIncrement() <<nl;
    vectorField wallRelax;
    if( varG ){
      wallRelax = mesh_rlx->wallRelaxation( mesh.boundaryMesh()[wallID], wWv, rlxTol);
    }
    else{
      wallRelax = mesh_rlx->wallRelaxation2( mesh.boundaryMesh()[wallID], wWv, rlxTol);
    }
    Info<<"Wall relaxation time: " << runTime.cpuTimeIncrement() << " s" << nl;
    
    mesh.movePoints( savedPointsAll );
    
    const vectorField wR1 = wallRelax/runTime.deltaTValue() 
            + wiEdgeRlx/runTime.deltaTValue()
            + woEdgeRlx/runTime.deltaTValue()
            + pointDispWall;
    
    vectorField vvff1 = wR1 * runTime.deltaTValue();
    const vectorField &vvff = vvff1;
    
    Info<<"Relaxing inlet..."<<nl;
    vectorField inlRelax = mesh_rlx->inletOutletRlx( mesh.boundaryMesh()[inletID], rlxTol, vvff);
    Info<<"Inlet relaxation time: " << runTime.cpuTimeIncrement() << " s" << nl;

    Info<<"Relaxing outlet..."<<nl;
    vectorField outRelax = mesh_rlx->inletOutletRlx( mesh.boundaryMesh()[outletID], rlxTol, vvff);
    Info<<"Outlet relaxation time: " << runTime.cpuTimeIncrement() << " s" << nl;
    
    mesh.movePoints( savedPointsAll );
    // *********************************************************************************
    

    // *********************************************************************************
    // Final mesh update. 3D
    wiEdgeRlx /= runTime.deltaTValue();
    woEdgeRlx /= runTime.deltaTValue();
    wallRelax /= runTime.deltaTValue();
    pointVelocity.boundaryField()[wallID] == wallRelax + pointDispWall + wiEdgeRlx + woEdgeRlx;

    inlRelax /= runTime.deltaTValue();
    pointVelocity.boundaryField()[inletID] == inlRelax + pointDispInlet;
    
    outRelax /= runTime.deltaTValue();
    pointVelocity.boundaryField()[outletID] == outRelax + pointDispOutlet;
    
    Info<<"Final mesh update"<<nl;
    
    mesh.update();
    // *********************************************************************************

    Info<<"Mesh update: ExecutionTime = "<<runTime.elapsedCpuTime()<<" s"
                  << "  ClockTime = " << runTime.elapsedClockTime()<<" s"<<nl<<nl;
    
    if( dissolDebug ){
      runTime.write();
      runTime++;
      return 0;  //std::exit(0);
    }
  
    run0timestep = true;
    Info<< "Process: Time = " << runTime.timeName() << nl << nl;
  }

  Info << "End" << nl;
  return 0;
}

// **************************** End of the solver ******************************** //
