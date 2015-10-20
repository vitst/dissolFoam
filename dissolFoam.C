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
#include "DissolMeshRlx.H"
#include "FieldOperations.H"

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
  bool dissolDebug(dissolProperties.lookupOrDefault<bool>("dissolDebug", false) );
  
  // if true the concentration for the points at the edge between wall and inlet
  // will be modified
  bool fixInletConcentration(readBool(dissolProperties.lookup("fixInletConcentration")));
  
  // if true it switches on the convection term in Navier-Stokes eqn
  bool NStokesInertia(dissolProperties.lookupOrDefault<bool>("NStokesInertia", false));
  
  // l_T=D/(k*h_0)
  scalar l_T( dissolProperties.lookupOrDefault<scalar>("lT", 1.0) );
  
  // Reynolds number
  scalar Re( dissolProperties.lookupOrDefault<scalar>("Re", 1.0/nu.value()) );
  
  // a tolerance for the relaxation cycles
  scalar rlxTol( dissolProperties.lookupOrDefault<scalar>("relaxationTolerance", 0.1) );
  
  // if true the grading in Z direction will change with time in accordance to the
  // formula G = inigradingZ/(timeCoef*t+1)
  bool varG( dissolProperties.lookupOrDefault<bool>("varG", false));
  scalar inigradingZ( dissolProperties.lookupOrDefault<scalar>("inigradingZ", 1.0) );
  scalar timeCoefZ( dissolProperties.lookupOrDefault<scalar>("timeCoefZ", 1.0) );
  int Nz( dissolProperties.lookupOrDefault<int>("numberOfCellsZ", 10) );
  
  
  Info << "*****************************************************************"<<nl;
  Info << "dissolFoamDict, fixInletConcentration:  " << fixInletConcentration <<nl;
  Info << "dissolFoamDict, rlxTol:  " << rlxTol <<nl;
  Info << "dissolFoamDict, inigradingZ:  " << inigradingZ <<nl;
  Info << "dissolFoamDict, timeCoefZ:  " << timeCoefZ <<nl;
  Info << "dissolFoamDict, numberOfCellsZ:  " << Nz <<nl;
  
  Info << "dissolFoamDict, NStokesInertia:  " << NStokesInertia <<nl;
  Info << "dissolFoamDict, Re:  " << Re <<nl;
  Info << "*****************************************************************"<<nl;
  
  // Get patch ID for boundaries we want to move ("walls" "inlet")
  label wallID  = mesh.boundaryMesh().findPatchID("walls");
  label inletID = mesh.boundaryMesh().findPatchID("inlet");
  label outletID = mesh.boundaryMesh().findPatchID("outlet");
  
  Info<< "Setup mesh relaxation class" << endl;
  DissolMeshRlx* mesh_rlx = new DissolMeshRlx(mesh, varG); // pointer to the mesh relaxation object
  Info<< "Setup field operation class" << endl;
  FieldOperations* fieldO = new FieldOperations(); // pointer to the mesh relaxation object
  
  // calculating initial area of the inlet in order to scale U later
  scalar areaCoef0 = fieldO->getInletArea(mesh);
  
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
  /*
   * A variable run0timestep is used to skip the Navier-Stokes and 
   * the convection-diffusion solvers if current time is not 0.
   * At the end of each cycle it is set to true.
   */
  bool run0timestep = false;
  if( runTime.value() == 0 ){  run0timestep = true; }
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  vectorFieldList wallWeightsV;
  if( !varG ){
    wallWeightsV = mesh_rlx->calc_weights2( mesh.boundaryMesh()[wallID]);
  }
    
  vectorFieldList inletWeightsV = mesh_rlx->calc_edge_weights( mesh.boundaryMesh()[inletID]);
  vectorFieldList outletWeightsV = mesh_rlx->calc_edge_weights( mesh.boundaryMesh()[outletID]);
  vectorFieldList& iwv = inletWeightsV;
  vectorFieldList& owv = outletWeightsV;
  
  
  /*
  forAll( mesh.boundaryMesh(), j1){
    Info<< mesh.boundaryMesh()[j1].type() <<nl;
  }
  */
  
  
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
      while ( simple.loop() ){
        if(NStokesInertia){
          #include "UEqnNavierStokesConv.H"
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
      scalar nU = fieldO->getConstFlowRateFactor(mesh, U, areaCoef0);
      U==U/nU;
      
// *********************************************************
// *    Danckwerts boundary condition loop, Convenction-Diffusion
// *********************************************************
      Info << "Steady-state convection-diffusion"<< endl;
      int inlet_count = 0;
      while ( true ){
        #include "ConvectionDiffusion.H"
        scalar tttol = fieldO->setConstFlowRateFactor(mesh, U, C, inletID, D.value());
        inlet_count+=1;
        Info<<"Inlet C iter "<<inlet_count<<"  tolerance: "<<tttol<< nl;
        if(tttol<convCritCD) break;
      }
      
// *********************************************************
// *    Write Output data
// *********************************************************
      (runTime.value()==0) ? runTime.writeNow() : runTime.write();
      Info<< "Write data, after conv-diff" << nl << nl;
    }
    else{
      Info<< "dissolFoam: Skip the Stokes and the convection-diffusion solvers."<<nl<<nl;
    }
    runTime++;
    
// *********************************************************
// *    Mesh motion & relaxation
// *********************************************************
    pointVectorField& pointVelocity = const_cast<pointVectorField&>(
      mesh.objectRegistry::lookupObject<pointVectorField>( "pointMotionU" )
    );

//  Mesh update 1: Calculate new displacements (interpolate C field to points)  
    vectorField pointDispWall = fieldO->getWallPointMotion(mesh, C, l_T, wallID);
    vectorField& pdw = pointDispWall;
    if(fixInletConcentration){
      Info << "Fix concentration on the edge between walls and inlet"<< nl;
      mesh_rlx->fixEdgeConcentration(pdw);
    }
    
//  Mesh update 2: boundary mesh relaxation
    if( varG ){
      Info<<nl<<"Calculating new Z grading...."<<nl;
      scalar Gz = inigradingZ / (timeCoefZ * runTime.value() + 1.0);
      scalar lambdaZ = 1/static_cast<double>(Nz-1) * std::log( Gz );
      wallWeightsV = mesh_rlx->calc_weights( mesh.boundaryMesh()[wallID], lambdaZ);
    }
    const vectorFieldList& wWv = wallWeightsV;
    
    pointField savedPointsAll = mesh.points();
    vectorField pointDispInlet = mesh_rlx->calculateInletDisplacement(pdw);
    vectorField pointDispOutlet = mesh_rlx->calculateOutletDisplacement(pdw);
    
//  Mesh update 3: boundary mesh relaxation
    
    const vectorField& wR = pointDispWall;
    pointField mpW = mesh_rlx->doWallDisplacement( wR * runTime.deltaTValue() );
    mesh.movePoints( mpW);
    
    const vectorField& iR = pointDispInlet;
    pointField mpI = mesh_rlx->doInletDisplacement( iR * runTime.deltaTValue() );
    mesh.movePoints( mpI );
    
    const vectorField& oR = pointDispOutlet;
    pointField mpO = mesh_rlx->doOutletDisplacement( oR * runTime.deltaTValue() );
    mesh.movePoints( mpO );
    
    if(dissolDebug){
      runTime.write();
      runTime++;
    }

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
    
    // if it is Debug mode finalize the run here
    if( dissolDebug ){
      runTime.write();
      runTime++;
      return 0;
    }
  
    run0timestep = true;
    Info<< "Process: Time = " << runTime.timeName() << nl << nl;
  }

  Info << "End" << nl;
  return 0;
}

// **************************** End of the solver ******************************** //
