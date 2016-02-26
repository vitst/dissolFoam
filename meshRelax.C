/* 
 * Surface relaxation for dissolFoam
 * Dec 09 2014
 * 
 */

#include "meshRelax.H"

#include "Pstream.H"
#include "pyramidPointFaceRef.H"
#include "triPointRef.H"

// mesh1 is the mesh at time 0
meshRelax::meshRelax(dynamicFvMesh& mesh, const argList& args)
:
  version("0.6"),
  date("Oct 2015"),
  mesh_(mesh)
{
  // get ID of each patch we need
  wallID   = mesh_.boundaryMesh().findPatchID("walls");
  inletID  = mesh_.boundaryMesh().findPatchID("inlet");
  outletID = mesh_.boundaryMesh().findPatchID("outlet");

  setUpLists();
  setUpPairsConc();
  
  Time& time = const_cast<Time&>(mesh_.time());
  deltaT = time.deltaTValue();
  
  // reading dissolFoam dictionary
  IOdictionary dissolProperties
  (
    IOobject
    (
      "dissolFoamDict",
      time.system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );
  
  // if true the solver will be stopped after a single iteration, writing to the disk
  // three times: overwriting 0, mesh with moved surface, relaxed mesh.
  dissolDebug = dissolProperties.lookupOrDefault<bool>("dissolDebug", false);
  
  // if true the concentration for the points at the edge between wall and inlet
  // will be modified
  fixInletWallEdgeDispl = readBool(dissolProperties.lookup("fixInletWallEdgeDispl"));
  // a tolerance for the relaxation cycles
  if( !dissolProperties.readIfPresent<scalar>("relaxationTolerance", rlxTol) ){
    SeriousErrorIn("meshRelax::meshRelax(dynamicFvMesh& mesh)")
            <<"There is no relaxationTolerance parameter in dissolProperties dictionary"
            <<exit(FatalError);
  }
  
  // if true the grading in Z direction will change with time in accordance to the
  // formula G = inigradingZ/(timeCoef*t+1)
  variableGrading = dissolProperties.lookupOrDefault<bool>("varG", false);
  inigradingZ = dissolProperties.lookupOrDefault<scalar>("inigradingZ", 1.0);
  timeCoefZ = dissolProperties.lookupOrDefault<scalar>("timeCoefZ", 1.0);
  Nz = dissolProperties.lookupOrDefault<int>("numberOfCellsZ", 10);
  
  // relaxation acceleration factors
  k_1 = dissolProperties.lookupOrDefault<scalar>("k_1", 1.0);
  k_2 = dissolProperties.lookupOrDefault<scalar>("k_2", 1.0);
  
  q_2 = dissolProperties.lookupOrDefault<int>("q_2", 1);
  q_norm_recalc = dissolProperties.lookupOrDefault<int>("q_norm_recalc", 1);
  
  k_1edge = dissolProperties.lookupOrDefault<scalar>("k_1edge", 1.0);
  k_2edge = dissolProperties.lookupOrDefault<scalar>("k_2edge", 1.0);
  q_2edge = dissolProperties.lookupOrDefault<int>("q_2edge", 1);
  q_edge_norm_recalc = dissolProperties.lookupOrDefault<int>("q_edge_norm_recalc", 1);
  
  Foam::Time timeTmp(Foam::Time::controlDictName, args);
  Foam::instantList timeDirs = Foam::timeSelector::select0(timeTmp, args);
  timeTmp.setTime(timeDirs[0], 0);

  // in case 0 time does not exist
  if( timeTmp.timeName()!="0" ){
    SeriousErrorIn("fieldOperations::getInletFlowRateT0")
            <<"There is no 0 time directory. Check your decomposition as well!"
            <<exit(FatalError);
  }
  
  Foam::fvMesh meshTmp
  (
      Foam::IOobject
      (
          Foam::fvMesh::defaultRegion,
          timeTmp.timeName(),
          timeTmp,
          Foam::IOobject::MUST_READ
      )
  );
  
  if( !variableGrading ){
    wallWeights = calc_weights2( meshTmp, meshTmp.boundaryMesh()[wallID]);
  }
  
  Info << "dissolFoamDict, dissolDebug:  " << dissolDebug << nl;
  Info << "dissolFoamDict, fixInletConcentration:  " 
       << fixInletWallEdgeDispl << nl;
  Info << "dissolFoamDict, rlxTol:  " << rlxTol << nl;
  Info << "dissolFoamDict, k_1:  " << k_1 << nl;
  Info << "dissolFoamDict, k_2:  " << k_2 << nl;
  Info << "dissolFoamDict, q_2:  " << q_2 << nl;
  Info << "dissolFoamDict, q_norm_recalc:  " << q_norm_recalc << nl;
  Info << "dissolFoamDict, k_1edge:  " << k_1edge << nl;
  Info << "dissolFoamDict, k_2edge:  " << k_2edge << nl;
  Info << "dissolFoamDict, q_2edge:  " << q_2edge << nl;
  Info << "dissolFoamDict, q_edge_norm_recalc:  " 
       << q_edge_norm_recalc << nl;
  
  if( variableGrading ){
    Info << "dissolFoamDict, inigradingZ:  " << inigradingZ << nl;
    Info << "dissolFoamDict, timeCoefZ:  " << timeCoefZ << nl;
    Info << "dissolFoamDict, numberOfCellsZ:  " << Nz << nl;
  }
  Info << "************************************************************"
       << nl << endl;
  
  // calculating weights for the edges inlet-wall, outlet-wall
  inletWallEdgeWeights  = calc_edge_weights( meshTmp, meshTmp.boundaryMesh()[inletID] );
  outletWallEdgeWeights = calc_edge_weights( meshTmp, meshTmp.boundaryMesh()[outletID]);
}

void meshRelax::meshUpdate(vectorField& pointDispWall, Time& time){
  pointVectorField& pointVelocity = const_cast<pointVectorField&>(
    mesh_.objectRegistry::lookupObject<pointVectorField>( "pointMotionU" )
  );
  
  if(fixInletWallEdgeDispl){
    Info << "Fix concentration on inlet-wall edge" << endl;
    fixIWEdgeDispl(pointDispWall);
  }
  
//  Mesh update 2: boundary mesh relaxation
  if( variableGrading ){
    Info << nl << "Calculating new Z grading...." << endl;
    scalar Gz = inigradingZ / (timeCoefZ * time.value() + 1.0);
    scalar lambdaZ = 1/static_cast<double>(Nz-1) * std::log( Gz );
    wallWeights = calc_weights( mesh_.boundaryMesh()[wallID], lambdaZ);
  }

  pointField savedPointsAll = mesh_.points();
  vectorField pointDispInlet = calculateInletDisplacement(pointDispWall);
  vectorField pointDispOutlet = calculateOutletDisplacement(pointDispWall);

//  Mesh update 3: boundary mesh relaxation
  pointField mpW = doWallDisplacement( pointDispWall * deltaT );
  mesh_.movePoints( mpW);

  pointField mpI = doInletDisplacement( pointDispInlet * deltaT );
  mesh_.movePoints( mpI );

  pointField mpO = doOutletDisplacement( pointDispOutlet * deltaT );
  mesh_.movePoints( mpO );

  if(dissolDebug){
    time.write();
    time++;
  }

// Relaxing edges. 1D
  Info << "Relaxing inlet-wall edge..." << endl;
  vectorField wiEdgeRlx = edgeRelaxation( mesh_.boundaryMesh()[inletID], inletWallEdgeWeights);
  pointField mpWIE = doWallDisplacement( wiEdgeRlx );
  mesh_.movePoints( mpWIE );

  Info << "Relaxing outlet-wall edge..." << endl;
  vectorField woEdgeRlx = edgeRelaxation( mesh_.boundaryMesh()[outletID], outletWallEdgeWeights);
  pointField mpWOE = doWallDisplacement( woEdgeRlx );
  mesh_.movePoints( mpWOE );
  Info << "Edge relaxation cpuTime: "
       << time.cpuTimeIncrement() << " s" << endl;

  if(dissolDebug){
    time.write();
    time++;
  }
// *********************************************************************
// Relaxing surfaces. 2D

  Info << "Relaxing the wall..." << endl;
  vectorField wallRelax = wallRelaxation(mesh_.boundaryMesh()[wallID], wallWeights);
  Info << "Wall relaxation cpuTime: "
       << time.cpuTimeIncrement() << " s" << endl;

  mesh_.movePoints( savedPointsAll );

  const vectorField vvff =  wallRelax 
                          + wiEdgeRlx
                          + woEdgeRlx
                          + pointDispWall * deltaT;

  Info << "Relaxing inlet..." << endl;
  vectorField inlRelax = inletOutletRlx( mesh_.boundaryMesh()[inletID], vvff);

  Info << "Relaxing outlet..." << endl;
  vectorField outRelax = inletOutletRlx( mesh_.boundaryMesh()[outletID], vvff);

  Info << "Inlet/Outlet relaxation cpuTime: "
       << time.cpuTimeIncrement() << " s" << endl;
  
  
  //mesh_.movePoints( savedPointsAll );
  // *********************************************************************************
  // Final mesh update. 3D
  wiEdgeRlx /= deltaT;
  woEdgeRlx /= deltaT;
  wallRelax /= deltaT;
  pointVelocity.boundaryField()[wallID] == wallRelax + pointDispWall + wiEdgeRlx + woEdgeRlx;
  
  inlRelax /= deltaT;
  pointVelocity.boundaryField()[inletID] == inlRelax + pointDispInlet;

  outRelax /= deltaT;
  pointVelocity.boundaryField()[outletID] == outRelax + pointDispOutlet;

  Info << "Final mesh update" << nl << endl;

  mesh_.update();
  
  // if it is Debug mode finalize the run here
  if( dissolDebug ){
    time.write();
    time++;
    std::exit(0);
  }
}














vectorField meshRelax::normalsOnTheEdge(){
  const List<face>& llf = mesh_.boundaryMesh()[wallID].localFaces();
  const pointField& boundaryPoints = mesh_.boundaryMesh()[wallID].localPoints();
  label NN = local_wall_WallsInletEdges.size();

  const labelListList& plistFaces = mesh_.boundaryMesh()[wallID].pointFaces();

  pointField faceCs = faceCentres(boundaryPoints, llf);
  vectorField faceNs = faceNormals(boundaryPoints, llf);

  vectorField pointNorm( NN, vector::zero );
  scalarList faceToPointSumWeights( NN, 0.0 );
  
  forAll(local_wall_WallsInletEdges, i){
    label&  curI = local_wall_WallsInletEdges[i];
    const point& curP = boundaryPoints[curI];
        
    const labelList& pFaces = plistFaces[curI];
    forAll(pFaces, j){
      label faceI = pFaces[j];
      point& faceC = faceCs[faceI];
      vector d = faceC - curP;
      scalar mag_d = mag(d);
          
      // this is for normal
      scalar nw = 1.0 / mag_d;
      pointNorm[i] += nw * faceNs[ faceI ];
      pointNorm[i].z()=0;
      faceToPointSumWeights[i] += nw;
    }
  }
  syncTools::syncPointList(mesh_, global_WallInletEdges, pointNorm, plusEqOp<vector>(), vector::zero);
  syncTools::syncPointList(mesh_, global_WallInletEdges, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
  forAll(pointNorm, i){
    pointNorm[i] /= mag( pointNorm[i] );
  }

  return pointNorm;
}


Foam::pointField meshRelax::faceCentres(const pointField& points, const List<face>& flist) const
{
  pointField fc( flist.size() );
  forAll(fc, facei)
  {
    fc[facei] = flist[facei].centre(points);
  }
  return fc;
}

Foam::vectorField meshRelax::faceNormals(const pointField& points, const List<face>& flist) const
{
  vectorField fn( flist.size() );
  forAll(fn, facei)
  {
    fn[facei] = flist[facei].normal(points);
    fn[facei] /= mag(fn[facei]) + VSMALL;
  }
  return fn;
}

Foam::vectorField meshRelax::localFaceToPointNormalInterpolate(const pointField& points,
        const pointField& faceCs,
        const vectorField& faceNs,
        const labelListList& pointFaces,
        const labelList& meshPoints
        
) const
{
  int N = points.size();
  scalarList faceToPointSumWeights( N );
  vectorField pointValue( N );
  
  forAll(pointFaces, pointi)
  {
      const labelList& curFaces = pointFaces[pointi];
      pointValue[pointi] = vector::zero;
      faceToPointSumWeights[pointi] = 0.0;
      forAll(curFaces, facei)
      {
        scalar pw = 1.0 / mag(faceCs[curFaces[facei]] - points[pointi]);
        pointValue[pointi] += pw * faceNs[ curFaces[facei] ];
        faceToPointSumWeights[pointi] += pw;
      }
  }

  // synchronization over coupled boundaries
  syncTools::syncPointList(mesh_, meshPoints, pointValue, plusEqOp<vector>(), vector::zero);
  syncTools::syncPointList(mesh_, meshPoints, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);

  // normalization
  forAll(pointValue, pointi){
    pointValue[pointi] /= faceToPointSumWeights[pointi];
  }

  return pointValue;
}


// ++
vectorField meshRelax::calculateInletDisplacement(vectorField& wallDispl){
  // currnt wall points
  const pointField& wallBP = mesh_.boundaryMesh()[wallID].localPoints();
  // list of neighbor faces
  //const labelListList& plistFaces = mesh_.boundaryMesh()[wallID].pointFaces();
  // list of faces
  const List<face>& llf = mesh_.boundaryMesh()[wallID].localFaces();
  // new points coordinates
  pointField curWallBP = wallBP + wallDispl;

  pointField curWallBP1 = wallBP + deltaT*wallDispl;

  // new faceCentres and faceNormals
  vectorField faceNs = faceNormals(curWallBP1, llf);
  coupledPatchInterpolation patchInterpolator( mesh_.boundaryMesh()[wallID], mesh_ );
  vectorField pointNs = patchInterpolator.faceToPointInterpolate(faceNs);
  
  // number of points on the edge
  int N = local_wall_WallsInletEdges.size();
  
  pointField edgePoints(N);
  vectorField edgeNorms(N);
  forAll(local_wall_WallsInletEdges, i){
    label pointI = local_wall_WallsInletEdges[i];
    edgePoints[i] = curWallBP[pointI];
    edgeNorms[i] = pointNs[pointI];
    edgeNorms[i]/=mag(edgeNorms[i]);
  }
  scalar maxZ = gMax( edgePoints.component(vector::Z) );
  // *************************************************************************************

  pointField inP = mesh_.boundaryMesh()[inletID].localPoints();
  
  vectorField pointDispInlet( inP.size(), vector::zero );
  forAll(inP,i){
    pointDispInlet[i].z() = maxZ - inP[i].z();
  }
  
  forAll(local_inlet_WallsInletEdges, i){
    pointDispInlet[ local_inlet_WallsInletEdges[i] ].z() = 0;
  }
  // *************************************************************************************
  vectorField displacement(N, vector::zero);
  
  scalar displ_tol = 1.0;
  int itt = 0;
  while( displ_tol>rlxTol ){
    forAll(local_wall_WallsInletEdges, i){
      displacement[i].z() = (maxZ - edgePoints[i].z());
    }
    vectorField proj_disp = transform(I - edgeNorms*edgeNorms, displacement);
    
    //displ_tol = gAverage( mag(proj_disp) );
    displ_tol = gAverage( mag(proj_disp - displacement) );
    
    scalar displ_tol_norm = gAverage( mag(proj_disp) );
      
    if(displ_tol_norm > SMALL){
      displ_tol /= displ_tol_norm;
    }
    else{
      displ_tol = 0.0;
    }
    
    // it should increase precision, points should end up on the surface
    if(displ_tol>rlxTol){
      edgePoints += proj_disp;
    }
    else{
      edgePoints += displacement;
    }
    
    if(itt%100==0){
      Info << "Wall-inlet  iter " << itt 
           << "  tolerance: "<< displ_tol << endl;
    }

    itt+=1;
  }
  
  forAll(local_wall_WallsInletEdges, i){
    label pointI = local_wall_WallsInletEdges[i];
    wallDispl[ pointI ] = edgePoints[i]-wallBP[pointI];
  }
  return pointDispInlet;
}

vectorField meshRelax::calculateOutletDisplacement(vectorField& wallDispl){
  // currnt wall points
  const pointField& wallBP = mesh_.boundaryMesh()[wallID].localPoints();
  // list of neighbor faces
  // list of faces
  const List<face>& llf = mesh_.boundaryMesh()[wallID].localFaces();
  // new points coordinates
  pointField curWallBP = wallBP+wallDispl;
  
  pointField curWallBP1 = wallBP + deltaT*wallDispl;

  vectorField faceNs = faceNormals(curWallBP1, llf);
  coupledPatchInterpolation patchInterpolator( mesh_.boundaryMesh()[wallID], mesh_ );
  vectorField pointNs = patchInterpolator.faceToPointInterpolate(faceNs);
  
  // number of points on the edge
  int N = local_wall_WallsOutletEdges.size();
  
  pointField edgePoints(N);
  vectorField edgeNorms(N);
  forAll(local_wall_WallsOutletEdges, i){
    label pointI = local_wall_WallsOutletEdges[i];
    edgePoints[i] = curWallBP[pointI];
    edgeNorms[i]  = pointNs[pointI] / mag( pointNs[pointI] );
  }

  scalar minZ = gMin( edgePoints.component(vector::Z) );
  
  // ***************************************************************************
  
  const pointField& outP = mesh_.boundaryMesh()[outletID].localPoints();
  
  vectorField pointDispOutlet( outP.size(), vector::zero );
  forAll(outP,i){
    pointDispOutlet[i].z() = minZ - outP[i].z();
  }
  
  forAll(local_outlet_WallsOutletEdges, i){
    pointDispOutlet[ local_outlet_WallsOutletEdges[i] ].z() = 0;
  }
  
  // ***************************************************************************
  // new faceCentres and faceNormals
  vectorField displacement(N, vector::zero);
  vectorField oldDisplacement(N, vector::zero);
  
  scalar displ_tol = 1.0;
  int itt = 0;
  while( displ_tol>rlxTol ){
    forAll(local_wall_WallsOutletEdges, i){
      displacement[i].z() = minZ - edgePoints[i].z();
    }
    vectorField proj_disp = transform(I - edgeNorms*edgeNorms, displacement);
    
    displ_tol = gAverage( mag(proj_disp - displacement) );
    scalar displ_tol_norm = gAverage( mag(proj_disp) );
      
    if(displ_tol_norm > SMALL){
      displ_tol /= displ_tol_norm;
    }
    else{
      displ_tol = 0.0;
    }
    
    if( displ_tol>rlxTol ){
      edgePoints += proj_disp;
    }
    else{
      edgePoints += displacement;
    }
    
    if(itt%100==0){
      Info << "Wall-outlet iter " << itt 
           << "  tolerance: "<< displ_tol << endl;
    }

    itt+=1;
  }
  
  forAll(local_wall_WallsOutletEdges, i){
    label pointI = local_wall_WallsOutletEdges[i];
    wallDispl[ pointI ] = edgePoints[i]-wallBP[pointI];
  }
  
  return pointDispOutlet;
}


pointField meshRelax::doWallDisplacement(const vectorField& wallDispl){
  pointField newPoints = mesh_.points();
  const labelList& wallsToAll  = mesh_.boundaryMesh()[wallID].meshPoints();
  forAll(wallDispl, i){
    label indx = wallsToAll[i];
    newPoints[indx] += wallDispl[i];
  }
  return newPoints;
}

pointField meshRelax::doInletDisplacement(const vectorField& inletDispl){
  pointField newPoints = mesh_.points();
  const labelList& inletToAll = mesh_.boundaryMesh()[inletID].meshPoints();
  forAll(inletDispl, i){
    label indx = inletToAll[i];
    newPoints[indx] += inletDispl[i];
  }
  return newPoints;
}

pointField meshRelax::doOutletDisplacement(const vectorField& outletDispl){
  pointField newPoints = mesh_.points();
  const labelList& outletToAll = mesh_.boundaryMesh()[outletID].meshPoints();
  forAll(outletDispl, i){
    label indx = outletToAll[i];
    newPoints[indx] += outletDispl[i];
  }
  return newPoints;
}


void meshRelax::fixIWEdgeDispl( vectorField& concNorm ){
  const pointField& boundaryPoints = mesh_.boundaryMesh()[wallID].localPoints();
  
  forAll( inletTriple, i ){
    const labelList& currentTriple = inletTriple[i];
    
    //vector cN0 = concNorm[ currentTriple[0] ];
    //scalar c0 = mag( cN0 );
    vector cN1 = concNorm[ currentTriple[1] ];
    scalar c1 = mag( cN1 );
    vector cN2 = concNorm[ currentTriple[2] ];
    scalar c2 = mag( cN2 );
    
    //scalar newC0 = extrapolateConcentrationLinear(boundaryPoints, c1, c2, currentTriple);
    //scalar newC0 = extrapolateConcentrationExp(boundaryPoints, c1, c2, currentTriple);
    //scalar newC0 = extrapolateConcentrationLinearZ(boundaryPoints, c1, c2, currentTriple);
    scalar newC0 = extrapolateConcentrationExpZ(boundaryPoints, c1, c2, currentTriple);
    vector newN0 = extrapolateVectorLinear(boundaryPoints, cN1, cN2, currentTriple);
    newN0 /= mag(newN0);
    
    concNorm[ currentTriple[0] ] = newC0 * newN0;
    //concNorm[ currentTriple[0] ] = newC0 * cN0 / c0;
  }
}



scalar meshRelax::extrapolateConcentrationExpZ(const pointField& loc_points,
                                scalar& c1, scalar& c2, 
                                const labelList& pnts){
  scalar r02 = loc_points[ pnts[0] ].z() - loc_points[ pnts[2] ].z();
  scalar r12 = loc_points[ pnts[1] ].z() - loc_points[ pnts[2] ].z();

  return c2 * std::pow( c1/c2, r02/r12 );
}
scalar meshRelax::extrapolateConcentrationExp(const pointField& loc_points,
                                scalar& c1, scalar& c2, 
                                const labelList& pnts){
  scalar r02 = mag( loc_points[ pnts[0] ] - loc_points[ pnts[2] ] );
  scalar r12 = mag( loc_points[ pnts[1] ] - loc_points[ pnts[2] ] );

  return c2 * std::pow( c1/c2, r02/r12 );
}
scalar meshRelax::extrapolateConcentrationLinear(const pointField& loc_points,
                                scalar& c1, scalar& c2, 
                                const labelList& pnt){
  scalar r02 = mag(  loc_points[pnt[0]] - loc_points[pnt[2]]  );
  scalar r12 = mag(  loc_points[pnt[1]] - loc_points[pnt[2]]  );

  return c2 - (c2-c1) * r02/r12;
}
vector meshRelax::extrapolateVectorLinear(const pointField& loc_points,
                                vector& r1, vector& r2, 
                                const labelList& pnt){
  scalar r02 = mag(  loc_points[pnt[0]] - loc_points[pnt[2]]  );
  scalar r12 = mag(  loc_points[pnt[1]] - loc_points[pnt[2]]  );

  return r2 - (r2-r1) * r02/r12;
}
scalar meshRelax::extrapolateConcentrationLinearZ(const pointField& loc_points,
                                scalar& c1, scalar& c2, 
                                const labelList& pnts){
  scalar r02 = loc_points[ pnts[0] ].z() - loc_points[ pnts[2] ].z();
  scalar r12 = loc_points[ pnts[1] ].z() - loc_points[ pnts[2] ].z();
  return c2 - (c2-c1) * r02/r12;
}


void meshRelax::setUpPairsConc()
{
  const labelListList& plistEdges = mesh_.boundaryMesh()[wallID].pointEdges();
  const edgeList& edgeL = mesh_.boundaryMesh()[wallID].edges();
  
  inletTriple.setSize(local_wall_WallsInletEdges.size());
  
  labelList secondLine(local_wall_WallsInletEdges.size());
  
  forAll(local_wall_WallsInletEdges, i)
  {
    inletTriple[i].setSize( 3 );
    inletTriple[i][0] = local_wall_WallsInletEdges[i];
    
    const labelList& currentNextEdges = plistEdges[local_wall_WallsInletEdges[i]];
    
    label secondLineP = -1;
    
    forAll(currentNextEdges, j)
    {
      edge edg = edgeL[ currentNextEdges[j] ];
      
      label edgePair = -1;
      
      if( edg.start() == local_wall_WallsInletEdges[i] ){
        edgePair = edg.end();
      }
      else{
        edgePair = edg.start();
      }
      
      if( findIndex(local_wall_WallsInletEdges, edgePair) == -1){
        secondLineP = edgePair;
        break;
      }
    }
    inletTriple[i][1] = secondLineP;
    secondLine[i] = secondLineP;
  }
  
  
  forAll(secondLine, i)
  {
    const labelList& currentNextEdges = plistEdges[ secondLine[i] ];
    
    label thirdLineP = -1;
    
    forAll(currentNextEdges, j)
    {
      edge edg = edgeL[ currentNextEdges[j] ];
      
      label edgePair = -1;
      
      if( edg.start() == secondLine[i] ){
        edgePair = edg.end();
      }
      else{
        edgePair = edg.start();
      }
      
      if( findIndex(secondLine, edgePair) == -1 && findIndex(local_wall_WallsInletEdges, edgePair) == -1){
        thirdLineP = edgePair;
        break;
      }
    }
    inletTriple[i][2] = thirdLineP;
  }
  
}

void meshRelax::commonPoints
(
  const labelList& list1, 
  const labelList& list2, 
  labelList& localList,
  labelList& globalList
)
{
  forAll(list1, loclLabel1){
    label globLabel1 = list1[loclLabel1];
    label loclLabel2 = findIndex(list2, globLabel1);
    if( loclLabel2 != -1){
      localList.append( loclLabel1 );
      globalList.append( globLabel1 );
    }
  }
}

void meshRelax::neighborListEdge
(
  const labelList& list,
  const edgeList& eList,
  const labelListList& pointEdgesLi,
  labelListList& neLi
)
{
  forAll(list, i){
    label ind = list[i];

    const labelList& lll = pointEdgesLi[ind];
    labelList nel;
    forAll(lll, ii){
      const edge& ed = eList[ lll[ii] ];
      label edb = ed.start();

      if( edb!=ind && findIndex(list, edb) != -1 ){
        nel.append(edb);
      }

      label ede = ed.end();
      if( ede!=ind && findIndex(list, ede) != -1){
        nel.append(ede);
      }
    }
    neLi.append(nel);
  }
}



void meshRelax::setUpLists()
{
  const labelList& wallsToAll  = mesh_.boundaryMesh()[wallID].meshPoints();
  const labelList& inletToAll  = mesh_.boundaryMesh()[inletID].meshPoints();
  const labelList& outletToAll = mesh_.boundaryMesh()[outletID].meshPoints();
  
  //commonPoints(wallsToAll, inletToAll, local_inlet_WallsInletEdges, global_WallEdges);
  
  
  forAll(wallsToAll, i){
    label lW = wallsToAll[i];
    
    label inlind = findIndex(inletToAll, lW);
    if( inlind != -1){
      local_wall_WallsInletEdges.append( i );
      global_WallInletEdges.append( lW );
      local_inlet_WallsInletEdges.append( inlind );
    }
    
    // create walls-outlet edge list of vertex IDs
    if( findIndex(outletToAll, lW) != -1){
      local_wall_WallsOutletEdges.append( i );
      global_WallOutletEdges.append( lW );
    }
  }
  
  // if cyclic boundary exists
  const label& cycID1 = mesh_.boundaryMesh().findPatchID("periodicx1");
  const label& cycID2 = mesh_.boundaryMesh().findPatchID("periodicx2");

  if(cycID1!=-1 && cycID2!=-1){
    // map vetex ID: patch to global
    const labelList& cycToAll1  = mesh_.boundaryMesh()[cycID1].meshPoints();
    const labelList& cycToAll2  = mesh_.boundaryMesh()[cycID2].meshPoints();
    
    forAll(wallsToAll, i)
    {
      label lW = wallsToAll[i];

      if( findIndex(cycToAll1, lW) != -1){
        local_wall_WallsCycEdges1.append( i );
      }
      if( findIndex(cycToAll2, lW) != -1){
        local_wall_WallsCycEdges2.append( i );
      }
    }
  }
  
}

word meshRelax::get_version() const{
  return version;
}
word meshRelax::get_date() const{
  return date;
}
