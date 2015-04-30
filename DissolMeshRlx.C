/* 
 * Surface relaxation for dissolFoam
 * Dec 09 2014
 * 
 */

#include "DissolMeshRlx.H"
#include <algorithm>

#include "Pstream.H"

#include "pyramidPointFaceRef.H"

#include "transformField.H"
#include "symmTransformField.H"

#include "triPointRef.H"


// mesh1 is the mesh at time 0
DissolMeshRlx::DissolMeshRlx( const fvMesh& mesh)
:
  version(0.4),
  mesh_(mesh)
{
  // get ID of each patch we need
  wallID   = mesh_.boundaryMesh().findPatchID("walls");
  inletID  = mesh_.boundaryMesh().findPatchID("inlet");
  outletID = mesh_.boundaryMesh().findPatchID("outlet");

  // map vetex ID: patch to global
  wallsToAll  = mesh_.boundaryMesh()[wallID].meshPoints();
  inletToAll  = mesh_.boundaryMesh()[inletID].meshPoints();
  outletToAll = mesh_.boundaryMesh()[outletID].meshPoints();

  setUpLists();
  setUpPairsConc();
  
  //scalar pn = 3*(1+Pstream::myProcNo());
  //Pout<< "proc: " << Pstream::myProcNo() << "   master?: "<< Pstream::master() << nl;
  //Pout<< "wallsToAll: " << wallsToAll.size() << "\n";
  //reduce(pn, maxOp<scalar>());
  //Pout<< "After pn: " << pn  << nl;
  //std::exit(0);
}

vectorField DissolMeshRlx::normalsOnTheEdge(){

  const List<face>& llf = mesh_.boundaryMesh()[wallID].localFaces();
  pointField boundaryPoints = mesh_.boundaryMesh()[wallID].localPoints();
  label NN = local_wall_WallsInletEdges.size();

  const labelListList& plistFaces = mesh_.boundaryMesh()[wallID].pointFaces();

  pointField faceCs = faceCentres(boundaryPoints, llf);
  vectorField faceNs = faceNormals(boundaryPoints, llf);

  vectorField pointNorm( NN, vector::zero );
  scalarList faceToPointSumWeights( NN, 0.0 );
  

  forAll(local_wall_WallsInletEdges, i){
    label  curI = local_wall_WallsInletEdges[i];
    point& curP = boundaryPoints[curI];
        
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

scalarField DissolMeshRlx::distanceToTheEdge(vectorField& nn){
  pointField ppoints = mesh_.boundaryMesh()[inletID].localPoints();
  pointField pfaces = mesh_.boundaryMesh()[inletID].faceCentres();
  
  scalarField distance( pfaces.size(), 0.0 );
  
  forAll(pfaces, fi){
    point& curF = pfaces[fi];
    
    scalarField alld( local_inlet_WallsInletEdges.size(), 0.0 );
    vectorField allvec( local_inlet_WallsInletEdges.size(), vector::zero );
    forAll(local_inlet_WallsInletEdges, pi){
      point& curP = ppoints[local_inlet_WallsInletEdges[pi]];
      vector dd = curF-curP;

      //alld[pi] = mag( curF-curP );
      //alld[pi] = mag( dd ^ nn[pi]);
      //alld[pi] = mag( dd );
      alld[pi] = mag(dd);
      allvec[pi] = dd;

      //alld[pi] = std::abs( mag( dd & nn[pi]) );
    }
    //distance[fi] = min( alld );


    label minid = findMin( alld );  
    distance[fi] = mag( allvec[minid] & nn[minid] );


  }
  
  return distance;
}


// ++
Foam::pointField DissolMeshRlx::faceCentres(const pointField& points, const List<face>& flist) const
{
  pointField fc( flist.size() );
  forAll(fc, facei)
  {
    fc[facei] = flist[facei].centre(points);
  }
  return fc;
}

Foam::vectorField DissolMeshRlx::faceNormals(const pointField& points, const List<face>& flist) const
{
  vectorField fn( flist.size() );
  forAll(fn, facei)
  {
    fn[facei] = flist[facei].normal(points);
    fn[facei] /= mag(fn[facei]) + VSMALL;
  }
  return fn;
}

Foam::vectorField DissolMeshRlx::localFaceToPointNormalInterpolate(const pointField& points,
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
vectorField DissolMeshRlx::calculateInletDisplacement(vectorField& wallDispl){

  scalar maxdZ = 0;
  forAll(local_wall_WallsInletEdges, i){
    vector A = wallDispl[ local_wall_WallsInletEdges[i] ];
    maxdZ = max( maxdZ, A.z() );    
  }

  reduce(maxdZ, maxOp<scalar>());
  
  vectorField pointDispInlet( mesh_.boundaryMesh()[inletID].localPoints().size(), vector::zero );
  pointDispInlet.replace( vector::Z, maxdZ);
  forAll(local_inlet_WallsInletEdges, i){
    pointDispInlet[ local_inlet_WallsInletEdges[i] ].z() = 0;
  }
  
  forAll(local_wall_WallsInletEdges, i){
    vector A = wallDispl[ local_wall_WallsInletEdges[i] ];
    vector dz(0.0, 0.0, maxdZ - A.z());
          
    wallDispl[ local_wall_WallsInletEdges[i] ] += dz;
  }

  return pointDispInlet;
}

vectorField DissolMeshRlx::calculateOutletDisplacement(vectorField& wallDispl){

  scalar maxdZ = 0;
  forAll(local_wall_WallsOutletEdges, i){
    vector A = wallDispl[ local_wall_WallsOutletEdges[i] ];
    if( mag(A.z()) > maxdZ ){
      maxdZ = A.z();
    }
  }

  reduce(maxdZ, maxOp<scalar>());
  
  vectorField pointDispOutlet( mesh_.boundaryMesh()[outletID].localPoints().size(), vector::zero );
  pointDispOutlet.replace( vector::Z, maxdZ);
  forAll(local_outlet_WallsOutletEdges, i){
    pointDispOutlet[ local_outlet_WallsOutletEdges[i] ].z() = 0;
  }

  forAll(local_wall_WallsOutletEdges, i){
    vector A = wallDispl[ local_wall_WallsOutletEdges[i] ];
    vector dz(0.0, 0.0, maxdZ - A.z());
          
    wallDispl[ local_wall_WallsOutletEdges[i] ] += dz;
  }

  return pointDispOutlet;
}


pointField DissolMeshRlx::doInletDisplacement(const vectorField& inletDispl){
  pointField newPoints = mesh_.points();
  forAll(inletDispl, i){
    label indx = inletToAll[i];
    newPoints[indx] += inletDispl[i];
  }
  return newPoints;
}

pointField DissolMeshRlx::doWallDisplacement(const vectorField& wallDispl){
  pointField newPoints = mesh_.points();
  forAll(wallDispl, i){
    label indx = wallsToAll[i];
    newPoints[indx] += wallDispl[i];
  }
  return newPoints;
}

pointField DissolMeshRlx::doOutletDisplacement(const vectorField& outletDispl){
  pointField newPoints = mesh_.points();
  forAll(outletDispl, i){
    label indx = outletToAll[i];
    newPoints[indx] += outletDispl[i];
  }
  return newPoints;
}


// ++
void DissolMeshRlx::fixEdgeConcentration( vectorField& concNorm ){
  const pointField& boundaryPoints = mesh_.boundaryMesh()[wallID].localPoints();
  
  forAll( inletTriple, i ){
    const labelList& currentTriple = inletTriple[i];
    
    vector cN0 = concNorm[ currentTriple[0] ];
    scalar c0 = mag( cN0 );
    vector cN1 = concNorm[ currentTriple[1] ];
    scalar c1 = mag( cN1 );
    vector cN2 = concNorm[ currentTriple[2] ];
    scalar c2 = mag( cN2 );
    
    //scalar newC0 = extrapolateConcentrationLinear(boundaryPoints, c1, c2, currentTriple);
    //scalar newC0 = extrapolateConcentrationExp(boundaryPoints, c1, c2, currentTriple);
    scalar newC0 = extrapolateConcentrationLinearZ(boundaryPoints, c1, c2, currentTriple);
    //scalar newC0 = extrapolateConcentrationExpZ(boundaryPoints, c1, c2, currentTriple);
    vector newN0 = extrapolateVectorLinear(boundaryPoints, cN1, cN2, currentTriple);
    newN0 /= mag(newN0);
    
    //concNorm[ currentTriple[0] ] = newC0 * newN0;
    concNorm[ currentTriple[0] ] = newC0 * cN0 / c0;
    
    //Pout << (cN0 / c0) << "   c0: "<< c0 << "     new: "<< newC0 << nl;
  }
}



scalar DissolMeshRlx::extrapolateConcentrationExpZ(const pointField& loc_points,
                                scalar& c1, scalar& c2, 
                                const labelList& pnts){
  scalar r02 = loc_points[ pnts[0] ].z() - loc_points[ pnts[2] ].z();
  scalar r12 = loc_points[ pnts[1] ].z() - loc_points[ pnts[2] ].z();

  return c2 * std::pow( c1/c2, r02/r12 );
}

scalar DissolMeshRlx::extrapolateConcentrationExp(const pointField& loc_points,
                                scalar& c1, scalar& c2, 
                                const labelList& pnts){
  scalar r02 = mag( loc_points[ pnts[0] ] - loc_points[ pnts[2] ] );
  scalar r12 = mag( loc_points[ pnts[1] ] - loc_points[ pnts[2] ] );

  return c2 * std::pow( c1/c2, r02/r12 );
}


scalar DissolMeshRlx::extrapolateConcentrationLinear(const pointField& loc_points,
                                scalar& c1, scalar& c2, 
                                const labelList& pnt){
  scalar r02 = mag(  loc_points[pnt[0]] - loc_points[pnt[2]]  );
  scalar r12 = mag(  loc_points[pnt[1]] - loc_points[pnt[2]]  );

  return c2 - (c2-c1) * r02/r12;
}

vector DissolMeshRlx::extrapolateVectorLinear(const pointField& loc_points,
                                vector& r1, vector& r2, 
                                const labelList& pnt){
  scalar r02 = mag(  loc_points[pnt[0]] - loc_points[pnt[2]]  );
  scalar r12 = mag(  loc_points[pnt[1]] - loc_points[pnt[2]]  );

  return r2 - (r2-r1) * r02/r12;
}


scalar DissolMeshRlx::extrapolateConcentrationLinearZ(const pointField& loc_points,
                                scalar& c1, scalar& c2, 
                                const labelList& pnts){
  scalar r02 = loc_points[ pnts[0] ].z() - loc_points[ pnts[2] ].z();
  scalar r12 = loc_points[ pnts[1] ].z() - loc_points[ pnts[2] ].z();
  return c2 - (c2-c1) * r02/r12;
}


void DissolMeshRlx::setUpPairsConc(){
  
  const labelListList& plistEdges = mesh_.boundaryMesh()[wallID].pointEdges();
  const edgeList& edgeL = mesh_.boundaryMesh()[wallID].edges();
  
  
  inletTriple.setSize(local_wall_WallsInletEdges.size());
  
  
  labelList secondLine(local_wall_WallsInletEdges.size());
  
  forAll(local_wall_WallsInletEdges, i){
    inletTriple[i].setSize( 3 );
    inletTriple[i][0] = local_wall_WallsInletEdges[i];
    
    const labelList& currentNextEdges = plistEdges[local_wall_WallsInletEdges[i]];
    
    label secondLineP = -1;
    
    forAll(currentNextEdges, j){
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
  
  
  forAll(secondLine, i){
    
    const labelList& currentNextEdges = plistEdges[ secondLine[i] ];
    
    label thirdLineP = -1;
    
    forAll(currentNextEdges, j){
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

void DissolMeshRlx::setUpLists(){

  forAll(wallsToAll, i){
    label lW = wallsToAll[i];
    
    label inlind = findIndex(inletToAll, lW);
    if( inlind != -1){
      local_wall_WallsInletEdges.append( i );
      global_WallInletEdges.append( lW );
      local_inlet_WallsInletEdges.append( inlind );
    }
    
    // create walls-outlet edge list of vertex IDs
    //std::map<int,int>::iterator iO = allToOutlet.find( lW );
    //if(iO != allToOutlet.end()){
    if( findIndex(outletToAll, lW) != -1){
      local_wall_WallsOutletEdges.append( i );
      global_WallOutletEdges.append( lW );
    }
  }
  
  scaleList.setSize( mesh_.boundaryMesh()[inletID].nPoints(), -1 );
  const pointField& pCoord = mesh_.boundaryMesh()[inletID].localPoints();
  
  forAll(scaleList, i){
    point pp = pCoord[ i ];
    forAll(local_inlet_WallsInletEdges, j){
      if( std::abs(pp.x()-pCoord[local_inlet_WallsInletEdges[j]].x())<0.1 
              &&  
          pp.y() * pCoord[local_inlet_WallsInletEdges[j]].y()>0
        ){
        scaleList[i] = j;
        break;
      }
    }
  }
  
  scaleListOutlet.setSize( mesh_.boundaryMesh()[outletID].nPoints(), -1 );
  const pointField& pCoordOut = mesh_.boundaryMesh()[outletID].localPoints();
  
  forAll(scaleListOutlet, i){
    point pp = pCoordOut[ i ];
    forAll(local_outlet_WallsOutletEdges, j){
      if( std::abs(pp.x()-pCoordOut[local_outlet_WallsOutletEdges[j]].x())<0.1 
              &&  
          pp.y() * pCoord[local_outlet_WallsOutletEdges[j]].y()>0
        ){
        scaleListOutlet[i] = j;
        break;
      }
    }
  }
}

float DissolMeshRlx::get_version(){
  return version;
}

