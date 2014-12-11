/* 
 * Surface relaxation for dissolFoam
 * Dec 09 2014
 * 
 */

#include "DissolMeshRlx.H"
#include <algorithm>

#include "Pstream.H"

#include "pyramidPointFaceRef.H"
#include "syncTools.H"

#include "transformField.H"
#include "symmTransformField.H"
#include "primitivePatchInterpolationSync.H"


// mesh1 is the mesh at time 0
DissolMeshRlx::DissolMeshRlx( const fvMesh& mesh)
:
  version(0.3),
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

// ++
vectorField DissolMeshRlx::wallRelaxation(){
  const labelList& meshPoints = mesh_.boundaryMesh()[wallID].meshPoints();

  const pointField& boundaryPoints = mesh_.boundaryMesh()[wallID].localPoints();
  const labelListList& plistFaces = mesh_.boundaryMesh()[wallID].pointFaces();
  const pointField& faceCs = mesh_.boundaryMesh()[wallID].faceCentres();
  
  vectorField displacement( boundaryPoints.size() );
  scalarField numberFaces( boundaryPoints.size() );
  
  forAll(boundaryPoints, i){
    point curP = boundaryPoints[i];
    const labelList& pFaces = plistFaces[i];
    
    numberFaces[i] = plistFaces[i].size();
    
    vector sumw(0,0,0);
    forAll(pFaces, j){
      label faceI = pFaces[j];
      point faceC = faceCs[faceI];
      
      vector d = faceC - curP;
      sumw += d;
    }
    // sum of all distances to face centres
    displacement[i] = sumw;
  }
  
  syncTools::syncPointList( mesh_, meshPoints, displacement, plusEqOp<vector>(), vector::zero);
  syncTools::syncPointList( mesh_, meshPoints, numberFaces, plusEqOp<scalar>(), 0.0);
  
  /*
  forAll(numberFaces, i){
    if(numberFaces[i]>2 && numberFaces[i]<4){
      Pout << numberFaces[i] << nl;
    }
  }
  std::exit(0);
   */

  forAll(displacement, i){
    displacement[i] /= numberFaces[i];
  }
  
  // getting correct point normal field
  vectorField fNorm = mesh_.boundaryMesh()[wallID].faceNormals();
  primitivePatchInterpolationSync patchInterpolator( mesh_.boundaryMesh()[wallID], mesh_ );
  vectorField pointNorm = patchInterpolator.faceToPointInterpolate(fNorm);
  
  vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);


  // project back inlet
  // TODO read from the dictionary
  vector ny(0,1,0);
  vector nz(0,0,1);
  
  vectorField inletEdgeD( local_wall_WallsInletEdges.size() );
  forAll(local_wall_WallsInletEdges, i){
    inletEdgeD[i] = projectedDisplacement[ local_wall_WallsInletEdges[i] ];
  }
  
  vectorField inletPrjY = transform(I - ny*ny, inletEdgeD);
  vectorField inletPrj = transform(I - nz*nz, inletPrjY);
  
  forAll(local_wall_WallsInletEdges, i){
    projectedDisplacement[ local_wall_WallsInletEdges[i] ] = inletPrj[i];
  }
  
  // project back outlet
  vectorField outletEdgeD( local_wall_WallsOutletEdges.size() );
  forAll(local_wall_WallsOutletEdges, i){
    outletEdgeD[i] = projectedDisplacement[ local_wall_WallsOutletEdges[i] ];
  }
  
  vectorField outletPrjY = transform(I - ny*ny, outletEdgeD);
  vectorField outletPrj = transform(I - nz*nz, outletPrjY);

  forAll(local_wall_WallsOutletEdges, i){
    projectedDisplacement[ local_wall_WallsOutletEdges[i] ] = outletPrj[i];
  }
  
  return projectedDisplacement;
}



// ++
vectorField DissolMeshRlx::calculateInletDisplacement(vectorField& wallDispl){

  scalar maxdZ = 0;
  forAll(local_wall_WallsInletEdges, i){
    vector A = wallDispl[ local_wall_WallsInletEdges[i] ];
    maxdZ = max( maxdZ, A.z() );    
  }

  reduce(maxdZ, maxOp<scalar>());

  forAll(local_wall_WallsInletEdges, i){
    vector A = wallDispl[ local_wall_WallsInletEdges[i] ];
    vector dz(0.0, 0.0, maxdZ - A.z());
          
    wallDispl[ local_wall_WallsInletEdges[i] ] += dz;
  }

  vectorField pointDispInlet( mesh_.boundaryMesh()[inletID].localPoints().size(), vector::zero );
      
  pointDispInlet.replace( vector::Z, maxdZ);

  return pointDispInlet;
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
    
    scalar newC0 = extrapolateConcentrationLinearZ(boundaryPoints, c1, c2, currentTriple);
    //scalar newC0 = extrapolateConcentrationExp(boundaryPoints, c1, c2, currentTriple) + 4.0;
    
    concNorm[ currentTriple[0] ] = newC0 * cN0 / c0;
    
    //Pout << (cN0 / c0) << "   c0: "<< c0 << "     new: "<< newC0 << nl;
  }
}



scalar DissolMeshRlx::extrapolateConcentrationExp(const pointField& loc_points,
                                scalar& c1, scalar& c2, 
                                const labelList& pnts){
  scalar r02 = loc_points[ pnts[0] ].z() - loc_points[ pnts[2] ].z();
  scalar r12 = loc_points[ pnts[1] ].z() - loc_points[ pnts[2] ].z();

  return c2 * std::pow( c1/c2, r02/r12 );
}

scalar DissolMeshRlx::extrapolateConcentrationLinear(const pointField& loc_points,
                                const scalarField& point_conc,
                                label pnt0, label pnt1, label pnt2){
  scalar r02 = mag(  loc_points[pnt0] - loc_points[pnt2]  );
  scalar r12 = mag(  loc_points[pnt1] - loc_points[pnt2]  );
  scalar c1 = point_conc[pnt1];
  scalar c2 = point_conc[pnt2];

  return c2 - (c2-c1) * r02/r12;
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
    
    if( findIndex(inletToAll, lW) != -1){
      local_wall_WallsInletEdges.append( i );
    }
    
    // create walls-outlet edge list of vertex IDs
    //std::map<int,int>::iterator iO = allToOutlet.find( lW );
    //if(iO != allToOutlet.end()){
    if( findIndex(outletToAll, lW) != -1){
      local_wall_WallsOutletEdges.append( i );
    }
  }
  
}

float DissolMeshRlx::get_version(){
  return version;
}

