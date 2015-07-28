/* 
 * Surface relaxation for dissolFoam
 * Dec 09 2014
 * 
 */

#include "DissolMeshRlx.H"
#include <algorithm>

#include "Pstream.H"

#include "pyramidPointFaceRef.H"

//#include "transformField.H"
//#include "symmTransformField.H"

#include "triPointRef.H"


// mesh1 is the mesh at time 0
DissolMeshRlx::DissolMeshRlx( const fvMesh& mesh)
:
  version(0.5),
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
  
  /*
    scalar pn = 3*(1+Pstream::myProcNo());
    Pout<< "proc: " << Pstream::myProcNo() << "   master?: "<< Pstream::master() << nl;
    Pout<< "wallsToAll: " << wallsToAll.size() << "\n";
    reduce(pn, maxOp<scalar>());
    Pout<< "After pn: " << pn  << nl;
    std::exit(0);
  */
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

      alld[pi] = mag(dd);
      allvec[pi] = dd;
    }
    label minid = findMin( alld );  
    distance[fi] = mag( allvec[minid] & nn[minid] );
  }
  return distance;
}

scalarField DissolMeshRlx::newC(scalarField& phph, scalar rlxTol){
  
  pointField boundaryPoints = mesh_.boundaryMesh()[inletID].localPoints();
  int N = boundaryPoints.size();
  
  const labelList& meshPoints = mesh_.boundaryMesh()[inletID].meshPoints();

  const labelListList& ppp = mesh_.boundaryMesh()[inletID].pointEdges();
  const edgeList& ee = mesh_.boundaryMesh()[inletID].edges();
  
  scalarList point_weight(N, 1.0);
  syncTools::syncPointList(mesh_, meshPoints, point_weight, plusEqOp<scalar>(), 0.0);
  forAll(point_weight, i){
    point_weight[i] = 1/point_weight[i];
  }

  const labelList& wallsTo  = mesh_.boundaryMesh()[wallID].meshPoints();
  labelList local_inlet_WallEdges;
  labelList global_InletEdges;
  
  forAll(meshPoints, i){
    label lW = meshPoints[i];

    label inlind = findIndex(wallsTo, lW);
    if( inlind != -1){
      local_inlet_WallEdges.append( i );
      global_InletEdges.append( lW );
    }
    
  }

  labelList local_neighbours(local_inlet_WallEdges.size(), -1);
  forAll(local_inlet_WallEdges, i){
    
    const labelList& currentNextEdges = ppp[local_inlet_WallEdges[i]];
    
    label secondLineP = -1;
    
    forAll(currentNextEdges, j){
      edge edg = ee[ currentNextEdges[j] ];
      
      label edgePair = -1;
      
      if( edg.start() == local_inlet_WallEdges[i] ){
        edgePair = edg.end();
      }
      else{
        edgePair = edg.start();
      }
      
      if( findIndex(local_inlet_WallEdges, edgePair) == -1){
        secondLineP = edgePair;
        break;
      }
    }
    local_neighbours[i] = secondLineP;
  }
  
  double displ_tol = 1.0;
  int itt = 0;
  scalarField nC(phph.size(), 1.0);
  scalar lp = 0.5;

  forAll(local_inlet_WallEdges, i){
    label ind_inl = local_inlet_WallEdges[i];
    label nei_inl = local_neighbours[i];

    vector del = boundaryPoints[ind_inl] - boundaryPoints[nei_inl];
    scalar kf = lp / (lp + mag(del));
    
    nC[ ind_inl ] = kf * nC[nei_inl];
  }
  
  while(displ_tol>rlxTol){
    
    scalarList sumWeights(N, 0.0);
    scalarField newCC(N, 0.0);
    
    forAll(newCC, i){
      point& curP = boundaryPoints[i];

      const labelList& lll = ppp[i];
      labelList nl(lll.size(), -1);
      forAll(lll, ii){
        const edge& ed = ee[lll[ii]];
        label edb = ed.start();
        label ede = ed.end();
        if(i != edb){
          nl[ii] = edb;
        }
        else{
          nl[ii] = ede;
        }
      }

      forAll(nl, j){
        label pI = nl[j];
        scalar cc = nC[pI];
        point& np = boundaryPoints[pI];

        vector wv = np - curP;
        scalar mag_wv = mag(wv);
        scalar nw = mag_wv;//1.0 / mag_wv;
        
        nw *= nw;
        nw *= nw;
        //nw *= nw;

        scalar weight = 1.0;
        
        if(point_weight[i]<1.0){
          weight = nw * point_weight[pI];
        }
        else{
          weight = nw;
        }
        scalar conc = weight * cc;
        newCC[i] += conc;
        sumWeights[i] += weight;
      }
    }

    syncTools::syncPointList(mesh_, meshPoints, newCC, plusEqOp<scalar>(), 0.0);
    syncTools::syncPointList(mesh_, meshPoints, sumWeights, plusEqOp<scalar>(), 0.0);
    // normalization
    forAll(newCC, pointi){
      newCC[pointi] /= sumWeights[pointi];
    }

    
    forAll(local_inlet_WallEdges, i){
      label ind_inl = local_inlet_WallEdges[i];
      label nei_inl = local_neighbours[i];

      vector del = boundaryPoints[ind_inl] - boundaryPoints[nei_inl];
      scalar kf = lp / (lp + mag(del));

      //newCC[ ind_inl ] = kf * kf * newCC[nei_inl];
      newCC[ ind_inl ] = 1.0102 * kf * newCC[nei_inl];
    }
    
    
    scalar totPhi = 0.0;
    scalar totCPhi = 0.0;
    forAll(phph, i){
      totPhi += phph[i];
      totCPhi += newCC[i] * phph[i];
    }

    scalar corrF = 0.0;
    if(totCPhi!=0.0){
      corrF = totPhi/totCPhi;
    }

    newCC *= corrF;
    

    scalarField aux_f = mag(newCC - nC);
    scalarField aux_f0(nC.size(), 0.0);

    if( aux_f == aux_f0 ){
      displ_tol = 0.0;
    }
    else{
      displ_tol = average( aux_f );
    }

    scalarField aux_fn = newCC;
    scalarField aux_fn0(nC.size(), 0.0);
    scalar displ_tol_norm;
    if( aux_fn == aux_fn0 ){
      displ_tol_norm = 0.0;
    }
    else{
      displ_tol_norm = average( aux_fn );
    }

    reduce(displ_tol, sumOp<scalar>());
    reduce(displ_tol_norm, sumOp<scalar>());
    displ_tol /= displ_tol_norm;

    nC = newCC;

    if(itt%100==0){
      Info<<" Inlet C: Iter "<<itt<<"  tolerance: "<<displ_tol<< nl;
    }
    itt+=1;
  }

  Info<< " Inlet C: converged in "<<itt<<" iterations. Tolerance: " << displ_tol<< nl;
  
  return nC;
}

scalarField DissolMeshRlx::newC1(scalarField& phph, scalar rlxTol){
  return phph;
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
  // currnt wall points
  const pointField& wallBP = mesh_.boundaryMesh()[wallID].localPoints();
  // list of neighbor faces
  const labelListList& plistFaces = mesh_.boundaryMesh()[wallID].pointFaces();
  // list of faces
  const List<face>& llf = mesh_.boundaryMesh()[wallID].localFaces();
  // new points coordinates
  pointField curWallBP = wallBP + wallDispl;

  // new faceCentres and faceNormals
  //pointField faceCs = faceCentres(curWallBP, llf);
  vectorField faceNs = faceNormals(curWallBP, llf);
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
  scalar maxZ = max( edgePoints.component(vector::Z) );
  reduce(maxZ, maxOp<scalar>());
  
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
  
  //const pointField& boundaryPoints = mesh_.boundaryMesh()[inletID].localPoints();
  scalar tol = 1.0;
  scalar tolerance = 0.00000000001;
  int itt = 0;
  while( tol>tolerance ){
    forAll(local_wall_WallsInletEdges, i){
      displacement[i].z() = (maxZ - edgePoints[i].z());
    }
    vectorField proj_disp = transform(I - edgeNorms*edgeNorms/8.0, displacement);
    //vectorField proj_disp = (displacement + transform(I - edgeNorms*edgeNorms*2.0, displacement))/2.0;
    
    /*
    if( Pstream::master() ){
      vector aa = displacement[0] - proj_disp[0];
      Info<<displacement[0]<<"  "
              <<proj_disp[0]<<"  "
              <<edgeNorms[0]<<"  aa "
              <<aa<<"  "
              <<nl;
    
    }
    */
    
    scalarField aux_f = mag(proj_disp);
    scalarField aux_f0(N, 0.0);

    if( aux_f == aux_f0 ){
      tol = 0.0;
    }
    else{
      tol = average( aux_f );
    }
    reduce(tol, sumOp<scalar>());
    tol /= static_cast<double>( Pstream::nProcs() );
    
    if(tol>tolerance){
      edgePoints += proj_disp;
    }
    else{
      edgePoints += displacement;
    }
    
    if(itt%100==0){
      Info<<"  Wall-inlet iter "<<itt<<"  tolerance: "<<tol<< nl;
    }

    itt+=1;
  }
  
  forAll(local_wall_WallsInletEdges, i){
    label pointI = local_wall_WallsInletEdges[i];
    wallDispl[ pointI ] = edgePoints[i]-wallBP[pointI];
  }
  
  /*
  forAll(local_wall_WallsInletEdges, i){
    vector A = wallDispl[ local_wall_WallsInletEdges[i] ];
    vector dz(0.0, 0.0, maxdZ - A.z());
          
    wallDispl[ local_wall_WallsInletEdges[i] ] += dz;
  }
  */

  return pointDispInlet;
}

vectorField DissolMeshRlx::calculateOutletDisplacement(vectorField& wallDispl){
  // currnt wall points
  const pointField& wallBP = mesh_.boundaryMesh()[wallID].localPoints();
  // list of neighbor faces
  const labelListList& plistFaces = mesh_.boundaryMesh()[wallID].pointFaces();
  // list of faces
  const List<face>& llf = mesh_.boundaryMesh()[wallID].localFaces();
  // new points coordinates
  pointField curWallBP = wallBP+wallDispl;
  

  //pointField faceCs = faceCentres(curWallBP, llf);
  vectorField faceNs = faceNormals(curWallBP, llf);
  coupledPatchInterpolation patchInterpolator( mesh_.boundaryMesh()[wallID], mesh_ );
  vectorField pointNs = patchInterpolator.faceToPointInterpolate(faceNs);
  
  // number of points on the edge
  int N = local_wall_WallsOutletEdges.size();
  
  pointField edgePoints(N);
  vectorField edgeNorms(N);
  forAll(local_wall_WallsOutletEdges, i){
    label pointI = local_wall_WallsOutletEdges[i];
    edgePoints[i] = curWallBP[pointI];
    edgeNorms[i] = pointNs[pointI];
    edgeNorms[i]/=mag(edgeNorms[i]);
  }

  scalar minZ = min( edgePoints.component(vector::Z) );
  reduce(minZ, minOp<scalar>());
  
  // ***************************************************************************
  
  pointField outP = mesh_.boundaryMesh()[outletID].localPoints();
  
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
  
  scalar tol = 1.0;
  scalar tolerance = 0.0000000001;
  int itt = 0;
  while( tol>tolerance ){
    forAll(local_wall_WallsOutletEdges, i){
      displacement[i].z() = minZ - edgePoints[i].z();
    }
    vectorField proj_disp = transform(I - edgeNorms*edgeNorms/8.0, displacement);
    
    scalarField aux_f = mag(proj_disp);
    scalarField aux_f0(N, 0.0);

    if( aux_f == aux_f0 ){
      tol = 0.0;
    }
    else{
      tol = average( aux_f );
    }
    reduce(tol, sumOp<scalar>());
    tol /= static_cast<double>( Pstream::nProcs() );
    
    if(tol>tolerance){
      edgePoints += proj_disp;
    }
    else{
      edgePoints += displacement;
    }
    
    if(itt%100==0){
      Info<<"  Wall-outlet iter "<<itt<<"  tolerance: "<<tol<< nl;
    }

    itt+=1;
  }
  
  forAll(local_wall_WallsOutletEdges, i){
    label pointI = local_wall_WallsOutletEdges[i];
    wallDispl[ pointI ] = edgePoints[i]-wallBP[pointI];
  }
  
  
  /*
  forAll(local_wall_WallsOutletEdges, i){
    vector A = wallDispl[ local_wall_WallsOutletEdges[i] ];
    vector dz(0.0, 0.0, maxdZ - A.z());
          
    wallDispl[ local_wall_WallsOutletEdges[i] ] += dz;
  }
  */

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
void DissolMeshRlx::fixEdgeConcentration1(vectorField& concNorm, scalarField& inlC, 
        vectorField& poNorm, scalarField& walC){
  
  //const List<face>& llf = mesh_.boundaryMesh()[wallID].localFaces();
  const pointField& wallBoundaryPoints = mesh_.boundaryMesh()[wallID].localPoints();
  const pointField& boundaryPoints = mesh_.boundaryMesh()[inletID].localPoints();
  const labelListList& wpfl = mesh_.boundaryMesh()[wallID].pointFaces();
  const labelListList& ipfl = mesh_.boundaryMesh()[inletID].pointFaces();
  const pointField& wfaceCs = mesh_.boundaryMesh()[wallID].faceCentres();
  const pointField& faceCs = mesh_.boundaryMesh()[inletID].faceCentres();
  

  int N = local_inlet_WallsInletEdges.size();
  scalarList faceToPointSumWeights(N, 0.0);
  scalarList wfaceToPointSumWeights(N, 0.0);
  scalarField conc(N, 0.0);
  scalarField iconc(N, 0.0);
  scalarField wconc(N, 0.0);
  scalarField idelt(N, 0.0);
  forAll(local_inlet_WallsInletEdges, i){
    label wallpi = local_wall_WallsInletEdges[i];
    label pointi = local_inlet_WallsInletEdges[i];
    
    const labelList& curFaces = ipfl[pointi];
    forAll(curFaces, facei)
    {
      scalar pw = 1.0 / mag(faceCs[curFaces[facei]] - boundaryPoints[pointi]);
      
      vector fwd = faceCs[curFaces[facei]] - boundaryPoints[pointi];
      
      scalar delta = mag(fwd & poNorm[wallpi]);
      
      //scalar factor = 0.5/(0.5+delta);
      //scalar nnC = (1-factor)/delta * inlC[curFaces[facei]];
      scalar nnC = inlC[curFaces[facei]];
      
      iconc[i] += pw * nnC;
      
      idelt[i] += pw * delta;
      
      faceToPointSumWeights[i] += pw;
    }
    
    const labelList& wCurFaces = wpfl[wallpi];
    forAll(wCurFaces, facei)
    {
      scalar pw = 1.0 / mag(wfaceCs[wCurFaces[facei]] - wallBoundaryPoints[wallpi]);
      scalar nnC = walC[wCurFaces[facei]];
      wconc[i]+= pw * nnC;
      wfaceToPointSumWeights[i] += pw;
    }

  
  }
  // synchronization over coupled boundaries
  //syncTools::syncPointList(mesh_, global_WallInletEdges, conc, plusEqOp<scalar>(), 0.0);
  //syncTools::syncPointList(mesh_, global_WallInletEdges, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
  
  syncTools::syncPointList(mesh_, global_WallInletEdges, iconc, plusEqOp<scalar>(), 0.0);
  syncTools::syncPointList(mesh_, global_WallInletEdges, idelt, plusEqOp<scalar>(), 0.0);
  syncTools::syncPointList(mesh_, global_WallInletEdges, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
  

  syncTools::syncPointList(mesh_, global_WallInletEdges, wconc, plusEqOp<scalar>(), 0.0);
  syncTools::syncPointList(mesh_, global_WallInletEdges, wfaceToPointSumWeights, plusEqOp<scalar>(), 0.0);
  
  // normalization
  forAll(iconc, i){
    //conc[i] /= faceToPointSumWeights[i];
    iconc[i] /= faceToPointSumWeights[i];
    idelt[i] /= faceToPointSumWeights[i];
    wconc[i] /= wfaceToPointSumWeights[i];
    
  }
  
  
  forAll( inletTriple, i ){
    const labelList& currentTriple = inletTriple[i];
    
    vector cN0 = concNorm[ currentTriple[0] ];
    scalar c0 = mag( cN0 );
    vector cN1 = concNorm[ currentTriple[1] ];
    scalar c1 = mag( cN1 );
    vector cN2 = concNorm[ currentTriple[2] ];
    scalar c2 = mag( cN2 );
    
    vector newN0 = extrapolateVectorLinear(wallBoundaryPoints, cN1, cN2, currentTriple);
    newN0 /= mag(newN0);
    scalar newC0 = extrapolateConcentrationLinearZ(wallBoundaryPoints, c1, c2, currentTriple);
    
    
    coupledPatchInterpolation patchInterpolator( mesh_.boundaryMesh()[wallID], mesh_ );
    scalarField powaC = patchInterpolator.faceToPointInterpolate(walC);
    scalar c11 = powaC[currentTriple[1]];
    scalar c22 = powaC[currentTriple[2]];
    scalar newC00 = extrapolateConcentrationLinearZ(wallBoundaryPoints, c11, c22, currentTriple);
    
    
    conc[i] = (iconc[i]-newC00)/idelt[i]/2.;
    
    Info<< i<< " AA  "<< conc[i] 
            <<"   "<< newC0
            <<"   "<< idelt[i]
            <<nl;
    
    //concNorm[ currentTriple[0] ] = newC0 * newN0;
    //concNorm[ currentTriple[0] ] = newC0 * cN0 / c0;
    //concNorm[ currentTriple[0] ] = newN0 * conc[i];
    concNorm[ currentTriple[0] ] = conc[i]* cN0 / c0;
  }
}



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
    //scalar newC0 = extrapolateConcentrationLinearZ(boundaryPoints, c1, c2, currentTriple);
    scalar newC0 = extrapolateConcentrationExpZ(boundaryPoints, c1, c2, currentTriple);
    vector newN0 = extrapolateVectorLinear(boundaryPoints, cN1, cN2, currentTriple);
    newN0 /= mag(newN0);
    
    concNorm[ currentTriple[0] ] = newC0 * newN0;
    //concNorm[ currentTriple[0] ] = newC0 * cN0 / c0;
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
    if( findIndex(outletToAll, lW) != -1){
      local_wall_WallsOutletEdges.append( i );
      global_WallOutletEdges.append( lW );
    }
  }
}

float DissolMeshRlx::get_version(){
  return version;
}

