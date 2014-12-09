/* 
 * Surface relaxation for dissolFoam
 * */

#include "DissolMeshRlx.H"
#include <algorithm>

#include "Pstream.H"

#include "pyramidPointFaceRef.H"
#include "syncTools.H"

#include "transformField.H"
#include "symmTransformField.H"
#include "primitivePatchInterpolationSync.H"


// mesh1 is the mesh at time 0
DissolMeshRlx::DissolMeshRlx( const fvMesh& mesh, const fvMesh& mesh1)
:
  version(0.2),
  mesh_(mesh),
  mesh1_(mesh1)
{
  /*
   *   TOOOOODOOOOOOO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   *   periodicx !!!!!!!!!!!!!!!!!!!
   * 
   * 
   * 
   * 
   */
  
  
  // get ID of each patch we need
  wallID   = mesh_.boundaryMesh().findPatchID("walls");
  inletID  = mesh_.boundaryMesh().findPatchID("inlet");
  outletID = mesh_.boundaryMesh().findPatchID("outlet");
  cyclicID1 = mesh_.boundaryMesh().findPatchID("periodicx1");
  cyclicID2 = mesh_.boundaryMesh().findPatchID("periodicx2");

  // map vetex ID: patch to global
  wallsToAll  = mesh_.boundaryMesh()[wallID].meshPoints();
  inletToAll  = mesh_.boundaryMesh()[inletID].meshPoints();
  outletToAll = mesh_.boundaryMesh()[outletID].meshPoints();
  cyclicToAll1 = mesh_.boundaryMesh()[cyclicID1].meshPoints();
  cyclicToAll2 = mesh_.boundaryMesh()[cyclicID2].meshPoints();

  // map vertex ID: global to patch
  forAll(wallsToAll,  i){  allToWalls[wallsToAll[i]]   = i;  }
  forAll(inletToAll,  i){  allToInlet[inletToAll[i]]   = i;  }
  forAll(outletToAll, i){  allToOutlet[outletToAll[i]] = i;  }
  forAll(cyclicToAll1, i){  allToCyclic1[cyclicToAll1[i]] = i;  }
  forAll(cyclicToAll2, i){  allToCyclic2[cyclicToAll2[i]] = i;  }

  setUpLists();
  setUpPairsRlx();
  setUpPairsConc();
  
  //scalar pn = 3*(1+Pstream::myProcNo());
  //Pout<< "proc: " << Pstream::myProcNo() << "   master?: "<< Pstream::master() << nl;
  //Pout<< "wallsToAll: " << wallsToAll.size() << "\n";
  //reduce(pn, maxOp<scalar>());
  //Pout<< "After pn: " << pn  << nl;
  //std::exit(0);
}


void DissolMeshRlx::checkFacePyramidsMy()
{
    const scalar minPyrVol = -1e-20;

    // check whether face area vector points to the cell with higher label
    const vectorField& ctrs = mesh_.cellCentres();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const faceList& f = mesh_.faces();
    const pointField& p = mesh_.points();

    label nErrorPyrs = 0;

    forAll (f, faceI){
        // Create the owner pyramid - it will have negative volume
        scalar pyrVol = pyramidPointFaceRef(f[faceI], ctrs[own[faceI]]).mag(p);

        if (pyrVol > -minPyrVol) nErrorPyrs++;

        if (mesh_.isInternalFace(faceI)){
            // Create the neighbour pyramid - it will have positive volume
            scalar pyrVol =
                pyramidPointFaceRef(f[faceI], ctrs[nei[faceI]]).mag(p);

            if (pyrVol < minPyrVol) nErrorPyrs++;
        }
    }

    reduce(nErrorPyrs, sumOp<label>());

    if (nErrorPyrs > 0)
      Info<< " ***Error in face pyramids: "
          << nErrorPyrs << " faces are incorrectly oriented."
          << endl;
    else
      Info<< "    Face pyramids OK." << endl;
}


void DissolMeshRlx::boundaryCheck(){
  const pointField& loc_points1 = mesh_.boundaryMesh()[cyclicID1].localPoints();
  const pointField& loc_points2 = mesh_.boundaryMesh()[cyclicID2].localPoints();
  
  for(std::map<int,int>::iterator itr  = cyclicWallToWall.begin();
                                  itr != cyclicWallToWall.end();
                                ++itr ){
    //Info << loc_points[pnt0] << "  : "<< aaa<<"   : "<<conc[pnt0]<< nl;
    scalar dist = std::sqrt( (loc_points1[itr->first].y()-loc_points2[itr->second].y())
                             *
                             (loc_points1[itr->first].y()-loc_points2[itr->second].y())
                             +
                             (loc_points1[itr->first].z()-loc_points2[itr->second].z())
                             *
                             (loc_points1[itr->first].z()-loc_points2[itr->second].z())
                            );
    
    if( dist>0.0000001 ){
      Info<<" r+["<<itr->first<<"]: "<< loc_points1[itr->first]
          <<"   r-["<<itr->second<<"]: "<< loc_points2[itr->second]
          << "  dist: "<< dist
          <<endl;
    }
  }
  //std::exit(0);
}

// @TODO it is wrong 
void DissolMeshRlx::boundaryFix( pointField& newPoints ){
  for(std::map<int,int>::iterator itr  = cyclicWallToWall.begin();
                                  itr != cyclicWallToWall.end();
                                ++itr ){
    point averPoint = (newPoints[itr->first]+newPoints[itr->second])/2.;
    newPoints[itr->first] = averPoint;
    newPoints[itr->second] = averPoint;
  }
}






vectorField DissolMeshRlx::wallRelaxation1(){
  const labelList& meshPoints = mesh_.boundaryMesh()[wallID].meshPoints();

  const pointField& boundaryPoints = mesh_.boundaryMesh()[wallID].localPoints();
  const labelListList& plistFaces = mesh_.boundaryMesh()[wallID].pointFaces();
  const pointField& faceCs = mesh_.boundaryMesh()[wallID].faceCentres();
  scalar alpha = 0.5;
  
  vectorField sumDistance( boundaryPoints.size() );
  scalarField numberNeighb( boundaryPoints.size() );
  scalarField shlp( boundaryPoints.size(), 1.0 );
  
  forAll(boundaryPoints, i){
    point curP = boundaryPoints[i];
    const labelList& pFaces = plistFaces[i];
    
    numberNeighb[i] = plistFaces[i].size();
    
    vector sumw(0,0,0);
    forAll(pFaces, j){
      label faceI = pFaces[j];
      point faceC = faceCs[faceI];
      
      vector d = faceC - curP;
      scalar dist = mag(d);
      
      sumw += d;
    }
    
    sumDistance[i] = sumw;
    
    shlp[i] = numberNeighb[i];
  }
  
  syncTools::syncPointList( mesh_, meshPoints, sumDistance, plusEqOp<vector>(), vector::zero);
  syncTools::syncPointList( mesh_, meshPoints, shlp, plusEqOp<scalar>(), 0.0);

  vectorField averDistance( boundaryPoints.size() );
  
  forAll(averDistance, i){
    averDistance[i] = sumDistance[i] / shlp[i];
  }
  
  /*
  forAll(shlp ,i){
    if( shlp[i]<4 ){
      point curP = boundaryPoints[i];
      Pout << shlp[i] << "  " << curP << nl;
    }
  }
  std::exit(0);
   */
   
  
  vectorField displacement( boundaryPoints.size() );
  
  forAll(displacement, i){
    displacement[i] =  averDistance[i];
  }
  
  vectorField fNorm = mesh_.boundaryMesh()[wallID].faceNormals();
  primitivePatchInterpolationSync patchInterpolator( mesh_.boundaryMesh()[wallID], mesh_ );
  vectorField pointNorm = patchInterpolator.faceToPointInterpolate(fNorm);
  vectorField ttt = transform(I - pointNorm*pointNorm, displacement);

  // project back inlet
  vector n(0,0,1);
  vectorField inletEdgeD( lcl_wall_list_wallsInletEdges.size() );
  for(std::list<int>::iterator iter  = lcl_wall_list_wallsInletEdges.begin();
                               iter != lcl_wall_list_wallsInletEdges.end();
                             ++iter )
  {
    label id = std::distance(lcl_wall_list_wallsInletEdges.begin(), iter);
    inletEdgeD[id] = ttt[*iter];
  }
  vectorField inletPrj = transform(I - n*n, inletEdgeD);
  for(std::list<int>::iterator iter  = lcl_wall_list_wallsInletEdges.begin();
                               iter != lcl_wall_list_wallsInletEdges.end();
                             ++iter )
  {
    label id = std::distance(lcl_wall_list_wallsInletEdges.begin(), iter);
    ttt[*iter] = inletPrj[id];
  }
  // project back outlet
  vectorField outletEdgeD( lcl_wall_list_wallsOutletEdges.size() );
  for(std::list<int>::iterator iter  = lcl_wall_list_wallsOutletEdges.begin();
                               iter != lcl_wall_list_wallsOutletEdges.end();
                             ++iter )
  {
    label id = std::distance(lcl_wall_list_wallsOutletEdges.begin(), iter);
    outletEdgeD[id] = ttt[*iter];
  }
  vectorField outletPrj = transform(I - n*n, outletEdgeD);
  for(std::list<int>::iterator iter  = lcl_wall_list_wallsOutletEdges.begin();
                               iter != lcl_wall_list_wallsOutletEdges.end();
                             ++iter )
  {
    label id = std::distance(lcl_wall_list_wallsOutletEdges.begin(), iter);
    ttt[*iter] = outletPrj[id];
  }
  
  return ttt; //displacement;
}



/*
vectorField DissolMeshRlx::wallRelaxation1(){
  const labelList& meshPoints = mesh_.boundaryMesh()[wallID].meshPoints();

  const pointField& boundaryPoints = mesh_.boundaryMesh()[wallID].localPoints();
  const labelListList& plistFaces = mesh_.boundaryMesh()[wallID].pointFaces();
  const pointField& faceCs = mesh_.boundaryMesh()[wallID].faceCentres();
  scalar alpha = 0.25;
  
  scalarField sumDistance( boundaryPoints.size() );
  scalarField numberNeighb( boundaryPoints.size() );
  
  forAll(boundaryPoints, i){
    point curP = boundaryPoints[i];
    const labelList& pFaces = plistFaces[i];
    
    numberNeighb[i] = plistFaces.size();
    
    scalar sumw = 0.0;
    forAll(pFaces, j){
      label faceI = pFaces[j];
      point faceC = faceCs[faceI];
      
      vector d = faceC - curP;
      scalar dist = mag(d);
      
      sumw += dist;
    }
    
    sumDistance[i] = sumw;
  }
  
  syncTools::syncPointList( mesh_, meshPoints, sumDistance, plusEqOp<scalar>(), 0.0);
  syncTools::syncPointList( mesh_, meshPoints, numberNeighb, plusEqOp<scalar>(), 0.0);

  scalarField averDistance( boundaryPoints.size() );
  
  forAll(averDistance, i){
    averDistance[i] = sumDistance[i] / numberNeighb[i];
  }
  
  vectorField displacement( boundaryPoints.size() );
  
  forAll(displacement, i){
    point curP = boundaryPoints[i];
    const labelList& pFaces = plistFaces[i];
    
    vector displ(0,0,0);
    forAll(pFaces, j){
      label faceI = pFaces[j];
      point faceC = faceCs[faceI];
      
      vector d = faceC - curP;
      scalar dist = mag(d);
      
      vector aux = d * (1 - averDistance[i] / dist);
      displ = (magSqr(aux)>=magSqr(displ) ? aux : displ);
    }
    
    displacement[i] = displ;
  }
  
  syncTools::syncPointList( mesh_, meshPoints, displacement, maxMagSqrEqOp<vector>(), vector::zero);

  forAll(displacement, i){
    displacement[i] *= alpha;
  }
  forAll(displacement, i){
    if( std::find(lcl_wall_list_wallsInletEdges.begin(), lcl_wall_list_wallsInletEdges.end(), i)
        != lcl_wall_list_wallsInletEdges.end()
            ||
        std::find(lcl_wall_list_wallsOutletEdges.begin(), lcl_wall_list_wallsOutletEdges.end(), i)
        != lcl_wall_list_wallsOutletEdges.end()
    ){
      displacement[i] = vector::zero;
    }
  }
  
  return displacement;
}
*/






vectorField DissolMeshRlx::wallRelaxation(int N){
  vectorField boundaryCopy = mesh_.boundaryMesh()[wallID].localPoints();
  
  for(int i=0; i<N; i++){
  
    vectorField boundaryRelax(mesh_.boundaryMesh()[wallID].localPoints().size(), vector::zero);

    for(std::map<int, std::pair<int, int> >::iterator iter  = neighbZrlx.begin();
                              iter != neighbZrlx.end();
                                ++iter){
      if(iter->second.first!=-1 && iter->second.second!=-1){
        point A0 = boundaryCopy[iter->first];
        point A1 = boundaryCopy[ allToWalls[iter->second.first] ];
        point A2 = boundaryCopy[ allToWalls[iter->second.second] ];
        boundaryRelax[iter->first] = vertexDispl(A0, A1, A2);
        
        /*
        if( A0.z()<24.1 && A0.z()>23.9 ){
          Info<< iter->first<< " " << A0 
                  << " | " << iter->second.first << "  " << iter->second.second << "  " 
                  << A1 << "  " << A2
                  << "  | bound relax " << boundaryRelax[iter->first]
                  << nl;
        }
        */
        
        
      }
    }
    for(std::map<int, std::pair<int, int> >::iterator iter  = neighbXrlx.begin();
                              iter != neighbXrlx.end();
                                ++iter){
      if(iter->second.first!=-1 && iter->second.second!=-1){
        point A0 = boundaryCopy[ iter->first ];
        point A1 = boundaryCopy[ allToWalls[iter->second.first] ];
        point A2 = boundaryCopy[ allToWalls[iter->second.second] ];
        boundaryRelax[iter->first] += vertexDispl(A0, A1, A2);
      }
    }
    
    //Foam::polyMesh& meshParent = refCast<Foam::polyMesh>(mesh_);
    const labelList& meshPoints = mesh_.boundaryMesh()[wallID].meshPoints();
    syncTools::syncPointList( mesh_, meshPoints, boundaryRelax, plusEqOp<vector>(), vector::zero);
      
    boundaryCopy += boundaryRelax;
  }
  
  vectorField finalDisplacements = boundaryCopy - mesh_.boundaryMesh()[wallID].localPoints();
  
  return finalDisplacements;
}



vectorField DissolMeshRlx::wallRelaxation(int N, const vectorField& wallPos){
  //vectorField boundaryCopy = mesh_.boundaryMesh()[wallID].localPoints();
  vectorField boundaryCopy = wallPos;
  
  for(int i=0; i<N; i++){
  
    vectorField boundaryRelax(wallPos.size(), vector::zero);

    for(std::map<int, std::pair<int, int> >::iterator iter  = neighbZrlx.begin();
                              iter != neighbZrlx.end();
                                ++iter){
      if(iter->second.first!=-1 && iter->second.second!=-1){
        point A0 = boundaryCopy[iter->first];
        point A1 = boundaryCopy[ allToWalls[iter->second.first] ];
        point A2 = boundaryCopy[ allToWalls[iter->second.second] ];
        boundaryRelax[iter->first] = vertexDispl(A0, A1, A2);
        
        /*
        if( A0.z()<24.1 && A0.z()>23.9 ){
          Info<< iter->first<< " " << A0 
                  << " | " << iter->second.first << "  " << iter->second.second << "  " 
                  << A1 << "  " << A2
                  << "  | bound relax " << boundaryRelax[iter->first]
                  << nl;
        }
        */
        
        
      }
    }
    for(std::map<int, std::pair<int, int> >::iterator iter  = neighbXrlx.begin();
                              iter != neighbXrlx.end();
                                ++iter){
      if(iter->second.first!=-1 && iter->second.second!=-1){
        point A0 = boundaryCopy[ iter->first ];
        point A1 = boundaryCopy[ allToWalls[iter->second.first] ];
        point A2 = boundaryCopy[ allToWalls[iter->second.second] ];
        boundaryRelax[iter->first] += vertexDispl(A0, A1, A2);
      }
    }
    
    boundaryCopy += boundaryRelax;
  }
  
  vectorField finalDisplacements = boundaryCopy - wallPos;
  
  return finalDisplacements;
}











/*#######################################################################################
 *  * A0 - boundary which has to be moved.
 *  * A1 - the point above A0
 *  * A2 - the point below A0
 *  * it returns a vector displacement
 *#######################################################################################*/
vector DissolMeshRlx::vertexDisplZ( point A0, point A1, point A2 ){
  scalar newZ = ( A1.z() + A2.z() ) / 2.0;
  newZ = (newZ + A0.z())/2.0;
  point refA;
  if( newZ < A0.z() ){
    refA = A2;
  }
  else{
    refA = A1;
  }
  
  vector A0I = refA - A0;
  
  scalar z0I = newZ - A0.z();
  
  //if(A0I.z() == 0){
  //  Info<< " WARNING!!! A0I.z()=0"<< nl;
  //}
  
  scalar z0IdivG = z0I / A0I.z();
  
  vector r0( A0I.x()*z0IdivG,
             A0I.y()*z0IdivG,
                     z0I);
  return r0;
}

/*#######################################################################################
 *  * the relaxation on the edge in order to keep uniform distribution over X direction
 *  * A0 - central point
 *  * A1 - point with x > x0
 *  * A2 - point with x < x0
 *#######################################################################################*/
vector DissolMeshRlx::vertexDisplX( point A0, point A1, point A2 ){
  scalar newX = ( A1.x() + A2.x() ) / 2.0;
  newX = (newX + A0.x())/2.0;
  
  point refA;
  if( newX < A0.x() ){
    refA = A2;
  }
  else{
    refA = A1;
  }
  
  vector A0I = refA - A0;
  
  scalar x0I = newX - A0.x();
  scalar x0IdivG = x0I / A0I.x(); 
  
  vector r0(             x0I,
             A0I.y()*x0IdivG,
             A0I.z()*x0IdivG);
  return r0;
}

/*#######################################################################################
 *  * the relaxation on the edge in order to keep uniform distribution
 *  * A0 - central point
 *  * A1, A2 - neighbors
 *#######################################################################################*/
vector DissolMeshRlx::vertexDispl( point A0, point A1, point A2 ){

  vector d1  = A1 - A0;
  vector d2  = A2 - A0;
  vector d12 = A1 - A2;
  scalar b   = (d1+d2) & d12;
  scalar a   = fabs(b)/((d12 & d12) + fabs(b));
  
  scalar alpha = 0.5;
  
  vector r0;
  if (b > 0)
    r0 = alpha * a*d1;
  else
    r0 = alpha * a*d2;
                
  return r0;
}


vectorField DissolMeshRlx::calculateInletDisplacement(vectorField& wallDispl){

  scalar maxdZ = 0;
  for(std::list<int>::iterator iter  = lcl_wall_list_wallsInletEdges.begin();
                               iter != lcl_wall_list_wallsInletEdges.end();
                             ++iter ){
    vector A = wallDispl[*iter];
    maxdZ = max( maxdZ, A.z() );    
  }

  reduce(maxdZ, maxOp<scalar>());

  for(std::list<int>::iterator iter  = lcl_wall_list_wallsInletEdges.begin();
                               iter != lcl_wall_list_wallsInletEdges.end();
                             ++iter ){
    vector A = wallDispl[*iter];
    vector dz(0.0, 0.0, maxdZ - A.z());
          
    wallDispl[*iter] += dz;
  }

  vectorField pointDispInlet( mesh_.boundaryMesh()[inletID].localPoints().size(), vector::zero );
      
  pointDispInlet.replace( vector::Z, maxdZ);
  pointDispInlet.replace( vector::Y, 0);
  pointDispInlet.replace( vector::X, 0);

  return pointDispInlet;
}

void DissolMeshRlx::fixWallDisplPeriodic( vectorField& pointDispWall ){
  for(std::map<int,int>::iterator itr  = wallCyclicEdgeToEdge.begin();
                                  itr != wallCyclicEdgeToEdge.end();
                                ++itr ){
    /*
    Info<<" r1["<<itr->first<<"]: "<< pointDispWall[itr->first]
        <<" r2["<<itr->second<<"]: "<< pointDispWall[itr->second]
        <<endl;
    */
    
    vector r_av = (pointDispWall[itr->first] + pointDispWall[itr->second]) * 0.5;
    r_av.x() = 0.0;
    pointDispWall[itr->first] = r_av;
    pointDispWall[itr->second] = r_av;
  }
}



void DissolMeshRlx::fixEdgeConcentration( scalarField& conc ){
  // exponential extrapolation of the concentration on the edge vertexes
  const pointField& loc_points = mesh_.boundaryMesh()[wallID].localPoints();
  // edgeConcentrationFixMap
  for(std::map<int, std::pair<int, int> >::iterator itr  = edgeConcentrationFixMap.begin();
                                                    itr != edgeConcentrationFixMap.end();
                                                  ++itr ){
    label pnt0 = itr->first;
    label pnt1 = mesh_.boundaryMesh()[wallID].whichPoint( (itr->second).first  );
    label pnt2 = mesh_.boundaryMesh()[wallID].whichPoint( (itr->second).second );
    
    //scalar aaa = conc[ pnt0 ];
    conc[ pnt0 ] = extrapolateConcentrationLinearZ(loc_points, conc, pnt0, pnt1, pnt2);
    //conc[ pnt0 ] = extrapolateConcentrationExp(loc_points, conc, pnt0, pnt1, pnt2);
    
    //Info << loc_points[pnt0] << "  : "<< aaa<<"   : "<<conc[pnt0]<< nl;
  }
  
/*
  for(std::map<int,int>::iterator itr  = wallCyclicEdgeToEdge.begin();
                                  itr != wallCyclicEdgeToEdge.end();
                                ++itr ){
    
    Info<<" c1["<<itr->first<<"]: "<< conc[itr->first]
        <<" c2["<<itr->second<<"]: "<< conc[itr->second]
        <<endl;
    
    scalar c_av = (conc[itr->first] + conc[itr->second]) * 0.5;
    conc[itr->first] = c_av;
    conc[itr->second] = c_av;
  }
*/
}


void DissolMeshRlx::fixCyclicNormal( vectorField& norm, scalarField& conc ){
  for(std::map<int,int>::iterator itr  = wallCyclicEdgeToEdge.begin();
                                  itr != wallCyclicEdgeToEdge.end();
                                ++itr ){
    /*
    vector av = norm[itr->first] + norm[itr->second];
    vector n_av = av / mag(av);
    norm[itr->first] = n_av;
    norm[itr->second] = n_av;
    
    scalar c_av = (conc[itr->first] + conc[itr->second]) * 0.5;
    conc[itr->first] = c_av;
    conc[itr->second] = c_av;
    */
    norm[itr->first] *= 0.5;
    norm[itr->second] *= 0.5;
    conc[itr->first] *= 0.5;
    conc[itr->second] *= 0.5;
  }

  /*
  forAll(cornerPoints, i){
    vector nnn = norm[cornerPoints[i]];
    vector& corner_n = norm[cornerPoints[i]];
    corner_n.z() = 0.0;
    corner_n.x() = 0.0;
    corner_n /= mag( corner_n );
    //Info << " Corner points! " << norm[cornerPoints[i]] << "  and before: " << nnn << nl;
  }
  */
  
  //std::exit(0);
}

scalar DissolMeshRlx::extrapolateConcentrationExp(const pointField& loc_points,
                                const scalarField& point_conc,
                                label pnt0, label pnt1, label pnt2){
  scalar r02 = mag(  loc_points[pnt0] - loc_points[pnt2]  );
  scalar r12 = mag(  loc_points[pnt1] - loc_points[pnt2]  );
  scalar c1 = point_conc[pnt1];
  scalar c2 = point_conc[pnt2];

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
                                const scalarField& point_conc,
                                label pnt0, label pnt1, label pnt2){
  scalar r02 = loc_points[pnt0].z() - loc_points[pnt2].z();
  scalar r12 = loc_points[pnt1].z() - loc_points[pnt2].z();
  scalar c1 = point_conc[pnt1];
  scalar c2 = point_conc[pnt2];
  return c2 - (c2-c1) * r02/r12;
}


void DissolMeshRlx::setUpPairsConc(){

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // this is the list for extrapolation of the C value to the vertices on the edge 
  // between side wall and inlet. Pair first is the next vertex in z direction, pair second
  // is the next to the next neighbor in z direction
  // @TODO move to function +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  forAll (pe, indp){
    // list of edges  
    labelList elist = pe[indp];
    // check if edge belong to the surface
    label neighb     = -1; // near neighbor in z
    label nextneighb = -1; // next near neighbor in z
    vector posP = mesh1_.boundaryMesh()[wallID].localPoints()[indp];// position of the vertex we check

    // this all just for first layer of vertices in z direction
    //if( posP.z()<0.005 ){
    //if( lcl_wall_list_wallsInletEdges.find(indp) != lcl_wall_list_wallsInletEdges.end() ){
    if( std::find(lcl_wall_list_wallsInletEdges.begin(), lcl_wall_list_wallsInletEdges.end(), indp)
        !=
        lcl_wall_list_wallsInletEdges.end()
    ){
      forAll (elist, ind){
        label segside = mesh1_.edges()[elist[ind]].start();
        if(segside==wallsToAll[indp]) segside = mesh1_.edges()[elist[ind]].end();
        
        label neibIDonBound = mesh1_.boundaryMesh()[wallID].whichPoint(segside);
        if ( neibIDonBound!=-1 ){
          // means segside is next neighbor of indp on the boundary    
          vector posN = mesh1_.points()[segside];// position of a neighbor
          
          // @TODO optimize, it may not work in general case
          if( posN.z()-posP.z()>0.01 ){
            neighb = segside;
            
            labelList nelist = pe[neibIDonBound];
            // @TODO the same for neighb, DO it recursive lazy boy!!
            forAll (nelist, ind11){
              label segside11 = mesh1_.edges()[nelist[ind11]].start();
              if(segside11==wallsToAll[neibIDonBound]) segside11 = mesh1_.edges()[nelist[ind11]].end();

              label neibIDonBound11 = mesh1_.boundaryMesh()[wallID].whichPoint(segside11);
              if ( neibIDonBound11!=-1 ){
                // means segside is next neighbor of neighb on the boundary
                vector posN11 = mesh1_.points()[segside11];// position of a neighbor

                // @TODO optimize, it may not work in general case
                if( posN11.z()-posN.z()>0.01 ){
                  nextneighb = segside11;
                  break;
                }
              }
            }
            
            break;
          }
        }
      }
      
      edgeConcentrationFixMap[indp] = std::make_pair(neighb, nextneighb);
    }
  }
  
  // to fix cyclic boundary
  
  for(std::list<int>::iterator iter  = lcl_wall_list_wallsCyclicEdges.begin();
                               iter != lcl_wall_list_wallsCyclicEdges.end();
                             ++iter ){
  
    for(std::list<int>::iterator iter1  = lcl_wall_list_wallsCyclicEdges.begin();
                                 iter1 != lcl_wall_list_wallsCyclicEdges.end();
                               ++iter1 ){

      if( iter != iter1 
          && ( std::abs(mesh_.boundaryMesh()[wallID].localPoints()[*iter].y() -
                  mesh_.boundaryMesh()[wallID].localPoints()[*iter1].y()) < 0.0001)
          && ( std::abs(mesh_.boundaryMesh()[wallID].localPoints()[*iter].z() -
                  mesh_.boundaryMesh()[wallID].localPoints()[*iter1].z()) < 0.0001)
          && search2IntMapByValue(wallCyclicEdgeToEdge, *iter) == wallCyclicEdgeToEdge.end()
          && search2IntMapByValue(wallCyclicEdgeToEdge, *iter1) == wallCyclicEdgeToEdge.end()
          && search2IntMapByKey(wallCyclicEdgeToEdge, *iter) == wallCyclicEdgeToEdge.end()
          && search2IntMapByKey(wallCyclicEdgeToEdge, *iter1) == wallCyclicEdgeToEdge.end()
        ){
        wallCyclicEdgeToEdge[*iter] = *iter1;
        break;
      }
      
    }
  }
  
  
  // fill the corner list at the outlet
  for(std::map<int,int>::iterator itr  = wallCyclicEdgeToEdge.begin();
                                  itr != wallCyclicEdgeToEdge.end();
                                ++itr ){
    label id1 = itr->first;
    label id2 = itr->second;

    label idAll1 = wallsToAll[id1];
    //label idAll2 = wallsToAll[id2];

    if( mesh1_.boundaryMesh()[outletID].whichPoint(idAll1)!=-1 ){
      cornerPoints.append( id1 );
      cornerPoints.append( id2 );
    }

  }
  

  
  
  
  pointField cyclic_points1 = mesh_.boundaryMesh()[cyclicID1].localPoints();
  pointField cyclic_points2 = mesh_.boundaryMesh()[cyclicID2].localPoints();
  forAll(cyclic_points1, i1){
    forAll(cyclic_points2, i2 ){
      if(
          ( std::abs(cyclic_points1[i1].y() -
                        cyclic_points2[i2].y()) < 0.0001)
          && ( std::abs(cyclic_points1[i1].z() -
                        cyclic_points2[i2].z()) < 0.0001)
          //&& search2IntMapByValue(cyclicWallToWall, i1) == cyclicWallToWall.end()
          && search2IntMapByValue(cyclicWallToWall, i2) == cyclicWallToWall.end()
          && search2IntMapByKey(cyclicWallToWall, i1) == cyclicWallToWall.end()
          //&& search2IntMapByKey(cyclicWallToWall, i2) == cyclicWallToWall.end()
              ){
        cyclicWallToWall[i1] = i2;
        break;
      }
    }
  }
  
}

std::map<int,int>::iterator DissolMeshRlx::search2IntMapByValue(std::map<int,int>& ii_map, int value){
  std::map<int,int>::iterator m_iter = ii_map.end();

  for (std::map<int,int>::iterator it = ii_map.begin(); it != ii_map.end(); ++it ){
    if (it->second == value){
       m_iter = it;
    }
  }
  
  return m_iter;
}

std::map<int,int>::iterator DissolMeshRlx::search2IntMapByKey(std::map<int,int>& ii_map, int key){
  return ii_map.find(key);
}


void DissolMeshRlx::setUpPairsRlx(){
  // global edge list
  const edgeList& e = mesh1_.edges();
  const polyBoundaryMesh& boundary_mesh_list = mesh1_.boundaryMesh();

  // set up storage for pointEdges
  List<SLList<label> > pointEdges(boundary_mesh_list[wallID].meshPoints().size());

  forAll (e, edgeI){
    label localID1 = boundary_mesh_list[wallID].whichPoint( e[edgeI].start() );
    label localID2 = boundary_mesh_list[wallID].whichPoint( e[edgeI].end() );
    // shift patch points
    if( localID1 >= 0 ) pointEdges[localID1].append(edgeI);
    if( localID2 >= 0 ) pointEdges[localID2].append(edgeI);
  }
 
  // sort out the list
  pe.setSize( pointEdges.size() );
  forAll (pointEdges, pointI){
    pe[pointI].setSize(pointEdges[pointI].size());

    label i = 0;
    for( SLList<label>::iterator curEdgesIter = pointEdges[pointI].begin();
         curEdgesIter != pointEdges[pointI].end();
         ++curEdgesIter, ++i ) {
      pe[pointI][i] = curEdgesIter();
    }
  }
  
  /*
  int poind = 0;
  if( Pstream::master() ){
    poind = 0;
  }
  else{
    poind = 1172;
  }
  forAll( pe[poind], i ){
    Pout<< "edges0: "<< e[pe[poind][i]]  <<nl;
    //Pout<< "edges0: "<< boundary_mesh_list[wallID].localPoints()[ e[pe[poind][i]][1] ] <<nl;
  }
  */

  forAll (pe, indp){
    // list of edges
    labelList el = pe[indp];
    // check if edge belong to the surface
    label intNeighb = -1; // internal neighbor
    label bouNeighb = -1; // boundary neighbor with z > z[indp]
    label bouNeighbLow = -1; // boundary neighbor with z < z[indp]
    
    // edge relaxation
    label nBX = -1; //neighbour with x > x0
    label nLX = -1; //neighbour with x < x0

    vector posP = boundary_mesh_list[wallID].localPoints()[indp];// position of the vertex we check

    label glob_lab = wallsToAll[indp];
    forAll (el, inde){
      label segside = mesh1_.edges()[el[inde]].start();
      if(segside==wallsToAll[indp]) segside = mesh1_.edges()[el[inde]].end();
      if ( boundary_mesh_list[wallID].whichPoint(segside)==-1 ){
        // means segside is next internal neighbor of indp
        intNeighb = segside;
      }
      else{
        vector posN = mesh1_.points()[segside];// position of a neighbor
        if( posN.z()-posP.z()>0.01){ //  && (posP.x()>-34.9) 
          bouNeighb = segside;
          //break;
        }
        if( posN.z()-posP.z()<-0.01){  // && (posP.x()>-34.9) 
            bouNeighbLow = segside;
        }
          
        // only for edge
        if(   ( allToWalls.find( glob_lab ) != allToWalls.end() )  
                &&
              ( std::abs(posN.z()-posP.z())<0.001 )
          ){
          if( posN.x()-posP.x()>0.001 ){
            nBX = segside;
          }
          if( posN.x()-posP.x()<-0.001 ){
            nLX = segside;
          }
        }
        
      }
    }
    if( allToInlet.find( glob_lab ) == allToInlet.end() && allToOutlet.find( glob_lab ) == allToOutlet.end() ){
      neighbZrlx[indp] = std::make_pair(bouNeighb, bouNeighbLow);
    }
    //if( allToCyclic.find( glob_lab ) == allToCyclic.end() ){
    if( 
            allToCyclic1.find( glob_lab ) == allToCyclic1.end()
            &&
            allToCyclic2.find( glob_lab ) == allToCyclic2.end()
      ){
      neighbXrlx[indp] = std::make_pair(nBX, nLX);
    }
    nBX = -1; nLX = -1;
  }
  
}

void DissolMeshRlx::setUpLists(){

  forAll(wallsToAll, ii){
    label lW = wallsToAll[ii];
    
    // create walls-inlet edge list of vertex IDs
    std::map<int,int>::iterator iI = allToInlet.find( lW );
    if(iI != allToInlet.end()){
      glb_list_wallsInletEdges.push_back( iI->second ); //element found;
      lcl_wall_list_wallsInletEdges.push_back( ii );
    }    
    
    // create walls-outlet edge list of vertex IDs
    std::map<int,int>::iterator iO = allToOutlet.find( lW );
    if(iO != allToOutlet.end()){
      glb_list_wallsOutletEdges.push_back( iO->second ); //element found;
      lcl_wall_list_wallsOutletEdges.push_back( ii );
    }

    // create walls-cyclic edge list of vertex IDs
    std::map<int,int>::iterator iC1 = allToCyclic1.find( lW );
    std::map<int,int>::iterator iC2 = allToCyclic2.find( lW );
    if( iC1 != allToCyclic1.end() ){
      glb_list_wallsCyclicEdges.push_back( iC1->second ); //element found;
      lcl_wall_list_wallsCyclicEdges.push_back( ii );
    }
    if( iC2 != allToCyclic2.end() ){
      glb_list_wallsCyclicEdges.push_back( iC2->second ); //element found;
      lcl_wall_list_wallsCyclicEdges.push_back( ii );
    }
    
    forAll(mesh1_.boundaryMesh(), patchI){
      if( mesh1_.boundaryMesh()[patchI].coupled() ){
        if( mesh1_.boundaryMesh()[patchI].whichPoint( lW ) != -1 ){
          edgePoints.append( ii );
        }
      }
    }

  }
  
}

float DissolMeshRlx::get_version(){
  return version;
}

