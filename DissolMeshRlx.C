/* 
 * Surface relaxation for dissolFoam
 * */

#include "DissolMeshRlx.H"
#include <algorithm>


// mesh1 is the mesh at time 0
DissolMeshRlx::DissolMeshRlx( const fvMesh& mesh, const fvMesh& mesh1)
:
  version(0.2),
  mesh_(mesh),
  mesh1_(mesh1)
{
          
  // get ID of each patch we need
  wallID   = mesh_.boundaryMesh().findPatchID("walls");
  inletID  = mesh_.boundaryMesh().findPatchID("inlet");
  outletID = mesh_.boundaryMesh().findPatchID("outlet");
  cyclicID = mesh_.boundaryMesh().findPatchID("periodicx");

  // map vetex ID: patch to global
  wallsToAll  = mesh_.boundaryMesh()[wallID].meshPoints();
  inletToAll  = mesh_.boundaryMesh()[inletID].meshPoints();
  outletToAll = mesh_.boundaryMesh()[outletID].meshPoints();
  cyclicToAll = mesh_.boundaryMesh()[cyclicID].meshPoints();

  // map vertex ID: global to patch
  forAll(wallsToAll,  i){  allToWalls[wallsToAll[i]]   = i;  }
  forAll(inletToAll,  i){  allToInlet[inletToAll[i]]   = i;  }
  forAll(outletToAll, i){  allToOutlet[outletToAll[i]] = i;  }
  forAll(cyclicToAll, i){  allToCyclic[cyclicToAll[i]] = i;  }

  setUpLists();
  setUpPairsRlx();
  setUpPairsConc();
  
  Pout<< "proc: " << Pstream::myProcNo() << "   master?: "<< Pstream::master() << "\n";
  Pout<< "wallsToAll: " << wallsToAll.size() << "\n";
}

void DissolMeshRlx::boundaryCheck(){
  const pointField& loc_points = mesh_.boundaryMesh()[cyclicID].localPoints();
  
  for(std::map<int,int>::iterator itr  = cyclicWallToWall.begin();
                                  itr != cyclicWallToWall.end();
                                ++itr ){
    //Info << loc_points[pnt0] << "  : "<< aaa<<"   : "<<conc[pnt0]<< nl;
    scalar dist = std::sqrt( (loc_points[itr->first].y()-loc_points[itr->second].y())
                             *
                             (loc_points[itr->first].y()-loc_points[itr->second].y())
                             +
                             (loc_points[itr->first].z()-loc_points[itr->second].z())
                             *
                             (loc_points[itr->first].z()-loc_points[itr->second].z())
                            );
    
    if( dist>0.0000001 ){
      Info<<" r+["<<itr->first<<"]: "<< loc_points[itr->first]
          <<"   r-["<<itr->second<<"]: "<< loc_points[itr->second]
          << "  dist: "<< dist
          <<endl;
    }
  }
  //std::exit(0);
}

void DissolMeshRlx::boundaryFix( pointField& newPoints ){
  for(std::map<int,int>::iterator itr  = cyclicWallToWall.begin();
                                  itr != cyclicWallToWall.end();
                                ++itr ){
    point averPoint = (newPoints[itr->first]+newPoints[itr->second])/2.;
    newPoints[itr->first] = averPoint;
    newPoints[itr->second] = averPoint;
  }
}


vectorField DissolMeshRlx::wallRelaxation(){
  vectorField boundaryRelax(mesh_.boundaryMesh()[wallID].localPoints().size(), vector::zero);

  for(std::map<int, std::pair<int, int> >::iterator iter  = neighbZrlx.begin();
						    iter != neighbZrlx.end();
					          ++iter){
    if(iter->second.first!=-1 && iter->second.second!=-1){
      point A0 = mesh_.boundaryMesh()[wallID].localPoints()[iter->first];
      point A1 = mesh_.points()[iter->second.first];
      point A2 = mesh_.points()[iter->second.second];
      boundaryRelax[iter->first] = vertexDisplZ(A0, A1, A2);
    }
  }

  for(std::map<int, std::pair<int, int> >::iterator iter  = neighbXrlx.begin();
						    iter != neighbXrlx.end();
					          ++iter){
    if(iter->second.first!=-1 && iter->second.second!=-1){
      point A0 = mesh_.boundaryMesh()[wallID].localPoints()[iter->first];
      point A1 = mesh_.points()[iter->second.first];
      point A2 = mesh_.points()[iter->second.second];
      boundaryRelax[iter->first] += vertexDisplX(A0, A1, A2);
    }
  }
  
  return boundaryRelax;
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



vectorField DissolMeshRlx::calculateInletDisplacement(vectorField& wallDispl){

  scalar maxdZ = 0;
  for(std::list<int>::iterator iter  = lcl_wall_list_wallsInletEdges.begin();
                               iter != lcl_wall_list_wallsInletEdges.end();
                             ++iter ){
    vector A = wallDispl[*iter];
    maxdZ = max( maxdZ, A.z() );    
  }


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
  
  pointField cyclic_points = mesh_.boundaryMesh()[cyclicID].localPoints();
  forAll(cyclic_points, i){
    forAll(cyclic_points, i1 ){
      if( i != i1
          && ( std::abs(cyclic_points[i].y() -
                        cyclic_points[i1].y()) < 0.0001)
          && ( std::abs(cyclic_points[i].z() -
                        cyclic_points[i1].z()) < 0.0001)
          && search2IntMapByValue(cyclicWallToWall, i) == cyclicWallToWall.end()
          && search2IntMapByValue(cyclicWallToWall, i1) == cyclicWallToWall.end()
          && search2IntMapByKey(cyclicWallToWall, i) == cyclicWallToWall.end()
          && search2IntMapByKey(cyclicWallToWall, i1) == cyclicWallToWall.end()
              ){
        cyclicWallToWall[i] = i1;
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
        if( posN.z()-posP.z()>0.01 ){
          bouNeighb = segside;
          //break;
        }
        if( posN.z()-posP.z()<-0.01 ){
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
    if( allToCyclic.find( glob_lab ) == allToCyclic.end() ){
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
    std::map<int,int>::iterator iC = allToCyclic.find( lW );
    if(iC != allToCyclic.end()){
      glb_list_wallsCyclicEdges.push_back( iC->second ); //element found;
      lcl_wall_list_wallsCyclicEdges.push_back( ii );
    }

  }

}

float DissolMeshRlx::get_version(){
  return version;
}

