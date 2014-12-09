/*
 *  This is a class which containes all the information about current dissolFoam system
 *  and can be used for different relaxation processes.
 * 
 * 
 */

#ifndef DissolMeshRlx_H
#define DissolMeshRlx_H

#include "fvCFD.H"
#include <list>
#include <map>

class DissolMeshRlx{

private:

  const float version;

  const fvMesh& mesh_;

  const fvMesh& mesh1_;
  
  // ID of each patch
  label wallID, inletID, outletID, cyclicID1, cyclicID2;

  // vetex ID maps: patch to global
  labelList wallsToAll, inletToAll, outletToAll, cyclicToAll1, cyclicToAll2;

  // oposite to previous lists. Map IDs global to patch
  std::map<int,int> allToInlet, allToOutlet, allToCyclic1, allToCyclic2, allToWalls;

  // list of global IDs for walls- inlet and outlet edges
  std::list<int> glb_list_wallsInletEdges;
  std::list<int> glb_list_wallsOutletEdges;
  std::list<int> glb_list_wallsCyclicEdges;
  // list of local (in terms of wall) IDs for walls- inlet and outlet edges
  std::list<int> lcl_wall_list_wallsInletEdges;
  std::list<int> lcl_wall_list_wallsOutletEdges;
  std::list<int> lcl_wall_list_wallsCyclicEdges;


  labelListList pe;

  labelList cornerPoints; // belong to outlet, wall and periodic. Stores wall point ID
  
  
  
  labelList edgePoints; // when a point belongs to the wall and any coupled patch. Stores wall point ID
  
  
  
  /*
   * main maps used for wall surface relaxation, approximatelly each vertex on the wall
   * initially looks like:
   *    |
   *   -*-
   *    |
   * It has 4 neighbors on the wall
   */   
  std::map<int, std::pair<int, int> > neighbZrlx; // map to top-bottom neighbor IDs for relaxation in Z direction
  std::map<int, std::pair<int, int> > neighbXrlx; // map to left-right neighbor IDs for relaxation in X direction

  // map for the extrapolation of the concentration on the edge
  // it fixes inlet boundary condition inconsistance
  std::map<int, std::pair<int, int> > edgeConcentrationFixMap;

  // map of corresponding cyclic-walls edges
  std::map<int, int> wallCyclicEdgeToEdge;
  std::map<int, int> cyclicWallToWall;
  
  std::map<int,int>::iterator search2IntMapByValue(std::map<int,int>&, int);
  std::map<int,int>::iterator search2IntMapByKey(std::map<int,int>&, int);

  scalar extrapolateConcentrationExp(const pointField&, const scalarField&, label, label, label);
  scalar extrapolateConcentrationLinear(const pointField&, const scalarField&, label, label, label);
  scalar extrapolateConcentrationLinearZ(const pointField&, const scalarField&, label, label, label);

public:
  
  // Constructors
  DissolMeshRlx(const fvMesh&, const fvMesh&);

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void fixEdgeConcentration( scalarField& );
  void fixWallDisplPeriodic( vectorField& );
  void fixCyclicNormal( vectorField&, scalarField& );
  
  //template<typename T> void averageCoupled( Field<T>& );
  template<typename T>
  void averageCoupled( Field<T>& field ){
    forAll(edgePoints, i){
      field[ edgePoints[i] ] *= 0.5;
    }
    Info<< "field size: " << field.size() << nl;
  }


  // fix edge displacement and calculate inlet displ.
  vectorField calculateInletDisplacement( vectorField& );

  vectorField wallRelaxation1();
  vectorField wallRelaxation(int);
  vectorField wallRelaxation(int, const vectorField&);

  vector vertexDisplZ(point, point, point);

  vector vertexDisplX(point, point, point);
  
  vector vertexDispl(point, point, point);
  
  void boundaryCheck();
  void boundaryFix( pointField& );
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // set up pairs for wall mesh relaxation: neighbZrlx, neighbXrlx
  void setUpPairsRlx();

  
  // set up pairs for concentration extrapollation on the edge
  void setUpPairsConc();

  // setup lists needed for mesh relaxation
  void setUpLists();

  // setup lists needed for mesh relaxation
  void checkFacePyramidsMy();
  
  
  // set/get methods
  float get_version();
};

#endif

