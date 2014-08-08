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
  label wallID, inletID, outletID, cyclicID;

  // vetex ID maps: patch to global
  labelList wallsToAll, inletToAll, outletToAll, cyclicToAll;

  // oposite to previous lists. Map IDs global to patch
  std::map<int,int> allToInlet, allToOutlet, allToCyclic, allToWalls;

  // list of global IDs for walls- inlet and outlet edges
  std::list<int> glb_list_wallsInletEdges;
  std::list<int> glb_list_wallsOutletEdges;
  std::list<int> glb_list_wallsCyclicEdges;
  // list of local (in terms of wall) IDs for walls- inlet and outlet edges
  std::list<int> lcl_wall_list_wallsInletEdges;
  std::list<int> lcl_wall_list_wallsOutletEdges;
  std::list<int> lcl_wall_list_wallsCyclicEdges;


  labelListList pe;


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
  //std::map<int, int> wallCyclicEdgeToEdge;

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

  // fix edge displacement and calculate inlet displ.
  vectorField calculateInletDisplacement( vectorField& );

  vectorField wallRelaxation();

  vector vertexDisplZ(point, point, point);

  vector vertexDisplX(point, point, point);
  
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // set up pairs for wall mesh relaxation: neighbZrlx, neighbXrlx
  void setUpPairsRlx();

  // set up pairs for concentration extrapollation on the edge
  void setUpPairsConc();

  // setup lists needed for mesh relaxation
  void setUpLists();

  // set/get methods
  float get_version();
};

#endif

