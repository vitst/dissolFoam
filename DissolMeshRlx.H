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

  // ID of each patch
  label wallID, inletID, outletID;

  // vetex ID maps: patch to global
  labelList wallsToAll, inletToAll, outletToAll;

  // list of local (in terms of wall) IDs for walls- inlet and outlet edges
  labelList local_wall_WallsInletEdges;
  labelList local_wall_WallsOutletEdges;


  
  labelListList inletTriple;
  // map for the extrapolation of the concentration on the edge
  // it fixes inlet boundary condition inconsistance
  std::map<int, std::pair<int, int> > edgeConcentrationFixMap;

  
  scalar extrapolateConcentrationExp(const pointField&, scalar&, scalar&, const labelList&);
  scalar extrapolateConcentrationLinear(const pointField&, const scalarField&, label, label, label);
  //scalar extrapolateConcentrationLinearZ(const pointField&, const scalarField&, label, label, label);
  scalar extrapolateConcentrationLinearZ(const pointField&, scalar&, scalar&, const labelList&);

public:
  
  // Constructors
  DissolMeshRlx(const fvMesh&);

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void fixEdgeConcentration( vectorField& );
  
  // fix edge displacement and calculate inlet displ.
  vectorField calculateInletDisplacement( vectorField& );

  vectorField wallRelaxation();

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // set up pairs for concentration extrapollation on the edge
  void setUpPairsConc();

  // setup lists needed for mesh relaxation
  void setUpLists();

  
  // set/get methods
  float get_version();
};

#endif

