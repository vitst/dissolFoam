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

#include "syncTools.H"
#include "primitivePatchInterpolationSync.H"


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
  void doInletDisplacement(vectorField&);
  
  template<class Patch>
  vectorField wallRelaxation(const Patch& patch_){
    const labelList& meshPoints = patch_.meshPoints();

    const pointField& boundaryPoints = patch_.localPoints();
    const labelListList& plistFaces = patch_.pointFaces();
    const pointField& faceCs = patch_.faceCentres();

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

    forAll(displacement, i)  displacement[i] /= numberFaces[i];

    // getting correct point normal field
    vectorField fNorm = patch_.faceNormals();
    primitivePatchInterpolationSync patchInterpolator( patch_, mesh_ );
    vectorField pointNorm = patchInterpolator.faceToPointInterpolate(fNorm);

    vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);

    // fixing edges
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
  

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // set up pairs for concentration extrapollation on the edge
  void setUpPairsConc();

  // setup lists needed for mesh relaxation
  void setUpLists();

  
  // set/get methods
  float get_version();
};

#endif

