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

  
  labelList global_WallInletEdges;
  labelList global_WallOutletEdges;
  
  
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

    vectorField inletPrj( global_WallInletEdges.size() );
    
    forAll(global_WallInletEdges, i){
      label ind = findIndex( meshPoints, global_WallInletEdges[i]);
      if( ind!=-1 ){
        inletPrj[i] = projectedDisplacement[ ind ];
      }
    }

    inletPrj = transform(I - ny*ny, inletPrj);
    inletPrj = transform(I - nz*nz, inletPrj);

    forAll(global_WallInletEdges, i){
      label ind = findIndex( meshPoints, global_WallInletEdges[i]);
      if( ind!=-1 ){
        projectedDisplacement[ ind ] = inletPrj[i];
      }
    }

    // project back outlet
    vectorField outletPrj( global_WallOutletEdges.size() );
    
    forAll(global_WallOutletEdges, i){
      label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
      if( ind!=-1 ){
        outletPrj[i] = projectedDisplacement[ ind ];
      }
    }

    outletPrj = transform(I - ny*ny, outletPrj);
    outletPrj = transform(I - nz*nz, outletPrj);

    forAll(global_WallOutletEdges, i){
      label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
      if( ind!=-1 ){
        projectedDisplacement[ ind ] = outletPrj[i];
      }
    }

    return projectedDisplacement;
  }
  
/*  
  template<class Patch>
  scalarField calc_weights(const Patch& patch_){

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

    
  }
  */

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // set up pairs for concentration extrapollation on the edge
  void setUpPairsConc();

  // setup lists needed for mesh relaxation
  void setUpLists();

  
  // set/get methods
  float get_version();
};

#endif

