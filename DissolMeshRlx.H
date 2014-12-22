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
  
  
  
  labelList local_inlet_WallsInletEdges;
  labelList scaleList;
  
  
  
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
  

/*  
  template<class Patch>
  vectorField inletRlx(const Patch& patch_, const vectorField& displ_){
    vectorField displ = displ_;
    vector nn(0,0,1);
    //displ = transform(I - nn*nn, displ);
    
    const labelList& meshPoints = patch_.meshPoints();

    const pointField& pCoord = patch_.localPoints();
    //const labelListList& plistFaces = patch_.pointFaces();
    const pointField& fCoord = patch_.faceCentres();
    
    vectorField pointDispl( pCoord.size() );
    vectorField faceDispl( fCoord.size() );
    
    primitivePatchInterpolationSync patchInterpolator( patch_, mesh_ );
    
    for(int iii=0; iii<10; iii++){
      
      forAll(global_WallInletEdges, i){
        label ind = findIndex( meshPoints, global_WallInletEdges[i]);
        if( ind!=-1 ){
          pointDispl[ind] = displ[ local_wall_WallsInletEdges[i] ];
        }
      }
      pointDispl = transform(I - nn*nn, pointDispl);
      
      faceDispl = patchInterpolator.pointToFaceInterpolate( pointDispl );
      //vectorField pointDispl1 = patchInterpolator.faceToPointInterpolate( faceDispl );
      pointDispl = patchInterpolator.faceToPointInterpolate( faceDispl );
      
    }
    return pointDispl;
  }
 */
  
  template<class Patch>
  vectorField inletRlx(const Patch& patch_, const vectorField& displ_){
    vectorField displ = displ_;
    vector nz(0,0,1);
    vector nx(1,0,0);
    //displ = transform(I - nn*nn, displ);
    
    const labelList& meshPoints = patch_.meshPoints();

    const pointField& pCoord = patch_.localPoints();
    //const labelListList& plistFaces = patch_.pointFaces();
    const pointField& fCoord = patch_.faceCentres();
    
    //vectorField pointDispl( pCoord.size() );
    vectorField pointDispl = pCoord;
    vectorField faceDispl( fCoord.size() );
    
    
    scalarField minScale( local_wall_WallsInletEdges.size() );
    
    forAll(global_WallInletEdges, i){
      label ind = findIndex( meshPoints, global_WallInletEdges[i]);
      if( ind!=-1 ){
        //pointDispl[ind] = displ[ local_wall_WallsInletEdges[i] ];
        //Pout<< pointDispl[ind].y()/pCoord[ind].y() <<nl;
        minScale[i] = displ[ local_wall_WallsInletEdges[i] ].y()/pCoord[ind].y();
      }
    }
    pointDispl = transform(I - nz*nz, pointDispl);
    pointDispl = transform(I - nx*nx, pointDispl);
    
    forAll(pointDispl, i){
      if( scaleList[i]!=-1 ){
        pointDispl[i].y() *= minScale[ scaleList[i] ];
      }
    }
    
    
    
    return pointDispl;
  }
  

  
  template<class Patch>
  vectorField wallRelaxation(const Patch& patch_, const scalarListList& weights){
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

      scalar ss =0;
      vector sumw(0,0,0);
      forAll(pFaces, j){
        label faceI = pFaces[j];
        point faceC = faceCs[faceI];
        //scalarList& pw = weights[i];

        vector d = faceC - curP;
        sumw += d * weights[i][j];
        ss += weights[i][j];
        
/*
        if( curP.x()>23.5 && curP.x()<24.5 && curP.y()>.4 ){
          Pout << curP << "   " << weights[i][j] << "   " << faceC << nl;
        }
 */
        
      }
      //Pout << curP << "   " << weights[i] << nl;
      // sum of all distances to face centres
      displacement[i] = sumw;
      //numberFaces[i] = ss;
    }

    syncTools::syncPointList( mesh_, meshPoints, displacement, plusEqOp<vector>(), vector::zero);
    syncTools::syncPointList( mesh_, meshPoints, numberFaces, plusEqOp<scalar>(), 0.0);

    //Pout << numberFaces << nl;
    //std::exit(0);
    
    //forAll(displacement, i)  displacement[i] /= numberFaces[i];
    
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
  
  template<class Patch>
  scalarListList calc_weights(const Patch& patch_){

    const labelList& meshPoints = patch_.meshPoints();

    const pointField& boundaryPoints = patch_.localPoints();
    const labelListList& plistFaces = patch_.pointFaces();
    const pointField& faceCs = patch_.faceCentres();

    scalarListList weights( boundaryPoints.size() );
    scalarField sumWeights( boundaryPoints.size() );

    forAll(boundaryPoints, i){
      point curP = boundaryPoints[i];
      const labelList& pFaces = plistFaces[i];

      scalarList& pw = weights[i];
      pw.setSize(pFaces.size());
      
      scalar sumw = 0.0;
      forAll(pFaces, j){
        label faceI = pFaces[j];
        point faceC = faceCs[faceI];

        vector d = faceC - curP;
        pw[j] = 1/mag(d);
        sumw += pw[j];
      }
      // sum of all distances to face centres
      sumWeights[i] = sumw;
    }

    syncTools::syncPointList( mesh_, meshPoints, sumWeights, plusEqOp<scalar>(), 0.0);

    forAll(weights, i){
      scalarList& pw = weights[i];
      forAll(pw, j){
        pw[j] /= sumWeights[i];
      }
    }

    return weights;
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

