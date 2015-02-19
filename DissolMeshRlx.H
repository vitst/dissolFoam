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
#include "dynamicFvMesh.H"

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
  
  labelList local_outlet_WallsOutletEdges;
  labelList scaleListOutlet;
  
  
  labelListList inletTriple;
  // map for the extrapolation of the concentration on the edge
  // it fixes inlet boundary condition inconsistance
  std::map<int, std::pair<int, int> > edgeConcentrationFixMap;

  
  scalar extrapolateConcentrationExp(const pointField&, scalar&, scalar&, const labelList&);
  scalar extrapolateConcentrationLinear(const pointField&, const scalarField&, label, label, label);
  //scalar extrapolateConcentrationLinearZ(const pointField&, const scalarField&, label, label, label);
  scalar extrapolateConcentrationLinearZ(const pointField&, scalar&, scalar&, const labelList&);

  Foam::pointField faceCentres(const pointField& points, const List<face>& flist) const;
  Foam::vectorField faceNormals(const pointField& points, const List<face>& flist) const;

  Foam::vectorField localFaceToPointNormalInterpolate(
          const pointField& points,
          const pointField& faceCs,
          const vectorField& faceNs,
          const labelListList& pointFaces,
          const labelList& meshPoints

  ) const;
  
public:
  
  // Constructors
  DissolMeshRlx(const fvMesh&);

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void fixEdgeConcentration( vectorField& );
  
  // fix edge displacement and calculate inlet displ.
  vectorField calculateInletDisplacement( vectorField& );
  
  vectorField calculateOutletDisplacement( vectorField& );
  
  pointField doInletDisplacement(const vectorField&);
  pointField doWallDisplacement(const vectorField&);
  
  
  template<class Patch>
  vectorField wallRelaxation(Patch& patch_, const scalarListList& weights, double rlxTol){
    pointField saveBoundaryPoints = patch_.localPoints();

    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    double displ_tol = 1.0;
    int itt = 0;
    while(displ_tol>rlxTol){
      
      // calculate current face centers
      pointField faceCs = faceCentres(boundaryPoints, llf);
      vectorField faceNs = faceNormals(boundaryPoints, llf);
      vectorField pointNorm( N, vector::zero );
      scalarList faceToPointSumWeights( N, 0.0 );

      // create a displacement field for points
      vectorField displacement(N, vector::zero);
      scalarField tol(N, 0.0);
    
      forAll(boundaryPoints, i){
        point& curP = boundaryPoints[i];
        const labelList& pFaces = plistFaces[i];
        const scalarList& pw = weights[i];

        forAll(pFaces, j){
          label faceI = pFaces[j];
          point& faceC = faceCs[faceI];

          vector d = faceC - curP;
          
          scalar mag_d = mag(d);
          
          displacement[i] += d * pw[j];
          tol[i] += mag_d * pw[j];
          
          // this is for normal
          scalar nw = 1.0 / mag_d;
          pointNorm[i] += nw * faceNs[ faceI ];
          faceToPointSumWeights[i] += nw;
        }
      }

      // synchronizing over the cyclic and processor boundaries
      syncTools::syncPointList(mesh_, meshPoints, displacement, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, meshPoints, tol, plusEqOp<scalar>(), 0.0);
      
      syncTools::syncPointList(mesh_, meshPoints, pointNorm, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, meshPoints, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
      // normalization
      forAll(pointNorm, pointi){
        pointNorm[pointi] /= faceToPointSumWeights[pointi];
      }

      // getting correct point normal field
      //vectorField faceNs = faceNormals(boundaryPoints, llf);
      //vectorField faceNs = patch_.faceNormals();
      
      //vectorField pointNorm = localFaceToPointNormalInterpolate(boundaryPoints, faceCs, faceNs, plistFaces, meshPoints);

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
      
      scalar factor = 6.0;
      if(itt%2==0) factor = 1.17;
      vectorField finalDisplacement = factor * projectedDisplacement;
      
      boundaryPoints += finalDisplacement;
      
      displ_tol = average( mag(finalDisplacement)/tol );
      reduce(displ_tol, sumOp<scalar>());
      displ_tol = displ_tol / static_cast<double>( Pstream::nProcs() );
      
      if(itt%1000==0){
        Info<<patch_.name()<<"  rlx iter "<<itt
                <<"  factor "<<factor<<"  tolerance: "<<displ_tol<< nl;
      }
      
      itt+=1;
    }
    Info<< patch_.name() << "  converged in "<< itt 
            << " iterations. Tolerance: " << displ_tol<< nl;
    
    return (boundaryPoints - patch_.localPoints());
  }






  template<class Patch>
  vectorField wallRelaxation1(Patch& patch_, const scalarListList& weights, double rlxTol){
    pointField saveBoundaryPoints = patch_.localPoints();

    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    double displ_tol = 1.0;
    int itt = 0;
    while(displ_tol>rlxTol){
      
      // calculate current face centers
      pointField faceCs = faceCentres(boundaryPoints, llf);
      vectorField faceNs = faceNormals(boundaryPoints, llf);
      vectorField pointNorm( N, vector::zero );
      scalarList faceToPointSumWeights( N, 0.0 );

      // create a displacement field for points
      vectorField displacement(N, vector::zero);
      scalarField tol(N, 0.0);
    
      forAll(boundaryPoints, i){
        point& curP = boundaryPoints[i];
        const labelList& pFaces = plistFaces[i];
        const scalarList& pw = weights[i];

        forAll(pFaces, j){
          label faceI = pFaces[j];
          point& faceC = faceCs[faceI];

          vector d = faceC - curP;
          
          scalar mag_d = mag(d);
          
          displacement[i] += d * pw[j];
          tol[i] += mag_d * pw[j];
          
          // this is for normal
          scalar nw = 1.0 / mag_d;
          pointNorm[i] += nw * faceNs[ faceI ];
          faceToPointSumWeights[i] += nw;
        }
      }

      // synchronizing over the cyclic and processor boundaries
      syncTools::syncPointList(mesh_, meshPoints, displacement, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, meshPoints, tol, plusEqOp<scalar>(), 0.0);
      
      syncTools::syncPointList(mesh_, meshPoints, pointNorm, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, meshPoints, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
      // normalization
      forAll(pointNorm, pointi){
        pointNorm[pointi] /= faceToPointSumWeights[pointi];
      }

      // getting correct point normal field
      //vectorField faceNs = faceNormals(boundaryPoints, llf);
      //vectorField faceNs = patch_.faceNormals();
      
      //vectorField pointNorm = localFaceToPointNormalInterpolate(boundaryPoints, faceCs, faceNs, plistFaces, meshPoints);

      //vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);
      
      // fixing edges
      vector ny(0,1,0);
      vector nz(0,0,1);
      vectorField projectedDisplacement = transform(I - nz*nz, displacement);

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
      
      scalar factor = 6.0;
      if(itt%2==0) factor = 1.17;
      vectorField finalDisplacement = factor * projectedDisplacement;
      
      boundaryPoints += finalDisplacement;
      
      displ_tol = average( mag(finalDisplacement)/tol );
      reduce(displ_tol, sumOp<scalar>());
      displ_tol = displ_tol / static_cast<double>( Pstream::nProcs() );
      
      if(itt%1000==0){
        Info<<patch_.name()<<"  rlx iter "<<itt
                <<"  factor "<<factor<<"  tolerance: "<<displ_tol<< nl;
      }
      
      itt+=1;
    }
    Info<< patch_.name() << "  converged in "<< itt 
            << " iterations. Tolerance: " << displ_tol<< nl;
    
    return (boundaryPoints - patch_.localPoints());
  }










  
  
  
  
  
  
  
  template<class Patch>
  scalarListList calc_weights(const Patch& patch_, scalar lambda, int dir){

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
        switch (dir){
          case 1:
          {
            if( std::abs(d.x()) < SMALL  ){
              Info<< "SMALLLLLLLLLLLLLL:   " << curP << "  " << faceC << endl;
              pw[j] = 1.0;
            }
            else{
              pw[j] = (1 - 0.5 * lambda * d.x()/std::abs(d.x()) * faceC.x()/std::abs(faceC.x()) );
            }
            break;
          }
          case 2:
          {
            if( std::abs(d.y()) < SMALL  ){
              Info<< "SMALLLLLLLLLLLLLL:   " << curP << "  " << faceC << endl;
              pw[j] = 1.0;
            }
            else{
              pw[j] = (1 - 0.5 * lambda * d.y()/std::abs(d.y()) * faceC.y()/std::abs(faceC.y()) );
            }
            break;
          }
          case 3:
          {
            if( std::abs(d.z()) < SMALL  ){
              Info<< "SMALLLLLLLLLLLLLL:   " << curP << "  " << faceC << endl;
              pw[j] = 1.0;
            }
            else{
              pw[j] = (1 - 0.5 * lambda * d.z()/std::abs(d.z()) * faceC.z()/std::abs(faceC.z()) );
            }
            break;
          }
          default:
          // TODO Error
            ;
        }
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

