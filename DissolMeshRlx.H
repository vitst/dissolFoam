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

// mesh search
#include "interpolation.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "meshSearch.H"

typedef List<vectorField> vectorFieldList;

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
  scalar extrapolateConcentrationExpZ(const pointField&, scalar&, scalar&, const labelList&);
  
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
  pointField doOutletDisplacement(const vectorField&);
  
  
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

      vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);
      
      // fixing edges
      vector ny(0,1,0); vector nz(0,0,1);

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
          
          vector disp;
          disp.x() = d.x() * 0.25;
          disp.y() = d.y() * 0.25;
          disp.z() = d.z() * pw[j];
          displacement[i] += disp;
          
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

      vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);
      //vectorField projectedDisplacement = displacement;
      //vectorField projectedDisplacement1 = transform(I - pointNorm*pointNorm, displacement);
      //vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, projectedDisplacement1);
      
      // fixing edges
      vector ny(0,1,0);
      vector nz(0,0,1);

      vectorField inletPrj( global_WallInletEdges.size() );
      vectorField inletPrjN( global_WallInletEdges.size() );
      forAll(global_WallInletEdges, i){
        label ind = findIndex( meshPoints, global_WallInletEdges[i]);
        if( ind!=-1 ){
          inletPrj[i] = projectedDisplacement[ ind ];
          inletPrjN[i] = pointNorm[ ind ];
        }
      }

      inletPrj = transform(I - ny*ny, inletPrj);
      //inletPrj = transform(I - inletPrjN*inletPrjN, inletPrj);
      inletPrj = transform(I - nz*nz, inletPrj);

      forAll(global_WallInletEdges, i){
        label ind = findIndex( meshPoints, global_WallInletEdges[i]);
        if( ind!=-1 ){
          projectedDisplacement[ ind ] = inletPrj[i];
        }
      }

      // project back outlet
      vectorField outletPrj( global_WallOutletEdges.size() );
      vectorField outletPrjN( global_WallOutletEdges.size() );

      forAll(global_WallOutletEdges, i){
        label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
        if( ind!=-1 ){
          outletPrj[i] = projectedDisplacement[ ind ];
          outletPrjN[i] = pointNorm[ ind ];
        }
      }

      outletPrj = transform(I - ny*ny, outletPrj);
      //outletPrj = transform(I - outletPrjN*outletPrjN, outletPrj);
      outletPrj = transform(I - nz*nz, outletPrj);

      forAll(global_WallOutletEdges, i){
        label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
        if( ind!=-1 ){
          projectedDisplacement[ ind ] = outletPrj[i];
        }
      }
      
      //scalar factor = 1;
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
    
    // project to the surface
    labelHashSet includePatches(1);
    includePatches.insert(wallID);
    triSurface wallTriSurface
    (
      triSurfaceTools::triangulate( mesh_.boundaryMesh(), includePatches )
    );
    const triSurfaceSearch querySurf(wallTriSurface);
    const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();
    bool ph = false;
    
    
    pointField faceCsFF = faceCentres(boundaryPoints, llf);
    vectorField faceNsFF = faceNormals(boundaryPoints, llf);
    vectorField pointNormFF( N, vector::zero );
    scalarList faceToPointSumWeightsFF( N, 0.0 );

    forAll(boundaryPoints, i){
      point& curP = boundaryPoints[i];
      const labelList& pFaces = plistFaces[i];

      forAll(pFaces, j){
        label faceI = pFaces[j];
        point& faceC = faceCsFF[faceI];
        vector d = faceC - curP;
        scalar mag_d = mag(d);
        // this is for normal
        scalar nw = 1.0 / mag_d;
        pointNormFF[i] += nw * faceNsFF[ faceI ];
        faceToPointSumWeightsFF[i] += nw;
      }
    }

    // synchronizing over the cyclic and processor boundaries
    syncTools::syncPointList(mesh_, meshPoints, pointNormFF, plusEqOp<vector>(), vector::zero);
    syncTools::syncPointList(mesh_, meshPoints, faceToPointSumWeightsFF, plusEqOp<scalar>(), 0.0);
    // normalization
    forAll(pointNormFF, pointi){
      pointNormFF[pointi] /= faceToPointSumWeightsFF[pointi];
    }
    
    

    forAll(boundaryPoints, i){
      point p = boundaryPoints[i];
      point searchStart = p + 5.0 * pointNormFF[i];
      //searchStart.y() = 0.0;
      
      /*
      if(searchStart.x()<=-25){
        //searchStart.z() += 0.01;
        searchStart.x() += 0.000001;
      }
      if(searchStart.x()>=25){
        //searchStart.z() += 0.01;
        searchStart.x() -= 0.000001;
      }
      */
      
      
      point searchEnd = p - 0.1 * pointNormFF[i];
      
      /*
      if(p.y()<0)
        searchStart.y() -= 10.0;
      else
        searchStart.y() += 10.0;
      */
        
      point hit = p;
      pointIndexHit pHit = tree.findLine(searchStart, searchEnd);
      if ( pHit.hit() )
      {
        hit =  pHit.hitPoint();
        ph = true;
      }
      else{
      }

      if(!ph || (hit.y() / std::abs(hit.y()) != p.y() / std::abs(p.y())) ){
        Info<<"WARNING! The ray did not find the surface for"
                << " boudaryPoint["<< i<< "]= "<< p 
                << " hitPoint= "<< hit 
                <<nl;
      }
      
      /*
      if( (p.x()<-24.9 || p.x()>24.9) && p.y()<0 && p.z()<4.5 && p.z()>4. ){
        Info << "Old point "
                << i << "  : "<< p
                << "  new point  : " << hit
                <<nl;
      }
      */
      
      boundaryPoints[i] = hit;
    }
    
    vectorField dd = (boundaryPoints - patch_.localPoints());
    scalarField cc(dd.size(), 1.0);
    syncTools::syncPointList(mesh_, meshPoints, dd, plusEqOp<vector>(), vector::zero);
    syncTools::syncPointList(mesh_, meshPoints, cc, plusEqOp<scalar>(), 0.0);
    
    forAll(dd, i){
      dd[i] /= cc[i];
    }
    
    return dd;
  }
  




  template<class Patch>
  vectorField wallRelaxation2(Patch& patch_, const vectorFieldList& weights, double rlxTol, const vectorField &vvff){
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
      
      
      forAll(global_WallInletEdges, i){
        label ind = findIndex( meshPoints, global_WallInletEdges[i]);
        label ind_wall = local_wall_WallsInletEdges[i];
        if( ind!=-1 ){
          displacement[ ind ] = vvff[ ind_wall ];
        }
      }
      
      
    
      forAll(boundaryPoints, i){
        point& curP = boundaryPoints[i];
        const labelList& pFaces = plistFaces[i];
        const vectorField& pw = weights[i];
        
        forAll(pFaces, j){
          label faceI = pFaces[j];
          point& faceC = faceCs[faceI];

          vector d = faceC - curP;
          
          scalar mag_d = mag(d);
          
          vector disp;
          disp.x() = d.x() * pw[j].x();
          disp.y() = d.y() * pw[j].y();
          disp.z() = d.z() * pw[j].z();
          displacement[i] += disp;
          
          tol[i] += mag( disp );
          
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

      vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);
      
      // fixing edges
      vector ny(0,1,0); vector nz(0,0,1);

      /*
      vectorField inletPrj( global_WallInletEdges.size() );
      forAll(global_WallInletEdges, i){
        label ind = findIndex( meshPoints, global_WallInletEdges[i]);
        if( ind!=-1 ){
          inletPrj[i] = projectedDisplacement[ ind ];
        }
      }
      inletPrj = transform(I - ny*ny, inletPrj);
      inletPrj = transform(I - nz*nz, inletPrj);
      */

      forAll(global_WallInletEdges, i){
        label ind = findIndex( meshPoints, global_WallInletEdges[i]);
        if( ind!=-1 ){
          //projectedDisplacement[ ind ] = inletPrj[i];
          projectedDisplacement[ ind ] = vector::zero;
        }
      }

      // project back outlet
      /*
      vectorField outletPrj( global_WallOutletEdges.size() );
      forAll(global_WallOutletEdges, i){
        label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
        if( ind!=-1 ){
          outletPrj[i] = projectedDisplacement[ ind ];
        }
      }
      outletPrj = transform(I - ny*ny, outletPrj);
      outletPrj = transform(I - nz*nz, outletPrj);
      */

      forAll(global_WallOutletEdges, i){
        label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
        if( ind!=-1 ){
          //projectedDisplacement[ ind ] = outletPrj[i];
          projectedDisplacement[ ind ] = vector::zero;
        }
      }
      
      scalar factor = 1.0;
      //scalar factor = 6.0;
      //if(itt%2==0) factor = 1.17;
      vectorField finalDisplacement = factor * projectedDisplacement;
      
      boundaryPoints += finalDisplacement;
      
      displ_tol = average( mag(finalDisplacement)/tol );
      reduce(displ_tol, sumOp<scalar>());
      displ_tol = 100000 * displ_tol / static_cast<double>( Pstream::nProcs() );
      
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
  vectorField wallRelaxation3(Patch& patch_, const vectorFieldList& weights, double rlxTol, const vectorField &vvff){
    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();

    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    double displ_tol = 1.0;
    int itt = 0;
    
    vectorField displacement(N, vector::zero);
    forAll(global_WallInletEdges, i){
      label ind = findIndex( meshPoints, global_WallInletEdges[i]);
      label ind_wall = local_wall_WallsInletEdges[i];
      if( ind!=-1 ){
        displacement[ ind ] = vvff[ ind_wall ];
      }
    }

    
    while(displ_tol>rlxTol){
      
      primitivePatchInterpolationSync patchInterpolator( patch_, mesh_ );
      vectorField faceDispl = patchInterpolator.pointToFaceInterpolate(displacement);
    
      vectorField newDisplacement = patchInterpolator.faceToPointInterpolate(faceDispl);
      forAll(global_WallInletEdges, i){
        label ind = findIndex( meshPoints, global_WallInletEdges[i]);
        label ind_wall = local_wall_WallsInletEdges[i];
        if( ind!=-1 ){
          newDisplacement[ ind ] = vvff[ ind_wall ];
        }
      }
      
      displ_tol = average( mag(newDisplacement - displacement) );
      scalar displ_tol_norm = average( mag(newDisplacement) );
      reduce(displ_tol, sumOp<scalar>());
      reduce(displ_tol_norm, sumOp<scalar>());
      displ_tol /= displ_tol_norm;
      displ_tol *= 10000;
      //displ_tol = displ_tol / static_cast<double>( Pstream::nProcs() );
      
      displacement = newDisplacement;
      
      if(itt%100==0){
        Info<<patch_.name()<<"  rlx iter "<<itt
                <<"  tolerance: "<<displ_tol<< nl;
      }
      
      itt+=1;
    }
    
    forAll(global_WallInletEdges, i){
      label ind = findIndex( meshPoints, global_WallInletEdges[i]);
      if( ind!=-1 ){
        displacement[ ind ] = vector::zero;
      }
    }
      
    Info<< patch_.name() << "  converged in "<< itt 
            << " iterations. Tolerance: " << displ_tol<< nl;
    
    return displacement;// (boundaryPoints - patch_.localPoints());
  }
  

  template<class Patch>
  vectorField wallRelaxation31(Patch& patch_, const vectorFieldList& weights, double rlxTol, const vectorField &vvff){
    pointField saveBoundaryPoints = patch_.localPoints();

    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    double displ_tol = 1.0;
    int itt = 0;
    
    
    vectorField displacement(N, vector::zero);
    forAll(global_WallOutletEdges, i){
      label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
      label ind_wall = local_wall_WallsOutletEdges[i];
      if( ind!=-1 ){
        displacement[ ind ] = vvff[ ind_wall ];
      }
    }

    
    while(displ_tol>rlxTol){
      
      // calculate current face centers
      pointField faceCs = faceCentres(boundaryPoints, llf);
      vectorField faceNs = faceNormals(boundaryPoints, llf);
      vectorField pointNorm( N, vector::zero );
      scalarList faceToPointSumWeights( N, 0.0 );

      // create a displacement field for points
      
      primitivePatchInterpolationSync patchInterpolator( patch_, mesh_ );
      vectorField faceDispl = patchInterpolator.pointToFaceInterpolate(displacement);
    
      vectorField newDisplacement = patchInterpolator.faceToPointInterpolate(faceDispl);
      forAll(global_WallOutletEdges, i){
        label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
        label ind_wall = local_wall_WallsOutletEdges[i];
        if( ind!=-1 ){
          newDisplacement[ ind ] = vvff[ ind_wall ];
        }
      }
      
      displ_tol = average( mag(newDisplacement - displacement) );
      scalar displ_tol_norm = average( mag(newDisplacement) );
      reduce(displ_tol, sumOp<scalar>());
      reduce(displ_tol_norm, sumOp<scalar>());
      displ_tol /= displ_tol_norm;
      displ_tol *= 10000;
      //displ_tol = displ_tol / static_cast<double>( Pstream::nProcs() );
      
      displacement = newDisplacement;
      
      if(itt%100==0){
        Info<<patch_.name()<<"  rlx iter "<<itt
                <<"  tolerance: "<<displ_tol<< nl;
      }
      
      itt+=1;
    }
    
    forAll(global_WallOutletEdges, i){
      label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
      if( ind!=-1 ){
        displacement[ ind ] = vector::zero;
      }
    }
      
    Info<< patch_.name() << "  converged in "<< itt 
            << " iterations. Tolerance: " << displ_tol<< nl;
    
    return displacement;// (boundaryPoints - patch_.localPoints());
  }
  
  
  

  
  template<class Patch>
  vectorField wallRelaxation4(Patch& patch_, const vectorFieldList& weights, double rlxTol, const vectorField &vvff){
    pointField saveBoundaryPoints = patch_.localPoints();

    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    double displ_tol = 1.0;
    int itt = 0;
    
    vectorField displacement(N, vector::zero);
    forAll(global_WallInletEdges, i){
      label ind = findIndex( meshPoints, global_WallInletEdges[i]);
      label ind_wall = local_wall_WallsInletEdges[i];
      if( ind!=-1 ){
        displacement[ ind ] = vvff[ ind_wall ];
      }
    }
    
    while(displ_tol>rlxTol){
      // create a displacement field for points
      vectorField faceDispl( llf.size(), vector::zero );
      forAll(faceDispl, facei)
      {
        const labelList& curPoints = llf[facei];
        forAll(curPoints, pointi)
        {
          faceDispl[facei] += displacement[curPoints[pointi]];
        }
        faceDispl[facei] /= curPoints.size();
      }
      
      vectorField newDisplacement(N, vector::zero);
      forAll(newDisplacement, i){
        const labelList& pFaces = plistFaces[i];
        const vectorField& pw = weights[i];
        
        forAll(pFaces, j){
          label faceI = pFaces[j];

          vector d = faceDispl[faceI];
          
          vector disp;
          disp.x() = d.x() * pw[j].x();
          disp.y() = d.y() * pw[j].y();
          //disp.z() = d.z() * pw[j].z();
          disp.z() = d.z() * 0.0;
          newDisplacement[i] += disp;
        }
      }
      
      syncTools::syncPointList(mesh_, meshPoints, newDisplacement, plusEqOp<vector>(), vector::zero);
      
      forAll(global_WallInletEdges, i){
        label ind = findIndex( meshPoints, global_WallInletEdges[i]);
        label ind_wall = local_wall_WallsInletEdges[i];
        if( ind!=-1 ){
          newDisplacement[ ind ] = vvff[ ind_wall ];
        }
      }
      
      displ_tol = average( mag(newDisplacement - displacement) );
      scalar displ_tol_norm = average( mag(newDisplacement) );
      reduce(displ_tol, sumOp<scalar>());
      reduce(displ_tol_norm, sumOp<scalar>());
      displ_tol /= displ_tol_norm;
      //displ_tol *= 10000;
      //displ_tol = displ_tol / static_cast<double>( Pstream::nProcs() );
      
      displacement = newDisplacement;
      
      if(itt%1000==0){
        Info<<patch_.name()<<"  rlx iter "<<itt
                <<"  tolerance: "<<displ_tol<< nl;
      }
      itt+=1;
    }
    
    forAll(global_WallInletEdges, i){
      label ind = findIndex( meshPoints, global_WallInletEdges[i]);
      if( ind!=-1 ){
        displacement[ ind ] = vector::zero;
      }
    }
      
    Info<< patch_.name() << "  converged in "<< itt 
            << " iterations. Tolerance: " << displ_tol<< nl;
    
    return displacement;// (boundaryPoints - patch_.localPoints());
  }
  

  
  
  
  template<class Patch>
  vectorField wallRelaxation41(Patch& patch_, const scalarListList& weights, double rlxTol, const vectorField &vvff){
    pointField saveBoundaryPoints = patch_.localPoints();

    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    const labelListList& ppp = patch_.pointEdges();
    
    const edgeList ee = patch_.edges();
    
    /*
    forAll(ppp, ii){
      if(ppp[ii].size()==2)
        Info << ii << "  size:  " << ppp[ii] <<"  "<< ee[ppp[ii][0]]<<"  "<< ee[ppp[ii][1]] << nl;
    }
    std::exit(0);
    */
    
    
    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    double displ_tol = 1.0;
    int itt = 0;
    
    vectorField displacement(N, vector::zero);
    forAll(global_WallInletEdges, i){
      label ind = findIndex( meshPoints, global_WallInletEdges[i]);
      label ind_wall = local_wall_WallsInletEdges[i];
      if( ind!=-1 ){
        displacement[ ind ] = vvff[ ind_wall ];
      }
    }
    
    while(displ_tol>rlxTol){
      // create a displacement field for points
      vectorField faceDispl( llf.size(), vector::zero );
      forAll(faceDispl, facei)
      {
        const labelList& curPoints = llf[facei];
        forAll(curPoints, pointi)
        {
          faceDispl[facei] += displacement[curPoints[pointi]];
        }
        faceDispl[facei] /= curPoints.size();
      }
      
      vectorField newDisplacement(N, vector::zero);
      forAll(newDisplacement, i){
        const labelList& pFaces = plistFaces[i];
        const scalarList& pw = weights[i];
        
        forAll(pFaces, j){
          label faceI = pFaces[j];

          vector d = faceDispl[faceI];
          
          vector disp;
          
          disp.x() = d.x() * 0.25;
          disp.y() = d.y() * pw[j];
          //disp.z() = d.z() * pw[j].z();
          disp.z() = d.z() * 0.0;
          
          newDisplacement[i] += disp;
        }
      }
      
      syncTools::syncPointList(mesh_, meshPoints, newDisplacement, plusEqOp<vector>(), vector::zero);
      
      forAll(global_WallInletEdges, i){
        label ind = findIndex( meshPoints, global_WallInletEdges[i]);
        label ind_wall = local_wall_WallsInletEdges[i];
        if( ind!=-1 ){
          newDisplacement[ ind ] = vvff[ ind_wall ];
        }
      }
      
      displ_tol = average( mag(newDisplacement - displacement) );
      scalar displ_tol_norm = average( mag(newDisplacement) );
      reduce(displ_tol, sumOp<scalar>());
      reduce(displ_tol_norm, sumOp<scalar>());
      displ_tol /= displ_tol_norm;
      //displ_tol *= 10000;
      //displ_tol = displ_tol / static_cast<double>( Pstream::nProcs() );
      
      displacement = newDisplacement;
      
      if(itt%1000==0){
        Info<<patch_.name()<<"  rlx iter "<<itt
                <<"  tolerance: "<<displ_tol<< nl;
      }
      itt+=1;
    }
    
    forAll(global_WallInletEdges, i){
      label ind = findIndex( meshPoints, global_WallInletEdges[i]);
      if( ind!=-1 ){
        displacement[ ind ] = vector::zero;
      }
    }
      
    Info<< patch_.name() << "  converged in "<< itt 
            << " iterations. Tolerance: " << displ_tol<< nl;
    
    return displacement;
  }
  

  
  
  template<class Patch>
  vectorField wallRelaxation42(Patch& patch_, const scalarListList& weights, double rlxTol, const vectorField &vvff){
    pointField saveBoundaryPoints = patch_.localPoints();

    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    const labelListList& ppp = patch_.pointEdges();
    const edgeList& ee = patch_.edges();


    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    double displ_tol = 1.0;
    int itt = 0;
    
    vectorField displacement(N, vector::zero);
    forAll(global_WallInletEdges, i){
      label ind = findIndex( meshPoints, global_WallInletEdges[i]);
      label ind_wall = local_wall_WallsInletEdges[i];
      if( ind!=-1 ){
        displacement[ ind ] = vvff[ ind_wall ];
      }
    }
    
    while(displ_tol>rlxTol){
      
      scalarList sumWeights( N, 0.0 );
      vectorField newDisplacement(N, vector::zero);
      forAll(newDisplacement, i){
        point& curP = boundaryPoints[i];
        const labelList& pFaces = plistFaces[i];
        const scalarList& pw = weights[i];
        
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
          vector d = displacement[pI];
          point& np = boundaryPoints[pI];
          
          vector wv = np - curP;
          scalar mag_wv = mag(wv);
          scalar nw = 1.0 / mag_wv;
          
          
          vector disp = nw * d;
          
          disp.z() = 0.0;
          
          newDisplacement[i] += disp;
          
          sumWeights[i] += nw;
        }
      }
      
      syncTools::syncPointList(mesh_, meshPoints, newDisplacement, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, meshPoints, sumWeights, plusEqOp<scalar>(), 0.0);
      // normalization
      forAll(newDisplacement, pointi){
        newDisplacement[pointi] /= sumWeights[pointi];
      }
      
      forAll(global_WallInletEdges, i){
        label ind = findIndex( meshPoints, global_WallInletEdges[i]);
        label ind_wall = local_wall_WallsInletEdges[i];
        if( ind!=-1 ){
          newDisplacement[ ind ] = vvff[ ind_wall ];
        }
      }
      
      displ_tol = average( mag(newDisplacement - displacement) );
      scalar displ_tol_norm = average( mag(newDisplacement) );
      reduce(displ_tol, sumOp<scalar>());
      reduce(displ_tol_norm, sumOp<scalar>());
      displ_tol /= displ_tol_norm;
      //displ_tol *= 10000;
      //displ_tol = displ_tol / static_cast<double>( Pstream::nProcs() );
      
      displacement = newDisplacement;
      
      if(itt%1000==0){
        Info<<patch_.name()<<"  rlx iter "<<itt
                <<"  tolerance: "<<displ_tol<< nl;
      }
      itt+=1;
    }
    
    forAll(global_WallInletEdges, i){
      label ind = findIndex( meshPoints, global_WallInletEdges[i]);
      if( ind!=-1 ){
        displacement[ ind ] = vector::zero;
      }
    }
      
    Info<< patch_.name() << "  converged in "<< itt 
            << " iterations. Tolerance: " << displ_tol<< nl;
    
    return displacement;
  }
  
  template<class Patch>
  vectorField wallRelaxation43(Patch& patch_, const scalarListList& weights, double rlxTol, const vectorField &vvff){
    pointField saveBoundaryPoints = patch_.localPoints();

    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    const labelListList& ppp = patch_.pointEdges();
    const edgeList& ee = patch_.edges();
    
    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    double displ_tol = 1.0;
    int itt = 0;
    
    vectorField displacement(N, vector::zero);
    forAll(global_WallOutletEdges, i){
      label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
      label ind_wall = local_wall_WallsOutletEdges[i];
      if( ind!=-1 ){
        displacement[ ind ] = vvff[ ind_wall ];
      }
    }
    
    while(displ_tol>rlxTol){
      
      scalarList sumWeights( N, 0.0 );
      vectorField newDisplacement(N, vector::zero);
      forAll(newDisplacement, i){
        point& curP = boundaryPoints[i];
        const labelList& pFaces = plistFaces[i];
        const scalarList& pw = weights[i];
        
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
          vector d = displacement[pI];
          point& np = boundaryPoints[pI];
          
          vector wv = np - curP;
          scalar mag_wv = mag(wv);
          scalar nw = 1.0 / mag_wv;
          
          vector disp = nw * d;
          
          disp.z() = 0.0;
          
          newDisplacement[i] += disp;
          
          sumWeights[i] += nw;
        }
      }
      
      syncTools::syncPointList(mesh_, meshPoints, newDisplacement, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, meshPoints, sumWeights, plusEqOp<scalar>(), 0.0);
      // normalization
      forAll(newDisplacement, pointi){
        newDisplacement[pointi] /= sumWeights[pointi];
      }
      
      forAll(global_WallOutletEdges, i){
        label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
        label ind_wall = local_wall_WallsOutletEdges[i];
        if( ind!=-1 ){
          newDisplacement[ ind ] = vvff[ ind_wall ];
        }
      }
      
      displ_tol = average( mag(newDisplacement - displacement) );
      scalar displ_tol_norm = average( mag(newDisplacement) );
      reduce(displ_tol, sumOp<scalar>());
      reduce(displ_tol_norm, sumOp<scalar>());
      displ_tol /= displ_tol_norm;
      
      displacement = newDisplacement;
      
      if(itt%1000==0){
        Info<<patch_.name()<<"  rlx iter "<<itt
                <<"  tolerance: "<<displ_tol<< nl;
      }
      itt+=1;
    }
    
    forAll(global_WallOutletEdges, i){
      label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
      if( ind!=-1 ){
        displacement[ ind ] = vector::zero;
      }
    }
      
    Info<< patch_.name() << "  converged in "<< itt 
            << " iterations. Tolerance: " << displ_tol<< nl;
    
    return displacement;
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
  
  

  template<class Patch>
  scalarListList calc_weights1(const Patch& patch_, scalar lambda, int dir){

    const labelList& meshPoints = patch_.meshPoints();

    const pointField& boundaryPoints = patch_.localPoints();
    const labelListList& plistFaces = patch_.pointFaces();
    const pointField& faceCs = patch_.faceCentres();
    
    scalarListList weights( boundaryPoints.size() );
    scalarList sumWeights( boundaryPoints.size() );

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
        
        pw[j] = 1/(mag(d) * mag(d) * mag(d) * mag(d) * mag(d) * mag(d));
        
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
  
  

  template<class Patch>
  vectorFieldList calc_weights2(const Patch& patch_, scalar lambda, int dir){

    const labelList& meshPoints = patch_.meshPoints();

    const pointField& boundaryPoints = patch_.localPoints();
    const labelListList& plistFaces = patch_.pointFaces();
    const pointField& faceCs = patch_.faceCentres();
    
    vectorFieldList weights( boundaryPoints.size() );
    vectorField sumWeights( boundaryPoints.size() );

    forAll(boundaryPoints, i){
      point curP = boundaryPoints[i];
      const labelList& pFaces = plistFaces[i];

      vectorField& pw = weights[i];
      pw.setSize(pFaces.size());
      
      vector sumw = vector::zero;
      forAll(pFaces, j){
        label faceI = pFaces[j];
        point faceC = faceCs[faceI];

        vector d = faceC - curP;
        
        vector w = vector::zero;
        
        for(int ii=0; ii<3; ii++){
          if(d[ii]!=0)
            w[ii] = 1/ std::abs( d[ii] );
          else
            w[ii] = 0.0;
        }
        
        //pw[j] = (1 - 0.5 * lambda * d.x()/std::abs(d.x()) * faceC.x()/std::abs(faceC.x()) );
        pw[j] = w;
        
        sumw += pw[j];
      }
      // sum of all distances to face centres
      sumWeights[i] = sumw;
    }

    syncTools::syncPointList( mesh_, meshPoints, sumWeights, plusEqOp<vector>(), vector::zero);
    
    forAll(weights, i){
      vectorField& pw = weights[i];
      vector &sw = sumWeights[i];
      //Info<< "sw  " << sw << "     weights: "<< pw << nl;
      forAll(pw, j){
        vector &ppw = pw[j];
        for(int ii=0; ii<3; ii++){
          if( sw[ii]!=0 ){
            ppw[ii] /= sw[ii];
          }
        }
        
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

