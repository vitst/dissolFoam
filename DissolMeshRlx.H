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
#include "coupledPatchInterpolation.H"
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
  
  scalar extrapolateConcentrationLinear(const pointField&, scalar&, scalar&, const labelList&);
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
  
  vector vertexDispl( point A0, point A1, point A2 );
  vector vertexDispl1( vector r1, vector r2 );
  
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
      
      scalarField aux_f = mag(finalDisplacement)/tol;
      scalarField aux_f0(finalDisplacement.size(), 0.0);
      
      if( aux_f == aux_f0 ){
        displ_tol = 0.0;
      }
      else{
        displ_tol = average( aux_f );
      }
      
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
      //vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, projectedDisplacement1);
      
      // fixing edges
      //vector ny(0,1,0);
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

      //inletPrj = transform(I - ny*ny, inletPrj);
      inletPrj = transform(I - inletPrjN*inletPrjN, inletPrj);
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

      //outletPrj = transform(I - ny*ny, outletPrj);
      outletPrj = transform(I - outletPrjN*outletPrjN, outletPrj);
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
      

      scalarField aux_f = mag(finalDisplacement)/tol;
      scalarField aux_f0(finalDisplacement.size(), 0.0);
      
      if( aux_f == aux_f0 ){
        displ_tol = 0.0;
      }
      else{
        displ_tol = average( aux_f );
      }

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
  vectorField wallRelaxationInletEdge(Patch& patch_, double rlxTol){
    pointField saveBoundaryPoints = patch_.localPoints();

    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    const labelListList& ppp = patch_.pointEdges();
    const edgeList& ee = patch_.edges();
    
    vectorField globalDisplacement(N, vector::zero);
    
    // relax edges
    // inlet
    label NN = local_wall_WallsInletEdges.size();
    labelList ne(NN, -1);
    labelListList nepe;
    forAll(local_wall_WallsInletEdges, i){
      label ind = local_wall_WallsInletEdges[i];

      const labelList& lll = ppp[ind];
      labelList nel;
      
      
      forAll(lll, ii){
        const edge& ed = ee[ lll[ii] ];
        label edb = ed.start();
        
        if( edb!=ind && findIndex(local_wall_WallsInletEdges, edb) != -1 ){
          nel.append(edb);
        }
        
        label ede = ed.end();
        if( ede!=ind && findIndex(local_wall_WallsInletEdges, ede) != -1){
          nel.append(ede);
        }
      }
      nepe.append(nel);
    }
    
    Info<<boundaryPoints.size()<<"    "<<nepe<<nl<<nl;
    
    
    double displ_tol = 1.0;
    int itt = 0;
    while(displ_tol>rlxTol){
      pointField maxNeib(NN, vector::zero);
      vectorField minNeib(NN, vector::zero);
      forAll(nepe, i){
        label  curI = local_wall_WallsInletEdges[i];
        point& curP = boundaryPoints[curI];
        const labelList& pNeib = nepe[i];
        point& neibP0 = boundaryPoints[pNeib[0]];
        vector d1 = curP - neibP0;
        maxNeib[i] = d1;
        minNeib[i] = d1;
        forAll(pNeib, ii){
          point& neibP = boundaryPoints[pNeib[ii]];
          vector d2 = curP - neibP;

          if( mag(maxNeib[i])>mag(d2) ){
            maxNeib[i] = d2;
          }
          if( mag(minNeib[i])<=mag(d2) ){
            minNeib[i] = d2;
          }
        }
      }

      syncTools::syncPointList(mesh_, global_WallInletEdges, maxNeib, maxEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, global_WallInletEdges, minNeib, minEqOp<vector>(), vector::zero);


      vectorField displacement(NN, vector::zero);
      scalarField tol(NN, 0.0);
      forAll(local_wall_WallsInletEdges, i){
        label ind = local_wall_WallsInletEdges[i];
        displacement[i] = vertexDispl1(maxNeib[i], minNeib[i]);
        boundaryPoints[ind] -= displacement[i];
        tol[i] = mag(maxNeib[i]) + mag(minNeib[i]);
      }
      
      scalarField aux_f = mag(displacement)/tol;
      scalarField aux_f0(NN, 0.0);

      if( aux_f == aux_f0 ){
        displ_tol = 0.0;
      }
      else{
        displ_tol = average( aux_f );
      }

      reduce(displ_tol, sumOp<scalar>());
      displ_tol = displ_tol / static_cast<double>( Pstream::nProcs() );
      
      if(itt%1==0){
        Info<<patch_.name()<<"  rlx iter "<<itt
                <<"  tolerance: "<<displ_tol<< nl;
      }
      
      itt+=1;
      
      //if(itt==10) break;
      
    }
    
    forAll(globalDisplacement, i){
      globalDisplacement[i] = patch_.localPoints()[i]-boundaryPoints[i];
    }
    
    return globalDisplacement;
  }

  
  template<class Patch>
  vectorField wallRelaxationInletEdge1(Patch& patch_, double rlxTol){
    pointField saveBoundaryPoints = patch_.localPoints();

    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    const labelListList& ppp = patch_.pointEdges();
    const edgeList& ee = patch_.edges();
    
    vectorField globalDisplacement(N, vector::zero);
    
    // relax edges
    // inlet
    label NN = local_wall_WallsInletEdges.size();
    labelList ne(NN, -1);
    labelListList nepe;
    forAll(local_wall_WallsInletEdges, i){
      label ind = local_wall_WallsInletEdges[i];

      const labelList& lll = ppp[ind];
      labelList nel;
      
      
      forAll(lll, ii){
        const edge& ed = ee[ lll[ii] ];
        label edb = ed.start();
        
        if( edb!=ind && findIndex(local_wall_WallsInletEdges, edb) != -1 ){
          nel.append(edb);
        }
        
        label ede = ed.end();
        if( ede!=ind && findIndex(local_wall_WallsInletEdges, ede) != -1){
          nel.append(ede);
        }
      }
      nepe.append(nel);
    }
    
    //Info<<boundaryPoints.size()<<"    "<<nepe<<nl<<nl;
    
    
    double displ_tol = 1.0;
    int itt = 0;
    while(displ_tol>rlxTol){
      
      pointField faceCs = faceCentres(boundaryPoints, llf);
      vectorField faceNs = faceNormals(boundaryPoints, llf);
      
      vectorField displacement(NN, vector::zero);
      scalarField tol(NN, 0.0);
      vectorField pointNorm( NN, vector::zero );
      scalarList sumWeights( NN, 0.0 );
      scalarList faceToPointSumWeights( NN, 0.0 );
      forAll(nepe, i){
        label  curI = local_wall_WallsInletEdges[i];
        point& curP = boundaryPoints[curI];
        const labelList& pNeib = nepe[i];
        
        
        forAll(pNeib, ii){
          label ind = pNeib[ii];
          point& neibP = boundaryPoints[ind];
          vector d2 = neibP - curP;

          scalar mag_d = mag(d2);
          
          vector disp;
          disp.x() = d2.x();// * 0.5;
          disp.y() = d2.y();// * 0.5;
          disp.z() = 0.0;
          displacement[i] += disp;
          
          tol[i] += mag_d; // * pw[j];
          
          // this is for normal
          scalar nw = 1.0 / mag_d;
          //pointNorm[i] += nz ^ (d2*nw);
          //pointNorm[i] += nw * faceNs[ faceI ];
          //sumWeights[i] += nw;
          sumWeights[i] += 1.0;
        }
        
        const labelList& pFaces = plistFaces[curI];
        forAll(pFaces, j){
          label faceI = pFaces[j];
          point& faceC = faceCs[faceI];
          vector d = faceC - curP;
          scalar mag_d = mag(d);
          
          // this is for normal
          scalar nw = 1.0 / mag_d;
          pointNorm[i] += nw * faceNs[ faceI ];
          pointNorm[i].z()=0;
          faceToPointSumWeights[i] += nw;
        }
      }

      syncTools::syncPointList(mesh_, global_WallInletEdges, displacement, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, global_WallInletEdges, tol, plusEqOp<scalar>(), 0.0);
      syncTools::syncPointList(mesh_, global_WallInletEdges, sumWeights, plusEqOp<scalar>(), 0.0);
      syncTools::syncPointList(mesh_, global_WallInletEdges, pointNorm, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, global_WallInletEdges, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
      // normalization
      forAll(pointNorm, i){
        displacement[i] /= sumWeights[i];
        pointNorm[i] /= mag( pointNorm[i] );
      }

      vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);
      
      forAll(local_wall_WallsInletEdges, i){
        label ind = local_wall_WallsInletEdges[i];
        boundaryPoints[ind] += 0.75*projectedDisplacement[i];
      }
      
      scalarField aux_f = mag(projectedDisplacement)/tol;
      scalarField aux_f0(NN, 0.0);

      if( aux_f == aux_f0 ){
        displ_tol = 0.0;
      }
      else{
        displ_tol = average( aux_f );
      }

      reduce(displ_tol, sumOp<scalar>());
      displ_tol = displ_tol / static_cast<double>( Pstream::nProcs() );
      
      if(itt%1000==0){
        Info<<patch_.name()<<"  rlx iter "<<itt
                <<"  tolerance: "<<displ_tol<< nl;
      }
      
      itt+=1;
    }
    
    forAll(globalDisplacement, i){
      globalDisplacement[i] = boundaryPoints[i]-patch_.localPoints()[i];
    }
    
    return globalDisplacement;
  }











  

  
  template<class Patch>
  vectorField wallRelaxation11(Patch& patch_, const scalarListList& weights, double rlxTol){
    pointField saveBoundaryPoints = patch_.localPoints();

    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    const labelListList& ppp = patch_.pointEdges();
    const edgeList& ee = patch_.edges();
    
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
        
        scalar pnd = mag(pointNorm[pointi]);
        pointNorm[pointi] /= pnd;
      }
      
      forAll(local_wall_WallsInletEdges, i){
        label ind = local_wall_WallsInletEdges[i];
        displacement[ind] = vector::zero;
      }
      forAll(local_wall_WallsOutletEdges, i){
        label ind = local_wall_WallsOutletEdges[i];
        displacement[ind] = vector::zero;
      }

      vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);
      //vectorField projectedDisplacement = displacement;
      //vectorField projectedDisplacement1 = transform(I - pointNorm*pointNorm, displacement);
      //vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, projectedDisplacement1);
      
      // fixing edges
      /*
      //vector ny(0,1,0);
      vector nz(0,0,1);

      vectorField inletPrj( global_WallInletEdges.size() );
      vectorField inletPrjN( global_WallInletEdges.size() );
      
      forAll(global_WallInletEdges, i){
        label ind = findIndex( meshPoints, global_WallInletEdges[i]);
        if( ind!=-1 ){
          inletPrj[i] = projectedDisplacement[ ind ];
          inletPrjN[i] = pointNorm[ ind ];
          inletPrjN[i].z() = 0.0;
          inletPrjN[i] /= mag(inletPrjN[i]);
        }
      }

      //inletPrj = transform(I - ny*ny, inletPrj);
      inletPrj = transform(I - inletPrjN*inletPrjN, inletPrj);
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
          outletPrjN[i].z() = 0.0;
          outletPrjN[i] /= mag(outletPrjN[i]);
        }
      }

      //outletPrj = transform(I - ny*ny, outletPrj);
      outletPrj = transform(I - outletPrjN*outletPrjN, outletPrj);
      outletPrj = transform(I - nz*nz, outletPrj);

      forAll(global_WallOutletEdges, i){
        label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
        if( ind!=-1 ){
          projectedDisplacement[ ind ] = outletPrj[i];
        }
      }
      */
      
      //scalar factor = 1;
      scalar factor = 5.5;
      if(itt%2==0) factor = 1.17;
      vectorField finalDisplacement = factor * projectedDisplacement;
      
      boundaryPoints += finalDisplacement;
      
      scalarField aux_f = mag(finalDisplacement)/tol;
      scalarField aux_f0(finalDisplacement.size(), 0.0);
      
      if( aux_f == aux_f0 ){
        displ_tol = 0.0;
      }
      else{
        displ_tol = average( aux_f );
      }
      
      reduce(displ_tol, sumOp<scalar>());
      displ_tol = displ_tol / static_cast<double>( Pstream::nProcs() );
      
      if(itt%1000==0){
        Info<<patch_.name()<<"  rlx iter "<<itt
                <<"  factor "<<factor<<"  tolerance: "<<displ_tol<< nl;
      }
      
      itt+=1;
      
      /*
      if(itt%16000==0){
        break;
      }
      */
      
    }
    Info<< patch_.name() << "  converged in "<< itt 
            << " iterations. Tolerance: " << displ_tol<< nl;
    
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
  vectorField wallRelaxation12(Patch& patch_, const scalarListList& weights, double rlxTol){
    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    const labelListList& ppp = patch_.pointEdges();
    const edgeList& ee = patch_.edges();
    
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
          point& np = boundaryPoints[pI];
          
          vector d = np - curP;
          
          vector disp;
          disp = d * pw[j];
          
          displacement[i] += disp;
          
          scalar mag_d = mag(d);
          tol[i] += mag_d * pw[j];
        }
        
        const labelList& pFaces = plistFaces[i];
        forAll(pFaces, j){
          label faceI = pFaces[j];
          point& faceC = faceCs[faceI];
          vector d = faceC - curP;
          scalar mag_d = mag(d);
          // this is for normal
          scalar nw = 1.0 / mag_d;
          pointNorm[i] += nw * faceNs[ faceI ];
          faceToPointSumWeights[i] += nw;
          
          /*
          if(i==8265){
            Info<<"nw  "<< nw << "   faceNs " << faceNs[ faceI ] <<nl<<nl;
          }
          */
          
        }
        
      }

      // synchronizing over the cyclic and processor boundaries
      syncTools::syncPointList(mesh_, meshPoints, displacement, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, meshPoints, tol, plusEqOp<scalar>(), 0.0);
      
      syncTools::syncPointList(mesh_, meshPoints, pointNorm, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, meshPoints, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
      
      // normalization
      forAll(pointNorm, pointi){
        /*
        if(pointi==8265){
          Info<<"AA "<< pointNorm[pointi] << "  " << faceToPointSumWeights[pointi] <<nl;
        }
        */
        pointNorm[pointi] /= faceToPointSumWeights[pointi];
        
        scalar pnd = mag(pointNorm[pointi]);
        
        pointNorm[pointi] /= pnd;
      }

      vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);
      
      //Info<< (projectedDisplacement - displacement) << nl;
      
      
      // fixing edges
      //vector ny(0,1,0);
      vector nz(0,0,1);

      vectorField inletPrj( global_WallInletEdges.size() );
      vectorField inletPrjN( global_WallInletEdges.size() );
      forAll(global_WallInletEdges, i){
        label ind = findIndex( meshPoints, global_WallInletEdges[i]);
        if( ind!=-1 ){
          inletPrj[i] = projectedDisplacement[ ind ];
          inletPrjN[i] = pointNorm[ ind ];
          inletPrjN[i].z() = 0.0;
          inletPrjN[i] /= mag(inletPrjN[i]);

          //Info << ind << "   " << pointNorm[ind] << "   " << mag(pointNorm[ind]) <<nl;
        }
      }
      //inletPrj = transform(I - ny*ny, inletPrj);
      inletPrj = transform(I - inletPrjN*inletPrjN, inletPrj);
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
          outletPrjN[i].z() = 0.0;
          outletPrjN[i] /= mag(outletPrjN[i]);
        }
      }

      //outletPrj = transform(I - ny*ny, outletPrj);
      outletPrj = transform(I - outletPrjN*outletPrjN, outletPrj);
      outletPrj = transform(I - nz*nz, outletPrj);

      forAll(global_WallOutletEdges, i){
        label ind = findIndex( meshPoints, global_WallOutletEdges[i]);
        if( ind!=-1 ){
          projectedDisplacement[ ind ] = outletPrj[i];
        }
      }
      
      scalar factor = 1;
      //scalar factor = 6.0;
      //if(itt%2==0) factor = 1.17;
      vectorField finalDisplacement = factor * projectedDisplacement;
      
      boundaryPoints += finalDisplacement;
      
      scalarField aux_f = mag(finalDisplacement)/tol;
      scalarField aux_f0(finalDisplacement.size(), 0.0);
      
      if( aux_f == aux_f0 ){
        displ_tol = 0.0;
      }
      else{
        displ_tol = average( aux_f );
      }
      
      reduce(displ_tol, sumOp<scalar>());
      displ_tol = displ_tol / static_cast<double>( Pstream::nProcs() );
      
      if(itt%1000==0) Info<<patch_.name()<<"  rlx iter "<<itt <<"  factor "<<factor<<"  tolerance: "<<displ_tol<< nl;
      
      itt+=1;
      
      if(itt%30000==0) break;
    }
    Info<< patch_.name() << "  converged in "<< itt 
            << " iterations. Tolerance: " << displ_tol<< nl;
    
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
  vectorField wallRelaxation42(Patch& patch_, const scalarListList& weights, double rlxTol, const vectorField &vvff){
    pointField saveBoundaryPoints = patch_.localPoints();

    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    const labelListList& ppp = patch_.pointEdges();
    const edgeList& ee = patch_.edges();
    
    const labelList& wallsTo  = mesh_.boundaryMesh()[wallID].meshPoints();
    
    labelList local_wall_WallEdges;
    labelList global_WallEdges;
    forAll(wallsTo, i){
      label lW = wallsTo[i];

      label inlind = findIndex(meshPoints, lW);
      if( inlind != -1){
        local_wall_WallEdges.append( i );
        global_WallEdges.append( lW );
      }
    }
    
    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    double displ_tol = 1.0;
    int itt = 0;
    
    vectorField displacement(N, vector::zero);
    forAll(global_WallEdges, i){
      label ind = findIndex( meshPoints, global_WallEdges[i]);
      label ind_wall = local_wall_WallEdges[i];
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
      
      forAll(global_WallEdges, i){
        label ind = findIndex( meshPoints, global_WallEdges[i]);
        label ind_wall = local_wall_WallEdges[i];
        if( ind!=-1 ){
          newDisplacement[ ind ] = vvff[ ind_wall ];
        }
      }
      
      scalarField aux_f = mag(newDisplacement - displacement);
      scalarField aux_f0(displacement.size(), 0.0);
      
      if( aux_f == aux_f0 ){
        displ_tol = 0.0;
      }
      else{
        displ_tol = average( aux_f );
      }
      
      scalarField aux_fn = mag(newDisplacement);
      scalarField aux_fn0(displacement.size(), 0.0);
      scalar displ_tol_norm;
      if( aux_fn == aux_fn0 ){
        displ_tol_norm = 0.0;
      }
      else{
        displ_tol_norm = average( aux_fn );
      }
      
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
      
      scalarField aux_f = mag(newDisplacement - displacement);
      scalarField aux_f0(displacement.size(), 0.0);
      
      if( aux_f == aux_f0 ){
        displ_tol = 0.0;
      }
      else{
        displ_tol = average( aux_f );
      }
      
      scalarField aux_fn = mag(newDisplacement);
      scalarField aux_fn0(displacement.size(), 0.0);
      scalar displ_tol_norm;
      if( aux_fn == aux_fn0 ){
        displ_tol_norm = 0.0;
      }
      else{
        displ_tol_norm = average( aux_fn );
      }
      
      
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
  scalarListList calc_weights0(const Patch& patch_, scalar lambda, int dir){

    const labelList& meshPoints = patch_.meshPoints();

    const pointField& boundaryPoints = patch_.localPoints();
    const labelListList& plistFaces = patch_.pointFaces();
    const pointField& faceCs = patch_.faceCentres();
    
    scalarListList weights( boundaryPoints.size() );
    scalarList sumWeights( boundaryPoints.size() );
    
    const labelListList& ppp = patch_.pointEdges();
    const edgeList& ee = patch_.edges();

    forAll(boundaryPoints, i){
      
      
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

      point curP = boundaryPoints[i];
      scalarList& pw = weights[i];
      pw.setSize(lll.size());
      
      scalar sumw = 0.0;
      forAll(nl, j){
        label pI = nl[j];
        const point& np = boundaryPoints[pI];

        vector d = np - curP;
        pw[j] = 1/(mag(d));
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
        
        pw[j] = 1/(mag(d));
        
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

