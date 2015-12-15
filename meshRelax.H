/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    auxiliary class for the dissolFoam solver.

Description
    It contains all mesh relaxation methods, mesh motion etc.
  
Contributors/Copyright:
    Vitaliy Starchenko (2015)

 
\*---------------------------------------------------------------------------*/

#ifndef meshRelax_H
#define meshRelax_H

#include "fvCFD.H"

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

class meshRelax{

private:

  const float version;
  const word date;

  //fvMesh& mesh_;
  dynamicFvMesh& mesh_;
  scalar deltaT;

  // ID of each patch
  label wallID, inletID, outletID;

  // vetex ID maps: patch to global
  labelList wallsToAll, inletToAll, outletToAll;
  
  // if cyclic is present
  label cycID1, cycID2;
  labelList cycToAll1, cycToAll2;
  labelList local_wall_WallsCycEdges1;
  labelList local_wall_WallsCycEdges2;

  // list of local (local for the wall) IDs for walls- inlet and outlet edges
  labelList local_wall_WallsInletEdges;
  labelList local_wall_WallsOutletEdges;

  // list of global IDs for walls- inlet and outlet edges
  labelList global_WallInletEdges;
  labelList global_WallOutletEdges;
  
  
  labelList local_inlet_WallsInletEdges;
  
  labelList local_outlet_WallsOutletEdges;
  
  // map for the extrapolation of the concentration on the edge
  // it fixes inlet boundary condition inconsistance
  labelListList inletTriple;
  
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  bool dissolDebug;
  bool fixInletWallEdgeDispl;
  scalar rlxTol;
  scalar k_1, k_2; // relaxation acceleration factors
  
  bool variableGrading;
  scalar inigradingZ;
  scalar timeCoefZ;
  int Nz;
  vectorFieldList wallWeights;
  vectorFieldList inletWeights;
  vectorFieldList outletWeights;
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  vector extrapolateVectorLinear(const pointField&, vector&, vector&, const labelList&);
  
  scalar extrapolateConcentrationExp(const pointField&, scalar&, scalar&, const labelList&);
  scalar extrapolateConcentrationExpZ(const pointField&, scalar&, scalar&, const labelList&);
  
  scalar extrapolateConcentrationLinear(const pointField&, scalar&, scalar&, const labelList&);
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
  meshRelax(dynamicFvMesh&);

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  void meshUpdate(vectorField&, Time&);

  void fixIWEdgeDispl( vectorField& );
  
  // fix edge displacement and calculate inlet displ.
  vectorField calculateInletDisplacement( vectorField& );
  
  vectorField calculateOutletDisplacement( vectorField& );
  
  pointField doInletDisplacement(const vectorField&);
  pointField doWallDisplacement(const vectorField&);
  pointField doOutletDisplacement(const vectorField&);
  
  vectorField normalsOnTheEdge();
  
  template<class Patch>
  vectorField edgeRelaxation(Patch& patchIO_, double rlxTol, vectorFieldList& wv){
    
    Patch& patch_ = mesh_.boundaryMesh()[wallID];
    const List<face>& llf = patch_.localFaces();
    const labelListList& plistFaces = patch_.pointFaces();

    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    const labelListList& ppp = patch_.pointEdges();
    const edgeList& ee = patch_.edges();
    
    vectorField globalDisplacement(N, vector::zero);
    
    const labelList& meshPoints = patchIO_.meshPoints();
    const labelList& wallsTo  = patch_.meshPoints();
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
    
    // relax edges
    label NN = local_wall_WallEdges.size();
    label globalNN = NN;
    reduce(globalNN, sumOp<label>());
    
    //Info<<"number of points: "<<globalNN<<nl;
    // @TODO check 2D <-> 3D
    //if(globalNN<5){
    //  return globalDisplacement;
    //}
    
    labelListList nepe;
    forAll(local_wall_WallEdges, i){
      label ind = local_wall_WallEdges[i];

      const labelList& lll = ppp[ind];
      labelList nel;
      forAll(lll, ii){
        const edge& ed = ee[ lll[ii] ];
        label edb = ed.start();
        
        if( edb!=ind && findIndex(local_wall_WallEdges, edb) != -1 ){
          nel.append(edb);
        }
        
        label ede = ed.end();
        if( ede!=ind && findIndex(local_wall_WallEdges, ede) != -1){
          nel.append(ede);
        }
      }
      nepe.append(nel);
    }
    
    double displ_tol = 1.0;
    int itt = 0;
    vector nCyc(1.0,0.0,0.0);
    
    while(displ_tol>rlxTol){
      pointField faceCs = faceCentres(boundaryPoints, llf);
      vectorField faceNs = faceNormals(boundaryPoints, llf);
      
      vectorField displacement(NN, vector::zero);
      scalarField tol(NN, 0.0);
      vectorField pointNorm( NN, vector::zero );
      scalarList sumWeights( NN, 0.0 );
      scalarList faceToPointSumWeights( NN, 0.0 );
      
      forAll(nepe, i){
        label  curI = local_wall_WallEdges[i];
        point& curP = boundaryPoints[curI];
        const labelList& pNeib = nepe[i];
        
        vectorField& curwv = wv[i];
        
        forAll(pNeib, ii){
          label ind = pNeib[ii];
          point& neibP = boundaryPoints[ind];
          vector d2 = neibP - curP;

          scalar mag_d = mag(d2);
          
          vector disp;
          disp.x() = d2.x() * curwv[ii].x(); // * nw; // TODO it is for the weighted 
          disp.y() = d2.y() * curwv[ii].y(); // * nw;
          disp.z() = 0.0;
          displacement[i] += disp;
          tol[i] += mag_d;
          sumWeights[i] += 1;
          //sumWeights[i] += nw;
        }

        
        // stick to cyclic boundary
        //if( findIndex(local_wall_WallsCycEdges1, curI) != -1 ||
        //        findIndex(local_wall_WallsCycEdges2, curI) != -1){
          //vector pdc = displacement[i];
          //displacement[i] = transform(I - nCyc*nCyc, pdc);
        //  displacement[i].x() = 0.0;
        //}
        
        
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
      
      syncTools::syncPointList(mesh_, global_WallEdges, displacement, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, global_WallEdges, tol, plusEqOp<scalar>(), 0.0);
      syncTools::syncPointList(mesh_, global_WallEdges, sumWeights, plusEqOp<scalar>(), 0.0);
      syncTools::syncPointList(mesh_, global_WallEdges, pointNorm, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, global_WallEdges, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
      // normalization
      forAll(pointNorm, i){
        displacement[i] /= sumWeights[i];
        pointNorm[i] /= mag( pointNorm[i] );
      }

      vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);
      
      forAll(local_wall_WallEdges, i){
        label ind = local_wall_WallEdges[i];
        
        boundaryPoints[ind] += projectedDisplacement[i];
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
        Info<<patch_.name()<<"  rlx iter "<<itt<<"  tolerance: "<<displ_tol<< nl;
      }
      
      itt+=1;
    }
    
    forAll(globalDisplacement, i){
      globalDisplacement[i] = boundaryPoints[i]-patch_.localPoints()[i];
    }
    
    return globalDisplacement;
  }

  
  template<class Patch>
  vectorField wallRelaxation(Patch& patch_, const vectorFieldList& weights, double rlxTol){
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
          
          tol[i] += mag_d * mag(pw[j]);
          
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
      
      scalar factor = (itt%2==0) ? k_2 : k_1;
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
      
      //if(itt%40000==0){
      //  break;
      //}
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
  vectorField wallRelaxation2(Patch& patch_, const vectorFieldList& weights, double rlxTol){
    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    double displ_tol = 1.0;
    int itt = 0;
    vector nCyc(1.0,0.0,0.0);
    
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
        const vectorField& pw = weights[i];

        forAll(pFaces, j){
          label faceI = pFaces[j];
          point& faceC = faceCs[faceI];

          vector d = faceC - curP;
          
          scalar mag_d = mag(d);
          
          vector disp;
          
          //if(global_WallInletEdges.size()<5){
            //disp.x() = d.x() * 0.0;
            //disp.y() = d.y() * 0.0;
          //}
          //else{
            disp.x() = d.x() * pw[j].x();
            disp.y() = d.y() * pw[j].y();
          //}
          //disp.x() = d.x() * 0.25;
          //disp.y() = d.y() * 0.25;
          disp.z() = d.z() * pw[j].z();
          displacement[i] += disp;
          
          tol[i] += mag_d * mag(pw[j]);
          
          // this is for normal
          scalar nw = 1.0 / mag_d;
          pointNorm[i] += nw * faceNs[ faceI ];
          faceToPointSumWeights[i] += nw;
        }
        
        // stick to cyclic boundary
        //if( findIndex(local_wall_WallsCycEdges1, i) != -1 ||
        //        findIndex(local_wall_WallsCycEdges2, i) != -1){
          //vector pdc = displacement[i];
          //displacement[i] = transform(I - nCyc*nCyc, pdc);
          //Info<<curP.x()<<nl;
        //  displacement[i].x() = 0.0;
        //}
        
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
      
      /*
      forAll(cycToAll1, i){
        label ind = cycToAll1[i];
        vector pdc = projectedDisplacement[ind];
        projectedDisplacement[ind] = transform(I - nCyc*nCyc, pdc);
      }
      forAll(cycToAll2, i){
        label ind = cycToAll2[i];
        vector pdc = projectedDisplacement[ind];
        projectedDisplacement[ind] = transform(I - nCyc*nCyc, pdc);
      }
      */
      
      scalar factor = (itt%2==0) ? k_2 : k_1;
      
      vectorField finalDisplacement = factor * projectedDisplacement;
      
      boundaryPoints += finalDisplacement;
      
      scalarField aux_f = mag(finalDisplacement/factor)/tol;
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
  vectorField wallRelaxationPointNeib(Patch& patch_, const scalarListList& weights, double rlxTol){
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

      vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);
      
      // fixing edges
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
  vectorField inletOutletRlx(Patch& patch_, double rlxTol, const vectorField &vvff){
    const labelList& meshPoints = patch_.meshPoints();

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
    
    vectorField displacement(N, vector::zero);
    forAll(global_WallEdges, i){
      label ind = findIndex( meshPoints, global_WallEdges[i]);
      label ind_wall = local_wall_WallEdges[i];
      if( ind!=-1 ){
        displacement[ ind ] = vvff[ ind_wall ];
      }
    }
    
    double displ_tol = 1.0;
    int itt = 0;
    while(displ_tol>rlxTol){
      
      scalarList sumWeights( N, 0.0 );
      vectorField newDisplacement(N, vector::zero);
      forAll(newDisplacement, i){
        point& curP = boundaryPoints[i];
        
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
      if(displ_tol_norm!=0.0){
        displ_tol /= displ_tol_norm;
      }
      else{
        displ_tol = 0.0;
      }
        
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
  vectorFieldList calc_weights(const Patch& patch_, scalar lambda){

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
          if( mag(d[ii]) < SMALL  ){
            w[ii] = 1.0;
          }
          else{
            if(ii==2)
              w[ii] = (1 - 0.5 * lambda * d.z()/mag(d.z()) * faceC.z()/mag(faceC.z()) );
            else
              w[ii] = 1.0 / mag( d[ii] );
          }
        }
        
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
  

  template<class Patch>
  scalarListList calc_weights0(const Patch& patch_){

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
  scalarListList calc_weights_faces(const Patch& patch_){
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




  template<class Patch>
  vectorFieldList calc_weights2(const Patch& patch_){

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
      
      //if(pFaces.size()<2){
      //  Info<<faceCs[pFaces[0]]<<"  "<<curP<<nl;
      //}
      
      vector sumw = vector::zero;
      forAll(pFaces, j){
        label faceI = pFaces[j];
        point faceC = faceCs[faceI];

        vector d = faceC - curP;
        
        vector w = vector::zero;
        
        for(int ii=0; ii<3; ii++){
          if(mag(d[ii])>SMALL){
            w[ii] = 1.0 / mag( d[ii] );
          }
          else{
            w[ii] = 0.0;
          }
        }
        
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
      forAll(pw, j){
        vector &ppw = pw[j];
        for(int ii=0; ii<3; ii++){
          if( mag(sw[ii])>SMALL ){
            ppw[ii] /= sw[ii];
          }
        }
        
      }
    }

    return weights;
  }



  template<class Patch>
  vectorFieldList calc_edge_weights(const Patch& patchIO_){
    
    const Patch& patch_ = mesh_.boundaryMesh()[wallID];

    pointField boundaryPoints = patch_.localPoints();
    
    const labelListList& ppp = patch_.pointEdges();
    const edgeList& ee = patch_.edges();
    
    const labelList& meshPoints = patchIO_.meshPoints();
    const labelList& wallsTo  = patch_.meshPoints();
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
    
    // relax edges
    label NN = local_wall_WallEdges.size();
    labelListList nepe;
    forAll(local_wall_WallEdges, i){
      label ind = local_wall_WallEdges[i];

      const labelList& lll = ppp[ind];
      labelList nel;
      forAll(lll, ii){
        const edge& ed = ee[ lll[ii] ];
        label edb = ed.start();
        
        if( edb!=ind && findIndex(local_wall_WallEdges, edb) != -1 ){
          nel.append(edb);
        }
        
        label ede = ed.end();
        if( ede!=ind && findIndex(local_wall_WallEdges, ede) != -1){
          nel.append(ede);
        }
      }
      nepe.append(nel);
    }
    
    vectorFieldList weights( NN );
    vectorField sumWeights( NN, vector::zero );

    forAll(nepe, i){
      label  curI = local_wall_WallEdges[i];
      point& curP = boundaryPoints[curI];
      const labelList& pNeib = nepe[i];

      vectorField& pw = weights[i];
      pw.setSize(pNeib.size());
      
      vector sumw = vector::zero;
      forAll(pNeib, ii){
        label ind = pNeib[ii];
        point& neibP = boundaryPoints[ind];
        vector d2 = neibP - curP;

        vector w = vector::zero;
        
        for(int j=0; j<3; j++){
          if(d2[j]!=0)
            w[j] = 1/ std::abs( d2[j] );
          else
            w[j] = 0.0;
        }
        pw[ii] = w;
        sumw += pw[ii];
      }
      sumWeights[i] = sumw;
    }
    
    syncTools::syncPointList(mesh_, global_WallEdges, sumWeights, plusEqOp<vector>(), vector::zero);
    
    forAll(weights, i){
      vectorField& pw = weights[i];
      vector &sw = sumWeights[i];
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

