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
typedef List<scalarField> scalarFieldList;

class meshRelax
{
public:
  
  // Constructors
  meshRelax(dynamicFvMesh&, const argList&);

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // public functions
  void meshUpdate(vectorField&, Time&);

  word get_version() const;
  word get_date() const;
  
private:
  // the version of meshRelax and date it was modified
  const word version;
  const word date;
  ////////////////////////////////////////////////////
  
  // dictionary variables
  bool dissolDebug;
  bool fixInletWallEdgeDispl;
  scalar rlxTol;
  
  // relaxation acceleration factors for the wall relaxation
  scalar k_1, k_2; 
  int q_2;
  int q_norm_recalc;
  
  // relaxation acceleration factors for the edge relaxation
  scalar k_1edge, k_2edge;
  int q_2edge;
  int q_edge_norm_recalc;
  
  bool variableGrading;
  scalar inigradingZ;
  scalar timeCoefZ;
  int Nz;
  ////////////////////////////////////////////////////
  dynamicFvMesh& mesh_;
  scalar deltaT;
  // ID of each patch
  label wallID, inletID, outletID;
  
  // list of local (local for the wall) IDs for walls- inlet and outlet edges
  labelList local_wall_WallsInletEdges;
  labelList local_wall_WallsOutletEdges;

  // list of global IDs for walls- inlet and outlet edges
  labelList global_WallInletEdges;
  labelList global_WallOutletEdges;
  
  labelList local_inlet_WallsInletEdges;
  labelList local_outlet_WallsOutletEdges;
  
  // if cyclic is present
  labelList local_wall_WallsCycEdges1;
  labelList local_wall_WallsCycEdges2;

  
  // map for the extrapolation of the concentration on the edge
  // it fixes inlet boundary condition inconsistance
  labelListList inletTriple;
  
  vectorFieldList wallWeights;
  scalarFieldList inletWeights;
  scalarFieldList outletWeights;
  
  ///////////////////////////////////////////////////////////////////////////
  // private functions
  ///////////////////////////////////////////////////////////////////////////
  vector extrapolateVectorLinear(const pointField&, vector&, vector&, const labelList&);
  
  scalar extrapolateConcentrationExp(const pointField&, scalar&, scalar&, const labelList&);
  scalar extrapolateConcentrationExpZ(const pointField&, scalar&, scalar&, const labelList&);
  
  scalar extrapolateConcentrationLinear(const pointField&, scalar&, scalar&, const labelList&);
  scalar extrapolateConcentrationLinearZ(const pointField&, scalar&, scalar&, const labelList&);

  pointField faceCentres(const pointField& points, const List<face>& flist) const;
  vectorField faceNormals(const pointField& points, const List<face>& flist) const;

  vectorField localFaceToPointNormalInterpolate(
          const pointField& points,
          const pointField& faceCs,
          const vectorField& faceNs,
          const labelListList& pointFaces,
          const labelList& meshPoints

  ) const;
  
  void fixIWEdgeDispl( vectorField& );
  // fix edge displacement and calculate inlet displ.
  vectorField calculateInletDisplacement( vectorField& );
  vectorField calculateOutletDisplacement( vectorField& );
  
  pointField doInletDisplacement(const vectorField&);
  pointField doWallDisplacement(const vectorField&);
  pointField doOutletDisplacement(const vectorField&);
  
  vectorField normalsOnTheEdge();
  
  // common points on label lists, result is two lists of local (for list1) and global labels
  void commonPoints(const labelList& list1, const labelList& list2, 
                    labelList& localList, labelList& globalList);

  void neighborListEdge(const labelList& list, const edgeList& eList, 
                        const labelListList& pointEdgesLi, labelListList& neLi);
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // set up pairs for concentration extrapollation on the edge
  void setUpPairsConc();

  // setup lists needed for mesh relaxation
  void setUpLists();

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  template<class Patch>
  vectorField edgeRelaxation(Patch& patchIO_, const scalarFieldList& wv)
  {
    Patch& wall_patch = mesh_.boundaryMesh()[wallID];
    const List<face>& llf = wall_patch.localFaces();
    const labelListList& plistFaces = wall_patch.pointFaces();

    pointField boundaryPoints = wall_patch.localPoints();
    
    const labelListList& pointEdgesLi = wall_patch.pointEdges();
    const edgeList& ee = wall_patch.edges();
    
    const labelList& meshPoints = patchIO_.meshPoints();
    const labelList& wallsTo  = wall_patch.meshPoints();
    labelList local_wall_WallEdges, global_WallEdges;
    commonPoints(wallsTo, meshPoints, local_wall_WallEdges, global_WallEdges);
    
    // relax edges
    label NN = local_wall_WallEdges.size();
    //label globalNN = NN;
    //reduce(globalNN, sumOp<label>());
    
    labelListList nepe;
    neighborListEdge(local_wall_WallEdges, ee, pointEdgesLi, nepe);
    
    double displ_tol = 1.0;
    int itt = 0;
    vectorField pointNorm( NN, vector::zero );
    scalarList faceToPointSumWeights( NN, 0.0 );
    while(displ_tol>rlxTol)
    {
      if(itt%q_edge_norm_recalc==0){
        forAll(pointNorm, pni)
        {
          pointNorm[pni] = vector::zero;
          faceToPointSumWeights[pni] = 0.0;
        }
      }
      
      pointField faceCs = faceCentres(boundaryPoints, llf);
      vectorField faceNs = faceNormals(boundaryPoints, llf);
      
      vectorField displacement(NN, vector::zero);
      scalarField tol(NN, 0.0);
      scalarList sumWeights( NN, 0.0 );
      
      forAll(nepe, i){
        label  curI = local_wall_WallEdges[i];
        point& curP = boundaryPoints[curI];
        const labelList& pNeib = nepe[i];
        
        const scalarField& curwv = wv[i];
        
        forAll(pNeib, ii){
          label ind = pNeib[ii];
          point& neibP = boundaryPoints[ind];
          vector d2 = neibP - curP;

          scalar mag_d = mag(d2);
          
          displacement[i] += curwv[ii] * d2;
          tol[i] += mag_d;
          sumWeights[i] += 1;
        }
        
        // stick to cyclic boundary
        if( findIndex(local_wall_WallsCycEdges1, curI) != -1 ||
            findIndex(local_wall_WallsCycEdges2, curI) != -1){
          displacement[i].x() = 0.0;
        }
        
        if(itt%q_edge_norm_recalc==0){
          const labelList& pFaces = plistFaces[curI];
          forAll(pFaces, j){
            label ind = pNeib[j];
            point& neibP = boundaryPoints[ind];
            vector d2 = neibP - curP;     // vector from the current point to its neighbor
            
            label faceI = pFaces[j]; 
            vector fnn = faceNs[ faceI ]; // face normal
            fnn.z() = 0.0;                // projection to the inlet surface ~ & (0,0,1)
            
            // correction of the normal, otherwise realN = fnn
            scalar middi = mag(d2)/2.0;   // half distance to neighbor
            point midpo = curP + d2/2.0;  // midpoint
            plane pll(midpo, d2);           // plane perpendicular to the edge via midpoint
            point endNorm = midpo + fnn;  // end of the face normal projected on the inlet
            point projp = pll.nearestPoint(endNorm); // projection of endNorm onto pll plane
            
            vector realN = projp - midpo; // normal to the edge between current point
                                          // and its neighbor in the inlet plane facing
                                          // outside the fracture (mag(realN)!=1 !!!)
            
            scalar nw = 1.0 / middi;
            pointNorm[i] += nw * realN;
            pointNorm[i].z()=0;
            faceToPointSumWeights[i] += nw;
          }
        }
      }
      
      syncTools::syncPointList(mesh_, global_WallEdges, displacement, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, global_WallEdges, tol, plusEqOp<scalar>(), 0.0);
      syncTools::syncPointList(mesh_, global_WallEdges, sumWeights, plusEqOp<scalar>(), 0.0);
      
      if(itt%q_edge_norm_recalc==0){
        syncTools::syncPointList(mesh_, global_WallEdges, pointNorm, plusEqOp<vector>(), vector::zero);
        syncTools::syncPointList(mesh_, global_WallEdges, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);
        // normalization
        forAll(pointNorm, i){
          displacement[i] /= sumWeights[i];
          pointNorm[i] /= mag( pointNorm[i] );
        }
      }

      vectorField projectedDisplacement = transform(I - pointNorm*pointNorm, displacement);
      
      scalar factor = (itt%q_2edge==0) ? k_2edge : k_1edge;
      projectedDisplacement*=factor;
      
      // stick to cyclic boundary
      forAll(local_wall_WallsCycEdges1, i){
        label ind = local_wall_WallsCycEdges1[i];
        label loc_ind = findIndex(local_wall_WallEdges, ind);
        if( loc_ind != -1){
          projectedDisplacement[loc_ind].x() = 0.0;
        }
      }
      forAll(local_wall_WallsCycEdges2, i){
        label ind = local_wall_WallsCycEdges2[i];
        label loc_ind = findIndex(local_wall_WallEdges, ind);
        if( loc_ind != -1){
          projectedDisplacement[loc_ind].x() = 0.0;
        }
      }
      
      forAll(local_wall_WallEdges, i){
        label ind = local_wall_WallEdges[i];
        boundaryPoints[ind] += projectedDisplacement[i];
      }
      
      displ_tol = gAverage( mag(projectedDisplacement)/tol );
      
      if(itt%1000==0){
        Info << patchIO_.name() << "  rlx iter " << itt
             << "  tolerance: " << displ_tol << endl;
      }
      
      itt++;
    }
    
    return (boundaryPoints-wall_patch.localPoints());
  }

  template<class Patch>
  vectorField wallRelaxation(Patch& patch_, const vectorFieldList& weights){
    const List<face>& llf = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();
    const labelListList& plistFaces = patch_.pointFaces();

    pointField boundaryPoints = patch_.localPoints();
    int N = boundaryPoints.size();
    
    double displ_tol = 1.0;
    int itt = 0;
    
    vectorField pointNorm( N, vector::zero );
    scalarList faceToPointSumWeights( N, 0.0 );
    while(displ_tol>rlxTol)
    {
      if(itt%q_norm_recalc==0){
        forAll(pointNorm, pni)
        {
          pointNorm[pni] = vector::zero;
          faceToPointSumWeights[pni] = 0.0;
        }
      }
      // calculate current face centers
      pointField faceCs = faceCentres(boundaryPoints, llf);
      vectorField faceNs = faceNormals(boundaryPoints, llf);

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
          if(itt%q_norm_recalc==0){
            scalar nw = 1.0 / mag_d;
            pointNorm[i] += nw * faceNs[ faceI ];
            faceToPointSumWeights[i] += nw;
          }
        }
        
        // stick to cyclic boundary
        if( findIndex(local_wall_WallsCycEdges1, i) != -1 ||
            findIndex(local_wall_WallsCycEdges2, i) != -1){
          displacement[i].x() = 0.0;
        }
        
      }
      
      // synchronizing over the cyclic and processor boundaries
      syncTools::syncPointList(mesh_, meshPoints, displacement, plusEqOp<vector>(), vector::zero);
      syncTools::syncPointList(mesh_, meshPoints, tol, plusEqOp<scalar>(), 0.0);
      
      if(itt%q_norm_recalc==0){
        syncTools::syncPointList(mesh_, meshPoints, pointNorm, plusEqOp<vector>(), vector::zero);
        syncTools::syncPointList(mesh_, meshPoints, faceToPointSumWeights, plusEqOp<scalar>(), 0.0);

        // normalization
        forAll(pointNorm, pointi){
          pointNorm[pointi] /= faceToPointSumWeights[pointi];
          pointNorm[pointi] /= mag(pointNorm[pointi]);
        }
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
      
      forAll(local_wall_WallsCycEdges1, i){
        label ind = local_wall_WallsCycEdges1[i];
        projectedDisplacement[ind].x() = 0.0;
      }
      forAll(local_wall_WallsCycEdges2, i){
        label ind = local_wall_WallsCycEdges2[i];
        projectedDisplacement[ind].x() = 0.0;
      }
      
      scalar factor = (itt%q_2==0) ? k_2 : k_1;
      
      vectorField finalDisplacement = factor * projectedDisplacement;
      
      boundaryPoints += finalDisplacement;
      
      displ_tol = gAverage( mag(finalDisplacement/factor)/tol );
      
      if(itt%1000==0){
        Info << patch_.name() << "  rlx iter " << itt
             << "  tolerance: " << displ_tol << endl;
      }
      
      itt++;
    }
    Info << patch_.name() << "  converged in " << itt 
         << " iterations. Tolerance: " << displ_tol<< endl;
    
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
  vectorField wallRelaxationOnPoints(Patch& patch_, const scalarListList& weights)
  {
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
      
      if(itt%1000==0){
        Info << patch_.name() << "  rlx iter " << itt
             <<"  tolerance: "<< displ_tol << endl;
      }
      
      itt+=1;
      
      if(itt%30000==0) break;
    }
    Info << patch_.name() << "  converged in "<< itt 
         << " iterations. Tolerance: " << displ_tol<< endl;
    
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
  vectorField inletOutletRlx(Patch& patch_, const vectorField &vvff)
  {
    const labelList& meshPoints = patch_.meshPoints();

    const labelListList& ppp = patch_.pointEdges();
    const edgeList& ee = patch_.edges();
    
    const labelList& wallsTo  = mesh_.boundaryMesh()[wallID].meshPoints();
    labelList local_wall_WallEdges, global_WallEdges;
    commonPoints(wallsTo, meshPoints, local_wall_WallEdges, global_WallEdges);
    
    const pointField& boundaryPoints = patch_.localPoints();
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
    while(displ_tol>rlxTol)
    {
      scalarList sumWeights( N, 0.0 );
      vectorField newDisplacement(N, vector::zero);
      forAll(newDisplacement, i){
        const point& curP = boundaryPoints[i];
        
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
          const point& np = boundaryPoints[pI];
          
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
      
      displ_tol = gAverage( mag(newDisplacement - displacement) );
      
      scalar displ_tol_norm = gAverage( mag(newDisplacement) );
      
      if(displ_tol_norm > SMALL){
        displ_tol /= displ_tol_norm;
      }
      else{
        displ_tol = 0.0;
      }
        
      displacement = newDisplacement;
      
      if(itt%1000==0){
        Info << patch_.name() << "  rlx iter " << itt
             << "  tolerance: " << displ_tol << endl;
      }
      itt+=1;
    }
    
    forAll(global_WallInletEdges, i){
      label ind = findIndex( meshPoints, global_WallInletEdges[i]);
      if( ind!=-1 ){
        displacement[ ind ] = vector::zero;
      }
    }
      
    Info << patch_.name() << "  converged in "<< itt 
         << " iterations. Tolerance: " << displ_tol<< endl;
    
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
  vectorFieldList calc_weights2(const fvMesh& meshTmp, const Patch& patch_){

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

    syncTools::syncPointList( meshTmp, meshPoints, sumWeights, plusEqOp<vector>(), vector::zero);
    
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
  scalarFieldList calc_edge_weights(const fvMesh& meshTmp, const Patch& patchIO_)
  {
    const Patch& patch_ = meshTmp.boundaryMesh()[wallID];
    pointField boundaryPoints = patch_.localPoints();
    
    const labelList& wallsTo  = patch_.meshPoints();
    const labelList& meshPoints = patchIO_.meshPoints();
    labelList local_wall_WallEdges, global_WallEdges;
    commonPoints(wallsTo, meshPoints, local_wall_WallEdges, global_WallEdges);
    
    const labelListList& pointEdgesLi = patch_.pointEdges();
    const edgeList& ee = patch_.edges();
    labelListList nepe;
    neighborListEdge(local_wall_WallEdges, ee, pointEdgesLi, nepe);
    
    label NN = local_wall_WallEdges.size();
    scalarFieldList weights( NN );
    scalarField sumWeights( NN, 0.0 );

    forAll(nepe, i){
      label  curI = local_wall_WallEdges[i];
      point& curP = boundaryPoints[curI];
      const labelList& pNeib = nepe[i];

      scalarField& pw = weights[i];
      pw.setSize(pNeib.size());
      
      scalar sumw = 0.0;
      forAll(pNeib, ii){
        label ind = pNeib[ii];
        point& neibP = boundaryPoints[ind];
        vector d2 = neibP - curP;
        scalar magd2 = mag(d2);

        scalar w = 0.0;
        if( magd2>SMALL) w = 1.0 / magd2;
        
        pw[ii] = w;
        sumw += pw[ii];
      }
      sumWeights[i] = sumw;
    }
    
    syncTools::syncPointList(meshTmp, global_WallEdges, sumWeights, plusEqOp<scalar>(), 0.0);
    
    forAll(weights, i){
      scalarField& pw = weights[i];
      scalar &sw = sumWeights[i];
      forAll(pw, j){
        if( mag(sw)>SMALL ) pw[j] /= sw;
      }
    }

    return weights;
  }

};

#endif

