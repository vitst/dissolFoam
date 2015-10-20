/* 
 * File:   FieldOperations.H
 * Author: vstar
 * Created on October 16, 2015
 */

#include "FieldOperations.H"

#include "coupledPatchInterpolation.H"

FieldOperations::FieldOperations()
{

}

scalar FieldOperations::getConstFlowRateFactor(const fvMesh& mesh, const volVectorField& U, scalar areaCoef0){

  scalar magU = mag( fvc::domainIntegrate( U ).value() );

  vectorField po = mesh.points();
  scalar maxZ = max( po.component(vector::Z) );
  reduce(maxZ, maxOp<scalar>());
  scalar minZ = min( po.component(vector::Z) );
  reduce(minZ, minOp<scalar>());

  scalar nU = magU / (maxZ-minZ) / areaCoef0;

  return nU;
}

scalar FieldOperations::setConstFlowRateFactor(const fvMesh& mesh, const volVectorField& U,
                                volScalarField& C,
                                label inletID,
                                scalar D
        ){
  // C field at the inlet
  const scalarField& oldC = C.boundaryField()[inletID];
  scalarField newC(oldC.size(), 0.0);

  const labelList& fc = mesh.boundaryMesh()[inletID].faceCells();
  const vectorField& fcr = mesh.boundaryMesh()[inletID].faceCentres();
  const vectorField& ccr = mesh.cellCentres();
  forAll(fc, ii){
    label fcL = fc[ii]; // label of the cell the face belong to
    vector vdel = fcr[ii]-ccr[fcL];
    scalar del = mag(vdel.z());

    scalar aa = D / ( mag(U[fcL].z()) * del);
    
    newC[ii] = (1.0+aa*C[fcL])/(1.0+aa);
  }
  scalarField diff = newC - oldC;

  C.boundaryField()[inletID]==newC;
  
  scalar ttt  = mag( gSum( diff ) );
  reduce(ttt, sumOp<scalar>());
  scalar tttn = mag( gSum( newC ) );
  reduce(tttn, sumOp<scalar>());

  scalar tttol = 0.0;
  if( tttn!=0.0 ) tttol = ttt/tttn;


        
  return tttol;
}


vectorField FieldOperations::getWallPointMotion(const fvMesh& mesh,const volScalarField& C,
                                               scalar l_T, label wallID
){
  // interpolate the concentration from cells to wall faces
  coupledPatchInterpolation patchInterpolator( mesh.boundaryMesh()[wallID], mesh );

  // concentration and normals on the faces
  scalarField pointCface = -C.boundaryField()[wallID].snGrad();
  vectorField pointNface = mesh.boundaryMesh()[wallID].faceNormals();

  scalarField motionC = patchInterpolator.faceToPointInterpolate(pointCface);
  vectorField motionN = patchInterpolator.faceToPointInterpolate(pointNface);
  // normalize point normals to 1
  forAll(motionN, ii) motionN[ii]/=mag(motionN[ii]);
  
  return (l_T*motionC*motionN);
}

scalar FieldOperations::getInletArea(const fvMesh& mesh){
  scalar areaCoef = 0.0;
  const surfaceScalarField& magSf2 = mesh.magSf();
  label inletID = mesh.boundaryMesh().findPatchID("inlet");
  forAll(magSf2.boundaryField()[inletID], facei){
    areaCoef += magSf2.boundaryField()[inletID][facei];
  }
  reduce(areaCoef, sumOp<scalar>());
  return areaCoef;
}
