/* 
 * File:   fieldOperations.H
 * Author: vstar
 * Created on October 16, 2015
 */

#include "fieldOperations.H"

#include "coupledPatchInterpolation.H"

#include "pointMesh.H"
#include "pointPatchField.H"
#include "volPointInterpolation.H"

fieldOperations::fieldOperations()
{

}

scalar fieldOperations::getConstFlowRateFactor(const fvMesh& mesh, 
        const volVectorField& U, 
        scalar areaCoef0)
{

  scalar magU = mag( fvc::domainIntegrate( U ).value() );

  vectorField po = mesh.points();
  scalar maxZ = max( po.component(vector::Z) );
  reduce(maxZ, maxOp<scalar>());
  scalar minZ = min( po.component(vector::Z) );
  reduce(minZ, minOp<scalar>());

  scalar nU = magU / ((maxZ-minZ) * areaCoef0);

  return nU;
}

scalar fieldOperations::calcDanckwerts(const fvMesh& mesh,
        const volVectorField& U,
        volScalarField& C,
        label inletID,
        scalar D)
{
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

    scalar aa = ( mag(U[fcL].z()) * del) / D;
    
    //newC[ii] = (1.0+aa*C[fcL])/(1.0+aa);
    newC[ii] = (aa+C[fcL])/(aa+1.0);
  }
  scalarField diff = newC - oldC;

  C.boundaryField()[inletID]==newC;
  
  scalar ttt  = mag( gSum( diff ) );
  //reduce(ttt, sumOp<scalar>());
  scalar tttn = mag( gSum( newC ) );
  //reduce(tttn, sumOp<scalar>());

  scalar tttol = 0.0;
  if( tttn!=0.0 ) tttol = ttt/tttn;
        
  return tttol;
}


vectorField fieldOperations::getWallPointMotion(const fvMesh& mesh,const volScalarField& C,
                                               scalar l_T, label wallID
){
  // interpolate the concentration from cells to wall faces
  coupledPatchInterpolation patchInterpolator( mesh.boundaryMesh()[wallID], mesh );

  // concentration and normals on the faces
  scalarField pointCface = -C.boundaryField()[wallID].snGrad();
  vectorField pointNface = mesh.boundaryMesh()[wallID].faceNormals();
  
  // ----------------------------------------------------------------
  
  scalarField motionC = patchInterpolator.faceToPointInterpolate(pointCface);
  vectorField motionN = patchInterpolator.faceToPointInterpolate(pointNface);
  
  // normalize point normals to 1
  forAll(motionN, ii) motionN[ii]/=mag(motionN[ii]);
  
  return (l_T*motionC*motionN);
}

scalar fieldOperations::getInletArea(const argList& args, bool constFlux){
  Foam::Time timeTmp(Foam::Time::controlDictName, args);
  if(!constFlux){
    Foam::instantList timeDirs = Foam::timeSelector::select0(timeTmp, args);
    timeTmp.setTime(timeDirs[0], 0);
  }
  
  Foam::fvMesh meshTmp
  (
      Foam::IOobject
      (
          Foam::fvMesh::defaultRegion,
          timeTmp.timeName(),
          timeTmp,
          Foam::IOobject::MUST_READ
      )
  );
  
  scalar areaCoef = 0.0;
  const surfaceScalarField& magSf2 = meshTmp.magSf();
  label inletID = meshTmp.boundaryMesh().findPatchID("inlet");
  forAll(magSf2.boundaryField()[inletID], facei){
    areaCoef += magSf2.boundaryField()[inletID][facei];
  }
  reduce(areaCoef, sumOp<scalar>());
  return areaCoef;
}
