/* 
 * File:   fieldOperations.H
 * Author: vstar
 * Created on October 16, 2015
 */

#include "fieldOperations.H"
#include "coupledPatchInterpolation.H"


fieldOperations::
fieldOperations(const argList& args, const label patchID):
args_(args),
scalingPatchID_(patchID)
{}


scalar fieldOperations::
getScalingArea(const fvMesh& mesh)
{
  const surfaceScalarField& magSf = mesh.magSf();
  return gSum( magSf.boundaryField()[scalingPatchID_] );
}


scalar fieldOperations::
getScalingFlowRate(const surfaceScalarField& phi)
{
  return mag( gSum(phi.boundaryField()[scalingPatchID_]) );
}


vectorField fieldOperations::
getWallPointMotion(const fvMesh& mesh, const volScalarField& C, 
                   const label movingPatchID)
{
  // interpolate the concentration from cells to wall faces
  coupledPatchInterpolation patchInterpolator
  (
    mesh.boundaryMesh()[movingPatchID], mesh 
  );

  // concentration and normals on the faces
  scalarField pointCface = -C.boundaryField()[movingPatchID].snGrad();
  vectorField pointNface = mesh.boundaryMesh()[movingPatchID].faceNormals();
  
  scalarField motionC = patchInterpolator.faceToPointInterpolate(pointCface);
  vectorField motionN = patchInterpolator.faceToPointInterpolate(pointNface);
  
  // normalize point normals to 1
  forAll(motionN, ii) motionN[ii]/=mag(motionN[ii]);
  
  return motionC*motionN;
}

