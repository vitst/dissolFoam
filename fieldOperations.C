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

fieldOperations::fieldOperations(const argList& args, const label patchID)
:
args_(args),
scalingPatchID_(patchID)
{
  scalingFlowRateT0calculated = false;
  scalingFlowRateT0 = 0.0;
}


scalar fieldOperations::getScalingFlowRateT0(const surfaceScalarField& phi)
{
  if( !scalingFlowRateT0calculated ){
    Foam::Time timeTmp(Foam::Time::controlDictName, args_);
    Foam::instantList timeDirs = Foam::timeSelector::select0(timeTmp, args_);
    timeTmp.setTime(timeDirs[0], 0);

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

    surfaceScalarField phiTmp
    (
        IOobject
        (
            "phi",
            timeTmp.timeName(),
            meshTmp,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        phi
    );
    
    // in case 0 time does not exist
    if( timeTmp.timeName()!="0" ){
      SeriousErrorIn("fieldOperations::getInletFlowRateT0")
              <<"There is no 0 time directory. Check your decomposition as well!"
              <<exit(FatalError);
    }
  
    scalingFlowRateT0 = mag( gSum( phiTmp.boundaryField()[scalingPatchID_] ) );
    scalingFlowRateT0calculated = true;
  }
  return scalingFlowRateT0;
}

scalar fieldOperations::getScalingFlowRate(const surfaceScalarField& phi, bool updateFlowRate)
{
  return ((updateFlowRate) ? mag( gSum(phi.boundaryField()[scalingPatchID_]) ):getScalingFlowRateT0(phi));
}


vectorField fieldOperations::getWallPointMotion(const fvMesh& mesh,const volScalarField& C,
                                               const scalar l_T, const label wallID)
{
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

scalar fieldOperations::getScalingAreaT0(){
  Foam::Time timeTmp(Foam::Time::controlDictName, args_);
  Foam::instantList timeDirs = Foam::timeSelector::select0(timeTmp, args_);
  timeTmp.setTime(timeDirs[0], 0);
  
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
  
  if( timeTmp.timeName()!="0" ){
    SeriousErrorIn("fieldOperations::getScalingAreaT0")
            <<"There is no 0 time directory. Check your decomposition as well!"
            <<exit(FatalError);
  }
  
  const surfaceScalarField& magSf2 = meshTmp.magSf();
  return gSum( magSf2.boundaryField()[scalingPatchID_] );
}
