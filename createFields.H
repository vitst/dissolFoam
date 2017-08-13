//  Create fields at runTime.timeName()

  Info << "Reading field p" << endl;

  volScalarField p
  (
    IOobject
    (
      "p",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
    ),
    mesh
  );


  Info << "Reading field U" << endl;

  volVectorField U
  (
    IOobject
    (
      "U",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
    ),
    mesh
  );


  Info << "Reading field C" << endl;

  volScalarField C
  (
    IOobject
    (
      "C",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
    ),
    mesh
  );


  Info << "Creating field fluxC" << endl;

  IOobject::writeOption writeFluxCOpt = 
          (writeFluxC) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE;

  surfaceScalarField fluxC
  (
    IOobject
    (
      "fluxC",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      writeFluxCOpt
    ),
    mesh,
    dimensionedScalar("fluxC", dimensionSet(0,1,-1,0,0,0,0), scalar(0))
  );

  # include "createPhi.H"

  label  pRefCell  = 0;
  scalar pRefValue = 0.0;
  setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue );

