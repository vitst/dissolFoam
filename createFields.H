  // reading dissolFoam and transportProperties dictionary
  const word dissolDictName("dissolFoamDict");
  IOdictionary dissolProperties
  (
    IOobject
    (
      dissolDictName,
      runTime.system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );
  
  IOdictionary transportProperties
  (
      IOobject
      (
          "transportProperties",
          runTime.constant(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      )
  );

  bool gradCwrite(dissolProperties.lookupOrDefault<bool>("gradCwrite", false));

  Info<< "Reading transportProperties\n" << nl;

  dimensionedScalar nu
  (
      transportProperties.lookup("nu")
  );
  dimensionedScalar D
  (
      transportProperties.lookup("D")
  );
  
  scalar l_T; // l_T=D/(k*h_0)
  if( !transportProperties.readIfPresent<scalar>("l_T", l_T) ){
    SeriousErrorIn("main")
            <<"There is no l_T parameter in transportProperties dictionary"
            <<exit(FatalError);
  }


  Info<< "Reading field p\n" << endl;

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


  Info<< "Reading field U\n" << endl;

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


  Info<< "Reading field C\n" << endl;

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


  Info<< "Creating field gradC\n" << endl;

  IOobject::writeOption gradCwriteOpt = (gradCwrite) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE;

  volScalarField maggradC
  (
      IOobject
      (
          "gradC",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          gradCwriteOpt
      ),
      mag(-fvc::grad(C))
  );


  # include "createPhi.H"

  label pRefCell = 0;
  scalar pRefValue = 0.0;
  setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue );


  /*
  singlePhaseTransportModel laminarTransport(U, phi);

  autoPtr<incompressible::RASModel> turbulence
  (
      incompressible::RASModel::New(U, phi, laminarTransport)
  );
  */
    
  Info<< "Reading dissolFoam flags\n" << nl;
  
  bool inertia;
  if( !dissolProperties.readIfPresent<bool>("inertia", inertia) ){
    SeriousErrorIn("main")
            <<"There is no inertia parameter in dissolFoam dictionary"
            <<exit(FatalError);
  }

  bool constFlux;
  if( !dissolProperties.readIfPresent<bool>("constFlux", constFlux) ){
    SeriousErrorIn("main")
            <<"There is no constFlux parameter in dissolFoam dictionary"
            <<exit(FatalError);
  }
