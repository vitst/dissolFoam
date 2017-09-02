  Info << "Reading dissolFoam flags" << nl;
  
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
  
  bool writeFluxC;
  if( !dissolProperties.readIfPresent<bool>("writeFluxC", writeFluxC) ){
    SeriousErrorIn("main")
          << "There is no writeFluxC parameter in dissolFoam dictionary"
          << exit(FatalError);
  }

  bool inertia;
  if( !dissolProperties.readIfPresent<bool>("inertia", inertia) ){
    SeriousErrorIn("main")
          << "There is no inertia parameter in dissolFoam dictionary"
          << exit(FatalError);
  }

  bool limitFlux;
  if( !dissolProperties.readIfPresent<bool>("limitFlux", limitFlux) ){
    SeriousErrorIn("main")
          << "There is no limitFlux parameter in dissolFoam dictionary"
          << exit(FatalError);
  }
  if(limitFlux)
    Info << "***************************************************" << nl
         << "'limitFlux' parameter is set to 'true'" << nl
         << "U & phi fields rescaled if outlet flux exceeds Qlim" << nl
         << "***************************************************" << nl
         << endl;

  bool constFlux;
  if( !dissolProperties.readIfPresent<bool>("constFlux", constFlux) ){
    SeriousErrorIn("main")
          << "There is no constFlux parameter in dissolFoam dictionary"
          << exit(FatalError);
  }
  if(constFlux)
    Info << "***************************************************" << nl
         << "'constFlux' parameter is set to 'true'" << nl
         << "U & phi fields always rescaled to Qlim" << nl
         << "***************************************************" << nl
         << endl;

  scalar Qlim = 1.0;
  if( !dissolProperties.readIfPresent<scalar>("Qlim", Qlim) ){
    SeriousErrorIn("main")
          << "There is no Qlim parameter in dissolFoam dictionary"
          << exit(FatalError);
  }


  Info << "Reading transportProperties" << nl;

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

  dimensionedScalar Kinv
  (
      transportProperties.lookup("Kinv")
  );

  dimensionedScalar nu
  (
      transportProperties.lookup("nu")
  );

  dimensionedScalar D
  (
      transportProperties.lookup("D")
  );

  dimensionedScalar lR
  (
      transportProperties.lookup("lR")
  );

  Info << "Reading convDiff subdirectory" << nl;

  dictionary convDiff = mesh.solutionDict().subDict("convDiff");

  double tolerance = 0;         // convergence criteria for convection-diffusion
  if( !convDiff.readIfPresent<double>("tolerance", tolerance) ){
    SeriousErrorIn("main")
        << "No tolerance parameter in convDiff subdictionary"
        << exit(FatalError);
  }
  int maxIter = 0;             // maximum number of iterations
  if( !convDiff.readIfPresent<int>("maxIter", maxIter) ){
    SeriousErrorIn("main")
        << "No maxIter parameter in convDiff subdictionary"
        << exit(FatalError);
  }


  Info << "***************************************************" << nl;
  Info << "transportProperties, Kinv:  " << Kinv                << nl;
  Info << "transportProperties, nu:    " << nu                  << nl;
  Info << "transportProperties, D:     " << D                   << nl;
  Info << "transportProperties, lR:    " << lR                  << nl;
  Info << "dissolFoamDict, inertia:    " << inertia             << nl;
  Info << "dissolFoamDict, limitFlux:  " << limitFlux           << nl;
  Info << "dissolFoamDict, constFlux:  " << constFlux           << nl;
  Info << "dissolFoamDict, Qlim:       " << Qlim                << nl;
  Info << "dissolFoamDict, writeFluxC: " << writeFluxC          << nl;
  Info << "fvSolutions, tolerance:     " << tolerance           << nl;
  Info << "fvSolutions, maxIter:       " << maxIter             << nl;
  Info << "***************************************************" << endl;

