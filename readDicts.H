  // reading dissolFoam dictionary

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
  
  bool inertia;
  if( !dissolProperties.readIfPresent<bool>("inertia", inertia) ){
    SeriousErrorIn("main")
          << "There is no inertia parameter in dissolFoam dictionary"
          << exit(FatalError);
  }

  bool constFlux;
  if( !dissolProperties.readIfPresent<bool>("constFlux", constFlux) ){
    SeriousErrorIn("main")
          << "There is no constFlux parameter in dissolFoam dictionary"
          << exit(FatalError);
  }

  bool gradCWrite;
  if( !dissolProperties.readIfPresent<bool>("gradCWrite", gradCWrite) ){
    SeriousErrorIn("main")
          << "There is no gradCWrite parameter in dissolFoam dictionary"
          << exit(FatalError);
  }
  
  bool limitFlux = false;
  scalar limitValue = 1.0;
  if( !constFlux ){
    if( !dissolProperties.readIfPresent<bool>("limitFlux", limitFlux) ){
      SeriousErrorIn("main")
            << "There is no limitFlux parameter in dissolFoam dictionary"
            << exit(FatalError);
    }
    
    if( !dissolProperties.readIfPresent<scalar>("limitValue", limitValue) ){
      SeriousErrorIn("main")
            << "There is no limitValue parameter in dissolFoam dictionary"
            << exit(FatalError);
    }
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

  dimensionedScalar nu
  (
      transportProperties.lookup("nu")
  );

  dimensionedScalar D
  (
      transportProperties.lookup("D")
  );

  scalar l_T;                                        // l_T=D/(k*l)
  if( !transportProperties.readIfPresent<scalar>("l_T", l_T) ){
    SeriousErrorIn("main")
        << "There is no l_T parameter in transportProperties dictionary"
        << exit(FatalError);
  }

  Info << "Reading confDiff subdirectory" << nl;

  dictionary convDiff = mesh.solutionDict().subDict("convDiff");

  double convCrit = 0;              // convergence criteria
  if( !convDiff.readIfPresent<double>("convergence", convCrit) ){
    SeriousErrorIn("main")
        << "No convergence parameter in convDiff subdictionary"
        << exit(FatalError);
  }
  int maxIter = 0;                  // maximum number of iterations
  if( !convDiff.readIfPresent<int>("maxIter", maxIter) ){
    SeriousErrorIn("main")
        << "No maxIter parameter in convDiff subdictionary"
        << exit(FatalError);
  }


  Info << "***************************************************" << nl;
  Info << "transportProperties, nu:    " << nu         << nl;
  Info << "transportProperties, D:     " << D          << nl;
  Info << "transportProperties, l_T:   " << l_T        << nl;
  Info << "dissolFoamDict, inertia:    " << inertia    << nl;
  Info << "dissolFoamDict, constFlux:  " << constFlux  << nl;
  Info << "dissolFoamDict, gradCWrite: " << gradCWrite << nl;
  Info << "dissolFoamDict, limitFlux:  " << limitFlux  << nl;
  Info << "dissolFoamDict, limitValue: " << limitValue << nl;
  Info << "fvSolutions, convCrit:      " << convCrit   << nl;
  Info << "fvSolutions, maxIter:       " << maxIter    << nl;
  Info << "***************************************************" << endl;

