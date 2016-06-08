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
  
  bool gradCWrite;
  if( !dissolProperties.readIfPresent<bool>("gradCWrite", gradCWrite) ){
    SeriousErrorIn("main")
          << "There is no gradCWrite parameter in dissolFoam dictionary"
          << exit(FatalError);
  }
  
  
  bool inertia;
  if( !dissolProperties.readIfPresent<bool>("inertia", inertia) ){
    SeriousErrorIn("main")
          << "There is no inertia parameter in dissolFoam dictionary"
          << exit(FatalError);
  }
  
  bool rescale;
  if( !dissolProperties.readIfPresent<bool>("rescale", rescale) ){
    SeriousErrorIn("main")
          << "There is no rescale parameter in dissolFoam dictionary"
          << exit(FatalError);
  }
  
  if(!rescale)
  {
    Info<<nl<<"**************************************************************"<<nl;
    Info<<nl<<"!!!! Warning !!!!"<<nl;
    WarningIn("readDicts.H")
            << "The 'rescale' parameter is set to 'false', which means\n"
            "The U and phi fields will not be rescaled. User should specify\n"
            "a pressure difference as 0 at the inlet and A*nu*L at the outlet,\n"
            "where nu is the kinematic viscosity from transportProperties file;\n"
            "L is a length of the fracture and A is a geometric parameter\n"
            "(12 in case of two parallel plates)."<<nl;
    Info<<"**************************************************************"<<nl<<nl;
  }
    
  scalar rescaleRelaxPar;
  if( !dissolProperties.readIfPresent<scalar>("rescaleRelaxPar", rescaleRelaxPar) ){
    SeriousErrorIn("main")
          << "There is no rescaleRelaxPar parameter in dissolFoam dictionary"
          << exit(FatalError);
  }
  

  bool constFlux;
  if( !dissolProperties.readIfPresent<bool>("constFlux", constFlux) ){
    SeriousErrorIn("main")
          << "There is no constFlux parameter in dissolFoam dictionary"
          << exit(FatalError);
  }

  scalar limitValue = 1.0;
  if( !dissolProperties.readIfPresent<scalar>("limitValue", limitValue) ){
    SeriousErrorIn("main")
          << "There is no limitValue parameter in dissolFoam dictionary"
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

  dimensionedScalar nu
  (
      transportProperties.lookup("nu")
  );

  dimensionedScalar D
  (
      transportProperties.lookup("D")
  );

  scalar l_T;                  // l_T=D/(k*h0)
  if( !transportProperties.readIfPresent<scalar>("l_T", l_T) ){
    SeriousErrorIn("main")
        << "There is no l_T parameter in transportProperties dictionary"
        << exit(FatalError);
  }

  Info << "Reading confDiff subdirectory" << nl;

  dictionary convDiff = mesh.solutionDict().subDict("convDiff");

  double convCrit = 0;         // convergence criteria for convection-diffusion
  if( !convDiff.readIfPresent<double>("convergence", convCrit) ){
    SeriousErrorIn("main")
        << "No convergence parameter in convDiff subdictionary"
        << exit(FatalError);
  }
  int maxIter = 0;             // maximum number of iterations
  if( !convDiff.readIfPresent<int>("maxIter", maxIter) ){
    SeriousErrorIn("main")
        << "No maxIter parameter in convDiff subdictionary"
        << exit(FatalError);
  }


  Info << "***************************************************" << nl;
  Info << "transportProperties, nu:    " << nu                  << nl;
  Info << "transportProperties, D:     " << D                   << nl;
  Info << "transportProperties, l_T:   " << l_T                 << nl;
  Info << "dissolFoamDict, inertia:    " << inertia             << nl;
  Info << "dissolFoamDict, rescale:    " << rescale             << nl;
  Info << "dissolFoamDict, scaleRlxr:  " << rescaleRelaxPar     << nl;
  Info << "dissolFoamDict, constFlux:  " << constFlux           << nl;
  Info << "dissolFoamDict, limitValue: " << limitValue          << nl;
  Info << "dissolFoamDict, gradCWrite: " << gradCWrite          << nl;
  Info << "fvSolutions, convCrit:      " << convCrit            << nl;
  Info << "fvSolutions, maxIter:       " << maxIter             << nl;
  Info << "***************************************************" << endl;

