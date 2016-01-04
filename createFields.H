    Info<< "Reading transportProperties\n" << nl;

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );
    dimensionedScalar D
    (
        transportProperties.lookup("D")
    );

    
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
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue    );
    
    
    /*
    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );
    */
    
    

