    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        //fvm::div(phi, U)
      - fvm::laplacian(nu, U)
      //+ turbulence->divDevReff(U)
      ==
        fvOptions(U)
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
