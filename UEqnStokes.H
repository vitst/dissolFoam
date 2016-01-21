    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
      -  fvm::laplacian(nu, U)
      == fvOptions(U)
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
