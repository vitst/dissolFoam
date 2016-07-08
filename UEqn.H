    // Momentum predictor

    //MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
     -  fvm::laplacian(nu, U)
    
//      + MRF.DDt(U)
//      + turbulence->divDevReff(U)
//     ==
//        fvOptions(U)
    );

    fvVectorMatrix& UEqn = tUEqn.ref();
		UEqn.relax();

//    fvOptions.constrain(UEqn());

    solve(UEqn == -fvc::grad(p));

//    fvOptions.correct(U);
