/*    
// Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
//        fvm::div(phi, U)
      - fvm::laplacian(nu, U)
      ==
        fvOptions(U)
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
*/

    tmp<fvVectorMatrix> UEqn
    (
//        fvm::div(phi, U)
      - fvm::laplacian(nu, U)
//      + turbulence->divDevReff(U)
    );

    UEqn().relax();

    eqnResidual = solve(  UEqn() == -fvc::grad(p)  ).initialResidual();
    
    //Info<<nl<< "UUUU residual" << eqnResidual << endl;
    
    // here initialResidual() returns only the biggest value for x,y,z 
    maxResidual = max(eqnResidual, maxResidual);
