// Momentum predictor
    tmp<fvVectorMatrix> UEqn
    (
//        fvm::div(phi, U)
      - fvm::laplacian(nu, U)
//      + turbulence->divDevReff(U)
    );

    UEqn().relax();

    eqnResidual = solve(  UEqn() == -fvc::grad(p)  ).initialResidual();
    
    // here initialResidual() returns only the biggest value for x,y,z 
    maxResidual = max(eqnResidual, maxResidual);
