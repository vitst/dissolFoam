// Momentum predictor. Stokes.

tmp<fvVectorMatrix> UEqn
(
  -  fvm::laplacian(nu, U)
);

UEqn().relax();

solve(UEqn() == -fvc::grad(p));