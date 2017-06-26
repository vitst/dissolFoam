// Momentum predictor. Stokes.

tmp<fvVectorMatrix> tUEqn
(
     fvm::Sp(Kinv*nu, U)
  -  fvm::laplacian(nu, U)
);

fvVectorMatrix& UEqn = tUEqn.ref();
UEqn.relax();

solve(UEqn == -fvc::grad(p));
