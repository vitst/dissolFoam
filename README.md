This is the `OpenFoam(R)` solver for reactive transport with dissolution.

It solves for steady flow (Stokes or inertial) and reactant transport.

Mesh motion is controlled by the `normalMotionSlip` [libsFoamAux](https://github.com/vitst/libsFoamAux). Parameters to the boundary condition include the name of the field that drives the dissolution and the magnitude of the motion. Examples can be found in `dissolCases`.

Reaction boundary conditions (linear/nonlinear/danckwerts) are coded boundary conditions that are located in the case directory in `constant/bcInclude` (template files in `libsFoamAux/boundaryConditions`) The moving interface in `polyMesh/boundary` should go first.

This version was developed with OpenFOAM-v1912. Check out the latest [releases](https://github.com/vitst/dissolFoam/releases).

Additional documentation about the libraries, solver and cases can be found in the .zip attachment to the each release.

If you use these codes for published work, please cite the following references:

1. V. Starchenko, C. J. Marra and A. J. C. Ladd, Three-dimensional simulations of fracture dissolution,
_J. Geophys. Res. Solid Earth_, __2016__, 121: 6421-6444, [DOI](https://doi.org/10.1002/2016JB013321)
2. V. Starchenko and A. J. C. Ladd, The development of wormholes in laboratory scale fractures: perspectives from three-dimensional simulations, _Water Resource Res._, __2018__, 54: 7946–7959, [DOI](https://doi.org/10.1029/2018WR022948) 
3. F. Dutka et al., Time-dependent shapes of a dissolving mineral grain: Comparisons of simulations with microfluidic experiments, _Chemical Geology_, __2020__, 540: 119459, [DOI](https://doi.org/10.1016/j.chemgeo.2019.119459)
