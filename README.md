This is the OpenFoam solver for reactive transport with dissolution.

It solves for steady flow (Stokes or inertial) and reactant transport

Mesh motion is controlled by the normalMotionSlip (libsFoamAux) Parameters to the boundary condition include the name of the field that drives the dissolution and the magnitude of the motion. Examples can be found in dissolCases.

Reaction boundary conditions (linear/nonlinear/danckwerts) are coded boundary conditions that are located in the case directory in constant/bcInclude (template files in libsFoamAux/boundaryConditions) The moving interface in polyMesh/boundary should go first.

This version was developed with OpenFOAM-v1706. Check out the latest [releases](https://github.com/vitst/dissolFoam/releases) compartible with [OpenFOAM-v1712](openfoam.com) and [OpenFOAM-6](openfoam.org):

Additional documentation about the libraries, solver and cases can be found in the file doc.pdf in the .zip attachment to the each release.

If you use these codes for published work, please cite the following references:

[1] V. Starchenko, C. J. Marra and A. J. C. Ladd, Three-dimensional simulations of fracture dissolution,
J. Geophys. Res. Solid Earth, 121:6421-6444, 2016.
<br>
[2] V. Starchenko and A. J. C. Ladd, The development of wormholes in laboratory scale fractures: perspectives from three-dimensional simulations, Water Resource Res., Submitted, 2018.
