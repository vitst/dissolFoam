This is the OpenFoam solver.

It solves for steady flow (Stokes or inertial) and reactant transport

Mesh motion is controlled by the normalMotionSlip condition (libsFoamAux)
Parameters to the boundary condition include the name of the field
that drives the dissolution and the magnitude of the motion. Examples
can be found in dissolFoamCases.

Reaction boundary conditions (linear/nonlinear/danckwerts) are coded
boundary conditions that are located in the case directory in
constant/bcInclude (template files in libsFoamAux/boundaryConditions)
The moving interface in polyMesh/boundary should go first.

This version is known to work with OpenFOAM-v1706. It should also work
with the official versions 4.x, 5.x and the v1712 extended release.

Additional documentation about the solver and associated codes /
cases can be found in the file doc.pdf

