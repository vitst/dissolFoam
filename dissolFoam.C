/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-201X OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    dissolFoam

Description
    Solves for steady flow (Stokes or inertial) and reactant transport
    normalMotionSlip boundary condition nmoves the mesh according to the
    reactant flux
 
    Works with official and extended versions of OpenFOAM
    Currently: 4.x and v1706    Mar 2 2018

\*---------------------------------------------------------------------------*/

// common and simpleFoam
#include "fvCFD.H"

// OF
#include "pointPatchField.H"
#include "dynamicFvMesh.H"
#include "syncTools.H"
#include "vectorIOList.H"

// dissolFoam project
#include "steadyStateControl.H"
#include "pointPatchField.H"


int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createDynamicFvMesh.H"

  #include "readDicts.H"
  #include "createFields.H"

  const label patchID = mesh.boundaryMesh().findPatchID("outlet");
  
// * * * * *   MAIN LOOP   * * * * * * * * * * * * * * * * * * * * * //

  runTime.functionObjects().execute();     // Execute cntlDict functions
  bool firstCycle = true;                  // Identify first cycle

  while (runTime.run())
  {
    if (firstCycle)  firstCycle = false;   // Skip runTime++ on first cycle
    else             runTime++;
    Info << "Begin cycle: Time = " << runTime.timeName() 
         << "    dt = " << runTime.deltaTValue()
         << nl << endl;
    
/*###############################################
 *   Steady-state flow solver
 *   Control parameters in fvSolution
 *###############################################*/

    steadyStateControl simple(mesh);
    while ( simple.loop() )
    {
      if(inertia)
      {
        #include "UEqn.H"        // Navier Stokes
        #include "pEqn.H"
      }
      else
      {
        #include "UEqnStokes.H"  // Stokes
        #include "pEqn.H"
      }
      
      // rescaling the flow velocities
      if(limitFlux)
      {
        scalar Q = mag( gSum(phi.boundaryField()[patchID]) );
        scalar compareTo = (constFlux) ? SMALL : Qlim;
        scalar scaleFactor = ( Q > compareTo ) ? Qlim / Q : 1.0;
                
        p == scaleFactor * p;
        U == scaleFactor * U;
      }
    }
    
    Info << "Final flow rate: " 
         << mag( gSum(phi.boundaryField()[patchID]) ) << endl;

    Info << "Flow solver: "
         << "ExecutionTime = " << runTime.elapsedCpuTime() << " s "
         << "ClockTime = "<< runTime.elapsedClockTime() << " s"
         << nl << endl;

/*##########################################
 *   Steady-state convection-diffusion solver
 *   Control parameters in fvSolution
 *##########################################*/

    Info << "Steady-state concentration solver"<< endl;

    int  iter = 0;
    while ( iter < maxIter )
    {
      fvScalarMatrix CEqn
      (
        fvm::div(phi, C) - fvm::laplacian(D, C)
      );
      CEqn.relax();
      double residual = solve( CEqn ).initialResidual();

      iter++;
      Info << " Step " << token::TAB << iter << token::TAB
           << " residual: "<< residual << " > " << tolerance << endl;

      if( residual < tolerance )
      {
        Info << "Convection-diffusion: "
             << "ExecutionTime = " << runTime.elapsedCpuTime() << " s "
             << "ClockTime = " << runTime.elapsedClockTime() << " s "
             << nl << "Converged in " << iter << " steps.  Residual = "
             << residual << nl << endl;
        break;                                      // Done
      }

      else if(iter >= maxIter)
      {
          Warning << "dissolFoam Runtime WARNING:"
                  << "Convection-diffusion solver did not converge." << nl
                  << "Maximum number of iterations"
                  << "  iter: "<< iter << exit(FatalError); // No convergence
      }
    }
    
// *********************************************************
// *    Write Output data
// *********************************************************

    if(writeFluxC)
        fluxC == phi/mag(mesh.Sf())*fvc::interpolate(C)
               - D*fvc::snGrad(C);
    
    Info << "Write fields: Time = " << runTime.timeName() << nl << endl;
    runTime.writeNow();

/*###############################################
 *    Mesh motion & relaxation
 *    Control parameters in dynamicMeshDict
 *###############################################*/

    mesh.update();
    Info << "Mesh update: ExecutionTime = " 
         << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl<< endl;
  }

  Info << "End" << endl;
  return 0;
}

// ********************* End of the solver ************************** //
