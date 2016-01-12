/* 
 * File:   convectionDiffusion.H
 * Author: vstar
 *
 * Created on October 14, 2015, 4:33 PM
 */

  /*############################################
  *   Steady-state convection-diffusion solver
  * ##########################################*/
  // define dictionary for convection diffusion solver
  dictionary conv_diff = mesh.solutionDict().subDict("CONVECTION_DIFFUSION");
  double convCritCD = 0;    // convergence criteria
  int maxNumIterCD = 0;     // maximum number of iterations
  // read values for convergence checks from the dictionary
  if( !conv_diff.readIfPresent<double>("convergence", convCritCD) ){
    SeriousErrorIn("main")
            <<"There is no convergence parameter in CONVECTION_DIFFUSION subdictionary"
            <<exit(FatalError);
  }
  if( !conv_diff.readIfPresent<int>("maxIter",     maxNumIterCD) ){
    SeriousErrorIn("main")
            <<"There is no maxIter parameter in CONVECTION_DIFFUSION subdictionary"
            <<exit(FatalError);
  }

  // Steady-state convection-diffusion solver main loop
  int counter = 0;
  while ( true ){
    counter++;

    double residual = solve
    (
      fvm::div(phi, C) - fvm::laplacian(D, C) == fvOptions(C)
    ).initialResidual();

    if( residual < convCritCD ){
      Info << "Convection-diffusion: "
           << "ExecutionTime = " << runTime.elapsedCpuTime() << " s "
           << "ClockTime = " << runTime.elapsedClockTime() << " s" <<nl 
           << "Converged in " << counter << " steps" << nl << endl;

      if(counter >= maxNumIterCD){
        Info << nl << "dissolFoam Runtime WARNING:"
             << "Convection-diffusion solver did not converge." << nl
             <<"Maximum number of iterations"
             <<"  counter: "<< counter << endl;
      }
      break;
    }
    else{
      Info << " Step " << counter
           << " residual: "<< residual << " > " << convCritCD << endl;
    }
  }
