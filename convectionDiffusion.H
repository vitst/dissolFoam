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
  conv_diff.readIfPresent("convergence", convCritCD);
  conv_diff.readIfPresent("maxIter",     maxNumIterCD);

  // Steady-state convection-diffusion solver main loop
  int counter = 0;
  while ( true ){
    counter++;

    double residual = solve
    (
      fvm::div(phi, C) - fvm::laplacian(D, C) == fvOptions(C)
    ).initialResidual();

    if( residual < convCritCD ){
      Info <<"Convection-diffusion: ExecutionTime = " << runTime.elapsedCpuTime() << " s "
           <<"ClockTime = "<<runTime.elapsedClockTime()<<" s"<<nl 
           <<"Converged in "<<counter<<" steps" 
           <<nl<<nl;
      if(counter >= maxNumIterCD){
        Info<<" dissolFoam Runtime WARNING:"
            <<" Steady-state convection-diffusion solver did not converge."
            <<nl
            <<"It reached maximum number of iterations"
            <<"  counter: "<<counter
            <<"  maxNumIterCD: " <<maxNumIterCD
            <<nl;
      }
      break;
    }
    else{
      Info<<nl<<" Step "<<counter<<" residual: "<<residual<<" > "<<convCritCD<<nl;
    }
  }
