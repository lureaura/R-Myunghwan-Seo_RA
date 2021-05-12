
##### The brief explanation of files in the dropbox #####

# 1. R files

 @@@ Functions_lapply.R @@@

   : This R file contains individual functions which are used in the simulation code.
     There are 10 functions in this file and there are three groups according to the role of the functions.

     The first group is simulation-related functions and there are 3 functions.
     The second group is random data generating functions and there are 2 functions.
     The third group is test-statistics and their p-value deriving functions and there are 5 functions.

  --------------------------------------------------------------------------------------------------------------
  1) Simulation procedure

    * simulation : This function is the Monte-Carlo simulation code. Set inputs of the function and 
                   then according to the "mod.type"(size or power), it derives the results of table3 
                   or table4 in the paper. That is, it calculates size or power of sup-wald, HW, ADF tests
                   for significant level 0.1 and 0.05 by using simu.size and simu.pow function below.
   
    * simu.size : This function is the Monte-Carlo simuation code for fixed ns(number of simulation),
                  and rho. It generates result of table3 in the paper. This function is included in 
                  simulation function above.

    * simu.pow : This function is the Monte-Carlo simuation code for fixed ns(number of simulation),
                  gamma and alpha. It generates result of table4 in the paper. This function is included in 
                  simulation function above.


  2) Random data generating procedure for Monte-Carlo simulation

    * dgp_linear : This function generates data x according to random bivariate standard normal data.
                   Data generating model is VAR model. It is used in the simu.size function.

    * dgp_band : This function generates data x according to random bivariate standard normal data.
                 Data generating model is band-TVECM model without lag-terms. 
                 It is used in the simu.pow function.

  3) Test-statistic and its p-value calculating procedure

    * boot.w : This function is a residual-based bootstrapping code. It calculates p-value of 
               bootstrapped joint and marginal sup-wald statistics. This function uses est.w and find.w
               functions below.

    * est.w : It is a function which calculates sup-wald statistics and minimum determinant of sigma 
              and each of gamma hat for them respectively for a given data. 
              It finds the argmax(for finding sup-wald) and argmin(for finding minimum value of 
              det. of sigma) gamma1 and gamma2 by iterating according to each pre-set grid. 
              It uses find.w function on its execution.

    * find.w : This function is the code which calculates determinant of sigma hat and 
               joint and marginal sup-wald statistics for a given gamma1, gamma2 and data.

    * hwcoint : It is a function which calculates HW test statistic for a given data.

    * adft : It is a function which calculates ADF test statistic for a given data.

 =================================================================================================================
 =================================================================================================================
 
 @@@ Simulation.R @@@

   : This R file which executes Monte-Carlo simulation. In this file, one can set some main input parameters
     for the simulation function and then get the result of table3 or table4 as explained above Functions_lapply.R file.
     

 =================================================================================================================
 =================================================================================================================
 
 @@@ Data_based result.R @@@

   : This R file has 3 individual functions which are "boot.w.detail", "est.w", "find.w" codes above.
     est.w and find.w functions are the same as in Function_lapply.R file. The difference is the boot.w.detail function.
     
     This function is the same function with boot.w function of Functions_lapply.R but it is different 
     on its output. This function gives more detail result in the bootstrapping procedure. While the boot.w
     function gives only p-value of bootstrap-based p-value of sup-wald, this boot.w.detail function gives
     sup-wald statistic value, gamma estimation values, coefficient estimation values, p-value, 
     minimum number of observation in each regime(m in the paper of 4-th chapter), and grid number for iterating.

 =================================================================================================================
 =================================================================================================================
 
 @@@ Se46.R @@@

   : This file is the same as Data_based result.R but this file has a given input parameter and data for
     the function boot.w.detail function. The data is from se46.csv. 