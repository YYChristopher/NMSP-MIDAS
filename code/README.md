# Code Organization

## simulation
The R programs in this folder are utilized to conduct simulation study of nonhomogeneous Markov-switching penalized mixed-frequency regression.

### DataGen.R
This R script is designed to generate the dataset utilized in the simulation study.

### tool.R
This R script contains utility functions required for the simulation study.

### update.cpp
This Rcpp script contains the MCMC algorithms for Bayesian inference of the NMSP-MIDAS model, including iterative update functions for various parameters as well as auxiliary functions for matrix and vector computations. 

### updateHMSP.cpp
This Rcpp script contains the MCMC algorithms for Bayesian inference of the HMSP-MIDAS model, including iterative update functions for various parameters as well as auxiliary functions for matrix and vector computations. 

### main.R
This R script serves as the main function for parameter estimation in the NMSP-MIDAS model. The estimation results will be saved in the `Paras_estimation` folder.

### mainHMSP.R
This R script serves as the main function for parameter estimation in the HMSP-MIDAS model. The estimation results will be saved in the `Paras_estimation_HMM` folder.

### Forecast.R
This R script is used for out-of-sample forecasting on the test dataset by NMSP-MIDAS model. The forecast results for various indicators will be saved in the `Forecast_results` folder.

### ForecastHMSP.R
This R script is used for out-of-sample forecasting on the test dataset by HMSP-MIDAS model. The forecast results for various indicators will be saved in the `Forecast_results_HMM` folder.

## realAnalysis
The R programs in this folder are utilized to conduct empirical analysis of American logarithmic GDP growth rate. 

### tool.R
This R script contains utility functions required for the empirical analysis.

### update.cpp
This Rcpp script implements MCMC algorithms for Bayesian estimation of the NMSP-MIDAS model into the empirical analysis, comprising parameter update routines and utility functions for matrix and vector operations.

### updateHMSP.cpp
This Rcpp script implements MCMC algorithms for Bayesian estimation of the HMSP-MIDAS model into the empirical analysis, comprising parameter update routines and utility functions for matrix and vector operations.

### RealMain.R
This R script serves as the main function for parameter estimation of the NMSP-MIDAS model into the empirical analysis. One-step forecasting with an expanding window is also implemented in this R file. In addition, the script can plot the unobservable states estimated by the NMSP-MIDAS model when h = 0.

### RealHMSP.R
This R script serves as the main function for parameter estimation of the HMSP-MIDAS model into the empirical analysis. One-step forecasting with an expanding window is also implemented in this R file. In addition, the script can plot the unobservable states estimated by the HMSP-MIDAS model when h = 0.

### BMIDAS.R
This file contains the functions for implementing the BMIDAS-AGL-SS and BMIDAS-AGL algorithms.

### BMIDAStrial.R
This function serves as the main function for conducting numerical simulation experiments using the BMIDAS method. It includes both the BMIDAS-AGL and BMIDAS-AGL-SS approaches, which can be manually selected by modifying the function. The resulting out-of-sample forecast evaluation results can be found in the corresponding folder 'Result_BMIDAS'.

# Workflow

## Simulation workflow

First, change the work directory into file 'DataGen.R', 'main.R' and 'Forecast.R'.

Second, running 'DataGen.R' to generate data utilized in simulation. The in-sample training data are stored into folder 'Gen_Data_Train_NS=2', the out-of-sample testing data are stored into folder 'Gen_Data_Test_NS=2', the real parameter information is stored into folder 'Paras_NS=2'.

Then, running 'main.R' to execute the main program file of NMSP-MIDAS, that is, use NMSP-MIDAS model to estimate the parameters of the in-sample data. The output estimate results are stored into folder 'Paras_estimation'. Running 'mainHMSP.R' to execute the main program file of HMSP-MIDAS, that is, use HMSP-MIDAS model to estimate the parameters of the in-sample data. The output estimate results are stored into folder 'Paras_estimation_HMM'. Running 'BMIDAStrial.R' to execute the main program file of BMIDAS-AGL and BMIDAS-AGL-SS, use these two models to estimate the parameters of the in-sample data. You can select two different functions into the file to use different methods. You can directly get the forecast results of out-of-sample testing data.

Last, running 'Forecast.R' to predict the out-of-sample testing data by NMSP-MIDAS and store the predicted evaluation results in the folder 'Forecast_results'. Running 'ForecastHMSP.R' to predict the out-of-sample testing data by HMSP-MIDAS and store the predicted evaluation results in the folder 'Forecast_results_HMM'. The predicted evaluation results of BMIDAS-AGL-SS and BMIDAS-AGL are in the folder 'Result_BMIDAS'.

For other comparing methods, we primarily rely on existing R packages for implementation. For method 'LASSO-MIDAS', 'MCP-MIDAS', 'SCAD-MIDAS', we use the classic penalized regression R package `ncvreg` for implementation. For method 'EN-MIDAS', we use the use the penalized regression R package `glmnet`. For main comparing model 'SG-LASSO-MIDAS', we use the R package `midasml`. For method 'BMA-MIDAS', we use the R package `BMA`. Last, for method 'factor-MIDAS', we use R package `psych` to reduce the dimensionality of high-frequency explanatory variables by 'PCA' method before performing regression.

## Real Analysis workflow

First, change the work directory into file 'RealMain.R' and 'RealHMSP.R'.

Second, running 'RealMain.R' to execute the main program file of NMSP-MIDAS, that is, use NMSP-MIDAS model to analysis real data. The (a) subfigure of figure 2 in the Paper which is stored into folder 'Plots'. You can also get the one-step forecasting results of the expanding window(You can store the result by yourself).

Then, running 'RealHMSP.R' to execute the main program file of HMSP-MIDAS, that is, use HMSP-MIDAS model to analysis real data. The (b) subfigure of figure 2 in the Paper which is stored into folder 'Plots'. You can also get the one-step forecasting results of the expanding window(You can store the result by yourself).

The main functions for the BMIDAS-AGL and BMIDAS-AGL-SS methods are provided in the folder `simulation`(BMIDAS.R). For the benchmark method AR(1), we mainly rely on the classic time series package `tseries` in R. The packages used for other comparison methods are consistent with those introduced in the simulation workflow section. 

## The utilized packages version:

ncvreg: >= 3.13.0; glmnet: >= 4.1.4; midasml: >= 0.1.9.1; BMA: >= 3.18.15; psych: >= 2.4.6.26; tseries: >= 0.10.50.

