# Code Organization

## simulation
The R programs in this folder are utilized to conduct simulation study of nonhomogeneous Markov-switching penalized mixed-frequency regression.

### DataGen.R
This R script is designed to generate the dataset utilized in the simulation study.

### tool.R
This R script contains utility functions required for the simulation study.

### update.cpp
This Rcpp script contains the MCMC algorithms for Bayesian inference of the NMSP-MIDAS model, including iterative update functions for various parameters as well as auxiliary functions for matrix and vector computations. 

### main.R
This R script serves as the main function for parameter estimation in the NMSP-MIDAS model. The estimation results will be saved in the `Paras_estimation` folder.

### Forecast.R
This R script is used for out-of-sample forecasting on the test dataset. The forecast results for various indicators will be saved in the `Forecast_results` folder.

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

