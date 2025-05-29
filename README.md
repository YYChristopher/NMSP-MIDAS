# Bayesian nonhomogeneous Markov-switching penalized mixed-frequency regression

Code to reproduce simulation results, figures, and real data analysis results from the paper "Bayesian nonhomogeneous Markov-switching penalized mixed-frequency regression" by Yue Yu, Xiangnan Feng, Kai Kang, and Xinyuan Song.

## Package need and Package version

Rcpp >= ‘1.0.8.3’

RcppArmadillo >= ‘0.11.0.0.0’

xts >= ‘0.12.1’

fbi >= ‘0.7.0’

gtools >= ‘3.9.2’

data.table >= ‘1.14.2’

dplyr >= ‘1.0.8’

midasr >= ‘0.8’

tseries >= ‘0.10.50’

MASS >= ‘7.3.51.6’

forecast >= ‘8.16’

LaplacesDemon >= ‘16.1.6’

ggplot2 >= ‘3.3.5’

ggExtra >= ‘0.10.1’

latex2exp >= ‘0.9.6’

## Organization

### code
Code of simulation and real analysis of the Paper.

#### realAnalysis
This folder contains the code for empirical analysis, encompassing data preprocessing procedures, a collection of tool functions, Rcpp update scripts, and the implementation of the NMSP-MIDAS and HMSP-MIDAS models.

#### simulation
This folder contains the code for the numerical simulation section, including data generation procedures, a suite of toolfunctions, Rcpp update scripts, the main experimental function for the NMSP-MIDAS model, and the code for out-of-sample forecasting. In addition, the folder also includes manually implemented code for the comparing methods, including estimation and forecasting function files for the HMSP-MIDAS, BMIDAS-AGL, and BMIDAS-AGL-SS methods.

### data
Data utilized into real analysis of the Paper.

### manuscript
The source files (often LaTeX or Rmd) for the manuscript, and introduction of functions utilized.

### output
The figure into the Paper plotting routines and corresponding output results used in the paper.

#### Description
This folder includes some figures plotting code in introduction and online appendix.

#### Simulation
This folder includes out-of-sample forecasting results and figures plotting code of the numerical simulations.

#### RealAnalysis
This folder includes several results and figures plotting code of the real analysis.






