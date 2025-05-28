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
This folder contains the code for the numerical simulation section, including data generation procedures, a suite of toolfunctions, Rcpp update scripts, the main experimental function for the NMSP-MIDAS model, and the code for out-of-sample forecasting.

### data
Data utilized into real analysis of the Paper.

### manuscript
The source files (often LaTeX or Rmd) for the manuscript, and introduction of functions utilized.

### output
The figure into the Paper plotting routines and corresponding output results used in the paper.

#### Simulation
This folder includes out-of-sample forecasting results and figures plotting code of the numerical simulations.

#### RealAnalysis
This folder includes several results figures plotting code of  the real analysis.

## Workflow

### Simulation workflow

First, change the work directory into file 'DataGen.R', 'main.R' and 'Forecast.R'.

Second, running 'DataGen.R' to generate data utilized in simulation. The in-sample training data are stored into folder 'Gen_Data_Train_NS=2', the out-of-sample testing data are stored into folder 'Gen_Data_Test_NS=2', the real parameter information is stored into folder 'Paras_NS=2'.

Then, running 'main.R' to execute the main program file, that is, use NMSP-MIDAS model to estimate the parameters of the in-sample data. The output estimate results are stored into folder 'Paras_estimation'.

Last, running 'Forecast.R' to predict the out-of-sample testing data and store the predicted evaluation results in the folder 'Forecast_results'.

### Real Analysis workflow

First, change the work directory into file 'RealMain.R' and 'RealHMSP.R'.

Second, running 'RealMain.R' to execute the main program file of NMSP-MIDAS, that is, use NMSP-MIDAS model to analysis real data. The (a) subfigure of figure 2 in the Paper which is stored into folder 'Plots'. You can also get the one-step forecasting results of the expanding window(You can store the result by yourself).

Then, running 'RealHMSP.R' to execute the main program file of HMSP-MIDAS, that is, use HMSP-MIDAS model to analysis real data. The (b) subfigure of figure 2 in the Paper which is stored into folder 'Plots'. You can also get the one-step forecasting results of the expanding window(You can store the result by yourself).

### Output workflow

#### Simulation Folder

First, change the work directory into file 'Analysis.R', 'DICPlot.R'(into the folder 'output/Simulation'). 

Then, running 'DICPlot.R' to plot the figure S3 of online appendix. The output figures are stored into folder 'Plots'.

Then, running 'Analysis.R' to plot the figure S4 of online appendix. The output figures are stored into folder 'Plots'.

#### RealAnalysis Folder

First, change the work directory into file 'VSPlot.R', 'OOS_NMSP.R' and 'OOS_HMSP.R'(into the folder 'output/RealAnalysis').

Then, running 'VSPlot.R' to plot the figure S6 and S7 of online appendix. The output figures are stored into folder 'Plots'.

Then, running 'OOS_NMSP.R' to plot the (a) subfigure of figure 3 in the Paper which is stored into folder 'Plots'.

Last, running 'OOS_HMSP.R' to plot the (b) subfigure of figure 3 in the Paper which is stored into folder 'Plots'.




