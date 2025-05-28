# rm(list=ls())
# library(BVAR)
library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(tseries)
library(ggplot2)
library(ggExtra)
library(latex2exp)
library(xts)
library(fbi)
library(dplyr)
library(midasr)
library(forecast)
library(data.table)
library(gtools)
library(latex2exp)

# enlarge memory for R
memory.limit(size = 200000)
#Work Space of the data
work_dir = 'F:/24561/Documents/ResearchWrite/NHMM-MIDAS/Submission Requirements/Code in Github'
#Set the Work Space of the data
setwd(work_dir)
source('./code/realAnalysis/tool.R')
sourceCpp("./code/realAnalysis/updateHMSP.cpp")

#################################### Data Process ########################################
# GDP Data to transform into GDP growth rate
# truncate data that we need
GDP = read.csv('./data/GDP2.csv', header = TRUE)
GDP_noindex <- GDP['GDP']
dates_gdp <- seq(as.Date("1947-03-01"), by = "quarter", length.out = 308)
GDP_xts <- xts(GDP_noindex, order.by = dates_gdp)
gdp <- GDP_xts['1959-03-01/2023-12-01']
gdp.train <- GDP_xts['1959-03-01/2011-03-01']
# log transform
len = as.numeric(length(gdp.train))
gdp_rate = rep(0, len - 1)
gdp_rate2 = rep(0, len - 1)
for (i in 1:(len - 1)) {
  gdp_rate[i] = 100 * log(as.numeric(gdp[i + 1]) / as.numeric(gdp[i]))
  gdp_rate2[i] = 100 * (log(as.numeric(gdp[i + 1])) / log(as.numeric(gdp[i])) - 1)
}
growth_dates = dates <- seq(as.Date("1959-03-01"), by = "quarter", length.out = len - 1)
gdp_growth = xts(gdp_rate, order.by = growth_dates)
gdp_growth_num = as.numeric(gdp_growth) 
gdp_growth.train <- gdp_growth["1958-12-01/2010-12-01"]
# hist(x=gdp_growth_num, freq = FALSE, breaks = 100, xlab = 'log US GDP', main = 'Histogram of log US GDP')
# hist(x=gdp_growth_num, freq = FALSE, breaks = 100, xlim = c(-0.1, 0.3), xlab = 'log US GDP growth rate', main = 'Histogram of log US GDP growth rate')
# lines(density(gdp_growth_num))

# Load Data and truncate data that we need
# FRED-MD Data (part of explanable variables, monthly)
filepath <- './data/fredmd.csv'
fred_md <- fredmd(filepath, date_start = NULL, date_end = NULL, transform = TRUE)
pm <- ncol(fred_md)
dates_fredmd <- seq(as.Date("1959-01-01"), by = "month", length.out = 780)
fredmd_xts <- xts(fred_md, order.by = dates_fredmd)[ ,-1]
fredmd <- as.data.frame(fredmd_xts['1959-01-01/2023-12-01'])
# fredmd_mat <- as.matrix(fredmd)
fredmd[] <- lapply(fredmd, function(x) as.numeric(as.character(x)))
# str(fredmd)

# divide month date into train set and test set
NT_month.train <- 624
NT_month.test <- length(dates_fredmd) - NT_month.train
dates_fredmd.train <- dates_fredmd[1:NT_month.train]
dates_fredmd.test <- dates_fredmd[(NT_month.train + 1):length(dates_fredmd)]
fredmd.train <- fredmd[1:NT_month.train, ]

# FF3 (Fama-French 3 factors model variables) monthly
FF3month <- as.data.frame(read.csv('./data/FF3month/ff3monthly.csv', header = TRUE))
ff3.MKTRF <- FF3month$Mkt.RF
ff3.SMB <- FF3month$SMB
ff3.HML <- FF3month$HML
ff3.RF <- FF3month$RF
dates_ff3 <- as.Date(paste(FF3month[ ,1], "01"), format = '%Y%m%d')
ff3MKTRF_xts <- xts(ff3.MKTRF, order.by = dates_ff3)
ff3SMB_xts <- xts(ff3.SMB, order.by = dates_ff3)
ff3HML_xts <- xts(ff3.HML, order.by = dates_ff3)
ff3RF_xts <- xts(ff3.RF, order.by = dates_ff3)
ff3MKTRF <- ff3MKTRF_xts['1959-01-01/2023-12-01']
ff3SMB <- ff3SMB_xts['1959-01-01/2023-12-01']
ff3HML <- ff3HML_xts['1959-01-01/2023-12-01']
ff3RF <- ff3RF_xts['1959-01-01/2023-12-01']
ff3 <- cbind(ff3MKTRF, ff3SMB, ff3HML, ff3RF)
ff3.train <- ff3[1:NT_month.train, ]

# ADS (Aruoba-Diebold-Scotti daily business conditions index) daily
ADS <- as.data.frame(read.csv('./data/ADS.csv', header = TRUE))
ads_noindex <- ADS['ADS_Index']
dates_ads <- as.Date(ADS[ ,1])
ads_xts <- xts(ads_noindex, order.by = dates_ads)
ads <- ads_xts['1960-03-01/2023-12-01']

# Federal FUNDS Data (monthly)
Fedfunds <- read.csv('./data/FEDFUNDS-month.csv', header = TRUE)
fedfunds_noindex <- Fedfunds['FEDFUNDS']
dates_fed <- seq(as.Date("1954-07-01"), by = "month", length.out = 835)
fedfunds_xts <- xts(fedfunds_noindex, order.by = dates_fed)
fedfunds <- fedfunds_xts['1958-12-01/2023-12-01']
diff_forward <- diff(fedfunds, lag = 1)
fedfunds <- diff_forward[-1]
# Data Process (NA filling)
# sum(is.na(ff3RF))
# sum(is.na(ff3MKTRF))
# sum(is.na(ff3HML))
# sum(is.na(ff3SMB))
# sum(is.na(ncfiInd))
# sum(is.na(ncfiRisk))
# sum(is.na(ncfiCredit))
# sum(is.na(ncfiLeverage))
# sum(is.na(ads))
# sum(is.na(tyff))
# sum(is.na(fedfunds))
# sum(is.na(gdp))
sum(is.na(fredmd.train))  # 387 NA need to be filled.
has_na <- colSums(is.na(fredmd.train)) > 0
na_columns <- names(fredmd.train)[has_na]
na_columns
na_counts <- colSums(is.na(fredmd.train[has_na]))
na_counts

nac1 <- as.numeric(fredmd.train$ACOGNO)
nac2 <- as.numeric(fredmd.train$TWEXAFEGSMTHx)
nac3 <- as.numeric(fredmd.train$UMCSENTx)
nac1_omit = na.omit(nac1)
nac2_omit = na.omit(nac2)
nac3_omit = na.omit(nac3)
# hist(x=nac1_omit, breaks = 100)
# lines(density(nac1_omit))
# hist(x=nac2_omit, breaks = 100)
# lines(density(nac2_omit))
# hist(x=nac3_omit, breaks = 100)
# lines(density(nac3_omit))
nac1_mean = mean(nac1_omit)
nac1_sd = sd(nac1_omit)
nac2_mean = mean(nac2_omit)
nac2_sd = sd(nac2_omit)
nac3_mean = mean(nac3_omit)
nac3_sd = sd(nac3_omit)

nac1_count <- sum(is.na(fredmd.train$ACOGNO))  
fill_values1 <- rnorm(nac1_count, mean = nac1_mean, sd = nac1_sd)
fredmd.train$ACOGNO[is.na(fredmd.train$ACOGNO)] <- fill_values1
nac2_count <- sum(is.na(fredmd.train$TWEXAFEGSMTHx)) 
fill_values2 <- rnorm(nac2_count, mean = nac2_mean, sd = nac2_sd)
fredmd.train$TWEXAFEGSMTHx[is.na(fredmd.train$TWEXAFEGSMTHx)] <- fill_values2
nac3_count <- sum(is.na(fredmd.train$UMCSENTx))  
fill_values3 <- rnorm(nac3_count, mean = nac3_mean, sd = nac3_sd)
fredmd.train$UMCSENTx[is.na(fredmd.train$UMCSENTx)] <- fill_values3

fredmd.train$`S&P div yield` <- as.numeric(fredmd.train$`S&P div yield`)
fredmd.train$CP3Mx <- as.numeric(fredmd.train$CP3Mx)
for(i in 1:ncol(fredmd.train)) {
  fredmd.train[ , i][is.na(fredmd.train[ , i])] <- median(fredmd.train[ , i], na.rm=TRUE)
}
fredmd.train <- as.data.frame(fredmd.train)

# fredmd_filled <- fredmd %>%
# mutate(across(where(is.numeric), ~if_else(is.na(.), median(., na.rm = TRUE), .)))
has_na <- colSums(is.na(fredmd.train)) > 0
na_columns <- names(fredmd.train)[has_na]
na_columns
na_counts <- colSums(is.na(fredmd.train[has_na]))
na_counts

# summary descriptive statistics
fredmd.summ <- sapply(fredmd.train, function(x) c(Mean = mean(x, na.rm = TRUE),
                                                  SD = sd(x, na.rm = TRUE),
                                                  Min = min(x, na.rm = TRUE),
                                                  Max = max(x, na.rm = TRUE)))


# Data Process (Transform to MIDAS Form)
# Constructing X matrices in the MIDAS form.
# extract data that we need (slope of the yield)
# dates_slope <- seq(as.Date("1960-01-01"), by = "month", length.out = 612)
# slope_yield <- xts(x = fredmd$GS10 - fredmd$TB3MS, order.by = dates_slope)
fredmd.train <- subset(fredmd.train, select = -c(VIXCLSx))
all.names <- c(names(fredmd.train))
fredmd.train <- subset(fredmd.train, select = -c(FEDFUNDS, CP3Mx, TB3MS, TB6MS, GS1, GS5, GS10, AAA, BAA))
fredmd.train[, "HWI"] <- fredmd.train[, "HWI"] / 100
fredmd.train[, c("CES0600000007", "AWHMAN")] <- fredmd.train[, c("CES0600000007", "AWHMAN")] / 10
ads.train.org <- ads_xts['1960-03-01/2010-12-01']
Xm.train.org <- as.data.frame(cbind(fredmd.train, ff3.train))
Xm.train <- as.matrix(scale(cbind(fredmd.train, ff3.train)))
Xm.mean <- colMeans(Xm.train.org)
Xm.std <- c(rep(0, ncol(Xm.train.org)))
for (m in 1:ncol(Xm.train.org)) {
  Xm.std[m] <- sd(Xm.train.org[ ,m])
}
ads.train <- scale(ads.train.org)
ads.mean <- mean(ads.train.org)
ads.std <- sd(ads.train.org)

# Utilizing scale
# fredmd.scale <- scale(fredmd)
# Xm <- as.matrix(fredmd.scale)
Nbvar_month = as.numeric(ncol(Xm.train))
Nbvar <- Nbvar_month
Spc = vector('list', 6)
names(Spc) = c('mseries', 'monthly', 'Km', 'kappam', 'sK', 'nbvar')
Spc$mseries = NULL
Spc$monthly = Nbvar
Spc$Km = 12 # lag orders. 
Spc$kappam = 3 
Spc$sK = rep(Spc$Km, Spc$monthly)
Spc$nbvar = Spc$monthly
K = Spc$Km

NT_gdp.train = as.numeric(nrow(gdp_growth.train))
K_mon = Spc$Km 
m_mon = as.numeric(nrow(fredmd.train) %/% nrow(gdp_growth.train))
# Monthly matrix
xmls.train <- matrix(0, nrow = NT_gdp.train, ncol = 0)
xx <- vector("list", length = Nbvar) # monthly is Nbvar
for (ii in 1:Nbvar) {
  xx[[ii]] <- Construct_DataMIDAS(g = gdp_growth.train, d = Xm.train[, ii], K = K_mon, m = m_mon)
  xmls.train <- cbind(xmls.train, xx[[ii]])
}
xd2 <- vector("list", length = 1)
m_day2 = as.numeric(length(ads.train) %/% nrow(gdp_growth)) # 89
K_day2 = m_day2
xd2[[1]] <- Construct_DataMIDAS(g = gdp_growth, d = ads.train, K = K_day2, m = m_day2)
xmls.train <- cbind(xmls.train, xd2[[1]][1:NT_gdp.train, ]) # dim(XX_all): 208 601
notuseindex <- 4
xmls.train <- xmls.train[-(1:notuseindex), ]


NT <- NT_train <- dim(xmls.train)[1]
NT_all <- length(gdp_growth) - notuseindex
Xused = xmls.train[1:NT_train, ]
y.train = gdp_growth.train[-(1:notuseindex)]
x_train = Xused
y_train = y.train
real_date <- growth_dates[-(1:notuseindex)]
real_date_train <- real_date[1:NT_train]

# ------------------------------------- Algorithm ------------------------------------------
# Algorithm basic settings
# set.seed(43)
Rep = 1
Y <- as.vector(y_train)
Y.pre <- as.vector(c(0, Y[-length(Y)]))
# Y.front <- as.vector(Y[-length(Y)])
# Y.back <- as.vector(Y[-1])
X <- x_train
N = 1
NS = 2

polydegree = 2# Almon lag degree
Ra = vector('list',2)
names(Ra) = c('fC','dfC')
Ra$fC = 1 # no tail restriction; if has set 1
Ra$dfC = 0 # no deravative restriction; if has set 1
Q = Almon_lag(polydegree = polydegree, C = Spc$Km, R = Ra)
Qrow = as.numeric(dim(Q)[1])
eg <- Qrow
NX_alm = Qrow * Nbvar
Qd = Almon_lag(polydegree = polydegree, C = 89, R = Ra)
Qdrow = as.numeric(dim(Qd)[1])
egd <- Qdrow
kdm <- as.numeric(dim(xd2[[1]])[2])
# QX = matrix(0, nrow = nrow(X), ncol = Nbvar * eg)
QX = matrix(0, nrow = nrow(X), ncol = Nbvar * eg + egd)
for (nb in 1:Nbvar) {
  QX[ ,((nb - 1)*eg + 1):(nb*eg)] = X[ ,((nb - 1)*K_mon + 1):(nb*K_mon)] %*% t(Q)
}
QX[ ,(Nbvar*eg + 1):(Nbvar*eg + egd)] =  X[ ,(Nbvar*K_mon+1):(Nbvar*K_mon + kdm)] %*% t(Qd)

set.seed(as.numeric(Sys.time()) + sample(1:10000, 1, replace = FALSE))
# some basic information for simulated data
p = eg * Nbvar + egd
NF <- p
ng <- Nbvar + 1  # number of groups. nrow(Q) is the polydegree + 1(contain 0).
nj <- c(rep(nrow(Q), Nbvar), egd)  # number of predictors in each group
groups <- matrix(0, nrow = ng, ncol = nrow(Q)) 
groups[1, ] <- seq_len(nj[1])
if (ng > 1){
  for (i in 2:ng) {
    groups[i, ] <- seq(sum(nj[1:(i - 1)]) + 1, sum(nj[1:(i - 1)]) + nj[i])
  }
}

# Standalization for data. 
x = array(0, dim = c(N, NT, NF))
y = array(0, dim = c(NT, N))
for (n in 1:N) {
  # result_x <- standardize_me2(QX)
  result_x <- standardize_me(QX)
  # result_y <- center_me(as.matrix(Y))
  result_y <- as.matrix(Y)
  x_scale = result_x$y
  mu_x = result_x$mu
  sig_x = result_x$sig
  y_scale = result_y
  x[n, , ] <- x_scale
  y[ ,n] <- y_scale
  # x[n, , ] <- QX
  # y[ ,n] <- Y
}

# sig_x for each hidden state
sig_x_NS = matrix(0, nrow = NS, ncol = ng * length(groups[1, ]))

sg <- numeric(ng) # sg stands for the dimension of each parameter theta
for (ig in 1:ng) {
  sg[ig] <- length(groups[ig, ])
}

# First Initialization
lambda2 <- matrix(rep(1, ng * NS), nrow = NS, ncol = ng) # lambda^2
intercept <- runif(NS, min = -1, max = 1)
sig2 <- runif(NS, min = 0, max = 1) # sigma^2
tau2 <- matrix(runif(NS * ng), nrow = NS, ncol = ng) # tau^2
beta <- matrix(0, nrow = NS, ncol = ng * sg[1]) 
betan <- matrix(0, nrow = NS, ncol = ng)
kappabar <- rep((1 + 1 / ng), NS)
u <- kappabar
aa <- kappabar * (ng^u) # stands for c in paper
bb <- rep(1, NS) # stands for d in paper
pi0 <- rbeta(NS, shape1 = aa, shape2 = bb, ncp = 0)
pi1 <- matrix(runif(NS * ng), nrow = NS, ncol = ng)
Z <- matrix(0, nrow = NS, ncol = ng)

x_NT=c()
for(i in 1:NT){
  x_NT=rbind(x_NT,x[,i,])
}
y_NT=c()
for(i in 1:NT){
  y_NT=c(y_NT,y[i, ])
}

fc=array(NA,dim=c(N,NT,Nbvar*eg+egd+1))
fc_NT=array(NA,dim=c(N*NT,Nbvar*eg+egd+1))
if(Nbvar!=0){
  fc_NT=cbind(x_NT[,1:((Nbvar+1)*eg)],rep(1,N*NT))
  for(i in 1:NT){
    fc[,i,]=fc_NT[(1+(i-1)*N):(i*N),]
  }
}else{
  fc_NT[,1]=rep(1,N*NT)
  for(i in 1:NT){
    fc[,i,]=fc_NT[(1+(i-1)*N):(i*N)]
  }
}

# tuning parameters for MH algorithm
c_lop = 15
c_zeta = 12
c_phi = 3
cal = rep(0, N)

## initial value of parameters
alpha_cpp = runif(NS * (eg * (Nbvar + 1) + 1), -1, 1)
lop_cpp = runif(NS - 1, -1, 1)
betan_cpp = as.vector(betan)
lambda2_cpp = as.vector(lambda2)
tau2_cpp = as.vector(tau2)
pi0_cpp = as.vector(pi0)
pi1_cpp = as.vector(pi1)
Z_cpp = as.vector(Z)
groups = as.vector(t(groups))
sigma_cpp = runif(NS, 0, 1)
s_cpp = sample(0:(NS - 1), N * NT, replace = TRUE)
y_cpp = y_NT
x_cpp = as.vector(x_NT)
fc_cpp = as.vector(fc_NT)
priorbeta = c(0, 0, 0)
ele = c(40, 40)

#iterations
iter = 1   ## number of replication
rep = 1    ## current replication
It = vector('list',3)
names(It) = c('nsave','nburn','nthin')
It$nsave = 10000 # total sample
It$nburn = 10000 # burn sample
It$nthin = 1
ntotal = It$nsave + It$nburn
ndraw = It$nsave / It$nthin 
# Update and Iter times for Gibbs Sampling
nsave <- It$nsave
nburn <- It$nburn
nthin <- It$nthin 
ntot <- nsave + nburn
ndraw <- nsave / nthin

seed1 = as.numeric(Sys.time()) + sample(1:10000, 1, replace = FALSE)
# seed1 = 1730975314
set.seed(seed1)
start.time = Sys.time()
result = mcmc_hmm(y_cpp, x_cpp, fc_cpp, alpha_cpp,
              sigma_cpp, betan_cpp, lambda2_cpp, tau2_cpp, pi0_cpp, 
              pi1_cpp, Z_cpp, sg, groups, s_cpp, priorbeta, 
              N, NT, NS, NF, iter, rep, ntot, ndraw, ng, eg, ele)
end.time = Sys.time()

# variables <- ls()
# sizes <- sapply(variables, function(x) object.size(get(x)))
# df <- data.frame(Variable = variables, Size = sizes / (1024 * 1024)) # 转换�? MB
# print(sum(df$Size))
# df <- df[order(df$Size), ]
# View(df)

## collect the estimates of MCMC iteration
sigma_result = array(result$sigma_str, dim = c(NS, ndraw))
alpha_result = array(0, dim = c(ndraw, NS, NF + 1))
for(i in 1:NS){
  for(j in 1:(NF + 1)){
    alpha_result[ ,i ,j] = result$alpha_str[(0:(ndraw - 1)) * NS * (NF + 1) +(i - 1) * (NF + 1) + j]
  }
}
p_result = array(0, dim = c(ndraw, NS, NS))
for(i in 1:NS){
  for(j in 1:NS){
    p_result[ ,i ,j] = result$p_str[(0:(ndraw - 1)) * NS * NS + (i - 1) * NS + j]
  }
}
state_result = array(0, dim = c(ndraw, N, NT))
for(i in 1:NT){
  for(j in 1:N){
    state_result[ ,j ,i] = result$s_str[(0:(ndraw - 1)) * NT * N + (i - 1) * N + j]
  }
}
pi1_result = array(0, dim = c(ndraw, NS, ng))
for(i in 1:NS){
  for(j in 1:ng){
    pi1_result[ ,i ,j] = result$pi1_str[(0:(ndraw - 1)) * NS * ng + (i - 1) * ng + j]
  }
}

## posterior mean
sigma.estimation = apply(sigma_result, 1, mean)
alpha.estimation = apply(alpha_result, c(2:3), mean)
intercept.estimation = alpha.estimation[ , ncol(alpha.estimation)]
pi1.estimation = apply(pi1_result, c(2:3), mean)
sigma.sd = apply(sigma_result, 1, sd)
alpha.sd = apply(alpha_result, c(2:3), sd)
intercept.sd = alpha.sd[ , ncol(alpha.estimation)]

# sigma.estimation2 = apply(sigma_result, 1, median)
# alpha.estimation2 = apply(alpha_result, c(2:3), median)
# intercept.estimation2 = alpha.estimation2[ , ncol(alpha.estimation)]

sigma.estimation
intercept.estimation
alpha.estimation
sigma.sd
alpha.sd

state.res = array(0, dim=c(N, NT))
temp.mat = matrix(0, nrow = NT, ncol = NS)
for(i in 1:N){
  for(j in 1:NT){
    temp = rep(0,NS)
    for(k in 1:NS){
      temp[k] = sum(state_result[ ,i,j] == (k - 1))
    }
    temp.mat[j, ] = temp
    # if(j == 73){
    #   set.seed(7)
    # }
    state.res[i, j] = sample(c(1, 2), size = 1, prob = temp)
  }
}

state.estimation = c(state.res)
state.estimation
sig_x_NS = matrix(0, nrow = NS, ncol = Nbvar * eg + egd)
mu_x_NS = matrix(0, nrow = NS, ncol = Nbvar * eg + egd)
beta.estimation = alpha.estimation[ ,-ncol(alpha.estimation)]
beta.real.estimation = matrix(0, nrow = NS, ncol = Nbvar + 1)
beta.real.estimation2 = matrix(0, nrow = NS, ncol = Nbvar + 1)
alpha.no.int = alpha.estimation[ ,1:(Nbvar * eg + egd)]
for (s in 1:NS) {
  # stdres = standardize_me2(QX[state == s, ])
  sig_x_NS[s, ] = standardize_me2(QX[state.estimation == s, ])$sig
  mu_x_NS[s, ] = standardize_me2(QX[state.estimation == s, ])$mu
  intercept.estimation[s] = intercept.estimation[s] + mean(Y)
  for (nb in 1:Nbvar) {
    beta.real.estimation[s, nb] = sum((alpha.estimation[s, ((nb - 1)*eg + 1):(nb*eg)] / sig_x_NS[s, ((nb - 1)*eg + 1):(nb*eg)]) %*% Q)
    beta.real.estimation2[s, nb] = sum((alpha.estimation[s, ((nb - 1)*eg + 1):(nb*eg)] / sig_x[((nb - 1)*eg + 1):(nb*eg)]) %*% Q)
  }
  beta.real.estimation[s, Nbvar + 1] = sum((alpha.estimation[s, (Nbvar*eg + 1):(Nbvar*eg + egd)] / sig_x_NS[s, (Nbvar*eg + 1):(Nbvar*eg + egd)]) %*% Qd)
  beta.real.estimation2[s, Nbvar + 1] = sum((alpha.estimation[s, (Nbvar*eg + 1):(Nbvar*eg + egd)] / sig_x[(Nbvar*eg + 1):(Nbvar*eg + egd)]) %*% Qd)
}
intercept.estimation
beta.real.estimation
beta.real.estimation2
which(abs(beta.real.estimation[1,])>0.1)
which(abs(beta.real.estimation[2,])>0.1)

apply(p_result, c(2, 3), mean)

beta_result = array(0, dim = c(ndraw, NS, Nbvar + 1))
for (i in 1:ndraw) {
  for (s in 1:NS) {
    for (nb in 1:Nbvar) {
      beta_result[i, s, nb] = sum((alpha_result[i, s, ((nb - 1)*eg + 1):(nb*eg)] / sig_x_NS[s, ((nb - 1)*eg + 1):(nb*eg)]) %*% Q)
    }
    beta_result[i, s, Nbvar + 1] = sum((alpha_result[i, s, (Nbvar*eg + 1):(Nbvar*eg + egd)] / sig_x_NS[s, (Nbvar*eg + 1):(Nbvar*eg + egd)]) %*% Qd)
  }
}
beta.sd = apply(beta_result, c(2:3), sd)


################################## Plot ########################################
NBER = read.csv('./data/USREC.csv', header = TRUE)
NBER_noindex <- NBER['USREC']
dates_nber_all <- seq(as.Date("1854-12-01"), by = "month", length.out = dim(NBER)[1])
NBER_xts <- xts(NBER_noindex, order.by = dates_nber_all)
dates_nber.quarter <- seq(as.Date("1960-01-01"), by = "quarter", length.out = NT_train)
dates_nber.month <- seq(as.Date("1960-01-01"), by = "month", length.out = NT_train * m_mon)
NBER_used <- NBER_xts['1960-01-01/2010-12-01']

df.nber <- data.frame(
  Date = dates_nber.month, 
  Recession = c(as.numeric(NBER_used))
)

recession_periods <- subset(df.nber, NBER == 1)
df.nber$Date <- as.Date(df.nber$Date)

recession_periods <- df.nber %>%
  mutate(Recession_Change = Recession != lag(Recession, default = 0)) %>%
  mutate(Period_ID = cumsum(Recession_Change)) %>%
  group_by(Period_ID) %>%
  filter(Recession == 1) %>%
  summarize(Start = min(Date), End = max(Date))

df.gdp.growth.color <- data.frame(
  Date = dates_nber.quarter,
  GDP.growth = c(Y),
  state = as.character(state.estimation)
)

plot.gdp.growth.color.hmsp <- ggplot() +
  geom_rect(data = recession_periods, aes(xmin = Start, xmax = End, ymin = -2, ymax = 6),
            fill = "grey", alpha = 0.5) + 
  geom_point(data = df.gdp.growth.color, aes(x = Date, y = GDP.growth, color = state)) +
  scale_color_manual(values = c("1" = "red", "2" = "blue")) +
  labs(y = TeX("$y_t"), x = "Date") + 
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 12),     
    plot.title = element_text(size = 18, face = "bold"),  
    legend.text = element_text(size = 14),   
    legend.title = element_text(size = 16)   
  ) +
  scale_x_date(date_breaks = "4 year", date_labels = "%Y") +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  theme(axis.text.y = element_text())  
show(plot.gdp.growth.color.hmsp)
ggsave("./code/realAnalysis/Plots/HMSP-MIDAS.png", plot = plot.gdp.growth.color.hmsp, width = 8, height = 6, dpi = 300)

################################## Forecast ########################################

######### Prepare Test Data ############
gdp <- GDP_xts['1959-03-01/2023-12-01']
gdp.test <- GDP_xts['2011-03-01/2023-12-01']
len.test = as.numeric(length(gdp.test))
gdp_rate_test = rep(0, len.test - 1)
for (i in 1:(len.test - 1)) {
  gdp_rate_test[i] = 100 * log(as.numeric(gdp.test[i + 1]) / as.numeric(gdp.test[i]))
}
growth_dates_test = dates_test <- seq(as.Date("2011-03-01"), by = "quarter", length.out = len.test - 1)
gdp_growth_test = xts(gdp_rate_test, order.by = growth_dates_test)
gdp_growth_test_num = as.numeric(gdp_growth_test) 
gdp_growth.test <- gdp_growth_test["2011-01-01/2023-12-01"]
Y.test <- gdp_growth.test

fredmd.test <- fredmd[(NT_month.train - 12 + 1):(NT_month.train + NT_month.test), ]
# process NA
sum(is.na(fredmd.test))  # 387 NA need to be filled.
has_na <- colSums(is.na(fredmd.test)) > 0
na_columns <- names(fredmd.test)[has_na]
na_columns
na_counts <- colSums(is.na(fredmd.test[has_na]))
na_counts

nac1 <- as.numeric(fredmd.test$`S&P div yield`)
nac2 <- as.numeric(fredmd.test$`S&P PE ratio`)
nac3 <- as.numeric(fredmd.test$CP3Mx)
nac1_omit = na.omit(nac1)
nac2_omit = na.omit(nac2)
nac3_omit = na.omit(nac3)
nac1_mean = mean(nac1_omit)
nac1_sd = sd(nac1_omit)
nac2_mean = mean(nac2_omit)
nac2_sd = sd(nac2_omit)
nac3_mean = mean(nac3_omit)
nac3_sd = sd(nac3_omit)
nac1_count <- sum(is.na(fredmd.test$`S&P div yield`))  
fill_values1 <- rnorm(nac1_count, mean = nac1_mean, sd = nac1_sd)
fredmd.test$`S&P div yield`[is.na(fredmd.test$`S&P div yield`)] <- fill_values1
nac2_count <- sum(is.na(fredmd.test$`S&P PE ratio`))  
fill_values2 <- rnorm(nac2_count, mean = nac2_mean, sd = nac2_sd)
fredmd.test$`S&P PE ratio`[is.na(fredmd.test$`S&P PE ratio`)] <- fill_values2
nac3_count <- sum(is.na(fredmd.test$CP3Mx))  
fill_values3 <- rnorm(nac3_count, mean = nac3_mean, sd = nac3_sd)
fredmd.test$CP3Mx[is.na(fredmd.test$CP3Mx)] <- fill_values3

fredmd.test$`S&P div yield` <- as.numeric(fredmd.test$`S&P div yield`)
fredmd.test$CP3Mx <- as.numeric(fredmd.test$CP3Mx)

for(i in 1:ncol(fredmd.test)) {
  fredmd.test[ , i][is.na(fredmd.test[ , i])] <- median(fredmd.test[ , i], na.rm=TRUE)
}
fredmd.test <- as.data.frame(fredmd.test)

has_na <- colSums(is.na(fredmd.test)) > 0
na_columns <- names(fredmd.test)[has_na]
na_columns
na_counts <- colSums(is.na(fredmd.test[has_na]))
na_counts

# summary descriptive statistics
fredmd.test.summ <- sapply(fredmd.test, function(x) c(Mean = mean(x, na.rm = TRUE),
                                                      SD = sd(x, na.rm = TRUE),
                                                      Min = min(x, na.rm = TRUE),
                                                      Max = max(x, na.rm = TRUE)))

fredmd.test <- subset(fredmd.test, select = -c(VIXCLSx))
all.names <- c(names(fredmd.test))
fredmd.test <- subset(fredmd.test, select = -c(FEDFUNDS, CP3Mx, TB3MS, TB6MS, GS1, GS5, GS10, AAA, BAA))
fredmd.test[, "HWI"] <- fredmd.test[, "HWI"] / 100
fredmd.test[, c("CES0600000007", "AWHMAN")] <- fredmd.test[, c("CES0600000007", "AWHMAN")] / 10

ff3.test <- ff3[(NT_month.train - 12 + 1):(NT_month.train + NT_month.test), ]
ads.test.org <- ads_xts['2010-12-01/2023-12-01']
Xm.test.org <- as.data.frame(cbind(fredmd.test, ff3.test))
Xm.test <- as.data.frame(matrix(0, nrow = nrow(Xm.test.org), ncol = ncol(Xm.test.org)))
for (m in 1:ncol(Xm.test)) {
  Xm.test[ ,m] <- (Xm.test.org[ ,m] - Xm.mean[m]) / Xm.std[m]
}
ads.test <- (ads.test.org - ads.mean) / ads.std

NT_gdp.test = as.numeric(nrow(gdp_growth.test))
NT_test <- NT_gdp.test
K_mon = Spc$Km 
m_mon = as.numeric(nrow(fredmd.test) %/% nrow(gdp_growth.test))
# Monthly matrix
xmls.test <- matrix(0, nrow = NT_gdp.test, ncol = 0)
xx.test <- vector("list", length = Nbvar) # monthly is Nbvar
for (ii in 1:Nbvar) {
  xx.test[[ii]] <- Construct_DataMIDAS(g = gdp_growth.test, d = Xm.test[, ii], K = K_mon, m = m_mon)
  xmls.test <- cbind(xmls.test, xx.test[[ii]])
}
xd2.test <- vector("list", length = 1)
m_day2.test = as.numeric(length(ads.test) %/% nrow(gdp_growth.test)) 
K_day2.test = m_day2.test
xd2.test[[1]] <- Construct_DataMIDAS(g = gdp_growth.test, d = ads.test, K = K_day2.test, m = m_day2.test)
xmls.test <- cbind(xmls.test, xd2.test[[1]][1:NT_gdp.test, ]) 
X.test <- xmls.test[1:NT_gdp.test, ]

polydegree = 2
Ra = vector('list',2)
names(Ra) = c('fC','dfC')
Ra$fC = 1 # no tail restriction; if has set 1
Ra$dfC = 0 # no deravative restriction; if has set 1
Q = Almon_lag(polydegree = polydegree, C = Spc$Km, R = Ra)
Qrow = as.numeric(dim(Q)[1])
eg <- Qrow
NX_alm = Qrow * Nbvar
Qd.test = Almon_lag(polydegree = polydegree, C = m_day2.test, R = Ra)
Qdrow.test = as.numeric(dim(Qd.test)[1])
egd.test <- Qdrow.test
kdm.test <- as.numeric(dim(xd2.test[[1]])[2])
# QX.test = matrix(0, nrow = nrow(X.test), ncol = Nbvar * eg)
QX.test = matrix(0, nrow = nrow(X.test), ncol = Nbvar * eg + egd.test)
for (nb in 1:Nbvar) {
  QX.test[ ,((nb - 1)*eg + 1):(nb*eg)] = X.test[ ,((nb - 1)*K_mon + 1):(nb*K_mon)] %*% t(Q)
}
QX.test[ ,(Nbvar*eg + 1):(Nbvar*eg + egd.test)] =  X.test[ ,(Nbvar*K_mon+1):(Nbvar*K_mon + kdm.test)] %*% t(Qd.test)

p = eg * Nbvar + egd.test
NF <- p
ng <- Nbvar + 1  # number of groups. nrow(Q) is the polydegree + 1(contain 0).
nj <- c(rep(nrow(Q), Nbvar), egd.test)  # number of predictors in each group
groups <- matrix(0, nrow = ng, ncol = nrow(Q)) 
groups[1, ] <- seq_len(nj[1])
if (ng > 1){
  for (i in 2:ng) {
    groups[i, ] <- seq(sum(nj[1:(i - 1)]) + 1, sum(nj[1:(i - 1)]) + nj[i])
  }
}

# Standalization for data. 
x.test = array(0, dim = c(N, NT_test, NF))
y.test = array(0, dim = c(NT_test, N))
for (n in 1:N) {
  # result_x <- standardize_me2(QX.test)
  result_x <- standardize_me(QX.test)
  # result_y <- center_me(as.matrix(Y.test))
  result_y <- as.matrix(Y.test)
  x_scale = result_x$y
  mu_x = result_x$mu
  sig_x = result_x$sig
  y_scale = result_y
  x.test[n, , ] <- x_scale
  y.test[ ,n] <- y_scale
}


# For Forecast with expanding window
forecast_prob <- matrix(0, nrow = NT_test, ncol = NS)
state_test_est <- array(0, dim = c(NT_test))
mean_test_fore <- array(dim = c(N, NT_test))
Y_test_fore <- array(dim = c(N, NT_test))
plast_prob <- temp.mat[nrow(temp.mat), ] / sum(temp.mat[nrow(temp.mat), ])
p_test_est <- apply(p_result, c(2, 3), mean)
tran_p_1 <- p_test_est
new_fore_prob <- c(rep(0, NS))
for (s1 in 1:NS) {
  for (s2 in 1:NS) {
    new_fore_prob[s1] = new_fore_prob[s1] + plast_prob[s2] * tran_p_1[s2, s1]
  }
}
forecast_prob[1, ] <- new_fore_prob
state_test_est[1] <- sample(c(1:NS), size = 1, replace = FALSE, prob = forecast_prob[1, ])
mean_test_fore[1, 1] <- t(x.test[1, 1, ]) %*% beta.estimation[state_test_est[1], ] + intercept.estimation[state_test_est[1]]
Y_test_fore[1, 1] <- rnorm(1, mean = mean_test_fore[1, 1], sd = sqrt(sigma.estimation[state_test_est[1]]))

