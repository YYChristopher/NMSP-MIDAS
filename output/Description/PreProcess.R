rm(list=ls())
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
sourceCpp("./code/realAnalysis/update.cpp")

#################################### Data Process ########################################
# GDP Data to transform into GDP growth rate
# truncate data that we need
GDP = read.csv('./data/GDP2.csv', header = TRUE)
GDP_noindex <- GDP['GDP']
dates_gdp <- seq(as.Date("1947-03-01"), by = "quarter", length.out = 308)
GDP_xts <- xts(GDP_noindex, order.by = dates_gdp)
gdp <- GDP_xts['1959-03-01/2023-12-01']
gdp.train <- GDP_xts['1959-03-01/2023-12-01']
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
gdp_growth.train <- gdp_growth["1958-12-01/2023-09-01"]
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
# Observe the data distribution of the three columns with more missing values
# (ACOGNO-255 TWEXAFEGSMTHx-25 UMCSENTx-86)
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

nac1_count <- sum(is.na(fredmd.train$ACOGNO))  # Count the number of NAs in column A
fill_values1 <- rnorm(nac1_count, mean = nac1_mean, sd = nac1_sd)
fredmd.train$ACOGNO[is.na(fredmd.train$ACOGNO)] <- fill_values1
nac2_count <- sum(is.na(fredmd.train$TWEXAFEGSMTHx)) 
fill_values2 <- rnorm(nac2_count, mean = nac2_mean, sd = nac2_sd)
fredmd.train$TWEXAFEGSMTHx[is.na(fredmd.train$TWEXAFEGSMTHx)] <- fill_values2
nac3_count <- sum(is.na(fredmd.train$UMCSENTx))  
fill_values3 <- rnorm(nac3_count, mean = nac3_mean, sd = nac3_sd)
fredmd.train$UMCSENTx[is.na(fredmd.train$UMCSENTx)] <- fill_values3

# Note that the following two columns of data are not in numeric format but in character format. 
# They need to be converted to numeric format.
fredmd.train$`S&P div yield` <- as.numeric(fredmd.train$`S&P div yield`)
fredmd.train$CP3Mx <- as.numeric(fredmd.train$CP3Mx)
# Median imputation method (for situations where the number of missing values is not large)
for(i in 1:ncol(fredmd.train)) {
  fredmd.train[ , i][is.na(fredmd.train[ , i])] <- median(fredmd.train[ , i], na.rm=TRUE)
}
fredmd.train <- as.data.frame(fredmd.train)

# fredmd_filled <- fredmd %>%
# mutate(across(where(is.numeric), ~if_else(is.na(.), median(., na.rm = TRUE), .)))
# Check for missing values
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
# Constructing X matrices in the MIDAS form
# extract data that we need (slope of the yield)
fredmd.train <- subset(fredmd.train, select = -c(VIXCLSx))
all.names <- c(names(fredmd.train))
fredmd.train <- subset(fredmd.train, select = -c(FEDFUNDS, CP3Mx, TB3MS, TB6MS, GS1, GS5, GS10, AAA, BAA))
fredmd.train[, "HWI"] <- fredmd.train[, "HWI"] / 100
fredmd.train[, c("CES0600000007", "AWHMAN")] <- fredmd.train[, c("CES0600000007", "AWHMAN")] / 10
ads.train.org <- ads_xts['1960-03-01/2023-12-01']
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
Nbvar_month = as.numeric(ncol(Xm.train))
Nbvar <- Nbvar_month
Spc = vector('list', 6)
names(Spc) = c('mseries', 'monthly', 'Km', 'kappam', 'sK', 'nbvar')
Spc$mseries = NULL
Spc$monthly = Nbvar
Spc$Km = 12 # lag orders. 
Spc$kappam = 3 # The ratio of high frequency to low frequency
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
# Response y
y.train = gdp_growth.train[-(1:notuseindex)]
x_train = Xused
y_train = y.train
real_date <- growth_dates[-(1:notuseindex)]
real_date_train <- real_date[1:NT_train]

# Data in transition model
filepath <- './data/fredqdADS.csv'
fred_qd <- fredqd(filepath, date_start = NULL, date_end = NULL, transform = TRUE)
pq <- ncol(fred_qd)
dates_fredqd <- seq(as.Date("1959-03-01"), by = "quarter", length.out = 260)
fredqd_xts <- xts(fred_qd, order.by = dates_fredqd)[ ,-1]
fredqd <- as.data.frame(fredqd_xts['1959-01-01/2023-12-01'])
fredqd[] <- lapply(fredqd, function(x) as.numeric(as.character(x)))
D <- as.matrix(cbind(fredqd$`S&P 500` * 100, fredqd$WPU0561 * 10, fredqd$PPIACO * 50, fredqd$ADS * 100))
D <- D[-(1:notuseindex), ]
NH <- ncol(D)
D_train <- D[1:NT_train, ]

# Real State
RealState <- read.csv('./output/Description/RealStateAll.csv', header = TRUE)
RealStatevec <- RealState['State']
dates_realState <- seq(as.Date("1959-03-01"), by = "quarter", length.out = 259)
realState_xts <- xts(RealStatevec, order.by = dates_realState)
realState.train <- realState_xts['1960-03-01/2023-12-01']
real_state_train <- as.vector(realState.train)

# set.seed(43)
Rep = 1
Y <- as.vector(y_train)
Y.pre <- as.vector(c(0, Y[-length(Y)]))
# Y.front <- as.vector(Y[-length(Y)])
# Y.back <- as.vector(Y[-1])
X <- x_train
N = 1
NS = 2

X.ads <- X[ ,1453:1541]
X.day <- matrix(0, nrow = NT)
for (d in 1:ncol(X.ads)) {
  if (d %% 7 == 0){
    X.day <- cbind(X.day, X.ads[ ,d])
  }
}
X.day <- X.day[ ,-1]
X.use <- cbind(X[ ,1:1452], X.day)

set.seed(as.numeric(Sys.time()) + sample(1:10000, 1, replace = FALSE))
# some basic information for simulated data
# p = eg * Nbvar
eg <- egd <- 2
p = eg * Nbvar + egd
NF <- p

# Standalization for data. 
x = array(0, dim = c(N, NT, NF * 12))
y = array(0, dim = c(NT, N))
d = array(0, dim = c(N, NT, NH))
for (n in 1:N) {
  # result_x <- standardize_me2(QX)
  result_x <- standardize_me(X.use)
  # result_y <- center_me(as.matrix(Y))
  result_y <- Y
  x_scale = result_x$y
  mu_x = result_x$mu
  sig_x = result_x$sig
  y_scale = result_y
  # y_scale = result_y$y
  # mu_y = result_y$mu
  x[n, , ] <- x_scale
  y[ ,n] <- y_scale
  # x[n, , ] <- QX
  # y[ ,n] <- Y
  d[n, , ] <- D_train
}
x <- x[1, , ]
