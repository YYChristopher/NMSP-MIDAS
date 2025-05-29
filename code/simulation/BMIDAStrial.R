library("tseries")
library("midasr")
library("ggplot2")
library(forecast)
library(verification)
library(data.table)
library(gtools)


#Work Space of the data
work_dir = 'F:/24561/Documents/ResearchWrite/NHMM-MIDAS/Submission Requirements/Code in Github'
#Set the Work Space of the data
setwd(work_dir)

source('./code/simulation/BMIDAS.R')

N = 1
NT = 500 
NT_oos = 50 
Rep = 1

# Total number of variables
Nbvar <- 30
# Number of relevant variables
Nbrelevant <- 5
# Number of irrelevant variables
Nbirrelevant <- Nbvar - Nbrelevant

# Basic specifications for explanatory variables (monthly)
Spc = vector('list', 6)
names(Spc) = c('mseries', 'monthly', 'Km', 'kappam', 'sK', 'nbvar')
Spc$mseries = NULL
Spc$monthly = Nbvar
Spc$Km = 24 # lag orders
Spc$kappam = 3
Spc$sK = rep(Spc$Km, Spc$monthly)
NX = Spc$Km * Nbvar + 1 

# construction of transform matrix Q
# Constrained Almon lag polynomial (monthly). 
polydegree = 3 # Almon lag degree
Ra = vector('list',2)
names(Ra) = c('fC','dfC')
Ra$fC = 1 # no tail restriction; if has set 1
Ra$dfC = 1 # no deravative restriction; if has set 1
Q = Almon_lag(polydegree = polydegree, C = Spc$Km, R = Ra)
Qrow = as.numeric(dim(Q)[1])
NX_alm = Qrow * Nbvar

It = vector('list',3)
names(It) = c('nsave','nburn','nthin')
It$nsave = 10000 # total sample
It$nburn = 10000 # burn sample
It$nthin = 1

p = NX_alm 
ng <- p / nrow(Q)  # number of groups. nrow(Q) is the polydegree + 1(contain 0).
nj <- rep(nrow(Q), ng)  # number of predictors in each group
groups <- matrix(0, nrow = ng, ncol = nrow(Q)) 
groups[1, ] <- seq_len(nj[1])
for (i in 2:ng) {
  groups[i, ] <- seq(sum(nj[1:(i-1)]) + 1, sum(nj[1:(i-1)]) + nj[i])
}

all_y_forecast <- array(0, dim = c(Rep, N, NT_oos))
all.beta.estimate <- array(0, dim = c(Rep, N, Nbvar))
all_y_rmse <- array(0, dim = Rep)
all_y_rsquare <- array(0, dim = Rep)
all_y_mae <- array(0, dim = Rep)
all_y_crps <- array(0, dim = Rep)

for (r in 1:Rep) {
  
  set.seed(29)
  directory_path <- "./code/simulation/Gen_Data_Train_NS=2"  # use offset data to illustrate that BMIDAS can't handle offset effect
  specific_number <- r
  # file_names <- list.files(path = directory_path, pattern = paste0("^.*", specific_number, "[^0-9].*\\.csv$"), full.names = TRUE)
  pattern <- sprintf("Train_rep%d_chain[0-9]+\\.csv", r)
  files <- list.files(directory_path, pattern = pattern, full.names = TRUE)
  file_names <- sort(files)
  file_names <- mixedsort(file_names, decreasing = FALSE) 
  # read files
  data_list <- lapply(file_names, read.csv)
  total_dim = as.numeric(dim(data_list[[1]])[2])
  D_train = array(0, dim = c(N, NT, total_dim))
  Y_train = array(0, dim = c(N, NT))
  X_train = array(0, dim = c(N, NT, NX))
  h_train = array(0, dim = c(N, NT, NH))
  X_train_alm = array(0, dim = c(N, NT, NX_alm))
  X_test_alm = array(0, dim = c(N, NT_oos, NX_alm))
  state_train = array(0, dim = c(N, NT))
  for (n in 1:N) {
    D_train[n, , ] <- matrix(unlist(data_list[[n]]), ncol = total_dim, byrow = FALSE)
    Y_train[n, ] <- D_train[n, , 1]
    X_train[n, , ] <- D_train[n, , 2:(NX + 1)]
    X_train_alm[n, , ] <- D_train[n, , (NX + 2):(NX + NX_alm + 1)]
    h_train[n, , ] <- D_train[n, , (NX + NX_alm + 2):(1 + NX + NX_alm + NH)]
    state_train[n, ] <- D_train[n, , total_dim]
  }
  
  D_train_mat = array(0, dim = c(N*NT, total_dim))
  # Y_train_mat = array(0, dim = c(N*NT))
  # X_train_mat = array(0, dim = c(N*NT, NX))
  # h_train_mat = array(0, dim = c(N*NT, NH))
  # X_train_alm_mat = array(0, dim = c(N*NT, NX_alm))
  # X_test_alm_mat = array(0, dim = c(N*NT_oos, NX_alm))
  # state_train_mat = array(0, dim = c(N*NT))
  for (n in 1:N) {
    D_train_mat[((n-1)*NT+1):(n*NT), ] <- matrix(unlist(data_list[[n]]), ncol = total_dim, byrow = FALSE) # 注意要先将list数据Unlist再填入矩???
  }
  Y_train_mat <- as.matrix(D_train_mat[ ,1])
  X_train_mat <- D_train_mat[ ,2:(NX + 1)]
  X_train_alm_mat <- D_train_mat[ ,(NX + 2):(NX + NX_alm + 1)]
  
  # # Standalization for data. 
  # x = array(0, dim = c(N, NT, NX_alm))
  # y = array(0, dim = c(N, NT))
  # for (n in 1:N) {
  #   result_x <- standardize_me(X_train_alm[n, , ])
  #   result_y <- center_me(as.matrix(Y_train[n, ]))
  #   # result_y <- standardize_me(as.matrix(Y_train[n, ]))
  #   x_scale = result_x$y
  #   mu_x = result_x$mu
  #   sig_x = result_x$sig
  #   y_scale = result_y$y
  #   mu_y = result_y$mu 
  #   x[n, , ] <- x_scale
  #   y[n, ] <- y_scale
  # }
  
  x = X_train_alm
  y = Y_train
  
  test_path <- "./code/simulation/Gen_Data_Test_NS=2"
  specific_number <- r
  # file_names <- list.files(path = directory_path, pattern = paste0("^.*", specific_number, "[^0-9].*\\.csv$"), full.names = TRUE)
  pattern <- sprintf("Test_rep%d_chain[0-9]+\\.csv", r)
  files <- list.files(test_path, pattern = pattern, full.names = TRUE)
  file_names <- sort(files)
  file_names <- mixedsort(file_names, decreasing = FALSE) 
  data_list <- lapply(file_names, read.csv)
  total_dim = as.numeric(dim(data_list[[1]])[2])
  
  D_test = array(0, dim = c(N, NT_oos, total_dim))
  Y_test = array(0, dim = c(N, NT_oos))
  X_test = array(0, dim = c(N, NT_oos, NX))
  X_test_alm = array(0, dim = c(N, NT_oos, NX_alm))
  for (n in 1:N) {
    D_test[n, , ] <- matrix(unlist(data_list[[n]]), ncol = total_dim, byrow = FALSE)
    Y_test[n, ] <- D_test[n, , 1]
    X_test[n, , ] <- D_test[n, , 2:(NX + 1)]
    X_test_alm[n, , ] <- D_test[n, , (NX + 2):(NX + NX_alm + 1)]
  }
  
  
  Y_test_fore <- array(dim = c(N, NT_oos))
  for (n in 1:N) {
    print(paste('sample chains', n, "\r"))
    X_reg = as.matrix(x[n, , ])
    Y_reg = as.matrix(y[n, ])
    X_for = as.matrix(X_test_alm[n, , ])
    result = BMIDAS_AGLasso_SS(X_reg, Y_reg, X_for, It, Q) # use BMIDAS-AGL-SS
    # result = BMIDAS_AGLasso(X_reg, Y_reg, X_for, It, Q) # use BMIDAS-AGL
    Y_test_fore[n, ] = result$y_fore[It$nsave, ]
    theta.est = apply(result$beta, 2, mean)
    beta.est = c(rep(0, Nbvar))
    for (nb in 1:Nbvar) {
      beta.est[nb] <- sum(theta.est[groups[nb, ]] %*% Q)
    }
    all.beta.estimate[r, n, ] = beta.est
  }
  
  Y_test_fore_vec = sample_merge_Ystate(Y_test_fore)
  Y_test_true_vec = sample_merge_Ystate(Y_test)
  y_rmse = rmse(Y_test_true_vec, Y_test_fore_vec)
  y_mae = mae(Y_test_true_vec, Y_test_fore_vec)
  mean_fore = mean(Y_test_fore_vec)
  std_fore = sd(Y_test_fore_vec)
  y_crps = crps2(Y_test_true_vec, Y_test_fore_vec)
  all_y_forecast[r, , ] = Y_test_fore
  all_y_rmse[r] = y_rmse
  all_y_mae[r] = y_mae
  all_y_crps[r] = y_crps
  print(paste("Trial times", r, "have completed...", "\r"))
}

# print(all_y_rmse)
# print(all_y_mae)
# print(all_y_crps)

y_rmse_res = mean(all_y_rmse)
y_mae_res = mean(all_y_mae)
y_crps_res = mean(all_y_crps)
print(y_rmse_res)
print(y_mae_res)
print(y_crps_res)
eva_list <- list(all_y_rmse = all_y_rmse, all_y_mae = all_y_mae, 
                 all_y_crps = all_y_crps)
char_eva_list <- lapply(eva_list, function(x) paste(x, collapse = ", "))
# If use BMIDAS-AGL-SS, write into csv
writeLines(unlist(char_eva_list), paste0("./code/simulation/Result_BMIDAS", "/", 
                                         "BMIDAS-AGL-SS", ".txt")) # If use BMIDAS-AGL-SS, write into csv
# writeLines(unlist(char_eva_list), paste0("./code/simulation/Result_BMIDAS", "/", 
#                                        "BMIDAS-AGL", ".txt")) # If use BMIDAS-AGL, write into csv
