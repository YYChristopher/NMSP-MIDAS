# import packages
library(LaplacesDemon)
library(gtools)
library(data.table)
library(ggplot2)
library(verification)

#Work Space of the data
work_dir = 'F:/24561/Documents/ResearchWrite/NHMM-MIDAS/Submission Requirements/Code in Github'
#Set the Work Space of the data
setwd(work_dir)
source('./code/simulation/tool.R')

Rep = 1
N = 1 # Number of subjects
NT = 500 # Number of in-sample
NT_oos = 50 # Number of out-of-sample
NS = 2
NH = 3
y <- rep(1, NT)
yoos <- rep(1, NT_oos)

# Total number of variables
Nbvar <- 30
# Nbvar <- 5
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

polydegree = 3 # Almon lag degree
Ra = vector('list',2)
names(Ra) = c('fC','dfC')
Ra$fC = 1 # no tail restriction; if has set 1
Ra$dfC = 1 # no deravative restriction; if has set 1
Q = Almon_lag(polydegree = polydegree, C = Spc$Km, R = Ra)
Qrow = as.numeric(dim(Q)[1])
NX_alm = Qrow * Nbvar

all_state_test_est <- array(0, dim = c(Rep, N, NT_oos))
all_y_forecast <- array(0, dim = c(Rep, N, NT_oos))
all_y_rmse <- array(0, dim = Rep)
all_y_rsquare <- array(0, dim = Rep)
all_y_mae <- array(0, dim = Rep)
all_y_crps <- array(0, dim = Rep)
all_conf_mat_test <- array(0, dim = c(Rep, NS, NS))

####### Forecast(Out-of-sample)
# read test data
set.seed(as.numeric(Sys.time()) + sample(1:10000, 1, replace = FALSE))
start_time = Sys.time()
for (r in 1:Rep) {
  
  directory_path <- "./code/simulation/Gen_Data_Train_NS=2"
  specific_number <- r
  # Regular expression matching files
  # file_names <- list.files(path = directory_path, pattern = paste0("^.*", specific_number, "[^0-9].*\\.csv$"), full.names = TRUE)
  pattern <- sprintf("Train_rep%d_chain[0-9]+\\.csv", r)
  # list all .csv file
  files <- list.files(directory_path, pattern = pattern, full.names = TRUE)
  # order the .csv file
  file_names <- sort(files)
  file_names <- mixedsort(file_names, decreasing = FALSE) 
  # read the file
  data_list <- lapply(file_names, read.csv)
  total_dim = as.numeric(dim(data_list[[1]])[2])
  # put list into array
  D_train = array(0, dim = c(N, NT, total_dim))
  Y_train = array(0, dim = c(N, NT))
  X_train = array(0, dim = c(N, NT, NX))
  h_train = array(0, dim = c(N, NT, NH))
  X_train_alm = array(0, dim = c(N, NT, NX_alm))
  X_test_alm = array(0, dim = c(N, NT_oos, NX_alm))
  state_train = array(0, dim = c(N, NT))
  for (n in 1:N) {
    D_train[n, , ] <- as.matrix(data_list[[n]][1:NT, ])
    Y_train[n, ] <- D_train[n, , 1]
    X_train[n, , ] <- D_train[n, , 2:(NX + 1)]
    X_train_alm[n, , ] <- D_train[n, , (NX + 2):(NX + NX_alm + 1)]
    h_train[n, , ] <- D_train[n, , (NX + NX_alm + 2):(1 + NX + NX_alm + NH)]
    state_train[n, ] <- D_train[n, , total_dim]
  }
  
  test_path <- "./code/simulation/Gen_Data_Test_NS=2"
  specific_number <- r
  # Regular expression matching files
  # file_names <- list.files(path = directory_path, pattern = paste0("^.*", specific_number, "[^0-9].*\\.csv$"), full.names = TRUE)
  pattern <- sprintf("Test_rep%d_chain[0-9]+\\.csv", r)
  # list all .csv file
  files <- list.files(test_path, pattern = pattern, full.names = TRUE)
  # order the .csv file
  file_names <- sort(files)
  file_names <- mixedsort(file_names, decreasing = FALSE) 
  # read files
  data_list <- lapply(file_names, read.csv)
  total_dim = as.numeric(dim(data_list[[1]])[2])
  
  D_test = array(0, dim = c(N, NT_oos, total_dim))
  Y_test = array(0, dim = c(N, NT_oos))
  X_test = array(0, dim = c(N, NT_oos, NX))
  h_test = array(0, dim = c(N, NT_oos, NH))
  X_test_alm = array(0, dim = c(N, NT_oos, NX_alm))
  state_test = array(0, dim = c(N, NT_oos))
  for (n in 1:N) {
    D_test[n, , ] <- as.matrix(data_list[[n]][1:NT_oos, ])
    Y_test[n, ] <- D_test[n, , 1]
    X_test[n, , ] <- D_test[n, , 2:(NX + 1)]
    X_test_alm[n, , ] <- D_test[n, , (NX + 2):(NX + NX_alm + 1)]
    h_test[n, , ] <- D_test[n, , (NX + NX_alm + 2):(1 + NX + NX_alm + NH)]
    state_test[n, ] <- D_test[n, , total_dim]
  }
  # state_begin = state_test[ ,1]
  
  D_test_mat = array(0, dim = c(N*NT, total_dim))
  for (n in 1:N) {
    D_test_mat[((n-1)*NT+1):(n*NT), ] <- matrix(unlist(data_list[[n]]), ncol = total_dim, byrow = FALSE) # 娉ㄦ剰瑕佸厛灏唋ist鏁版嵁Unlist鍐嶅～鍏ョ煩闃?
  }
  Y_test_mat <- as.matrix(D_test_mat[ ,1])
  X_test_mat <- D_test_mat[ ,2:(NX + 1)]
  X_test_alm_mat <- D_test_mat[ ,(NX + 2):(NX + NX_alm + 1)]
  h_test_mat <- D_test_mat[ ,(NX + NX_alm + 2):(1 + NX + NX_alm + NH)]
  state_test_mat <- as.matrix(D_test_mat[ ,total_dim])
  
  # Load Parameters Estimations
  parameters_path <- "./code/simulation/Paras_estimation_HMM"
  para_names <- sprintf("%s/Paras_est_rep%d.txt", parameters_path, r)
  para_names <- mixedsort(para_names, decreasing = FALSE) 
  para_lines <- readLines(para_names)
  # analyse each data
  parsed_data <- lapply(para_lines, function(line) unlist(strsplit(line, ", ")))
  beta_est <- matrix(as.numeric(parsed_data[[1]]), nrow = NS, ncol = Nbvar) 
  beta_org_est <- matrix(as.numeric(parsed_data[[2]]), nrow = NS, ncol = NX_alm)
  intercept_est <- as.numeric(parsed_data[[3]]) 
  sigma_est <- as.numeric(parsed_data[[4]])
  p_est <- matrix(as.numeric(parsed_data[[5]]), nrow = NS, ncol = NS)
  plast_prob <- as.numeric(parsed_data[[6]])
  
  # standardize test set of x
  x = array(0, dim = c(N, NT_oos, NX_alm))
  y = array(0, dim = c(N, NT_oos))
  for (n in 1:N) {
    result_x <- standardize_me(X_train_alm[n, , ])
    # result_y <- center_me(as.matrix(Y_train[n, ]))
    # result_y <- as.matrix(Y_train[n, ])
    x_scale = result_x$y
    mu_x = result_x$mu
    sig_x = result_x$sig
    x_scale <- (X_test_alm[n, , ] - matrix(rep(mu_x, each = nrow(X_test_alm[n, , ])), nrow = nrow(X_test_alm[n, , ]), byrow = FALSE)) /
      matrix(rep(sig_x, each = nrow(X_test_alm[n, , ])), nrow = nrow(X_test_alm[n, , ]), byrow = FALSE)
    y_scale <- Y_test[n, ]
    x[n, , ] <- x_scale
    y[n, ] <- y_scale
  }
  
  trial_count = 10000
  Y_test_fore_total = array(0, dim = c(trial_count, N, NT_oos))
  for (tr in 1:trial_count) {
    # Calculate the initial probability distribution
    p0_est <- plast_prob
    initial_state <- sample(c(1:NS), size = 1, prob = p0_est)
    # Calculate state sequence
    state_test_est <- tran_state_hmm(p_est, initial_state, NT_oos)
    # state_test_est <- tran_state_new(tran_p_test_est, state_begin)
    # Assume all_t_state has been initialized and has the appropriate dimensions
    all_state_test_est[r, , ] <- state_test_est
    state_test_est_vec = sample_merge_Ystate(state_test_est)
    state_test_vec = sample_merge_Ystate(state_test)
    # conf_mat_test = confuse_matrix(state_test_est_vec, state_test_vec, NS)
    # print(sum(state_test_est == state_test) / (N * NT_oos))
    
    location1_test <- which(state_test_est == 1, arr.ind = TRUE)
    location2_test <- which(state_test_est == 2, arr.ind = TRUE) 
    
    mean_test_fore <- array(dim = c(N, NT_oos))
    Y_test_fore <- array(dim = c(N, NT_oos))
    for (i in 1:nrow(location1_test)) {
      idx <- location1_test[i,]
      mean_test_fore[idx[1], idx[2]] <- t(x[idx[1], idx[2], ]) %*% beta_org_est[1, ] + intercept_est[1]
      Y_test_fore[idx[1], idx[2]] <- rnorm(1, mean = mean_test_fore[idx[1], idx[2]], sd = sigma_est[1])
    }
    
    for (i in 1:nrow(location2_test)) {
      idx <- location2_test[i,]
      mean_test_fore[idx[1], idx[2]] <- t(x[idx[1], idx[2], ]) %*% beta_org_est[2, ] + intercept_est[2]
      Y_test_fore[idx[1], idx[2]] <- rnorm(1, mean = mean_test_fore[idx[1], idx[2]], sd = sigma_est[2])
    }
    
    Y_test_fore_total[tr, , ] = Y_test_fore
  }
  Y_test_fore_res = apply(Y_test_fore_total, c(2, 3), mean)
  
  Y_test_fore_vec = sample_merge_Ystate(Y_test_fore_res)
  Y_test_true_vec = sample_merge_Ystate(Y_test)
  y_rmse = rmse(Y_test_true_vec, Y_test_fore_vec)
  y_mae = mae(Y_test_true_vec, Y_test_fore_vec)
  mean_fore = mean(Y_test_fore_vec)
  std_fore = sd(Y_test_fore_vec)
  y_crps = crps2(Y_test_true_vec, Y_test_fore_vec)
  # all_conf_mat_test[r, , ] = confuse_matrix_NS2(state_test_est_vec, state_test_vec, NS)
  all_y_forecast[r, , ] = Y_test_fore
  all_y_rmse[r] = y_rmse
  all_y_mae[r] = y_mae
  all_y_crps[r] = y_crps
  print(paste('Replication', r, 'has completed'))
}
end_time = Sys.time()
conf_mat_mean = apply(all_conf_mat_test, c(2, 3), mean)

# RMSFE
all_y_rmse
# MAFE
all_y_mae
# CRPS
all_y_crps

for (r in 1:Rep) {
  forecast_est_list <- list(y_rmse = all_y_rmse[r], y_mae = all_y_mae[r], y_crps = all_y_crps[r])
  char_est_list <- lapply(forecast_est_list, function(x) paste(x, collapse = ", "))
  writeLines(unlist(char_est_list), paste0("./code/simulation/Forecast_results_HMM", "/", "forecast_res", r, ".txt"))
}