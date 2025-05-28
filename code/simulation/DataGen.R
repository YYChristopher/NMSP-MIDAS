setwd('F:/24561/Documents/ResearchWrite/NHMM-MIDAS/Submission Requirements/Code in Github')
# rm(list=ls())
source('./code/simulation/tool.R')
# Generation: The weight function corresponding to X in different hidden states is different

N = 1 # number of subjects
NT = 500 # number of training set
NT_oos = 50 # number of out-of-sample
NT_all = NT + NT_oos
y <- rep(1, NT)
yoos <- rep(1, NT_oos)

NS = 2 # number of hidden states
NH = 3 # external variables in transition probability model
Rep = 1

# Total number of variables
Nbvar <- 30
# Number of relevant variables
Nbrelevant <- 5
# Number of irrelevant variables
Nbirrelevant <- Nbvar - Nbrelevant

beta_real <- matrix(0, nrow = NS, ncol = Nbvar, byrow = TRUE)
beta_real[1, ] = beta <- c(0, 1, -1, -0.8, 1, -1, 2, 0, -2, rep(0, Nbirrelevant - 4))
beta_real[2, ] = beta <- c(-1, 0, -1, 1, 1, 1, -2, 0, 2, rep(0, Nbirrelevant - 4))
# Intercept
b0_real <- c(-2, 2)

# Basic specifications for explanatory variables (monthly)
Spc = vector('list', 6)
names(Spc) = c('mseries', 'monthly', 'Km', 'kappam', 'sK', 'nbvar')
Spc$mseries = NULL
Spc$monthly = Nbvar
Spc$Km = 24 # lag orders. 
Spc$kappam = 3 # correspond to m_k
Spc$sK = rep(Spc$Km, Spc$monthly)

# Weighting scheme and parameters under the DGP MIDAS
Spc$nbvar = Spc$monthly
# slow decaying
theta_ele_11 = rep(0.0007, Spc$nbvar)
theta_ele_12 = rep(-0.006, Spc$nbvar)
theta_ele_21 = rep(0.0007, Spc$nbvar)
theta_ele_22 = rep(-0.006, Spc$nbvar)
theta_NS = matrix(c(0.0007, -0.006, 0.0007, -0.006), nrow = NS, ncol = 2, byrow = TRUE)

# Linearizing the MIDAS weighting scheme
W_NS <- matrix(0, nrow = NS, ncol = Nbvar * Spc$Km + 1)
for (s in 1:NS) {
  theta = c(rep(theta_NS[s, 1], Spc$nbvar), rep(theta_NS[s, 2], Spc$nbvar))
  result_WM <- weights_midas(theta, beta_real[s, ], Spc)
  WM <- result_WM$WM
  weiii <- result_WM$ww
  W <- c(WM, b0_real[s])
  W_NS[s, ] <- W
}

# Constrained Almon lag polynomial (monthly)
polydegree = 3 # Almon lag orders
Ra = vector('list',2)
names(Ra) = c('fC','dfC')
Ra$fC = 1 # no tail restriction; if has set 1
Ra$dfC = 1 # no deravative restriction; if has set 1
Q = Almon_lag(polydegree = polydegree, C = Spc$Km, R = Ra)
Qrow = as.numeric(dim(Q)[1])

#---------------------Real Parameters------------------------
#---------------------transition model------------------------
lop_real <- array(0, dim = c(Rep, NS - 1))
phi_real <- array(0, dim = c(Rep, NH))
zeta_real <- array(0, dim = c(Rep, NS, NS))
#---------------------Emission model------------------------
betas_real <- array(0, dim = c(Rep, NS, 1 + Nbvar))
sigmas_real <- array(0, dim = c(Rep, NS))

#--------------------------About data storage----------------------------------
all_t_state_train <- array(0, dim = c(Rep, N, NT))
all_t_state_test <- array(0, dim = c(Rep, N, NT_oos))
all_x_train <- array(0, dim = c(Rep, N, NT, (Spc$Km * Nbvar) + 1))
all_x_test <- array(0, dim = c(Rep, N, NT_oos, (Spc$Km * Nbvar) + 1))
all_x_train_alm <- array(0, dim = c(Rep, N, NT, Spc$nbvar * Qrow))
all_x_test_alm <- array(0, dim = c(Rep, N, NT_oos, Spc$nbvar * Qrow))
all_h_train <- array(0, dim = c(Rep, N, NT, NH))
all_h_test <- array(0, dim = c(Rep, N, NT_oos, NH))
all_y_train <- array(0, dim = c(Rep, N, NT))
all_y_test <- array(0, dim = c(Rep, N, NT_oos))

set.seed(as.numeric(Sys.time()) + sample(1:10000, 1, replace = FALSE))
for (r in 1:Rep) {
  
  for (n in 1:N) {
    Sigma <- matrix(1, nrow = Nbvar, ncol = Nbvar)
    SS <- 0.5
    
    for (ii in 1:Nbvar) {
      for (iii in 1:Nbvar) {
        if (ii != iii) {
          Sigma[ii, iii] <- SS^(abs(ii - iii))
        }
      }
    }
    
    eps <- mvrnorm(2 * (NT + NT_oos) * Spc$kappam, mu = rep(0, Nbvar), Sigma = Sigma)
    Xm <- matrix(0, nrow = 2 * (NT + NT_oos) * Spc$kappam, ncol = Nbvar)
    
    for (iii in 1:Nbvar) {
      delta <- 0.5
      for (ii in 2:(2 * (NT + NT_oos) * Spc$kappam)) {
        Xm[ii, iii] <- delta * Xm[ii - 1, iii] + eps[ii, iii]
      }
    }
    
    Xm <- Xm[11:(2 * (NT + NT_oos) * Spc$kappam), ]
    
    # Constructing X matrices in the MIDAS form. XX_Reg对应训练集，XX_For对应测试???
    XX_Reg <- matrix(0, nrow = NT, ncol = 0)
    XX_For <- matrix(0, nrow = NT_oos, ncol = 0)
    
    # Monthly matrix
    x <- vector("list", length = Spc$monthly)
    for (ii in 1:Spc$monthly) {
      x[[ii]] <- Construct_DataMIDAS(g = c(y, yoos), d = Xm[, ii], K = Spc$Km, m = Spc$kappam)
      XX_Reg <- cbind(XX_Reg, x[[ii]][1:NT, ])
      XX_For <- cbind(XX_For, x[[ii]][(NT + 1):(NT + NT_oos), ])
    }
    
    # Adding a column of one in the covariates matrix
    if (!is.null(b0_real)) {
      XX_Reg <- cbind(XX_Reg, rep(1, NT))
      XX_For <- cbind(XX_For, rep(1, NT_oos))
    }
    
    all_x_train[r, n, , ] <- XX_Reg
    all_x_test[r, n, , ] <- XX_For
    
    x <- vector("list", length = Spc$monthly)
    
    # Prepare data for X with transform matrix Q. 
    XX_Reg_alm <- matrix(0, nrow = dim(Q)[1] * Nbvar, ncol = NT)
    XX_For_alm <- matrix(0, nrow = dim(Q)[1] * Nbvar, ncol = NT_oos)
    
    for (ii in 1:Spc$monthly) {
      x[[ii]] <- Construct_DataMIDAS(g = c(y, yoos), d = Xm[, ii], K = Spc$Km, m = Spc$kappam)
      
      # Data for MIDAS with Almon Lag
      start_row <- dim(Q)[1] * (ii - 1) + 1
      end_row <- dim(Q)[1] * ii
      
      XX_Reg_alm[start_row:end_row, ] <- Q %*% t(x[[ii]][1:NT, ])
      XX_For_alm[start_row:end_row, ] <- Q %*% t(x[[ii]][(NT + 1):(NT + NT_oos), ])
    }
    
    XX_Reg_alm = t(XX_Reg_alm)
    XX_For_alm = t(XX_For_alm)
    all_x_train_alm[r, n, , ] <- XX_Reg_alm
    all_x_test_alm[r, n, , ] <- XX_For_alm
  }
  
  # h in transition probability model
  h_train <- array(0, dim = c(N, NT, NH))
  h_train[ , ,1] <- matrix(rnorm(N * NT, mean = 0, sd = 1), nrow = N, ncol = NT)
  h_train[ , ,2] <- matrix(runif(N * NT, min = 0, max = 1), nrow = N, ncol = NT)
  h_train[ , ,3] <- matrix(rbinom(N * NT, size = 1, prob = 0.9), nrow = N, ncol = NT)
  all_h_train[r, , , ] <- h_train
  
  # t_lop <- rep(-1, NS - 1)
  t_lop <- c(0)
  t_p0 <- p0_gen(t_lop)
  # define t_zeta和t_phi
  t_zeta <- matrix(c(-1, 0, 1, 0), nrow = NS, ncol = NS, byrow = TRUE) 
  t_phi <- c(1, -1, 0.5)
  t_tran_train <- tran(t_zeta, t_phi, h_train, NS)
  t_tran_p_train <- tran_p(t_tran_train, NS)
  t_state_train <- tran_state(t_tran_p_train, t_p0)
  all_t_state_train[r, , ] <- t_state_train
  state1_pro_train <- sum(t_state_train == 1) / (nrow(t_state_train) * ncol(t_state_train))
  state2_pro_train <- sum(t_state_train == 2) / (nrow(t_state_train) * ncol(t_state_train))
  Y_train <- array(dim = c(N, NT))
  location1_train <- which(t_state_train == 1, arr.ind = TRUE)
  location2_train <- which(t_state_train == 2, arr.ind = TRUE) 
  t_mean_train <- array(dim = c(N, NT))
  # DGP from Gaussian noise. Set the signal-to-noise ratio. 
  SNR <- 5
  sigma_NS <- rep(0, NS)
  noise <- rnorm(NT)
  # signal_NS <- XX_Reg %*% t(W_NS)
  # for (s in 1:NS) {
  #   sigma_NS[s] <- sqrt(var(signal_NS[s, ]) / (SNR * var(noise)))
  # }
  
  for (i in 1:nrow(location1_train)) {
    idx <- location1_train[i,]
    t_mean_train[idx[1], idx[2]] <- all_x_train[r, idx[1], idx[2], ] %*% W_NS[1, ]
  }
  
  for (i in 1:nrow(location2_train)) {
    idx <- location2_train[i,]
    t_mean_train[idx[1], idx[2]] <- all_x_train[r, idx[1], idx[2], ] %*% W_NS[2, ]
  }
  
  t_mean_train_vec = sample_merge_Ystate(t_mean_train)
  t_state_train_vec = sample_merge_Ystate(t_state_train)
  for (s in 1:NS) {
    sigma_NS[s] <- sqrt(var(t_mean_train_vec[t_state_train_vec == s]) / (SNR * var(noise)))
  }
  
  for (i in 1:nrow(location1_train)) {
    idx <- location1_train[i,]
    Y_train[idx[1], idx[2]] <- rnorm(1, mean = t_mean_train[idx[1], idx[2]], sd = sigma_NS[1])
  }
  
  for (i in 1:nrow(location2_train)) {
    idx <- location2_train[i,]
    Y_train[idx[1], idx[2]] <- rnorm(1, mean = t_mean_train[idx[1], idx[2]], sd = sigma_NS[2])
  }
  
  all_y_train[r, , ] <- Y_train
  
  
  # testing set
  h_test <- array(0, dim = c(N, NT_oos, NH))
  h_test[ , ,1] <- matrix(rnorm(N * NT_oos, mean = 0, sd = 1), nrow = N, ncol = NT_oos)
  h_test[ , ,2] <- matrix(runif(N * NT_oos, min = 0, max = 1), nrow = N, ncol = NT_oos)
  h_test[ , ,3] <- matrix(rbinom(N * NT_oos, size = 1, prob = 0.9), nrow = N, ncol = NT_oos)
  all_h_test[r, , , ] <- h_test
  
  t_tran_test <- tran_acc(t_zeta, t_phi, h_test, NS)
  t_tran_p_test <- tran_p(t_tran_test, NS)
  t_state_test <- tran_state(t_tran_p_test, t_p0)
  all_t_state_test[r, , ] <- t_state_test
  state1_pro_test <- sum(t_state_test == 1) / (nrow(t_state_test) * ncol(t_state_test))
  state2_pro_test <- sum(t_state_test == 2) / (nrow(t_state_test) * ncol(t_state_test))
  
  Y_test <- array(dim = c(N, NT_oos))
  location1_test <- which(t_state_test == 1, arr.ind = TRUE)
  location2_test <- which(t_state_test == 2, arr.ind = TRUE)
  t_mean_test <- array(dim = c(N, NT_oos))
  
  for (i in 1:nrow(location1_test)) {
    idx <- location1_test[i,]
    t_mean_test[idx[1], idx[2]] <- all_x_test[r, idx[1], idx[2], ] %*% W_NS[1, ]
    Y_test[idx[1], idx[2]] <- rnorm(1, mean = t_mean_test[idx[1], idx[2]], sd = sigma_NS[1])
  }
  
  for (i in 1:nrow(location2_test)) {
    idx <- location2_test[i,]
    t_mean_test[idx[1], idx[2]] <- all_x_test[r, idx[1], idx[2], ] %*% W_NS[2, ]
    Y_test[idx[1], idx[2]] <- rnorm(1, mean = t_mean_test[idx[1], idx[2]], sd = sigma_NS[2])
  }
  
  all_y_test[r, , ] <- Y_test
  
  for (s in 1:NS) {
    betas_real[r, s, ] <- c(b0_real[s], beta_real[s, ])
    sigmas_real[r, ] <- sigma_NS
  }
  lop_real[r, ] <- t_lop
  phi_real[r, ] <- t_phi
  zeta_real[r, , ] <- t_zeta
  
}

# save data
for (r in 1:Rep) {
  for (n in 1:N) {
    D_train <- cbind(all_y_train[r, n, ], all_x_train[r, n, , ], all_x_train_alm[r, n, , ], all_h_train[r, n, , ], all_t_state_train[r, n, ])
    write.csv(D_train, file = paste0("./code/simulation/Gen_Data_Train_NS=2", "/", "Train_rep", r, "_", "chain", n, ".csv"), row.names = FALSE)
  }
}

for (r in 1:Rep) {
  for (n in 1:N) {
    D_test <- cbind(all_y_test[r, n, ], all_x_test[r, n, , ], all_x_test_alm[r, n, , ], all_h_test[r, n, , ], all_t_state_test[r, n, ])
    write.csv(D_test, file = paste0("./code/simulation/Gen_Data_Test_NS=2", "/", "Test_rep", r, "_", "chain", n, ".csv"), row.names = FALSE)
  }
}

for (r in 1:Rep) {
  paras_list <- list(betas_real = betas_real[r, , ], sigmas_real = (sigmas_real[r, ])^2, lop_real = lop_real[r, ],
                     phi_real = phi_real[r, ], zeta_real = zeta_real[r, , ])
  char_list <- lapply(paras_list, function(x) paste(x, collapse = ", "))
  writeLines(unlist(char_list), paste0("./code/simulation/Paras_NS=2", "/", "Paras_rep", r, "_.txt"))
}

# print(c(state1_pro_train, state2_pro_train))
