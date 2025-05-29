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

# enlarge memory for R
memory.limit(size = 200000)
#Work Space of the data
work_dir = 'F:/24561/Documents/ResearchWrite/NHMM-MIDAS/Submission Requirements/Code in Github'
#Set the Work Space of the data
setwd(work_dir)
source('./code/simulation/tool.R')
sourceCpp("./code/simulation/updateHMSP.cpp")

N = 1 # number of subjects
NT = 500 # number of in-sample
NT_oos = 50 # number of out-of-sample
NT_all = NT + NT_oos
y <- rep(1, NT)  
yoos <- rep(1, NT_oos)

NS = 2 # number of hidden states
Rep = 1

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

# Total number of variables
Nbvar <- 30
# Number of relevant variables
Nbrelevant <- 5
# Number of irrelevant variables
Nbirrelevant <- Nbvar - Nbrelevant

beta_real <- matrix(0, nrow = NS, ncol = Nbvar, byrow = TRUE)
beta_real[1, ] = beta <- c(0, 1, -1, -0.8, 1, -1, 2, 0, -2, rep(0, Nbirrelevant - 4))
beta_real[2, ] = beta <- c(-1, 0, -1, 1, 1, 1, -2, 0, 2, rep(0, Nbirrelevant - 4))
b0_real <- c(-2, 2)

Spc = vector('list', 6)
names(Spc) = c('mseries', 'monthly', 'Km', 'kappam', 'sK', 'nbvar')
Spc$mseries = NULL
Spc$monthly = Nbvar
Spc$Km = 24 # lag orders
Spc$kappam = 3
Spc$sK = rep(Spc$Km, Spc$monthly)
NX = Spc$Km * Nbvar + 1 

polydegree = 3 # Almon lag orders
Ra = vector('list',2)
names(Ra) = c('fC','dfC')
Ra$fC = 1 # no tail restriction; if has set 1
Ra$dfC = 1 # no deravative restriction; if has set 1
Q = Almon_lag(polydegree = polydegree, C = Spc$Km, R = Ra)
Qrow = as.numeric(dim(Q)[1])
eg <- Qrow
NX_alm = Qrow * Nbvar

p = eg * Nbvar
NF <- p
ng <- Nbvar
nj <- c(rep(nrow(Q), Nbvar))
groups <- matrix(0, nrow = ng, ncol = nrow(Q)) 
groups[1, ] <- seq_len(nj[1])
if (ng > 1){
  for (i in 2:ng) {
    groups[i, ] <- seq(sum(nj[1:(i - 1)]) + 1, sum(nj[1:(i - 1)]) + nj[i])
  }
}
ggnum <- length(groups[1, ])

all_intercept <- array(0, dim = c(Rep, NS))
all_beta <- array(0, dim = c(Rep, NS, Nbvar))
all_alpha <- array(0, dim = c(Rep, NS, p))
all_sigma <- array(0, dim = c(Rep, NS))
all_state <- array(0, dim = c(Rep, NT))
all_p <- array(0, dim = c(Rep, NS, NS))

for (r in 1:Rep) {
  directory_path <- "./code/simulation/Gen_Data_Train_NS=2"
  specific_number <- r
  pattern <- sprintf("Train_rep%d_chain[0-9]+\\.csv", r)
  files <- list.files(directory_path, pattern = pattern, full.names = TRUE)
  file_names <- sort(files)
  file_names <- mixedsort(file_names, decreasing = FALSE) 
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
    D_train[n, , ] <- as.matrix(data_list[[n]][1:NT, ])
    Y_train[n, ] <- D_train[n, , 1]
    X_train[n, , ] <- D_train[n, , 2:(NX + 1)]
    X_train_alm[n, , ] <- D_train[n, , (NX + 2):(NX + NX_alm + 1)]
    h_train[n, , ] <- D_train[n, , (NX + NX_alm + 2):(1 + NX + NX_alm + NH)]
    state_train[n, ] <- D_train[n, , total_dim]
  }
  
  D_train_mat = array(0, dim = c(N*NT, total_dim))
  for (n in 1:N) {
    D_train_mat[((n-1)*NT+1):(n*NT), ] <- D_train[n, , ] 
  }
  Y_train_mat <- as.matrix(D_train_mat[ ,1])
  X_train_mat <- D_train_mat[ ,2:(NX + 1)]
  X_train_alm_mat <- D_train_mat[ ,(NX + 2):(NX + NX_alm + 1)]
  h_train_mat <- D_train_mat[ ,(NX + NX_alm + 2):(1 + NX + NX_alm + NH)]
  state_train_mat <- as.matrix(D_train_mat[ ,total_dim])
  
  test_path <- "./code/simulation/Gen_Data_Test_NS=2"
  specific_number <- r
  # Regular expression matching files
  # file_names <- list.files(path = directory_path, pattern = paste0("^.*", specific_number, "[^0-9].*\\.csv$"), full.names = TRUE)
  pattern <- sprintf("Test_rep%d_chain[0-9]+\\.csv", r)
  # list all .csv file
  files.test <- list.files(test_path, pattern = pattern, full.names = TRUE)
  # order the .csv file
  file_names.test <- sort(files.test)
  file_names.test <- mixedsort(file_names.test, decreasing = FALSE) 
  # read files
  data_list.test <- lapply(file_names.test, read.csv)
  
  D_test = array(0, dim = c(N, NT_oos, total_dim))
  Y_test = array(0, dim = c(N, NT_oos))
  X_test = array(0, dim = c(N, NT_oos, NX))
  h_test = array(0, dim = c(N, NT_oos, NH))
  X_test_alm = array(0, dim = c(N, NT_oos, NX_alm))
  state_test = array(0, dim = c(N, NT_oos))
  for (n in 1:N) {
    D_test[n, , ] <- as.matrix(data_list.test[[n]][1:NT_oos, ])
    Y_test[n, ] <- D_test[n, , 1]
    X_test[n, , ] <- D_test[n, , 2:(NX + 1)]
    X_test_alm[n, , ] <- D_test[n, , (NX + 2):(NX + NX_alm + 1)]
    h_test[n, , ] <- D_test[n, , (NX + NX_alm + 2):(1 + NX + NX_alm + NH)]
    state_test[n, ] <- D_test[n, , total_dim]
  }
  # state_begin = state_test[ ,1]
  
  D_test_mat = array(0, dim = c(N*NT_oos, total_dim))
  for (n in 1:N) {
    D_test_mat[((n-1)*NT+1):(n*NT_oos), ] <- matrix(unlist(data_list.test[[n]]), ncol = total_dim, byrow = FALSE) # 娉ㄦ剰瑕佸厛灏唋ist鏁版嵁Unlist鍐嶅～鍏ョ煩闃?
  }
  Y_test_mat <- as.matrix(D_test_mat[ ,1])
  X_test_mat <- D_test_mat[ ,2:(NX + 1)]
  X_test_alm_mat <- D_test_mat[ ,(NX + 2):(NX + NX_alm + 1)]
  h_test_mat <- D_test_mat[ ,(NX + NX_alm + 2):(1 + NX + NX_alm + NH)]
  state_test_mat <- as.matrix(D_test_mat[ ,total_dim])
  
  
  # Load Real Parameters
  parameters_path <- "./code/simulation/Paras_NS=2"
  para_names <- sprintf("%s/Paras_rep%d_.txt", parameters_path, r)
  para_names <- mixedsort(para_names, decreasing = FALSE) 
  para_lines <- readLines(para_names)
  # analysis each data
  parsed_data <- lapply(para_lines, function(line) unlist(strsplit(line, ", ")))
  betas_int_real <- matrix(as.numeric(parsed_data[[1]]), nrow = NS, ncol = 1 + Nbvar) 
  betas_real <- betas_int_real[ ,2:(Nbvar + 1)] 
  intercept_real <- betas_int_real[ ,1] 
  sigmas_real <- as.numeric(parsed_data[[2]])
  lop_real <- as.numeric(parsed_data[[3]])
  phi_real <- matrix(as.numeric(parsed_data[[4]]), nrow = 1, ncol = NH)
  zeta_real <- matrix(as.numeric(parsed_data[[5]]), nrow = NS, ncol = NS)
  
  # Standalization for data. 
  x = array(0, dim = c(N, NT, NF))
  y = array(0, dim = c(NT, N))
  for (n in 1:N) {
    # result_x <- standardize_me2(QX)
    result_x <- standardize_me(X_train_alm_mat)
    # result_y <- center_me(as.matrix(Y_train_mat))
    result_y <- Y_train_mat
    x_scale = result_x$y
    mu_x = result_x$mu
    sig_x = result_x$sig
    y_scale = result_y
    # y_scale = result_y$y
    # mu_y = result_y$mu
    x[n, , ] <- x_scale
    y[ ,n] <- y_scale
  }
  
  sig_x_NS = matrix(0, nrow = NS, ncol = ng * ggnum)
  sg <- numeric(ng) # sg stands for the dimension of each parameter theta
  for (ig in 1:ng) {
    sg[ig] <- ggnum
  }
  
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
  
  fc=array(NA,dim=c(N,NT,Nbvar*eg+1))
  fc_NT=array(NA,dim=c(N*NT,Nbvar*eg+1))
  if(Nbvar!=0){
    fc_NT=cbind(x_NT[,1:((Nbvar)*eg)],rep(1,N*NT))
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
  # cal = rep(0, N)
  
  error_count <- 0
  max_errors <- 10
  skip_r <- FALSE
  
  seed1 = as.numeric(Sys.time()) + sample(1:10000, 1, replace = FALSE)
  # seed1 = 1730975314
  set.seed(seed1)
  repeat{
    tryCatch({
      ## initial value of parameters
      alpha_cpp = runif(NS * (eg * (Nbvar) + 1), 0, 2)
      # phi_cpp = as.vector(runif(NS*NH, -1, 1))
      betan_cpp = as.vector(betan)
      lambda2_cpp = as.vector(lambda2)
      tau2_cpp = as.vector(tau2)
      pi0_cpp = as.vector(pi0)
      pi1_cpp = as.vector(pi1)
      Z_cpp = as.vector(Z)
      groups = as.vector(t(groups))
      sigma_cpp = runif(NS, 0, 1)
      s_cpp = sample(0:(NS - 1), N * NT, replace = TRUE)
      # s_cpp = as.vector(real_state) - 1
      y_cpp = y_NT
      x_cpp = as.vector(x_NT)
      fc_cpp = as.vector(fc_NT)
      priorbeta = c(0, 0)
      ele = c(40, 40)
      
      start.time = Sys.time()
      result = mcmc_hmm(y_cpp, x_cpp, fc_cpp, alpha_cpp,
                        sigma_cpp, betan_cpp, lambda2_cpp, tau2_cpp, pi0_cpp, 
                        pi1_cpp, Z_cpp, sg, groups, s_cpp, priorbeta, 
                        N, NT, NS, NF, iter, rep, ntot, ndraw, ng, eg, ele)
      end.time = Sys.time()
      
      # Result Analysis
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
      p.estimation = apply(p_result, c(2:3), mean)
      sigma.sd = apply(sigma_result, 1, sd)
      alpha.sd = apply(alpha_result, c(2:3), sd)
      intercept.sd = alpha.sd[ , ncol(alpha.estimation)]
      
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
      # sum(c(state_train) == state.estimation)
      # state.estimation
      sig_x_NS = matrix(0, nrow = NS, ncol = Nbvar * eg)
      mu_x_NS = matrix(0, nrow = NS, ncol = Nbvar * eg)
      beta.estimation = alpha.estimation[ ,-ncol(alpha.estimation)]
      beta.real.estimation = matrix(0, nrow = NS, ncol = Nbvar)
      beta.real.estimation2 = matrix(0, nrow = NS, ncol = Nbvar)
      alpha.no.int = alpha.estimation[ ,1:(Nbvar * eg)]
      plast_prob <- temp.mat[nrow(temp.mat), ] / sum(temp.mat[nrow(temp.mat), ])
      for (s in 1:NS) {
        sig_x_NS[s, ] = standardize_me2(X_train_alm_mat[state.estimation == s, ])$sig
        # mu_x_NS[s, ] = standardize_me2(X_train_alm_mat[state.estimation == s, ])$mu
        for (nb in 1:Nbvar) {
          beta.real.estimation[s, nb] = sum((alpha.estimation[s, ((nb - 1)*eg + 1):(nb*eg)] / sig_x_NS[s, ((nb - 1)*eg + 1):(nb*eg)]) %*% Q)
          beta.real.estimation2[s, nb] = sum((alpha.estimation[s, ((nb - 1)*eg + 1):(nb*eg)] / sig_x[((nb - 1)*eg + 1):(nb*eg)]) %*% Q)
        }
      }
      
      all_intercept[r, ] <- intercept.estimation
      all_beta[r, , ] <- beta.real.estimation
      all_alpha[r, , ] <- alpha.estimation[ ,-ncol(alpha.estimation)]
      all_sigma[r, ] <- sigma.estimation
      all_p[r, , ] <- p.estimation
      all_state[r, ] <- state.estimation
      print(paste0("\r Trial times", r, "Running MCMC, ", "has completed..."))
      
      specific_number <- r
      paras_est_list <- list(beta_est = all_beta[r, , ], beta_org_est = all_alpha[r, , ],
                             intercept_est = all_intercept[r, ], sigma_est = all_sigma[r, ],
                             p_est = all_p[r, , ], last_prob = plast_prob)
      char_est_list <- lapply(paras_est_list, function(x) paste(x, collapse = ", "))
      writeLines(unlist(char_est_list), paste0("./code/simulation/Paras_estimation_HMM", "/", "Paras_est_rep", specific_number, ".txt"))
      
      break
    }, error = function(e) {
      # change seed
      message("error catching: ", e$message)
      error_count <<-  error_count + 1
      print(error_count)
      new_seed <- as.numeric(Sys.time()) + sample(1:10000, 1, replace = FALSE)
      set.seed(new_seed)
      message("change seed", new_seed, " Retry ...")
      # repeat 
      
      if (error_count >= max_errors) {
        skip_r <<- TRUE
        # cat("  error has get limited, skip r =", r, "\n\n")
        # break
      }
    })
    
    if (skip_r) {
      cat("r =", r, " error has get limited then skip this cycle\n\n")
      break  # break out of tryCatch
    }
  }
  
}

