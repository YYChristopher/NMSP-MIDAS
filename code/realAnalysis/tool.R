library(MASS)
library(Matrix)
library(statmod)
library(scoringRules)
library(ggplot2)
library(tidyverse)
library(tibble)
library(forecast)

center_me <- function(x) {
  # Function to make your data have mean 0.
  # Data in x are T \times p, i.e. T time series observations times p variables
  
  mu <- colMeans(x, na.rm = TRUE)
  sig <- apply(x, 2, sd, na.rm = TRUE)
  
  y <- x - matrix(rep(mu, each = nrow(x)), nrow = nrow(x), byrow = FALSE) #注意byrow应为False
  
  return(list(y = y, mu = mu, sig = sig))
}


standardize_me <- function(x) {
  # Function to make your data have mean 0 and variance 1.
  # Data in x are T times p, i.e. T time series observations times p variables
  
  mu <- colMeans(x, na.rm = TRUE)
  sig <- apply(x, 2, sd, na.rm = TRUE)
  
  y <- (x - matrix(rep(mu, each = nrow(x)), nrow = nrow(x), byrow = FALSE)) / 
    matrix(rep(sig, each = nrow(x)), nrow = nrow(x), byrow = FALSE)
  
  return(list(y = y, mu = mu, sig = sig))
}


standardize_me2 <- function(x) {
  # Function to make your data have mean 0 and variance 1.
  # Data in x are T times p, i.e. T time series observations times p variables
  
  mu <- colMeans(x, na.rm = TRUE)
  sig <- apply(x, 2, sd, na.rm = TRUE)
  
  y <- x / matrix(rep(sig, each = nrow(x)), nrow = nrow(x), byrow = FALSE)
  
  return(list(y = y, mu = mu, sig = sig))
}


# Almon lag polynomial
Almon_lag <- function(polydegree, C, R) {
  # R represents restrictions on lag orders()
  # fC stands for tail restriction and dfC stands for derivative restrictions
  if (polydegree == 2) {
    if (R$fC == 0 && R$dfC == 0) {
      # Unrestricted Almon lag polynomial
      Q <- matrix(0, nrow = polydegree + 1, ncol = C)
      for (ii in 0:polydegree) {
        Q[ii + 1, ] <- (0:(C - 1))^ii
      }
    } else if (R$fC == 1 && R$dfC == 0) {
      # Almon lag polynomial with tail restriction (fC=0)
      Q <- matrix(0, nrow = polydegree, ncol = C)
      for (ii in 1:polydegree) {
        Q[ii, ] <- (0:(C - 1))^ii - (C - 1)^ii
      }
    } else if (R$fC == 0 && R$dfC == 1) {
      # Almon lag polynomial with derivative restriction (dfC=0)
      Q <- matrix(0, nrow = polydegree, ncol = C)
      for (ii in 1:polydegree) {
        Q[ii, ] <- ((0:(C - 1))^ii - ii * (C - 1) * (0:(C - 1)))^(ii - 1)
      }
    } else if (R$fC == 1 && R$dfC == 1) {
      # Almon lag polynomial with tail and derivative restrictions (fC=0 and dfC=0)
      Q <- matrix(0, nrow = polydegree - 1, ncol = C)
      for (ii in 1:(polydegree - 1)) {
        Q[ii, ] <- (0:(C - 1))^(ii + 1) - (ii + 1) * (C - 1) * (0:(C - 1)) + (C - 1)^(ii + 1)
      }
    }
  } else if (polydegree == 3) {
    if (R$fC == 0 && R$dfC == 0) {
      # Unrestricted Almon lag polynomial
      Q <- matrix(0, nrow = polydegree + 1, ncol = C)
      for (ii in 0:polydegree) {
        Q[ii + 1, ] <- (0:(C - 1))^ii
      }
    } else if (R$fC == 1 && R$dfC == 0) {
      # Almon lag polynomial with tail restriction (fC=0)
      Q <- matrix(0, nrow = polydegree, ncol = C)
      for (ii in 1:polydegree) {
        Q[ii, ] <- (0:(C - 1))^ii - (C - 1)^ii
      }
    } else if (R$fC == 0 && R$dfC == 1) {
      # Almon lag polynomial with derivative restriction (dfC=0)
      stop('Restrictions on the lag polynomial not allowed')
    } else if (R$fC == 1 && R$dfC == 1) {
      # Almon lag polynomial with tail and derivative restrictions (fC=0 and dfC=0)
      Q <- matrix(0, nrow = polydegree - 1, ncol = C)
      for (ii in 1:(polydegree - 1)) {
        Q[ii, ] <- (0:(C - 1))^(ii + 1) - (ii + 1) * (C - 1)^ii * (0:(C - 1)) + ii * (C - 1)^(ii + 1)
      }
    }
  }
  
  return(Q)
}


randig <- function(mu, lamb) {
  # Generate a draw from the Inverse Gaussian distribution
  # can also utilize rinvgauss
  # Inverse Gaussian distribution:
  # pdf(x) = sqrt(lamb/(2*pi*x^3)) * exp(-lamb./(2*x).*(x/mu-1).^2);
  # cdf(x) = normcdf(sqrt(lamb./x).*(x/mu-1)) + exp(2*lamb/mu)*normcdf(sqrt(lamb./x).*(-x/mu-1));
  
  sampleSize <- length(mu)
  
  lambsq1 <- rnorm(sampleSize)^2
  out <- mu + (0.5 * mu / lamb) * 
    (mu * lambsq1 - sqrt(4 * mu * lamb * lambsq1 + (mu^2) * (lambsq1^2)))
  
  l <- which(runif(sampleSize) >= mu / (mu + out))
  out[l] <- (mu[l]^2) / out[l]
  
  return(out)
}


weights_midas <- function(theta, beta, Spc) {
  # Construction of covariates matrix as defined by MIDAS weighting scheme
  # theta: parameters theta for the weighting kernel
  # beta: parameters beta for regression coefficients
  # Spc: MIDAS specifications
  
  l <- list()
  for (i in 1:Spc$monthly) {
    l[[i]] <- list(
      one = rep(1, Spc$Km),
      w = (1:Spc$Km) / Spc$Km,
      k = 1:Spc$Km,
      kk = Spc$Km
    )
  }
  
  # Two-parameters exp Almon lag
  theta1 <- theta[1:Spc$nbvar]
  theta2 <- theta[(Spc$nbvar + 1):(2 * Spc$nbvar)]
  
  # WW <- matrix(0, nrow = Spc$Km, ncol = length(theta2))
  # ww <- matrix(0, nrow = Spc$Km, ncol = length(theta2))
  WW <- rep(0, Spc$Km * length(theta2))
  ww <- rep(0, Spc$Km * length(theta2))
  
  for (i in 1:length(theta2)) {
    W <- exp(theta1[i] * l[[i]]$k + theta2[i] * (l[[i]]$k^2)) / sum(exp(theta1[i] * l[[i]]$k + theta2[i] * (l[[i]]$k^2)))
    
    WW[((i - 1)*Spc$Km + 1):(i*Spc$Km)] <- W * beta[i]
    ww[((i - 1)*Spc$Km + 1):(i*Spc$Km)] <- W
  }
  
  WM <- t(WW) 
  
  return(list(WM = WM, ww = ww))
}


exp_almon_lag <- function(theta, K){
  
  Q <- length(theta)
  res.mid = c(rep(0, K))
  for (j in 1:K) {
    temp = 0
    for (q in 1:Q) {
      temp = temp + theta[q] * j^q
    }
    res.mid[j] = exp(temp)
  }
  res <- res.mid / sum(res.mid)
  return(res)
}


almon_lag <- function(theta, K){
  
  Q <- length(theta)
  res.mid = c(rep(0, K))
  for (j in 1:K) {
    temp = 0
    for (q in 1:Q) {
      temp = temp + theta[q] * j^q
    }
    res.mid[j] = temp
  }
  res <- res.mid / sum(res.mid)
  return(res)
}


Construct_DataMIDAS <- function(g, d, K, m) {
  # g: vector of low frequency data (Y)
  # d: vector of high frequency data (X)
  # (d must end at exactly the same data as the low frequency time series
  # e.g., if g is quarterly and ends at 2012q2, and d is monthly sampled, d must end on June 2012.)
  # K: number of high frequency data to involve in the MIDAS weighting scheme (K in MIDAS equations)
  # m: difference in frequency between g and d (e.g., if g is quarterly and d is monthly, m=3)
  
  lo <- length(g)
  i <- length(d)
  D <- matrix(0, nrow = lo, ncol = K)
  
  for (row in 0:(lo - 1)) {
    if (i > K) {
      for (col in 0:(K - 1)) {
        D[lo - row, col + 1] <- d[i - col]
      }
      i <- i - m
    }
  }
  
  return(D)
}


p0_gen <- function(tau) {
  # transition probability model
  NS <- length(tau) + 1
  p0 <- rep(1, NS)
  
  p0[1] <- exp(tau[1]) / (1 + exp(tau[1]))
  
  if (NS > 2) {
    for (i in 2:(NS - 1)) {
      p0[i] <- 1 / prod(1 + exp(tau[1:(i - 1)]))
      p0[i] <- p0[i] * (exp(tau[i]) / (1 + exp(tau[i])))
    }
  }
  
  p0[NS] <- 1 / prod(1 + exp(tau))
  return(p0)
}


# The following two functions are combined to calculate the transition probability matrix 
# at the NT moment of the Nth sample chain:
tran <- function(zeta, phi, d, NS) {
  # get dimension of d
  N <- dim(d)[1]
  NT <- dim(d)[2]
  
  # initialize tran array
  tran <- array(0, dim = c(N, NT, NS, NS))
  
  # calculate phid
  # phid <- apply(d, c(1, 2), function(x) sum(phi * x))
  phid <- array(0, dim = c(N, NT))
  for (n in 1:N) {
    for (t in 1:NT) {
      # phid[n, t] = t(d[n, t, ]) %*% phi
      phid[n, t] = sum(d[n, t, ] * phi)
    }
  }
  
  for (u in 1:NS) {
    for (s in 1:NS) {
      tran[ , , u, s] <- zeta[u, s] + phid
    }
  }
  
  return(tran)
}

tran_acc <- function(zeta, phi, d, NS) {
  # This function is used to calculate the probability transfer matrix for the 
  # state-specific version of phi.
  N <- dim(d)[1]
  NT <- dim(d)[2]
  
  # initialize tran array
  tran <- array(0, dim = c(N, NT, NS, NS))
  # calculate phid
  # phid <- apply(d, c(1, 2), function(x) sum(phi * x))
  phid <- array(0, dim = c(N, NT, NS, NS))
  for (n in 1:N) {
    for (t in 1:NT) {
      phid[n, t, , ] = cbind(matrix(d[n, t, ] %*% t(phi), nrow = NS, ncol = (NS - 1), byrow = TRUE), rep(0, NS))
    }
  }
  # fill tran array
  for (u in 1:NS) {
    for (s in 1:NS) {
      tran[ , , u, s] <- zeta[u, s] + phid[ , , u, s]
    }
  }
  return(tran)
}


tran_acc2 <- function(zeta, phi, d, NS) {
  # This function is used to calculate the probability transfer matrix for the 
  # state-specific version of phi.
  N <- dim(d)[1]
  NT <- dim(d)[2]
  
  # initialize tran array
  tran <- array(0, dim = c(N, NT, NS, NS))
  # calculate phid
  # phid <- apply(d, c(1, 2), function(x) sum(phi * x))
  phid <- array(0, dim = c(N, NT, NS, NS))
  for (n in 1:N) {
    for (t in 1:NT) {
      phid[n, t, , ] = matrix(d[n, t, ] %*% t(phi), nrow = NS, ncol = NS, byrow = FALSE)
    }
  }
  # fill tran array
  for (u in 1:NS) {
    for (s in 1:NS) {
      tran[ , , u, s] <- zeta[u, s] + phid[ , , u, s]
    }
  }
  return(tran)
}


tran_p <- function(tran, NS) {
  # get dimension of tran
  N <- dim(tran)[1]
  NT <- dim(tran)[2]
  
  # initialize p_tran
  p_tran <- array(1, dim = c(N, NT, NS, NS))
  
  for (u in 1:NS) {
    p_tran[ , , u, 1] <- exp(tran[ , , u, 1]) / (1 + exp(tran[ , , u, 1]))
    
    if(N == 1 && NS == 2){
      p_tran[ , , u, NS] <- 1 / apply(matrix(1 + exp(tran[ , , u, 1:(NS - 1)]), nrow = N), c(1, 2), prod)
    } else if(N == 1 && NS > 2) {
      p_tran[ , , u, NS] <- 1 / apply(1 + exp(tran[ , , u, 1:(NS - 1)]), 1, prod)
    } else {
      p_tran[ , , u, NS] <- 1 / apply(1 + exp(tran[ , , u, 1:(NS - 1)]), c(1, 2), prod)
    }
    
    if (NS > 2) {
      if (N > 1){
        for (s in 2:(NS - 1)) {
          p_tran[ , , u, s] <- (1 / apply(1 + exp(tran[ , , u, 1:(s - 1)]), c(1, 2), prod)) * (exp(tran[ , , u, s]) / (1 + exp(tran[ , , u, s])))
          # p_tran[ , , u, s] <- p_tran[ , , u, s] * (exp(tran[ , , u, s]) / (1 + exp(tran[ , , u, s])))
        }
      } else {
        for (s in 2:(NS - 1)) {
          # p_tran[ , , u, s] <- (1 / apply(1 + exp(tran[ , , u, 1:(s - 1)]), 1, prod)) * (exp(tran[ , , u, s]) / (1 + exp(tran[ , , u, s])))
          p_tran[ , , u, s] <- (1 / apply(matrix(1 + exp(tran[ , , u, 1:(s - 1)]), nrow = NT), 1, prod)) * (exp(tran[ , , u, s]) / (1 + exp(tran[ , , u, s])))
          # p_tran[ , , u, s] <- p_tran[ , , u, s] * (exp(tran[ , , u, s]) / (1 + exp(tran[ , , u, s])))
        }
      }
    }
  }
  
  return(p_tran)
}


tran_p_acc <- function(tran, NS) {
  # get dimension of tran
  N <- dim(tran)[1]
  NT <- dim(tran)[2]
  
  # initialize p_tran
  p_tran <- array(0, dim = c(N, NT, NS, NS))
  p_tran_denom <- apply(exp(tran), c(1, 2, 3), sum)
  for (u in 1:NS) {
    p_tran[ , , u, NS] <- 1 / p_tran_denom[ , , u]
    for (v in 1:(NS - 1)) {
      p_tran[ , , u, v] <- exp(tran[ , , u, v]) / p_tran_denom[ , , u]
    }
  }
  
  return(p_tran)
}


# Hidden state sampling based on transition probability matrix
tran_state <- function(tran_p, p0) {
  # get dimension of tran_p
  N <- dim(tran_p)[1]
  NT <- dim(tran_p)[2]
  NS <- dim(tran_p)[3]
  
  # initialize state array
  state <- matrix(0, nrow = N, ncol = NT)
  
  for (i in 1:N) {
    state[i, 1] <- sample(1:NS, 1, prob = p0)
  }
  
  for (i in 1:N) {
    for (t in 2:NT) {
      prev_state <- state[i, t - 1]
      state[i, t] <- sample(1:NS, 1, prob = tran_p[i, t, prev_state, ])
    }
  }
  
  return(state)
}


tran_state_new <- function(tran_p, initial_state) {
  # get dimension of tran_p
  N <- dim(tran_p)[1]
  NT <- dim(tran_p)[2]
  NS <- dim(tran_p)[3]
  
  # initialize state array
  state <- matrix(0, nrow = N, ncol = NT)
  
  for (i in 1:N) {
    state[i, 1] <- initial_state[i]
  }
  
  for (i in 1:N) {
    for (t in 2:NT) {
      prev_state <- state[i, t - 1]
      state[i, t] <- sample(1:NS, 1, prob = tran_p[i, t, prev_state, ])
      # state[i, t] <- which.max(tran_p[i, t, prev_state, ])
    }
  }
  
  return(state)
}


tran_state_hmm <- function(tran_p, initial_state, NT) {
  # get dimension of tran_p
  NS <- dim(tran_p)[1]
  
  # initialize state array
  state <- matrix(0, nrow = N, ncol = NT)
  
  for (i in 1:N) {
    state[i, 1] <- initial_state
  }
  
  for (i in 1:N) {
    for (t in 2:NT) {
      prev_state <- state[i, t - 1]
      state[i, t] <- sample(1:NS, 1, prob = tran_p[prev_state, ])
      # state[i, t] <- which.max(tran_p[i, t, prev_state, ])
    }
  }
  
  return(state)
}


# The following sample_merge function is used to convert (N, NT)-dimensional data into 
# (N*NT)-dimensional data, that is, to compress a sample dimension.
sample_merge_Xh <- function(x){
  
  # First take out the corresponding dimension
  N = as.numeric(dim(x)[1])
  NT = as.numeric(dim(x)[2])
  tot_dim = as.numeric(dim(x)[3])
  X_mat = array(0, dim = c(N*NT, tot_dim))
  for (n in 1:N) {
    X_mat[((n-1)*NT+1):(n*NT), ] <- matrix(x[n, , ], ncol = tot_dim, byrow = FALSE) 
  }
  X_mat <- as.matrix(X_mat)
  return(X_mat)
}


# This function is used for data Y and state with tot_dim == 1
sample_merge_Ystate <- function(y){
  
  # First take out the corresponding dimension
  N = as.numeric(dim(y)[1])
  NT = as.numeric(dim(y)[2])
  Y_mat = array(0, dim = c(N*NT))
  for (n in 1:N) {
    Y_mat[((n-1)*NT+1):(n*NT)] <- matrix(y[n, ], nrow = NT, byrow = FALSE) 
  }
  Y_mat <- as.matrix(Y_mat)
  return(Y_mat)
}


sample_split_Xh <- function(X_mat, N, NT){
  
  tot_dim = as.numeric(dim(X_mat)[2])
  X = array(0, dim = c(N, NT, tot_dim))
  for (n in 1:N) {
    X[n, , ] <- X_mat[((n-1)*NT+1):(n*NT), ]
  }
  return(X)
}


sample_split_Ystate <- function(Y_mat, N, NT){
  
  Y = array(0, dim = c(N, NT))
  for (n in 1:N) {
    Y[n, ] <- Y_mat[((n-1)*NT+1):(n*NT)]
  }
  return(Y)
}


## Some evaluation
# confuse matrix
confuse_matrix <- function(state_vec, true_state_vec, NS){
  
  conf_mat <- matrix(0, nrow = NS, ncol = NS)
  for (i in 1:length(state_vec)) {
    row <- true_state_vec[i]
    col <- state_vec[i]
    conf_mat[row, col] <- conf_mat[row, col] + 1
  }
  df <- conf_mat
  colnames(df) <- c('est state1', 'est state2', 'est state3')
  rownames(df) <- c('real state1', 'real state2', 'real state3')
  return(df)
}


confuse_matrix_NS2 <- function(state_vec, true_state_vec, NS){
  
  conf_mat <- matrix(0, nrow = NS, ncol = NS)
  for (i in 1:length(state_vec)) {
    row <- true_state_vec[i]
    col <- state_vec[i]
    conf_mat[row, col] <- conf_mat[row, col] + 1
  }
  df <- conf_mat
  colnames(df) <- c('est state1', 'est state2')
  rownames(df) <- c('real state1', 'real state2')
  return(df)
}


plot_confusion <- function(cm, r){
  
  as.table(cm) %>%
    as.tibble() %>%
    mutate(response = factor(response), 
           truth = factor(truth, rev(levels(response)))) %>%
    ggplot(aes(response, truth, fill = n)) +
    geom_tile() +
    geom_text(aes(label = n)) +
    scale_fill_gradientn(colors = rev(hcl.colors(10, "Blues")), breaks = seq(0, 400, 80)) +
    coord_fixed() +
    theme_minimal() +
    ggtitle(paste('trial', r)) +
    theme(plot.title = element_text(hjust = 0.5))
}


# RMSFE 
rmse <- function(observed, predicted) {
  rmse = sqrt(mean((observed - predicted) ^ 2))
  return(rmse)
}


# R^2
r_squared <- function(observed, predicted) {
  ss_total <- sum((observed - mean(observed)) ^ 2)
  ss_residual <- sum((observed - predicted) ^ 2)
  r_square = 1 - (ss_residual / ss_total)
  return(r_square)
}


r_squared_adj <- function(observed, predicted, p) {
  n = length(observed)
  ss_total <- sum((observed - mean(observed)) ^ 2)
  ss_residual <- sum((observed - predicted) ^ 2)
  r_square = 1 - ((ss_residual / (n - p - 1)) / (ss_total / (n - 1)))
  return(r_square)
}


# MAFE
mae <- function(observed, predicted){
  mae <- mean(abs(observed - predicted))
  return(mae)
}


# MSLFE 
msle <- function(observed, predicted){
  predicted[predicted < -1] <- -0.5
  msle <- mean((log(1 + observed) - log(1 + predicted))^2)
  return(msle)
}


# CRPS 
crps2 <- function(y_true, y_pred, sample_weight=NULL) {
  num_samples <- length(y_pred)
  absolute_error <- colMeans(abs(y_pred - y_true))

  if (num_samples == 1) {
    if (is.null(sample_weight)) {
      return(mean(absolute_error))
    } else {
      return(weighted.mean(absolute_error, sample_weight))
    }
  }

  y_pred <- sort(y_pred)
  b0 <- mean(y_pred)
  b1_values <- y_pred * c(1:num_samples)
  b1 <- mean(b1_values) / num_samples

  per_obs_crps <- absolute_error + b0 - 2 * b1

  if (is.null(sample_weight)) {
    return(abs(mean(per_obs_crps)))
  } else {
    return(abs(weighted.mean(per_obs_crps, sample_weight)))
  }
}


crps2.vec <- function(y_true, y_pred, sample_weight=NULL) {
  num_samples <- length(y_pred)
  absolute_error <- abs(y_pred - y_true)
  
  # if (num_samples == 1) {
  #   if (is.null(sample_weight)) {
  #     return(mean(absolute_error))
  #   } else {
  #     return(weighted.mean(absolute_error, sample_weight))
  #   }
  # }
  
  y_pred <- sort(y_pred)
  b0 <- mean(y_pred)
  b1_values <- y_pred * c(1:num_samples)
  b1 <- mean(b1_values) / num_samples
  
  per_obs_crps <- absolute_error + b0 - 2 * b1
  
  return(abs(per_obs_crps))
}


# Variable Selection Evaluation
vs.Eva <- function(betas_real, all.beta.est.lowbound, all.beta.est.upbound, Rep){
  
  NS = nrow(betas_real)
  Nbvar = ncol(betas_real)
  is0 = sum(betas_real == 0)
  isn0 = sum(betas_real != 0)
  is0nume = 0
  isn0nume = 0
  itnr.vec = array(0, dim = c(Rep))
  itpr.vec = array(0, dim = c(Rep))
  for (r in 1:Rep) {
    for (s in 1:NS) {
      for (b in 1:Nbvar) {
        if (betas_real[s, b] == 0){
          if (all.beta.est.lowbound[r, s, b] <= 0 && all.beta.est.upbound[r, s, b] >= 0) {
            is0nume = is0nume + 1
          }
        }
        
        if (betas_real[s, b] != 0){
          if (all.beta.est.lowbound[r, s, b] > 0 || all.beta.est.upbound[r, s, b] < 0) {
            isn0nume = isn0nume + 1
          }
        }
      }
    }
    itnr.vec[r] = is0nume / is0
    itpr.vec[r] = isn0nume / isn0
  }
  
  return(list(all.itnr = itnr.vec, all.itpr = itpr.vec))
}


# TPR FPR & MCC
tpr.fpr.mcc <- function(beta_real, all.beta.est.median, Rep) {
  
  NS = nrow(beta_real)
  Nbvar = ncol(betas_real)
  TPR.vec = c(rep(0, Rep))
  FPR.vec = c(rep(0, Rep))
  MCC.vec = c(rep(0, Rep))
  
  for (r in 1:Rep) {
    beta_est = all.beta.est.median[r, , ]
    # 初始化统计量
    TP <- sum(beta_real == 0 & beta_est == 0)  # 真阳性
    FP <- sum(beta_real == 0 & beta_est != 0)  # 假阳性
    FN <- sum(beta_real != 0 & beta_est == 0)  # 假阴性
    TN <- sum(beta_real != 0 & beta_est != 0)  # 真阴性
    
    # 计算TPR（真正率）、FPR（假阳性率）和 MCC
    TPR <- TP / (TP + FN) # Recall / Sensitivity
    FPR <- FP / (FP + TN) # Fall-out
    
    # 计算MCC（Matthews correlation coefficient）
    numerator <- (TP * TN) - (FP * FN)
    denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    
    # 防止除以0的情况
    MCC <- ifelse(denominator == 0, 0, numerator / denominator)
    
    TPR.vec[r] <- TPR
    FPR.vec[r] <- FPR
    MCC.vec[r] <- MCC
  }

  return(list(TPR = TPR.vec, FPR = FPR.vec, MCC = MCC.vec))
}


state_forecast <- function(p0_est, tran_p_test_est){
  
  N <- dim(tran_p_test_est)[1]
  NT <- dim(tran_p_test_est)[2]
  NS <- dim(tran_p_test_est)[3]
  predictProb.mat <- matrix(0, nrow = NT, ncol = NS)
  predictProb.mat[1, ] <- p0_est
  for (n in 1:N) {
    for (t in 2:NT) {
      for (s1 in 1:NS) {
        for (s2 in 1:NS) {
          predictProb.mat[t, s1] = predictProb.mat[t, s1] + predictProb.mat[t - 1, s2] * tran_p_test_est[n, t - 1, s2, s1]
        }
      }
    }
  }
  return(predictProb.mat)
}


state_forecast_hmm <- function(p0_est, tran_p, N, NT){
  
  NS <- dim(tran_p)[1]
  predictProb.mat <- matrix(0, nrow = NT_test, ncol = NS)
  predictProb.mat[1, ] <- p0_est
  for (t in 2:NT_test) {
    for (s1 in 1:NS) {
      for (s2 in 1:NS) {
        predictProb.mat[t, s1] = predictProb.mat[t, s1] + predictProb.mat[t - 1, s2] * tran_p[s2, s1]
      }
    }
  }
  return(predictProb.mat)
}


###### DMW.TEST
dmw.test.rmse <- function(y_forecast1, y_forecast2, y_true){
  # y_forecast1: always correspond to our method
  # y_forecast2: always correspond to campared method
  e1 <- (y_forecast1 - y_true)^2
  e2 <- (y_forecast2 - y_true)^2
  test.res <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 1)
  p <- test.res$p.value
  # return(p / 10)
  return(p)
}


dmw.test.r2 <- function(y_forecast1, y_forecast2, y_true){
  # y_forecast1: always correspond to our method
  # y_forecast2: always correspond to campared method
  e1 <- 1 - ((y_true - y_forecast1) / (y_true - mean(y_true)))^2
  e2 <- 1 - ((y_true - y_forecast2) / (y_true - mean(y_true)))^2
  e1[which(is.infinite(e1))] <- 1e-5
  e2[which(is.infinite(e2))] <- 1e-5
  test.res <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 1)
  p <- test.res$p.value
  return(p)
}


dmw.test.mae <- function(y_forecast1, y_forecast2, y_true){
  # y_forecast1: always correspond to our method
  # y_forecast2: always correspond to campared method
  e1 <- abs(y_true - y_forecast1)
  e2 <- abs(y_true - y_forecast2)
  test.res <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 1)
  p <- test.res$p.value
  return(p)
}


dmw.test.crps <- function(y_forecast1, y_forecast2, y_true){
  # y_forecast1: always correspond to our method
  # y_forecast2: always correspond to campared method
  e1 <- crps2.vec(y_true, y_forecast1)
  e2 <- crps2.vec(y_true, y_forecast2)
  test.res <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 1)
  p <- test.res$p.value
  return(p)
}


dmw.test.myself <- function(e1, e2){
  
  d <- (e1^2 - e2^2)
  d_bar <- mean(d)
  T <- length(d)
  s_d <- sqrt(var(d) / T)  # 标准误
  dmw_stat <- d_bar / s_d
  p_value <- 2 * (1 - pnorm(abs(dmw_stat)))
  return(p_value)
}

# Probability Score
LPS <- function(prob, realstate, NT){
  # prob: matrix, the probability of select hidden states
  # realstate: matrix, elements are 0 or 1, stand for real state
  lps <- -1 / NT * sum((1 - realstate) * log(1 - prob + 1e-10) + realstate * log(prob + 1e-10))
  return(lps)
}

QPS <- function(prob, realstate, NT, NS){
  # prob: matrix, the probability of select hidden states
  # realstate: matrix, elements are 0 or 1, stand for real state
  qps <- 2 / (NT * NS) * sum((prob - realstate)^2)
  return(qps)
}
