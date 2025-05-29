source('./code/simulation/tool.R')

BMIDAS_AGLasso <- function(x, y, X_fore, It, Q){
  # Main Function for BMIDAS-AGL
  # solve function can be replaced by cholesky decomposition
  
  T <- nrow(x)
  T_oos <- nrow(X_fore)
  p <- ncol(x)
  ng <- p / nrow(Q)  # number of groups == Nbvar. nrow(Q) is the polydegree + 1(contain 0).
  nj <- rep(nrow(Q), ng)  # number of predictors in each group
  groups <- vector("list", ng)
  groups[[1]] <- seq_len(nj[1])
  for (i in 2:ng) {
    groups[[i]] <- seq(sum(nj[1:(i-1)]) + 1, sum(nj[1:(i-1)]) + nj[i])
  }
  
  result_x <- standardize_me(x)
  result_y <- center_me(y)
  
  x = result_x$y
  mu_x = result_x$mu
  sig_x = result_x$sig
  y = result_y$y
  mu_y = result_y$mu
  # sig_y = result_y$sig
  # X_fore_star <- (X_fore - mu_x) / sig_x
  X_fore_star <- (X_fore - matrix(rep(mu_x, each = nrow(X_fore)), nrow = nrow(X_fore), byrow = FALSE)) / 
    matrix(rep(sig_x, each = nrow(X_fore)), nrow = nrow(X_fore), byrow = FALSE)
  
  xg <- vector("list", ng) # xg stands for expalanable variables with the index in groups 
  sg <- numeric(ng) # sg stands for the dimension of each parameter theta
  for (ig in 1:ng) {
    xg[[ig]] <- x[ ,groups[[ig]]]
    sg[ig] <- length(groups[[ig]])
  }
  
  # Gibbs Sampling
  nsave <- It$nsave
  nburn <- It$nburn
  nthin <- It$nthin # The parameter thin allows the user to specify if and how much the MCMC chains should be thinned out before storing them.
  ntot <- nsave + nburn
  ndraw <- nsave / nthin
  
  beta_draws <- matrix(0, nrow = ndraw, ncol = p + 1)
  tau2_draws <- matrix(0, nrow = ndraw, ncol = ng)
  lambda2_draws <- matrix(0, nrow = ndraw, ncol = ng)
  sig2_draws <- numeric(ndraw)
  y_fore_draws <- matrix(0, nrow = ndraw, ncol = T_oos)
  
  # Initialization
  XgtXg <- vector("list", ng)
  for (ig in 1:ng) {
    XgtXg[[ig]] <- t(x[, groups[[ig]]]) %*% x[, groups[[ig]]]
  }
  
  # Priors
  # Parameters initialization
  lambda2 <- rep(1, ng) # lambda^2
  sig2 <- 1 # sigma^2
  tau2 <- numeric(ng) # tau^2
  beta <- vector("list", ng)
  betan <- numeric(ng) # which is used to store the norm of beta_j in each group
  
  for (ig in 1:ng) { # ng: number of groups
    tau2[ig] <- (1 / rgamma(1, shape = (sg[ig] + 1) / 2, rate = lambda2[ig] / 2)) + 1e-10
    cov <- diag(rep(sig2 * tau2[ig], sg[ig]))
    # initialization for beta
    beta[[ig]] <- rep(0, sg[ig]) + rnorm(sg[ig]) %*% chol(cov)
  }
  
  beta <- do.call(cbind, beta)
  
  # Gibbs Iteration
  start_time = Sys.time()
  for (irep in 1:ntot) {
    # 1. Update beta and tau2
    for (ig in 1:ng) { # iteration with each group one by one
      D <- (1 / tau2[ig]) * diag(sg[ig])
      A <- XgtXg[[ig]] + D # calculate A_j in paper
      AA <- (A + t(A)) / 2
      beta[groups[[ig]]] <- numeric(sg[ig]) # helpful to calculate y - Z_{\j}theta_{\j} in paper
      betan[ig] <- 0
      b <- y - x %*% t(beta)
      xb <- t(x[, groups[[ig]]]) %*% b # calculate C_j in paper
      
      beta[groups[[ig]]] <- t(t(solve(AA) %*% xb) + rnorm(sg[ig]) %*% chol(sig2 * solve(AA)))
      betan[ig] <- norm(beta[groups[[ig]]], type = "2")
      a1 <- sqrt(lambda2[ig] * sig2) / betan[ig]
      a2 <- lambda2[ig]
      tau2[ig] <- (1 / rinvgauss(1, mean = a1, shape = a2)) + 1e-10
    }
    
    # 2. Update sigma2 from Inverse Gamma
    shape <- (T - 1) / 2 + p / 2 + 1e-3
    rate <- 0.5 * (crossprod(y - x %*% t(beta)) + sum((betan^2) / tau2)) + 1e-3
    sig2 <- (1 / rgamma(1, shape = shape, rate = rate)) + 1e-10
    
    # 3. Update lambda2_j (utilize Stabilization algorithm(MCEM algorithm)). zeta, kappa, nu are used in this algorithm.
    for (ig in 1:ng) {
      lambda2[ig] = rgamma(1, shape = (length(groups[[ig]])) + 1 / 2 + 1e-3, rate = tau2[ig] / 2 + 1e-3)
    }
    # print(lambda2)
    # print(eps_zeta)
    
    if ((irep > nburn) && (irep - nburn) %% nthin == 0) {
      intercept <- mu_y + sqrt(sig2 / T) * rnorm(1)
      beta_draws[(irep - nburn) / nthin, ] <- c(intercept, beta)
      tau2_draws[(irep - nburn) / nthin, ] <- tau2
      sig2_draws[(irep - nburn) / nthin] <- sig2
      lambda2_draws[(irep - nburn) / nthin, ] <- lambda2
      # y_fore_draws[(irep - nburn) / nthin, ] <- cbind(1, X_fore_star) %*% c(intercept, beta)
      y_fore_draws[(irep - nburn) / nthin, ] <- cbind(1, X_fore_star) %*% c(intercept, beta) +rnorm(1, mean = 0, sd = sqrt(sig2))
    }
    
    if (irep %% 10 == 0) {
      cat(paste("Running MCMC, ", 100 * (irep / ntot), "% completed...", "\r"))
    }
  }
    
    end_time = Sys.time()
    
    # Output
    out <- list()
    out$beta0 <- beta_draws[ ,1]
    out$original_beta <- beta_draws[ ,2:(p + 1)]
    temp <- lapply(1:ng, function(k) (t(t(beta_draws[, groups[[k]] + 1]) / sig_x[groups[[k]]]))[, , drop = FALSE]) # 为什么要除方�?
    beta.org <- do.call(cbind, temp)
    beta_true <- matrix(0, nrow = It$nsave, ncol = Nbvar)
    for (nb in 1:Nbvar) {
      beta_true[ ,nb] <- (beta.org[ ,groups[[nb]]] %*% Q) %*% as.vector(rep(1, ncol(Q)))
    }
    out$beta <- beta_true
    out$lambda2 <- lambda2_draws
    out$tau2 <- tau2_draws
    out$sigma2 <- sig2_draws
    out$y_fore <- y_fore_draws
    out$y_fore_median <- median(y_fore_draws)
    out$MCMC_time <- end_time - start_time
    return(out)
}


BMIDAS_AGLasso_SS <- function(x, y, X_fore, It, Q){
  # Main Function for BMIDAS-AGL-SS
  # solve function can be replaced by cholesky decomposition
  
  T <- nrow(x)
  T_oos <- nrow(X_fore)
  p <- ncol(x)
  ng <- p / nrow(Q)  # number of groups == Nbvar. nrow(Q) is the polydegree + 1(contain 0).
  nj <- rep(nrow(Q), ng)  # number of predictors in each group
  groups <- vector("list", ng)
  groups[[1]] <- seq_len(nj[1])
  for (i in 2:ng) {
    groups[[i]] <- seq(sum(nj[1:(i-1)]) + 1, sum(nj[1:(i-1)]) + nj[i])
  }
  
  result_x <- standardize_me(x)
  result_y <- center_me(y)
  
  x = result_x$y
  mu_x = result_x$mu
  sig_x = result_x$sig
  y = result_y$y
  mu_y = result_y$mu
  # sig_y = result_y$sig
  # X_fore_star <- (X_fore - mu_x) / sig_x
  X_fore_star <- (X_fore - matrix(rep(mu_x, each = nrow(X_fore)), nrow = nrow(X_fore), byrow = FALSE)) / 
    matrix(rep(sig_x, each = nrow(X_fore)), nrow = nrow(X_fore), byrow = FALSE)
  
  xg <- vector("list", ng) # xg stands for expalanable variables with the index in groups 
  sg <- numeric(ng) # sg stands for the dimension of each parameter theta
  for (ig in 1:ng) {
    xg[[ig]] <- x[ ,groups[[ig]]]
    sg[ig] <- length(groups[[ig]])
  }
  
  # Gibbs Sampling
  nsave <- It$nsave
  nburn <- It$nburn
  nthin <- It$nthin # The parameter thin allows the user to specify if and how much the MCMC chains should be thinned out before storing them.
  ntot <- nsave + nburn
  ndraw <- nsave / nthin
  
  beta_draws <- matrix(0, nrow = ndraw, ncol = p + 1)
  tau2_draws <- matrix(0, nrow = ndraw, ncol = ng)
  lambda2_draws <- matrix(0, nrow = ndraw, ncol = ng)
  sig2_draws <- numeric(ndraw)
  y_fore_draws <- matrix(0, nrow = ndraw, ncol = T_oos)
  pi0_draws <- numeric(ndraw)
  pi1_draws <- matrix(0, nrow = ndraw, ncol = ng)
  Z_draws <- matrix(0, nrow = ndraw, ncol = ng)
  
  # initial values for stabilization algorithm
  zeta <- rep(1, ng) # zeta is used to calculate the monotone non-increasing sequences e^{(s)}
  kappa <- 0
  nu <- 0
  
  logmax <- log(3)
  logmin <- log(1)
  diff_ <- (logmin - logmax) / (ntot - 1)
  range <- exp(seq(logmax, logmin, by = diff_))
  eps_zeta <- exp(logmax) * rep(1, ng)
  
  # Initialization
  XgtXg <- vector("list", ng)
  for (ig in 1:ng) {
    XgtXg[[ig]] <- t(x[, groups[[ig]]]) %*% x[, groups[[ig]]]
  }
  
  # Priors
  # Parameters initialization
  lambda2 <- rep(1, ng) # lambda^2
  sig2 <- 1 # sigma^2
  tau2 <- numeric(ng) # tau^2
  beta <- vector("list", ng)
  betan <- numeric(ng) # which is used to store the norm of beta_j in each group
  
  for (ig in 1:ng) { # ng: number of groups
    tau2[ig] <- (1 / rgamma(1, shape = (sg[ig] + 1) / 2, rate = lambda2[ig] / 2)) + 1e-10
    cov <- diag(rep(sig2 * tau2[ig], sg[ig]))
    # initialization for beta
    beta[[ig]] <- rep(0, sg[ig]) + rnorm(sg[ig]) %*% chol(cov)
  }
  
  beta <- do.call(cbind, beta)
  
  kappabar <- (1 + 1 / ng)
  u <- kappabar
  aa <- kappabar * (ng^u) # stands for c in paper
  bb <- 1 # stands for d in paper
  # g1 <- rgamma(1, shape = aa, rate = 1)
  # g2 <- rgamma(1, shape = bb, rate = 1)
  # pi0 <- g1 / (g1 + g2)
  pi0 <- rbeta(1, shape1 = aa, shape2 = bb, ncp = 0)
  pi1 <- numeric(ng)
  Z <- numeric(ng) # which is used to record whether theta_j == 0. If Z[ig] == 0 then theta_j == 0. Correspond to gamma_j in paper.
  
  # Gibbs Iteration
  start_time = Sys.time()
  for (irep in 1:ntot) {
    # 1. Update beta and tau2
    u_rand <- runif(ng, min = 0, max = 1)
    for (ig in 1:ng) { # iteration with each group one by one
      D <- (1 / tau2[ig]) * diag(sg[ig])
      A <- XgtXg[[ig]] + D # calculate A_j in paper
      AA <- (A + t(A)) / 2
      beta[groups[[ig]]] <- numeric(sg[ig]) # helpful to calculate y - Z_{\j}theta_{\j} in paper
      Z[ig] <- 0
      betan[ig] <- 0
      b <- y - x %*% t(beta)
      xb <- t(x[, groups[[ig]]]) %*% b # calculate C_j in paper
      
      # calculate pi1
      maxAA <- max(AA)
      L <- (-sg[ig] / 2) * log(tau2[ig]) + (-0.5) * log(det(AA / maxAA)) + 
        (-dim(AA)[1] / 2) * log(maxAA) + (0.5 / sig2) * (t(xb) %*% solve(AA) %*% xb)
      # L <- (-sg[ig] / 2) * log(tau2[ig]) + (-0.5) * log(det(AA)) + (0.5 / sig2) * (t(xb) %*% solve(AA) %*% xb)
      pi1[ig] <- pi0 / (pi0 + (1 - pi0) * exp(L))
      
      # calculate tau
      if (u_rand[ig] >= pi1[ig]) {
        beta[groups[[ig]]] <- t(t(solve(AA) %*% xb) + rnorm(sg[ig]) %*% chol(sig2 * solve(AA)))
        Z[ig] <- 1
        betan[ig] <- norm(beta[groups[[ig]]], type = "2")
        a1 <- sqrt(lambda2[ig] * sig2) / betan[ig]
        a2 <- lambda2[ig]
        tau2[ig] <- (1 / rinvgauss(1, mean = a1, shape = a2)) + 1e-10
        # tau2[ig] <- (1 / randig(a1, a2)) + 1e-10
      }
    }
    tau2[Z == 0] <- (1 / rgamma(ng - sum(Z), shape = (sg[Z == 0] + 1) / 2, rate = lambda2[Z == 0] / 2)) + 1e-10
    
    # 2. Update sigma2 from Inverse Gamma
    shape <- (T - 1) / 2 + sum(Z * sg) / 2 + 1e-3
    rate <- 0.5 * (crossprod(y - x %*% t(beta)) + sum((betan^2) / tau2)) + 1e-3
    sig2 <- (1 / rgamma(1, shape = shape, rate = rate)) + 1e-10
    # sig2 <- (1 / rgamma(1, shape = shape, scale = rate)) + 1e-10
    
    # 3. Update pi from Beta
    # g1 <- rgamma(1, shape = aa + ng - sum(Z), rate = 1)
    # g2 <- rgamma(1, shape = bb + sum(Z), rate = 1)
    # pi0 <- g1 / (g1 + g2)
    pi0 <- rbeta(1, shape1 = aa + ng - sum(Z), shape2 = bb + sum(Z))
    
    # 4. Update lambda2_j (utilize Stabilization algorithm(MCEM algorithm)). zeta, kappa, nu are used in this algorithm.
    a_n <- 1 / (zeta^0.8) # step-size for gradient update
    s_i_1 <- log(sqrt(lambda2)) # transformation of lambda to omega in paper
    s_i <- s_i_1 + a_n * ((sg + 1) - exp(2 * s_i_1) * tau2) # update of omega
    lambda2 <- exp(2 * s_i) # transformation of omega to lambda in paper
    
    # Stabilization algorithm (Core, using gradient algorithm, Algorithm 1 in appendix A) 
    # c(pmax(-5, -kappa - 1), kappa + 1) is the K^(s) in appendix A which is a monotone increasing sequence
    Kappa <- matrix(c(pmax(-5, -kappa - 1), kappa + 1), nrow = ng, ncol = 2, byrow = TRUE) # set c=5
    Delta <- abs(s_i - s_i_1) # change between two iterations
    
    if (mean(s_i >= Kappa[, 1]) == 1 && mean(s_i <= Kappa[, 2]) == 1 && mean(Delta <= eps_zeta) == 1) { # if judge condition in Algorithm 1
      zeta[Z == 1] <- zeta[Z == 1] + 1
      nu <- nu + 1
      eps_zeta[Z == 1] <- exp(log(eps_zeta[Z == 1]) + diff_)
    } else {
      ii <- 1
      while ((mean(s_i < Kappa[, 1]) > 0 || mean(s_i > Kappa[, 2]) > 0) || mean(Delta > eps_zeta) > 0) {
        l1 <- pmin(s_i_1[s_i >= Kappa[, 2]], Kappa[s_i >= Kappa[, 2], 2])
        u1 <- pmax(s_i_1[s_i >= Kappa[, 2]], Kappa[s_i >= Kappa[, 2], 2]) # Correspond to \bar{omega} >= Kappa_{u}^{(s-1)}
        s_i[s_i >= Kappa[, 2]] <- l1 + (u1 - l1) * runif(sum(s_i >= Kappa[, 2]))
        
        l2 <- pmin(s_i_1[s_i <= Kappa[, 1]], Kappa[s_i <= Kappa[, 1], 1])
        u2 <- pmax(s_i_1[s_i <= Kappa[, 1]], Kappa[s_i <= Kappa[, 1], 1]) # Correspond to \bar{omega} < Kappa_{l}^{(s-1)}
        s_i[s_i <= Kappa[, 1]] <- l2 + (u2 - l2) * runif(sum(s_i <= Kappa[, 1]))
        
        lambda2 <- exp(2 * s_i)
        Delta <- abs(s_i - s_i_1)
        tau2 <- (1 / rgamma(ng, shape = (sg + 1) / 2, rate = lambda2 / 2)) + 1e-10
        ii <- ii + 1
        if (ii >= 100) {
          break
        }
      }
      zeta[Z == 1] <- zeta[Z == 1] + 1
      kappa <- kappa + 1
      nu <- 0
      eps_zeta[Z == 1] <- exp(log(eps_zeta[Z == 1]) + diff_)
    }
    # print(lambda2)
    # print(eps_zeta)
    
    if ((irep > nburn) && (irep - nburn) %% nthin == 0) {
      intercept <- mu_y + sqrt(sig2 / T) * rnorm(1)
      # intercept <- mean(y - x %*% t(beta)) + sqrt(sig2 / T) * rnorm(1)
      beta_draws[(irep - nburn) / nthin, ] <- c(intercept, beta)
      tau2_draws[(irep - nburn) / nthin, ] <- tau2
      sig2_draws[(irep - nburn) / nthin] <- sig2
      lambda2_draws[(irep - nburn) / nthin, ] <- lambda2
      pi0_draws[(irep - nburn) / nthin] <- pi0
      Z_draws[(irep - nburn) / nthin, ] <- Z
      pi1_draws[(irep - nburn) / nthin, ] <- pi1
      y_fore_draws[(irep - nburn) / nthin, ] <- cbind(1, X_fore_star) %*% c(intercept, beta) + sqrt(sig2)  * rnorm(1)
      # y_fore_draws[(irep - nburn) / nthin, ] <- cbind(1, X_fore_star) %*% c(intercept, beta)
    }
    
    if (irep %% 10 == 0) {
      cat(paste("Running MCMC, ", 100 * (irep / ntot), "% completed...", "\r"))
    }
  }
  
  end_time = Sys.time()
  
  # Output
  out <- list()
  out$beta0 <- beta_draws[ ,1]
  out$original_beta <- beta_draws[ ,2:(p + 1)]
  out$pi0 <- pi0_draws
  out$pi1 <- pi1_draws
  temp <- lapply(1:ng, function(k) (t(t(beta_draws[, groups[[k]] + 1]) / sig_x[groups[[k]]]))[, , drop = FALSE]) # 为什么要除方�?
  beta.org <- do.call(cbind, temp)
  beta_true <- matrix(0, nrow = It$nsave, ncol = Nbvar)
  for (nb in 1:Nbvar) {
    beta_true[ ,nb] <- (beta.org[ ,groups[[nb]]] %*% Q) %*% as.vector(rep(1, ncol(Q)))
  }
  out$beta <- beta_true
  out$lambda2 <- lambda2_draws
  out$tau2 <- tau2_draws
  out$sigma2 <- sig2_draws
  out$Z <- apply(Z_draws, 2, median)
  out$y_fore <- y_fore_draws
  out$y_fore_median <- median(y_fore_draws)
  out$MCMC_time <- end_time - start_time
  return(out)
}