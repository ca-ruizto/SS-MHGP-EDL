

physical_properties_funs <- function(mu_gp, var_gp, n_funs, tension = "T"){
  #coefs gp are log strain, log stress, other coefs
  
  if(is.null(dim(mu_gp))) mu_gp <- matrix(mu_gp, nrow = 1)
  if(is.null(dim(var_gp))) var_gp <- matrix(var_gp, nrow = 1)
  
  ln_mu <- exp(mu_gp[,1:2] + var_gp[,1:2]/2)
  ln_var <- (exp(var_gp[,1:2]) - 1) * ln_mu^2
  
  
  ## Young's module
  ratio_mu <- sapply(1:(n_funs-1), function(i) mu_gp[,1] - i*mu_gp[,2])
  ratio_var <- sapply(1:(n_funs-1),
                      function(i) var_gp[,1] + i^2*var_gp[,2])
  
  lnr_mu <- exp(ratio_mu + ratio_var/2)
  lnr_var <- (exp(ratio_var) - 1) * lnr_mu ^ 2
  
  der_v <- 0:(n_funs - 1) * 0.01^(c(0,0:(n_funs - 2)))
  der_v <- t(compute_coefs_leg_bf(n_funs)) %*% der_v
  
  mult_mu <- lnr_mu * mu_gp[,-c(1:3)] 
  mult_var <- lnr_var* var_gp[,-c(1:3)] +
    lnr_var * mu_gp[,-c(1:3)]^2 +
    var_gp[,-c(1:3)] * lnr_mu^2
  ym_mean <- mult_mu%*% der_v[-1]
  ym_var <- mult_var %*% (der_v[-1] ^ 2)
  
  ## thoughness
  
  #int_b <- 1/(1:n_funs)
  #int_b <- int_b %*% compute_coefs_leg_bf(n_funs)
  # for shifted legendre poly, int pm = 0 for m > 1, so only 
  # the first coefs matters
  
  th_mean <- th_var <- rep(0, nrow(mu_gp))
  
  if(tension){
    
    th_mean <- lnr_mu[,1]  * mu_gp[,3]
    th_var <- lnr_var[,1]* var_gp[,3] +
      lnr_var[,1] * mu_gp[,3]^2 +
      var_gp[,3] * lnr_mu[,3]^2
  }
  
  out_mean <- tibble(frac_stress = ln_mu[,1],
                     perc_elong_frac = ln_mu[,2],
                     thoughness = th_mean,
                     young_mod = ym_mean)
  
  out_var <- tibble(frac_stress = ln_var[,1],
                    perc_elong_frac = ln_var[,2],
                    thoughness = th_var,
                    young_mod = ym_var)
  
  return(list(mean = out_mean,
              var = out_var))
  
}


physical_prop_uns_funs <- function(mu_gp, var_gp, n_funs, tension = "T"){
  #coefs gp are log strain, other coefs
  
  if(is.null(dim(mu_gp))) mu_gp <- matrix(mu_gp, nrow = 1)
  if(is.null(dim(var_gp))) var_gp <- matrix(var_gp, nrow = 1)
  
  
  # % Elongation
  ln_mu <- exp(mu_gp[,1] + var_gp[,1]/2)
  ln_var <- (exp(var_gp[,1]) - 1) * ln_mu^2
  
  
  ## Young's module
  ratio_mu <- sapply(1:(n_funs-1), function(i) - i*mu_gp[,1])
  ratio_var <- sapply(1:(n_funs-1),
                      function(i) i^2*var_gp[,1])
  
  lnr_mu <- exp(ratio_mu + ratio_var/2)
  lnr_var <- (exp(ratio_var) - 1) * lnr_mu ^ 2
  
  der_v <- 0:(n_funs - 1) * 0.01^(c(0,0:(n_funs - 2)))
  b_coefs <- t(compute_coefs_leg_bf (n_funs))
  #der_v <- b_coefs %*% der_v
  
  mu_gpB <- mu_gp[,-1] %*% b_coefs
  var_gpB <- var_gp[,-1]%*% (b_coefs^2)
  mult_mu <- lnr_mu * mu_gpB[,-1] 
  mult_var <- lnr_var* var_gpB[,-1] +
    lnr_var * mu_gpB[,-1]^2 +
    var_gpB[,-1] * lnr_mu^2
  ym_mean <- mult_mu%*% der_v[-1]
  ym_var <- mult_var %*% (der_v[-1] ^ 2)
  
  ## Frac strenght
  fs_v <- b_coefs %*% rep(1,n_funs)
  fs_mu <- mu_gp[,-1]%*%fs_v
  fs_var <- var_gp[,-1]%*% (fs_v^2)
  
  ## thoughness
  
  #int_b <- 1/(1:n_funs)
  #int_b <- int_b %*% compute_coefs_leg_bf(n_funs)
  # for shifted legendre poly, int pm = 0 for m > 1, so only 
  # the first coefs matters
  
  th_mean <- th_var <- rep(0, nrow(mu_gp))
  
  if(tension){
    
    th_mean <- ln_mu  * mu_gp[,2]
    th_var <- ln_var* var_gp[,2] +
      ln_var * mu_gp[,2]^2 +
      var_gp[,2] * ln_mu^2
  }
  
  out_mean <- tibble(frac_stress = as.vector(fs_mu),
                     perc_elong_frac = as.vector(ln_mu),
                     thoughness = as.vector(th_mean),
                     young_mod = as.vector(ym_mean))
  
  out_var <- tibble(frac_stress = as.vector(fs_var),
                    perc_elong_frac = as.vector(ln_var),
                    thoughness = as.vector(th_var),
                    young_mod = as.vector(ym_var))
  
  return(list(mean = out_mean,
              var = out_var))
  
}

physical_prop_sampling_funs <- function(mu_gp, var_gp, n_funs, n_samples = 1e3,
                                        cof_level= 0.05, 
                                        rseed = 1){
  #coefs gp are log strain, other coefs
  
  if(is.null(dim(mu_gp))) mu_gp <- matrix(mu_gp, nrow = 1)
  if(is.null(dim(var_gp))) var_gp <- matrix(var_gp, nrow = 1)
  
  sd_gp <- sqrt(var_gp)
  
  # % Elongation
  zq <- qnorm(1 - cof_level/2)
  el_mu <- exp(mu_gp[,1] + var_gp[,1]/2)
  el_ub <- exp(mu_gp[,1] +sd_gp[,1] *zq)
  el_lb <- exp(mu_gp[,1] - sd_gp[,1] * zq)
  
  
  ## max/fracture Strength
  b_coefs <- t(compute_coefs_leg_bf(n_funs))
  fs_v <- b_coefs %*% rep(1,n_funs)
  fs_mu <- mu_gp[,-1]%*%fs_v
  fs_sd <- sqrt(var_gp[,-1]%*% (fs_v^2))
  fs_ub <- fs_mu + fs_sd * zq
  fs_lb <- fs_mu - fs_sd * zq
  
  
  ## Young's module and thoughness
  zs <- rnorm(n_samples * ncol(mu_gp))
  zs <- matrix(zs, nrow = n_samples,
               ncol = ncol(mu_gp))
  
  n_alloys <- nrow(mu_gp)
  
  th_mu <- th_lb <- th_ub <- ym_mu <- ym_lb <- ym_ub <- rep(0, n_alloys)
  
  der_v <- 0:(n_funs - 1) * 0.01^(c(0,0:(n_funs - 2)))
  
  
  for( i in 1:n_alloys){
    X <- (zs %r*% sd_gp[i,]) %r+% mu_gp[i,] 
    
    X[,1] <- exp(X[,1])
    xig <- X[,2]*X[,1]
    #xig <- xig[xig > 0]
    
    th_mu[i] <- mean(xig)
    thq <- quantile(xig, probs= c(cof_level/2,
                                  1- cof_level/2))
    th_lb[i] <- thq[1]
    th_ub[i] <- thq[2]
    
    X[,-1] <- X[,-1] %*% b_coefs
    xig <- 0
    for(j in 2:n_funs){
      xig <- xig + (X[,j+1]/(X[,1]^(j-1)))*der_v[j]
    }
    #xig <- xig[xig > 0]
    xig <- xig[!is.na(xig)]
    ym_mu[i] <- mean(xig)
    ymq <- quantile(xig, probs= c(cof_level/2, 
                                  1- cof_level/2))
    ym_lb[i] <- ymq[1]
    ym_ub[i] <- ymq[2]
  }
  
  
  
  out_mean <- tibble(ultimate_stress = as.vector(fs_mu),
                     ultimate_strain = as.vector(el_mu),
                     toughness = as.vector(th_mu),
                     young_mod = as.vector(ym_mu))
  
  out_lb <- tibble(ultimate_stress = as.vector(fs_lb),
                   ultimate_strain = as.vector(el_lb),
                   toughness = as.vector(th_lb),
                   young_mod = as.vector(ym_lb))
  
  out_ub <- tibble(ultimate_stress = as.vector(fs_ub),
                   ultimate_strain = as.vector(el_ub),
                   toughness = as.vector(th_ub),
                   young_mod = as.vector(ym_ub))
  
  return(list(mean = out_mean,
              lb = out_lb,
              ub = out_ub))
  
}


ev_moment_fun <- function( m = 1, mu_v, var_v,log_s, probs){
  
  n_alloys <- length(mu_v)
  e_mu <- exp(m*mu_v + m^2 * var_v/2)
  sd_v <- sqrt(var_v)
  
  ep <- sapply(1:n_alloys, function(j) e_mu[j] * pnorm(q = (log_s - mu_v[j])/sd_v[j] -
                                                         m*sd_v[j],
                                                       lower.tail = T)
  )
  
  return(ep / probs)
}

predict_ss_curve_fun <- function(mu_gp, var_gp, n_points = 1e2,
                                 max_strain = 0.2){
  
  n_funs <- ncol(mu_gp) - 1
  n_alloys <- nrow(mu_gp)
  xt <- seq(0, 1, length = n_points)
  xt <- xt[-1]
  xs <- xt * max_strain
  log_xs <- log(xs)
  
  sd_gp <- sqrt(var_gp[,1])
  prob_ms <- sapply(1:n_alloys, function(i) pnorm(q = log_xs,
                                                  mean =  mu_gp[i,1],
                                                  sd = sd_gp[i],
                                                  lower.tail = F)
  )
  
  xMat <- cbind(1, stats::poly(xs, n_funs - 1, raw=TRUE))
  b_coefs <- t(compute_coefs_leg_bf (n_funs))
  
  
  mu_gpB <- mu_gp[,-1] %*% b_coefs
  var_gpB <- var_gp[,-1]%*% b_coefs^2
  
  
  out_mean <- tcrossprod(xMat[,1], mu_gpB[,1])
  out_var <- tcrossprod(xMat[,1]^2, var_gpB[,1])
  
  for(i in 1:(n_funs - 1)){
    m1 <- ev_moment_fun(m = 1,
                        mu_v = -i*mu_gp[,1],
                        var_v = i^2 * var_gp[,1],
                        log_s = -i*log_xs,
                        probs = prob_ms)
    
    m2 <- ev_moment_fun(m = 2,
                        mu_v = -i*mu_gp[,1],
                        var_v = i^2 * var_gp[,1],
                        log_s = -i*log_xs,
                        probs = prob_ms)
    
    mv <- m2 - m1^2
    
    tm <- tcrossprod(xMat[,i + 1], mu_gpB[,i + 1])
    temp_mean <- tm * m1
    tv <- tcrossprod(xMat[,i+1]^2, var_gpB[,i+1])
    temp_var <- mv*tv + mv*tm^2 + tv * m1^2
    
    out_mean <- out_mean + temp_mean
    out_var <- out_var + temp_var
  }
  
  out_mean <- out_mean * prob_ms
  out_var <- out_var * prob_ms^2
  
  return(list(mean = out_mean,
              sd = sqrt(out_var),
              strain = xs))
  
}


predict_ss_curve_sampling_fun <- function(mu_gp, var_gp,
                                          n_points = 1e2,
                                          n_log_s = 1e2,
                                          n_cond_s = 1e2,
                                          a_prob = 0.05,
                                          max_strain = 0.2){
  
  n_funs <- ncol(mu_gp) - 1
  n_alloys <- nrow(mu_gp)
  xt <- seq(0, 1, length = n_points  +1)
  xt <- xt[-1]
  xs <- xt * max_strain
  log_xs <- log(xs)
  
  sd_gp <- sqrt(var_gp[,1])
  prob_ms <- sapply(1:n_alloys, function(i) pnorm(q = log_xs,
                                                  mean =  mu_gp[i,1],
                                                  sd = sd_gp[i],
                                                  lower.tail = F)
  )
  
  out_lb <- out_ub <- out_mean <- matrix(0,nrow = n_points ,
                                         ncol = n_alloys)
  
  sd_gp <- sqrt(var_gp)
  
  n_ts <- n_log_s * n_cond_s
  b_coefs <- t(compute_coefs_leg_bf (n_funs))
  
  
  mu_gpB <- mu_gp[,-1] %*% b_coefs
  var_gpB <- var_gp[,-1]%*% b_coefs^2
  
  for(i in 1:n_alloys){
    max_s_samples <- exp(rnorm(n_log_s,
                               mean = mu_gp[i,1],
                               sd = sd_gp[i,1]))
    
    sample_stress <- matrix(NA, nrow = n_ts,
                            ncol = n_points)
    
    counter_idx <- rep(0, n_points)
    for(j in 1:n_log_s){
      
      x_temp <- xs / max_s_samples[j]
      idx_add <- which(x_temp <= 1)
      x_temp <- x_temp[idx_add]
      counter_idx[idx_add] <- counter_idx[idx_add] + 1
      
      xMat <- cbind(1, stats::poly(x_temp, n_funs - 1, raw=TRUE))
      
      sample_coefs <- MASS::mvrnorm(n = n_cond_s,
                                    mu = mu_gpB[i,],
                                    Sigma = diag(var_gpB[i,]))
      
      temp <- tcrossprod(sample_coefs, xMat)
      
      sample_stress[1:n_cond_s + (j-1)*n_cond_s, idx_add] <- temp
      
    }
    
    out_mean[,i] <- colSums(sample_stress, na.rm = T)/n_ts
    
    for(k in 1:n_points){
      sq <- quantile(sample_stress[,k], 
                     probs = c(a_prob/2, 1 - a_prob/2),
                     na.rm = TRUE)
      out_lb[k,i] <- sq[1]*counter_idx[k]/n_log_s
      out_ub[k,i] <- sq[2]*counter_idx[k]/n_log_s
      
    }
    
  }
  
  
  return(list(mean = out_mean,
              lb = out_lb,
              ub = out_ub,
              c_idx = counter_idx,
              strain = xs))
  
}


####### MV Sampling #############

rmv_simple <- function(n_samples, mu, sig){
  n_dim <- length(mu)
  n_obs <- n_samples * n_dim
  mat <- matrix(rnorm(n_obs),
                ncol = n_dim,
                nrow = n_obs)
  return((mat %r*% sig) %r+% mu)
  
}

runif_strat <- function(n_samples, n_strat){
  u <- NULL
  
  for(i in 1:n_strat) u <- c(u,
                             runif(n_samples/n_strat,
                                   min = (i-1)/n_strat,
                                   max = i/n_strat))
  
  return(u)
}

rmv_strat <- function(n_samples, mu, sig, n_strat = 10){
  n_dim <- length(mu)
  n_obs <- n_samples * n_dim
  
  breaks <- sqrt(-2 * log(1 - 1:(n_strat-1)/n_strat))
  breaks <- c(0,breaks)
  
  rads <-  sqrt(-2 * log(1 - runif_strat(n_samples = n_samples,
                                         n_strat = n_strat)))
  
  angles <- replicate(n_dim - 1 , runif_strat(n_samples = n_samples,
                                              n_strat = n_strat))
  angles <- 2*pi*angles
  
  mat <- hsc_to_cart_mat_fun(angles = angles,
                             rad = rads)
  return((mat %r*% sig) %r+% mu)
}

compute_funs_sampling <- function(n_samples, mu, sig, sampling_fun,
                                  list_funs,
                                  fun_names = NULL){
  
  # each function outputs a single number
  
  sample_coefs <- sampling_fun(n_samples = n_samples,
                               mu = mu,
                               sig = sig)
  
  n_funs <- length(list_funs)
  
  if(is.null(fun_names)) fun_names <- paste("func",
                                            1:n_funs,
                                            sep = "_")
  
  out <- vector(mode = "list",
                length = n_funs)
  for(i in 1:n_funs){
    temp_fun <- list_funs[[i]]
    
    vals <- sapply(1:n_samples,
                   function(j) temp_fun(sample_coefs[j,]) )
    
    out[[i]] <- tibble(val = vals,
                       func = fun_names[i],
                       sample = 1:n_samples)
  }
  return(rbind(out))
}

