### Initial Solution ####

uni_gp_fit_fun <- function(vs, log_mult, th0,
                           D,Dx, cor_fun,
                           nlopt_ctrl = list("algorithm" = "NLOPT_LN_BOBYQA",
                                             "maxeval" = 1e5,
                                             "xtol_abs" = 1e-8,
                                             "ftol_abs" = 1e-8)){
  
  n_prop <- ncol(vs)
  min_fun <- function(par){
    ep <- exp(par[-1])
    
    K <- cor_fun(D_mat = D, 
                 theta = th0 + par[1] )
    kx <- cor_fun(D_mat = Dx, 
                  theta = th0 + par[1])
    
    tau <- bet <- rep(0,n_prop)
    out <- 0
    
    for(i in 1:n_prop){
      A <- chol(ep[i]*K + diag(nrow(K)))
      
      Ai <- chol2inv(A)
      bet[i] <- sum(Ai%*%vs[,i])/sum(Ai)
      res <- vs[,i] - bet[i]
      tau[i] <- crossprod(res,Ai)%*% res/nrow(vs)
      out <- out + 2*sum(log(diag(A)))
    }
    
    out <- out + nrow(vs)*sum(log(tau))
    
    attr(out,"tau") <- tau
    attr(out, "bet") <- bet
    attr(out, "K") <- K
    attr(out, "kx") <- kx
    return(out)
  }
  
  opt_sol <- nloptr(x0 = c(log_mult[2],rep(0,n_prop)),
                    lb = c(log_mult[1],rep(-log(1e4),n_prop)),
                    ub = c(log_mult[3],rep(log(1e4),n_prop)),
                    eval_f = min_fun,
                    opts = nlopt_ctrl)
  
  f_eval <- min_fun(opt_sol$solution)
  bet <- attr(f_eval, "bet")
  tau <- attr(f_eval, "tau")
  K <- attr(f_eval, "K")
  kx <- attr(f_eval, "kx")
  
  Delta <- matrix(0, nrow = nrow(kx),
                  ncol = n_prop)
  es <- exp(opt_sol$solution)
  
  for(i in 1:n_prop){
    cK <- chol(es[i+1]*K + diag(nrow(K)))
    KK <-es[i+1]*kx%*%chol2inv(cK)
    
    Delta[,i] <- bet[i] + KK%*%(vs[,i] -bet[i])
  }
  
  
  return(list(Delta = Delta,
              tau = tau,
              g = es[-1],
              thm = opt_sol$solution[1]))
  
}

compute_mean_var_internal <- function(Y,a_mult,n_prop,n_unique){
  n_last <- 0
  idx_use <- which(a_mult > 1)
  y_bar <- rep(0, n_prop * n_unique)
  var_mat <- NULL
  
  j <- 0
  keep_v <- NULL
  for(i in 1:n_unique) {
    idx_p <- 1:n_prop + (i-1)*n_prop
    if(a_mult[i] == 1){
      y_bar[idx_p] <- Y[n_last + 1:a_mult[i],]
    }else{
      j <- j + 1
      id_y <- n_last + 1:a_mult[i]
      y_bar[idx_p] <- colMeans(Y[id_y,])
      v_temp <- rep(0,n_prop)
      for( k in 1:n_prop) v_temp[k] <- var(Y[id_y, k])
      
      if(min(v_temp)>0){
        keep_v <- c(keep_v,
                    j)
        if(j == 1){
          var_mat <- v_temp
        }else{
          var_mat <- rbind(var_mat,
                           v_temp)
        }
      }
    }
    
    n_last <- a_mult[i] + n_last
    
  }
  return(list(y_bar = y_bar,
              idx = idx_use[keep_v],
              var_mat = var_mat))
}

compute_mean_var_fun <- function(Y, a_mult, rescale_y = "after"){
  #re-scale = 0, means no rescaling
  #=1 to rescale before calculating means
  #=2 for after
  n_unique <- length(a_mult)
  
  n_prop <- ncol(Y)
  
  rs_mean <- rs_mult <- 1
  
  if(rescale_y == "before"){
    rs_mean <- colMeans(Y)
    rs_mult <- sapply(1:ncol(Y), function(i) sd(Y[,i]) )
    Y <- (Y %r-% rs_mean) %r/% rs_mult
  } 
  
  out <- compute_mean_var_internal(Y=Y,
                                   a_mult = a_mult,
                                   n_prop =n_prop,
                                   n_unique = n_unique)
  
  
  Yb <- matrix(out$y_bar,
               nrow = n_unique,
               byrow = T)
  if(rescale_y == "after"){
    rs_mean <- colMeans(Yb)
    rs_mult <- sapply(1:ncol(Yb), function(i) sd(Yb[,i]) )
    Y <- (Y %r-% rs_mean) %r/% rs_mult
  } 
  
  
  out <- compute_mean_var_internal(Y=Y,
                                   a_mult = a_mult,
                                   n_prop =n_prop,
                                   n_unique = n_unique)
  
  Yb <- matrix(out$y_bar,
               nrow = n_unique,
               byrow = T)
  out$var_y <- diag(cov(Yb))
  
  out$var_rep <- diag(cov(out$var_mat))
  out$rs_mean = rs_mean
  out$rs_mult = rs_mult
  out$Y <- Y
  
  return(out)
}

compute_mean_var_new_sample_fun <- function(Y, a_mult, rs_mean,rs_mult){
  #re-scale = 0, means no rescaling
  #=1 to rescale before calculating means
  #=2 for after
  n_unique <- length(a_mult)
  
  n_prop <- ncol(Y)

  Y <- (Y %r-% rs_mean) %r/% rs_mult
  
  
  out <- compute_mean_var_internal(Y=Y,
                                   a_mult = a_mult,
                                   n_prop =n_prop,
                                   n_unique = n_unique)
  
  Yb <- matrix(out$y_bar,
               nrow = n_unique,
               byrow = T)
  out$var_y <- diag(cov(Yb))
  
  out$var_rep <- diag(cov(out$var_mat))
  out$rs_mean = rs_mean
  out$rs_mult = rs_mult
  out$Y <- Y
  
  return(out)
}



initial_delta_reg_fun <- function(var_mat, X, 
                                  edge_locs= NULL, pX, 
                                  idx_use, cor_fun, theta,
                                  log_mult){
  
  #dX_mat <- mv_dist_fun(pX)
  X0 <- X[idx_use,]
  if(!is.null(edge_locs)) X0 <- rbind(X0, edge_locs)
  dComb_mat <- mv_dist_fun(X2=X0,X=pX)
  X_mat <- mv_dist_fun(X0)
  
  
  vs <- log(var_mat)
  if(!is.null(edge_locs)){
    vs <- rbind(vs,
                matrix(-25,nrow = nrow(edge_locs),
                       ncol = ncol(vs)))
  } 
  result <- uni_gp_fit_fun(vs = vs,
                           log_mult = log_mult,
                           th0 = theta,
                           D = X_mat,
                           Dx = dComb_mat,
                           cor_fun = cor_fun)
  
  return(result)
}

ini_pars_hom_fun <- function(rs = NULL, n_dim, prev_pars = NULL,
                             log_bounds = log(c(1e-2,2,50))){
  
  if(is.null(prev_pars)){
    tmp <- log(rs$var_y[1])
    g <- log(rs$var_rep) - tmp
    s <- log(rs$var_y[-1]) - tmp
  }else{
    n_prop <- length(prev_pars) - n_dim
    n_prop <- ceiling(n_prop/2)
    g <- prev_pars[1:n_prop + n_dim]
    s <- prev_pars[1:(n_prop-1) + n_dim + n_prop]
  }
  
  
  ini_sol <- c(rep(log(0.5),n_dim),
               s,g)
  
  lb <- c(rep(log_bounds[1], n_dim),
         s - log_bounds[3],
         g - log_bounds[3])
  ub <- c(rep(log_bounds[2], n_dim),
          s + log_bounds[3],
          g + log_bounds[3])
  
  return(list(ini = ini_sol,
              lb = lb,
              ub = ub))
  
}


update_ini_pars_het <- function(prev_pars, n_dim, n_prop,
                                log_bounds = log(c(pi/10,2*pi,50))){
  id_ths <- c(1:n_dim, 1:n_dim + n_dim + n_prop -1)
  het_pars <- prev_pars
  low_b <- prev_pars - log_bounds[3]
  upp_b <- prev_pars + log_bounds[3]
  low_b[id_ths] <- log_bounds[1]
  upp_b[id_ths] <- log_bounds[2]
  
  return(list(ini=het_pars,
              ub = upp_b,
              lb = low_b))
}

update_ini_pars_fun <- function(new_ini, ini_sol, log_bounds = log(2)){
  ini_sol$ini = new_ini
  
  idx <- which(new_ini >= ini_sol$ub)
  if(length(idx) > 0) ini_sol$ub[idx] <- ini_sol$ub[idx] + log_bounds
  
  idx <- which(new_ini <= ini_sol$lb)
  if(length(idx) > 0) ini_sol$lb[idx] <- ini_sol$lb[idx] - log_bounds
  
  return(ini_sol)
}


pars_hom_to_het_reg <- function(hom_pars,
                            n_properties, 
                            var_mat,X,pX,idx_use,
                            cor_fun,
                            phi_h,
                            edge_locs = NULL,
                            #n_edges = 0,
                            log_mult = log(0),
                            log_bounds = log(c(pi/10,2*pi,50))){
  
  #length scale and cov of mean process
  n_dim <- ncol(X)
  hp1 <- hom_pars[1:n_dim]
  hp2 <- hom_pars[1:(n_properties-1) + n_dim]
  low_b <- c(rep(log_bounds[1],n_dim),
             hp2 - log_bounds[3])
  upp_b <- c(rep(log_bounds[2],n_dim),
             hp2 + log_bounds[3])
  het_pars <- c(hp1, hp2)
  
  #reg solution
  #check if phi should be multiplied or divided
  vm <- var_mat/phi_h
  
  result <- initial_delta_reg_fun(var_mat = vm,
                                  edge_locs = edge_locs,
                                  X = X, pX = pX,
                                  idx_use = idx_use,
                                  cor_fun = cor_fun,
                                  theta = hp1,
                                  log_mult = log_mult)
  
  tau_t <- result$tau*result$g
  g_t <- result$tau
  hp3 <- log(c(tau_t,g_t)) #thg
  
  het_pars <- c(het_pars,
                hp1 + result$thm,
                hp3)
                
  low_b <- c(low_b,
             rep(log_bounds[1],n_dim),
             hp3 - log_bounds[3])
  
  upp_b <- c(upp_b,
             rep(log_bounds[2],n_dim),
             hp3 + log_bounds[3])
  
  #Delta
  delta_ini <- as.vector(t(result$Delta))
  low_b <- c(low_b,
          delta_ini - abs(delta_ini))
  upp_b <- c(upp_b,
             delta_ini + abs(delta_ini))
  
  
  het_pars <- c(het_pars , 
                delta_ini)
  
  return(list(ini=het_pars,
              ub = upp_b,
              lb = low_b))
}



log_chol_fun <- function(par,n_dim){
  #returns a correlation matrix under log-chol parameterization
  L <- Diagonal(x = exp(par[1:n_dim]))
  
  L[lower.tri(L)] <- par[-(1:n_dim)]
  return(tcrossprod(L))
}

log_chol_diag_spa_fun <- function(par)  Diagonal(x=exp(par))

log_chol_diag_fun <- function(par) diag(x = exp(par))

fast_diagMK_fun <- function(M, K){
  
  MK <- tcrossprod(M, K)
  
  out <- sapply(1:nrow(K), function(i) K[i,] %*% MK[,i] )
  
  return(as.vector(out))
}


####### Multivariate MVRE ###########

compute_beta_fun <- function(D,K,v){
  DK <- crossprod(D,K)
  solve(DK%*% D,
        DK %*% v)
}


log_like_mv_het_spa_gp_fun <- function(D0, X_pred,X_pred_all,
                                       y_bar, y, Delta, a_mult, 
                                   n_properties, 
                               theta, theta_cov, 
                               mult_theta_sv=NULL, theta_sv = NULL,
                               theta_cov_sv = NULL,
                               g, dX = NULL, a_mult_pX= NULL,
                               n_edge = 0,#id_train_edge = NULL,
                               cor_fun_gp, cor_fun_gp_sv = NULL, cov_fun, dComb = NULL,
                      eps = sqrt(.Machine$double.eps),
                      penalty = F){
  #sv stands for simulation variability
  
  #setups
  n_obs <- nrow(D0[[1]])
  N_obs <- length(y)/n_properties
  N_tobs <- length(y)
  
  if(is.null(theta_sv) ) theta_sv <- k_theta_sv * theta_sv
  if(is.null(cor_fun_gp_sv) ) cor_fun_gp_sv <- cor_fun_gp
  A_mult <- rep(a_mult,each = n_properties)
  
  
  # Temporarily store Cholesky transform of K in Ki
  
  Sig_sv <- cov_fun(theta_cov_sv)
  
  
  if(is.null(dX)){
    G <- rep(g, n_obs)
    
      C_sv <- cor_fun_gp_sv(D_mat = D0, 
                         theta = theta_sv)
      CS_sv <- kronecker(C_sv, Sig_sv)
      
      CS_sv <- CS_sv + Diagonal(x = eps + G/A_mult)
      
    K_sv_c <- chol(CS_sv)
    K_svi <- chol2inv(K_sv_c)
    
    Di_sv <- kronecker(rep(1,n_obs),
                       Diagonal(n_properties))
    beta_sv <- compute_beta_fun(D= Di_sv,
                                K= CS_sv,
                                v = Delta)
    
    nmean <- rep(beta_sv, n_obs)
    M <- CS_sv %*% (K_svi %*% (Delta - nmean))
  }else{
    if(n_edge > 0) Delta <- c(Delta, rep(-25, n_properties*n_edge))
    
    G <- rep(g, nrow(dX[[1]]))
    A_mult_px <- rep(a_mult_pX,each = n_properties)
    
      C_sv <- cor_fun_gp_sv(D_mat = dX, 
                            theta = theta_sv)
      CS_sv <- kronecker(C_sv, Sig_sv)
      CS_sv <- CS_sv + Diagonal(x = eps + G/A_mult_px)
    
    K_sv_c <- chol(CS_sv)
    K_svi <- chol2inv(K_sv_c)
    
    kg <- cor_fun_gp_sv(D_mat = dComb, 
                     theta = theta_sv)
    kg <- kronecker(kg,Sig_sv)
    
    Di_sv <- kronecker(rep(1,nrow(dX[[1]])),
                         Diagonal(n_properties))
    
    beta_sv <- compute_beta_fun(D= Di_sv,
                                K= K_svi,
                                v = Delta)
    
    #nmean <- drop(rowSums(K_svi) %*% Delta / sum(K_svi)) ## ordinary kriging mean
    nmean <- rep(beta_sv, nrow(dX[[1]]))
    
    M <- kg %*% (K_svi %*% (Delta - nmean))
    #nmean <- 0
    #M <- kg %*% (K_svi %*% Delta )
    
  }
  
  log_lam <- drop(rep(beta_sv, n_obs) + M)
  Lambda <- exp(log_lam)
  
  LambdaN <- NULL
  
  for(i in 1:n_obs) {
    idx <- 1:n_properties + (i-1)*n_properties
    LambdaN <- c(LambdaN,
                 rep(Lambda[idx],times = a_mult[i]))
  }
  
  
  
  # Temporarily store Cholesky transform of K in Ki
  C <- cor_fun_gp(D_mat = D0, 
                  theta = theta)
  Sig <- cov_fun(theta_cov)
  CS <- kronecker(C,Sig)
  CS_lam <- CS + Diagonal(x = Lambda/A_mult + eps)
  Ki <- chol(CS_lam)
  ldetKi <- - 2 * sum(log(diag(Ki))) # log determinant from Cholesky
  Ki <- chol2inv(Ki)
  
  
  beta0 <- compute_beta_fun(D= X_pred,
                            K= Ki,
                            v = y_bar)
  B <- X_pred %*% beta0
  ybB <- y_bar - B
  psi_0 <- drop(crossprod(ybB, Ki) %*% ybB)
  
  B2 <- X_pred_all %*% beta0
  
  psi_1 <- crossprod((y - B2)/LambdaN, y - B2) - 
    crossprod((y_bar - B) * A_mult/Lambda, y_bar - B) 
  
  psi <- (psi_0 + psi_1)/N_tobs
  
  loglik <- - N_tobs* log(psi) + ldetKi - 
            sum((A_mult - 1) * log_lam + log(A_mult)) -
            N_tobs*(log(2*pi)+1)
  
  #n_tobs <- n_obs * n_properties
  if(penalty){
    nu_hat_var <- drop(crossprod(Delta - nmean, K_svi) %*% (Delta - nmean))/length(Delta)
    
    ## To avoid 0 variance, e.g., when Delta = nmean
    if(nu_hat_var < eps){
      return(loglik)
    } 
    
    # if(hardpenalty)
    #   return(loglik + min(0, - n/2 * log(nu_hat_var) - sum(log(diag(Kg_c))) - n/2*log(2*pi) - n/2))
    G <- rep(g, nrow(D0[[1]]))
    
    C_svN <- cor_fun_gp_sv(D_mat = D0, 
                            theta = theta_sv)
    CS_svN <- kronecker(C_svN, Sig_sv)
    CS_svN <- CS_svN + diag(x = eps + G/A_mult)
    
    K_svN <- chol(CS_svN)
    
    n_tobs <- n_obs * n_properties
    pen <- - n_tobs * log(nu_hat_var) -
            2*sum(log(diag(K_svN))) -
            n_tobs*log(2*pi) - n_tobs
    
    
    return(loglik + pen)
  }
  return(as.numeric(loglik))
}

log_like_mv_hom_gp_fun <- function(D0, y_bar, X_pred = NULL, X_pred_all = NULL,
                                   y, a_mult, n_properties,
                               theta, theta_cov, 
                               g,  beta0 = NULL, 
                               cor_fun_gp,  cov_fun, 
                               eps = sqrt(.Machine$double.eps)){
  #sv stands for simulation variability
  
  #setups
  n_obs <- nrow(D0[[1]])
  N_obs <- length(y)/n_properties
  N_tobs <- length(y)
  
  A_mult <- rep(a_mult,each = n_properties)
  
  # Temporarily store Cholesky transform of K in Ki
  #print(theta)
  #print(theta_cov)
  #print(g)
  
  C <- cor_fun_gp(D_mat = D0, 
                  theta = theta)
  Sig <- cov_fun(theta_cov)
  CS <- kronecker(C,Sig)
  G <- rep(g, n_obs)
  GN <- rep(g, N_obs)
  CS_lam <- CS + Diagonal(x = G/A_mult + eps)
  Ki <- chol(CS_lam)
  ldetKi <- - 2*sum(log(diag(Ki))) # log determinant from Cholesky
  Ki <- chol2inv(Ki)

  
  if(is.null(beta0))
    beta0 <- solve(t(X_pred)%*% Ki %*% X_pred,
                   t(X_pred)%*% Ki %*% y_bar)
  B <- X_pred %*% beta0
  psi_0 <- drop(crossprod(y_bar - B, Ki) %*% (y_bar - B))
  
  B2 <- X_pred_all %*% beta0
  
  #psi <- 1/N * ((crossprod(Z - beta0) - crossprod((Z0 - beta0) * mult, Z0 - beta0))/g + psi_0)
  psi <-  psi_0 + sum((y - B2)^2/GN) - sum((y_bar - B)^2 * A_mult/G) 
  psi <- psi/N_tobs
  
  loglik <- - N_tobs* log(psi) + ldetKi -
            sum((A_mult - 1) * log(G)) - sum(log(A_mult)) -
            N_tobs*(1 + log(2*pi))
  
  attr(loglik, "psi") <- psi
  
  return(loglik)
}

fit_mvgp_sv_fun <- function(Y,X,
                            X_pred = NULL, X_pred_all = NULL,
                            a_mult,
                            pX = NULL,a_mult_px = NULL, 
                            n_edge = 0,
                            cor_fun,cor_sv_fun = NULL,
                            cov_fun = log_chol_diag_spa_fun,
                            het_sv = T,
                            ini_pars , ub = NULL,lb = NULL,
                            max_it = 1,verbose = TRUE,ftol_it = 1e-3,
                            rescale_y = "after",rs = NULL,
                            opt = T,
                            eps = sqrt(.Machine$double.eps),
                            nlopt_ctrl = list("algorithm" = "NLOPT_LN_BOBYQA", #NLOPT_LN_BOBYQA
                                              "maxeval" = 1e5,
                                              "xtol_abs" = 1e-6,
                                              "ftol_abs" = 1e-6),
                            ...){
  
  #parameters are covariance kernel parameters (ncols X) and parameters of Sigma 
  #kernel_function is the kernel used to compute the correlation between alloys
  # ind_corr = T constructs Omega as a diagonal funciton.
  # ini_pars = initial estimate, usually a vector of zeros
  # Y has multiple replications, ordered as  candidate, replication, 1 property per column
  
  ## Set up
  
  D_mat <- mv_dist_fun(X)
  dX_mat <- dComb_mat <- NULL
  if(!is.null(pX)){
    dX_mat <- mv_dist_fun(pX)
    dComb_mat <- mv_dist_fun(X=X,X2=pX)
  }
  
  n_dim <- ncol(X)
  n_X <- nrow(X)
  n_properties <- ncol(Y)
  
  if(is.null(rs)){
    rs <- compute_mean_var_fun(Y=Y,a_mult = a_mult,
                               rescale_y = rescale_y)
    Y <- rs$Y
  }else{
    Y <- (Y %r-% rs$rs_mean) %r/% rs$rs_mult
    rs_new <- compute_mean_var_internal(Y = Y,
                                    a_mult = a_mult,
                                    n_prop = n_properties,
                                    n_unique = length(a_mult))
    rs$y_bar <- rs_new$y_bar
  }
  
  
  y <- as.vector(t(Y))
  y_bar <- rs$y_bar
  
  if(is.null(X_pred)) X_pred <-  kronecker(rep(1,n_X),
                                           Diagonal(n_properties))
  
  if(is.null(X_pred_all)) X_pred_all <- kronecker(rep(1,sum(a_mult)),
                                                  Diagonal(n_properties))
  
  if(verbose){
    message("Starting")
    message(Sys.time())
  }
  
  
  ## negative log-likelihood function
  
  nll_fun <- function(pars){
    
    th <- pars[1:n_dim]
    th_cov <- c(0,pars[n_dim + 1:(n_properties -1)])
    th_cov_sv <- th_sv <- mts <- NULL
    
    if(het_sv){
      idx <- n_dim + n_properties -1
        
      th_sv <- pars[1:n_dim + idx]
      th_cov_sv <- pars[1:n_properties + 2*n_dim + n_properties -1 ]
      th_g <- pars[ 1:(n_properties)+ 2*n_dim + 2*n_properties -1]
      last_idx <- 3*n_properties + 2*n_dim -1
      th_delta <- pars[-c(1:last_idx)]
      
      ll <- log_like_mv_het_spa_gp_fun(D0 = D_mat,
                                       X_pred = X_pred,
                                       X_pred_all = X_pred_all,
                                     y_bar = y_bar,
                                     y = y,
                                     Delta = th_delta, 
                                     a_mult = a_mult,
                                     n_edge = n_edge,
                                     n_properties = n_properties,
                                     theta = th, 
                                     theta_cov = th_cov, 
                                     theta_cov_sv = th_cov_sv,
                                     mult_theta_sv=mts, theta_sv = th_sv,
                                     g = exp(th_g),
                                     dX = dX_mat,
                                     a_mult_pX = a_mult_px,
                                     cor_fun_gp = cor_fun, 
                                     cor_fun_gp_sv = cor_sv_fun,
                                     cov_fun = cov_fun, 
                                     dComb = dComb_mat,
                                     eps = eps)
      
      
    }else{
      th_g <- pars[ 1:(n_properties)+ n_dim + n_properties -1]
      
      ll <- log_like_mv_hom_gp_fun(D0 = D_mat,
                                   X_pred = X_pred,
                                   X_pred_all = X_pred_all,
                                   y_bar = y_bar,
                                   y = y,
                                   a_mult = a_mult,
                                   n_properties = n_properties,
                                   theta = th, 
                                   theta_cov = th_cov, 
                                   g = exp(th_g),
                                   cor_fun_gp = cor_fun, 
                                   cov_fun = cov_fun)
      
    }
    
    return(-ll)
  }
  
  if(!opt) return(nll_fun(ini_pars))
  
  if(length(lb) != length(ini_pars)) lb <- rep(lb, length(ini_pars))
  if(length(ub) != length(ini_pars)) ub <- rep(ub,length(ini_pars))
  
  over_ub <- which(ini_pars > ub)
  if(length(over_ub) > 0) ini_pars[over_ub] <- ub[over_ub] - eps
  under_lb <- which(ini_pars < lb)
  if(length(under_lb) > 0) ini_pars[under_lb] <- lb[under_lb] + eps
  
  opt_sol <- nloptr(x0 = ini_pars,
                    eval_f = nll_fun,
                    lb = lb,
                    ub = ub,
                    opts = nlopt_ctrl)
  
  if(max_it > 1){
    
    sol_hist <- vector(mode = "list",
                       length = max_it)
    sol_hist[[1]] <- opt_sol
    
    continue <- TRUE
    iteration <- 1

    while(continue){
      iteration <- iteration + 1
      
      if (verbose){
        message(gettextf("iteration %d", iteration))
        message(Sys.time())
      }
      
      prev_sol <- opt_sol
      
      new_ini_sol <- prev_sol$solution
      
      over_ub <- which(new_ini_sol > ub)
      if(length(over_ub) > 0) new_ini_sol[over_ub] <- ub[over_ub] - eps
      under_lb <- which(new_ini_sol < lb)
      if(length(under_lb) > 0) new_ini_sol[under_lb] <- lb[under_lb] + eps
      
      opt_sol <- nloptr(x0 = new_ini_sol,
                        eval_f = nll_fun,
                        lb = lb,
                        ub = ub,
                        opts = nlopt_ctrl)
      sol_hist[[iteration]] <- opt_sol
      
      
      
      continue <- (prev_sol$objective - opt_sol$objective) > abs(prev_sol$objective) * ftol_it
      continue <- iteration < max_it & continue
      #continue <- length(under_lb)  +length(under_lb)  == 0  & continue
      
    }
    
    return(list(opt_sol = opt_sol,
                sol_hist = sol_hist,
                Y = Y,
                y_bar = y_bar,
                rs_y = rescale_y,
                rs_mean = rs$rs_mean,
                rs_mult = rs$rs_mult))
     
  }else{
    return(list(opt_sol = opt_sol,
                Y = Y,
                y_bar = y_bar,
                rs_y = rescale_y,
                rs_mean = rs$rs_mean,
                rs_mult = rs$rs_mult))
  }
  
}

### Prediction Objects #####

construct_pred_hom_obj_fun <- function(pars,X, Y_mat, y_bar, 
                                       X_pred = NULL, X_pred_all = NULL,
                                   n_properties, a_mult,
                                   cor_fun_gp,
                                   cor_fun_gp_sv = NULL,
                                   joint_cov= T,
                                   cov_type= "Matern5_2",
                                   rescale_y=F,
                                   rs_mean = NULL,
                                   rs_mult = NULL,
                                   eps = sqrt(.Machine$double.eps)){
  
  #Preliminaries
  y <- as.vector(t(Y_mat))
  N_obs <- length(y)/n_properties
  N_tobs <- length(y)
  
  if(is.null(cor_fun_gp_sv) ) cor_fun_gp_sv <- cor_fun_gp
  
  if(is.null(X_pred)) X_pred <-  kronecker(rep(1,nrow(X)),
                                           Diagonal(n_properties))
  
  if(is.null(X_pred_all)) X_pred_all <- kronecker(rep(1,sum(a_mult)),
                                                  Diagonal(n_properties))

  D0 <- mv_dist_fun(X)
  
  n_dim <- ncol(X)
  n_obs <- nrow(X)
  Di_ones <- kronecker(rep(1,nrow(X)),
                       Diagonal(n_properties))
  
  
  #Extract parameters
  
  
  theta <- pars[1:n_dim]
  theta_cov <- c(0,pars[n_dim + 1:(n_properties-1)])
  g <- pars[ 1:n_properties + n_dim + n_properties-1]
  g <- exp(g)
  
  A_mult <- rep(a_mult,each = n_properties)
  
  
  Lambda <- rep(g, n_obs)
  
  LambdaN <- rep(g, N_obs)
  
  C <- cor_fun_gp(D_mat = D0, 
                  theta = theta)
  Sig <- log_chol_diag_fun(theta_cov)
  CS <- kronecker(C,Sig)
  CS_lam <- CS + Diagonal(x = Lambda/A_mult + eps)
  Ki <- chol2inv(chol(CS_lam))
  
  beta0 <- compute_beta_fun(D= X_pred,
                            K= Ki,
                            v = y_bar)
  
  B <- X_pred %*% beta0
  KiYB <- Ki %*% (y_bar - B)
  psi_0 <- drop(crossprod(y_bar - B, Ki) %*% (y_bar - B))
  
  B2 <- X_pred_all%*%beta0
  
  tau_hat <-  psi_0 + crossprod((y - B2)/LambdaN, y - B2) - 
    crossprod((y_bar - B) * A_mult/Lambda, y_bar - B) 
  tau_hat <- tau_hat/N_tobs
  
  out <- list(y_bar = y_bar,
              X0 = X,
              a_mult = a_mult,
              g = g,
              theta = theta,
              Sig = Sig,
              Ki = Ki,
              KiYB = KiYB,
              beta0 = beta0,
              tau = tau_hat[1],
              cor_fun_gp_sv = cor_fun_gp_sv,
              cor_fun_gp = cor_fun_gp,
              eps = eps,
              rescale_y = rescale_y,
              rs_mean = rs_mean,
              rs_mult = rs_mult,
              n_properties = n_properties,
              cov_type = cov_type)
}

predict_mv_hom_gp_fun <- function(newX,newX_pred = NULL,
                                  pred_obj, lambda_only = F){
  
  if(is.null(dim(newX))) newX <- matrix(newX, nrow=1)
  
  if(is.null(newX_pred)) newX_pred <- kronecker(rep(1,nrow(newX)),
                                                Diagonal(pred_obj$n_properties))
  
  n_locs <- nrow(newX)
  
  
  dX <- mv_dist_fun(X = newX, 
                    X2 = pred_obj$X0)
  
  kx <- pred_obj$cor_fun_gp(D_mat = dX,
                            theta = pred_obj$theta)
  kx <- kronecker(kx, pred_obj$Sig)
  
  Lambda <- rep(pred_obj$g,n_locs)
  if(lambda_only) return(Lambda)
  
  nugs <- pred_obj$tau * Lambda

  
  sd2 <- rep(diag(pred_obj$Sig),n_locs) -
    fast_diagMK_fun(pred_obj$Ki, kx)
  
  sd2 <- pred_obj$tau * sd2
  
  overall_mean <- drop( newX_pred %*% pred_obj$beta0 +
                          kx %*% pred_obj$KiYB)
  
  if(pred_obj$rescale_y != "none" ){
    om <- matrix(as.vector(overall_mean), 
                 ncol = pred_obj$n_properties,
                 byrow = T)
    om <- (om %r*% pred_obj$rs_mult) %r+% pred_obj$rs_mean
    
    overall_mean <- as.vector(t(om))
    
    sd2 <- matrix(as.vector(sd2),
                  ncol = pred_obj$n_properties,
                  byrow = T)
    sd2 <- sd2 %r*% pred_obj$rs_mult^2
    sd2 <- as.vector(t(sd2))

    
    nugs <- matrix(as.vector(nugs),
                   ncol = pred_obj$n_properties,
                   byrow = T)
    nugs <- nugs %r*% pred_obj$rs_mult^2
    nugs <- as.vector(t(nugs))
  }
  
  out <- list(mean = overall_mean,
              sd2 = sd2,
              nugs = nugs)
  return(out)
}


construct_pred_obj_fun <- function(pars,X, 
                                   X_pred = NULL, X_pred_all = NULL,
                                   Y_mat,y_bar, pX = NULL,
                                   n_properties, a_mult,a_mult_pX,
                                   cor_fun_gp,
                                   n_edge =0 ,
                                   cor_fun_gp_sv = NULL,
                                   cov_type= "Matern5_2",
                                   rescale_y=F,
                                   rs_mean = NULL,
                                   rs_mult = NULL,
                                   eps = sqrt(.Machine$double.eps)){
  
  #Preliminaries
  y <- as.vector(t(Y_mat))
  N_obs <- length(y)/n_properties
  N_tobs <- length(y)
  
  if(is.null(X_pred)) X_pred <-  kronecker(rep(1,length(a_mult)),
                                           Diagonal(n_properties))
  
  if(is.null(X_pred_all)) X_pred_all <- kronecker(rep(1,sum(a_mult)),
                                                  Diagonal(n_properties))
  
  if(is.null(cor_fun_gp_sv) ) cor_fun_gp_sv <- cor_fun_gp
  
  D0 <- mv_dist_fun(X)
  
  n_dim <- ncol(X)
  n_obs <- nrow(X)
  Di_ones <- kronecker(rep(1,nrow(X)),
                       Diagonal(n_properties))
  
  
  #Extract parameters
  theta <- pars[1:n_dim]
  theta_cov <- c(0,pars[n_dim + 1:(n_properties-1)])
  theta_cov_sv <- theta_sv <- NULL
  theta_sv <- pars[1:n_dim + n_dim + n_properties -1]
  theta_cov_sv <- pars[1:n_properties + 2*n_dim + n_properties -1]
  g <- pars[ 1:(n_properties)+ 2*n_dim + 2*n_properties -1]
  last_idx <- 3*n_properties + 2*n_dim -1
  Delta <- pars[-c(1:last_idx)]
  g <- exp(g)
  
  A_mult <- rep(a_mult,each = n_properties)
  
  
  # Temporarily store Cholesky transform of K in Ki
  
  Sig_sv <- log_chol_diag_spa_fun(theta_cov_sv)
  
  
  if(is.null(pX)){
    G <- rep(g, n_obs)
    
    C_sv <- cor_fun_gp_sv(D_mat = D0, 
                       theta = theta_sv)
    CS_sv <- kronecker(C_sv, Sig_sv)
    CS_sv <- CS_sv + Diagonal(x = eps + G/A_mult)
    
    K_sv_c <- chol(CS_sv)
    K_svi <- chol2inv(K_sv_c)
    
    Di_sv <- kronecker(rep(1,n_obs),
                       Diagonal(n_properties))
    
    beta_sv <- compute_beta_fun(D= Di_sv,
                                K= K_svi,
                                v = Delta)
    kg <- NULL
    
    nmean <- Di_sv%*%beta_sv
    KsiDn <- K_svi %*% (Delta - nmean)
    M <- CS_sv %*% KsiDn
    
  }else{
    
    if(n_edge > 0) Delta <- c(Delta, rep(-25, n_properties*n_edge))
    dX <- mv_dist_fun(pX)
    dComb <- mv_dist_fun(X=X,X2=pX)
    
    G <- rep(g, nrow(dX[[1]]))
    A_mult_px <- rep(a_mult_pX,each = n_properties)
    
    C_sv <- cor_fun_gp_sv(D_mat = dX, 
                          theta = theta_sv)
    CS_sv <- kronecker(C_sv, Sig_sv)
    CS_sv <- CS_sv + Diagonal(x = eps + G/A_mult_px)
    
    K_sv_c <- chol(CS_sv)
    K_svi <- chol2inv(K_sv_c)
    
    kg <- cor_fun_gp_sv(D_mat = dComb, 
                        theta = theta_sv)
    kg <- kronecker(kg,Sig_sv)
    
    Di_sv <- kronecker(rep(1,nrow(dX[[1]])),
                       Diagonal(n_properties))
    beta_sv <- compute_beta_fun(D= Di_sv,
                                K= K_svi,
                                v = Delta)
    
    nmean <- Di_sv%*%beta_sv
    
    KsiDn <- K_svi %*% (Delta - nmean)
    M <- kg %*% KsiDn
  }
  
  Lambda <- exp(drop(Di_ones%*%beta_sv + M))
  
  LambdaN <- NULL
  
  for(i in 1:n_obs) {
    idx <- 1:n_properties + (i-1)*n_properties
    LambdaN <- c(LambdaN,
                 rep(Lambda[idx],times = a_mult[i]))
  }
  
  C <- cor_fun_gp(D_mat = D0, 
                  theta = theta)
  Sig <- log_chol_diag_spa_fun(theta_cov)
  CS <- kronecker(C,Sig)
  CS_lam <- CS + Diagonal(x = Lambda/A_mult + eps)
  Ki <- chol2inv(chol(CS_lam))
  
  
  beta0 <- compute_beta_fun(D= X_pred,
                            K= Ki,
                            v = y_bar)
  
  B <- X_pred %*% beta0
  KiYB <- Ki %*% (y_bar - B)
  psi_0 <- drop(crossprod(y_bar - B, Ki) %*% (y_bar - B))
  
  B2 <- X_pred_all%*%beta0
  
  tau_hat <-  psi_0 + crossprod((y - B2)/LambdaN, y - B2) - 
    crossprod((y_bar - B) * A_mult/Lambda, y_bar - B) 
  tau_hat <- tau_hat/N_tobs
  
  
  tau_hat_var <- drop(crossprod(Delta - nmean, K_svi) %*% (Delta - nmean))/length(Delta)
  tau_hat_var <- max(eps, tau_hat_var)
  
  out <- list(y_bar = y_bar,
              K_svi = K_svi,
              X0 = X,
              pX = pX,
              a_mult = a_mult,
              g = g,
              theta = theta,
              theta_sv = theta_sv,
              Lambda = Lambda,
              Sig = Sig,
              Sig_sv = Sig_sv,
              beta_sv = beta_sv,
              Ki = Ki,
              KiYB = KiYB,
              beta0 = beta0,
              tau = tau_hat[1],
              tau_hat_var = tau_hat_var[1],
              KsiDn= KsiDn,
              cor_fun_gp_sv = cor_fun_gp_sv,
              cor_fun_gp = cor_fun_gp,
              eps = eps,
              rescale_y = rescale_y,
              rs_mean = rs_mean,
              rs_mult = rs_mult,
              n_properties = n_properties,
              cov_type = cov_type)
}

predict_mv_het_gp_fun <- function(newX, newX_pred = NULL,
                                  pred_obj, lambda_only = F){
  
  if(is.null(dim(newX))) newX <- matrix(newX, nrow=1)
  
  if(is.null(newX_pred)) newX_pred <- kronecker(rep(1,nrow(newX)),
                                                Diagonal(pred_obj$n_properties))
  
  n_locs <- nrow(newX)
  

  dX <- mv_dist_fun(X = newX, 
                    X2 = pred_obj$X0)
  
  kx <- pred_obj$cor_fun_gp(D_mat = dX,
                            theta = pred_obj$theta)
  kx <- kronecker(kx, pred_obj$Sig)
  
  dpX <- mv_dist_fun(X = newX, 
                    X2 = pred_obj$pX)
  kx_sv <- pred_obj$cor_fun_gp(D_mat = dpX,
                               theta = pred_obj$theta_sv)
  kx_sv <- kronecker(kx_sv, pred_obj$Sig_sv)
  
  M <- kx_sv %*% pred_obj$KsiDn
  Lambda <- exp(rep(pred_obj$beta_sv,n_locs) + M)
  if(lambda_only) return(Lambda)
  
  nugs <- pred_obj$tau * Lambda
  
  
  sd2var <- rep(diag(pred_obj$Sig_sv),n_locs) -
                fast_diagMK_fun(pred_obj$K_svi, kx_sv) #+ 
                   #(1 - tcrossprod(rowSums(pred_obj$K_svi), kx_sv))^2/sum(pred_obj$K_svi))

  sd2var <- pred_obj$tau * pred_obj$tau_hat_var* sd2var
  
  
  sd2 <- rep(diag(pred_obj$Sig),n_locs) -
              fast_diagMK_fun(pred_obj$Ki, kx)
  
  sd2 <- pred_obj$tau * sd2
  
  overall_mean <- drop( newX_pred %*% pred_obj$beta0+
                          kx %*% pred_obj$KiYB)
  
  if(pred_obj$rescale_y != "none"){
    om <- matrix(as.vector(overall_mean), 
                 ncol = pred_obj$n_properties,
                    byrow = T)
    om <- (om %r*% pred_obj$rs_mult) %r+% pred_obj$rs_mean
    
    overall_mean <- as.vector(t(om))
    
    sd2 <- matrix(as.vector(sd2),
                  ncol = pred_obj$n_properties,
                  byrow = T)
    sd2 <- sd2 %r*% pred_obj$rs_mult^2
    sd2 <- as.vector(t(sd2))
    
    sd2var <- matrix(as.vector(sd2var),
                  ncol = pred_obj$n_properties,
                  byrow = T)
    sd2var <- sd2var %r*% pred_obj$rs_mult^2
    sd2var <- as.vector(t(sd2var))
    
    nugs <- matrix(as.vector(nugs),
                  ncol = pred_obj$n_properties,
                  byrow = T)
    nugs <- nugs %r*% pred_obj$rs_mult^2
    nugs <- as.vector(t(nugs))
  }
  
  out <- list(mean = overall_mean,
              sd2 = sd2,
              nugs = nugs,
              sd2var = sd2var)
  return(out)
}


import_predict_mod_het_gp_fun <- function(newX,Xdesign,newX_pred = NULL,file_path, 
                                          mod_id = 1, step_id=1,
                                          n_funs = 4,
                                          prop_names){
  load(file_path)
  
  pred <- predict_mv_het_gp_fun(newX = newX,
                                newX_pred = newX_pred,
                                pred_obj = pol[[mod_id]])
  imspe <- IMSPE_MV(pred_obj = pol[[mod_id]])
  
  gc()
  
  gp_mu <- matrix(pred$mean, 
                  ncol = length(prop_names),
                  byrow = T)
  gp_nug <- matrix(pred$nugs,
                   ncol = length(prop_names),
                   byrow = T)
  
  gp_var <- matrix(pred$sd2 + pred$nugs,
                   ncol = length(prop_names),
                   byrow = T)
  
  
  
  mp_preds <- physical_prop_uns_funs(mu_gp = gp_mu,
                                       var_gp = gp_var,
                                       n_funs = n_funs)
  colnames(gp_mu) <- prop_names
  colnames(gp_var) <- prop_names
  colnames(gp_nug) <- prop_names
  gp_mu <- as_tibble(gp_mu)%>%
    mutate(id = row_number(),
           Step = step_id)
  gp_nug <- as_tibble(gp_nug)%>%
    mutate(id = row_number(),
           Step = step_id)
  gp_var <- as_tibble(gp_var)%>%
    mutate(id = row_number(),
           Step = step_id)
  
  mp_mean <- as_tibble(Xdesign)%>%
    bind_cols(mp_preds$mean)%>%
    mutate(Type = "Mean",
           Step = step_id,
           id_obs = row_number())
  
  mp_var <- as_tibble(Xdesign)%>%
    bind_cols(mp_preds$var)%>%
    mutate(Type = "Var",
           Step = step_id,
           id_obs = row_number())
  
  return(list(gp_mu = gp_mu,
              gp_var = gp_var,
              gp_nug = gp_nug,
              mp_preds = mp_preds,
              mp_mean = mp_mean,
              mp_var  = mp_var,
              raw_pred = pred,
              imspe = imspe))
}


import_predict_mod_gp_sampling_fun <- function(newX,Xdesign,
                                                   newX_pred = NULL,
                                                   file_path = NULL,
                                                   pol = NULL,
                                               gp_type = "Het",
                                          mod_id = 1, step_id=1,
                                          n_funs = 4,
                                          n_samples = 1e3, cof_level= 0.05,
                                          rseed = 1,
                                          prop_names){
  if(!is.null(file_path)) load(file_path)
  if(is.null(nrow(newX))) newX <- matrix(newX,nrow =1)
  
  if(gp_type == "Het"){
    pred <- predict_mv_het_gp_fun(newX = newX,
                                  newX_pred = newX_pred,
                                  pred_obj = pol[[mod_id]])
  }else{
    pred <- predict_mv_hom_gp_fun(newX = newX,
                                  newX_pred = newX_pred,
                                  pred_obj = pol[[mod_id]])
  }
  
  imspe <- IMSPE_MV(pred_obj = pol[[mod_id]])
  tau <- pol[[mod_id]]$tau
  sum_S <- sum(diag(pol[[mod_id]]$Sig))
  
  
  gc()
  
  gp_mu <- matrix(pred$mean, 
                  ncol = length(prop_names),
                  byrow = T)
  gp_nug <- matrix(pred$nugs,
                   ncol = length(prop_names),
                   byrow = T)
  
  gp_var <- matrix(pred$sd2 + pred$nugs,
                   ncol = length(prop_names),
                   byrow = T)
  
  
  
  mp_preds <- physical_prop_sampling_funs(mu_gp = gp_mu,
                                     var_gp = gp_var,
                                     n_funs = n_funs, 
                                     n_samples = n_samples,
                                     cof_level= cof_level, 
                                     rseed = rseed)
  colnames(gp_mu) <- prop_names
  colnames(gp_var) <- prop_names
  colnames(gp_nug) <- prop_names
  
  test_type <- mod_id
  
  gp_mu <- as_tibble(gp_mu)%>%
    mutate(id = row_number(),
           Step = step_id,
           Test = test_type)
  gp_nug <- as_tibble(gp_nug)%>%
    mutate(id = row_number(),
           Step = step_id,
           Test = test_type)
  gp_var <- as_tibble(gp_var)%>%
    mutate(id = row_number(),
           Step = step_id,
           Test = test_type)
  
  mp_mean <- as_tibble(Xdesign)%>%
    bind_cols(mp_preds$mean)%>%
    mutate(Type = "Mean",
           Step = step_id,
           id_obs = row_number(),
           Test = test_type)
  
  mp_lb <- as_tibble(Xdesign)%>%
    bind_cols(mp_preds$lb)%>%
    mutate(Type = "LB",
           Step = step_id,
           id_obs = row_number(),
           Test = test_type)
  
  mp_ub <- as_tibble(Xdesign)%>%
    bind_cols(mp_preds$ub)%>%
    mutate(Type = "UB",
           Step = step_id,
           id_obs = row_number(),
           Test = test_type)
  
  return(list(gp_mu = gp_mu,
              gp_var = gp_var,
              gp_nug = gp_nug,
              mp_preds = mp_preds,
              mp_mean = mp_mean,
              mp_lb  = mp_lb,
              mp_ub = mp_ub,
              raw_pred = pred,
              imspe = imspe,
              tau = tau,
              theta = pol[[mod_id]]$theta,
              sumS = sum_S))
}



predict_mod_gp_sampling_fun <- function(Xdesign,
                                            gp_mu, gp_var,
                                            n_funs = 4,
                                            test_type = "Tension", step_id = 1,
                                            n_samples = 1e3, cof_level= 0.05,
                                            rseed = 1,prop_names){
  
  mp_preds <- physical_prop_sampling_funs(mu_gp = gp_mu,
                                          var_gp = gp_var,
                                          n_funs = n_funs, 
                                          n_samples = n_samples,
                                          cof_level= cof_level, 
                                          rseed = rseed)
  colnames(gp_mu) <- prop_names
  colnames(gp_var) <- prop_names
  
  
  gp_mu <- as_tibble(gp_mu)%>%
    mutate(id = row_number(),
           Step = step_id,
           Test = test_type)
  gp_var <- as_tibble(gp_var)%>%
    mutate(id = row_number(),
           Step = step_id,
           Test = test_type)
  
  mp_mean <- as_tibble(Xdesign)%>%
    bind_cols(mp_preds$mean)%>%
    mutate(Type = "Mean",
           Step = step_id,
           id_obs = row_number(),
           Test = test_type)
  
  mp_lb <- as_tibble(Xdesign)%>%
    bind_cols(mp_preds$lb)%>%
    mutate(Type = "LB",
           Step = step_id,
           id_obs = row_number(),
           Test = test_type)
  
  mp_ub <- as_tibble(Xdesign)%>%
    bind_cols(mp_preds$ub)%>%
    mutate(Type = "UB",
           Step = step_id,
           id_obs = row_number(),
           Test = test_type)
  
  return(list(gp_mu = gp_mu,
              gp_var = gp_var,
              mp_preds = mp_preds,
              mp_mean = mp_mean,
              mp_lb  = mp_lb,
              mp_ub = mp_ub))
}






