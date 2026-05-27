### Multivariate IMSPE #####

Wij_int <- function(X, X2 = NULL,
                    theta, cov_fun){
  
  if(is.null(X2)){
    n_row <- length(X)
    out <- matrix(0, nrow = n_row, ncol = n_row)
    for(i in 1:n_row){
      for(j in i:n_row){
        int_f <- function(x) cov_fun( abs(x-X[i])/theta ) * cov_fun(abs(x-X[j])/theta)
        int_val <- cubature::hcubature(int_f, lowerLimit = 0, upperLimit = 1)
        out[i,j] <- out[j,i] <- int_val$integral
      }
    }
  }else{
    n_row <- length(X)
    n_col <- length(X2)
    out <- matrix(0, nrow = n_row, ncol = n_col)
    for(i in 1:n_row){
      for(j in 1:n_col){
        int_f <- function(x) cov_fun( abs(x-X[i])/theta ) * cov_fun(abs(x-X2[j])/theta)
        int_val <- cubature::hcubature(int_f, lowerLimit = 0, upperLimit = 1)
        out[i,j] <- int_val$integral
      }
    }
  }
  
  return(out)
}

Wij_int_mv <- function(X, X2 = NULL, theta, cov_fun){
  n_dim <- ncol(X)
  W <- 1
  if(is.null(X2)){
    for(i in 1:n_dim) W <- W * Wij_int(X[,i], theta = theta[i], cov_fun = cov_fun)
  }else{
    for(i in 1:n_dim) W <- W * Wij_int(X[,i], X2[,j],
                                       theta = theta[i], cov_fun = cov_fun)
  }
  
  return(W)
}

Wij_mv <- function(X, X2 = NULL, theta, 
                   scale_factor  = 1, 
                   cov_type="Matern5_2"){
  #scale_factor to transform x to [0,1] range
  n_dim <- ncol(X)
  th <- exp(theta)/scale_factor
  
  X_in <- X/scale_factor
  X2_in <- X_in
  if(!is.null(X2)) X2_in <- X2 / scale_factor
  
  wij <- 1
  
  if(cov_type == "Gaussian") th <- th ^ 2
  
  for(i in 1:n_dim){
    wij <- wij * hetGP::Wij(mu1 = matrix(X_in[,i] , ncol = 1),
                            mu2 = matrix(X2_in[,i] , ncol =1),
                            theta = th[i], 
                            type = cov_type)  
  }
  
  return(wij)
}

Wij_mv2 <- function(X, X2 = NULL, theta, 
                    scale_factor  = 1, 
                    cov_type="Matern5_2"){
  #scale_factor to transform x to [0,1] range
  n_dim <- ncol(X)
  th <- exp(theta)/scale_factor 
  
  X_in <- X/scale_factor
  X2_in <- NULL
  if(!is.null(X2)) X2_in <- X2 / scale_factor
  
  
  hetGP::Wij(mu1 = X_in,
             mu2 = X2_in,
             theta = th, 
             type = cov_type) 
}

trace_inv_mv <- function(Ki, W, n_properties){
  kw <- diag(Ki %*% W)
  M <- kronecker(rep(1, length(kw)/n_properties),
                 diag(n_properties))
  
  as.vector(kw%*%M)
  
}


IMSPE_MV <- function(pred_obj, scale_factor = 1){
  
  Wijs <- Wij_mv2(X = pred_obj$X0,
                 theta = pred_obj$theta,
                 scale_factor = scale_factor,
                 cov_type = pred_obj$cov_type)
  
  Wijs_mv <- kronecker(Wijs,pred_obj$Sig^2)
  
  sfn <- scale_factor^pred_obj$n_properties
  
  result <- diag(pred_obj$Sig) - trace_inv_mv(Ki = pred_obj$Ki,
                                              W=Wijs_mv,
                                              n_properties = pred_obj$n_properties)
  
  pred_obj$tau * sum(result) *sfn
}

IMSPE_MV2 <- function(pred_obj, scale_factor = 1){
  
  Wijs <- Wij_mv(X = pred_obj$X0,
                 theta = pred_obj$theta,
                 scale_factor = scale_factor,
                 cov_type = pred_obj$cov_type)
  
  Wijs_mv <- kronecker(Wijs,pred_obj$Sig %*% pred_obj$Sig)
  sfn <- scale_factor^pred_obj$n_properties
  #Km <- matrixcalc::commutation.matrix(ncol(pred_obj$Sig), ncol(Wijs))
  #Km <- Matrix(Km, sparse = T)
  
  # result <- diag(pred_obj$Sig) - trace_inv_mv(Ki = pred_obj$Ki %*% Km,
  #                                             W=Wijs_mv,
  #                                             n_properties = pred_obj$n_properties)
  
  result <- sum(diag(pred_obj$Sig)) - sum(diag(pred_obj$Ki %*% Wijs_mv))
  
  pred_obj$tau * result *sfn
}


IMSPE_MV_NI <- function(pred_obj, max_eval = 1e5, return_object = F ){
  n_dim <- pred_obj$theta |> length()
  
  sd2_fun <- function(x){
    predict_mv_het_gp_fun(newX = x,
                                 pred_obj = pred_obj,
                                 sd2_only = T) |>
    sum()

  }  
  
  int_result <- cubature::hcubature(f = sd2_fun,
                                    lowerLimit = rep(0, n_dim),
                                    upperLimit = rep(1, n_dim),
                                    vectorInterface = FALSE,
                                    maxEval = max_eval,
                                    absError = 1e-2)
  if(return_object) return(int_result)
  
  return(int_result$integral )
}

IMSPE_MV_Adjust <- function(pred_obj, accuracy = 0.01){
  imspe <- IMSPE_MV(pred_obj = pred_obj) # analytical
  
  a <- IMSPE_MV_NI(pred_obj = pred_obj, #rough integral
                           max_eval = 1e2,
                   return_object = T)
  
  if(abs(imspe - a$integral) < a$error) return(imspe)
  
  b <- IMSPE_MV_NI(pred_obj = pred_obj, #improved estimate, assuming log reduction in error
                   max_eval = 1e3,
                   return_object = T)
  
  if(b$error < accuracy) return(b$integral)
  
  aa <- (log(b$error) - log(a$error) ) / (log10(b$functionEvaluations) - log10(a$functionEvaluations))
  nn <- (log(accuracy) - log(a$error)) / aa + log10(a$functionEvaluations)
  
  b <- IMSPE_MV_NI(pred_obj = pred_obj, #final estimate
                   max_eval = min(10^nn, 10^5),
                   return_object = F)
  
  return(b)
}

cIMSPE_MV <- function(x, pred_obj, 
                      id = NULL, reps =1,
                      scale_factor = 1,
                      Wijs = NULL,
                      cv = T){
  #id indicates design to replicate
  
  ## Test to add an extra summation
  
  ## Precalculations
  if(is.null(Wijs)) Wijs <- Wij_mv(X = pred_obj$X0,
                                   theta = pred_obj$theta,
                                   scale_factor = scale_factor,
                                   cov_type = pred_obj$cov_type)
  
  ## Evaluate for each property, multiply imspe / mean
  
  n_properties <- pred_obj$n_properties
  
  Wijs_mv <- kronecker(Wijs,pred_obj$Sig %*% pred_obj$Sig)
  Ki <- pred_obj$Ki 
  sfn <- scale_factor^pred_obj$n_properties
  if(!is.null(id)){
    idx <- 1:n_properties + (id-1)*n_properties
    tmp <- pred_obj$a_mult[id]*(pred_obj$a_mult[id] + reps)/reps
    tmp <- tmp/as.vector(pred_obj$Lambda[idx])
    tmp <- tmp - diag(Ki[idx, idx])
    result <- diag(pred_obj$Sig) - 
      trace_inv_mv(Ki = Ki,W=Wijs_mv,
                   n_properties = n_properties) -
      diag(crossprod(Ki[,idx],Wijs_mv %*% Ki[,idx]))/tmp
    
    if(cv) return(pred_obj$tau * sum(result/as.vector(pred_obj$beta0^2)) )
    
    return(pred_obj$tau * sum(result)*sfn)
  }
  
  if(is.null(dim(x))) x <- matrix(x, nrow=1)
  
  newWijs <- Wij_mv(X = x,
                    X2 = pred_obj$X0, 
                    theta = pred_obj$theta,
                    scale_factor = scale_factor,
                    cov_type = pred_obj$cov_type)
  newWijs <- kronecker(newWijs,
                       pred_obj$Sig %*% pred_obj$Sig)
  
  W11 <- Wij_mv(X2 = x, X = x, 
                theta = pred_obj$theta, 
                scale_factor = scale_factor,
                cov_type = pred_obj$cov_type)
  W11 <- kronecker(W11, pred_obj$Sig %*% pred_obj$Sig)
  
  D <- mv_dist_fun(X=x,X2=pred_obj$X0)
  
  kn1 <- pred_obj$cor_fun_gp(D,
                             theta = pred_obj$theta)
  kn1 <- kronecker(kn1, pred_obj$Sig)
  
  new_lambda <- predict_mv_het_gp_fun(pred_obj = pred_obj, 
                                      newX = x, 
                                      lambda_only = TRUE)
  n_d <- length(new_lambda)
  new_lambda <- diag(as.vector(new_lambda))
  
  vn <- drop(pred_obj$Sig - kn1 %*% tcrossprod(Ki, kn1)) +
    new_lambda + pred_obj$eps * diag(n_d)
  gn <- - tcrossprod(Ki, kn1)%*%solve(vn)
  
  result <- diag(pred_obj$Sig) - (trace_inv_mv(Ki = Ki,W=Wijs_mv,
                                               n_properties = n_properties) + 
                                    diag(crossprod(gn, Wijs_mv %*% gn)%*%vn +
                                           2 * newWijs %*% gn + W11 %*% solve(vn)))
  
  if(cv) return(sfn * pred_obj$tau * sum(result/as.vector(pred_obj$beta0^2)) )
  
  return(pred_obj$tau * sum(result)*sfn )
}


cIMSPE_MV_list <- function(x, pred_obj, 
                           id = NULL, reps =1,
                           scale_factor = 1,
                           Wijs = NULL,
                           cv = T){
  #id indicates design to replicate
  
  out <- 0
  
  for( i in 1:length(pred_obj)){
    t_tau <- pred_obj[[i]]$tau
    
    out <- out + cIMSPE_MV(x = x,
                           pred_obj = pred_obj[[i]],
                           reps = reps,
                           scale_factor = scale_factor,
                           Wijs = Wijs[[i]],
                           id = id,
                           cv = cv)
  }
  return(out)
  
}

cIMSPE_MV_new_X_batch <- function(x, pred_obj,
                                  reps = 2,
                                  scale_factor = 1,
                                  Wijs = NULL,
                                  cv = T){
  
  if(is.null(dim(x))) x<- matrix(x, nrow = 1)
  
  pred_obj_new <- pred_obj
  
  #update Wijs
  if(is.null(Wijs)) Wijs <- Wij_mv(X = pred_obj$X0,
                                   theta = pred_obj$theta,
                                   scale_factor = scale_factor,
                                   cov_type = pred_obj$cov_type)
  newWijs <- Wij_mv(X = x,
                    X2 = pred_obj$X0, 
                    theta = pred_obj$theta,
                    scale_factor = scale_factor,
                    cov_type = pred_obj$cov_type)
  
  w11 <- Wij_mv(X2 = x, X = x, 
                theta = pred_obj$theta, 
                scale_factor = scale_factor,
                cov_type = pred_obj$cov_type)
  
  n_w <- nrow(Wijs)
  Wijs_new <- diag( n_w+ 1)
  Wijs_new[1:n_w, 1:n_w] <- Wijs
  Wijs_new[n_w+1,] <- c(as.vector(newWijs), w11[1])
  Wijs_new[,n_w + 1] <- c(as.vector(newWijs), w11[1])
  
  
  #update X0
  
  pred_obj_new$X0 <- rbind(pred_obj$X0,x)
  pred_obj_new$a_mult <- c(pred_obj$a_mult, 1)
  
  new_lambda <- predict_mv_het_gp_fun(pred_obj = pred_obj, 
                                      newX = x, 
                                      lambda_only = TRUE)
  pred_obj_new$Lambda <- c(pred_obj$Lambda , 
                           as.vector(new_lambda))
  
  A_mult <- rep(pred_obj_new$a_mult,
                each = pred_obj$n_properties)
  D0 <- mv_dist_fun(X =pred_obj_new$X0)
  C <- pred_obj_new$cor_fun_gp(D_mat = D0, 
                               theta = pred_obj_new$theta)
  CS <- kronecker(C,pred_obj_new$Sig)
  CS_lam <- CS + diag(pred_obj_new$Lambda/A_mult +
                        pred_obj_new$eps)
  
  pred_obj_new$Ki <- chol2inv(chol(CS_lam))
  
  return(cIMSPE_MV(x= NULL,
                   id = nrow(pred_obj_new$X0),
                   pred_obj = pred_obj_new,
                   reps = reps - 1,
                   Wijs = Wijs_new,
                   cv = cv))
  
}


cIMSPE_MV_new_X_batch_list <- function(x, pred_obj,
                                       reps = 2,
                                       scale_factor = 1,
                                       Wijs = NULL,
                                       cv = T){
  out <- 0
  
  for( i in 1:length(pred_obj)){
    out <- out + cIMSPE_MV_new_X_batch(x = x,
                                       pred_obj = pred_obj[[i]],
                                       reps = reps,
                                       scale_factor = scale_factor,
                                       Wijs = Wijs[[i]],
                                       cv = cv)
  }
  return(out)
  
}

IMSPE_MV_search <- function(pred_obj, replicate = FALSE, Xcand = NULL, 
                            nlopt_ctrl = list("algorithm" = "NLOPT_LN_BOBYQA",
                                              "maxeval" = 1e5,
                                              "xtol_abs" = 1e-8,
                                              "ftol_abs" = 1e-8),
                            maximim = F,
                            n_reps = 1,
                            tol_diff =1e-6, tol_dist =1e-6,
                            n_intial_points = 50,
                            scale_factor  = pi/2,
                            Wijs = NULL, seed = NULL,
                            n_cores = 1){
  # Only search on existing designs
  if(replicate){
    ## Discrete optimization
    res <- unlist(mclapply(1:nrow(model$X0), 
                           cIMSPE_MV,
                           Wijs = Wijs, 
                           pred_obj = pred_obj,
                           reps = n_reps,
                           x = NULL, mc.cores = ncores))
    
    return(list(par = model$X0[which.min(res),,drop = FALSE],
                value = min(res), 
                new = FALSE,
                id = which.min(res)))
  }
  
  
  n_dim <- ncol(pred_obj$X0)
  
  ## Precalculations
  if(is.null(Wijs)) Wijs <- Wij_mv(X = pred_obj$X0,
                                   theta = pred_obj$theta,
                                   cov_type = pred_obj$cov_type)
  
  
  ## Optimization
  if(is.null(Xcand)){
    ## Continuous optimization
    if(is.null(seed)) seed <- sample(1:2^10, 1) ## To be changed?
    
    Xstart <- DiceDesign::lhsDesign(n_intial_points, 
                                    n_dim, 
                                    seed = seed)$design
    if(maximin) Xstart <- DiceDesign::maximimSA_LHS(XStart)$design
    
    res <- list(par = NA, 
                value = Inf, 
                new = NA)
    
    local_opt_fun <- function(i){
      out <- try(nloptr(x0 = Xstart[i,, drop = FALSE],
                        eval_f = cIMSPE_MV, 
                        lb = rep(0, n_dim),
                        ub = rep(scale_factor, n_dim),
                        Wijs = Wijs, 
                        pred_obj = pred_obj,
                        opts = nlopt_ctrl))
      if(is(out, "try-error")) return(NULL)
      return(out)
    }
    
    all_res <- mclapply(1:nrow(Xstart), 
                        local_opt_fun,
                        mc.cores = n_cores)
    
    res_min <- which.min(Reduce(c, lapply(all_res, function(x) x$objective)))
    
    res <- list(par = apply(all_res[[res_min]]$solution,
                            c(1, 2),
                            function(x) max(min(x , 1), 0)),
                value = all_res[[res_min]]$objective,
                new = TRUE,
                id = NULL)
    
    if(control$tol_dist > 0 && control$tol_diff > 0){
      ## Check if new design is not to close to existing design
      dists <- dist(res$par, pred_obj$X0)
      if(min(dists) < tol_dist){
        res <- list(par = pred_obj$X0[which.min(dists),,drop = FALSE],
                    value = cIMSPE_MV(x = pred_obj$X0[which.min(dists),, drop = F],
                                      pred_obj = pred_obj, 
                                      id = which.min(dists),
                                      Wijs = Wijs),
                    new = FALSE, id = which.min(dists))
      }else{
        ## Check if IMSPE difference between replication and new design is significative
        id_closest <- which.min(dists) # closest point to new design
        imspe_rep <- cIMSPE_MV(pred_obj = pred_obj, 
                               id = id_closest, 
                               Wijs = Wijs)
        
        if((imspe_rep - res$objective)/res$objective < tol_diff){ #(1 - sum(model$Ki * Wijs)) < ){
          res <- list(par = pred_obj$X0[which.min(dists),,drop = FALSE],
                      value = imspe_rep,
                      new = FALSE, id = which.min(dists))
        }
      }
    }
    
    
    return(res)
    
    
  }else{
    ## Discrete optimization
    
    # redefine crit_IMSPE for use with mclapply
    cIMSPE_mcl <- function(i, pred_obj, Wijs, Xcand){
      cIMSPE_MV(x = Xcand[i,,drop = F], 
                Wijs = Wijs, 
                pred_obj = pred_obj)
    }
    # res <- unlist(apply(Xcand, 1, crit_IMSPE, Wijs = Wijs, model = model))
    res <- unlist(mclapply(1:nrow(Xcand),
                           cIMSPE_mcl,
                           Xcand = Xcand, 
                           Wijs = Wijs, 
                           pred_obj = pred_obj,
                           mc.cores = n_cores))
    
    tmp <- which(duplicated(rbind(pred_obj$X0, 
                                  Xcand[which.min(res),,drop = FALSE]),
                            fromLast = TRUE))
    
    if(length(tmp) > 0) return(list(par = Xcand[which.min(res),,drop = FALSE],
                                    value = min(res),
                                    new = FALSE, id = tmp))
    return(list(par = Xcand[which.min(res),,drop = FALSE], 
                value = min(res), new = TRUE, id = NULL))
  }
  
  
}


## Create a new one that updates Wijs for list and updates pred_objs list.

IMSPE_MV_batch_search <- function(pred_obj,  Xcand = NULL, 
                                  nlopt_ctrl = list("algorithm" = "NLOPT_LN_BOBYQA",
                                                    "maxeval" = 1e5,
                                                    "xtol_abs" = 1e-8,
                                                    "ftol_abs" = 1e-8),
                                  maximin = F,
                                  n_reps_batch = 3,
                                  tol_diff =1e-6, tol_dist =1e-6,
                                  n_intial_points = 50,
                                  scale_factor  = pi/2,
                                  Wijs = NULL,
                                  seed = NULL,
                                  Xstart = NULL,
                                  n_cores = 1,
                                  n_block = 1,
                                  n_dim = 4,
                                  cv = T){
  
  
  ## Precalculations
  
  if(is.null(Wijs)){
    if(n_block > 1){
      Wijs <- vector(mode = "list",
                     length = n_block)
      
      for(i in 1:n_block) Wijs[[i]] <- Wij_mv(X = pred_obj[[i]]$X0,
                                              scale_factor = scale_factor,
                                              theta = pred_obj[[i]]$theta,
                                              cov_type = pred_obj[[i]]$cov_type)
      
    }else{
      Wijs <- Wij_mv(X = pred_obj$X0,
                     scale_factor = scale_factor,
                     theta = pred_obj$theta,
                     cov_type = pred_obj$cov_type)
    }
  }
  
  # Add n_batch for existing designs
  
  if(n_block > 1){
    cIMSPE_MV_fun <- cIMSPE_MV_list
    n_obs <- nrow(pred_obj[[1]]$X0)
  }else{
    cIMSPE_MV_fun <- cIMSPE_MV
    n_obs <- nrow(pred_obj$X0)
  }
  
  ## Discrete optimization
  res_replicate <- unlist(mclapply(1:n_obs, 
                                   cIMSPE_MV_fun,
                                   Wijs = Wijs, 
                                   pred_obj = pred_obj,
                                   reps = n_reps_batch,
                                   scale_factor = scale_factor,
                                   cv = cv,
                                   x = NULL,
                                   mc.cores = 1))
  
  best_obj <- list(value = min(res_replicate), 
                   new = FALSE,
                   id = which.min(res_replicate),
                   id_cand = NULL)
  
  ## Optimization new batch
  cIMSPE_MV_new_X_batch_f <- cIMSPE_MV_new_X_batch
  
  if(n_block > 1) cIMSPE_MV_new_X_batch_f <- cIMSPE_MV_new_X_batch_list
  
  if(is.null(Xcand)){
    ## Continuous optimization
    
    if(is.null(Xstart)){
      if(is.null(seed)) seed <- sample(1:2^10, 1) ## To be changed?
      
      Xstart <- DiceDesign::lhsDesign(n_intial_points, 
                                      n_dim, 
                                      seed = seed)$design
      if(maximin) Xstart <- DiceDesign::maximinSA_LHS(Xstart)$design
      
      Xstart <- Xstart * scale_factor
    }
    
    local_opt_fun <- function(i){
      out <- try(nloptr(x0 = Xstart[i,, drop = FALSE],
                        eval_f = cIMSPE_MV_new_X_batch_f, 
                        lb = rep(0, n_dim),
                        ub = rep(scale_factor, n_dim),
                        Wijs = Wijs, 
                        scale_factor = scale_factor,
                        reps = n_reps_batch,
                        pred_obj = pred_obj,
                        cv = cv,
                        opts = nlopt_ctrl))
      if(is(out, "try-error")) return(NULL)
      return(out)
    }
    
    
    all_res <- mclapply(1:nrow(Xstart), 
                        local_opt_fun,
                        mc.cores = 1)
    
    
    res_X <- Reduce(rbind, lapply(all_res, function(x) x$solution))
    res_obj <- Reduce(c, lapply(all_res, function(x) x$objective))
    
    better_new <- which(res_obj < best_obj$value)
    
    if(best_obj$value - min(res_obj) > tol_diff ){
      res_min <- which.min(res_obj)
      best_obj <- list(par = all_res[[res_min]]$solution,
                       value = all_res[[res_min]]$objective,
                       new = TRUE,
                       id = NULL,
                       id_cand = NULL)
    }
    
    
    
    return(list(best_obj = best_obj,
                res_obj = res_obj,
                res_X = res_X,
                res_rep = res_replicate))
    
    # if(control$tol_dist > 0 && control$tol_diff > 0){
    #   ## Check if new design is not to close to existing design
    #   dists <- rdist::cdist(res_X,
    #                         pred_obj$X0)
    #   r_mins <- matrixStats::rowMins(dists)
    #   r_mins <- which(r_mins > tol_dist)
    #   
    #   if(min(dists) < tol_dist){
    #     res <- list(par = pred_obj$X0[which.min(dists),,drop = FALSE],
    #                 value = cIMSPE_MV(x = pred_obj$X0[which.min(dists),, drop = F],
    #                                   pred_obj = pred_obj, 
    #                                    id = which.min(dists),
    #                                    Wijs = Wijs),
    #                 new = FALSE, id = which.min(dists))
    #   }else{
    #     ## Check if IMSPE difference between replication and new design is significative
    #     id_closest <- which.min(dists) # closest point to new design
    #     imspe_rep <- cIMSPE_MV(pred_obj = pred_obj, 
    #                            id = id_closest, 
    #                            Wijs = Wijs)
    #     
    #     if((imspe_rep - res$objective)/res$objective < tol_diff){ #(1 - sum(model$Ki * Wijs)) < ){
    #       res <- list(par = pred_obj$X0[which.min(dists),,drop = FALSE],
    #                   value = imspe_rep,
    #                   new = FALSE, id = which.min(dists))
    #     }
    #   }
    # }
    
    
  }else{
    ## Discrete optimization
    
    # redefine crit_IMSPE for use with mclapply
    cIMSPE_mcl <- function(i, pred_obj, Wijs, Xcand){
      cIMSPE_MV_new_X_batch_f(x = Xcand[i,,drop = F], 
                              Wijs = Wijs, 
                              scale_factor = scale_factor,
                              reps = n_reps_batch,
                              cv = cv,
                              pred_obj = pred_obj)
    }
    # res <- unlist(apply(Xcand, 1, crit_IMSPE, Wijs = Wijs, model = model))
    res_new <- unlist(mclapply(1:nrow(Xcand),
                               cIMSPE_mcl,
                               Xcand = Xcand, 
                               Wijs = Wijs, 
                               pred_obj = pred_obj,
                               mc.cores = 1))
    
    best_obj <- list(par = Xcand[which.min(res_new),,drop = FALSE], 
                     value = min(res_new), 
                     new = TRUE,
                     id = NULL,
                     id_cand = which.min(res_new))
    
    return(list(best_obj = best_obj,
                res_obj = res_new,
                res_rep = res_replicate))
  }
  
  
}





IMSPE_max_batch<- function(pred_obj, h = 1, batch_size = 3,
                           Xcand = NULL, 
                           control = list(tol_dist = 1e-6, 
                                          tol_diff = 1e-6, 
                                          multi_start = 20,
                                          maxit = 100),
                           nlopt_ctrl = list("algorithm" = "NLOPT_LN_BOBYQA",
                                             "maxeval" = 1e5,
                                             "xtol_abs" = 1e-8,
                                             "ftol_abs" = 1e-8),
                           Wijs = NULL,
                           scale_factor = 1,
                           seed = NULL, n_cores = 1){
  ## Assume that the new design will be simulated batch_size times.
  
  d <- ncol(model$X0)
  
  ## Precalculations
  if(is.null(Wijs)) Wijs <- Wij_mv(X = pred_obj$X0,
                                   theta = pred_obj$theta,
                                   cov_type = pred_obj$cov_type)
  
  
  ## A) Setting to beat: first new point then replicate h times
  IMSPE_A <- IMSPE_MV_search(pred_obj = pred_obj, 
                             Xcand = Xcand, 
                             nlopt_ctrl = nlopt_ctrl,
                             n_reps = 1,
                             tol_diff =control$tol_diff,
                             tol_dist =control$tol_dist,
                             n_intial_points = control$multi_start,
                             scale_factor  = scale_factor,
                             Wijs = Wijs, 
                             seed = seed,
                             n_cores = n_cores)
  new_designA <- IMSPE_A$par ## store first considered design to be added
  path_A <- list(IMSPE_A)
  
  if(h > 0){
    newmodelA <- model
    
    if(IMSPE_A$new){
      newWijs <- Wij(mu1 = model$X0,
                     mu2 = new_designA,
                     theta = model$theta, 
                     type = model$covtype)
      WijsA <- cbind(Wijs, newWijs)
      WijsA <- rbind(WijsA, c(newWijs, Wij(IMSPE_A$par,
                                           IMSPE_A$par,
                                           theta = model$theta,
                                           type = model$covtype)))
    }else{
      WijsA <- Wijs
    }
    
    for(i in 1:h){
      newmodelA <- update(object = newmodelA, 
                          Xnew = IMSPE_A$par,
                          Znew = NA, 
                          maxit = 0)
      IMSPE_A <- IMSPE.search(model = newmodelA,
                              replicate = TRUE,
                              control = control,
                              Wijs = WijsA, ncores = ncores)
      path_A <- c(path_A, list(IMSPE_A))
    }
    
  }
  
  if(h == -1) return(list(par = new_designA, value = IMSPE_A$value, path = path_A)) 
  
  ## B) Now compare with waiting to add new point
  newmodelB <- model
  
  if(h == 0){
    IMSPE_B <- IMSPE.search(model = newmodelB, replicate = TRUE, control = control, Wijs = Wijs, ncores = ncores)
    new_designB <- IMSPE_B$par ## store considered design to be added
    
    # search from best replicate
    if(is.null(Xcand)){
      IMSPE_C <- IMSPE.search(model = newmodelB, Wijs = Wijs,
                              control = list(Xstart = IMSPE_B$par, maxit = control$maxit,
                                             tol_dist = control$tol_dist, tol_diff = control$tol_diff), ncores = ncores)
    }else{
      IMSPE_C <- IMSPE_B
    }
    
    if(IMSPE_C$value < min(IMSPE_A$value, IMSPE_B$value)) return(list(par = IMSPE_C$par, value = IMSPE_C$value, path = list(IMSPE_C)))
    
    if(IMSPE_B$value < IMSPE_A$value){
      return(list(par = IMSPE_B$par, value = IMSPE_B$value, path = list(IMSPE_B)))
    } 
  }else{
    for(i in 1:h){
      ## Add new replicate
      IMSPE_B <- IMSPE.search(model = newmodelB, replicate = TRUE, control = control, Wijs = Wijs, ncores = ncores)
      
      if(i == 1){
        new_designB <- matrix(IMSPE_B$par, nrow = 1) ##store first considered design to add
        path_B <- list()
      } 
      
      path_B <- c(path_B, list(IMSPE_B))
      newmodelB <- update(object = newmodelB, Xnew = IMSPE_B$par, Znew = NA, maxit = 0)
      
      ## Add new design
      IMSPE_C <- IMSPE.search(model = newmodelB, control = control, Xcand = Xcand, Wijs = Wijs, ncores = ncores)
      path_C <- list(IMSPE_C)
      
      if(i < h){
        newmodelC <- newmodelB
        
        if(!any(duplicated(rbind(model$X0, IMSPE_C$par)))){ 
          newWijs <- Wij(mu1 = model$X0, mu2 = IMSPE_C$par, theta = model$theta, type = model$covtype)
          WijsC <- cbind(Wijs, newWijs)
          WijsC <- rbind(WijsC, c(newWijs, Wij(mu1 = IMSPE_C$par, theta = model$theta, type = model$covtype)))
        }else{
          WijsC <- Wijs
        }
        
        for(j in i:(h-1)){
          ## Add remaining replicates
          newmodelC <- update(object = newmodelC, Xnew = IMSPE_C$par, Znew = NA, maxit = 0)
          IMSPE_C <- IMSPE.search(model = newmodelC, replicate = TRUE, control = control, Wijs = WijsC, ncores = ncores)
          path_C <- c(path_C, list(IMSPE_C))
        }
      }
      
      if(IMSPE_C$value < IMSPE_A$value) return(list(par = new_designB, value = IMSPE_C$value, path = c(path_B, path_C)))
    }
  }
  
  return(list(par = new_designA, value = IMSPE_A$value, path = path_A))
  
}
