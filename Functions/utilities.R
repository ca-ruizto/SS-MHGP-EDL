mod_coefs_to_df_fun <- function(mod,
                                ss_dat_sort,
                                intercept = T){
  
  #returns a dataframe of values
  # first log(max_strain)
  # then model coefficients
  
  mod_coefs <- coef(mod)
  alloy_names <- levels(ss_dat_sort$alloy)
  
  max_strain_dat <- ss_dat_sort%>%
                          arrange(alloy)%>%
                          group_by(alloy)%>%
                          slice_head(n=1)
  
  n_basis = length(mod_coefs)/length(alloy_names)
  
  if(intercept){
    coef_names_v <- paste("bc", 1:(n_basis-1), sep = "_")
    coef_names_v <- c("intercept", coef_names_v)
  }else{
    coef_names_v <- paste("bc", 1:n_basis, sep = "_")
  }
  
  vals <- c(log(max_strain_dat$max_strain),
            unname(mod_coefs))
  coef_names_v <- c("log_max_strain", coef_names_v)
  
  return(tibble(val = vals,
                coef_names =  rep(coef_names_v, 
                                  each = length(alloy_names)),
                alloy = rep(alloy_names,
                            length(coef_names_v))))
  
}

mod_all_coefs_to_df_fun <- function(mod,
                                ss_dat_sort,
                                intercept = T){
  
  #returns a dataframe of values
  # first log(max_strain)
  # then model coefficients
  
  mod_coefs <- coef(mod)
  alloy_names <- levels(ss_dat_sort$alloy_rep_id)
  
  max_strain_dat <- ss_dat_sort%>%
    arrange(alloy_rep_id)%>%
    group_by(alloy_rep_id)%>%
    slice_head(n=1)
  
  n_basis = length(mod_coefs)/length(alloy_names)
  
  if(intercept){
    coef_names_v <- paste("bc", 1:(n_basis-1), sep = "_")
    coef_names_v <- c("intercept", coef_names_v)
  }else{
    coef_names_v <- paste("bc", 1:n_basis, sep = "_")
  }
  
  vals <- c(log(max_strain_dat$max_strain),
            unname(mod_coefs))
  coef_names_v <- c("log_max_strain", coef_names_v)
  
  return(tibble(val = vals,
                coef_names =  rep(coef_names_v, 
                                  each = length(alloy_names)),
                alloy_rep_id = rep(alloy_names,
                            length(coef_names_v))))
  
}

create_leg_basis_fun <- function(pts,n_funs = 3, intercept = T){
  #modified from fdapace::CreateBasis
  
  if (n_funs == 1) {
    res <- matrix(1, length(pts), 1)
  } else if (n_funs > 1) {
    coefMat <- sapply(seq_len(n_funs), function(n) {
      coef <- rep(0, n_funs)
      # # coef[1] <- (-1)^(n - 1)
      for (k in seq_len(n)) {
        coef[k] <- (-1)^(n - k) * choose(n - 1, k - 1) * 
          choose(n + k - 2, k - 1)
      }
      coef * sqrt(2 * n - 1)
    })
    # 
    # if(intercept){
    #   xMat <- cbind(1, stats::poly(pts, n_funs - 1, raw=TRUE))
    # }else{
    #   xMat <- stats::poly(pts, n_funs - 1, raw=TRUE)
    #   coefMat <- coefMat[,-1]
    # }
    xMat <- cbind(1, stats::poly(pts, n_funs - 1, raw=TRUE))
    res <- xMat %*% coefMat
    
    if(! intercept) res <- res[,-1]
  }
    
    return(res)
}


fit_leg_basis_mod_fun <- function(ss_dat,
                                    n_funs = 5,
                                  lam = 1,
                                  intercept = T,
                                  return_pred = F){
  #fits a linear model using orthogonal legendred basis
  #fits for every group
  
 
  alr_id <- levels(ss_dat$alloy_rep_id)
  

  
  
  if(intercept){
    fun_names <- paste("bc", 1:n_funs,
                       sep = "_")
    k_v <-  0:(n_funs -1)
  }else{
    fun_names <- paste("bc", 1:(n_funs-1),
                       sep = "_")
    k_v <- 1:(n_funs-1) 
  }
  
  k_v <- lam * (2*k_v + 1)*(k_v - 1) *(k_v)*(k_v + 1)*(k_v + 2)

  
  out <- out_pred <-  vector(mode = "list",
                length = length(alr_id))
  
  for(i in 1:length(alr_id)){
    temp_dat <- ss_dat%>%
      filter(alloy_rep_id==alr_id[i])
    
    basis_mat <- create_leg_basis_fun(n_funs = n_funs,
                                pts = temp_dat$scaled_strain,
                                intercept = intercept)
    colnames(basis_mat) <- fun_names
    
    reg_dat <- as_tibble(basis_mat)%>%
      mutate(scaled_stress = temp_dat$scaled_stress)
    
    reg_mod <- lm(scaled_stress ~ -1 + .,
                  data = reg_dat)
    
    #thoughness
    
    pred_stress <- fitted.values(reg_mod)
    ds <- diff(c(0,temp_dat$strain))
    ds[c(1,length(ds))] <- ds[c(1,length(ds))]/2
    th <- sum(pred_stress*ds)
    
    #young modulus
    idx <- which.min(abs(temp_dat$strain - 0.01))
    idx <- c(idx- 1, idx + 1)
    dstrain <- temp_dat$strain[idx]
    dstress <- temp_dat$stress[idx]
    
    ym <- (pred_stress[idx[2]] - pred_stress[idx[1]])/(dstrain[2] - dstrain[1])
    
    n_obs <- nrow(temp_dat)
    kt <- 1 +  k_v / n_obs
    edf <- sum(1/kt)
    gcv_val <- sqrt(sum(residuals(reg_mod)^2))
    #gcv_val <- gcv_val * n_obs / (n_obs - edf)
    
    #return value
    coefs_v <- c(temp_dat$max_stress[1],
                 temp_dat$max_strain[1])
    coefs_v <- c(log(coefs_v),
                 unname(coef(reg_mod)),
                 th, ym,pred_stress[length(pred_stress)],
                 gcv_val)
    
    out[[i]] <- tibble(val = coefs_v,
                       coef_names = c("log_max_stress",
                                      "log_max_strain", 
                                      fun_names,
                                      "toughness",
                                      "young_mod",
                                      "frac_stress",
                                      "gcv"),
                       alloy_rep_id = alr_id[i])
    
    out_pred[[i]] <- temp_dat%>%
                      mutate(pred_stress = pred_stress)
  }
  
  if(return_pred) return(list(coefs = bind_rows(out),
                              preds = bind_rows(out_pred)))
  
  return(bind_rows(out))
  
}





