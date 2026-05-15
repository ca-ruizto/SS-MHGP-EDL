### Legendre Polynomials ####

compute_coefs_leg_bf <- function(n_funs = 3){
  coefMat <- sapply(seq_len(n_funs), function(n) {
    coef <- rep(0, n_funs)
    # # coef[1] <- (-1)^(n - 1)
    for (k in seq_len(n)) {
      coef[k] <- (-1)^(n - k) * choose(n - 1, k - 1) * 
        choose(n + k - 2, k - 1)
    }
    coef * sqrt(2 * n - 1)
  })
  
  return(coefMat)
}

create_leg_basis_fun <- function(pts,n_funs = 3, intercept = T){
  #modified from fdapace::CreateBasis
  
  if (n_funs == 1) {
    res <- matrix(1, length(pts), 1)
  } else if (n_funs > 1) {
    coefMat <- compute_coefs_leg_bf(n_funs)
    xMat <- cbind(1, stats::poly(pts, n_funs - 1, raw=TRUE))
    res <- xMat %*% coefMat
    
    if(! intercept) res <- res[,-1]
  }
  
  return(res)
}


### Other basis #####


create_univ_basis_fun <- function(x, nf = 3,
                                  knots = NULL, type = c("bs,ns"), 
                                  intercept = F,
                                  sparse_mat = T){
  
  if(is.null(knots)) knots <- seq(min(x), max(x), length.out = nf)
  
  n_b <- length(knots)
  int_knots <- knots[-c(1,n_b)]
  bknots <- knots[c(1,n_b)]
  
  if(type == "bs"){
    out <- splines::bs(x, knots = int_knots,
                       intercept = intercept,
                       Boundary.knots = bknots)
  }else{
    out <- splines::ns(x, knots = int_knots, 
                       intercept = intercept,
                       Boundary.knots = bknots)
  }
  
  attr(out, "class") <- "matrix"
  
  if(sparse_mat) out <- as(out, "sparseMatrix")
  return(out)
}

mgcv_tps_fun <- function(obs, n_out,
                         prop_names = paste("p", 1:n_out, sep="_"),
                         k_s = 5, max_val = 1){
  var_names <- colnames(obs)
  
  k_list <- vector(mode = "list")
  
  for(i in 1:length(var_names) ) k_list[[var_names[i]]] <- c(0,max_val)
  
  s_eval <- mgcv::s(var_names,
                    bs = "tp",
                    k = k_s)
  
  s_eval$term <- var_names
  s_eval$dim <- length(var_names)
  
  XX <- mgcv::smooth.construct2(
    s_eval,
    knots = k_list,
    data = obs
  )
  
  out <- XX$X
  colnames(out) <- paste("bf", 1:ncol(out),
                         sep = "_")
  
  if(n_out > 1) {
    c_names <- colnames(out)
    out <- kronecker(out, Diagonal(n_out))
    colnames(out) <- as.vector(sapply(c_names, 
                                      function(x) paste(prop_names,
                                                        x,sep = "_")))
  }
  return(out)
}

create_tensor_bf <- function(X, n_out = 2,
                             prop_names = paste("p", 1:n_out, sep="_"),
                             var_names = paste("x",1:ncol(X), sep="_"),
                             knots_list, 
                             type = "bs",
                             intercept = TRUE,
                             sparse_mat = T){
  n_vars <- ncol(X)
  if(length(type) == 1) type <- rep(type, n_vars)
  
  int <- kronecker(rep(1,nrow(X)),
                   Diagonal(n_out))
  
  colnames(int) <- paste(prop_names, "intercept", sep = "_")
  
  out <- create_univ_basis_fun(X[,1],
                               knots = knots_list[[1]],
                               type = type[1],
                               intercept = intercept,
                               sparse_mat =  sparse_mat)
  
  if(intercept) out <- cbind(1, out)
  
  for(i in 2:n_vars){
    new_b <- create_univ_basis_fun(X[,i],
                                   knots = knots_list[[i]],
                                   type = type[i],
                                   sparse_mat =  sparse_mat)
    
    if(intercept) new_b <- cbind(1, new_b)
    nc_out <- ncol(out)
    
    out <- kronecker( t(rep(1, ncol(new_b))), out)
    new_b <- kronecker(new_b, t(rep(1,nc_out)))
    out <- out * new_b

    
  }
  
  colnames(out) <- paste("bf", 1:ncol(out),
                         sep = "_")
  
  return(out)
}

create_mv_basis_fun <- function(X, n_out,
                                prop_names = paste("p", 1:n_out, sep="_"),
                                var_names = paste("x",1:ncol(X), sep="_"),
                                knots_list,
                                type = "bs",
                                func_intercept = T,
                                omit_last_var = T,
                                sparse_mat = T){
  
  n_vars <- ncol(X) - omit_last_var
  if(length(type) == 1) type <- rep(type, n_vars)
  
  out <- kronecker(rep(1,nrow(X)),
                   Diagonal(n_out))
  
  colnames(out) <- paste(prop_names, "intercept", sep = "_")
  
  
  for(i in 1:n_vars){
    new_b <- create_univ_basis_fun(X[,i],
                                   knots = knots_list[[i]],
                                   type = type[i],
                                   intercept = func_intercept,
                                   sparse_mat =  sparse_mat)
    
    c_names <- paste(var_names[i],
                     1:ncol(new_b), sep = "_")
    
    new_b <- kronecker(new_b, Diagonal(n_out))
    
    colnames(new_b) <- as.vector(sapply(c_names, 
                                      function(x) paste(prop_names,
                                                        x,sep = "_")))
    
    out <- cbind(out, new_b)

  }
  
  return(out)
  
}




create_poly_basis <- function(X, n_out = 2, 
                              prop_names = paste("p", 1:n_out, sep="_"),
                              degree = 1, 
                              omit_last_var = T,
                              sparse_mat = T,
                              var_names = paste("x",1:ncol(X), sep="_")){
  
  int <- kronecker(rep(1,nrow(X)),
                   Diagonal(n_out))
  
  colnames(int) <- paste(prop_names,
                         "intercept", 
                         sep = "_")
  
  out <- stats::poly(X[,1],degree, raw = T)
  colnames(out) <- paste(var_names[1],
                         1:degree, sep = "_")
  
  n_vars <- ncol(X) - omit_last_var
  
  for(i in 2:n_vars){
    temp <- stats::poly(X[,i],degree, raw = T)
    colnames(temp) <- paste(var_names[i],
                           1:degree, sep = "_")
    
    out <- cbind(out, temp)
  }
  
  c_names <- colnames(out)
  out <- kronecker(out, Diagonal(n_out))
  colnames(out) <- as.vector(sapply(c_names, 
                                    function(x) paste(prop_names,
                                                      x,sep = "_")))
  
  out <- cbind(int , out)
  
  if(sparse_mat) out <- as(out, "sparseMatrix")
  
  return(out)
}

