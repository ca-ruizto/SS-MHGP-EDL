#### Kernels ######


matern_5_2_fun <- function(D){
  Dt <- D*sqrt(5)
  exp(-Dt)*(1 + Dt + (Dt^2)/3)
} 

matern_3_2_fun <- function(D){
  Dt <- D*sqrt(3)
  exp(-Dt)*(1 + Dt )
} 

exp_ker_fun <- function(D) exp(-D)

gaus_ker_fun <- function(D) exp(-D^2 / 2)

### Distances functions #####

mv_dist_fun <- function(X, X2 = NULL){
  #compute the distance between points in two matrices for each column
  n_c <- ncol(X)
  
  out <- vector(mode = "list",
                length = n_c)
  
  if(is.null(X2)) X2 <- X
  
  for( i in 1:n_c){
    out[[i]] <- rdist::cdist(X[,i],X2[,i])
  }
  
  return(out)
}

mv_dist_vec_fun <- function(x, x2 = NULL){
  n_c <- length(x)
  
  out <- vector(mode = "list",
                length = n_c)
  
  if(is.null(x2)) x2 <- x
  
  for( i in 1:n_c){
    out[[i]] <- rdist::cdist(x[i],x2[i])
  }
  
  return(out)
}


general_cor_fun <- function(theta,n_block = 1,D_mat,kernel_fun){
  #computes the correlation matrix for the paper (C=K and Sig = Omega in the paper notation)
  n_dim <- length(theta)/n_block
  exp_th <- exp(theta)
  C <- kernel_fun(D_mat[[1]]/exp_th[1])
  for(i in 2:n_dim) C <-C * kernel_fun(D_mat[[i]]/exp_th[i])
  
  if(n_block == 1) return(C)
  
  C_out <- vector(mode = "list",
                  length = n_block)
  
  C_out[[1]] <- C
  
  for(j in 2:n_block){
    idx <- (j-1)*n_dim
    
    C <- kernel_fun(D_mat[[1]]/exp_th[1 + idx])
    for(i in 2:n_dim) C <-C * kernel_fun(D_mat[[i]]/exp_th[i + idx])
    
    C_out[[j]] <- C
  }
  
  return(C_out)
}



additive_cor_fun <- function(theta,D_mat,n_block=1,kernel_fun){
  #computes the correlation matrix for the paper (C=K and Sig = Omega in the paper notation)
  n_dim <- length(theta)/n_block
  exp_th <- exp(theta)
  
  D <- D_mat[[1]]/exp_th[1]
  for(i in 2:n_dim) D <- D + D_mat[[i]] / exp_th[i]
  
  return(kernel_fun(D))
}
