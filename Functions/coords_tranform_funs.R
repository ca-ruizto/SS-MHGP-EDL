# 1. Helper: Create the Orthonormal ILR Basis Matrix
# create_ilr_basis <- function(D) {
#   V <- matrix(0, nrow = D, ncol = D - 1)
#   for (k in 1:(D - 1)) {
#     V[1:k, k] <- -1 / sqrt(k * (k + 1))
#     V[k + 1, k] <- sqrt(k / (k + 1))
#   }
#   return(V)
# }

create_ilr_basis <- function(D, method = "basic"){
  compositions::ilrBase(D = D, method = method)
}

# 1. Additive Zero Replacement (Epsilon Smoothing)
additive_zero_replace <- function(x, epsilon) {
  D <- ncol(x)
  # Add epsilon to all elements and re-normalize
  x_new <- (x + epsilon) / (1 + D * epsilon)
  return(x_new)
}

ilr_bounded_minmax_additive <- function(x, b = 0, epsilon = 1e-4, 
                                        z_min = NULL, z_max = NULL, method = "balanced") {
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  D <- ncol(x)
  
  # Apply Additive Replacement uniformly
  x_replaced <- additive_zero_replace(x, epsilon)
  
  # Shift and scale to standard simplex w
  x_shifted <- x_replaced - b
  if (any(x_shifted <= 0)) {
    stop("epsilon is too small relative to lower bound 'b'. Ensure epsilon / (1 + D*epsilon) > b")
  }
  
  free_space <- 1 - (D * b)
  w <- x_shifted / free_space
  
  # Apply standard ILR
  V <- create_ilr_basis(D, method)
  z <- log(w) %*% V
  
  # Min-Max Scaling
  if (is.null(z_min)) z_min <- apply(z, 2, min)
  if (is.null(z_max)) z_max <- apply(z, 2, max)
  
  range_z <- z_max - z_min
  range_z[range_z == 0] <- 1 
  
  u <- sweep(z, 2, z_min, "-")
  u <- sweep(u, 2, range_z, "/")
  
  return(list(
    u = u,
    z_min = z_min,
    z_max = z_max,
    replaced_x = x_replaced 
  ))
}

# 4. Backward Transformation: [0,1] Hypercube -> Bounded Simplex (with Additive Inverse)
inv_ilr_bounded_minmax_additive <- function(u, b = 0, epsilon = 1e-4,
                                            z_min, z_max) {
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  D <- ncol(u) + 1 
  free_space <- 1 - (D * b)
  
  # 1. Un-scale Min-Max to get back to unconstrained ILR space (z)
  range_z <- z_max - z_min
  range_z[range_z == 0] <- 1
  
  z <- sweep(u, 2, range_z, "*")
  z <- sweep(z, 2, z_min, "+")
  
  # 2. Inverse ILR to standard simplex (w)
  V <- create_ilr_basis(D)
  log_w <- z %*% t(V)
  w_unnormalized <- exp(log_w)
  w <- w_unnormalized / rowSums(w_unnormalized)
  
  # 3. Scale and shift back to the restricted bounded simplex (x_replaced)
  x_replaced <- (w * free_space) + b
  
  # 4. Invert the additive epsilon smoothing
  x <- (x_replaced * (1 + D * epsilon)) - epsilon
  
  return(x)
}