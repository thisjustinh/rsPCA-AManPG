set.seed(1)
z <- normalize(mvrnorm(n=30, mu=rep(0, 10), Sigma=sigma))
l1 <- 0.1
k <- 2
tol=1e-5
maxiter=1000
fmax=1e6
type=0 
gamma=0.5
verbose=TRUE

rspca.amanpg <- function(z, l1, k, tol=1e-5, maxiter=1000, fmax=1e6, type=0, 
                         gamma=0.5, normalize=TRUE, verbose=FALSE) {
  start <- Sys.time()
  
  if (normalize) {
    z <- normalize(z)
  }
  
  dims <- dim(z)
  n <- dims[1]
  p <- dims[2]
  
  if (p < n * 2) {
    z <- t(z) %*% z
    type <- 1
  }
  
  ### Initialized values ###
  init_svd <- svd(z, nu=k, nv=k)
  a <- normalize(init_svd$u, center=FALSE)
  b <- init_svd$d * init_svd$v
  
  ata <- t(a) %*% a
  btb <- t(b) %*% b
   
  # Get Lipschitz
  # TODO: Need to be different depending on what matrix I use?
  stepsize_a <- 1 / (2 * svd(ata)$d[1])
  stepsize_b <- 1 / (2 * svd(btb)$d[1])
  
  linesearch_flag <- 1
  total_linesearch <- 0
  min_step <- 0
  alpha <- 100 / p
  
  ### Main Loop ###
  for (iter in 1:maxiter) {
    if (verbose) {
      iter_start <- Sys.time()
      print("=================")
      print(paste("On iteration", iter))
    }
    
    ### Update A ###
    # s = 2*t*(z-a %*% t(b)) %*% b
    grad_a <- 2 * stepsize_a * (z %*% b - a %*% btb)  # well it's negative grad a times 2
    da <- grad_a - (a %*% t(a) %*% grad_a + a %*% t(grad_a) %*% a) / 2
    print(dim(da))
    
    ## Retract ##
    if (!linesearch_flag) alpha <- alpha * 1.1
    if (min_step) alpha <- 1 / d
    linesearch_flag <- 0
    min_step <- 0
    
    xi <- alpha * da
    retracted_a <- polar_decomp(a, xi)
    
    while (norm(z-retracted_a %*% t(b))^2 - norm(z-a %*% t(b))^2> -(alpha/(2*stepsize_a)) * norm(da)^2) {
      alpha <- gamma * alpha
      if (alpha < 1e-5 / d) {
        min_step <- 1
        break
      }
      
      linesearch_flag <- 1
      total_linesearch <- total_linesearch + 1
      xi <- alpha * da
      retracted_a <- polar_decomp(a, xi)
    }
    
    a <- retracted_a
    ata <- t(a) %*% a
    
    ## Update B ##
    c <- b + 2 * stepsize_b * (t(z) %*% a - b %*% ata)
    db <- c - b + stepsize_b * sign(c) * matrix(as.numeric(abs(c) > l1), nrow=d, ncol=k)  # done element-wise
    b <- db
    btb <- t(b) %*% b
    
    check <- norm(da/stepsize_a)^2 + norm(db/stepsize_b)^2
    
    if (verbose) {
      print(paste("sum:", check))
      print(paste("time:", difftime(Sys.time(), iter_start)))
    }
    
    # Check for convergence
    if (check <= tol^2 && check < fmax) {
      if (verbose) print(paste("Final difference": check))
      break
    }
  }
}

polar_decomp <- function(x, xi) {
  eigendecomp_xi <- eigen(diag(1, dim(xi)[2]) + t(xi) %*% xi, 
                          symmetric = TRUE)
  u <- eigendecomp_xi$vectors
  half_power <- u %*% diag(sqrt(1 / eigendecomp_xi$values)) %*% t(u)
  return((x + xi) %*% half_power)
}

normalize <- function(x, center=TRUE, scale=TRUE) {
  if (center) {
    x <- scale(x, scale=FALSE)
  }
  
  if (scale) {
    if (is.matrix(x)) {
      dims <- dim(x)
      x <- x / matrix(sqrt(rowSums(x^2)), nrow=dims[1], ncol=dims[2])
    } else if (is.numeric(x)) {
      x <- x / sqrt(sum(x^2))
    } else {
      stop("'x' must be a numeric vector or matrix")
    }
  }
  
  return(x)
}
