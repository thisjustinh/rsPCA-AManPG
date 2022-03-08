set.seed(1)
z <- normalize(mvrnorm(n=30, mu=rep(0, 10), Sigma=sigma))
l1 <- 1
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
  a <- init_svd$u / norm(init_svd$u, 'F')
  b <- init_svd$d[1] * init_svd$v

  ata <- t(a) %*% a
  btb <- t(b) %*% b
   
  # Get Lipschitz
  # stepsize_a <- 0.01
  stepsize_a <- 1 / (2 * svd(ata)$d[1])
  stepsize_b <- 1 / (2 * svd(btb)$d[1])
  
  linesearch_flag <- 1
  total_linesearch <- 0
  min_step <- 0
  alpha <- 100 / p
  
  grad_a <- 2 * stepsize_a * (z - a %*% t(b)) %*% b

  f = norm(z - a %*% t(b), 'F')^2 + l1 * sum(abs(b))
  
  ### Main Loop ###
  for (iter in 2:maxiter) {
    if (verbose) {
      iter_start <- Sys.time()
      print("=================")
      print(paste("On iteration", iter))
    }
    
    ### Update A ###
    # s = 2*t*(z-a %*% t(b)) %*% b
    gtx <- t(grad_a) %*% a
    da <- grad_a - a %*% (gtx + t(gtx)) / 2

    ## Retract ##
    if (!linesearch_flag) alpha <- alpha * 1.1
    if (min_step) alpha <- 1 / p
    linesearch_flag <- 0
    min_step <- 0
    
    xi <- alpha * da
    retracted_a <- polar_decomp(a, xi)
    norm_retract <- norm(z - retracted_a %*% t(b), 'F')^2
        
    while (norm_retract > norm(z - a %*% t(b), 'F')^2 - (alpha / (2*stepsize_a)) * norm(da, 'F')^2) {
      alpha <- gamma * alpha
      if (alpha < 1e-5 / p) {
        min_step <- 1
        break
      }

      linesearch_flag <- 1
      total_linesearch <- total_linesearch + 1
      xi <- alpha * da
      retracted_a <- polar_decomp(a, xi)
      norm_retract <- norm(z - retracted_a %*% t(b), 'F')^2 
    }
    
    a <- retracted_a
    ata <- t(a) %*% a
    
    ## Update B ##
    grad_b <- 2 * stepsize_b * t(z - a %*%  t(b)) %*% a
    c <- b + grad_b
    tmp <- abs(c) - as.vector(l1 * stepsize_b)
    if (k < 15) act_set <- as.numeric(tmp > 0) else act_set <- tmp > 0
    db <-  tmp * act_set * sign(c) - b
    # db <- matrix(1, nrow=p, ncol=k)
    # for (i in 1:p) {
    #   for (j in 1:k) {
    #     db[i,j] <- stepsize_b * max(abs(c[i,j])-l1, 0) * sign(c[i,j])-b[i,j]
    #   }
    # }
    # db <- c - b - stepsize_b * sign(c) * matrix(as.numeric(abs(c) > l1), nrow=p, ncol=k)  # done element-wise
    
    #tmp <- abs(c) - as.vector(l1)
    #if (k < 15) act_set <- as.numeric(tmp > 0) else act_set <- tmp > 0
    #db <- act_set * sign(c) * stepsize_b
    
    b <- b + db
    btb <- t(b) %*% b
    
    ### Check for convergence ###
    # check <- norm(da/stepsize_a)^2 + norm(db/stepsize_b)^2
    prev_grad <- -grad_a/stepsize_a - grad_b/stepsize_b
    # grad_a <- 2 * stepsize_a * (z %*% b - a %*% btb)
    # grad_b <- 2 * stepsize_b * (t(z) %*% a - b %*% ata)
    grad_a <- 2 * stepsize_a * (z - a %*% t(b)) %*% b
    grad_b <- 2 * stepsize_b * t(z - a %*%  t(b)) %*% a
    # check <- norm(-grad_a/stepsize_a - grad_b/stepsize_b, 'F')^2 # - norm(prev_grad, 'F')
    
    f <- c(f, norm(z - a %*% t(b), 'F')^2 + l1 * sum(abs(b)))
    check <- abs(f[iter] - f[iter-1])
    
    if (verbose) {
      print(paste("difference:", check))
      print(paste("time:", difftime(Sys.time(), iter_start)))
      print(paste("loss:", f[iter]))
    }
    
    # if (abs(check) <= tol || check > fmax) {
    #   if (verbose) print(paste("Final difference:", check))
    #   break
    # }
    if(check <= tol || check > fmax) {
      if (verbose) print(paste("Final difference: ", check))
      break
    }
  }
  
  b_norm <- sqrt(colSums(b^2))
  b_norm[b_norm == 0] <- 1
  
  return(list(
    iterations=iter,
    scores=a,
    b=b,
    loadings=b / matrix(1, p, 1) %*% b_norm,
    sparsity=sum(b == 0) / (p * k),
    time=difftime(Sys.time(), start)
  ))
}

# Polar decomposition as retraction method
# Pretty sure this is correct
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

# sprout <- rspca.amanpg(x, 1, 2, verbose=TRUE, maxiter=2000)
# sprout
