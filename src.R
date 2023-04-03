rspca.amanpg <- function(z, l1, k, tol=1e-5, maxiter=1000, fmax=1e6, type=0, 
                         gamma=0.5, normalize=TRUE, verbose=FALSE) {
  start <- Sys.time()
  
  if (normalize) {
    z <- normalize(z)
  }
  
  dims <- dim(z)
  n <- dims[1]
  p <- dims[2]
  
  if (p < n * 2 && n > 100) {
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
  stepsize_a <- 1 / (2 * svd(btb)$d[1])
  stepsize_b <- 1 / (2 * svd(ata)$d[1])
  # print(stepsize_a)
  # print(stepsize_b)
  stepsize_a <- 1.2
  stepsize_b <- 0.05

  linesearch_flag <- 1
  total_linesearch <- 0
  min_step <- 0
  alpha <- 100 / p
  
  # grad_a <- 2 * stepsize_a * (z - a %*% t(b)) %*% b
  # grad_b <- 2 * stepsize_b * t(z - a %*%  t(b)) %*% a
  grad_a <- 2 * stepsize_a * (z %*% b - a %*% btb)
  grad_b <- 2 * stepsize_b * (t(z) %*% a - b %*% ata)

  f = norm(z - a %*% t(b), 'F')^2 + l1 * sum(abs(b))

  ### Main Loop ###
  for (iter in 2:maxiter) {
    if (verbose) {
      iter_start <- Sys.time()
      print("=================")
      print(paste("On iteration", iter))
    }
    
    ### Update A ###
    # s = 2*stepsize_a*grad_a
    # da <- s - (a*t(a)*s+a*t(s)*a) / 2
    gtx <- t(grad_a) %*% a
    da <- grad_a - a %*% (gtx + t(gtx)) / 2
    # da <- grad_a - (a %*% t(a) %*% grad_a + a %*% t(grad_a) %*% a) / 2

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
    c <- b + grad_b
    # tmp <- abs(c) - as.vector(l1 * stepsize_b)
    tmp <- abs(c) - as.vector(l1)
    if (k < 15) act_set <- as.numeric(tmp > 0) else act_set <- tmp > 0
    db <-  tmp * act_set * sign(c) - b
    
    b <- b + db
  
    btb <- t(b) %*% b
    
    ### Check for convergence ###
    grad_a <- 2 * stepsize_a * (z %*% b - a %*% btb)
    grad_b <- 2 * stepsize_b * (t(z) %*% a - b %*% ata)

    f <- c(f, norm(z - a %*% t(b), 'F')^2 + l1 * sum(abs(b)))
    check <- abs(f[iter] - f[iter-1])
    # check <- norm(grad_a / stepsize_a + grad_b / stepsize_b)
    # check <- max(1 - diag(abs(crossprod(b, oldb))))
    
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
      # if (verbose) print(paste("Final difference: ", check))
      break
    }
  }
  
  b_norm <- sqrt(colSums(b^2))
  b_norm[b_norm == 0] <- 1
  loadings <- b / matrix(1, p, 1) %*% b_norm
  
  # # Calculate PEV and CPEV
  # total <- var()
  # pev <- 
  # 
  # scores <- 
    
  return(list(
    iterations=iter,
    scores=a,
    b=b,
    loadings=loadings,
    sparsity=sum(b == 0) / (p * k),
    time=difftime(Sys.time(), start),
    loss=f
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

scree.plot <- function(x, data, normalized=TRUE) {
  if (!normalized) {
    sum <- sum(apply(data, 2, var))
  } else {
    sum <- 1
  }
  
  ev <- apply(x$scores, 2, var)
  pev <- ev / sum
  
  qplot(c(1:ncol(x$loadings)), pev) + 
    geom_line() +
    ylim(0,1) +
    xlab("PC") +
    ylab("Proportion of Explained Variance") +
    ggtitle("Scree Plot")
}

