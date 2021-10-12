rspca.amanpg <- function(z, l1, k, tol=1e-5, maxiter=1000, fmax=1e6, type=0, 
                         gamma=0.5, normalize=TRUE, verbose=TRUE) {
  start <- Sys.time()
  
  if (normalize) {
    z <- normalize(z)
  }
  
  dims <- dim(z)
  m <- dims[1]
  d <- dims[2]
  
  if (d < m * 2) {
    z <- Conj(t(z)) %*% z
    type <- 1
  }
  
  ### Initialized values ###
  svds <- svd(z, nu=k, nv=k)
  a <- normalize(svds$u, center=FALSE)
  b <- svds$d * svds$v
  
  # Get Lipschitz
  if (!type) {
    t <- 1 / (2 * svds$d[1]^2)  # TODO: Why is svd involved? can I use norm(x)^2?
  } else {
    t <- 1 / (2 * svds$d[1])
  }
  
  total_linesearch <- 0
  min_step <- 0
  alpha <- 100 / d
  
  ### Main Loop ###
  for (iter in 1:maxiter) {
    if (verbose) {
      iter_start <- Sys.time()
      print("=================")
      print(paste("On iteration", iter))
    }
    
    ### Update A ###
    # TODO: Solve for A subproblem
    da = 0
    
    ## Retract ##
    if (!linesearch_flag) alpha <- alpha * 1.1
    if (min_step == 1) alpha <- 1 / d
    linesearch_flag <- 0
    min_step <- 0
    
    xi <- alpha * da
    retracted_a <- (a + xi) %*% (diag(r) + t(xi) %*% xi) ^ (-1/2)
    
    while (sum(2*x*(a - retracted_a)) > -(a/(2*t)) * norm(da)^2) {
      alpha <- gamma * alpha
      if (alpha < 1e-5 / d) {
        min_step <- 1
        break
      }
      
      linesearch_flag <- 1
      xi <- alpha * da
      retracted_a <- (a + xi) %*% (diag(r) + t(xi) %*% xi) ^ (-1/2)
    }
    
    a <- retracted_a
    
    ## Update B ##
    c <- b + 2 * t * t(z - a %*% t(b)) %*% a
    db <- c - b + t * sign(c) * matrix(as.numeric(abs(mat) > l1), nrow=d, ncol=k)  # done element-wise
    b <- db
    
    check <- norm(da/t)^2 + norm(db/t)^2
    
    if (verbose) {
      print(paste("da:", da))
      print(paste("db:", db))
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

normalize <- function(x, center=TRUE, scale=TRUE) {
  if (center) {
    x <- scale(x, scale=FALSE)
  }
  
  if (scale) {
    dims <- dim(x)
    x <- x / matrix(sqrt(rowSums(x^2)), nrow=dims[1], ncol=dims[2])
  }
  
  return(x)
}