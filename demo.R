library(MASS)
library(ltsspca)
library(tidyverse)

### Functions ###
synthetic_tests <- function(v1, v2, eig, seed, iters=50, n=30, error_tol=1e-2) {
  v1 <- normalize(v1, center=FALSE)
  v2 <- normalize(v2, center=FALSE)
  
  # result data structs
  cnames <- c('v1_angle', 
              'v2_angle', 
              'time', 
              'v1_accuracy', 
              'v1_sparse_accuracy', 
              'v1_sparse_fnr',
              'v2_accuracy',
              'v2_sparse_accuracy',
              'v2_sparse_fnr')
  df_amanpg <- data.frame(matrix(nrow=0, ncol=9))
  df_rspca <- data.frame(matrix(nrow=0, ncol=9))
  colnames(df_amanpg) <- cnames
  colnames(df_rspca) <- cnames
  
  v1_gt <- abs(v1)
  v2_gt <- abs(v2)
  
  for (i in 1:iters) {
    v3 <- numeric(0)
    set.seed(seed+i)
    for (i in 1:8) {
      v3 <- c(v3, runif(10))
    }
    v3 <- normalize(v3, center=FALSE)
    vstar <- matrix(c(v1, v2, v3), nrow=10, ncol=10)
    v <- qr.Q(qr(vstar))
    c <- diag(eig)  # eigenvalues
    sigma <- v %*% c %*% t(v)
    
    # Generated data
    set.seed(1)
    x <- mvrnorm(n=30, mu=rep(0, 10), Sigma=sigma)
    
    sprout <- rspca.amanpg(x, 0.455, 2, verbose=FALSE)
    # amanpg <- spca.amanpg(x, 1, 0.1, k=2)
    v1_pred <- abs(sprout$loadings[, 1])
    v2_pred <- abs(sprout$loadings[, 2])
    
    amanpg_error_v1 <- calc_pc_error(v1_pred, v1_gt, tol=error_tol)
    amanpg_error_v2 <- calc_pc_error(v2_pred, v2_gt, tol=error_tol)
    
    amanpg_row <- data.frame(
      v1_angle=acos(sum(v1_gt*v1_pred) / (sqrt(sum(v1_gt*v1_gt)) * sqrt(sum(v1_pred*v1_pred)))),  # v1 angle
      v2_angle=acos(sum(v2_gt*v2_pred) / (sqrt(sum(v2_gt*v2_gt)) * sqrt(sum(v2_pred*v2_pred)))),  # v2 angle
      time=sprout$time,  # time
      v1_accuracy=amanpg_error_v1$accuracy,  # v1_accuracy
      v1_sparse_accuracy=amanpg_error_v1$sparse_accuracy,  #v1_sparse_accuracy
      v1_sparse_fnr=amanpg_error_v1$sparse_fnr,  # v1_sparse_fnr
      v2_accuracy=amanpg_error_v2$accuracy,  # v2_accuracy
      v2_sparse_accuracy=amanpg_error_v2$sparse_accuracy,  #v2_sparse_accuracy
      v2_sparse_fnr=amanpg_error_v2$sparse_fnr  # v2_sparse_fnr
    )
    
    rsprout_start <- Sys.time()
    rsprout <- sPCA_rSVD(x, 2, center=T)
    # store runtime of regularized method
    time_rspca <- difftime(Sys.time(), rsprout_start)
    # calculate angles for regularized method
    v1_pred <- abs(rsprout$loadings[, 1])
    v2_pred <- abs(rsprout$loadings[, 2])
    
    rspca_error_v1 <- calc_pc_error(v1_pred, v1_gt, tol=error_tol)
    rspca_error_v2 <- calc_pc_error(v2_pred, v2_gt, tol=error_tol)
    
    rspca_row <- data.frame(
      v1_angle=acos(sum(v1_gt*v1_pred) / (sqrt(sum(v1_gt*v1_gt)) * sqrt(sum(v1_pred*v1_pred)))),  # v1 angle
      v2_angle=acos(sum(v2_gt*v2_pred) / (sqrt(sum(v2_gt*v2_gt)) * sqrt(sum(v2_pred*v2_pred)))),  # v2 angle
      time=time_rspca,  # time
      v1_accuracy=rspca_error_v1$accuracy,  # v1_accuracy
      v1_sparse_accuracy=rspca_error_v1$sparse_accuracy,  #v1_sparse_accuracy
      v1_sparse_fnr=rspca_error_v1$sparse_fnr,  # v1_sparse_fnr
      v2_accuracy=rspca_error_v2$accuracy,  # v2_accuracy
      v2_sparse_accuracy=rspca_error_v2$sparse_accuracy,  #v2_sparse_accuracy
      v2_sparse_fnr=rspca_error_v2$sparse_fnr  # v2_sparse_fnr
    )
    
    df_amanpg <- rbind(df_amanpg, amanpg_row)
    df_rspca <- rbind(df_rspca, rspca_row)
    # scree.plot(rsprout, x, normalized=FALSE)
  }
  
  return(list(
    df_amanpg,
    df_rspca
  ))
}

calc_pc_error <- function(pred, gt, tol=1e-3) {
  sparse_gt_idx <- which(gt == 0)
  sparse_pred_idx <- which(pred == 0)
  common <- intersect(sparse_gt_idx, sparse_pred_idx)
  
  sparse_accuracy <- length(common) / length(sparse_gt_idx)
  sparse_fnr <- max(0, length(sparse_pred_idx) - length(common))  / length(pred)
  
  diff <- ifelse(abs(pred - gt) < tol, 1, 0)
  accuracy <- sum(diff) / length(gt)
  
  return(as.data.frame(cbind(
    accuracy, 
    sparse_accuracy, 
    sparse_fnr
  )))
}

### Synthetic Data ###

# Covariance matrix using sparse & uniformly-generated eigenvectors
v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
eig <- c(200,100,5,5,6,5,4,3,2,1)  # eigenvalues

synth_test_results <- synthetic_tests(v1, v2, eig, 1, iters=10)
synth_amanpg <- synth_test_results[[1]]
synth_rspca <- synth_test_results[[2]]

# TODO: more generated data results, check how others measure performance
# TODO: check how scores are restored to nxk
# TODO: optimize code

### UCI Breast Cancer Dataset ###

# cancer <- read.table("wdbc.data", sep=",") %>% select(-c(1, 2))  # drop ID and diagnosis
# cancer <- normalize(cancer)
# system.time(sprout2 <- rspca.amanpg(as.matrix(cancer), 1, 2, maxiter=2000, tol=1e-4))
# amanpg2 <- spca.amanpg(as.matrix(cancer), 1, 0.1, k=2, verbose=TRUE, maxiter=2000, tol=1e-4)
# system.time(rsprout2 <- sPCA_rSVD(x, 2, center=T))

