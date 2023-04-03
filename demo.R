library(MASS)
# library(ltsspca)
# library(PMA)
library(irlba)
library(tidyverse)

### Functions ###
synthetic_tests <- function(v1, v2, eig, seed, rank=2, iters=50, n=30, error_tol=1e-2) {
  v1 <- normalize(v1, center=FALSE)
  v2 <- normalize(v2, center=FALSE)
  
  # result data structs
  cnames <- c('v1_angle', 
              'v2_angle', 
              'time',
              'v1_sparse_accuracy', 
              'v1_sparse_fnr',
              'v2_sparse_accuracy',
              'v2_sparse_fnr',
              'iter',
              'final_loss',
              'v1_SSE',
              'v2_SSE')
  df_amanpg <- data.frame(matrix(nrow=0, ncol=13))
  df_rspca <- data.frame(matrix(nrow=0, ncol=13))
  df_pca <- data.frame(matrix(nrow=0, ncol=13))
  df_pmd <- data.frame(matrix(nrow=0, ncol=13))
  colnames(df_amanpg) <- cnames
  colnames(df_rspca) <- cnames
  colnames(df_pmd) <- cnames
  colnames(df_pca) <- cnames
  
  v1_gt <- abs(v1)
  v2_gt <- abs(v2)
  v1_amanpg <- numeric(0)
  v2_amanpg <- numeric(0)
  v1_rsvd <- numeric(0)
  v2_rsvd <- numeric(0)
  
  p = length(eig)
  
  for (i in 1:iters) {
    print(i)
    v3 <- numeric(0)
    set.seed(seed+i)
    for (i in 1:(p-2)) {
      v3 <- c(v3, runif(p))
    }
    v3 <- normalize(v3, center=FALSE)
    vstar <- matrix(c(v1, v2, v3), nrow=p, ncol=p)
    v <- qr.Q(qr(vstar))
    c <- diag(eig)  # eigenvalues
    sigma <- v %*% c %*% t(v)
    
    # Generated data
    set.seed(1)
    x <- mvrnorm(n=n, mu=rep(0, p), Sigma=sigma)
    x <- normalize(x)

    # start = Sys.time()
    t1 <- system.time(sprout <- rspca.amanpg(x, 1, gamma=0.5, rank, verbose=FALSE, normalize=FALSE, maxiter=1000))
    # sprout <- rspca.amanpg(x, 1, 5, verbose=FALSE)
    # t1 = difftime(Sys.time(), start)
    # pmd <- SPC(x, niter=20, K=2, center=FALSE)
    # amanpg <- spca.amanpg(x, 1, 0.1, k=2)
    v1_pred <- abs(sprout$loadings[, 1])
    v2_pred <- abs(sprout$loadings[, 2])
    
    amanpg_error_v1 <- calc_pc_error(v1_pred, v1_gt, tol=error_tol)
    amanpg_error_v2 <- calc_pc_error(v2_pred, v2_gt, tol=error_tol)
    
    amanpg_row <- data.frame(
      v1_angle=acos(sum(v1_gt*v1_pred) / (sqrt(sum(v1_gt*v1_gt)) * sqrt(sum(v1_pred*v1_pred)))),  # v1 angle
      v2_angle=acos(sum(v2_gt*v2_pred) / (sqrt(sum(v2_gt*v2_gt)) * sqrt(sum(v2_pred*v2_pred)))),  # v2 angle
      time=t1[3],  # time
      v1_sparse_accuracy=amanpg_error_v1$sparse_accuracy,  #v1_sparse_accuracy
      v1_sparse_fnr=amanpg_error_v1$sparse_fnr,  # v1_sparse_fnr
      v2_sparse_accuracy=amanpg_error_v2$sparse_accuracy,  #v2_sparse_accuracy
      v2_sparse_fnr=amanpg_error_v2$sparse_fnr,  # v2_sparse_fnr
      iter=sprout$iterations,
      final_loss=sprout$loss[length(sprout$loss)] - sprout$loss[length(sprout$loss)-1],
      v1_SSE=sum((v1_gt-v1_pred)^2),
      v2_SSE=sum((v2_gt-v2_pred)^2)
    )
    
    print(sprout$iterations)
    
    v1_amanpg <- v1_pred
    v2_amanpg <- v2_pred
    
    t2 <- system.time(rsprout <- spca.rsvd(x, k=rank, tol=1e-5, n=40, lambda=0.25, maxit=1000))
    loadings <- rsprout$v
    # rsprout <- sPCA_rSVD(x, 2, center=T)
    # calculate angles for regularized method
    v1_pred <- abs(loadings[, 1])
    v2_pred <- abs(loadings[, 2])
    
    print(rsprout$lambda)
    
    rspca_error_v1 <- calc_pc_error(v1_pred, v1_gt, tol=error_tol)
    rspca_error_v2 <- calc_pc_error(v2_pred, v2_gt, tol=error_tol)
    
    rspca_row <- data.frame(
      v1_angle=acos(sum(v1_gt*v1_pred) / (sqrt(sum(v1_gt*v1_gt)) * sqrt(sum(v1_pred*v1_pred)))),  # v1 angle
      v2_angle=acos(sum(v2_gt*v2_pred) / (sqrt(sum(v2_gt*v2_gt)) * sqrt(sum(v2_pred*v2_pred)))),  # v2 angle
      time=t2[3],  # time
      v1_sparse_accuracy=rspca_error_v1$sparse_accuracy,  #v1_sparse_accuracy
      v1_sparse_fnr=rspca_error_v1$sparse_fnr,  # v1_sparse_fnr
      v2_sparse_accuracy=rspca_error_v2$sparse_accuracy,  #v2_sparse_accuracy
      v2_sparse_fnr=rspca_error_v2$sparse_fnr,  # v2_sparse_fnr
      iter=rsprout$iter,
      final_loss=tail(rsprout$loss, n=1),
      v1_SSE=sum((v1_gt-v1_pred)^2),
      v2_SSE=sum((v2_gt-v2_pred)^2)
    )
    
    v1_rsvd <- v1_pred
    v2_rsvd <- v2_pred
    
    t3 <- system.time(pr.out <- prcomp(x))
    loadings <- pr.out$rotation
    v1_pred <- abs(loadings[, 1])
    v2_pred <- abs(loadings[, 2])
    
    pca_error_v1 <- calc_pc_error(v1_pred, v1_gt, tol=error_tol)
    pca_error_v2 <- calc_pc_error(v2_pred, v2_gt, tol=error_tol)
    
    pca_row <- data.frame(
      v1_angle=acos(sum(v1_gt*v1_pred) / (sqrt(sum(v1_gt*v1_gt)) * sqrt(sum(v1_pred*v1_pred)))),  # v1 angle
      v2_angle=acos(sum(v2_gt*v2_pred) / (sqrt(sum(v2_gt*v2_gt)) * sqrt(sum(v2_pred*v2_pred)))),  # v2 angle
      time=t3[3],  # time
      v1_sparse_accuracy=pca_error_v1$sparse_accuracy,  #v1_sparse_accuracy
      v1_sparse_fnr=pca_error_v1$sparse_fnr,  # v1_sparse_fnr
      v2_sparse_accuracy=pca_error_v2$sparse_accuracy,  #v2_sparse_accuracy
      v2_sparse_fnr=pca_error_v2$sparse_fnr,  # v2_sparse_fnr
      iter=NA,
      final_loss=NA,
      v1_SSE=sum((v1_gt-v1_pred)^2),
      v2_SSE=sum((v2_gt-v2_pred)^2)
    )
    
    # t3 <- system.time(pmd <- pmd.spc(x, sumabsv=1, niter=2000, K=5, center=FALSE, compute.pve=FALSE))
    # 
    # v1_pred <- abs(pmd$v[, 1])
    # v2_pred <- abs(pmd$v[, 2])
    # 
    # pmd_error_v1 <- calc_pc_error(v1_pred, v1_gt, tol=error_tol)
    # pmd_error_v2 <- calc_pc_error(v2_pred, v2_gt, tol=error_tol)
    # 
    # pmd_row <- data.frame(
    #   v1_angle=acos(sum(v1_gt*v1_pred) / (sqrt(sum(v1_gt*v1_gt)) * sqrt(sum(v1_pred*v1_pred)))),  # v1 angle
    #   v2_angle=acos(sum(v2_gt*v2_pred) / (sqrt(sum(v2_gt*v2_gt)) * sqrt(sum(v2_pred*v2_pred)))),  # v2 angle
    #   time=t3[3],  # time
    #   v1_accuracy=pmd_error_v1$accuracy,  # v1_accuracy
    #   v1_sparse_accuracy=pmd_error_v1$sparse_accuracy,  #v1_sparse_accuracy
    #   v1_sparse_fnr=pmd_error_v1$sparse_fnr,  # v1_sparse_fnr
    #   v2_accuracy=pmd_error_v2$accuracy,  # v2_accuracy
    #   v2_sparse_accuracy=pmd_error_v2$sparse_accuracy,  #v2_sparse_accuracy
    #   v2_sparse_fnr=pmd_error_v2$sparse_fnr,  # v2_sparse_fnr
    #   iter=NA,
    #   final_loss=NA
    # )
    
    df_amanpg <- rbind(df_amanpg, amanpg_row)
    df_rspca <- rbind(df_rspca, rspca_row)
    df_pca <- rbind(df_pca, pca_row)
    # df_pmd <- rbind(df_pmd, pmd_row)
    # scree.plot(rsprout, x, normalized=FALSE)
  }
  
  pcs <- rbind(data.frame(type="GT", v1=v1_gt, v2=v2_gt),
               data.frame(type="Ours", v1=v1_amanpg, v2=v2_amanpg),
               data.frame(type="sPCA-rSVD", v1=v1_rsvd, v2=v2_rsvd))
  
  return(list(
    amanpg=df_amanpg,
    rspca=df_rspca,
    pca=df_pca,
    # pmd=df_pmd,
    pcs=pcs,
    amanpg_data=sprout,
    x=x
  ))
}

calc_pc_error <- function(pred, gt, tol=1e-2) {
  sparse_gt_idx <- which(gt == 0)
  sparse_pred_idx <- which(pred == 0)
  common <- intersect(sparse_gt_idx, sparse_pred_idx)
  
  sparse_accuracy <- length(common) / length(sparse_gt_idx)
  sparse_fnr <- max(0, length(sparse_pred_idx) - length(common))  / length(pred)
  
  return(as.data.frame(cbind(
    sparse_accuracy, 
    sparse_fnr
  )))
}

### Synthetic Data ###

# Covariance matrix using sparse & uniformly-generated eigenvectors
set.seed(5105)
v1 <- sample(c(1, 0.5), replace=TRUE, size=100)
v2 <- sample(c(1, 0.5), replace=TRUE, size=100)
v1_indices <- sample(1:100, 40)
v1[v1_indices] <- 0
v2_indices <- sample(1:100, 40)
v2[v2_indices] <- 0
# v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
# v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
eig <- c(1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, rep.int(1, 100-15))  # eigenvalues
# eig <- c(200, 100, 5, 6, 5, 4, 3, 4, 2, 2)

synth_test_results <- synthetic_tests(v1, v2, eig, 5105, rank=15, iters=5, n=1000)
synth_amanpg <- synth_test_results$amanpg
synth_rspca <- synth_test_results$rspca
# synth_pmd <- synth_test_results$pmd

synth_data <- synth_amanpg %>%
  mutate(type="Ours") %>%
  rbind(synth_rspca %>% mutate(type="sPCA-rSVD")) %>%
  rbind(synth_test_results$pca %>% mutate(type="PCA"))

time_plot <- ggplot(data=synth_data, aes(x=type, y=time)) +
  geom_boxplot() + 
  # ylim(0, 10) + 
  ggtitle("CPU Time versus Algorithm") +
  xlab("Sparse PCA Algorithm") +
  ylab("CPU Time (s)")
time_plot
# ggsave("~/Downloads/time_plot.png", plot=time_plot, width=10, height=4, dpi=300)

accuracies <- synth_data %>%
  pivot_longer(cols=c(v1_sparse_accuracy, v2_sparse_accuracy),
               names_to="pc",
               values_to="sparse_accuracy") %>%
  mutate(pc=ifelse(pc=="v1_sparse_accuracy", 1, 2))

acc_plot <- ggplot(data=accuracies, aes(x=type, y=sparse_accuracy, fill=pc)) +
  geom_boxplot() +
  ggtitle("Principal Component Sparse Accuracy by Algorithm") +
  xlab("Sparse PCA Algorithm") +
  ylab("Correct Sparse Loadings (%)") +
  labs(fill="PC")
acc_plot
ggsave("~/Downloads/acc_plot.png", plot=acc_plot, width=10, height=4, dpi=300)

pcs <- synth_test_results$pcs
  
ggplot(data=pcs[pcs$type != "GT",], aes(x=v1, y=v2)) +
  geom_point(aes(colour=type), pch=21, fill=NA, size=3, stroke=1, alpha=0.5) +
  geom_point(data=pcs[pcs$type == "GT",]) +
  ggtitle("Predicted Loadings versus Ground Truth") +
  xlab("PC1") +
  ylab("PC2") +
  labs(colour="Method")
  # annotate("text", x = 0.05, y = 0.11, label = "Black points are GT", hjust=1, size = 3)

synth_data %>%
  mutate(type = as.factor(type)) %>%
  group_by(type) %>%
  summarize_at(1:9, median)

resid_plot <- ggplot(data=synth_data, aes(x=v1_SSE+v2_SSE, fill=type)) +
  geom_boxplot() +
  xlab("Sum of squared error of residuals") +
  ggtitle("Error with 100x500 input matrix") +
  theme(axis.text.y=element_blank())
resid_plot
ggsave("~/Downloads/resid_plot_100x500.png", plot=resid_plot, width=10, height=4, dpi=300)
