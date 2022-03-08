library(MASS)
library(ltsspca)
library(tidyverse)

### Generated Data

# Covariance matrix using sparse & uniformly-generated eigenvectors
v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)

v1 <- normalize(v1, center=FALSE)
v2 <- normalize(v2, center=FALSE)

v3 <- numeric(0)
set.seed(1)
for (i in 1:8) {
  v3 <- c(v3, runif(10))
}
v3 <- normalize(v3, center=FALSE)
vstar <- matrix(c(v1, v2, v3), nrow=10, ncol=10)
v <- qr.Q(qr(vstar))
c <- diag(c(200,100,5,5,6,5,4,3,2,1))  # eigenvalues
sigma <- v %*% c %*% t(v)

# Generated data
set.seed(1)
x <- mvrnorm(n=30, mu=rep(0, 10), Sigma=sigma)

sprout <- rspca.amanpg(x, 1, 2)
# amanpg <- spca.amanpg(x, 1, 0.1, k=2)
rsprout <- sPCA_rSVD(x, 2, center=T)

### UCI Breast Cancer Dataset

cancer <- read.table("wdbc.data", sep=",") %>% select(-c(1, 2))  # drop ID and diagnosis
cancer <- normalize(cancer)
sprout2 <- rspca.amanpg(as.matrix(cancer), 1, 2, verbose=TRUE, maxiter=2000, tol=1e-4)
amanpg2 <- spca.amanpg(as.matrix(cancer), 1, 0.1, k=2, verbose=TRUE, maxiter=2000, tol=1e-4)

# TODO: more generated data results, check how others measure performance
# TODO: check how scores are restored to nxk
# TODO: optimize code



