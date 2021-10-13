library(MASS)

v1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
v2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)

v1 <- normalize(v1, center=FALSE)
v2 <- normalize(v2, center=FALSE)

v3 <- numeric(0)
for (i in 1:8) {
  set.seed(1)
  v3 <- c(v3, runif(10))
}
v3 <- normalize(v3, center=FALSE)
vstar <- matrix(c(v1, v2, v3), nrow=10, ncol=10)
v <- qr.Q(qr(vstar))
c <- diag(c(200,100,50,50,6,5,4,3,2,1))
sigma <- v %*% c %*% t(v)

set.seed(1)
x <- mvrnorm(n=30, mu=rep(0, 10), Sigma=sigma)

rsprout <- rspca.amanpg(x, 0.1, 2, verbose=TRUE)
