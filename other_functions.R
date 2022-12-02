library(irlba)

spca.rsvd <- function (x, k = 1, n = 2, maxit = 500, tol = 0.001, lambda = 1, center = FALSE, 
          scale. = FALSE, alpha = 0, tsvd = NULL, ...) 
{
  if (alpha < 0 || alpha >= 1) 
    stop("0 <= alpha < 1")
  if (is.logical(center) && center) 
    center <- colMeans(x)
  if (is.logical(scale.)) {
    if (scale.) {
      if (is.numeric(center)) {
        f <- function(i) sqrt(sum((x[, i] - center[i])^2)/(nrow(x) - 
                                                             1L))
        scale. <- vapply(seq(ncol(x)), f, pi, USE.NAMES = FALSE)
      }
      else scale. <- apply(x, 2L, function(v) sqrt(sum(v^2)/max(1, 
                                                                length(v) - 1L)))
    }
  }
  if (all(n > ncol(x) - 1)) {
    warning("no sparsity constraints specified")
    return(irlba(x, k, ...))
  }
  n <- ncol(x) - n
  if (length(n) != k) 
    n <- rep(n, length.out = k)
  s <- tsvd
  if (is.null(tsvd)) 
    s <- irlba(x, k, scale = scale., center = center, ...)
  lambda <- c()
  soft <- function(x, u, p) {
    y <- crossprod(x, u)
    if (is.numeric(center)) 
      y <- y - sum(u) * center
    if (is.numeric(scale.)) 
      y <- y/scale.
    a <- abs(y)
    z <- apply(a, 2, sort)
    lambda <<- vapply(seq(length(p)), function(j) (1 - alpha) * 
                        z[p[j], j] + alpha * z[p[j] + 1, j], pi, USE.NAMES = FALSE)
    sign(y) * pmax(sweep(a, 2, lambda, `-`), 0)
    # sign(y) * pmax(y - lambda, 0)
  }
  s$v <- s$d * s$v
  iter <- 0
  delta_u <- Inf
  loss <- c()
  
  while (delta_u > tol && iter < maxit) {
    u <- s$u
    s$v <- soft(x, s$u, n)
    if (is.numeric(scale.)) 
      s$v <- s$v/scale.
    if (is.numeric(center)) {
      xsv <- x %*% s$v - drop(crossprod(center, s$v))
      s$u <- qr.Q(qr(xsv))
    }
    else {
      xsv <- x %*% s$v
      s$u <- qr.Q(qr(xsv))
    }
    s$u <- sweep(s$u, 2, apply(xsv, 2, function(x) sign(head(x[x != 
                                                                 0], 1)))/apply(s$u, 2, function(x) sign(head(x[x != 
                                                                                                                  0], 1))), `*`)
    delta_u <- max(1 - diag(abs(crossprod(u, s$u))))
    iter <- iter + 1
    loss <- c(loss, delta_u)
  }
  if (iter >= maxit) 
    warning("Maximum number of iterations reached before convergence: solution may not be optimal. Consider increasing 'maxit'.")
  s$v <- s$v %*% diag(1/sqrt(apply(s$v, 2, crossprod)), ncol(s$v), 
                      ncol(s$v))
  d <- s$v
  if (is.numeric(scale.)) 
    d <- d/scale.
  d1 <- x %*% d
  if (is.numeric(center)) 
    d1 <- d1 - drop(crossprod(center, d))
  d <- crossprod(s$u, d1)
  list(u = s$u, v = s$v, d = d, iter = iter, lambda = lambda, 
       center = center, scale = scale., n = n, alpha = alpha, 
       loss=loss)
}

pmd.spc <- function (x, sumabsv = 4, niter = 20, K = 1, orth = FALSE, trace = TRUE, 
                     v = NULL, center = TRUE, cnames = NULL, vpos = FALSE, vneg = FALSE, 
                     compute.pve = TRUE) 
{
  if (vpos && vneg) 
    stop("Cannot constrain elements to be positive AND negative.")
  out <- PMDL1L1(x, sumabsu = sqrt(nrow(x)), sumabsv = sumabsv, 
                 niter = niter, K = K, orth = orth, trace = trace, v = v, 
                 center = center, cnames = cnames, upos = FALSE, uneg = FALSE, 
                 vpos = vpos, vneg = vneg)
  if (compute.pve) {
    v <- matrix(out$v, ncol = K)
    ve <- NULL
    xfill <- x
    # if (center) 
    #   xfill <- x - mean.na(x)
    # xfill[is.na(x)] <- mean.na(xfill)
    for (k in 1:K) {
      vk <- matrix(v[, 1:k], ncol = k)
      xk <- xfill %*% vk %*% solve(t(vk) %*% vk) %*% t(vk)
      svdxk <- svd(xk)
      ve <- c(ve, sum(svdxk$d^2))
    }
    pve <- ve/sum(svd(xfill)$d^2)
    out$prop.var.explained <- pve
  }
  out$vpos <- vpos
  out$vneg <- vneg
  class(out) <- "SPC"
  return(out)
}


PMDL1L1 <- function(x,sumabs=.4,sumabsu=NULL,sumabsv=NULL,niter=20,K=1,v=NULL, trace=TRUE, orth=FALSE, center=TRUE, rnames=NULL, cnames=NULL, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg){
  if(center){
    # meanx <- mean.na(x)
    # x <- x-meanx
  } else {
    meanx <- NULL
  }
  if(orth){
    if(is.null(sumabsu) || sumabsu < sqrt(nrow(x))){
      orth <- FALSE
      warning("Orth option ignored because sparse PCA results only when sumabsu equals sqrt(nrow(x))")
    }
  }
  #  if((is.null(sumabsu) && !is.null(sumabsv)) || (is.null(sumabsv) && !is.null(sumabsu))) warning("Sumabsu and sumabsv BOTH must be input when type=standard. Since only one was given, it was ignored and sumabs was used instead.")
  if(is.null(sumabsu) || is.null(sumabsv)){
    sumabsu <- sqrt(nrow(x))*sumabs
    sumabsv <- sqrt(ncol(x))*sumabs
  }
  call <-  match.call()
  # if(trace && abs(mean.na(x)) > 1e-15) warning("PMDL1L1 was run without first subtracting out the mean of x.")
  if(!is.null(sumabsu) && (sumabsu<1 || sumabsu>sqrt(nrow(x)))) stop("sumabsu must be between 1 and sqrt(n)")
  if(!is.null(sumabsv) && (sumabsv<1 || sumabsv>sqrt(ncol(x)))) stop("sumabsv must be between 1 and sqrt(p)")
  v <- CheckPMDV(v,x,K)
  if(K>1 && !orth) out <- (MultiSMD(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,K=K, trace=trace, v=v, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg))
  if(K>1 && orth) out <- MultiSMDOrth(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,K=K, trace=trace, v=v,  vpos=vpos, vneg=vneg)
  if(K==1) out <- SMD(x,sumabsu=sumabsu,sumabsv=sumabsv,niter=niter, trace=trace, v=v, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
  obj <- (list(u=out$u,v=out$v, d=out$d, v.init=out$v.init, call=call, meanx=meanx,sumabsu=sumabsu, sumabsv=sumabsv, rnames=rnames, cnames=cnames, K=K, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg))
  class(obj) <- "PMDL1L1"
  return(obj)
}

CheckPMDV <-  function(v,x,K){
  if(!is.null(v) && is.matrix(v) && ncol(v)>=K){
    v <- matrix(v[,1:K], ncol=K)
  } else if(ncol(x)>nrow(x)){
    # x[is.na(x)] <- mean.na(x)
    v <- matrix(t(x)%*%(safesvd(x%*%t(x))$v[,1:K]),ncol=K)
    if(sum(is.na(v))>0) v <- matrix(safesvd(x)$v[,1:K], ncol=K)
    v <- sweep(v,2,apply(v, 2, l2n), "/")
    if(sum(is.na(v))>0) stop("some are NA")
  } else if (ncol(x)<=nrow(x)){
    # x[is.na(x)] <- mean.na(x)
    v <- matrix(safesvd(t(x)%*%x)$v[,1:K],ncol=K)
  }
  return(v)
}

MultiSMD <- function(x, sumabsu, sumabsv, K=3, niter=20,v, trace=TRUE, upos, uneg, vpos, vneg){
  nas <- is.na(x)
  v.init <- v
  xuse <- x
  ds <- numeric(K)
  us <- matrix(0,nrow=nrow(x),ncol=K)
  vs <- matrix(0,nrow=ncol(x),ncol=K)
  for(k in 1:K){
    out <- SMD(xuse, sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,v=matrix(v[,k],ncol=1), trace=trace, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
    us[,k] <- out$u
    vs[,k] <- out$v
    ds[k] <- out$d
    res <- xuse - out$d*out$u%*%t(out$v)
    xuse[!nas] <- res[!nas] # [!nas] is new on July 24 2009
  }
  return(list(u=us,v=vs,d=ds, v.init=v.init))
}

SMD <- function(x, sumabsu, sumabsv, niter=20,trace=TRUE, tol=1e-5, v, upos, uneg, vpos, vneg){
  # This gets a single factor. Do MultiSMD to get multiple factors.
  nas <- is.na(x)
  v.init <- v
  xoo <- x
  if(sum(nas)>0) xoo[nas] <- mean(x[!nas])
  oldv <- rnorm(ncol(x))
  delta_u <- Inf
  iter <- 0
  loss <- c()
  
  for(i in 1:niter) {
    if(sum(abs(oldv-v))>tol){
      oldv <- v
      if(trace) cat(iter,fill=F)
      # update u #
      argu <- xoo%*%v
      if(upos) argu <- pmax(argu,0)
      if(uneg) argu <- pmin(argu,0)
      lamu <- BinarySearch(argu,sumabsu)
      su <- soft(argu,lamu)
      u <- matrix(su/l2n(su),ncol=1)
      # done updating u #
      # update v #
      argv <- t(u)%*%xoo
      if(vpos) argv <- pmax(argv,0)
      if(vneg) argv <- pmin(argv,0)
      lamv <- BinarySearch(argv, sumabsv)
      sv <- soft(argv,lamv)
      v <- matrix(sv/l2n(sv),ncol=1)
      # done updating v #
    }
  }
  d <- as.numeric(t(u)%*%(xoo%*%v))
  if(trace) cat(fill=TRUE)
  return(list(d=d, u=u, v=v, v.init=v.init, iter=iter, loss=loss))
}

soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}

mean.na <- function(vec){
  return(mean(vec[!is.na(vec)]))
}

safesvd <- function(x){
  i <- 1
  out <- try(svd(x), silent=TRUE)
  while(i<10 && class(out)=="try-error"){
    out <- try(svd(x), silent=TRUE)
    i <- i+1
  }
  if(class(out)=="try-error") out <- svd(matrix(rnorm(nrow(x)*ncol(x)), ncol=ncol(x)))
  return(out)
}

BinarySearch <- function(argu,sumabs){
  if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter < 150){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    if((lam2-lam1)<1e-6) return((lam1+lam2)/2)
    iter <- iter+1
  }
  warning("Didn't quite converge")
  return((lam1+lam2)/2)
}

l2n <- function(vec){
  a <- sqrt(sum(vec^2))
  if(a==0) a <- .05
  return(a)
}
