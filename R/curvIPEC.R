curvIPEC <-
function( expr, theta, x, y, tol = 1e-16, alpha=0.05, method="Richardson",                            
              method.args=list(eps=1e-4, d=0.11, zero.tol=sqrt(.Machine$double.eps/7e-7), 
              r=6, v=2, show.details=FALSE), side=NULL ){    

  x     <- rbind( x )
  y     <- as.vector(y)
  if( min(dim(x))[1] ==1 )    x <- cbind( x )
  if( nrow(x) != length(y) )  x <- t(x)
  n     <- nrow(x)
  p     <- length(theta)
  yhat  <- expr(theta, x)
  sigma <- sqrt( sum((y-yhat)^2)/(n-p) )
  sp    <- sigma * sqrt(p)
  cc    <- 1/( sqrt(qf(1-alpha, p, n-p)) )
  v1    <- matrix(NA, nrow=n, ncol=p)
  v2    <- array(NA, c(p, p, n))
  for(k in 1L:n){
    resu      <- derivIPEC(expr, theta, x[k, ], method=method, 
                     method.args=method.args, side=side)
    v1[k,]    <- resu$Jacobian
    v2[, , k] <- resu$Hessian
  }  
  a <- array( NA, c(n, p, p) )
  for(k in 1L:p){
    temp.mat  <- matrix(NA, nrow=1, ncol=p) 
    for(q in 1L:n) temp.mat <- rbind(temp.mat, v2[k, , q])
    temp.mat  <- temp.mat[-1, ]
    a[ , , k] <- temp.mat
  }
  D <- v1
  for (j in 1L:p){
    D <- cbind(D, a[, 1L:p, j])
  }
  qrd <- qr(D, tol)
  Q   <- qr.Q(qrd)
  rnk <- qrd$rank
  if (rnk <= p){ 
    warning("regression apparently linear")
  }
  Q1 <- Q[, 1L:rnk]
  C  <- array(0, c(rnk, p, p))
  for (j in 1L:p){
    C[, , j] <- crossprod(Q1, a[, , j])
  }
  C    <- aperm(C, c(2, 3, 1))
  r11i <- ginv( qr.R(qrd)[1L:p, 1L:p], tol=tol )
  ct   <- 0
  for (j in 1L:p) {
    C[, , j] <- crossprod(r11i, C[, , j]) %*% r11i * sp
    ct       <- ct + 2 * sum(C[, , j]^2) + sum(diag(C[, , j]))^2
  }
  ci <- 0
  for (j in (p + 1):rnk) {
    C[, , j] <- crossprod(r11i, C[, , j]) %*% r11i * sp
    ci <- ci + 2 * sum(C[, , j]^2) + sum(diag(C[, , j]))^2
  }
  ct <- sqrt(ct/(p * (p + 2)))
  ci <- sqrt(ci/(p * (p + 2)))
  # pe <- ct / cc
  # ic <- ci / cc
  list(rms.ic = ci, rms.pec = ct, critical.c = cc)
}






