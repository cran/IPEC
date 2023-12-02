skewIPEC <- function( expr, theta, x, y, tol=sqrt(.Machine$double.eps), method = "Richardson", 
        method.args = list(eps = 1e-04, d = 0.11, zero.tol = sqrt(.Machine$double.eps/7e-07), 
        r = 6, v = 2, show.details = FALSE), side = NULL ) 
{
    x     <- rbind( x )
    y     <- as.vector(y)
    if( min(dim(x))[1] ==1 )    x <- cbind( x )
    if( nrow(x) != length(y) )  x <- t(x)
    n     <- nrow(x)
    p     <- length(theta)
    yhat  <- expr(theta, x)
    sigma <- sqrt(sum((y - yhat)^2)/(n - p))
    v1    <- matrix(NA, nrow = n, ncol = p)
    v2    <- array(NA, c(p, p, n))
    for (k in 1L:n) {
        resu <- derivIPEC(expr, theta, x[k, ], method = method, 
            method.args = method.args, side = side)
        v1[k, ]   <- resu$Jacobian
        v2[, , k] <- resu$Hessian
    }
  
    # L is a matrix with p rows and p columns
    L <- t(v1) %*% v1   
    L <- ginv(L, tol=tol)

    W <- array(NA, c(p, p, p))
    for (k1 in 1L:p) {
      for (k2 in 1L:p){
        for (k3 in 1L:p){
          temp.sum <- 0
          for(k4 in 1L:n){           
            temp.val <- sum( v1[k4, k1] * v2[k2, k3, k4] )   
            temp.sum <- temp.sum + temp.val         
          }
          W[k2, k3, k1] <- temp.sum
        }
      }
    }

    numer <- c()
    for (k4 in 1:p){
      mysum <- 0
      for (k1 in 1:p){
        for (k2 in 1:p){
          for (k3 in 1:p){
            mysum <- mysum + L[k4, k1] * L[k4, k2] * L[k4, k3] * (
                       W[k2,k3,k1] + W[k1,k3,k2] + W[k2,k1,k3])
          }     
        }    
      } 
      numer[k4] <- -sigma^4 * mysum    
    }

    denomi <- ( sigma^2 * diag(L) )^(3/2)
    g      <- numer / denomi
    list(skewness = g)
}
