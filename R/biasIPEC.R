biasIPEC <-
function( expr, theta, x, y, tol = .Machine$double.eps, method="Richardson",
                  method.args=list(eps=1e-4, d=0.11, zero.tol=sqrt(.Machine$double.eps/7e-7), 
                  r=6, v=2, show.details=FALSE), side=NULL ){

  x     <- rbind( x )
  y     <- as.vector(y)
  if( min(dim(x))[1] ==1 )    x <- cbind( x )
  if( nrow(x) != length(y) )  x <- t(x)
  n     <- nrow(x)
  p     <- length(theta)
  yhat  <-  expr(theta, x)
  sigma <- sqrt( sum((y-yhat)^2)/(n-p) )
  temp1 <- 0
  for(k in 1:n){
    resu  <- derivIPEC(expr, theta, x[k, ], method=method, method.args=method.args, side=side)
    Fu    <- resu$Jacobian
    temp1 <- temp1 + Fu %*% t(Fu)    
  }
  temp1 <- solve(temp1, tol=tol)
  temp2 <- 0
  for(k in 1:n){    
    resu  <- derivIPEC(expr, theta, x[k, ], method=method, method.args=method.args, side=side)    
    Fu    <- resu$Jacobian
    Ht    <- resu$Hessian 
    temp2 <- temp2 + Fu * sum(diag( temp1 %*% Ht ))    
  }  
  bias         <- array( -sigma^2/2*temp1 %*% temp2 )
  percent.bias <- bias/theta * 100
  list(bias = bias, percent.bias = paste(round(percent.bias,5), "%", sep=""))  
}
