derivIPEC <-
function( expr, theta, z, method="Richardson", 
              method.args=list(eps=1e-4, d=0.11, zero.tol=sqrt(.Machine$double.eps/7e-7), 
              r=6, v=2, show.details=FALSE), side=NULL ){

  z <- rbind( z )
  if(nrow(z) !=  1) stop("'z' should be a vector") 
  myfun <- function(P){
      expr(P, z)
  }
  Jacobian <- jacobian(myfun, x=theta, method=method, side=side, method.args=method.args)
  Jacobian <- array( Jacobian )
  Hessian  <- hessian(func=myfun, x=theta, method=method, method.args=method.args)
  list(Jacobian=Jacobian, Hessian = Hessian)

}
