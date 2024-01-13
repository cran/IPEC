parinfo <- function(object, x, CI = 0.95, method = "Richardson", 
             method.args = list(eps = 1e-04, d = 0.11, 
             zero.tol = sqrt(.Machine$double.eps/7e-07), r = 6, 
             v = 2, show.details = FALSE), side = NULL){
    if( ("sample.size" %in% names(object)) & ("n" %in% names(object)) )
        n <- object$sample.size
    if( ("sample.size" %in% names(object)) & !("n" %in% names(object)) )
        n <- object$sample.size
    if( !("sample.size" %in% names(object)) & ("n" %in% names(object)) )
        n <- object$n
    expr <- object$expr
    RSS  <- object$RSS
    par  <- object$par
    p    <- length( object$par )
    D    <- matrix(NA, nrow=n, ncol=p)
    x    <- rbind(x)
    if(nrow(x)==1){
      x <- matrix(x, ncol=1)
    }
    for(j in 1:nrow(x)){
      temp <- derivIPEC(expr, theta=par, z=x[j, ], method=method,
                        method.args=method.args, side=side)
      D[j,] <- temp$Jacobian
    }
    par.cov <- RSS/(n-p)*ginv( t(D) %*% D )
    par.se  <- c()
    for(j in 1:p){
      par.se <- c(par.se, sqrt(par.cov[j,j]))
    }
    alpha <- 1-CI
    tc    <- qt(1-alpha/2, df=n-p)
    lci   <- par - tc * par.se
    uci   <- par + tc * par.se
    par.tab <- data.frame(Estimate=par, SE=par.se, LCI=lci, UCI=uci)
    list(D=D, partab=par.tab, covmat=par.cov)
}

