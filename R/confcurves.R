confcurves <- function(expr, x, y, ini.val, weights = NULL, control=list(), 
               fig.opt = TRUE, fold = 5, np = 20, alpha = seq(1, 0.001, by=-0.001), 
               show.CI = NULL, 
               method = "Richardson", method.args = list(eps = 1e-04, d = 0.11, 
               zero.tol = sqrt(.Machine$double.eps/7e-07), r = 6, 
               v = 2, show.details = FALSE), side = NULL){ 

  if(!is.null(weights) && !is.numeric(weights)) 
      stop("'weights' must be a numeric vector!")
  if(!is.null(weights) && length(weights) != length(y))
      stop("The length of 'weights' must be the same as that of 'y'!")

  res1 <- fitIPEC(expr=expr, x=x, y=y, ini.val=ini.val, 
                  weights=weights, control=control, fig.opt=FALSE)
  n    <- res1$n
  p    <- length(res1$par)
  df   <- n-p
  tc   <- qt(1-alpha/2, df=df)
  RSS  <- res1$RSS
  MSE  <- RSS/df   
  res2 <- parinfo(res1, x=x, CI = 0.95, method=method, 
                  method.args=method.args, side=side)
  theta.hat <- res1$par
  theta.SE  <- res2$partab[,2]  
  parlist   <- list() 
  for(j in 1:p){
    Pj.hat <- theta.hat[j]
    Pj.SE  <- theta.SE[j]
    wald.lci.yj <- Pj.hat - tc * Pj.SE
    wald.uci.yj <- Pj.hat + tc * Pj.SE
    WaldCI      <- cbind(tc, wald.lci.yj, wald.uci.yj)
    colnames(WaldCI) <- c("tc", "LCI", "UCI")    
    ylim0 <- c(min(c(wald.lci.yj, wald.uci.yj)), max(c(wald.lci.yj, wald.uci.yj)))  
    delta <- fold * Pj.SE / np
    loop1.Pj.hat <- Pj.hat
    loop2.Pj.hat <- Pj.hat
    RSS.lower    <- c()
    RSS.upper    <- c()
    y.lower      <- c()
    y.upper      <- c()   
    for(q in 1:np){
      if(p==1){ 
        expr.lower0 <- function(x){   
          expr(loop1.Pj.hat-delta, x)  
        }
        expr.upper0 <- function(theta, x){      
          expr(loop2.Pj.hat+delta, x)  
        }
      }
      if(j==1 & p > 1){
        expr.lower <- function(theta, x){      
         expr(c(loop1.Pj.hat-delta, theta), x)  
        }
        expr.upper <- function(theta, x){      
          expr(c(loop2.Pj.hat+delta, theta), x)  
        }
      }
      if(j==p & p > 1){
        expr.lower <- function(theta, x){      
          expr(c(theta, loop1.Pj.hat-delta), x)  
        }
        expr.upper <- function(theta, x){      
          expr(c(theta, loop2.Pj.hat+delta), x)  
        }
      }
      if(j!=1 & j!=p & p > 1){
        expr.lower <- function(theta, x){  
          theta0 <- rep(NA, len=p)
          for(k in 1:p){
            if(k < j)
              theta0[k] <- theta[k]
            if(k == j)
              theta0[k] <- loop1.Pj.hat-delta
            if(k > j)
              theta0[k] <- theta[k-1] 
          }
          expr(theta0, x)  
        }
       expr.upper <- function(theta, x){ 
          theta0 <- rep(NA, len=p)
          for(k in 1:p){
            if(k < j)
              theta0[k] <- theta[k]
            if(k == j)
              theta0[k] <- loop2.Pj.hat+delta
            if(k > j)
              theta0[k] <- theta[k-1] 
          }
          expr(theta0, x)  
        }
      }
      if(p > 1){
        ini.val2  <- ini.val[-j]
        res.lower <- NULL
        res.upper <- NULL        
        try( res.lower <- fitIPEC(expr.lower, x=x, y=y, ini.val=ini.val2, 
             weights=weights, control=control, fig.opt=FALSE), silent=TRUE )
        try( res.upper <- fitIPEC(expr.upper, x=x, y=y, ini.val=ini.val2, 
             weights=weights, control=control, fig.opt=FALSE), silent=TRUE )
        if(!is.null(res.lower))
           RSS1 <- res.lower$RSS
        if(!is.null(res.upper))
           RSS2 <- res.upper$RSS
        if(is.null(res.lower))
           RSS1 <- NA
        if(is.null(res.upper))
           RSS2 <- NA
       }
      if(p == 1){
        RSS1 <- expr.lower0(x)
        RSS2 <- expr.upper0(x)
      }
      RSS.lower <- c(RSS.lower, RSS1)
      RSS.upper <- c(RSS.upper, RSS2)
      loop1.Pj.hat <- loop1.Pj.hat - delta
      loop2.Pj.hat <- loop2.Pj.hat + delta
      y.lower <- c(y.lower, loop1.Pj.hat)
      y.upper <- c(y.upper, loop2.Pj.hat)
    }
    temp1 <- cbind(RSS.lower, RSS.upper)
    temp1 <- na.omit(temp1)
    temp2 <- as.numeric(attr(temp1,"na.action"))
    RSS.lower <- temp1[,1]
    RSS.upper <- temp1[,2]
    if(length(temp2) > 0){
      y.lower <- y.lower[-temp2]
      y.upper <- y.upper[-temp2]
    } 
    x.lower  <- sqrt( abs(RSS.lower-RSS)/MSE ) 
    x.upper  <- sqrt( abs(RSS.upper-RSS)/MSE ) 
    x2.lower <- c(0, x.lower) 
    x2.upper <- c(0, x.upper)
    y2.lower <- c(Pj.hat, y.lower)
    y2.upper <- c(Pj.hat, y.upper)
    RSS2.lower <- c(RSS, RSS.lower)
    RSS2.upper <- c(RSS, RSS.upper)
    lhCI       <- cbind(x2.lower, y2.lower, x2.upper, y2.upper, RSS2.lower, RSS2.upper)
    
    colnames(lhCI) <- c("x.lower", "lhLCI", "x.upper", "lhUCI", "RSS.lower", "RSS.upper")
    if( (fig.opt=="TRUE" | fig.opt=="T" | fig.opt=="True")){ 
      dev.new()
      par(family="serif")
      par(mar=c(5,5,2,2))
      plot(tc, wald.lci.yj, cex.lab=1.5, cex.axis=1.5, type="l", 
           col=1, ylim=ylim0, 
           xlab=expression(paste(italic(t), {""}[paste(alpha, "/2", sep="")], sep="")), 
           ylab=bquote(paste(theta, {""}[.(j)], sep="")))
      lines(tc, wald.uci.yj, col=1)
      points(x2.lower, y2.lower, pch=1, col=4, cex=1.5)
      points(x2.upper, y2.upper, pch=1, col=2, cex=1.5)
      if(!is.null(show.CI)){
        alpha2 <- 1 - show.CI
        cols   <- hcl.colors(length(alpha2))
        tc2    <- qt(1-alpha2/2, df=df)
        abline(v=tc2, col=cols)
        for(i in 1:length(show.CI)){
          text(tc2[i], Pj.hat, bquote(paste(.(show.CI[i]*100), 
            "%", sep="")), cex=1.5, col=cols[i], srt=90)
        }
      }
    }  
    parlist[[j]] <- list(pari=theta.hat[j], WaldCI=WaldCI, lhCI=lhCI)
  }
  list(partab=res2$partab, parlist=parlist)
} 
