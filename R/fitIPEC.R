fitIPEC <-
function( expr, x, y, ini.val, weights=NULL, control=list(), fig.opt=TRUE, 
              xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL ){

  x <- rbind( x )
  y <- as.vector(y)
  if( min(dim(x))[1] == 1 )   x <- cbind( x )
  if( nrow(x) != length(y) )  x <- t(x)

  if(!is.null(weights) && !is.numeric(weights)) 
      stop("'weights' must be a numeric vector!")
  if(!is.null(weights) && length(weights) != length(y))
      stop("The length of 'weights' must be the same as that of 'y'!")

  object.fun <- function(P){
    y.theor <- expr(P, x)
    if(is.null(weights)){
        temp <- sum( (y-y.theor)^2 )
    }
    if(!is.null(weights)){
        temp <- sum(weights*(y-y.theor)^2)
    }
    return(temp)
  }
  ini.val <- as.list(ini.val)  
  p       <- length(ini.val)
  r       <- 1
  for(i in 1:p){
    r <- r * length( ini.val[[i]] )
  }
  ini.val <- expand.grid(ini.val)
  mat     <- matrix(NA, nrow=r, ncol=(p+1))
  for(i in 1:nrow(ini.val)){    
    res      <- optim( ini.val[i, ], object.fun, control=control ) 
    mat[i, ] <- c(res$par, res$val)
  }   

  Names <- rep(NA, len=p)
  for(k in 1:p){
    Names[k] <- paste("theta[", k, "]", sep="")
  }
  colnames(mat) <- c(Names, "RSS")
  ind    <- which( mat[, p+1] == min(mat[, p+1])[1] )[1]
  par    <- as.vector( mat[ind, 1:p] )        
  yhat   <- expr(par, x)
  if(is.null(weights)){
    RSS    <- sum((yhat - y)^2)
    R.sq   <- 1 - RSS/sum((y-mean(y))^2)
  }
  if(!is.null(weights)){
    RSS    <- sum(weights*(yhat - y)^2)
    R.sq   <- 1 - RSS/sum(weights*(y-mean(y))^2)
  }

  if( (fig.opt=="TRUE" | fig.opt=="T" | fig.opt=="True") & ncol(x) == 1 ){ 
    dev.new()
    par(family="serif")
    par(mar=c(5,5,2,2))
    if( is.null(xlim) ){
        xmin  <- min(x)[1]
        xmax  <- max(x)[1]
        x.int <- (xmax - xmin)/8
        x2    <- seq(xmin-x.int, xmax+x.int, len=2000)
    }
    else{
      if(xlim[1] >= xlim[2]) stop("The first element should be less than the second in 'xlim'")
      x2 <- seq(xlim[1], xlim[2], len=2000)
    }
    y2 <- expr(par, x2) 
    if( is.null(ylim) ){
        if( is.null(xlab) & is.null(ylab) )
            plot(x2, y2, cex.lab=1.5, cex.axis=1.5, type="n", xlim=xlim, xlab="x", ylab="y") 
        if( !is.null(xlab) & is.null(ylab) )
            plot(x2, y2, cex.lab=1.5, cex.axis=1.5, type="n", xlim=xlim, xlab=xlab, ylab="y")
        if( is.null(xlab) & !is.null(ylab) )
            plot(x2, y2, cex.lab=1.5, cex.axis=1.5, type="n", xlim=xlim, xlab="x", ylab=ylab)
        if( !is.null(xlab) & !is.null(ylab) )
            plot(x2, y2, cex.lab=1.5, cex.axis=1.5, type="n", xlim=xlim, xlab=xlab, ylab=ylab)
    }
    if( !is.null(ylim) ){
        if(ylim[1] >= ylim[2]) 
            stop("The first element should be less than the second in 'ylim'")
        if( is.null(xlab) & is.null(ylab) )
            plot(x2, y2, cex.lab=1.5, cex.axis=1.5, type="n", xlim=xlim, ylim=ylim, xlab="x", ylab="y") 
        if( !is.null(xlab) & is.null(ylab) )
            plot(x2, y2, cex.lab=1.5, cex.axis=1.5, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab="y")
        if( is.null(xlab) & !is.null(ylab) )
            plot(x2, y2, cex.lab=1.5, cex.axis=1.5, type="n", xlim=xlim, ylim=ylim, xlab="x", ylab=ylab)
        if( !is.null(xlab) & !is.null(ylab) )
            plot(x2, y2, cex.lab=1.5, cex.axis=1.5, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
    }
    lines(x2, y2, col=4, lwd=1)
    points(x, y, pch=1, cex=1.5, col=2) 
  }

  if( (fig.opt=="TRUE" | fig.opt=="T" | fig.opt=="True") & ncol(x) > 1 ){ 
    dev.new()
    par(family="serif")
    par(mar=c(5,5,2,2))
    if( is.null(xlim) ){
      x.int <- (max(y)[1] - min(y)[1])/8
      xlim  <- c(min(y)[1] - x.int, max(y)[1] + x.int)
    }
    else{
      if(xlim[1] >= xlim[2]) stop("The first element should be less than the second in 'xlim'")      
    }
    if( is.null(ylim) ){
        ylim <- xlim
        if( is.null(xlab) & is.null(ylab) )
            plot(y, yhat, type="n", cex.lab=1.5, cex.axis=1.5, xlim=xlim, xlab="Observations of response variable", ylab="Predicted values of response variable") 
        if( !is.null(xlab) & is.null(ylab) )
            plot(y, yhat, type="n", cex.lab=1.5, cex.axis=1.5, xlim=xlim, xlab=xlab, ylab="Predicted values of response variable")
        if( is.null(xlab) & !is.null(ylab) )
            plot(y, yhat, type="n", cex.lab=1.5, cex.axis=1.5, xlim=xlim, xlab="Observations of response variable", ylab=ylab)
        if( !is.null(xlab) & !is.null(ylab) )
            plot(y, yhat, type="n", cex.lab=1.5, cex.axis=1.5, xlim=xlim, xlab=xlab, ylab=ylab)
    }
    if( !is.null(ylim) ){
        if(ylim[1] >= ylim[2]) 
            stop("The first element should be less than the second in 'ylim'")
        if( is.null(xlab) & is.null(ylab) )
            plot(y, yhat, type="n", cex.lab=1.5, cex.axis=1.5, xlim=xlim, ylim=ylim, xlab="Observations of response variable", ylab="Predicted values of response variable") 
        if( !is.null(xlab) & is.null(ylab) )
            plot(y, yhat, type="n", cex.lab=1.5, cex.axis=1.5, xlim=xlim, ylim=ylim, xlab=xlab, ylab="Predicted values of response variable")
        if( is.null(xlab) & !is.null(ylab) )
            plot(y, yhat, type="n", cex.lab=1.5, cex.axis=1.5, xlim=xlim, ylim=ylim, xlab="Observations of response variable", ylab=ylab)
        if( !is.null(xlab) & !is.null(ylab) )
            plot(y, yhat, type="n", cex.lab=1.5, cex.axis=1.5, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
    }
    abline(0, 1, col=4, lwd=1)
    points(y, yhat, pch=1, cex=1.5, col=2)
  }

  list(expr=expr, mat=mat, par=par, RSS=RSS, R.sq=R.sq, n=length(y))

}
