bootIPEC <-
function( expr, x, y, ini.val, target.fun = "RSS", 
                  control=list(), nboot=200, alpha=0.05, fig.opt=TRUE, fold=3.5, 
                  seed=NULL, unique.num=2, prog.opt=TRUE){

  x       <- rbind( x )
  y       <- as.vector(y)
  if( min(dim(x))[1] == 1 )   x <- cbind( x )
  if( nrow(x) != length(y) )  x <- t(x)
  ini.val <- as.list(ini.val) 
  res1    <- fitIPEC(expr, x, y, ini.val, target.fun, control, fig.opt=FALSE)
  n       <- nrow(x)
  p       <- length(ini.val)
  M0      <- c(res1$par, res1$RSS, res1$chi.sq, res1$R.sq)
  M       <- matrix( NA, nrow=nboot, ncol=(p+3) )
  if(!is.null(seed)){
      if(!is.numeric(seed)) stop("'seed' should be an integer")
  }
  OldSeed <- set.seed(seed)
  for(i in 1L:nboot){
    if(prog.opt=="TRUE" | prog.opt=="T"){
        # Sys.sleep(0.0005)
        # cat(i, paste(" of ", nboot, "\r", sep = ""))
        # flush.console()
        if(i %% 50 == 0){
          print(paste("The current running progress is ", i, "/", nboot, sep=""))
          cat("\n")
        }
        if (i %% nboot == 0)    cat("\n")
    }
    ind   <- sample( 1:n, n, replace=TRUE )
    xboot <- x[ind, ]
    yboot <- y[ind]
    # To set the permitted least non-overlapped number of 
    #     sampled data points for nonlinear regression
    if( nrow( unique(cbind(xboot, yboot)) ) < unique.num ){
        M[i,] <- NA
    }
    else{
        res2  <- fitIPEC( expr, xboot, yboot, ini.val=res1$par, target.fun=target.fun, 
                          control=control, fig.opt=FALSE ) 
        M[i,] <- c(res2$par, res2$RSS, res2$chi.sq, res2$R.sq)
    }   
  }
  M <- na.omit(M)

  #### To drop the extreme points ##################################
  inde <- c()
  for(j in 1:ncol(M)){
      v    <- M[, j]
      cl   <- fold * ( quantile(v, 0.75)[[1]] - quantile(v, 0.25)[[1]] )  
      inde <- c(inde, which( abs(v-median(v)) >= cl ))
   }
   inde <- sort( unique(inde) )   
   if(length(inde) > 0){
      M <- M[-inde,]
   }
  ##################################################################
  perc.ci.mat <- matrix(NA, nrow=p, ncol=6)



  # print( M )  


  
  for(j in 1:p){
    z                <- M[,j]
    lower            <- quantile(z, c(alpha, 1 - alpha))[[1]]
    upper            <- quantile(z, c(alpha, 1 - alpha))[[2]]
    perc.ci.mat[j, ] <- c(res1$par[j], sd(z), median(z), mean(z), lower, upper)  

    if( fig.opt=="TRUE" | fig.opt=="T" ){    
      dev.new()
      par(mar=c(5, 5, 2, 2))      
      z.int   <- ( max(z)[1] - min(z)[1] )/10
      z.range <- seq( min(z)[1]-z.int, max(z)[1]+z.int, len=2000 )
      den     <- dnorm(z.range, mean=mean(z), sd=sd(z))
      max.den <- max( c(hist(z, freq=FALSE)$density, den) )[1]
      e       <- bquote( expression(hat(theta)[.(j)]) )
      hist( z, freq=FALSE, cex.lab=1.5, cex.axis=1.5, xlab=eval(e), 
          main="", col="grey90", ylim=c(0, max.den*1.2))
      lines(z.range, den, col=2, lwd=2)
      abline(v=mean(z), lty=2, col=4, lwd=1)
      box()
    }
  }
  if( nboot >=2 ){  
    covar.mat <- matrix(NA, nrow=p, ncol=p)
    cor.mat   <- covar.mat
    for(i in 1L:p){
      z1 <- M[,i]
      e1 <- bquote( expression(hat(theta)[.(i)]) )
      for(j in 1L:p){
        z2 <- M[,j]
        e2 <- bquote( expression(hat(theta)[.(j)]) )
        covar.mat[i, j] <- sum( (z1-mean(z1)) * (z2-mean(z2)) ) / (nboot-1)
        cor.mat[i, j]   <- cor(z1, z2)
        if(j > i & fig.opt=="TRUE" | fig.opt=="T"){           
          dev.new()
          par(mar=c(5, 5, 2, 2))
          plot( z1, z2, pch=1, cex=1.5, cex.lab=1.5, cex.axis=1.5,
              xlab=eval(e1), ylab=eval(e2) )
          abline(v=res1$par[i], lty=2, col=3)
          abline(h=res1$par[j], lty=2, col=3)
        }
      }
    }
  }

  Names <- rep(NA, len=p)
  for(k in 1:p){
    Names[k] <- paste("theta[", k, "]", sep="")
  }
  rownames(perc.ci.mat) <- Names
  colnames(perc.ci.mat) <- c("Estimate", "SD", "Median", "Mean", "perc LCI", "perc UCI")
  colnames(M)           <- c(Names, "RSS", "chi.sq", "R.sq") 

  # Calculate the lower and upper limits of confidence intervals based on the BCa method
  number <- 0
  for (k in 1:p){
    number[k] <- sum(M[, k] < M0[k])
    }
  z0 <- qnorm(number/nrow(M))
  M1 <- matrix(NA, nrow=n, ncol=p)
  for (i in 1L:n){
    index3 <- 1:n
    xone   <- x[index3 != i, ]
    yone   <- y[index3 != i]
    res3   <- fitIPEC( expr, xone, yone, ini.val=res1$par, target.fun=target.fun, 
                       control=control, fig.opt=FALSE )   
    M1[i,] <- res3$par
  }
  a <- c()
  for (k in 1:p){
    a[k] <- sum((mean(M1[,k])-M1[,k])^3)/6/sum((M1[,k]-mean(M1[,k]))^2)^(3/2)
    }  
  ci.adj <- matrix(NA, nrow=p, ncol=2)
  alpha1 <- c()
  alpha2 <- c()
  for(k in 1:p){
    alpha1.temp  <- pnorm(z0[k]+(z0[k]+qnorm(alpha))/(1-a[k]*(z0[k]+qnorm(alpha))))
    alpha2.temp  <- pnorm(z0[k]+(z0[k]+qnorm(1-alpha))/(1-a[k]*(z0[k]+qnorm(1-alpha))))
    alpha1       <- c(alpha1, alpha1.temp)
    alpha2       <- c(alpha2, alpha2.temp)
    lower.temp   <- quantile(M[,k], c(alpha1.temp, alpha2.temp))[[1]]
    upper.temp   <- quantile(M[,k], c(alpha1.temp, alpha2.temp))[[2]]
    ci.adj[k,]   <- c(lower.temp, upper.temp)
    }

  bca.ci.mat      <- perc.ci.mat
  bca.ci.mat[, 5] <- ci.adj[,1]
  bca.ci.mat[, 6] <- ci.adj[,2] 

  rownames(bca.ci.mat) <- Names
  colnames(bca.ci.mat) <- c("Estimate", "SD", "Median", "Mean", "bca LCI", "bca UCI")

  if(nboot == 1){
      covar.mat <- NA; cor.mat <- NA
  }
  on.exit( set.seed(OldSeed) )
  list(M=M, perc.ci.mat=perc.ci.mat, bca.ci.mat=bca.ci.mat, covar.mat=covar.mat, cor.mat=cor.mat)
}
