\name{crops}
\alias{crops}
\title{
Whole-plant biomass Data of 12 Species of Crops
}
\description{
The whole-plant biomass data of 12 species of crops growing in northern China in 2011.
}
\usage{
data(crops)
}
\details{
In the data set, there are six columns: \code{Code}, \code{CommonName}, \code{Date},  
\code{Time}, \code{FM}, and \code{DM}. 
\code{Code} is used to save the codes of crops; 
\code{CommonName} is used to save the common names of crops; 
\code{Date} is used to save the investigation date; 
\code{Time } is used to save the ages of crops from the sowing date (27 June, 2011) in days; 
\code{FM} is used to save the whole-plant fresh mass of crops in g; 
\code{DM} is used to save the whole-plant dry mass of crops in g. 
 
\code{Code = 1} represents sunflowers; 

\code{Code = 2} represents peanuts; 

\code{Code = 3} represents black soybeans; 

\code{Code = 4} represents soybeans; 

\code{Code = 5} represents kidney beans; 

\code{Code = 6} represents garden peas; 

\code{Code = 7} represents adzuki beans; 

\code{Code = 8} represents mungbeans; 

\code{Code = 9} represents cottons; 

\code{Code = 10} represents sweet sorghums; 

\code{Code = 11} represents corns; 

\code{Code = 12} represents Mexican corns. 
}
\references{
Shi, P., Fan, M., Ratkowsky, D.A., Huang, J., Wu, H., Chen, L., Fang, S. and 
Zhang, C. (2017) Comparison of two ontogenetic growth equations for animals and plants. 
\emph{Ecol. Model.} 349, 1\eqn{-}10. \doi{10.1016/j.ecolmodel.2017.01.012}

Shi, P., Men, X., Sandhu, H.S., Chakraborty, A., Li, B., Ouyang, F., Sun, Y., Ge, F. 
(2013) The "general" ontogenetic growth model is inapplicable to crop growth. 
\emph{Ecol. Model.} 266, 1\eqn{-}9. \doi{10.1016/j.ecolmodel.2013.06.025}
}
\examples{
data(crops)
ind   <- 6
xv    <- crops$Time[crops$Code == ind]
yv    <- crops$DM[crops$Code == ind] 
xlab0 <- "Time (d)"
ylab0 <- "Dry mass (g)"

dev.new()
plot(xv, yv, cex=1.5, cex.lab=1.5, cex.axis=1.5, xlab=xlab0, ylab=ylab0)

# Define the beta sigmoid model (bsm)
bsm <- function(P, x){
  P  <- cbind(P)
  if(length(P) !=4 ) {stop("The number of parameters should be 4!")}
  ropt <- P[1]
  topt <- P[2]
  tmin <- P[3]
  tmax <- P[4]
  tailor.fun <- function(x){
    x[x < tmin] <- tmin
    x[x > tmax] <- tmax
    return(x)
  }
  x <- tailor.fun(x)   
  ropt*(x-tmin)*(x-2*tmax+topt)/(topt+tmin-2*tmax)*(
       (x-tmin)/(topt-tmin) )^((topt-tmin)/(tmax-topt))   
}

# For the original beta sigmoid model
ini.val0  <- c(60, 30, seq(0, 10, 20), 100)
fit1 <- fitIPEC( bsm, x=xv, y=yv, ini.val=ini.val0, xlim=NULL, ylim=NULL, 
                 xlab=xlab0, ylab=ylab0, fig.opt=TRUE, 
                 control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
fit1$par

w    <- rep(1/as.numeric(tapply(yv, xv, var)), tapply(yv, xv, length))
fit2 <- fitIPEC( bsm, x=xv, y=yv, ini.val=ini.val0, weights=w, xlim=NULL,  
                 ylim=NULL, xlab=xlab0, ylab=ylab0, fig.opt=TRUE, 
                 control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
fit2$par

dev.new()
xp  <- seq(0, 120, len=2000)
yp  <- bsm(P=fit2$par, x=xp)
xv2 <- as.numeric(tapply(xv, xv, mean))
yv2 <- as.numeric(tapply(yv, xv, mean))
sd2 <- as.numeric(tapply(yv, xv, sd))
Up  <- yv2+sd2
Low <- yv2-sd2
plot( xv2, yv2, xlab=xlab0, ylab=ylab0, cex.lab=1.5, 
      cex.axis=1.5, xlim=c(0,120), ylim=c(-5, 100), type="n" )
lines( xp, yp, col=4 )
points( xv2, yv2, pch=1, cex=1.5, col=2 )
for(i in 1:length(Up)){
  lines(c(xv2[i], xv2[i]), c(Low[i], Up[i]), col=6)
}  
}
