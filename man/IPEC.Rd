\name{IPEC}
\alias{IPEC}
\docType{package}
\title{
Root Mean Square Curvature Calculation
}
\description{
Calculates the RMS intrinsic and parameter-effects curvatures of a nonlinear regression model.
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{IPEC}
\packageIndices{IPEC}

}
\note{
We are deeply thankful to Drs. Paul Gilbert and Jinlong Zhang for their invaluable help 
during creating this package.
}
\author{
Peijian Shi, Peter Ridland, David A. Ratkowsky

Maintainer: Peijian Shi <peijianshi@gmail.com>
}
\references{
Bates, D.M and Watts, D.G. (1988) \emph{Nonlinear Regression Analysis and its Applications}. Wiley, New York.

Ratkowsky, D.A. (1983) \emph{Nonlinear Regression Modeling: A Unified Practical Approach}. Marcel Dekker, New York.

Ratkowsky, D.A. (1990) \emph{Handbook of Nonlinear Regression Models}, Marcel Dekker, New York.
}

\keyword{ package }
\seealso{ \code{\link[numDeriv]{hessian}} in package \pkg{numDeriv}, \code{\link[numDeriv]{jacobian}} 
in package \pkg{numDeriv}, \code{\link[MASS]{rms.curv}} in package \pkg{MASS}}
\examples{
#### Example 1 ##################################################################################
graphics.off()
# The velocity of the reaction (counts/min^2) under different substrate concentrations 
#   in parts per million (ppm) (Page 269 of Bates and Watts 1988)

x1 <- c(0.02, 0.02, 0.06, 0.06, 0.11, 0.11, 0.22, 0.22, 0.56, 0.56, 1.10, 1.10)
y1 <- c(76, 47, 97, 107, 123, 139, 159, 152, 191, 201, 207, 200)

# Define the Michaelis-Menten model
MM <- function(theta, x){
    theta[1]*x / ( theta[2] + x )    
}

res0 <- fitIPEC( MM, x=x1, y=y1, ini.val=c(200, 0.05), 
                 xlim=c( 0, 1.5 ), ylim=c(0, 250), fig.opt=TRUE )
par1 <- res0$par
par1

res1 <- derivIPEC( MM, theta=par1, z=x1[1], method="Richardson",
                   method.args=list(eps=1e-4, d=0.11, 
                   zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) )
res1

# To calculate curvatures
res2 <- curvIPEC( MM, theta=par1, x=x1, y=y1, alpha=0.05, method="Richardson",
                  method.args=list(eps=1e-4, d=0.11, 
                  zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
res2

# To calculate bias
res3 <- biasIPEC(MM, theta=par1, x=x1, y=y1, tol= 1e-20)
res3

\dontrun{
  res4 <- bootIPEC( MM, x=x1, y=y1, ini.val=par1, target.fun = "RSS", 
                    control=list(reltol=1e-20, maxit=40000), 
                    nboot=2000, alpha=0.05, fig.opt=TRUE, seed=123 )
  res4
}

# To calculate skewness
res5 <- skewIPEC(MM, theta=par1, x=x1, y=y1, tol= 1e-20)
res5
#################################################################################################


\dontrun{
#### Example 2 ##################################################################################
graphics.off()
# Development data of female pupae of cotton bollworm (Wu et al. 2009)
# References:
#   Ratkowsky, D.A. and Reddy, G.V.P. (2017) Empirical model with excellent statistical 
#       properties for describing temperature-dependent developmental rates of insects  
#       and mites. Ann. Entomol. Soc. Am. 110, 302-309.
#   Wu, K.-J., Gong, P.-Y. and Ruan, Y.-M. (2009) Estimating developmental rates of 
#       Helicoverpa armigera (Lepidoptera: Noctuidae) pupae at constant and
#       alternating temperature by nonlinear models. Acta Entomol. Sin. 52, 640-650.

# 'x2' is the vector of temperature (in degrees Celsius)
# 'D2' is the vector of developmental duration (in d)
# 'y2' is the vector of the square root of developmental rate (in 1/d)

x2 <- seq(15, 37, by=1)
D2 <- c(41.24,37.16,32.47,26.22,22.71,19.01,16.79,15.63,14.27,12.48,
       11.3,10.56,9.69,9.14,8.24,8.02,7.43,7.27,7.35,7.49,7.63,7.9,10.03)
y2 <- 1/D2
y2 <- sqrt( y2 )

ini.val1 <- c(0.14, 30, 10, 40)

# Define the square root function of the Lobry-Rosso-Flandrois (LRF) model
sqrt.LRF <- function(P, x){
  ropt <- P[1]
  Topt <- P[2]
  Tmin <- P[3]
  Tmax <- P[4]
  fun0 <- function(z){
    z[z < Tmin] <- Tmin
    z[z > Tmax] <- Tmax
    return(z)
  }
  x <- fun0(x)
  sqrt( ropt*(x-Tmax)*(x-Tmin)^2/((Topt-Tmin)*((Topt-Tmin
      )*(x-Topt)-(Topt-Tmax)*(Topt+Tmin-2*x))) )
}

myfun <- sqrt.LRF
xlab1 <- expression( paste("Temperature (", degree, "C)", sep="" ) )
ylab1 <- expression( paste("Developmental rate"^{1/2}, " (", d^{"-1"}, ")", sep="") )
resu0 <- fitIPEC( myfun, x=x2, y=y2, ini.val=ini.val1, xlim=NULL, ylim=NULL, 
                  xlab=xlab1, ylab=ylab1, fig.opt=TRUE, 
                  control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
par2  <- resu0$par
par2

resu1 <- derivIPEC( myfun, theta=par2, z=x2[1], method="Richardson", 
                    method.args=list(eps=1e-4, d=0.11, 
                    zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) )
resu1

# To calculate curvatures
resu2 <- curvIPEC( myfun, theta=par2, x=x2, y=y2, alpha=0.05, method="Richardson", 
                   method.args=list(eps=1e-4, d=0.11, 
                   zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
resu2

# To calculate bias
resu3 <- biasIPEC(myfun, theta=par2, x=x2, y=y2, tol= 1e-20)
resu3

resu4 <- bootIPEC( myfun, x=x2, y=y2, ini.val=ini.val1, target.fun = "RSS", 
                   nboot=2000, alpha=0.05, fig.opt=TRUE, seed=123 )
resu4

# To calculate skewness
resu5 <- skewIPEC(myfun, theta=par2, x=x2, y=y2, tol= 1e-20)
resu5
#################################################################################################
}


#### Example 3 ##################################################################################
graphics.off()
# Height growth data of four species of bamboos (Gramineae: Bambusoideae)
# Reference(s):
# Shi, P.-J., Fan, M.-L., Ratkowsky, R.A., Huang, J.-G., Wu, H.-I, Chen, L., Fang, S.-Y. and 
#     Zhang, C.-X. (2017) Comparison of two ontogenetic growth equations for animals and plants. 
#     Ecol. Model. 349, 1-10.

data(shoots)
attach(shoots)
# Choose a species
# 1: Phyllostachys iridescens; 2: Phyllostachys mannii; 
# 3: Sinobambusa tootsik; 4: Pleioblastus maculatus
# 'x3' is the vector of the observation times from a specific starting time of growth
# 'y3' is the vector of the aboveground height values of bamboo shoots at 'x3' 

ind <- 3
x3  <- time[code == ind]
y3  <- height[code == ind] 

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

# Define the simplified beta sigmoid model (simp.bsm)
simp.bsm <- function(P, x, tmin=0){
  P  <- cbind(P)  
  ropt  <- P[1]
  topt  <- P[2]
  tmax  <- P[3]
  tailor.fun <- function(x){
    x[x < tmin] <- tmin
    x[x > tmax] <- tmax
    return(x)
  }
  x <- tailor.fun(x)   
  ropt*(x-tmin)*(x-2*tmax+topt)/(topt+tmin-2*tmax)*
       ((x-tmin)/(topt-tmin))^((topt-tmin)/(tmax-topt))   
}

# For the original beta sigmoid model
ini.val2 <- c(40, 30, 5, 50)
xlab2    <- "Time (d)"
ylab2    <- "Height (cm)"

re0 <- fitIPEC( bsm, x=x3, y=y3, ini.val=ini.val2, xlim=NULL, ylim=NULL, 
                xlab=xlab2, ylab=ylab2, fig.opt=TRUE, 
                control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
par3 <- re0$par
par3

re1 <- derivIPEC( bsm, theta=par3, x3[15], method="Richardson", 
                  method.args=list(eps=1e-4, d=0.11, zero.tol=
                  sqrt(.Machine$double.eps/7e-7), r=6, v=2) )
re1

re2 <- curvIPEC( bsm, theta=par3, x=x3, y=y3, alpha=0.05, method="Richardson", 
                 method.args=list(eps=1e-4, d=0.11, zero.tol=
                 sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
re2

re3 <- biasIPEC( bsm, theta=par3, x=x3, y=y3, tol= 1e-20 )
re3

\dontrun{
  re4 <- bootIPEC( bsm, x=x3, y=y3, ini.val=ini.val2, target.fun = "RSS", 
                   control=list(trace=FALSE, reltol=1e-20, maxit=50000),
                   nboot=2000, alpha=0.05, fig.opt=TRUE, fold=3.5, seed=123 )
  re4
}

re5 <- skewIPEC( bsm, theta=par3, x=x3, y=y3, tol= 1e-20 )
re5

# For the simplified beta sigmoid model 
#  (in comparison with the original beta sigmoid model)
ini.val7 <- c(40, 30, 50)

RESU0 <- fitIPEC( simp.bsm, x=x3, y=y3, ini.val=ini.val7, 
                  xlim=NULL, ylim=NULL, xlab=xlab2, ylab=ylab2, 
                  fig.opt=TRUE, control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
par7  <- RESU0$par
par7

RESU1 <- derivIPEC( simp.bsm, theta=par7, x3[15], method="Richardson", 
                    method.args=list(eps=1e-4, d=0.11, 
                    zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) )
RESU1

RESU2 <- curvIPEC( simp.bsm, theta=par7, x=x3, y=y3, alpha=0.05, method="Richardson", 
                   method.args=list(eps=1e-4, d=0.11, 
                   zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
RESU2

RESU3 <- biasIPEC( simp.bsm, theta=par7, x=x3, y=y3, tol= 1e-20 )
RESU3

\dontrun{
  RESU4 <- bootIPEC( simp.bsm, x=x3, y=y3, ini.val=ini.val7, target.fun = "RSS", 
                     control=list(trace=FALSE, reltol=1e-20, maxit=50000),
                     nboot=2000, alpha=0.05, fig.opt=TRUE, fold=3.5, seed=123 )
  RESU4
}

RESU5 <- skewIPEC( simp.bsm, theta=par7, x=x3, y=y3, tol= 1e-20 )
RESU5
##################################################################################################


#### Example 4 ###################################################################################
# Data on biochemical oxygen demand (BOD; Marske 1967)
# References:
# Pages 56, 255 and 271 in Bates and Watts (1988)
# Carr, N.L. (1960) Kinetics of catalytic isomerization of n-pentane. Ind. Eng. Chem.
#     52, 391-396.   

graphics.off()
data(isom)
Y <- isom[,1]
X <- isom[,2:4]

# There are three independent variables saved in matrix 'X' and one response variable (Y)
# The first column of 'X' is the vector of partial pressure of hydrogen
# The second column of 'X' is the vector of partial pressure of n-pentane
# The third column of 'X' is the vector of partial pressure of isopentane
# Y is the vector of experimental reaction rate (in 1/hr)

isom.fun <- function(theta, x){
  x1     <- x[,1]
  x2     <- x[,2]
  x3     <- x[,3]
  theta1 <- theta[1]
  theta2 <- theta[2]
  theta3 <- theta[3]
  theta4 <- theta[4]
  theta1*theta3*(x2-x3/1.632) / ( 1 + theta2*x1 + theta3*x2 + theta4*x3 )
}

ini.val8 <- c(35, 0.1, 0.05, 0.2)
cons1    <- fitIPEC( isom.fun, x=X, y=Y, ini.val=ini.val8, control=list(
                     trace=FALSE, reltol=1e-20, maxit=50000) )
par8     <- cons1$par 
cons2    <- curvIPEC( isom.fun, theta=par8, x=X, y=Y, alpha=0.05, method="Richardson", 
                      method.args=list(eps=1e-4, d=0.11, 
                      zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2)) 
cons2
cons3    <- biasIPEC( isom.fun, theta=par8, x=X, y=Y, tol= 1e-20 )
cons3

\dontrun{
  cons4 <- bootIPEC( isom.fun, x=X, y=Y, ini.val=ini.val8, target.fun = "RSS", 
                     control=list(trace=FALSE, reltol=1e-20, maxit=50000),
                     nboot=2000, alpha=0.05, fig.opt=TRUE, fold=10000, seed=123 )
  cons4
}

cons5    <- skewIPEC( isom.fun, theta=par8, x=X, y=Y, tol= 1e-20 )
cons5
##################################################################################################
}
