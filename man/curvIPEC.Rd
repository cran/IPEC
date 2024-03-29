\name{curvIPEC}
\alias{curvIPEC}
\title{
RMS Curvature Calculation Function
}
\description{
Calculates the root mean square curvatures (intrinsic and parameter-effects curvatures) 
of a nonlinear regression model.
}
\usage{
curvIPEC(expr, theta, x, y, tol = 1e-16, alpha = 0.05, method = "Richardson", 
         method.args = list(eps = 1e-04, d = 0.11, 
         zero.tol = sqrt(.Machine$double.eps/7e-07), 
         r = 6, v = 2, show.details = FALSE), side = NULL)
}
\arguments{
  \item{expr}{A given parametric model}
  \item{theta}{A vector of parameters of the model}
  \item{x}{A vector or matrix of observations of independent variable(s)}
  \item{y}{A vector of observations of response variable}
  \item{tol}{The tolerance for detecting linear dependencies in the columns of a matrix in the QR decomposition. 
See the input argument of \code{tol} of the \code{\link{qr}} function in package \pkg{base} }
  \item{alpha}{Parameter controlling the significance level for testing the significance of a curvature}
  \item{method}{It is the same as the input argument of \code{method} of the \code{\link[numDeriv]{hessian}} function in package \pkg{numDeriv}}
  \item{method.args}{It is the same as the input argument of \code{method.args} of the \code{\link[numDeriv]{hessian}} function in package \pkg{numDeriv}}
  \item{side}{It is the same as the input argument of \code{side} of the \code{\link[numDeriv]{jacobian}} function in package \pkg{numDeriv}}
}
\details{
This function was built based on the \code{\link{hessian}} and \code{\link{jacobian}} functions in package \pkg{numDeriv}, 
with reference to the \code{\link{rms.curv}} function in package \pkg{MASS}. 
However, it is more general without being limited by the \code{\link{deriv3}} function in package \pkg{stats} 
and \code{\link{nls}} class like the \code{\link{rms.curv}} function in package \pkg{MASS}. It mainly relies on package \pkg{numDeriv}.
The users only need provide the defined model, the fitted parameter vector, and the observations 
of independent and response variables, they will obtain the curvatures. The input argument \code{theta} can be obtained using the \code{\link{fitIPEC}} 
function in the current package, and it also can be obtained using the other nonlinear regression functions.
}
\value{
  \item{rms.ic}{The root mean square intrinsic curvature}
  \item{rms.pec}{The root mean square parameter-effects curvature}
  \item{critical.c}{The critical curvature value}
}
\references{
Bates, D.M and Watts, D.G. (1988) \emph{Nonlinear Regression Analysis and its Applications}. Wiley, New York. \doi{10.1002/9780470316757}

Gebremariam, B. (2014) Is nonlinear regression throwing you a curve? New diagnostic and inference tools in the NLIN Procedure. 
Paper SAS384-2014. \url{http://support.sas.com/resources/papers/proceedings14/SAS384-2014.pdf}

Ratkowsky, D.A. (1983) \emph{Nonlinear Regression Modeling: A Unified Practical Approach}. Marcel Dekker, New York.

Ratkowsky, D.A. (1990) \emph{Handbook of Nonlinear Regression Models}, Marcel Dekker, New York.

Ratkowsky, D.A. & Reddy, G.V.P. (2017) Empirical model with excellent statistical properties for describing temperature-dependent 
developmental rates of insects and mites. \emph{Ann. Entomol. Soc. Am.} 110, 302\eqn{-}309. \doi{10.1093/aesa/saw098}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Peter M. Ridland \email{p.ridland@unimelb.edu.au}, 
David A. Ratkowsky \email{d.ratkowsky@utas.edu.au}, Yang Li \email{yangli@fau.edu}.
}
\note{The calculation precision of curvature mainly depends on the setting of \code{method.args}. 
The two important default values in the list of \code{method.args} are \code{d = 0.11}, and \code{r = 6}.

This function cannot be used to calculate the maximum intrinsic and parameter-effects curvatures.
}
\seealso{
\code{\link{derivIPEC}}, \code{\link[numDeriv]{hessian}} in package \pkg{numDeriv}, 
\code{\link[numDeriv]{jacobian}} in package \pkg{numDeriv}, \code{\link{rms.curv}} in package \pkg{MASS}
}
\examples{
#### Example 1 ##################################################################################
# The velocity of the reaction (counts/min^2) under different substrate concentrations 
#   in parts per million (ppm) (Pages 255 and 269 of Bates and Watts 1988)

x1 <- c(0.02, 0.02, 0.06, 0.06, 0.11, 0.11, 0.22, 0.22, 0.56, 0.56, 1.10, 1.10)
y1 <- c(76, 47, 97, 107, 123, 139, 159, 152, 191, 201, 207, 200)

# Define the Michaelis-Menten model
MM <- function(theta, x){
    theta[1]*x / ( theta[2] + x )    
}

par1 <- c(212.68490865, 0.06412421)
# To calculate curvatures
res2 <- curvIPEC(MM, theta=par1, x=x1, y=y1, alpha=0.05, method="Richardson",
            method.args=list(eps=1e-4, d=0.11, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2)) 
res2
##################################################################################################


#### Example 2 ###################################################################################
# Development data of female pupae of cotton bollworm (Wu et al. 2009)
# References:
#   Ratkowsky, D.A. and Reddy, G.V.P. (2017) Empirical model with excellent statistical 
#       properties for describing temperature-dependent developmental rates of insects  
#       and mites. Ann. Entomol. Soc. Am. 110, 302-309.
#   Wu, K., Gong, P. and Ruan, Y. (2009) Estimating developmental rates of 
#       Helicoverpa armigera (Lepidoptera: Noctuidae) pupae at constant and
#       alternating temperature by nonlinear models. Acta Entomol. Sin. 52, 640-650.

# 'x2' is the vector of temperature (in degrees Celsius)
# 'D2' is the vector of developmental duration (in d)
# 'y2' is the vector of the square root of developmental rate (in 1/d)

x2 <- seq(15, 37, by=1)
D2 <- c( 41.24,37.16,32.47,26.22,22.71,19.01,16.79,15.63,14.27,12.48,
         11.3,10.56,9.69,9.14,8.24,8.02,7.43,7.27,7.35,7.49,7.63,7.9,10.03 )
y2 <- 1/D2
y2 <- sqrt( y2 )

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
  if (Tmin >= Tmax | ropt <= 0 | Topt <= Tmin | Topt >= Tmax) 
    temp <- Inf
  if (Tmax > Tmin & ropt > 0 & Topt > Tmin & Topt < Tmax){
    temp <- sqrt( ropt*(x-Tmax)*(x-Tmin)^2/((Topt-Tmin)*((Topt-Tmin
      )*(x-Topt)-(Topt-Tmax)*(Topt+Tmin-2*x))) )  
  }
  return( temp )
}

myfun <- sqrt.LRF
par2  <- c(0.1382926, 33.4575663, 5.5841244, 38.8282021)

# To calculate curvatures
resu2 <- curvIPEC( myfun, theta=par2, x=x2, y=y2, alpha=0.05, method="Richardson", 
                   method.args=list(eps=1e-4, d=0.11, 
                   zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
resu2
##################################################################################################


#### Example 3 ###################################################################################
# Height growth data of four species of bamboo (Gramineae: Bambusoideae)
# Reference(s):
# Shi, P., Fan, M., Ratkowsky, D.A., Huang, J., Wu, H., Chen, L., Fang, S. and  
#     Zhang, C. (2017) Comparison of two ontogenetic growth equations for animals and plants. 
#     Ecol. Model. 349, 1-10.

data(shoots)
# Choose a species
# 1: Phyllostachys iridescens; 2: Phyllostachys mannii; 
# 3: Pleioblastus maculatus; 4: Sinobambusa tootsik. 
# 'x3' is the vector of the investigation times (in d) from a specific starting time of growth
# 'y3' is the vector of the aboveground height values (in cm) of bamboo shoots at 'x3' 

ind <- 4
x3  <- shoots$x[shoots$Code == ind]
y3  <- shoots$y[shoots$Code == ind] 

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
     (x-tmin)/(topt-tmin))^((topt-tmin)/(tmax-topt))   
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
  ropt*(x-tmin)*(x-2*tmax+topt)/(topt+tmin-2*tmax)*(
  (x-tmin)/(topt-tmin))^((topt-tmin)/(tmax-topt))   
}

# For the original beta sigmoid model 
ini.val2 <- c(40, 30, 5, 50)
xlab2    <- "Time (d)"
ylab2    <- "Height (cm)"
re0      <- fitIPEC( bsm, x=x3, y=y3, ini.val=ini.val2, 
                     xlim=NULL, ylim=NULL, xlab=xlab2, ylab=ylab2, 
                     fig.opt=TRUE, control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
par3  <- re0$par
par3
re1   <- derivIPEC( bsm, theta=par3, x3[20], method="Richardson", 
                    method.args=list(eps=1e-4, d=0.11, 
                    zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) )
re1
re2   <- curvIPEC( bsm, theta=par3, x=x3, y=y3, alpha=0.05, method="Richardson",                    
                   method.args=list(eps=1e-4, d=0.11, 
                   zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
re2

# For the simplified beta sigmoid model (in comparison with the original beta sigmoid model)
ini.val7 <- c(40, 30, 50)

RESU0 <- fitIPEC( simp.bsm, x=x3, y=y3, ini.val=ini.val7, 
                  xlim=NULL, ylim=NULL, xlab=xlab2, ylab=ylab2, 
                  fig.opt=TRUE, control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
par7  <- RESU0$par
par7

RESU2 <- curvIPEC( simp.bsm, theta=par7, x=x3, y=y3, alpha=0.05, method="Richardson",             
                   method.args=list(eps=1e-4, d=0.11, 
                   zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
RESU2
##################################################################################################


#### Example 4 ###################################################################################
# Weight of cut grass data (Pattinson 1981)
# References:
#   Clarke, G.P.Y. (1987) Approximate confidence limits for a parameter function in nonlinear 
#       regression. J. Am. Stat. Assoc. 82, 221-230.
#   Gebremariam, B. (2014) Is nonlinear regression throwing you a curve? 
#       New diagnostic and inference tools in the NLIN Procedure. Paper SAS384-2014.
#       http://support.sas.com/resources/papers/proceedings14/SAS384-2014.pdf
#   Pattinson, N.B. (1981) Dry Matter Intake: An Estimate of the Animal
#       Response to Herbage on Offer. unpublished M.Sc. thesis, University
#       of Natal, Pietermaritzburg, South Africa, Department of Grassland Science.

# 'x4' is the vector of weeks after commencement of grazing in a pasture
# 'y4' is the vector of weight of cut grass from 10 randomly sited quadrants

x4 <- 1:13
y4 <- c(3.183, 3.059, 2.871, 2.622, 2.541, 2.184, 
        2.110, 2.075, 2.018, 1.903, 1.770, 1.762, 1.550)

# Define the first case of Mitscherlich equation
MitA <- function(P1, x){
    P1[3] + P1[2]*exp(P1[1]*x)
}

# Define the second case of Mitscherlich equation
MitB <- function(P2, x){
    log( P2[3] ) + exp(P2[2] + P2[1]*x)
}

# Define the third case of Mitscherlich equation
MitC <- function(P3, x, x1=1, x2=13){
    theta1 <- P3[1]
    beta2  <- P3[2]
    beta3  <- P3[3]
    theta2 <- (beta3 - beta2)/(exp(theta1*x2)-exp(theta1*x1))
    theta3 <- beta2/(1-exp(theta1*(x1-x2))) - beta3/(exp(theta1*(x2-x1))-1)
    theta3 + theta2*exp(theta1*x)
}

ini.val3 <- c(-0.1, 2.5, 1)
r0       <- fitIPEC( MitA, x=x4, y=y4, ini.val=ini.val3, xlim=NULL, ylim=NULL,  
                     fig.opt=TRUE, control=list(
                     trace=FALSE, reltol=1e-20, maxit=50000) )
parA     <- r0$par
parA
r2 <- curvIPEC( MitA, theta=parA, x=x4, y=y4, alpha=0.05, method="Richardson", 
                method.args=list(eps=1e-4, d=0.11, 
                zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
r2

ini.val4 <- c(exp(-0.1), log(2.5), 1)

R0       <- fitIPEC( MitB, x=x4, y=y4, ini.val=ini.val3, xlim=NULL, ylim=NULL,  
                     fig.opt=TRUE, control=list(
                     trace=FALSE, reltol=1e-20, maxit=50000) )
parB     <- R0$par
parB
R2       <- curvIPEC( MitB, theta=parB, x=x4, y=y4, alpha=0.05, method="Richardson", 
                      method.args=list(eps=1e-4, d=0.11, 
                      zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
R2

ini.val6 <- c(-0.15, 2.52, 1.09)
RES0     <- fitIPEC( MitC, x=x4, y=y4, ini.val=ini.val6, xlim=NULL, ylim=NULL,  
                     fig.opt=TRUE, control=list(trace=FALSE, 
                     reltol=1e-20, maxit=50000) )
parC     <- RES0$par
parC
RES2     <- curvIPEC( MitC, theta=parC, x=x4, y=y4, 
                      tol=1e-20, alpha=0.05, method="Richardson", 
                      method.args=list(eps=1e-4, d=0.11, 
                      zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
RES2
##################################################################################################


#### Example 5 ###################################################################################
# Conductance of a thermistor (y5) as a function of temperature (x5) (Meyer and Roth, 1972)
# References:
#   Page 120 in Ratkowsky (1983)
#   Meyer, R.R. and Roth P.M. (1972) Modified damped least squares:
#       A algorithm for non-linear estimation. J. Inst. Math. Appl. 9, 218-233.

x5 <- seq(50, 125, by=5)
y5 <- c( 34780, 28610, 23650, 19630, 16370, 13720, 11540, 9744, 
         8261, 7030, 6005, 5147, 4427, 3820, 3307, 2872 )
y5 <- log(y5)

conduct.fun <- function(P, x){
-P[1]+P[2]/(x+P[3])
}

ini.val5 <- c(5, 10^4, 0.5*10^3)
RE0      <- fitIPEC( conduct.fun, x=x5, y=y5, ini.val=ini.val5, xlim=NULL, ylim=NULL,  
                     fig.opt=TRUE, control=list(
                     trace=FALSE, reltol=1e-20, maxit=50000) )
par5     <- RE0$par
par5
RE2      <- curvIPEC( conduct.fun, theta=par5, x=x5, y=y5, alpha=0.05, method="Richardson", 
                      method.args=list(eps=1e-4, d=0.11, 
                      zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
RE2
##################################################################################################


#### Example 6 ###################################################################################
# Data on biochemical oxygen demand (BOD; Marske 1967)
# References
# Pages 255 and 270 in Bates and Watts (1988)
# Marske, D. (1967) Biochemical oxygen demand data interpretation using sum of squares surface.
#     M.Sc. Thesis, University of Wisconsin-Madison.

# 'x6' is a vector of time (in d)
# 'y6' is a vector of biochemical oxygen demand (mg/l)

x6 <- c(1, 2, 3, 4, 5, 7)
y6 <- c(8.3, 10.3, 19.0, 16.0, 15.6, 19.8)

BOD.fun <- function(P, x){
  P[1]*(1-exp(P[2]*x))
}

ini.val7 <- c(210, 0.06)
consq0   <- fitIPEC( BOD.fun, x=x6, y=y6, ini.val=ini.val7, xlim=NULL, ylim=NULL,  
                     fig.opt=TRUE, control=list(
                     trace=FALSE, reltol=1e-20, maxit=50000) )
par7     <- consq0$par
par7
consq2   <- curvIPEC( BOD.fun, theta=par7, x=x6, y=y6, alpha=0.05, method="Richardson", 
                      method.args=list(eps=1e-4, d=0.11, 
                      zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
consq2
##################################################################################################


#### Example 7 ###################################################################################
# Data on biochemical oxygen demand (BOD; Marske 1967)
# References:
# Pages 56, 255 and 271 in Bates and Watts (1988)
# Carr, N.L. (1960) Kinetics of catalytic isomerization of n-pentane. Ind. Eng. Chem.
#     52, 391-396.   

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

par8  <- c(35.92831619, 0.07084811, 0.03772270, 0.16718384) 
cons2 <- curvIPEC( isom.fun, theta=par8, x=X, y=Y, alpha=0.05, method="Richardson", 
                   method.args=list(eps=1e-4, d=0.11, 
                   zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) ) 
cons2
##################################################################################################
}
