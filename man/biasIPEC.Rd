\name{biasIPEC}
\alias{biasIPEC}
\title{
Bias Calculation Function
}
\description{
Calculates the bias in the estimates of the parameters of a given nonlinear regression model.
}
\usage{
biasIPEC(expr, theta, x, y, tol = 1e-16, method = "Richardson", 
         method.args = list(eps = 1e-04, d = 0.11, 
         zero.tol = sqrt(.Machine$double.eps/7e-07), r = 6, v = 2, 
         show.details = FALSE), side = NULL)
}
\arguments{
  \item{expr}{A given nonlinear regression model}
  \item{theta}{A vector of parameters of the model}
  \item{x}{A vector or matrix of observations of independent variable(s)}
  \item{y}{A vector of observations of response variable}
  \item{tol}{The tolerance for detecting linear dependencies in the columns of a matrix for calculating its inverse. 
See the input argument of \code{tol} of the \code{\link{solve}} function in package \pkg{base} }
  \item{method}{It is the same as the input argument of \code{method} of the \code{\link[numDeriv]{hessian}} function in package \pkg{numDeriv}}
  \item{method.args}{It is the same as the input argument of \code{method.args} of the \code{\link[numDeriv]{hessian}} function in package \pkg{numDeriv}}
  \item{side}{It is the same as the input argument of \code{side} of the \code{\link[numDeriv]{jacobian}} function in package \pkg{numDeriv}}
}
\details{
The defined model should have two input arguments: a parameter vector and an independent variable vector or matrix, 
e.g. \code{myfun <- function(P, x)\{...\}}, where \code{P} represents the parameter vector and \code{x} represents the independent variable vector or matrix. 

An absolute value of \code{percent.bias} (see below) in excess of \strong{1\%} appears to be a good rule of thumb for indicating nonlinear behavior (Ratkowsky 1983).
}
\value{
  \item{bias}{The bias}
  \item{percent.bias}{The percentage bias that is equal to \code{bias/estimate * 100\%}}
}
\references{
Box, M.J. (1971) Bias in nonlinear estimation. \emph{J. R. Statist. Soc., Ser. B} 33, 171\eqn{-}201. 
\doi{10.1111/j.2517-6161.1971.tb00871.x}

Ratkowsky, D.A. (1983) \emph{Nonlinear Regression Modeling: A Unified Practical Approach}. Marcel Dekker, New York.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Peter M. Ridland \email{p.ridland@unimelb.edu.au}, 
David A. Ratkowsky \email{d.ratkowsky@utas.edu.au}, Yang Li \email{yangli@fau.edu}.
}
\note{
The current function can be applicable to nonlinear models with multiple independent variables.
}
\seealso{
\code{\link{derivIPEC}}, \code{\link[numDeriv]{hessian}} 
in package \pkg{numDeriv}, \code{\link[numDeriv]{jacobian}} in package \pkg{numDeriv}
}
\examples{
#### Example 1 #################################################################################
# The velocity of the reaction (counts/min^2) under different substrate concentrations 
#   in parts per million (ppm) (Page 269 of Bates and Watts 1988)
x1 <- c(0.02, 0.02, 0.06, 0.06, 0.11, 0.11, 0.22, 0.22, 0.56, 0.56, 1.10, 1.10)
y1 <- c(76, 47, 97, 107, 123, 139, 159, 152, 191, 201, 207, 200)

# Define the Michaelis-Menten (MM) model
MM <- function(theta, x){
    theta[1]*x / ( theta[2] + x )    
}

par1 <- c(212.68490865, 0.06412421)
res3 <- biasIPEC(MM, theta=par1, x=x1, y=y1, tol= 1e-20)
res3
#################################################################################################


#### Example 2 ##################################################################################
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
D2 <- c(41.24,37.16,32.47,26.22,22.71,19.01,16.79,15.63,14.27,12.48,
       11.3,10.56,9.69,9.14,8.24,8.02,7.43,7.27,7.35,7.49,7.63,7.9,10.03)
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

# To calculate bias
resu3 <- biasIPEC(myfun, theta=par2, x=x2, y=y2, tol= 1e-20)
resu3
#################################################################################################


#### Example 3 ##################################################################################
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
y4 <- c(3.183, 3.059, 2.871, 2.622, 2.541, 2.184, 2.110, 2.075, 2.018, 1.903, 1.770, 1.762, 1.550)

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
                     fig.opt=TRUE, control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
parA     <- r0$par
parA
r3       <- biasIPEC( MitA, theta=parA, x=x4, y=y4, tol=1e-20 ) 
r3

ini.val4 <- c(exp(-0.1), log(2.5), 1)
R0       <- fitIPEC( MitB, x=x4, y=y4, ini.val=ini.val3, xlim=NULL, ylim=NULL,  
                     fig.opt=TRUE, control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
parB     <- R0$par
parB
R3       <- biasIPEC( MitB, theta=parB, x=x4, y=y4, tol=1e-20 ) 
R3

ini.val6 <- c(-0.15, 2.52, 1.09)
RES0     <- fitIPEC( MitC, x=x4, y=y4, ini.val=ini.val6, xlim=NULL, ylim=NULL,  
                     fig.opt=TRUE, control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
parC     <- RES0$par
parC
RES3     <- biasIPEC(MitC, theta=parC, x=x4, y=y4, tol=1e-20) 
RES3
#################################################################################################


#### Example 4 ###################################################################################
# Data on biochemical oxygen demand (BOD; Marske 1967)
# References
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
cons3 <- biasIPEC( isom.fun, theta=par8, x=X, y=Y, tol= 1e-20 )
cons3
#################################################################################################
}
