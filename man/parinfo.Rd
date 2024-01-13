\name{parinfo}
\alias{parinfo}
\title{
Detailed Information of Estimated Model Parameters 
}
\description{
Provides the estimates, standard errors, confidence intervals, 
Jacobian matrix, and the covariance matrix of model parameters.
}
\usage{
parinfo(object, x, CI = 0.95, method = "Richardson", 
        method.args = list(eps = 1e-04, d = 0.11, 
        zero.tol = sqrt(.Machine$double.eps/7e-07), r = 6, 
        v = 2, show.details = FALSE), side = NULL)
}
\arguments{
  \item{object}{A fitted model object for which there exist the model expression(\code{expr}), the sample size (\code{sample.size} or \code{n}), estimate(s) of model parameter(s) (\code{par}), and residual sum of squares (\code{RSS})}
  \item{x}{A vector or a matrix of observations of independent variable(s)}
  \item{CI}{The confidence level(s) of the required interval(s)}
  \item{method}{It is the same as the input argument of \code{method} of the \code{\link[numDeriv]{hessian}} function in package \pkg{numDeriv}}
  \item{method.args}{It is the same as the input argument of \code{method.args} of the \code{\link[numDeriv]{hessian}} function in package \pkg{numDeriv}}
  \item{side}{It is the same as the input argument of \code{side} of the \code{\link[numDeriv]{jacobian}} function in package \pkg{numDeriv}}
}
\details{
The \code{object} argument cannot be a list. It is a fitted model object from using the \code{\link{fitIPEC}} function.

}
\value{
  \item{D}{The Jacobian matrix of model parameters at all the \eqn{x} observations}
  \item{partab}{The estimates, standard errors and confidence intervals of model parameters}
  \item{covmat}{The covariance matrix of model parameters}
}
\references{
Bates, D.M and Watts, D.G. (1988) \emph{Nonlinear Regression Analysis and its Applications}. Wiley, New York. \doi{10.1002/9780470316757}

Ratkowsky, D.A. (1983) \emph{Nonlinear Regression Modeling: A Unified Practical Approach}. Marcel Dekker, New York.

Ratkowsky, D.A. (1990) \emph{Handbook of Nonlinear Regression Models}, Marcel Dekker, New York.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Peter M. Ridland \email{p.ridland@unimelb.edu.au}, 
David A. Ratkowsky \email{d.ratkowsky@utas.edu.au}, Yang Li \email{yangli@fau.edu}.
}
\note{
When there are \code{sample.size} and \code{n} in \code{object} at the same time, the default of 
the sample size is \code{sample.size}, which is superior to \code{n}. 
}
\seealso{
\code{\link{biasIPEC}}, \code{\link{confcurves}}, \code{\link{curvIPEC}}, \code{\link{skewIPEC}}, 
\code{\link[numDeriv]{hessian}} in package \pkg{numDeriv}, 
\code{\link[numDeriv]{jacobian}} in package \pkg{numDeriv}
}
\examples{
#### Example 1 ###################################################################################
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
r1       <- fitIPEC( MitA, x=x4, y=y4, ini.val=ini.val3, xlim=NULL, ylim=NULL,  
                     fig.opt=TRUE, control=list(
                     trace=FALSE, reltol=1e-20, maxit=50000) )
parA     <- r1$par
parA
result1  <- parinfo(r1, x=x4, CI=0.95)
result1

ini.val4 <- c(-0.10, 0.90, 2.5)
R0       <- fitIPEC( MitB, x=x4, y=y4, ini.val=ini.val4, xlim=NULL, ylim=NULL,  
                     fig.opt=TRUE, control=list(
                     trace=FALSE, reltol=1e-20, maxit=50000) )
parB     <- R0$par
parB
result2  <- parinfo(R0, x=x4, CI=0.95)
result2

ini.val6 <- c(-0.15, 2.52, 1.09)
RES0     <- fitIPEC( MitC, x=x4, y=y4, ini.val=ini.val6, xlim=NULL, ylim=NULL,  
                     fig.opt=TRUE, control=list(trace=FALSE, 
                     reltol=1e-20, maxit=50000) )
parC     <- RES0$par
parC
result3  <- parinfo(RES0, x=x4, CI=0.95)
result3
##################################################################################################


#### Example 2 ###################################################################################
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

ini.val8 <- c(35, 0.1, 0.05, 0.2)
cons1    <- fitIPEC( isom.fun, x=X, y=Y, ini.val=ini.val8, control=list(
                     trace=FALSE, reltol=1e-20, maxit=50000) )
par8     <- cons1$par 
result2  <- parinfo(cons1, x=X, CI=0.95)
result2
##################################################################################################

graphics.off()
}
