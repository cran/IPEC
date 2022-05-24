\name{derivIPEC}
\alias{derivIPEC}
\title{
Derivative Calculation Function
}
\description{
Calculates the Jacobian and Hessian matrices of model parameters at a vector \code{z}.
}
\usage{
derivIPEC(expr, theta, z, method = "Richardson", 
          method.args = list(eps = 1e-04, d = 0.11, 
          zero.tol = sqrt(.Machine$double.eps/7e-07), r = 6, v = 2, 
          show.details = FALSE), side = NULL)
}
\arguments{
  \item{expr}{A given parametric model}
  \item{theta}{A vector of parameters of the model}
  \item{z}{A vector where the derivatives are calculated}
  \item{method}{It is the same as the input argument of \code{method} of the \code{\link[numDeriv]{hessian}} function in package \pkg{numDeriv}}
  \item{method.args}{It is the same as the input argument of \code{method.args} of the \code{\link[numDeriv]{hessian}} function in package \pkg{numDeriv}}
  \item{side}{It is the same as the input argument of \code{side} of the \code{\link[numDeriv]{jacobian}} function in package \pkg{numDeriv}}
}
\details{
The Hessian and Jacobian matrices are calculated at a vector \code{z}, 
which represents a value of a single independent variable 
or a combination of different values of multiple independent variables.
}
\value{
  \item{Jacobian}{The Jacobian matrix of parameters at \code{z}}
  \item{Hessian}{The Hessian matrix of parameters at \code{z}}
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
\seealso{
\code{\link{biasIPEC}}, \code{\link{skewIPEC}}, \code{\link{curvIPEC}}, 
\code{\link[numDeriv]{hessian}} in package \pkg{numDeriv}, 
\code{\link[numDeriv]{jacobian}} in package \pkg{numDeriv}
}
\examples{
#### Example 1 #####################################################################################
# Define the Michaelis-Menten model
MM <- function(theta, x){
    theta[1]*x / ( theta[2] + x )    
}

par1 <- c(212.68490865, 0.06412421)
res1 <- derivIPEC(MM, theta=par1, z=0.02, method="Richardson",
            method.args=list(eps=1e-4, d=0.11, 
            zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2))
res1
####################################################################################################


#### Example 2 #####################################################################################
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
resu1 <- derivIPEC( myfun, theta=par2, z=15, method="Richardson", 
            method.args=list(eps=1e-4, d=0.11, 
            zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) )
resu1
####################################################################################################


#### Example 3 #####################################################################################
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

# Define the third case of Mitscherlich equation
MitC <- function(P3, x){
    theta1 <- P3[1]
    beta2  <- P3[2]
    beta3  <- P3[3]
    x1     <- 1
    x2     <- 13
    theta2 <- (beta3 - beta2)/(exp(theta1*x2)-exp(theta1*x1))
    theta3 <- beta2/(1-exp(theta1*(x1-x2))) - beta3/(exp(theta1*(x2-x1))-1)
    theta3 + theta2*exp(theta1*x)
}

ini.val6 <- c(-0.15, 2.52, 1.09)
RES0     <- fitIPEC( MitC, x=x4, y=y4, ini.val=ini.val6, xlim=NULL, ylim=NULL,  
                     fig.opt=TRUE, control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
parC     <- RES0$par
parC
RES1     <- derivIPEC( MitC, theta=parC, z=2, method="Richardson", 
                       method.args=list(eps=1e-4, d=0.11, 
                       zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2) )
RES1
#################################################################################################
}
