\name{fitIPEC}
\alias{fitIPEC}
\title{
Nonlinear Fitting Function
}
\description{
Estimates the parameters of a given parametric model 
using the \code{\link{optim}} function in package \pkg{stats}.
}
\usage{
fitIPEC( expr, x, y, ini.val, target.fun = "RSS", control = list(), 
         fig.opt = TRUE, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL )
}
\arguments{
  \item{expr}{A given parametric model}
  \item{x}{Vector or matrix of observations of independent variable(s)}
  \item{y}{Vector of observations of response variable}
  \item{ini.val}{A vector or list of initial values of model parameters}
  \item{target.fun}{Objective function, \dQuote{RSS} 
and \dQuote{chi.sq}, for the optimization}
  \item{control}{A list of control parameters for using 
the \code{\link{optim}} function in package \pkg{stats}}
  \item{fig.opt}{An option to determine whether to draw the fitted curve}
  \item{xlim}{The shown range of the \emph{x}-axis}
  \item{ylim}{The shown range of the \emph{y}-axis}
  \item{xlab}{The label of the \emph{x}-axis}
  \item{ylab}{The label of the \emph{y}-axis}
}
\details{
The Nelder-Mead algorithm is the default in the \code{\link{optim}} function in package \pkg{stats}. 
The user can accurately estimate the model parameters by setting smaller relative convergence 
tolerance and larger maximum number of iterations in the input argument of \code{control}, 

e.g. \code{control=list(trace=FALSE, reltol=1e-20, maxit=50000)}, 

at the expense of the running speed.

\code{ini.val} can be a vector or a list that has saved initial values for model parameters, 

e.g. \code{y = beta0 + beta1 * x + beta2 * x^2}, 

\code{ini.val = list(beta0=seq(5, 15, len=2), beta1=seq(0.1, 1, len=9), 
beta2=seq(0.01, 0.05, len=5))}, which is similar to the usage of the 
input argument of \code{start} of \code{\link{nls}} in package \pkg{stats}.
}
\value{
  \item{expr}{The formula used}
  \item{par}{Vector of estimates of parameters}
  \item{RSS}{Calculated residual sum of squares}
  \item{chi.sq}{Calculated \eqn{\chi^2}}
  \item{R.sq}{Calculated \eqn{R^2}}
}
\references{
Nelder, J.A. and Mead, R. (1965) A simplex algorithm for function minimization. 
\emph{Comput. J.} 7, 308-313.
}
\author{
Peijian Shi, Peter Ridland, David A. Ratkowsky
}
\note{
This function can be applicable to a nonlinear parametric model with 
a single independent variable or with multiple independent variables. 

\code{R.sq} is only used to help users intuitively judge whether the fitted curve seriously 
deviates from the actual observations. However, it should NOT be used to decide which of 
several competing models is the most appropriate (Pages 44-45 in Ratkowsky 1990). 
\code{RSS} and curvatures are among the suitable candidates to answer such a question.
}
\seealso{
\code{\link{bootIPEC}}, \code{\link{optim}} in package \pkg{stats}
}
\examples{
#### Example 1 ###################################################################################
graphics.off()
# The velocity of the reaction (counts/min^2) under different substrate concentrations 
#   in parts per million (ppm) (Page 269 of Bates and Watts 1988)

x1 <- c(0.02, 0.02, 0.06, 0.06, 0.11, 0.11, 0.22, 0.22, 0.56, 0.56, 1.10, 1.10)
y1 <- c(76, 47, 97, 107, 123, 139, 159, 152, 191, 201, 207, 200)

# Define the Michaelis-Menten model
MM <- function(theta, x){
    theta[1]*x / ( theta[2] + x )    
}

res0 <- fitIPEC(MM, x=x1, y=y1, ini.val=c(200, 0.05), 
                xlim=c(0, 1.5), ylim=c(0, 250), fig.opt=TRUE)
par1 <- res0$par
par1
res0

# The input names of parameters will not affect the fitted results.
# We can use other names to replace theta1 and theta2.   
iv.list1 <- list( theta1=seq(100, 300, by=50), theta2=seq(10, 100, by=10) )  
result0  <- fitIPEC( MM, x=x1, y=y1, ini.val=iv.list1, xlim=c(0, 1.5), ylim=c(0, 250), 
                     fig.opt=FALSE, control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
param1   <- result0$par
param1
##################################################################################################


\dontrun{
#### Example 2 ###################################################################################
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
ylab1 <- expression( paste("Developmental rate"^{1/2}, 
                     " (", d^{"-1"}, ")", sep="") )
resu0 <- fitIPEC( myfun, x=x2, y=y2, ini.val=ini.val1, xlim=NULL, 
                  ylim=NULL, xlab=xlab1, ylab=ylab1, fig.opt=TRUE, 
				  control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
par2  <- resu0$par
par2
resu0
##################################################################################################
}


#### Example 3 ###################################################################################
graphics.off()
# Height growth data of four species of bamboos (Gramineae: Bambusoideae)
# Reference(s):
# Shi, P.-J., Fan, M.-L., Ratkowsky, D.A., Huang, J.-G., Wu, H.-I, Chen, L.,   
#     Fang, S.-Y. and Zhang, C.-X. (2017) Comparison of two ontogenetic  
#     growth equations for animals and plants. Ecol. Model. 349, 1-10.

data(shoots)
attach(shoots)
# Choose a species
# 1: Phyllostachys iridescens; 2: Phyllostachys mannii; 
# 3: Sinobambusa tootsik; 4: Pleioblastus maculatus
# 'x3' is the vector of the observation times from a specific starting time of growth
# 'y3' is the vector of the aboveground height values of bamboo shoots at 'x3' 
ind <- 2
x3  <- time[code == ind]
y3  <- height[code == ind] 

# Define the beta sigmoid model (bsm)
bsm <- function(P, x){
  P  <- cbind(P)
  if(length(P) !=4 ) {stop(" The number of parameters should be 4!")}
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

ini.val2 <- c(40, 30, 5, 50)
xlab2    <- "Time (d)"
ylab2    <- "Height (cm)"

re0  <- fitIPEC( bsm, x=x3, y=y3, ini.val=ini.val2, 
                 xlim=NULL, ylim=NULL, xlab=xlab2, ylab=ylab2, 
                 fig.opt=TRUE, control=list(trace=FALSE, reltol=1e-20, maxit=50000) )
par3 <- re0$par
par3
##################################################################################################


#### Example 4 ###################################################################################
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
##################################################################################################
}
