\name{bootIPEC}
\alias{bootIPEC}
\title{Bootstrap Function for Nonlinear Regression}
\description{
Generates the density distributions, standard deviations, confidence intervals, 
covariance matrices and correlation matrices of parameters 
based on bootstrap replications.
}
\usage{
bootIPEC( expr, x, y, ini.val, weights = NULL, control = list(), 
          nboot = 200, CI = 0.95, fig.opt = TRUE, fold = 3.5, 
          unique.num = 2, prog.opt = TRUE )
}
\arguments{
  \item{expr}{A given parametric model}
  \item{x}{A vector or matrix of observations of independent variable(s)}
  \item{y}{A vector of observations of response variable}
  \item{ini.val}{A vector or list of initial values of model parameters}
  \item{weights}{An optional vector of weights to be used in the fitting process. 
\code{weights} should be NULL or a numeric vector. If non-NULL, weighted least squares is used with weights \code{weights}; 
otherwise ordinary least squares is used.}
  \item{control}{A list of control parameters for using the \code{\link{optim}} 
function in package \pkg{stats}}
  \item{nboot}{The number of bootstrap replications}
  \item{CI}{The confidence level(s) of the required interval(s)}
  \item{fig.opt}{An option of drawing figures of the distributions of bootstrap values of 
parameters and figures of pairwise comparisons of bootstrap values}
  \item{fold}{A parameter removing the extreme bootstrap values of parameters}
  \item{unique.num}{The least number of sampled non-overlapping data points 
for carrying out a bootstrap nonlinear regression}
  \item{prog.opt}{An option of showing the running progress of bootstrap}
}
\details{
\code{ini.val} can be a vector or a list that has saved initial values for model parameters, 

e.g. \code{y = beta0 + beta1 * x + beta2 * x^2}, 

\code{ini.val = list(beta0=seq(5, 15, len=2), beta1=seq(0.1, 1, len=9), 
beta2=seq(0.01, 0.05, len=5))}, which is similar to the usage of the 
input argument of \code{start} of \code{\link{nls}} in package \pkg{stats}. 

In the \code{weights} argument option, the default is \code{weights = NULL}. 
In that case, ordinary least squares is used. 
The residual sum of squares (RSS) between the observed and predicted \eqn{y} values 
is minimized to estimate a model's parameters, i.e., 
\deqn{\mbox{RSS} = \sum_{i=1}^{n}\left(y_i-\hat{y}_i\right)^{2}}
where \eqn{y_i} and \eqn{\hat{y}_i} represent the observed and predicted \eqn{y} values, respectively; 
and \eqn{n} represents the sample size. If \code{weights} is a numeric vector, 
the weighted residual sum of squares is minimized, i.e.,   
\deqn{\mbox{RSS} = \sum_{i=1}^{n}w_i\left(y_i-\hat{y}_i\right)^{2}} 
where \eqn{w_i} is the \eqn{i} elements of \code{weights}.

\code{CI} determines the width of confidence intervals. 

\code{fold} is used to delete the data whose differences from the median exceed 
a certain \code{fold} of the difference between 3/4 and 1/4 quantiles of the 
bootstrap values of a model parameter. 

The default of \code{unique.num} is 2. That is, at least two non-overlapping data 
points randomly sampled from \code{(x,y)} are needed for carrying out 
a bootstrap nonlinear regression.
}
\value{
  \item{M}{The matrix saving the fitted results of all \code{nboot} bootstrap 
values of model parameters and goodness of fit}
  \item{perc.ci.mat}{The matrix saving the estimate, standard deviation, median, mean, 
and the calculated lower and upper limits of confidence interval based on 
the bootstrap percentile method}
  \item{bca.ci.mat}{The matrix saving the estimate, standard deviation, median, mean, 
and the calculated lower and upper limits of confidence interval based on 
the bootstrap \eqn{BC_a}{BC[a]} method}
  \item{covar.mat}{The covariance matrix of parameters based on the bootstrap 
values when \code{nboot > 1}}
  \item{cor.mat}{The correlation matrix of parameters based on the bootstrap 
values when \code{nboot > 1}}
}
\references{
Efron, B. and Tibshirani, R.J. (1993) \emph{An Introduction to the Bootstrap}. 
Chapman and Hall (CRC), New York. \doi{10.2307/2532810}

Sandhu, H.S., Shi, P., Kuang, X., Xue, F. and Ge, F. (2011) Applications of the bootstrap to 
insect physiology. \emph{Fla. Entomol.} 94, 1036\eqn{-}1041. \doi{10.1653/024.094.0442}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Peter M. Ridland \email{p.ridland@unimelb.edu.au}, 
David A. Ratkowsky \email{d.ratkowsky@utas.edu.au}, Yang Li \email{yangli@fau.edu}.
}
\note{
To obtain reliable confidence intervals of model parameters, more than \strong{2000} 
bootstrap replications are recommended; whereas to obtain a reliable standard deviation of the estimate 
of a parameter, more than \strong{30} bootstrap replications are sufficient (Efron and Tibshirani 1993). 
\code{bca.ci.mat} is recommended to show better confidence intervals of parameters than 
those in \code{perc.ci.mat}. 

The outputs of model parameters will all be represented by \eqn{\theta_i}{theta[i]}, \code{i} from 1 
to \code{p}, where \code{p} represents the number of model parameters. The letters of model 
parameters defined by users such as \eqn{\beta_i}{beta[i]} will be automatically replaced by \eqn{\theta_i}{theta[i]}.
}
\seealso{
\code{\link{fitIPEC}}
}
\examples{
#### Example 1 #################################################################################
graphics.off()
# The velocity of the reaction (counts/min^2) under different substrate concentrations 
#   in parts per million (ppm) (Page 269 of Bates and Watts 1988)

x1 <- c(0.02, 0.02, 0.06, 0.06, 0.11, 0.11, 0.22, 0.22, 0.56, 0.56, 1.10, 1.10)
y1 <- c(76, 47, 97, 107, 123, 139, 159, 152, 191, 201, 207, 200)

# Define the Michaelis-Menten (MM) model
MM <- function(theta, x){
    theta[1]*x / ( theta[2] + x )    
}

\donttest{
  set.seed(123)
  res4 <- bootIPEC( MM, x=x1, y=y1, ini.val=c(200, 0.05), 
                    control=list(reltol=1e-20, maxit=40000), nboot=2000, CI=0.95, 
                    fig.opt=TRUE )
  res4
  set.seed(NULL)
}
#################################################################################################


#### Example 2 ##################################################################################
graphics.off()
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
  if (Tmin >= Tmax | ropt <= 0 | Topt <= Tmin | Topt >= Tmax) 
    temp <- Inf
  if (Tmax > Tmin & ropt > 0 & Topt > Tmin & Topt < Tmax){
    temp <- sqrt( ropt*(x-Tmax)*(x-Tmin)^2/((Topt-Tmin)*((Topt-Tmin
      )*(x-Topt)-(Topt-Tmax)*(Topt+Tmin-2*x))) )  
  }
  return( temp )
}

myfun <- sqrt.LRF
\donttest{
  set.seed(123)
  resu4 <- bootIPEC( myfun, x=x2, y=y2, ini.val=ini.val1, 
                     nboot=2000, CI=0.95, fig.opt=TRUE )
  resu4
  set.seed(NULL)
}
#################################################################################################


#### Example 3 ##################################################################################
graphics.off()
# Height growth data of four species of bamboo (Gramineae: Bambusoideae)
# Reference(s):
# Shi, P., Fan, M., Ratkowsky, D.A., Huang, J., Wu, H., Chen, L., Fang, S. and  
#     Zhang, C. (2017) Comparison of two ontogenetic growth equations for animals and plants. 
#     Ecol. Model. 349, 1-10.

data(shoots)
# Choose a species
# 1: Phyllostachys iridescens; 2: Phyllostachys mannii; 
# 3: Pleioblastus maculatus; 4: Sinobambusa tootsik.
# 'x3' is the vector of the observation times from a specific starting time of growth
# 'y3' is the vector of the aboveground height values of bamboo shoots at 'x3' 

ind <- 4
x3  <- shoots$x[shoots$Code == ind]
y3  <- shoots$y[shoots$Code == ind] 

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
  return(ropt*(x-tmin)*(x-2*tmax+topt)/(topt+tmin-
         2*tmax)*( (x-tmin)/(topt-tmin) )^((topt-tmin)/(tmax-topt)))   
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
  return(ropt*(x-tmin)*(x-2*tmax+topt)/(topt+tmin-
         2*tmax)*((x-tmin)/(topt-tmin) )^((topt-tmin)/(tmax-topt)))   
}

# For the original beta sigmoid model
ini.val2 <- c(40, 30, 5, 50)
xlab2    <- "Time (d)"
ylab2    <- "Height (cm)"

\donttest{
  set.seed(123)
  re4 <- bootIPEC( bsm, x=x3, y=y3, ini.val=ini.val2,    
                   control=list(trace=FALSE, reltol=1e-20, maxit=50000),
                   nboot=2000, CI=0.95, fig.opt=TRUE, fold=10 )
  re4
  set.seed(NULL)
}

# For the simplified beta sigmoid model (in comparison with the original beta sigmoid model)
ini.val7 <- c(40, 30, 50)

\donttest{
  set.seed(123)
  RESU4 <- bootIPEC( simp.bsm, x=x3, y=y3, ini.val=ini.val7,   
                     control=list(trace=FALSE, reltol=1e-20, maxit=50000),
                     nboot=2000, CI=0.95, fig.opt=TRUE, fold=10 )
  RESU4
  set.seed(NULL)
}
#################################################################################################


#### Example 4 ##################################################################################
graphics.off()
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
y4 <- c( 3.183, 3.059, 2.871, 2.622, 2.541, 2.184, 
         2.110, 2.075, 2.018, 1.903, 1.770, 1.762, 1.550 )

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

\donttest{
  set.seed(123)
  ini.val3 <- c(-0.1, 2.5, 1.0)
  r4       <- bootIPEC( MitA, x=x4, y=y4, ini.val=ini.val3,    
                        nboot=2000, CI=0.95, fig.opt=TRUE )
  r4

  ini.val4 <- c(exp(-0.1), log(2.5), 1)
  R4       <- bootIPEC( MitB, x=x4, y=y4, ini.val=ini.val4, 
                        nboot=2000, CI=0.95, fig.opt=TRUE )
  R4

  # ini.val6 <- c(-0.15, 2.52, 1.09)
  iv.list2 <- list(seq(-2, -0.05, len=5), seq(1, 4, len=8), seq(0.05, 3, by=0.5))
  RES0 <- fitIPEC( MitC, x=x4, y=y4, ini.val=iv.list2,    
                   control=list(trace=FALSE, reltol=1e-10, maxit=5000) )
  RES0$par
  RES4 <- bootIPEC( MitC, x=x4, y=y4, ini.val=iv.list2, 
                    control=list(trace=FALSE, reltol=1e-10, maxit=5000), 
                    nboot=5000, CI=0.95, fig.opt=TRUE, fold=3.5, unique.num=2 )
  RES4
  set.seed(NULL)
}
#################################################################################################
}
