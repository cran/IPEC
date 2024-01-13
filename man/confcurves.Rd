\name{confcurves}
\alias{confcurves}
\title{
Wald Confidence Curves and the Likelihood Confidence Curves
}
\description{
Calculates the Wald confidence curves and the likelihood confidence curves of model parameters.
}
\usage{
confcurves( expr, x, y, ini.val, weights = NULL, control=list(), 
            fig.opt = TRUE, fold = 5, np = 20, alpha = seq(1, 0.001, by=-0.001), 
            show.CI = NULL, method = "Richardson", method.args = 
            list(eps = 1e-04, d = 0.11, zero.tol = sqrt(.Machine$double.eps/7e-07), 
            r = 6, v = 2, show.details = FALSE), side = NULL )
}
\arguments{
  \item{expr}{A given parametric model}
  \item{x}{A vector or matrix of observations of independent variable(s)}
  \item{y}{A vector of observations of response variable}
  \item{ini.val}{A vector or list of initial values of model parameters}
  \item{weights}{An optional vector of weights to be used in the fitting process. 
\code{weights} should be NULL or a numeric vector. If non-NULL, weighted least squares is used with weights \code{weights}; 
otherwise ordinary least squares is used.}
  \item{control}{A list of control parameters for using 
the \code{\link{optim}} function in package \pkg{stats}}
  \item{fig.opt}{An option to determine whether to draw the confidence curves of each parameter}
  \item{fold}{The fold of \eqn{ \mbox{SE}\left(\hat{\theta_{i}}\right) } for controlling the width of the confidence interval of \eqn{\hat{\theta_{i}}} that represents the estimate of the \eqn{i}th parameter}
  \item{np}{The number of data points for forming the lower or upper bounds of a likelihood confidence interval of \eqn{ \hat{\theta_{i}} }, which controlls the step size (i.e., \eqn{\delta}) in the \eqn{y} coordinates of the likelihood confidence curves}
  \item{alpha}{The significance level(s) for calculating the \eqn{x} coordinate(s) of the \eqn{(1-\alpha)100\%} Wald confidence curves, which equals to 
\eqn{ t_{\alpha/2}\left(n-p\right) }}
  \item{show.CI}{The \eqn{ t_{\alpha/2}\left(n-p\right) } value(s) associated with the confidence level(s) of each parameter to be showed, i.e., \code{c(0.80, 0.90, 0.95, 0.99)}}
  \item{method}{It is the same as the input argument of \code{method} of the \code{\link[numDeriv]{hessian}} function in package \pkg{numDeriv}}
  \item{method.args}{It is the same as the input argument of \code{method.args} of the \code{\link[numDeriv]{hessian}} function in package \pkg{numDeriv}}
  \item{side}{It is the same as the input argument of \code{side} of the \code{\link[numDeriv]{jacobian}} function in package \pkg{numDeriv}}
}
\details{
For the \eqn{(1-\alpha)100\%} Wald confidence curves, the corresponding \eqn{x} and \eqn{y} coordinates are:
\deqn{ x = t_{\alpha/2}\left(n-p\right), }
and
\deqn{ y = \hat{\theta_{i}} \pm t_{\alpha/2}\left(n-p\right)\,\mbox{SE}\left(\hat{\theta_{i}}\right), }
where \eqn{n} denotes the number of the observations, \eqn{p} denotes the number of model parameters, and 
\eqn{\mbox{SE}\left(\hat{\theta_{i}}\right)} denotes the standard error of the \eqn{i}th model parameter's estimate. 

\eqn{\quad} For the likelihood confidence curves (Cook and Weisberg, 1990), the corresponding \eqn{x} and \eqn{y} coordinates are:
\deqn{ x = \sqrt{\frac{\mbox{RSS}\left(\hat{\theta}^{\,\left(-i\right)}\right)-\mbox{RSS}\left(\hat{\theta}\right)}{\mbox{RSS}\left(\hat{\theta}\right)/(n-p)}}, }
where \eqn{\mbox{RSS}\left(\hat{\theta}\right)} represents the residual sum of squares for fitting the model with all model parameters; 
\eqn{\mbox{RSS}\left(\hat{\theta}^{\,\left(-i\right)}\right)} represents the residual sum of squares for fitting the model with the \eqn{i}th model parameter 
being fixed to be \eqn{\hat{\theta_{i}} \pm k\,\delta}. Here, \eqn{k} denotes the \eqn{k}th iteration, and \eqn{\delta} denotes the step size, 
which equals  
\deqn{ \delta = \frac{\hat{\theta_{i}} \pm \mbox{fold}\,\mbox{SE}\left(\hat{\theta_{i}}\right)}{\mbox{np}}. }
\deqn{ y = \hat{\theta_{i}} \pm k\,\delta. }
Here, \code{fold} and \code{np} are defined by the user in the arguments. 

\eqn{\quad} For other arguments, please see the \code{\link{fitIPEC}} and \code{\link{parinfo}} functions for details.
}
\value{
  \item{partab}{The estimates, standard errors and confidence intervals of model parameters; also see the \code{\link{parinfo}} function}
  \item{parlist}{A list for the estimate, Wald interval curves and likelihood interval curves of each model parameter.}
}
\references{
Cook, R.D. and Weisberg, S. (1990) Confidence curves in nonlinear regression. 
\emph{J. Am. Statist. Assoc.} 82, 221\eqn{-}230. \doi{10.1080/01621459.1990.10476233}

Nelder, J.A. and Mead, R. (1965) A simplex method for function minimization. 
\emph{Comput. J.} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Ratkowsky, D.A. (1990) \emph{Handbook of Nonlinear Regression Models}, Marcel Dekker, New York.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Peter M. Ridland \email{p.ridland@unimelb.edu.au}, 
David A. Ratkowsky \email{d.ratkowsky@utas.edu.au}, Yang Li \email{yangli@fau.edu}.
}
\note{
In the value of \code{parlist}, there are the estimate (\code{pari}), the Wald interval curves 
(\code{WaldCI}), and the likelihood interval curves (\code{lhCI}) of the \eqn{i}th model parameter. 
In \code{WaldCI}, there are three columns. The first column, \code{tc}, 
represents \eqn{ t_{\alpha/2}\left(n-p\right) }, the second and third columns, \code{LCI} and \code{UCI}, 
represent the lower and upper bounds of the \eqn{(1-\alpha)100\%} Wald confidence intervals, respectively. 
In \code{lhCI}, there are six columns. The first and second columns, \code{x.lower} and \code{lhLCI}, 
represent the lower bounds of the likelihood confidence intervals the and corresponding \eqn{x} values, 
respectively; the third and fourth columns, \code{x.upper} and \code{lhUCI}, represent the upper 
bounds of the likelihood confidence intervals and the corresponding \eqn{x} values, respectively; 
the fifth and sixth columns, \code{RSS.lower} and \code{RSS.upper}, represent 
the values of the residual sum of squares of the lower bounds 
and those of the upper bounds, respectively. Please see Cook and Weisberg (1990) for details.

}
\seealso{
\code{\link{parinfo}}, \code{\link{fitIPEC}}, \code{\link{optim}} in package \pkg{stats}
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
MitA <- function(P, x){    
    P[3] + P[2]*exp(P[1]*x)
}

# Define the second case of Mitscherlich equation
MitB <- function(P2, x){
    if(P2[3] <= 0)
        temp <- mean(y4)
    if(P2[3] > 0)
       temp <- log(P2[3]) + exp(P2[2] + P2[1]*x)
    return( temp )
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
RESU1    <- confcurves( MitA, x=x4, y=y4, ini.val=ini.val3, fig.opt = TRUE, 
                        fold=5, np=20, alpha=seq(1, 0.001, by=-0.001), 
                        show.CI=c(0.8, 0.9, 0.95, 0.99) )

ini.val4 <- c(-0.10, 0.90, 2.5)
RESU2    <- confcurves( MitB, x=x4, y=y4, ini.val=ini.val4, fig.opt = TRUE, 
                        fold=5, np=200, alpha=seq(1, 0.001, by=-0.001), 
                        show.CI=c(0.8, 0.9, 0.95, 0.99) )

ini.val6 <- c(-0.15, 2.5, 1)
RESU3    <- confcurves( MitC, x=x4, y=y4, ini.val=ini.val6, fig.opt = TRUE, 
                        fold=5, np=20, alpha=seq(1, 0.001, by=-0.001), 
                        show.CI=c(0.8, 0.9, 0.95, 0.99) )
##################################################################################################

graphics.off()
}
