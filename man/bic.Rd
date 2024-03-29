\name{bic}
\alias{bic}
\title{
Bayesian Information Criterion (BIC) Calculation Function
}
\description{
Calculates the BIC value(s) of the object(s) obtained from 
using the \code{\link{fitIPEC}} function.
}
\usage{
bic( object, ... )
}
\arguments{
  \item{object}{A fitted model object for which there exists the sample size (\code{sample.size} or \code{n}), 
estimate(s) of model parameter(s) (\code{par}), and residual sum of squares (\code{RSS})}
  \item{...}{Optionally more fitted model objects}
}
\details{

BIC = \emph{p} ln(\emph{n}) - 2 ln(\emph{L}), where \emph{p} represents 
the number of model parameter(s) plus 1 for the error, \emph{n} represents the sample size, 
and ln(\emph{L}) represents the maximum log-likelihood of the estimated model (Spiess and Neumeyer, 2010).

}
\value{
There is a BIC value corresponding to one object, 
and there is a vector of BIC values corresponding to the multiple objects.
}
\references{
Spiess, A-N and Neumeyer, N. (2010) An evaluation of R squared as an inadequate measure 
for nonlinear models in pharmacological and biochemical research: a Monte Carlo approach. 
\emph{BMC Pharmacol.} 10, 6. \doi{10.1186/1471-2210-10-6}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Peter M. Ridland \email{p.ridland@unimelb.edu.au}, 
David A. Ratkowsky \email{d.ratkowsky@utas.edu.au}, Yang Li \email{yangli@fau.edu}.
}
\note{
When there are \code{sample.size} and \code{n} in \code{object} at the same time, the default of 
the sample size is \code{sample.size}, which is superior to \code{n}. 
The BIC gives a higher penalty on the number of model parameters than the AIC.
}
\seealso{
\code{\link{aic}}, \code{\link{AIC}} in package \pkg{stats}, and \code{\link{BIC}} in package \pkg{stats}
}
\examples{
#### Example #####################################################################################
data(leaves)
attach(leaves)
# Choose a geographical population (see Table S1 in Wang et al. [2018] for details)
# Wang, P., Ratkowsky, D.A., Xiao, X., Yu, X., Su, J., Zhang, L. and Shi, P. 
#   (2018) Taylor's power law for leaf bilateral symmetry. Forests 9, 500. doi: 10.3390/f9080500
# 1: AJ; 2: HN; 3: HW; 4: HZ; 5: JD; 
# 6: JS; 7: SC; 8: TC; 9: TT; 10: TX
ind <- 1
L   <- Length[PopuCode == ind]
W   <- Width[PopuCode == ind] 
A   <- Area[PopuCode == ind]

# Define a model y = a*(x1*x2), where a is a parameter to be estimated
propor <- function(theta, x){
    a  <- theta[1]
    x1 <- x[,1]
    x2 <- x[,2]
    a*x1*x2
}

# Define a model y = a*(x1^b)*(x2^c), where a, b and c are parameters to be estimated    
threepar <- function(theta, x){
    a  <- theta[1]
    b  <- theta[2]
    c  <- theta[3]
    x1 <- x[,1]
    x2 <- x[,2]
    a*x1^b*x2^c
}

# Define a model y = a*x^b, where a and b are parameters to be estimated    
twopar <- function(theta, x){
    a  <- theta[1]
    b  <- theta[2]
    a*x^b
}

\donttest{
  A1 <- fitIPEC(propor, x=cbind(L, W), y=A, fig.opt=FALSE,
            ini.val=list(seq(0.1, 1.5, by=0.1)))
  B1 <- curvIPEC(propor, theta=A1$par, x=cbind(L, W), y=A)    
  A2 <- fitIPEC(threepar, x=cbind(L, W), y=A, fig.opt=FALSE,
            ini.val=list(A1$par, seq(0.5, 1.5, by=0.1), seq(0.5, 1.5, by=0.1)))    
  B2 <- curvIPEC(threepar, theta=A2$par, x=cbind(L, W), y=A)
  A3 <- fitIPEC(twopar, x=L, y=A, fig.opt=FALSE,
                ini.val=list(1, seq(0.5, 1.5, by=0.05)))    
  B3 <- curvIPEC(twopar, theta=A3$par, x=L, y=A)
  A4 <- fitIPEC(twopar, x=W, y=A, fig.opt=FALSE,
                ini.val=list(1, seq(0.5, 1.5, by=0.05)))    
  B4 <- curvIPEC(twopar, theta=A4$par, x=W, y=A)
  aic(A1, A2, A3, A4)
  bic(A1, A2, A3, A4)
}
##################################################################################################
}