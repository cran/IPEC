\name{isom}
\alias{isom}
\title{
Data on biochemical oxygen demand
}
\description{
Data on the reaction rate of the catalytic isomerization of \emph{n}-pentane to isopentane 
versus the partial pressures of hydrogen, \emph{n}-pentane, and isopentane.
}
\usage{
data(isom)
}
\details{
There are four columns in the data set:

'y' is the vector of experimental reaction rate (in 1/hr);

'x1' is the vector of partial pressure of hydrogen;

'x2' is the vector of partial pressure of \emph{n}-pentane;

'x3' is the vector of partial pressure of isopentane.
}
\note{
There were errors about the definitions of 'x2' and 'x3' in page 272 in Bates and Watts (1988). 
Here, we redefined them according to the paper of Carr (1960).
}
\references{
Bates, D.M and Watts, D.G. (1988) \emph{Nonlinear Regression Analysis and its Applications}. Wiley, New York.

Carr, N.L. (1960) Kinetics of catalytic isomerization of \emph{n}-pentane. \emph{Ind. Eng. Chem.} 52, 391-396. 
}
\examples{
data(isom)
isom
Y <- isom[,1]
X <- isom[,2:4]
X
Y
}
