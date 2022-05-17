\name{shoots}
\alias{shoots}
\title{
Height Growth Data of Bamboo Shoots
}
\description{
The height growth data of four species of bamboo at the Nanjing Forestry University campus in 2016.
}
\usage{
data(shoots)
}
\details{
In the data set, there are four variables: \code{Code}, \code{LatinName}, \code{x} and \code{y}. 
\code{Code} is used to save the number codes of different bamboo species; 
\code{LatinName} is used to save the Latin names of different bamboo species; 
\code{x} is used to save the investigation times (d) from a specific starting time of growth, 
and every bamboo has a different starting time of growth;  
\code{y} is used to save the measured aboveground height values (cm).
 
\code{Code = 1} represents \emph{Phyllostachys iridescens}, and the starting time (namely \code{x = 0}) was defined as 12:00, 3rd April, 2016; 
 
\code{Code = 2} represents \emph{Phyllostachys mannii}, and the starting time (namely \code{x = 0}) was defined as 12:00, 4th April, 2016; 

\code{Code = 3} represents \emph{Pleioblastus maculatus}, and the starting time (namely \code{x = 0}) was defined as 12:00, 29th April, 2016; 

\code{Code = 4} represents \emph{Sinobambusa tootsik}, and the starting time (namely \code{x = 0}) was defined as 12:00, 18th April, 2016.
}
\references{
Shi, P., Fan, M., Ratkowsky, D.A., Huang, J., Wu, H., Chen, L., Fang, S. and 
Zhang, C. (2017) Comparison of two ontogenetic growth equations for animals and plants. 
\emph{Ecol. Model.} 349, 1\eqn{-}10. \doi{10.1016/j.ecolmodel.2017.01.012}
}
\examples{
data(shoots)
# Choose a species
# 1: Phyllostachys iridescens; 2: Phyllostachys mannii; 
# 3: Pleioblastus maculatus; 4: Sinobambusa tootsik.
ind <- 4
x3  <- shoots$x[shoots$Code == ind]
y3  <- shoots$y[shoots$Code == ind] 
dev.new()
plot(x3, y3, cex=1.5, cex.lab=1.5, cex.axis=1.5, xlab="Time (d)", ylab="Height (cm)")
}
