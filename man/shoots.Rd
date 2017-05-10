\name{shoots}
\alias{shoots}
\title{
Height Growth Data of Bamboo Shoots
}
\description{
The height growth data of four species of bamboos in the Nanjing Forestry University campus in 2016.
}
\usage{
data(shoots)
}
\details{
In the data set, there are four variables: \code{code}, \code{species}, \code{time} and \code{height}. 
\code{code} is used to save the number codes of different bamboo species; 
\code{species} is used to save the Latin names of different bamboo species; 
\code{time} is used to save the observation times (d) from a specific starting time of growth, 
and every bamboo has a different starting time of growth;  
\code{height} is used to save the measured aboveground height values (cm).
 
\code{code = 1} represents \emph{Phyllostachys iridescens}, and the starting time (namely \code{time = 0}) was defined as 12:00, 3rd April, 2016; 
 
\code{code = 2} represents \emph{Phyllostachys mannii}, and the starting time (namely \code{time = 0}) was defined as 12:00, 4th April, 2016; 

\code{code = 3} represents \emph{Sinobambusa tootsik}, and the starting time (namely \code{time = 0}) was defined as 12:00, 18th April, 2016; 

\code{code = 4} represents \emph{Pleioblastus maculatus}, and the starting time (namely \code{time = 0}) was defined as 12:00, 29th April, 2016.
}
\references{
Shi, P.-J., Fan, M.-L., Ratkowsky, D.A., Huang, J.-G., Wu, H.-I, Chen, L., Fang, S.-Y. and 
Zhang, C.-X. (2017) Comparison of two ontogenetic growth equations for animals and plants. \emph{Ecol. Model.} 349, 1-10.
}
\examples{
data(shoots)
attach(shoots)
# Choose a species
# 1: Phyllostachys iridescens; 2: Phyllostachys mannii; 
# 3: Sinobambusa tootsik; 4: Pleioblastus maculatus
ind <- 1
x3  <- time[code == ind]
y3  <- height[code == ind] 
dev.new()
par(mar=c(5,5,2,2))
plot(x3, y3, cex=1.5, cex.lab=1.5, cex.axis=1.5, xlab="Time (d)", ylab="Height (cm)")
}
