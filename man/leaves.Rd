\name{leaves}
\alias{leaves}
\title{
Leaf Data of \emph{Parrotia subaequalis} (Hamamelidaceae)
}
\description{
The data consist of the area, length and width of the leaves of 10 geographical populations of \emph{P. subaequalis} 
collected in Southern China from July to September, 2016.
}
\usage{
data(leaves)
}
\details{
In the data set, there are five variables: \code{PopuCode}, \code{LeafCode}, \code{Length}, \code{Width} and \code{Area}. 
\code{PopuCode} is used to save the number codes of different geographical populations; 
\code{LeafCode} is used to save the number codes of the leaves of a geographical population; 
\code{Length} is used to save the scanned leaf length data (cm);
\code{Width} is used to save the scanned leaf width data (cm);   
\code{Area} is used to save the scanned leaf area data (cm squared). 
}
\references{
Wang, P., Ratkowsky, D.A., Xiao, X., Yu, X., Su, J., Zhang, L. and Shi, P. 
(2018) Taylor's power law for leaf bilateral symmetry. \emph{Forests} 9, 500. doi: 10.3390/f9080500 
}
\examples{
data(leaves)
attach(leaves)
# Choose a geographical population (see Table S1 in Wang et al. [2018] for details)
# 1: AJ; 2: HN; 3: HW; 4: HZ; 5: JD; 
# 6: JS; 7: SC; 8: TC; 9: TT; 10: TX
ind <- 1
L   <- Length[PopuCode == ind]
W   <- Width[PopuCode == ind] 
A   <- Area[PopuCode == ind]
x   <- L*W
fit <- lm(A ~ x-1)
summary(fit)

# Show the leaf areas of the 10 geographical populations
dev.new()
par(mar=c(5,6,2,2))
boxplot(Area~PopuCode, cex=1.5, cex.lab=1.5, cex.axis=1.5, 
        col="grey70", xlab=expression(bold("Population code")), 
        ylab=expression(bold(paste("Leaf area (cm", ""^{"2"}, ")", sep=""))),
        ylim=c(0, 50), xaxs="i", yaxs="i", las=1)
}
