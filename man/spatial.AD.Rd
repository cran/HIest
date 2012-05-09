\name{spatial.AD}
\alias{spatial.AD}
\title{
Simulate admixture in continuous, 2-dimensional space
}
\description{
Simulate admixture dynamics with or without selection on a few loci. Although space matters for local density-dependent competition, mating and dispersal are random (uniform) with respect to space.
}
\usage{
spatial.AD(minX ,minY, maxX, maxY, XY, Genotypes, 
						beta=0,sel=0, mid=0,  
						h=0, DM = matrix(0,ncol=3,nrow=3),
						sigmac, R, M, gens,
						immigrants=FALSE,plotgrowth=FALSE,m=0.10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{minX, minY, maxX, maxY}{Limits of the model space in x and y dimensions.}
  \item{XY}{Matrix of initial x,y coordinates of individual organinisms.}
  \item{Genotypes}{Matrix of genotypes of initial organisms. Each genotype should be coded as 0, 0.5, or 1 for the frequency of alleles derived from one parental population. Rows are organisms, columns are unlinked loci. The first four loci can cause fitness variation.}
  \item{beta}{Steepness of an environmental gradient affecting the first locus.}
  \item{sel}{Strength of environmental selection affecting the first locus.}
  \item{mid}{Midpoint of the environmental gradient affecting the first locus.}
  \item{h}{Selection on heterozygotes at the second locus.}
  \item{DM}{Matrix of 2-locus fitness values for the 3rd and 4th loci (see \code{details}).}
  \item{sigmac}{Local competition parameter: Standard deviation of Gaussian competition function.}
  \item{R}{Instantaneous growth rate of the Beverton-Holt model.}
  \item{M}{Determines the local carrying capacity of the Beverton-Holt (K = (R-1)*M).}
  \item{gens}{Number of generations to simulate.}
  \item{immigrants}{If FALSE, the model space is closed to immigration and all boundaries are reflecting. If TRUE, the model is open to immigrants from pure parental populations at each edge of the x-dimension. If TRUE, \code{m} controls the edge dynamic (see below).}
  \item{plotgrowth}{If TRUE, the population size at each generation will be plotted.}
  \item{m}{Immigration parameter. If \code{immigrants=TRUE}, any individual within m/2 of each edge in the x-dimension will be replaced by a pure parental genotype.}
}
\details{
For the DM incompatibility, the matrix of fitnesses is 3x3, with rows corresponding to the first DM locus and columns correspinding to the second DM locus. Entries are W[i,j], where i and j index genotypes 0, 1, and 2 at the first and second locus, respectively. See example.
}
\value{
A list with
\item{XY }{The x,y coordinates of the diploid individuals in the final generation.}
\item{Genotypes }{The genotypes of the diploid individuals (rows) in the final generation.}
\item{mothers}{The genotypes of the successful mothers in the next-to-last generation (roughly, an "after selection" sample from that generation).}
}

\references{
Fitzpatrick, B. M. Alternative forms for genomic clines. In review
}
\author{
Benjamin M. Fitzpatrick}

\seealso{
See \code{\link{spatial.HZ}} for a version with limited dispersal and mating distances. The simulated data can be analyzed with \code{\link{Cline.fit}}, but the genotypes must be multiplied by 2.}
\examples{
\dontrun{	
# define space:
minX <- minY <- -3
maxX <- maxY <- 3
# 100 individuals randomly placed:
XY <- cbind(runif(100,minX,maxX),runif(100,minY,maxY)) 
# simulate secondary contact by sorting aling the x dimension and assigning parental genotypes on each side of the centerline:
XY <- XY[order(XY[,1]),]
Genotypes <- rbind(matrix(0,nrow=sum(XY[,1]<=0),ncol=10),matrix(1,nrow=sum(XY[,1]>0),ncol=10))
# competition parameters:
 sigmac <- 0.2; R <- 1.75; M <- 5

# selection, inlcuding heterozygote disadvantage at locus 2 and a DM incompatibility between 3 and 4:
beta <- 0
sel <- 0
mid <- 0
h <- 0.4
DM <- rbind(
c(0,0.2,0.4),
c(0,0.0,0.2),
c(0,0.0,0.0))

# simulate 10 generations, open to immigration:
G10 <- spatial.AD(minX,minY,maxX,maxY,XY,Genotypes,beta,sel,mid,h,DM,sigmac,R,M,gens=10,immigrants=TRUE)

}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
