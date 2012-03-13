\name{HIsurf}
\alias{HIsurf}

\title{
Describe the joint likelihood surface of ancestry and heterozygosity for a hybrid genotype.
}
\description{
\code{HIsurf} calculates the log likelihood of points on a bivariate grid to describe the joint likelihood surface of ancestry and interclass heterozygosity for a genotype given parental allele frequencies.
}
\usage{
HIsurf(G, P, type, size)
}
\arguments{
  \item{G}{
An individual diploid genotype as matrix \code{G}. If \code{type == "codominant"}, \code{G} must be a two-row matrix with one column for each locus, as in STRUCTURE 

(http://pritch.bsd.uchicago.edu/structure.html). If \code{type == "dominant"}, \code{G} is a vector of 0,1 for absence,presence of the dominant allele. If \code{type == "allele.count"}, \code{G} must be a vector of genotypes coded as 0,1,2 for the number of "j" alleles. That is, genotype 2 is homozygous for allele j, genotype 1 is heterozygous, and genotype 0 has no j alleles.
}
  \item{P}{
Parental allele frequencies. A matrix or data frame with the following columns (order is important!): Locus name, Allele name, P1 allele frequency, P2 allele frequency. For \code{type="dominant"} or \code{type="allele.count"}, there should be one row per locus, giving the frequencies of the dominant or "1" allele. For \code{type="codominant"} there should be a separate row for each allele AND the Allele names should match the data in \code{G}.
}
  \item{type}{
A string representing the data type. The options are \code{"codominant"}, \code{"dominant"}, and \code{"allele.count"}.
}
  \item{size}{
An integer giving the desired number of gridlines in each direction. The function will calculate the likelihood for all \code{(size^2)/2} combinations of S and H that fall within the triangular sample space.
}
}
\details{
Given two ancestral species or parental populations (P1 and P2), the ancestry index (S) is the proportion of an individual's alleles descending from alleles in the P1 population and the interclass heterozygosity (H) is the proportion of an individual's loci that have one allele from each ancestral population (Lynch 1991). The likelihood functions are described in Fitzpatrick (2012).
}
\value{
A \code{size x size} matrix of log likelihoods for all combinations of ancestry (S) and interclass heterozygosity (H). Rows correspond to the \code{size} values of S, and columns the \code{size} values of H. For impossible combinations (\code{H > min(2*S,2-2*S)}), NA's are returned.
}
\references{
Fitzpatrick, B. M. 2008. Hybrid dysfunction: Population genetic and quantitative genetic perspectives. American Naturalist 171:491-198.

Fitzpatrick, B. M. 2012. Estimating ancestry and heterozygosity of hybrids using molecular markers.

Lynch, M. 1991. The genetic interpretation of inbreeding depression and outbreeding depression. Evolution 45:622-629.
}
\author{
Ben Fitzpatrick
}

\seealso{
\code{\link{HIest}} for maximum likelihood estimation of S and H,  \code{\link{HIclass}} for likelihoods of early generation hybrid classes, \code{\link{HItest}} to compare the classification to the maximum likelihood, \code{\link{HILL}} for the basic likelihood function.
}
\examples{
data(Bluestone)
Bluestone <- replace(Bluestone,is.na(Bluestone),-9)
# parental allele frequencies (assumed diagnostic)
BS.P <- data.frame(Locus=names(Bluestone),Allele="BTS",P1=1,P2=0)

# a small surface to view in the console
BS.surf.5 <- HIsurf(Bluestone[21,],BS.P,type="allele.count",size=5)
BS.surf.5 # the maximum likelihood is very near the center (S ~ H ~ 0.5)

# # a more finely sampled surface to visualize with image
# BS.surf <- HIsurf(Bluestone[21,],BS.P,type="allele.count",size=99)

# image(-BS.surf,col=gray(seq(from=0,to=1,length.out=6)),
	# breaks=seq(from=min(-BS.surf,na.rm=TRUE),by=2,length.out=7),
	# cex.axis=1.5,bty="n",xaxs="r",yaxs="r")
# # the breaks option is set so that each level of shading corresponds to 2 log-likelihood
# # units (for one unit increments, set by=1).
# # now make it pretty:
# image(is.na(BS.surf),col="light blue",breaks=c(.5,1),add=TRUE)
# axis(1,labels=FALSE,lwd=2); axis(2,labels=FALSE,lwd=2)
# title(xlab=expression(italic(S)),ylab=expression(italic(H[I])),cex.lab=1.5)
# lines(c(0,.5,1,0),c(0,1,0,0),lty=2,lwd=2)

## The function is currently defined as
function(G,P,type,size){
	S <- H <- seq(from=0,to=1,length.out=size)
	surf <- matrix(NA,ncol=size,nrow=size)
	for(i in 1:size){
		for(j in 1:size){
			surf[i,j] <- HILL(c(S[i],H[j]),G,P,type)
			}}
	for(i in 1:size){
		for(j in 1:size){
			if(H[j]>min(2*S[i],2-2*S[i])) surf[i,j] <- NA
		}}
surf
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
