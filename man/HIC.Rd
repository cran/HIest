\name{HIC}
\alias{HIC}
\title{
Closed form maximum-likelihood estimates of ancestry and heterozygosity for diagnostic markers
}
\description{
For hybrid individuals genotyped with diagnostic markers (one allele fixed in each parental lineage), the ancestry index and interclass heterozygosity are calculated using closed form estimators using allele counts.
}
\usage{
HIC(G)
}
\arguments{
  \item{G}{
Data matrix of individual genotypes (individuals in rows, markers in columns), coded as 0, 1, or 2 for the number of alleles inherited from parental lineage 1. 
}
}
\details{
Data must be coded as allele counts and markers must be assumed diagnostic. The MLE of the ancestry index is \eqn{S=\sum(x_i)/(2n)}, where \eqn{x_i} is the allele count for locus \eqn{i}. The MLE of interclass heterozygosity is simply the observed fraction of markers that are heterozygous.
}
\value{
A matrix with three named columns is returned.
  \item{S }{The ancestry index for each individual}
  \item{H }{The interclass heterozygosity for each individual}
  \item{logLik }{The log-likelihood of the joint estimate, given the data for each individual}
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
\code{\link{HIest}} finds maximum likelihood estimates for non-diagnostic markers.
}
\examples{
## A random set of allele counts for 10 loci and 10 individuals
G <- matrix(rbinom(100,2,.5),nrow=10)
HIC(G)

## The function is currently defined as
function(G){
	n <- nrow(G)
	k <- ncol(G)
	if(is.null(names(G))){Locus <- 1:k}else{Locus=names(G)}
	P <- data.frame(Locus=Locus,Allele="1",P1=1,P2=0)
	G <- replace(G,G<0,NA)
	out <- matrix(nrow=n,ncol=3)
	out[,1] <- rowMeans(G,na.rm=TRUE)/2
	out[,2] <- rowMeans(G==1,na.rm=TRUE)
	G <- replace(G,is.na(G),-9)
	for(i in 1:n){
		out[i,3] <- HILL(out[i,1:2],G[i,],P,type="allele.count")
	}
	colnames(out) <- c("S","H","logLik")
	out
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
