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

Fitzpatrick, B. M. 2012. Estimating ancestry and heterozygosity of hybrids using molecular markers. BMC Evolutionary Biology 12:131. http://www.biomedcentral.com/1471-2148/12/131

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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
