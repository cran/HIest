\name{HILL}
\alias{HILL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hybrid Index Log Likelihood
}
\description{
\code{HILL} calulates the log-likelihood of a hybrid genotype given parental allele frequencies, an ancestry index (S), and interclass heterozygosity (H). This function is used by \code{\link{HIest}} and \code{\link{HIclass}} for data analysis.
}
\usage{
HILL(par = c(S, H), G, P, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
A numeric vector including a value of S and a value of H
}
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
}
\details{
Given two ancestral species or parental populations (P1 and P2), the ancestry index (S) is the proportion of an individual's alleles descending from alleles in the P1 population and the interclass heterozygosity (H) is the proportion of an individual's loci that have one allele from each ancestral population (Lynch 1991). The likelihood functions are described in Fitzpatrick (2012).
}
\value{
The natural log of the likelihood of the data given S, H, and P.
}
\references{
Fitzpatrick, B. M. 2008. Hybrid dysfunction: Population genetic and quantitative genetic perspectives. American Naturalist 171:491-198.

Fitzpatrick, B. M. 2012. Estimating ancestry and heterozygosity of hybrids using molecular markers. BMC Evolutionary Biology 12:131. http://www.biomedcentral.com/1471-2148/12/131

Lynch, M. 1991. The genetic interpretation of inbreeding depression and outbreeding depression. Evolution 45:622-629.
}
\author{
Ben Fitzpatrick
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{HIest}} for maximum likelihood estimation of S and H, \code{\link{HIsurf}} for a likelihood surface, \code{\link{HIclass}} for likelihoods of early generation hybrid classes, \code{\link{HItest}} to compare the classification to the maximum likelihood.
}
\examples{
##-- A random codominant data set of 1 individual with 5 markers with three possible alleles each
P <- data.frame(Locus=rep(1:5,each=3),Allele=rep(1:3,5),
	P1=as.vector(rmultinom(5,10,c(.7,.2,.1)))/10, 
	P2=as.vector(rmultinom(5,10,c(.1,.2,.7)))/10)
G <- matrix(nrow=2,ncol=5)
for(i in 1:5){
	G[,i] <- sample(c(1,2,3),size=2,replace=TRUE,prob=rowMeans(P[P$Locus==i,3:4]))
	}
	
HILL(par=c(0.5,0.5),G,P,type="codominant")

##-- Make it into allele count data (count "3" alleles)
P.c <- P[seq(from=3,to=15,by=3),]
G.c <- colSums(G==3)

HILL(par=c(0.5,0.5),G.c,P.c,type="allele.count")

##-- Make it into dominant data where allele 3 is dominant
G.d <- replace(G.c,G.c==2,1)

HILL(par=c(0.5,0.5),G.d,P.c,type="dominant")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
