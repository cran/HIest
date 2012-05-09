\name{thirdclass}
\alias{thirdclass}
\title{
Calculate likelihoods for early generation hybrid genotype classes in a three-way hybrid zone.
}
\description{
\code{thirdclass} uses genetic marker data and parental allele frequencies to calculate the likelihoods for each of the 15 diploid genotype classes possible in the first two generations of admixture (each parental, F1, F2, and each backcross) for each pair of parental lineages.
}
\usage{
thirdclass(G, P, type = "codominant")
}
\arguments{
  \item{G}{
A matrix or data frame of genetic marker data. Each column is a locus. For \code{type="dominant"}, there should be one row per individual. For \code{type="codominant"}, each individual is to be represented in consecutive rows (one for each allele).
}
  \item{P}{
A matrix or data frame with the following columns (order is important!): Locus name, Allele name, P1 allele frequency, P2 allele frequency, P3 allele frequency. For \code{type="dominant"}, there should be one row per locus, giving the frequencies of the dominant or "1" allele. For \code{type="codominant"} there should be a separate row for each allele AND the Allele names should match the data in \code{G}.
}
  \item{type}{
A string representing the data type. The options are \code{"codominant"} and \code{"dominant"}.
}
}
\details{
The function evaluates genotype classes including only two parental lineages in the context of the possibility of a third parental lineage. It would be straightforward but tedious to extend it to include complex classes such as the offspring of a 1x2 F1 with a 2x3 F1. However, it is not entirely clear how to account for the additional complexity introduced by allowing more classes. 
}
\value{
A data matrix with the following columns (one row per individual)
\item{c100000}{log-likelihood for genotype class expected for pure parental P1}
\item{c010000}{log-likelihood for genotype class expected for P2}
\item{c001000}{log-likelihood for genotype class expected for P3}
\item{c000100}{log-likelihood for genotype class expected for 1x2 F1 hybrids}
\item{c110200}{log-likelihood for genotype class expected for 1x2 F2 hybrids}
\item{c100100}{log-likelihood for genotype class expected for 1x1x2 backcrosses}
\item{c010100}{log-likelihood for genotype class expected for 1x2x2 backcrosses}
\item{c000010}{log-likelihood for genotype class expected for 1x3 F1 hybrids}
\item{c101020}{log-likelihood for genotype class expected for 1x3 F2 hybrids}
\item{c100010}{log-likelihood for genotype class expected for 1x1x3 backcrosses}
\item{c001010}{log-likelihood for genotype class expected for 1x3x3 backcrosses}
\item{c000001}{log-likelihood for genotype class expected for 2x3 F1 hybrids}
\item{c011002}{log-likelihood for genotype class expected for 2x3 F2 hybrids}
\item{c010001}{log-likelihood for genotype class expected for 2x2x3 backcrosses}
\item{c001001}{log-likelihood for genotype class expected for 2x3x3 backcrosses}
\item{Best}{The class with the highest likelihood of the 15}
\item{LLD}{The difference in log-likelihood between the best and second-best fit class. This can be used as a criterion for deciding whether the best fit class is enough better to reject the alternatives.}
}
\references{
Fitzpatrick, B. M. 2012. Estimating ancestry and heterozygosity of hybrids using molecular markers. BMC Evolutionary Biology 12:131. http://www.biomedcentral.com/1471-2148/12/131
}
\author{
Ben Fitzpatrick
}

\seealso{
\code{\link{threeway}} finds maximum likelihood estimates for non-diagnostic marker data. \code{\link{HItest3}} compares classification to the continuous model MLE. For conventional two-way hybrid zone analyses, see \code{\link{HIC}, \link{HIest}, \link{HIclass}, \link{HItest}}.
}
\examples{
	\dontrun{
	## all possible 2-way crosses after 2 generations
G <- rbind(
rep(1,12),rep(1,12),               # parental 1
rep(2,12),rep(2,12),               # parental 2
rep(3,12),rep(3,12),               # parental 3
rep(1,12),rep(2,12),               # 1 x 2 F1
rep(1:2,each=6),rep(1:2,6),        # 1 x 2 F2
rep(1,12),rep(1:2,6),              # 1 x 1 x 2 BC
rep(2,12),rep(1:2,6),              # 1 x 2 x 2 BC
rep(1,12),rep(3,12),               # 1 x 3 F1
rep(c(1,3),each=6),rep(c(1,3),6),  # 1 x 3 F2
rep(1,12),rep(c(1,3),6),           # 1 x 1 x 3 BC
rep(3,12),rep(c(1,3),6),           # 1 x 3 x 3 BC
rep(2,12),rep(3,12),               # 2 x 3 F1
rep(2:3,each=6),rep(2:3,6),        # 2 x 3 F2
rep(3,12),rep(2:3,6),              # 2 x 3 x 3 BC
rep(2,12),rep(2:3,6)               # 2 x 2 x 3 BC
)

P <- data.frame(Locus=rep(1:12,each=3),allele=rep(1:3,12),P1=rep(c(1,0,0),12),
	P2=rep(c(0,1,0),12),P3=rep(c(0,0,1),12))

Est <- HIC3(G,P)
Class <- thirdclass(G,P)
HItest3(Class,Est)

## now for some three-way mixes
G3 <- matrix(1+rbinom(200,2,.5),ncol=10)
Est3 <- HIC3(G3,P)
Class3 <- thirdclass(G3,P)
HItest3(Class3,Est3)   	# usually all classifications will be impossible because all 
						# individuals will have nonzero contributions from each of the 
						# three parentals

## bias toward parental 1
G3 <- matrix(1+rbinom(200,2,.25),ncol=10)
Est3 <- HIC3(G3,P)
Class3 <- thirdclass(G3,P)
HItest3(Class3,Est3) 	# now you might have a few that look like F2's 
						# between 1 and 2 (c110200)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
