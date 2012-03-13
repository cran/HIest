\name{HIC3}
\alias{HIC3}
\title{
Closed form maximum-likelihood estimates of ancestry and heterozygosity for diagnostic markers in a three-way hybrid zone
}
\description{
For hybrids with up to three parental lineages with diagnostic molecular markers (one allele fixed in each parental lineage). The function computes genomic proportions and ancestry indices using closed form estimators.
}
\usage{
HIC3(G, P)
}
\arguments{
  \item{G}{
A matrix or data frame of diploid genotypes. Each individual occupies two consecutive rows (one for each allele). Each marker is a column. For heterozygotes, it does not matter which allele is entered first.
}
  \item{P}{
A matrix or data frame identifying which alleles are fixed in each parental lineage. The first column stores locus identifiers and the second allele identifiers (one per parental lineage). The value of the third column is 1 for each allele fixed in the first parental lineage, and 0 otherwise. The value of the fourth column is 1 for each allele fixed in the second parental lineage, and 0 otherwise. The value of the fifth column is 1 for each allele fixed in the third parental lineage, and 0 otherwise. 
}
}
\value{
A matrix with 10 named columns for each individual
  \item{p11 }{Proportion of markers homozygous for lineage 1 alleles}
  \item{p22 }{Proportion of markers homozygous for lineage 2 alleles}
  \item{p33 }{Proportion of markers homozygous for lineage 3 alleles}
  \item{p12 }{Proportion of markers heterozygous for lineage 1 and 2 alleles}
  \item{p13 }{Proportion of markers heterozygous for lineage 1 and 3 alleles}
  \item{p23 }{Proportion of markers heterozygous for lineage 2 and 3 alleles}
  \item{S1 }{Lineage 1 ancestry index: proportion of alleles derived from parental lineage 1}
  \item{S2 }{Lineage 2 ancestry index: proportion of alleles derived from parental lineage 2}
  \item{S3 }{Lineage 3 ancestry index: proportion of alleles derived from parental lineage 3}
  \item{logLik }{log-likelihood of the genomic proportions given the individual marker data}

}
\references{
Fitzpatrick, B. M. 2012. Estimating ancestry and heterozygosity of hybrids using molecular markers.
}
\author{
Ben Fitzpatrick
}

\seealso{
\code{\link{threeway}} finds maximum likelihood estimates for non-diagnostic marker data. \code{\link{thirdclass}} and \code{\link{HItest3}} evaluate simple classification of three-way hybrids into parental, F1, F2, and backcross categories. For conventional two-way hybrid zone analyses, see \code{\link{HIC}, \link{HIest},link{HIclass},link{HItest}}.
}
\examples{
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

HIC3(G,P)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
