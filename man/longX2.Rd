\name{longX2}
\alias{longX2}
\title{
Longs heterogeneity test for admixture proportions in a single population.
}
\description{
For hybrid individuals genotyped with diagnostic markers (one allele fixed in each parental lineage), the null hypothesis that variation in allele frequencies is explained by sampling and drift (no selection) is evaluated with a chi-squared test statistic.
}
\usage{
longX2(Freqs)
}
\arguments{
  \item{Freqs}{
A vector of allele frequencies (one estimate per locus in a single population).
}
}
\value{
A list with the following items is returned.
  \item{test }{a list including the global test statistic, degrees of freedom, and p-value.}
  \item{chisq.res }{Chi-squared residuals for each marker}
}
\references{
Long, J. C. (1991). The genetic structure of admixed populations. Genetics, 127:417-428.
Fitzpatrick, B. M., Johnson, J. R., Kump, D. K., Shaffer, H. B., Smith, J. J., and Voss, S. R. (2009). Rapid fixation of non-native alleles revealed by genome-wide snp analysis of hybrid tiger salamanders. BMC Evolutionary Biology, 9:176.
}
\author{
Ben Fitzpatrick
}

\examples{
longX2(c(0.95,rbeta(20,2,5)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
