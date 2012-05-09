\name{HIclass}
\alias{HIclass}
\title{
Calculate likelihoods for early generation hybrid genotype classes
}
\description{
\code{HIclass} uses genetic marker data and parental allele frequencies to calculate the likelihoods for each of the six diploid genotype classes possible in the first two generations of admixture (each parental, F1, F2, and each backcross)
}
\usage{
HIclass(G, P, type)
}
\arguments{
  \item{G}{
A matrix or data frame of genetic marker data. Each column is a locus. For \code{type="dominant"} or \code{type="allele.count"}, there should be one row per individual. For \code{type="codominant"}, each individual is to be represented in consecutive rows (one for each allele).
}
  \item{P}{
A matrix or data frame with the following columns (order is important!): Locus name, Allele name, P1 allele frequency, P2 allele frequency. For \code{type="dominant"} or \code{type="allele.count"}, there should be one row per locus, giving the frequencies of the dominant or "1" allele. For \code{type="codominant"} there should be a separate row for each allele AND the Allele names should match the data in \code{G}.
}
  \item{type}{
A string representing the data type. The options are \code{"codominant"}, \code{"dominant"}, and \code{"allele.count"}.
}
}
\details{
Classification should be used only for systems for which two generations of admixture are credible. \code{\link{HIest}} can be used to find the joint maximum likelihood estimates of ancestry (S) and interclass heterozygosity (H), which offer a more nuanced summary of hybrid genotypes.
}
\value{
A data frame with the following columns (one row per individual)
\item{class100}{log-likelihood for genotype class expected for pure parental P2 (P = 0, H = 0)}
\item{class010}{log-likelihood for genotype class expected for F1 hybrids (P = 0.5, H = 1)}
\item{class121}{log-likelihood for genotype class expected for F2 hybrids (P = 0.5, H = 0.5)}
\item{class110}{log-likelihood for genotype class expected for backcrosses to parental 2 (P = 0.25, H = 0.5)}
\item{class011}{log-likelihood for genotype class expected for backcrosses to parental 1 (P = 0.75, H = 0.5)}
\item{class001}{log-likelihood for genotype class expected for pure parental P1 (P = 1, H = 0)}
\item{Best}{The class with the highest likelihood of the six}
\item{LLD}{The difference in log-likelihood between the best and second-best fit class. This can be used as a criterion for deciding whether the best fit class is enough better to reject the alternatives.}
}
\references{
Fitzpatrick, B. M. 2008. Hybrid dysfunction: Population genetic and quantitative genetic perspectives. American Naturalist 171:491-198.

Fitzpatrick, B. M. 2012. Estimating ancestry and heterozygosity of hybrids using molecular markers. BMC Evolutionary Biology 12:131. http://www.biomedcentral.com/1471-2148/12/131

Lynch, M. 1991. The genetic interpretation of inbreeding depression and outbreeding depression. Evolution 45:622-629.
}
\author{
Ben Fitzpatrick}

\seealso{
\code{\link{HIest}} for maximum likelihood estimation of S and H, \code{\link{HIsurf}} for a likelihood surface, \code{\link{HItest}} to compare the classification to the maximum likelihood, \code{\link{HILL}} for the basic likelihood function.
}
\examples{
data(Bluestone)
Bluestone <- replace(Bluestone,is.na(Bluestone),-9)
# parental allele frequencies (assumed diagnostic)
BS.P <- data.frame(Locus=names(Bluestone),Allele="BTS",P1=1,P2=0)

# estimate ancestry and heterozygosity
# BS.est <-HIest(Bluestone,BS.P,type="allele.count")

# shortcut for diagnostic markers
BS.est <- HIC(Bluestone)

# calculate likelihoods for early generation hybrid classes
BS.class <- HIclass(Bluestone,BS.P,type="allele.count")

# compare classification with maximum likelihood estimates
BS.test <- HItest(BS.class,BS.est,thresholds=c(2,2))

table(BS.test$c1)
# all 41 are TRUE, meaning the best classification is at least 2 log-likelihood units
# better than the next best

table(BS.test$c2)
# 2 are TRUE, meaning the MLE S and H are within 2 log-likelihood units of the best
# classification, i.e., the simple classification is rejected in all but 2 cases

table(BS.test$Best.class,BS.test$c2)
# individuals were classified as F2-like (class 3) or backcross to CTS (class 4),
# but only two of the F2's were credible 

BS.test[BS.test$c2,]
# in only one case was the F2 classification a better fit (based on AIC) than the
# continuous model.

# equivalent to the AIC criterion:
BS.test <- HItest(BS.class,BS.est,thresholds=c(2,1))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
