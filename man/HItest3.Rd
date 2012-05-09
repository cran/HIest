\name{HItest3}
\alias{HItest3}
\title{
Compare the likelihood of hybrid classification to MLE estimates of ancestry and heterozygosity for three-way hybrid zones.
}
\description{
\code{HItest3} compares the best fit of fifteen early generation diploid hybrid genotypes (parental, F1, F2, backcross between all three pairs of parental lineages) to the maximum likelihood genotype proportions estimated by \code{threeway} or \code{HIC3}. 
}
\usage{
HItest3(class, MLE, thresholds = c(2, 8))
}
\arguments{
  \item{class}{
Object containing output from \code{thirdclass}
}
  \item{MLE}{
Object containing output from \code{threeway} or \code{HIC3}.
}
  \item{thresholds}{
Vector of length 2, containing criteria for classification. The first criterion
 (\code{thresholds[1]}) is a cutoff for the difference in log-likelihood for the best vs. second best genotype class. The second criterion (\code{thresholds[2]}) is a cutoff for the difference in log-likelihood for the best genotype class vs. the MLE. For the three-way model, \code{thresholds[2] = 8} corresponds to the criterion that the best classification must have a lower AIC than the MLE to favor classification over the continuous model.
}
}
\details{
The AIC for the continuous model accounts for k = 5 estimated paramaters. For the classification model, on a per individual basis, there is only one estimated parameter. Classification could be expanded to include more complex genotypic combinations, but it is not clear how (or whether) to account for additional complexity in model comparison. It is probably advisable to consider testing only classifications reflecting a clear biological question, such as whether F1 hybrids are completely sterile.
}
\value{
A data matrix with 7 named columns and one row per individual.
  \item{Best.class }{Most likely classification of the individual given the data. Classes are indicated by their expected genomic proportions in order: \code{p11, p22, p33, p12, p13, p23}. For example, a pure parental 1 would be \code{c100000}, an F1 between parentals 1 and 3 would be \code{c000010}, an F2 between parentals 2 and 3 would be \code{c011002}, etc. NOTE, if two or more classes are equally likely, only the first will be reported here, but \code{LLD.class} will be 0 or \code{NaN}.
  	}
  \item{LL.class }{The log-likelihood of the best classification given the data. This will be \code{-Inf} if no classification is theoretically possible.}
  \item{LLD.class }{The difference in log-likelihood between the best and second best classification.}
  \item{LL.max }{The log-liklihood of the MLE genomic proportions, not constrained to fit a class (continuous model).}
  \item{dAIC }{Difference in AIC between the MLE and best classification. This will be negative when the AIC for classification (1 parameter) is lower than the AIC for the continuous model MLE (5 parameters).}
\item{c1 }{Logical: \code{TRUE} if the best-fit class is supported by more than \code{thresholds[1]} log-likelihood units over the second best.}
\item{c2 }{Logical: \code{TRUE} if the best-fit class is WITHIN \code{thresholds[2]} log-likelihood units of the continuous model MLE.}
}
\references{
Fitzpatrick, B. M. 2012. Estimating ancestry and heterozygosity of hybrids using molecular markers. BMC Evolutionary Biology 12:131. http://www.biomedcentral.com/1471-2148/12/131
}
\author{
Ben Fitzpatrick
}
\seealso{
\code{\link{threeway}} finds maximum likelihood estimates for non-diagnostic marker data. \code{\link{thirdclass}} evaluates simple classification of three-way hybrids into parental, F1, F2, and backcross categories. For conventional two-way hybrid zone analyses, see \code{\link{HIC}, \link{HIest}, \link{HIclass}, \link{HItest}}.
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
HItest3(Class3,Est3)    # usually all classifications will be impossible because all
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
