\name{HIest-package}
\alias{HIest-package}
\docType{package}
\title{
Functions for estimating hybrid indices (ancestry and heterozygosity) and fitting genomic clines.
}
\description{
Uses likelihood to estimate ancestry, heterozygosity, and genomic cline parameters. Also evaluates simple hybrid classifications (parentals, F1, F2, backcrosses).
}
\details{
\tabular{ll}{
Package: \tab HIest\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-02-13\cr
License: \tab GPL 3.0 \cr
LazyLoad: \tab yes\cr
}
}
\author{
Benjamin M. Fitzpatrick
Maintainer: <benfitz@utk.edu>
}
\references{
Fitzpatrick, B. M. 2008. Hybrid dysfunction: Population genetic and quantitative genetic perspectives. American Naturalist 171:491-198.

Fitzpatrick, B. M. 2012. Estimating ancestry and heterozygosity of hybrids using molecular markers.

Fitzpatrick, B. M. 2012. Alternative forms for genomic clines.

Lynch, M. 1991. The genetic interpretation of inbreeding depression and outbreeding depression. Evolution 45:622-629.
}
\keyword{ package }
%\seealso{}
\examples{
data(Bluestone)

######################
# Fit genomic clines #
######################

BS.fit <- Cline.fit(Bluestone)
Cline.plot(BS.fit)

########################################
# Estimate ancestry and heterozygosity #
########################################

Bluestone <- replace(Bluestone,is.na(Bluestone),-9)

# parental allele frequencies (assumed diagnostic)
BS.P <- data.frame(Locus=names(Bluestone),Allele="BTS",P1=1,P2=0)

# estimate ancestry and heterozygosity
BS.est <-HIC(Bluestone)

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
# individuals were classified as F2-like (class 3) or backcross to CTS (class 4), but
# only two of the F2's were credible 

BS.test[BS.test$c2,]
# in only one case was the F2 classification a better fit (based on AIC) than the
# continuous model.

# equivalent to the AIC criterion:
BS.test <- HItest(BS.class,BS.est,thresholds=c(2,1))

#########################
# three-way hybrid zone #
#########################

# for example: make each parental, F1, F2, and backcross
G <- rbind(
rep(1,12),rep(1,12),
rep(2,12),rep(2,12),
rep(3,12),rep(3,12),
rep(1,12),rep(2,12),
rep(1:2,each=6),rep(1:2,6),
rep(1,12),rep(1:2,6),
rep(2,12),rep(1:2,6),
rep(1,12),rep(3,12),
rep(c(1,3),each=6),rep(c(1,3),6),
rep(1,12),rep(c(1,3),6),
rep(3,12),rep(c(1,3),6),
rep(2,12),rep(3,12),
rep(2:3,each=6),rep(2:3,6),
rep(3,12),rep(2:3,6),
rep(2,12),rep(2:3,6)
)

# 12 diagnostic markers
P <- data.frame(Locus=rep(1:12,each=3), allele=rep(1:3,12), P1=rep(c(1,0,0),12), 
P2=rep(c(0,1,0),12), P3=rep(c(0,0,1),12))

# find MLE with simulated annealing ... takes a few minutes with default iterations
# Est <- threeway(G,P,method="SANN",surf=FALSE)

# shortcut for diagnostic markers
Est <- HIC3(G,P)
CL <- thirdclass(G,P)

}
