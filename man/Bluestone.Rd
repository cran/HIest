\name{Bluestone}
\alias{Bluestone}
\docType{data}
\title{Original marker data for hybrid tiger salamander larvae collected from Bluestone Quarry pond.
}
\description{Each row is an individual, each column is a marker. This is an example type="allele.count". Genotypes are 0 (homozygous for native California Tiger Salamander allele), 1 (heterozygous), or 2 (homozygous for introduced Barred Tiger Salamander allele). There are NA's.
}
\usage{data(Bluestone)}
\format{
  A data frame with 41 observations on 64 markers.
}
\source{
Fitzpatrick, B. M., J. R. Johnson, D. K. Kump, H. B. Shaffer, J. J. Smith, and S. R. Voss. 2009. Rapid fixation of non-native alleles revealed by genome-wide SNP analysis of hybrid tiger salamanders. BMC Evolutionary Biology 9:176. http://www.biomedcentral.com/1471-2148/9/176
}
\examples{
	\dontrun{
data(Bluestone)
BS.fit <- Cline.fit(Bluestone[,1:12], model = c("logit.logistic", "Barton"))
Cline.plot(BS.fit)

# # parental allele frequencies (assumed diagnostic)
BS.P <- data.frame(Locus=names(Bluestone),Allele="BTS",P1=1,P2=0)

# # estimate ancestry and heterozygosity
BS.est <-HIest(Bluestone,BS.P,type="allele.count")

# shortcut for diagnostic markers and allele count data:
BS.est <- HIC(Bluestone) 

# # calculate likelihoods for early generation hybrid classes
BS.class <- HIclass(Bluestone,BS.P,type="allele.count")

# # compare classification with maximum likelihood estimates
BS.test <- HItest(BS.class,BS.est)

table(BS.test$c1)
# # all 41 are TRUE, meaning the best classification is at least 2 log-likelihood units
# # better than the next best

table(BS.test$c2)
# # 2 are TRUE, meaning the MLE S and H are within 2 log-likelihood units of the best
# # classification, i.e., the simple classification is rejected in all but 2 cases

table(BS.test$Best.class,BS.test$c2)
# # individuals were classified as F2-like (class 3) or backcross to CTS (class 4),
# # but only two of the F2's were credible 

BS.test[BS.test$c2,]
# # in only one case was the F2 classification a better fit (based on AIC) than the
# # continuous model.
}
}
\keyword{datasets}
