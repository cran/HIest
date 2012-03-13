\name{HItest}
\alias{HItest}
\title{
Compare the likelihood of hybrid classification to MLE estimates of ancestry and heterozygosity.
}
\description{
\code{HItest} compares the best fit of six early generation diploid hybrid genotypes (parental, F1, F2, backcross) to the maximum likelihood genotype decribed by ancestry (S) and interclass heterozygosity (H). 
}
\usage{
HItest(class, MLE, thresholds = c(2, 1))
}
\arguments{
  \item{class}{
Output from \code{HIclass}: a data frame summarizing the fit of each individual to the six genotype classes.
}
  \item{MLE}{
Output from \code{HIest}: a data frame giving the MLE S and H and associated log-likelihood.
}
  \item{thresholds}{
Criteria for classification. The first criterion (\code{thresholds[1]}) is a cutoff for the difference in log-likelihood for the best vs. second best genotype class. The second criterion (\code{thresholds[2]}) is a cutoff for the difference in log-likelihood for the best genotype class vs. the MLE.
}
}
\details{
As a quick-and-dirty rule of thumb, one might accept a putative classification as credible if the log-likelihood of the best-fit class was over 2 units greater than the log-likelihood of the second best-fit class AND within 2 units of the maximum log-likelihood. The first criterion is based on the approximate equivalence of a 2 x log-likelihood interval to a 95 percent confidence interval for some distributions (Hudson 1971; Hillborn and Mangel 1997). The second is based on the conventional penalty of two log-likelihood units for an additional estimated parameter in model selection (Edwards 1972; Burnaham and Anderson 2004). The classification model can be viewed as having one free parameter (once the best-fit class is set to "chosen", the other five are constrained to "not chosen"), while the continuous model has two (S and H). This approach has the disadvantage of effectively treating the classification as a null model, which is not biologically justified. 

A better approach might be to accept the classification only if its AIC is lower than the AIC of the MLE, i.e., if \code{dAIC} is negative (Fitzpatrick 2012). Note that \code{dAIC} cannot be less than -2 (the case where the MLE is identical to the expectation for a class).
}
\value{
A data frame with one row per individual. Columns are:
\item{S }{The maximim likelihood estimate of the ancestry index from \code{HIest}.}
\item{H }{The maximum likelihood estimate of the interclass heterozygosity from \code{HIest}.}
\item{Best.class}{The class with the highest likelihood of the six from \code{HIclass}.}
\item{LL.class }{The log-likelihood of the data for the best-fit class from \code{HIclass}.}
\item{LLD.class}{The difference in log-likelihood between the best and second-best fit class from \code{HIclass}.}
\item{LL.max}{The maximum log-likelihood from \code{HIest}.}
\item{dAIC }{The difference in AIC between the continuous model MLE (2 estimated parameters) and the best-fit class (1 estimated parameter).}
\item{c1 }{Logical: \code{TRUE} if the best-fit class is supported by more than \code{thresholds[1]} log-likelihood units over the second best.}
\item{c2 }{Logical: \code{TRUE} if the best-fit class is WITHIN \code{thresholds[2]} log-likelihood units of the MLE.}
}
\references{
Burnham, K. P., and D. R. Anderson. 2004. Multimodel inference: understanding AIC and BIC in model selection. Sociological Methods and Research 33:261-304.

Edwards, A. W. F. 1972. Likelihood. Cambridge University Press, Cambridge.

Fitzpatrick, B. M. 2008. Hybrid dysfunction: Population genetic and quantitative genetic perspectives. American Naturalist 171:491-198.

Fitzpatrick, B. M. 2012. Estimating ancestry and heterozygosity of hybrids using molecular markers.

Hilborn, R., and M. Mangel. 1997. The ecological detective: Confronting models with data. Princeton University Press, New Jersey.

Hudson, D. J. 1971. Interval estimation from the likelihood function. Journal of the Royal Statistical Society, Series B 33: 256-262.

Lynch, M. 1991. The genetic interpretation of inbreeding depression and outbreeding depression. Evolution 45:622-629.
}
\author{
Ben Fitzpatrick
}

\seealso{
\code{\link{HIest}} for maximum likelihood estimation of S and H, \code{\link{HIsurf}} for a likelihood surface, \code{\link{HIclass}} for likelihoods of early generation hybrid classes, \code{\link{HILL}} for the basic likelihood function.
}
\examples{
data(Bluestone)
Bluestone <- replace(Bluestone,is.na(Bluestone),-9)
# parental allele frequencies (assumed diagnostic)
BS.P <- data.frame(Locus=names(Bluestone),Allele="BTS",P1=1,P2=0)

# estimate ancestry and heterozygosity
BS.est <-HIC(Bluestone)

# calculate likelihoods for early generation hybrid classes
BS.class <- HIclass(Bluestone,BS.P,type="allele.count")

# compare classification with maximum likelihood estimates
BS.test <- HItest(BS.class,BS.est)

table(BS.test$c1)
# all 41 are TRUE, meaning the best classification is at least 2 log-likelihood units
# better than the next best

table(BS.test$c2)
# 2 are TRUE, meaning the MLE S and H are within 2 log-likelihood units of the best
# classification, i.e., the simple classification is rejected in all but 2 cases.

table(BS.test$Best.class,BS.test$c2)
# individuals were classified as F2-like (class 3) or backcross to CTS (class 4), but
# only two of the F2's were credible 

BS.test[BS.test$c2,]
# in only one case was the F2 classification a better fit (based on AIC) than the
# continuous model.

## The function is currently defined as
function(class,MLE,thresholds=c(2,2)){
	LLC <- apply(class[,1:6],1,max)
	LLM <- MLE[,3]
	dAIC <- (2*1-2*LLC) - (2*2-2*LLM)
	c1 <- class$LLD>2
	c2 <- (LLM-LLC)<2
	data.frame(S=MLE[,1],H=MLE[,2],Best.class = class$Best,LL.class=LLC,
		LLD.class=class$LLD,LL.max=MLE[,3],dAIC,c1,c2)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
