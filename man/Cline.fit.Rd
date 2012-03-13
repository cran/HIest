\name{Cline.fit}
\alias{Cline.fit}
\title{
Fit alternative cline functions to data from one or more genetic loci.
}
\description{
This function takes individual or population data for multiple genetic loci and fits genomic clines and functions describing interclass heterozygosity. This implementation assumes markers are diagnostic.
}
\usage{
Cline.fit(Data, By = NULL, S = NULL,
 models = c(bart = TRUE, beta = TRUE, rich = TRUE, het1 = TRUE, het2 = TRUE),
  Start = NULL, Methods = NULL, iterations = 99, SD = NULL, headstart = TRUE,
   Grid = FALSE, ploidy = 2, trim = 0, include = 1:ncol(Data))
}
\arguments{
  \item{Data}{
Data matrix of allele counts (allele diagnostic for one species desgnated 1, all others 0) by individual or sample (rows) and marker (columns).
}
  \item{By}{
Optional factor defining aggregation of individuals into samples (e.g., sites or populations). Length must be equal to the number of rows in \code{Data}.
}
  \item{S}{
Optional alternative values for the genome wide hybrid index. By default, the mean ancestry across all loci in \code{Data} will be used. If \code{S} is specified, the values given by the user trump the default. 
}
  \item{models}{
Vector of five logical TRUE/FALSE indicators of which models to consider in the follwoing order: Barton, Beta, Richards, het1, and het2. het1 is a symmetrical power function for interclass heterozygosity \eqn{E(h)=(d2S(1-S))^z}, and het2 is an asymmetrical alternative \eqn{E(h)=d2S^{z_1}(1-S)^{z_2}}.Note that the asymmetrical heterozygosity curves are not necessarily tied down at zero and one.
}
  \item{Start}{
Optional starting values for optimization. If \code{NULL}, start values at the 1:1 line are provided by default.
}
  \item{Methods}{
Optional vector with three strings indicating the optimization method for each cline model (Barton, Beta, and Richards). Choices are \code{"L-BFGS-B"} (the default for all three), \code{"SANN"}, and \code{"mcmc"}.
}
  \item{iterations}{
The number of MCMC generations to use if \code{"mcmc"} is used for any \code{Methods}. 
}
  \item{SD}{
Dispersion parameters for the \code{"mcmc"} and \code{"SANN"} methods. In these methods, new parameter values are proposed by drawing values from normal distributions centered on the current value and with standard deviations from \code{SD}.
}
  \item{headstart}{
Logical: if TRUE and \code{method="mcmc"} or \code{"SANN"}, starting values will be found by first using \code{optim} with \code{"L-BFGS-B"}.
}
  \item{Grid}{Logical: if TRUE and \code{method="mcmc"} and \code{model="Beta"}, starting values for the Markov Chain will be found by finding the highest likelihood on a 100 x 100 grid made by \code{mu <- seq(from=0.02,to=0.90,length.out=10);	nu <- 2^(0:9)/10}.}
  \item{ploidy}{Ploidy of the data, i.e., an integer (1 or 2) to multiply the sample size (number of individuals).}
  \item{trim}{Optional fraction of extreme values to omit from calculation of \eqn{S} (see \code{\link{mean}}).}
  \item{include}{Optional vector of column indices to include in calculation of \eqn{S} - can be used to omit biological outliers.}
}

\value{
	A list including

  \item{Models }{Strings indicating which cline models were fit.}
  \item{Comparison }{Data frame summarizing the model fits. For each marker, the locus name is given, followed by \eqn{\Delta}AIC for the null model and each fitted model, and a column stating which model was best according to AIC.}
  \item{Barton }{Data frame giving the fitted coefficients of the Barton cline, maximum likelihood, likelihood of the null (1:1 line), likelihood ratio p-value, and AICc.}
  \item{Beta }{As obove for the Beta cline}
  \item{Richards }{As above for the Richards cline}
  \item{Heterozygosity }{As above for the symmetrical heterozygosity curve, with the null being the Hardy-Weinberg reference: \eqn{2S(1-S)}}
  \item{Hbeta }{As above for the asymmetrical heterozygosity curve}

}
\references{
Fitzpatrick, B. M. 2012. Alternative forms for genomic clines. In prep
}
\author{
Benjamin M. Fitzpatrick}

\seealso{
See \code{\link{gcline.fn}} for the basic fitting function. \code{\link{Cline.plot}} provides an easy way to visualize the output of \code{Cline.fit}.}
\examples{
	
data(Bluestone)
BS.fit <- Cline.fit(Bluestone)
#par(mfrow=c(2,3))
Cline.plot(BS.fit)

# pretend this was actually 10 sites with different mean ancestry
S <- apply(Bluestone,1,FUN=mean,na.rm=TRUE)/2
BS10 <- Bluestone[order(S),]
Sites <- c(1,rep(1:10,each=4))

BS.fitS <- Cline.fit(BS10,By=Sites)
Cline.plot(BS.fitS)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
