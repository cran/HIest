\name{Cline.fit}
\alias{Cline.fit}
\title{
Fit alternative cline functions to data from one or more genetic loci.
}
\description{
This function takes individual or population data for multiple genetic loci and fits genomic clines. This implementation assumes markers are diagnostic.
}
\usage{
Cline.fit(Data, By = NULL, S = NULL, model,
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
  \item{model}{
Character vector specifying which models to fit. Valid choices are "multinom","binom","logit.logistic","Barton","Beta", and "Richards".
}
  \item{Start}{
Optional starting values for optimization. If \code{NULL}, start values at the 1:1 line are provided by default.
}
  \item{Methods}{
Optional named list of strings indicating the optimization method for each cline model (logit-logistic, Barton, Beta, and Richards). Choices are \code{"L-BFGS-B"}, \code{"SANN"}, and \code{"mcmc"}. Default is \code{Methods=list(logit.logistic="L-BFGS-B",Barton="L-BFGS-B",Beta="L-BFGS-B",Richards="L-BFGS-B")}.
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
A list including a named matrix for each fitted model. Each such matrix includes the parameter estimates, maximum log-likelihood, AICc, squared Mahalanobis distance \code{D2}, P-value for Mahalanobis-based outlier detection, and TRUE/FALSE declaration of whether a marker is an outlier based on a Bonferroni-adjusted critical P-value. If more than one model was fit, an additional data frame \code{best.fit} is included, giving the model with lowest AIC for each marker.
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
BS.fit <- Cline.fit(Bluestone[,1:12],model=c("logit.logistic","Barton"))
Cline.plot(BS.fit)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
