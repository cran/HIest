\name{Cline.plot}
\alias{Cline.plot}
\title{
Graph genomic clines fit to data
}
\description{
A quick way to summarize and view the results of \code{\link{Cline.fit}}
}
\usage{
Cline.plot(cfit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cfit}{
An object: the output of \code{\link{Cline.fit}}.
}
}
\details{
Creates up to sixteen graphs, depending on which models are fit to the data. For each fitted model, fitted clines are graphed first. For two-parameter models, a bivariate graph of parameter estimates is plotted next (this is left blank for four-parameter models). Finally, a Q-Q plot of squared Mahalanobis distances vs. CHi-squared quantiles is graphed. Putative outlier markers are indicated in red.
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
