\name{Cline.plot}
\alias{Cline.plot}
\title{
Graph genomic clines fit to data
}
\description{
A quick way to summarize and view the results of \code{\link{Cline.fit}}
}
\usage{
Cline.plot(cfit, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cfit}{
An object: the output of \code{\link{Cline.fit}}.
}
  \item{\dots}{
Optional arguments to be passed to \code{\link{plot}}.
}
}
\details{
Creates up to six graphs, depending on which models are fit to the data. If all models are fit, the last panel shows the frequency with which the Barton, Beta, and Richards clines were best fit according to AIC. Lighter shading in this final barplot represents the fraction of markers where the null (1:1 line) could be rejected according to the minimal criterion that any fitted model had lower AIC. 
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
Cline.plot(BS.fit)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
