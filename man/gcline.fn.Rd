\name{gcline.fn}
\alias{gcline.fn}
\title{
Fit a genomic cline using maximum likelihood}
\description{
Fit a genomic cline and compare it to a null expectation. Choices are Barton cline, Beta cline, and Richards cline. This function is used by \code{\link{Cline.fit}} to compare clines fit to a data set.}
\usage{
gcline.fn(x, n, y, start, model = "Beta", method = "L-BFGS-B", iterations = 99,
 SD = rep(0.01, length(start)), headstart = FALSE,Grid = TRUE)
}
\arguments{
  \item{x}{
A numeric vector of genome-wide mean ancestry (or any independent variable on the unit interval).}
  \item{n}{
A numeric vector of sample sizes for each value in \code{x}}
  \item{y}{
A numeric vector containing the dependent variable: usually an allele frequency for each \code{x}
}
  \item{start}{
A vector of starting values: \eqn{\mu} and \eqn{\nu} for the Beta cline, \eqn{a} and \eqn{b} for the Barton cline, and \eqn{U}, \eqn{L}, \eqn{m}, and \eqn{b} for the Richards cline.
}
  \item{model}{
Character string indicating which cline function to fit: \code{"Beta"}, \code{"Barton"}, or \code{"Richards"}
}
  \item{method}{
Character string indicating which algorithm to use to find the MLE: \code{"L-BFGS-B"} and \code{"SANN"} are used by the native optimization function \code{\link{optim}}, \code{"mcmc"} is a Markov Chain Monte Carlo using Metropolis-Hastings sampling. If \code{"mcmc"} is used, the following four items are also used.
}
  \item{iterations}{
The desired number of MCMC generations. The larger this number is, the greater the chance that the chain will find the maximum likelihood.
}
  \item{SD}{
Dispersion parameters for the \code{"mcmc"} and \code{"SANN"} methods. In these methods, new parameter values are proposed by drawing values from normal distributions centered on the current value and with standard deviations from \code{SD}.
}
  \item{headstart}{
Logical: if TRUE and \code{method="mcmc"} or \code{"SANN"}, starting values will be found by first using \code{optim} with \code{"L-BFGS-B"}.
}
  \item{Grid}{
  	Logical: if TRUE and \code{method="mcmc"} and \code{model="Beta"}, starting values for the Markov Chain will be found by finding the highest likelihood on a 100 x 100 grid made by \code{mu <- seq(from=0.02,to=0.90,length.out=10);	nu <- 2^(0:9)/10}.
  }
}
\value{
A list:
	\item{model }{The cline function used}
    \item{method }{The optimization method used}
    \item{estimates}{Maximum likelihood parameter estimates}
    \item{lnL}{The maximum likelihood and the likelihood of the data given the null model \code{E(y)=x}}
    \item{k}{The number of fitted coefficients}
    \item{AICc}{Akaike's information criterion with sample size correction
    	}
    \item{convergence}{From \code{\link{optim}}: 0 means the algorithm thinks it did a good job, otherwise the MLE may be incorrect. If \code{mcmc} was used, this will be the full Markov Chain, which could be analyzed for convergence (e.g., see \code{coda}).}
}
\references{
Fitzpatrick, B. M. 2012. Alternative forms for genomic clines. In prep
}
\author{
Benjamin M. Fitzpatrick}

\seealso{
See \code{\link{Cline.fit}} for application to multilocus data sets and analysis of interclass heterozygosity.
}
\examples{
x <- 0:50/50          # hypothetical genomic ancestry proportions
y <- rbinom(50,2,x)   # random diploid genotypes for a diagnostic marker
n=rep(2,50)			  # sample size is two alleles per diploid individual

gcline.fn(x=x,n=n,y=y,start=c(.5,2),model="Beta")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
