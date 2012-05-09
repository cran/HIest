\name{threeway}
\alias{threeway}
\title{
Find joint maximum-likelihood estimates of ancestry and heterozygosity for a sample of hybrids in a three-way hybrid zone
}
\description{
For hybrids with up to three parental lineages, the function estimates genomic proportions and ancestry indices by numerically searching likelihood space.
}
\usage{
threeway(G, P, type = "codominant", surf = TRUE, method = "SANN", iterations = 500, 
	start = rep(1/6, 6), Cscale = NULL, props = NULL, 
	control = list(fnscale = -1, maxit = iterations))
}
\arguments{
  \item{G}{
A matrix or data frame of genetic marker data. Each column is a locus. For \code{type="dominant"}, there should be one row per individual. For \code{type="codominant"}, each individual is to be represented in consecutive rows (one for each allele).
}
  \item{P}{
A matrix or data frame with the following columns (order is important!): Locus name, Allele name, P1 allele frequency, P2 allele frequency, P3 allele frequency. For \code{type="dominant"}, there should be one row per locus, giving the frequencies of the dominant or "1" allele. For \code{type="codominant"} there should be a separate row for each allele AND the Allele names should match the data in \code{G}.
}
  \item{type}{
A string representing the data type. The options are \code{"codominant"} or \code{"dominant"}.
}
  \item{surf}{
Logical: should the function find starting values by evaluating likelihoods on a grid?
}
  \item{method}{
  	Optimization method to search for maximum likelihood estimates of ancestry and heterozygosity. Alternatives are \code{"SANN"},  \code{"L-BFGS-B"}, \code{"surf"}, and \code{"mcmc"}. See details.
}
  \item{iterations}{
The desired number of MCMC steps to perform when \code{method="mcmc" or "SANN"}.
}
  \item{start}{
A vector including the starting values of the six genomic proportions. See details.
}
  \item{Cscale}{
An integer, controlling the the proposal distribution for \code{method = "SANN"} or \code{"mcmc"}. Smaller values will cause the algorithm to search more broadly, but could make the search inefficient. See details.
}
  \item{props}{
Optional: a matrix of genomic proportions to evaluate if \code{surf = TRUE} or \code{method = "surf"}. Columns correspond to the six genomic proportions in order: p11, p22, p33, p12, p13, p23.
}
  \item{control}{
A list of options to be passed to \code{control} in the \code{\link{optim}} function. Whatever else is chosen, be sure \code{fnscale} is negative to make \code{optim} search for a maximum rather than a minimum. Specifying \code{maxit} will override \code{iterations} for \code{method = "SANN"}.
}
}
\details{
Given three ancestral species or parental populations (P1, P2, and P3), the genome of a hybrid can be described by six genomic proportions: 
\eqn{p_{11}} = proportion of one's genome that is homozygous for alleles inherited from P1,
\eqn{p_{22}} = proportion of one's genome that is homozygous for alleles inherited from P2,
\eqn{p_{33}} = proportion of one's genome that is homozygous for alleles inherited from P3,
\eqn{p_{12}} = proportion of one's genome that is heterozygous for alleles inhertied from P1 and P2,
\eqn{p_{13}} = proportion of one's genome that is heterozygous for alleles inhertied from P1 and P3,
\eqn{p_{23}} = proportion of one's genome that is heterozygous for alleles inhertied from P2 and P3.

Currently, the function provides four methods for searching the likelihood space (a 5-simplex with vertices wherever each genomic proportion is equal to 1.0). \code{method = "SANN"} is probably the best; it uses the general purpose optimization function \code{\link{optim}} with its simulated annealing algorithm. For this estimation problem, a custom proposal function is passed to the option \code{gr}. This proposal function draws new genomic proportions from a 6-dimensional Dirichlet distribution centered on the old genomic proportions. The concentration of the proposal distribution is controlled by \code{Cscale}; the larger this value, the more the proposal distribution is concentrated near the current state.

\code{method = "mcmc"} uses a Markov-Chain Monte Carlo with Metropolis-Hastings sampling to explore the likelihood space. It also uses the Dirichlet proposal distribution, and could be useful (with some modification of the code) for generating posterior distributions. \code{method = "SANN"} is probaby superior for simply finding the MLE.

\code{method = "L-BFGS-B"} also uses \code{\link{optim}}, but with a quasi-Newton likelihood search algorithm to look for the maximum likelihood. This method is relatively fast, but it can miss the MLE if it is near an edge of the sample space.

\code{"surf"}, finds all likelihoods on a 6-dimensional grid defined by \code{props} and chooses the maximum. By default, \code{props} is all possible combinations of 10 equally spaced proportions (from 0 to 1), subject to the constraint that they add to 1. By itself, this method is not going to find the MLE unless the MLE happens to be one of the grid points. However, using the option \code{surf = TRUE} with the \code{SANN} or \code{mcmc} methods can improve efficiency by initiating the search at the grid point nearest the MLE.

}
\value{
A matrix with 10 named columns for each individual, containing estimated genomic proportions, ancestry indices, and the log-likelihood:
  \item{p11 }{Proportion of markers homozygous for lineage 1 alleles}
  \item{p22 }{Proportion of markers homozygous for lineage 2 alleles}
  \item{p33 }{Proportion of markers homozygous for lineage 3 alleles}
  \item{p12 }{Proportion of markers heterozygous for lineage 1 and 2 alleles}
  \item{p13 }{Proportion of markers heterozygous for lineage 1 and 3 alleles}
  \item{p23 }{Proportion of markers heterozygous for lineage 2 and 3 alleles}
  \item{S1 }{Lineage 1 ancestry index: proportion of alleles derived from parental lineage 1}
  \item{S2 }{Lineage 2 ancestry index: proportion of alleles derived from parental lineage 2}
  \item{S3 }{Lineage 3 ancestry index: proportion of alleles derived from parental lineage 3}
  \item{logLik }{log-likelihood of the genomic proportions given the individual marker data}

}
\references{
Fitzpatrick, B. M. 2012. Estimating ancestry and heterozygosity of hybrids using molecular markers. BMC Evolutionary Biology 12:131. http://www.biomedcentral.com/1471-2148/12/131
}
\author{
Ben Fitzpatrick
}

\seealso{
\code{\link{HIC3}} calculates closed-form maximum likelihood estimates for diagnostic marker data. \code{\link{thirdclass}} and \code{\link{HItest3}} evaluate simple classification of three-way hybrids into parental, F1, F2, and backcross categories. For conventional two-way hybrid zone analyses, see \code{\link{HIC}, \link{HIest}, \link{HIclass}, \link{HItest}}.
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

P <- data.frame(Locus=rep(1:12,each=3),allele=rep(1:3,12),
	P1=rep(c(1,0,0),12),P2=rep(c(0,1,0),12),P3=rep(c(0,0,1),12))

mle.o <- threeway(G,P,surf=FALSE,iterations=99)
mle.c <- HIC3(G,P)

# compare the optimization (mle.o) to the closed-form (mle.c): 
# 99 iterations is not enough to converge on the known true values. 
# Try setting surf=TRUE and/or increasing iterations. 
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
