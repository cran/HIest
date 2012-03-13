\name{HIest}
\alias{HIest}
\title{
Find the joint maximum likelihood estimates of ancestry and interclass heterozygosity in a sample of hybrids.
}
\description{
\code{HIest} provides approaches to find maximum likelihood estimates of S and H, using the likelihood functions described by Fitzpatrick (2012).
}
\usage{
HIest(G, P, type, method = "SANN", iterations = 1000, Cscale = NULL, surf = FALSE,
 startgrid = 99, start = c(.5,.5), control = list(fnscale = -1, maxit = iterations))
}
\arguments{
  \item{G}{
A matrix or data frame of genetic marker data. Each column is a locus. For \code{type="dominant"} or \code{type="allele.count"}, there should be one row per individual. For \code{type="codominant"}, each individual is to be represented in consecutive rows (one for each allele).
}
  \item{P}{
A matrix or data frame with the following columns (order is important!): Locus name, Allele name, P1 allele frequency, P2 allele frequency. For \code{type="dominant"} or \code{type="allele.count"}, there should be one row per locus, giving the frequencies of the dominant or "1" allele. For \code{type="codominant"} there should be a separate row for each allele AND the Allele names should match the data in \code{G}.
}
  \item{type}{
A string representing the data type. The options are \code{"codominant"}, \code{"dominant"}, and \code{"allele.count"}.
}
  \item{method}{
  	Optimization method to search for maximum likelihood estimates of ancestry and heterozygosity. Alternatives are \code{"SANN"},  \code{"L-BFGS-B"}, \code{"surf"}, and \code{"mcmc"}. See details.
}
  \item{iterations}{
The desired number of MCMC steps to perform when \code{method="mcmc" or "SANN"}.
}
  \item{Cscale}{
An integer, controlling the the proposal distribution for \code{method = "SANN"} or \code{"mcmc"}. The default value is 100. Smaller values will cause the algorithm to search more broadly, but could make the search inefficient. See details.
}
  \item{surf}{Logical: Should the function find starting values by evaluating likelihoods on a grid?}
  \item{startgrid}{ Integer.
This controls the size of the grid for \code{method="surf"}. It is the same as the argument \code{size} in the function \code{\link{HIsurf}}.
}
  \item{start}{
A vector including the starting values of S and H.
}
  \item{control}{
A list of options to be passed to \code{control} in the \code{\link{optim}} function. Whatever else is chosen, be sure \code{fnscale} is negative to make \code{optim} search for a maximum rather than a minimum.
}
}
\details{
Given two ancestral species or parental populations (P1 and P2), the ancestry index (S) is the proportion of an individual's alleles descending from alleles in the P1 population and the interclass heterozygosity (H) is the proportion of an individual's loci that have one allele from each ancestral population (Lynch 1991). The likelihood functions are described in Fitzpatrick (2012). The likelihood functions take advantage of the correspondence between ancestry and heterozygosity and the genomic proportions: \eqn{p_{11}} = proportion of one's genome that is homozygous for alleles inherited from parental lineage 1, \eqn{p_{12} = H} = proportion one's genome that is heterozygous for alleles inhertied from each parental lineage, and \eqn{p_{22}} = proportion of one's genome that is homozygous for alleles inherited from parental lineage 2.

Currently, the function provides four methods for searching the likelihood surface. \code{method = "SANN"} is probably the best; it uses the general purpose optimization function \code{\link{optim}} with its simulated annealing algorithm. For this estimation problem, a custom proposal function is passed to the option \code{gr}. This proposal function draws new genomic proportions (\eqn{p_{11},p_{12},p_{22}}) from a three dimensional Dirichlet distribution centered on the old genomic proportions. The concentration of the proposal distribution is controlled by \code{Cscale}; the larger this value, the more the proposal distribution is concentrated near the current state.

\code{method = "mcmc"} uses a Markov-Chain Monte Carlo with Metropolis-Hastings sampling to explore the likelihood surface. It also uses the Dirichlet proposal distribution, and could be useful (with some modification of the code) for generating posterior distributions. \code{method = "SANN"} is probaby superior for simply finding the MLE.

\code{method = "L-BFGS-B"} also uses \code{\link{optim}}, but with a quasi-Newton likelihood search algorithm to look for the maximum likelihood. This method is relatively fast, but it can miss the MLE if it is near the edge of the sample space.

\code{"surf"}, finds all likelihoods on a grid defined by \code{startgrid} and chooses the maximum. This is not going to find the MLE unless the MLE happens to be one of the grid points. However, using the option \code{surf = TRUE} with the \code{SANN} or \code{mcmc} methods can improve efficiency by initiating the search at the grid point nearest the MLE.

}
\value{
A data frame with one row for each individual and three columns:
\item{S }{The maximim likelihood estimate of the ancestry index}
\item{H }{The maximum likelihood estimate of the interclass heterozygosity}
\item{logLik}{The maximum log-likelihood}
}
\references{
Fitzpatrick, B. M. 2008. Hybrid dysfunction: Population genetic and quantitative genetic perspectives. American Naturalist 171:491-198.

Fitzpatrick, B. M. 2012. Estimating ancestry and heterozygosity of hybrids using molecular markers.

Fitzpatrick, B. M., J. R. Johnson, D. K. Kump, H. B. Shaffer, J. J. Smith, and S. R. Voss. 2009. Rapid fixation of non-native alleles revealed by genome-wide SNP analysis of hybrid tiger salamanders. BMC Evolutionary Biology 9:176. http://www.biomedcentral.com/1471-2148/9/176

Lynch, M. 1991. The genetic interpretation of inbreeding depression and outbreeding depression. Evolution 45:622-629.
}
\author{
Ben Fitzpatrick
}

\seealso{
\code{\link{HIsurf}} for a likelihood surface, \code{\link{HIclass}} for likelihoods of early generation hybrid classes, \code{\link{HItest}} to compare the classification to the maximum likelihood, \code{\link{HILL}} for the basic likelihood function.
}
\examples{
##-- A random codominant data set of 5 individuals and 5 markers with three alleles each
L <- 10
P <- data.frame(Locus=rep(1:L,each=3),Allele=rep(1:3,L),
	P1=as.vector(rmultinom(L,10,c(.7,.2,.1)))/10, 
	P2=as.vector(rmultinom(L,10,c(.1,.2,.7)))/10)
G <- matrix(nrow=10,ncol=L)
for(i in 1:L){
	G[,i] <- sample(c(1,2,3),size=10,replace=TRUE,prob=rowMeans(P[P$Locus==i,3:4]))
	}
	
HI.cod <- HIest(G,P,type="codominant",iterations=99,surf=TRUE,startgrid=20)

# this is unlikely to converge on the MLE: increase iterations and/r startgrid.

# # optional plot
# plot(c(0,.5,1,0),c(0,1,0,0),type="l",xlab=expression(italic(S)),
	# ylab=expression(italic(H[I])),lwd=2,cex.lab=1.5,cex.axis=1.5,bty="n")
# points(HI.cod$S,HI.cod$H,cex=1.5,lwd=2)
# axis(1,labels=FALSE,lwd=2);axis(2,labels=FALSE,lwd=2)

# # other examples

# ##-- Make it into allele count data (count "3" alleles)
# P.c <- P[seq(from=3,to=dim(P)[2],by=3),]
# G.c <- matrix(nrow=5,ncol=L)
# for(i in 1:5){
	# G.c[i,] <- colSums(G[c(i*2-1,i*2),]==3)
	# }

# HI.ac <- HIest(G.c,P.c,type="allele.count",iterations=500,surf=TRUE,startgrid=50)

# ##-- Make it into dominant data where allele 3 is dominant
# G.d <- replace(G.c,G.c==2,1)

# HI.dom <- HIest(G.d,P.c,type="dominant",iterations=500,surf=TRUE,startgrid=50)

# ## -- A real dataset (Fitzpatrick et al. 2009)
# data(Bluestone)
# Bluestone <- replace(Bluestone,is.na(Bluestone),-9)
# # parental allele frequencies (assumed diagnostic)
# BS.P <- data.frame(Locus=names(Bluestone),Allele="BTS",P1=1,P2=0)

# # estimate ancestry and heterozygosity
# # BS.est <-HIest(Bluestone,BS.P,type="allele.count")
# ## shortcut for diagnostic markers
# BS.est <- HIC(Bluestone)

# # calculate likelihoods for early generation hybrid classes
# BS.class <- HIclass(Bluestone,BS.P,type="allele.count")

# # compare classification with maximum likelihood estimates
# BS.test <- HItest(BS.class,BS.est)

# table(BS.test$c1)
# # all 41 are TRUE, meaning the best classification is at least 2 log-likelihood units
# # better than the next best

# table(BS.test$c2)
# # 2 are TRUE, meaning the MLE S and H are within 2 log-likelihood units of the best
# # classification, i.e., the simple classification is rejected in all but 2 cases

# table(BS.test$Best.class,BS.test$c2)
# # individuals were classified as F2-like (class 3) or backcross to CTS (class 4), but
# # only two of the F2's were credible 

# BS.test[BS.test$c2,]
# # in only one case was the F2 classification a better fit (based on AIC) than the
# # continuous model.


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
