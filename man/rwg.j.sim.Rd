\name{rwg.j.sim}
\alias{rwg.j.sim}
\title{Simulate rwg(j) values from a random null distribution}
\description{
This function is based on the work of Cohen, Doveh and Eick (2001) and Cohen, Doveh and Nahum-Shani
(2009).  The function draws data from a random uniform null distribution, and calculates the James,
Demaree and Wolf (1984) within group agreement measure rwg(j) for multiple item scales.  By
repeatedly drawing random samples, a distribution of the rwg(j) is generated.  The sampling
distribution can be used to calculate confidence intervals for different combinations of group
sizes and number of items (J).  Users provide the number of scale response options (A) and the
number of random samples.  By default, items (J) drawn in the simulation are independent
(non-correlated); however, an optional argument (itemcors) allows the user to specify a correlation
matrix with relationships among items.  Cohen et al. (2001) show that values of rwg(j) are primarily
a function of the number of items and the group size and are not strongly influenced by
correlations among items; nonetheless, assuming correlations among items is more realistic and
thereby is a preferred model (see Cohen et al., 2009).  If item correlations are provided, the MASS
library also needs to be attached. 
}
\usage{
rwg.j.sim(gsize, nitems, nresp, itemcors=NULL, nrep)
}
\arguments{
 \item{gsize}{Group size used in the rwg(j) simulation.}
 \item{nitems}{The number of items (J) in the multi-item scale on which to base the simulation.  If
itemcors are provided, this is an optional argument as nitems will be calculated from the
correlation matrix.}
 \item{nresp}{The number of response options for the J items in the simulation (e.g., there would be 5
response options if using Strongly Disagree, Disagree, Neither, Agree, Strongly Agree).}
 \item{itemcors}{An optional argument containing a correlation matrix with the item correlations.}
 \item{nrep}{The number of rwg(j) values to simulate.  This will generally be 10,000 or more, but
only 1,000 are used in the examples to increase the speed.}

}

\value{
	\item{rwg.j}{rwg(j) value from each of the nrep simulations.}
	\item{gsize}{Simulation group size.}
	\item{nresp}{Simulated number of response options.}
	\item{nitems}{Number of items in the multiple item scale. Either provided in the call or
calculated from the correlation matrix, if given.}
	\item{rwg.j.95}{95 percent confidence interval estimate associated with a p-value of .05.  Values
greater than or equal to the rwg.j.95 value are considered significant.}
}


\author{
Paul Bliese
\email{paul.bliese@moore.sc.edu}
}
\references{
Cohen, A., Doveh, E., & Nahum-Shani, I. (2009).  Testing agreement for multi-item scales
with the indices rwg(j) and adm(j).  Organizational Research Methods, 12, 148-164.

Cohen, A., Doveh, E., & Eick, U. (2001).  Statistical properties of the rwg(j)
index of agreement.  Psychological Methods, 6, 297-310.

James, L.R., Demaree, R.G., & Wolf, G. (1984).  Estimating within-group
interrater reliability with and without response bias.  Journal of Applied 
Psychology, 69, 85-98.
}

\seealso{
 \code{\link{rwg.j}}
 \code{\link{rwg}}
 \code{\link{rwg.sim}}
 \code{\link{rwg.j.lindell}}
 \code{\link{rgr.agree}}

}

\examples{
#An example assuming independent items
RWG.J.OUT<-rwg.j.sim(gsize=10,nitems=6,nresp=5,nrep=1000)
summary(RWG.J.OUT)
quantile(RWG.J.OUT, c(.95,.99))

#A more realistic example assuming correlated items.  The 
#estimate in Cohen et al. (2006) is .61.

data(lq2002)
library(MASS)
RWG.J.OUT<-rwg.j.sim(gsize=10,nresp=5,
  itemcors=cor(lq2002[,c("TSIG01","TSIG02","TSIG03")]),
  nrep=1000)
summary(RWG.J.OUT)
quantile(RWG.J.OUT,c(.95,.99))
}
\keyword{attribute}