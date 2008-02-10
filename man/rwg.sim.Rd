\name{rwg.sim}
\alias{rwg.sim}
\title{Simulate rwg values from a random null distribution}
\description{
This function is based on the work of Dunlap, Burke & Smith-Crowe (2003).  The function draws data
from a random uniform null distribution, and calculates the within group agreement measure
rwg for single item measures as described in James, Demaree & Wolf (1984).  By repeatedly
drawing random samples, a distribution of the rwg is generated.  The sampling
distribution can be used to calculate confidence intervals for different combinations of group
sizes and number of response options (A).
}
\usage{
rwg.sim(gsize, nresp, nrep)
}
\arguments{
 \item{gsize}{Group size upon which to base the rwg simulation.}
 \item{nresp}{The number of response options (e.g., there would be 5
response options if using Strongly Disagree, Disagree, Neither, Agree, Strongly Agree).}
 \item{nrep}{The number of rwg values to simulate.  This will generally be 10,000 or more,
although the examples use nrep of 1000 to make the calculations fast.}

}

\value{
	\item{rwg}{rwg value from each simulation.}
	\item{gsize}{Group size used in the rwg simulation.}
	\item{nresp}{Simulated number of response options.}
	\item{nitems}{Will always be 1 for an rwg estimate.}
	\item{rwg.95}{Estimated 95 percent confidence interval. Values
greater than or equal to rwg.95 are considered significant, p<.05.}

}

\author{
Paul Bliese
\email{paul.bliese@us.army.mil}
}
\references{
Cohen, A., Doveh, E., & Eick, U. (2001).  Statistical properties of the rwg(j)
index of agreement.  Psychological Methods, 6, 297-310.

Dunlap, W. P., Burke, M. J., & Smith-Crowe, K. (2003).  Accurate tests of statistical
significance for rwg and average deviation interrater agreement indices.  Journal of
Applied Psychology, 88, 356-362.

James, L.R., Demaree, R.G., & Wolf, G. (1984).  Estimating within-group
interrater reliability with and without response bias.  Journal of Applied 
Psychology, 69, 85-98.
}

\seealso{
 \code{\link{ad.m}}
 \code{\link{rwg.j}}
 \code{\link{rwg}}
 \code{\link{rwg.j.sim}}
 \code{\link{rgr.agree}}

}

\examples{
#An example from Dunlap et al. (2003).  The estimate from Dunlap
#et al. Table 2 is 0.53 (p=.05)
RWG.OUT<-rwg.sim(gsize=10,nresp=5,nrep=1000)
summary(RWG.OUT)
quantile(RWG.OUT, c(.95,.99))
}
\keyword{attribute}