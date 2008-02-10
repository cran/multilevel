\name{ad.m.sim}
\alias{ad.m.sim}
\title{Simulate significance of average deviation around mean or median}
\description{
This function uses procedures detailed in Dunlap, Burke, and Smith-Crowe (2003) and Cohen,
Doveh, and Nahum-Shani (2009) to estimate the significance of the average deviation of the
mean or median (AD.M).  Dunlap et al. proposed a strategy to use Monte Carlo techniques to
estimate the significane of single item AD.M measures.  Cohen et al., (2009) expanded these ideas
to cover multiple item scales, ADM(J), and account for correlations among items.  The ad.m.sim
function is flexible and covers single item or multiple item measures.  In the case of multiple
item measures, correlations among items can be included (preferred method) or excluded.  If item
correlations are provided, the MASS library must also be attached.  In the Monte Carlo
simulations conducted by both Dunlap et al. (2003) and Cohen et al., (2009), 100,000 repetitions
were used.  In practice, it will require considerable time to perform 100,000 repititions and
in most cases 10,000 should suffice.  The examples use 1,000 repetitions simply to speed up the
process.
}
\usage{
ad.m.sim(gsize, nitems=1, nresp, itemcors=NULL, type="mean",nrep)
}
\arguments{
 \item{gsize}{Simulated group size.}
 \item{nitems}{Number of items to simulate.  The default is 1 for single item measures.  If
itemcors are provided, this is an optional argument as nitems will be calculated from the
correlation matrix, thus it is only necessary for multiple item scales where no correlation
matrix is provided.}
 \item{nresp}{The number of response options on the items.  For instance, nresp would equal
5 for a 5-point response option of strongly disagree, disagree, neither, agree,
strongly agree.}
 \item{itemcors}{An optional matrix providing correlations among items.}
 \item{type}{A character string with either "mean" or "median".}
 \item{nrep}{The number of simulation repetitions.}
}

\value{
	\item{ad.m}{Simulated estimates of AD.M values for each of the nrep runs.}
	\item{gsize}{Simulated group size.}
	\item{nresp}{Simulated number of response options.}
	\item{nitems}{Number of items. Either provided in the call (default of 1) or
calculated from the correlation matrix, if given.}
	\item{ad.m.05}{Estimated p=.05 value.  Observed values equal to or smaller
than this value are considered significant.}
	\item{pract.sig}{Estimate of practical significance calculated as nresp/6 (see ad.m).}
}

\author{
Paul Bliese
\email{paul.bliese@us.army.mil}
}
\references{
Cohen, A., Doveh, E., & Nahum-Shani, I. (2009).  Testing agreement for multi-item scales
with the indices rwg(j) and adm(j).  Organizational Research Methods, 12, 148-164.
 
Dunlap, W. P., Burke, M. J., & Smith-Crowe, K. (2003).  Accurate tests of statistical
significance for rwg and average deviation interrater agreement indices.  Journal of
Applied Psychology, 88, 356-362.
}

\seealso{
 \code{\link{ad.m}}
 \code{\link{rgr.agree}}
 \code{\link{rwg.sim}}
 \code{\link{rwg.j.sim}}
}

\examples{
#Example from Dunlap et al. (2003), Table 3.  The listed significance
#value (p=.05) for a group of size 5 with a 7-item response format is
#0.64 or less
 
SIMOUT<-ad.m.sim(gsize=5, nitems=1, nresp=7, itemcors=NULL,
        type="mean", nrep=1000) 
summary(SIMOUT)

#Example with a multiple item scale basing item correlations on observed
#correlations among 11 leadership items in the lq2002 data set.  Estimate
#in Cohen et al., (2009) is 0.99

library(MASS)
data(lq2002)
SIMOUT<-ad.m.sim(gsize=10, nresp=5, itemcors=cor(lq2002[,3:13]), 
        type="mean", nrep=1000) 
summary(SIMOUT)
quantile(SIMOUT,c(.05,.10))
}

\keyword{attribute}