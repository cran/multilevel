\name{rgr.agree}
\alias{rgr.agree}
\title{Random Group Resampling for Within-group Agreement}
\description{
This function uses random group resampling (RGR) to estimate within group
agreement.  RGR agreement compares within group variances from
actual groups to within group variances from pseudo groups.  If variances
on the construct of interest in the actual groups are significantly smaller than variances from
pseudo groups, one can concludes that there is significant agreement. 
}
\usage{
rgr.agree(x, grpid, nrangrps)
}
\arguments{
 \item{x}{The vector representing the construct of interest.  This is the
variable about which one is interested in determining whether or not
group members agree.}
 \item{grpid}{A vector identifying the groups from which x originated (actual
group membership).}
 \item{nrangrps}{A number representing the number of random groups that one is
interested in generating.  Note that the number of random groups
created must be directly divisible by the number of actual groups
to ensure that group sizes of pseudo groups and actual groups are
identical.  The rgr.agree routine will generate the number of 
pseudo groups that most closely approximates nrangrps given the
group size characteristics of one's data.}
}
\value{An object of class 'rgr.agree' with the following components:

	\item{NRanGrp}{The number of random groups created.}
	\item{AvRGRVar}{The average within-group variance of the random groups.}
	\item{SDRGRVar}{The standard deviation of the random group variances used in the z-score estimate.}
	\item{zvalue}{The z-score for the difference between the actual group and random group variances.}
	\item{RGRVARS}{The random group variances.} 
}

\author{
Paul Bliese
\email{paul.bliese@us.army.mil}
}
\references{
Bliese, P. D., & Halverson, R. R. (2002).  Using random group resampling in multilevel research.  Leadership Quarterly, 13, 53-68.

Bliese, P.D., Halverson, R. R., & Rothberg, J. (2000).  Using random
group resampling (RGR) to estimate within-group agreement with examples
using the statistical language R.  Walter Reed Army Institute of Research.
}
\seealso{
 \code{\link{rwg}}
 \code{\link{rwg.j}}
}
\examples{

data(bh1996)
RGROUT<-rgr.agree(bh1996$HRS,bh1996$GRP,1000)
summary(RGROUT)
}
\keyword{attribute}