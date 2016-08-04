\name{rwg}
\alias{rwg}
\title{James et al., (1984) agreement index for single item measures}
\description{
This function calculates the within group agreement measure rwg for single
item measures as described in James, Demaree and Wolf (1984).  The rwg is calculated
as rwg = 1-(Observed Group Variance/Expected Random Variance).  James et al. (1984)
recommend truncating the Observed Group Variance to the Expected Random Variance in
cases where the Observed Group Variance was larger than the Expected Random Variance.
This truncation results in an rwg value of 0 (no agreement) for groups with large
variances.
}
\usage{
rwg(x, grpid, ranvar=2)
}
\arguments{
 \item{x}{A vector representing the construct on which to estimate agreement.}
 \item{grpid}{A vector identifying the groups from which x originated.}
 \item{ranvar}{The random variance to which actual group variances are compared.
The value of 2 represents the variance from a rectangular
distribution in the case where there are 5 response options (e.g.,
Strongly Disagree, Disagree, Neither, Agree, Strongly Agree).
In cases where there are not 5 response options, the rectangular
distribution is estimated using the formula
\eqn{\mathtt{ranvar}=(A^{2}-1)/12}{ranvar=(A^2-1)/12} where
A is the number of response options.  While the rectangular distribution is
widely used, other random values may be more appropriate.}
}
\value{
	\item{grpid}{The group identifier}
	\item{rwg}{The rwg estimate for the group}
	\item{gsize}{The group size}
}     

\author{Paul Bliese
\email{paul.bliese@moore.sc.edu}
}

\references{
Bliese, P. D. (2000).  Within-group agreement, non-independence, and reliability:
 Implications for data aggregation and analysis.  In K. J. Klein & S. W. Kozlowski (Eds.),
Multilevel Theory, Research, and Methods in Organizations (pp. 349-381). San Francisco, CA:  Jossey-Bass, Inc.

James, L.R., Demaree, R.G., & Wolf, G. (1984).  Estimating within-group
interrater reliability with and without response bias.  Journal of Applied 
Psychology, 69, 85-98.
}
\seealso{
 \code{\link{ad.m}}
 \code{\link{rwg.j}}
 \code{\link{rwg.sim}}
 \code{\link{rgr.agree}}
 \code{\link{rwg.j.lindell}}
}
\examples{
data(lq2002)
RWGOUT<-rwg(lq2002$LEAD,lq2002$COMPID)
RWGOUT[1:10,]
summary(RWGOUT)

}
\keyword{attribute}