\name{rwg.j.lindell}
\alias{rwg.j.lindell}
\title{Lindell et al. r*wg(j) agreement index for multi-item scales}
\description{
Calculates the Lindell et al r*wg(j) within-group agreement index for multiple
item measures. The r*wg(j) is similar to the James, Demaree and Wolf (1984) rwg and rwg(j) indices.  The
r*wg(j) index is calculated by taking the average item variability as the Observed Group Variance, and using the
average item variability in the numerator of the rwg formula (rwg=1-(Observed Group Variance/
Expected Random Variance)). In practice, this means that the r*wg(j) does not increase as the number of items in the 
scale increases as does the rwg(j). Additionally, the r*wg(j) allows Observed Group Variances to be
larger than Expected Random Variances. In practice this means that r*wg(j) values can be negative.
}
\usage{
rwg.j.lindell(x, grpid, ranvar=2)
}
\arguments{
 \item{x}{A matrix representing the scale of interest upon which one is interested
in estimating agreement. Each column of the matrix represents a separate scale item,
and each row represents an individual respondent.}
 \item{grpid}{A vector identifying the groups from which x originated.}
 \item{ranvar}{The random variance to which actual group variances are compared.
The value of 2 represents the variance from a rectangular
distribution in the case where there are 5 response options (e.g.,
Strongly Disagree, Disagree, Neither, Agree, Strongly Agree).
In cases where there are not 5 response options, the rectangular
distribution is estimated using the formula
\eqn{\mathtt{ranvar}=(A^{2}-1)/12}{ranvar=(A^2-1)/12}
where
A is the number of response options.  Note that one is not limited
to the rectangular distribution; rather, one can include any 
appropriate random value for ranvar.}
}

\value{
	\item{grpid}{The group identifier}
	\item{rwg.lindell}{The r*wg(j) estimate for the group}
	\item{gsize}{The group size}
}


\author{
Paul Bliese
\email{pdbliese@gmail.com}
}

\references{
James, L.R., Demaree, R.G., & Wolf, G. (1984).  Estimating within-group
interrater reliability with and without response bias.  Journal of Applied 
Psychology, 69, 85-98.

Lindell, M. K. & Brandt, C. J. (1999).  Assessing interrater agreement on 
the job relevance of a test: A comparison of CVI, T, rWG(J), and r*WG(J) indexes.
Journal of Applied Psychology, 84, 640-647.
}

\seealso{
\code{\link{ad.m}}
\code{\link{awg}}
\code{\link{rwg}}
\code{\link{rwg.j}}
\code{\link{rgr.agree}}
}

\examples{
data(lq2002)
RWGOUT<-rwg.j.lindell(lq2002[,3:13],lq2002$COMPID)
RWGOUT[1:10,]
summary(RWGOUT)
}

\keyword{attribute}