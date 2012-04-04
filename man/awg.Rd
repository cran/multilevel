\name{awg}
\alias{awg}
\title{Brown and Hauenstein (2005) awg agreement index}
\description{
This function calculates the awg index proposed by Brown and Hauenstein (2005).  The awg
agreement index can be applied to either a single item vector or a multiple item matrix
representing a scale.  The awg is an analogue to Cohen's kappa. Brown
and Hauenstein (pages 177-178) recommend interpreting the awg similarly to how the rwg
(James et al., 1984) is commonly interpreted with values of .70 indicating acceptable
agreement; values between .60 and .69 as reasonable agreement, and values less than .60
as unacceptable levels of agreement.
}
\usage{
awg(x, grpid, range=c(1,5))
}
\arguments{
 \item{x}{A vector representing a single item or a matrix representing a scale of interest.
If a matrix, each column of the matrix represents a scale item, and each row represents
an individual respondent.}
 \item{grpid}{A vector identifying the groups from which x originated.}
 \item{range}{A vector with the lower and upper response options (e.g., c(1,5)) for
              a five-point scale from strongly disagree to strongly agree.}
}

\value{
	\item{grpid}{The group identifier.}
	\item{a.wg}{The awg estimate for each group.}
	\item{nitems}{The number of scale items when x is a matrix or dataframe representing
                      a multi-item scale.  This value is not returned when x is a vector.}
	\item{nraters}{The number of raters.  Given that the awg estimate is based on the
                       sample estimate of variance with N-1 in the denominator,
                       Brown and Hauenstein (2005) contend that awg can be estimated on as
                       few as A-1 raters where A represents the number of response options
                       specified by the range option (5 as the default).  Note that
                       in many situations nraters will correspond to group size.}

}

\author{
Paul Bliese
\email{paul.bliese@us.army.mil}
}

\references{
Brown, R. D. & Hauenstein, N. M. A. (2005). Interrater Agreement Reconsidered: An Alternative to the
rwg Indices.  Organizational Research Methods, 8, 165-184.

Wagner, S. M., Rau, C., & Lindemann, E. (2010).  Multiple informant methodology:  A critical review
and recommendations.  Sociological Methods and Research, 38, 582-618.
}

\seealso{
 \code{\link{rwg}}
 \code{\link{rwg.j}}
 \code{\link{ad.m}}
}

\examples{
data(lq2002)

#Examples for multiple item scales
awg.out<-awg(lq2002[,3:13],lq2002$COMPID,range=c(1,5))
summary(awg.out)

#Example for single item measure
awg.out<-awg(lq2002$LEAD05,lq2002$COMPID,range=c(1,5))
summary(awg.out)
}

\keyword{attribute}