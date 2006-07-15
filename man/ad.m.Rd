\name{ad.m}
\alias{ad.m}
\title{Average deviation around mean or median}
\description{
This function calculates the average deviation around the mean or median as a measure
of within-group agreement as proposed by Burke, Finkelstein and Dusig (1999).  A basic
rule for interpreting whether or not the results display practically significant levels
of agreement is whether the AD value is smaller than A/6 where A represents the number of
response options.  For instance, A would be 5 on a five-point response option format ranging
from strongly disagree, disagree, neither, agree, strongly agree (see Dunlap, Burke & Smith-Crowe,
2003).  To estimate statistical significance see the ad.m.sim function and help files.
}
\usage{
ad.m(x, grpid, type="mean")
}
\arguments{
 \item{x}{A vector representing a single item or a matrix representing a scale of interest
upon which one is interested in estimating agreement.  If a matrix, each column of the matrix
represents a separate scale item,and each row represents an individual respondent.}
 \item{grpid}{A vector identifying the groups from which x originated.}
 \item{type}{A character string for either the mean or median.}
}

\value{
	\item{grpid}{The group identifier}
	\item{AD.M}{The average deviation around the mean or median for each group}
	\item{gsize}{Group size}
}

}
\author{
Paul Bliese
\email{paul.bliese@us.army.mil}
}
\references{
Burke, M. J., Finkelstein, L. M., & Dusig, M. S. (1999).  On average deviation indices
for estimating interrater agreement.  Organizational Research Methods, 2, 49-68.

Dunlap, W. P., Burke, M. J., & Smith-Crowe, K. (2003).  Accurate tests of statistical
significance for rwg and average deviation interrater agreement indices.  Journal of
Applied Psychology, 88, 356-362.
}

\seealso{
 \code{\link{ad.m.sim}}
 \code{\link{rwg}}
 \code{\link{rwg.j}}
 \code{\link{rgr.agree}}
 \code{\link{rwg.sim}}
 \code{\link{rwg.j.sim}}
}

\examples{
data(bhr2000)

#Examples for multiple item scales
AD.VAL<-ad.m(bhr2000[,2:12],bhr2000$GRP)
AD.VAL[1:5,]
summary(AD.VAL)
summary(ad.m(bhr2000[,2:12],bhr2000$GRP,type="median"))

#Example for single item measure
summary(ad.m(bhr2000$HRS,bhr2000$GRP))
}

\keyword{attribute}