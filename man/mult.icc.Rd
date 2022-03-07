\name{mult.icc}
\alias{mult.icc}
\title{Multiple ICCs from a dataset}
\description{Estimates ICC(1) and ICC(2) values for each column given a data frame and a group identifier.
Uses a mixed-effects estimate of the ICC, so ICC values cannot be negative. In cases
where ICC values are assumed to be zero or negative, the ANOVA-based formulas should be used (see the ICC1
and ICC2 functions). The mult.icc function only works with one level of nesting.
}
\usage{
mult.icc(x, grpid)
}
\arguments{
 \item{x}{A data frame containing the variables of interest in each column.}
 \item{grpid}{A vector identifying the groups from which the variables originated.}
}
\value{
	\item{Variable}{The variable name.}
	\item{ICC1}{Intraclass correlation coefficient 1.}
	\item{ICC2}{Group mean reliability or intraclass correlation coefficient 2.}
}
\author{
 Paul Bliese
 \email{pdbliese@gmail.com}
}
\references{
Bartko, J.J. (1976).  On various intraclass correlation reliability
coefficients.  Psychological Bulletin, 83, 762-765.

Bliese, P. D. (2000).  Within-group agreement, non-independence, and
reliability:  Implications for data aggregation and Analysis.  In K. J. 
Klein & S. W. Kozlowski (Eds.), Multilevel Theory, Research, and Methods 
in Organizations (pp. 349-381). San Francisco, CA:  Jossey-Bass, Inc.

Bliese, P. D., Maltarich, M. A., Hendricks, J. L., Hofmann, D. A., & Adler, A. B. (2019).
Improving the measurement of group-level constructs by optimizing between-group
differentiation. Journal of Applied Psychology, 104, 293-302.
}
\seealso{
 \code{\link{ICC2}}
 \code{\link{ICC1}}
 \code{\link{sim.icc}}
}
\examples{
library(nlme)
data(bh1996)
mult.icc(bh1996[,c("HRS","LEAD","COHES")],grpid=bh1996$GRP)

}
\keyword{attribute}