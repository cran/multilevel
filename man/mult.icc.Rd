\name{mult.icc}
\alias{mult.icc}
\title{Multiple ICCs from a dataset}
\description{
Given a data frame and a group identifier, this function will estimate ICC(1)
and ICC(2) values for each column in the dataframe.  Note that this function depends upon the nlme 
package, and it only works with one level of nesting (e.g., students within schools).  The dependent
variable is assumed to be gaussian.
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
	\item{ICC1}{The intraclass correlation coefficient 1.}
	\item{ICC2}{The group mean reliability or intraclass correlation coefficient 2.}
}
\author{
 Paul Bliese
 \email{paul.bliese@us.army.mil}
}
\references{
Bliese, P. D. (2000).  Within-group agreement, non-independence, and
reliability:  Implications for data aggregation and Analysis.  In K. J. 
Klein & S. W. Kozlowski (Eds.), Multilevel Theory, Research, and Methods 
in Organizations (pp. 349-381). San Francisco, CA:  Jossey-Bass, Inc.

Bartko, J.J. (1976).  On various intraclass correlation reliability
coefficients.  Psychological Bulletin, 83, 762-765.
}
\seealso{
 \code{\link{ICC2}}
 \code{\link{ICC1}}
}
\examples{
library(nlme)
data(bh1996)
mult.icc(bh1996[,c("HRS","LEAD","COHES")],grpid=bh1996$GRP)

}
\keyword{attribute}