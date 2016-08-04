\name{ICC2}
\alias{ICC2}
\title{Function to Estimate Intraclass Correlation Coefficient 2 or ICC(2) from an aov model}
\description{This function calculates the Intraclass Correlation Coefficient 2 or ICC(2)
from an ANOVA model.  This value represents the reliability of the group
means.
}
\usage{
ICC2(object)
}
\arguments{
  \item{object}{An ANOVA (aov) object from an one-way analysis of variance.}
}
\value{Provides an estimate of ICC(1) for the sample.}
\references{ Bliese, P. D. (2000).  Within-group agreement, non-independence, and
reliability:  Implications for data aggregation and Analysis.  In K. J. 
Klein & S. W. Kozlowski (Eds.), Multilevel Theory, Research, and Methods 
in Organizations (pp. 349-381). San Francisco, CA:  Jossey-Bass, Inc.

Bartko, J.J. (1976).  On various intraclass correlation reliability
coefficients.  Psychological Bulletin, 83, 762-765.}

\author{ Paul Bliese
 \email{paul.bliese@moore.sc.edu}
}

\seealso{  \code{\link{ICC1}}
 \code{\link{aov}}
}
\examples{
   data(bh1996)
   hrs.mod<-aov(HRS~as.factor(GRP),data=bh1996)
   ICC2(hrs.mod)
}
\keyword{attribute}
