\name{ICC1}
\alias{ICC1}
\title{Estimate Intraclass Correlation Coefficient 1 or ICC(1) from an aov model}
\description{Calculates the Intraclass Correlation Coefficient 1 or ICC(1)
from an ANOVA model. This value is equivalent to the ICC discussed in the
random coefficient modeling literature, and represents the amount of 
individual-level variance that can be "explained" by group membership.
}
\usage{
ICC1(object)
}
\arguments{
  \item{object}{An ANOVA (aov) object from an one-way analysis of variance.}
}
\value{Provides an estimate of ICC(1) for the sample.}
\references{

Bliese, P. D. (2000).  Within-group agreement, non-independence, and
reliability:  Implications for data aggregation and Analysis.  In K. J. 
Klein & S. W. Kozlowski (Eds.), Multilevel Theory, Research, and Methods 
in Organizations (pp. 349-381). San Francisco, CA:  Jossey-Bass, Inc.

Bartko, J.J. (1976).  On various intraclass correlation reliability
coefficients.  Psychological Bulletin, 83, 762-765.}

\author{ Paul Bliese
 \email{pdbliese@gmail.com}
}

\seealso{  \code{\link{ICC2}}
 \code{\link{aov}}
 \code{\link{sim.icc}}
 \code{\link{sim.mlcor}}
}
\examples{
   data(bh1996)
   hrs.mod<-aov(HRS~as.factor(GRP),data=bh1996)
   ICC1(hrs.mod)
}
\keyword{attribute}
