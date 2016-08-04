\name{boot.icc}
\alias{boot.icc}
\title{Bootstrap ICC values in 2-level data}
\description{
Implements a 2-level bootstrap.  The bootstrap first draws a sample of level-2 units with replacement,
and in a second stage draws a sample of level-1 observations with replacement from the level-2 units.
Following each bootstrap replication, the Intraclass Correlation Coefficient 1 is estimated using
the lme function.
}
\usage{
boot.icc(x, grpid, nboot, aov.est=FALSE)
}
\arguments{
 \item{x}{A vector representing the variable upon which to estimate the ICC values.}
 \item{grpid}{A vector representing the level-2 unit identifier.}
 \item{nboot}{The number of bootstrap iterations.  Computational demands underlying
a 2-level bootstrap are heavy, so the examples use 100; however, the number of interations
should generally be 10,000.}
 \item{aov.est}{An option to estimate the ICC using aov.}
}
\value{Provides ICC(1) estimates for each bootstrap draw.}

\author{
Paul Bliese
\email{paul.bliese@moore.sc.edu}
}
\references{
Bliese, P. D. (2000).  Within-group agreement, non-independence, and reliability: 
Implications for data aggregation and analysis.  In K. J. Klein & S. W. Kozlowski (Eds.),
Multilevel Theory, Research, and Methods in Organizations (pp. 349-381). San Francisco, CA:  Jossey-Bass, Inc.}

\seealso{
 \code{\link{ICC1}}
 \code{\link{ICC2}}

}

\examples{
\dontrun{
data(bh1996)
ICC.OUT<-boot.icc(bh1996$WBEING,bh1996$GRP,100)
quantile(ICC.OUT,c(.025,.975))
}
}
\keyword{attribute}