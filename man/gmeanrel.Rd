\name{gmeanrel}
\alias{gmeanrel}
\title{Group Mean Reliability from an lme model (nlme package)}
\description{Calculates the group-mean reliability from a linear mixed effects (lme) model.
If group sizes are identical, the group-mean reliabilty estimate equals the ICC(2) estimate from an
ANOVA model. When group sizes differ, however, a group-mean reliability estimate
is calculated for each group based on the Spearman-Brown formula, the overall ICC, and
group size for each group.
}
\usage{
gmeanrel(object)
}
\arguments{
 \item{object}{A Linear Mixed Effect (lme) object.}
}
\value{
	\item{ICC}{Intraclass Correlation Coefficient}
	\item{Group}{A vector containing all the group names.}
	\item{GrpSize}{A vector containing all the group sizes.}
	\item{MeanRel}{A vector containing the group-mean reliability estimate for each group.}
}
     
\author{
 Paul Bliese
 \email{pdbliese@gmail.com}
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
 \code{\link{ICC1}}
 \code{\link{ICC2}}
 \code{\link{lme}}
}
\examples{
data(bh1996)
library(nlme)
tmod<-lme(WBEING~1,random=~1|GRP,data=bh1996)
gmeanrel(tmod)
}
\keyword{attribute}