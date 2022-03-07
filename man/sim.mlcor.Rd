\name{sim.mlcor}
\alias{sim.mlcor}
\title{Simulate a multilevel correlation}
\description{In multilevel or hierarchical nested data, correlations at the within
and between levels often differ in magnitude and/or sign. For instance, Bliese
and Halverson (1996) showed that the within correlation between individual reports of work
hours and individual well-being was -.11. When these same data were mean-aggregated to
the group-level, the between correlation based on group means was -.71. A necessary, but not
sufficient, condition for differences across levels is a non-zero ICC1 value for
both variables (Bliese, 2000). This simulation creates a dataset with a group ID and an X and Y
variable for any combination of group size, number of groups, within and between correlations,
ICC1 values and reliability (alpha).  
}
\usage{
sim.mlcor(gsize,ngrp,gcor,wcor,icc1x,icc1y,alphax=1,alphay=1)
}
\arguments{
	\item{gsize}{The simulated group size.}
	\item{ngrp}{The simulated number of groups.}
	\item{gcor}{The simulated between group correlation.}
	\item{wcor}{The simulated within group correlation.}
	\item{icc1x}{The simulated ICC1 value for X.}
	\item{icc1y}{The simulated ICC1 value for Y.}
	\item{alphax}{The reliability (alpha) of X.}
	\item{alphay}{The reliability (alpha) of Y.}
}
\value{
	\item{GRP}{The grouping designator.}
        \item{X}{The simulated X value.}
	\item{Y}{The simulated Y value.}
}

\author{
Paul Bliese
\email{pdbliese@gmail.com}
}

\references{
Bliese, P. D. (2000).  Within-group agreement, non-independence, and reliability:
Implications for data aggregation and analysis.  In K. J. Klein & S. W. Kozlowski (Eds.),
Multilevel Theory, Research, and Methods in Organizations (pp. 349-381). San
Francisco, CA:  Jossey-Bass, Inc.

Bliese, P. D. & Halverson, R. R. (1996).  Individual and nomothetic models of job
stress:  An examination of work hours, cohesion, and well-being.  Journal of 
Applied Social Psychology, 26, 1171-1189.

Bliese, P. D., Maltarich, M. A., Hendricks, J. L., Hofmann, D. A., & Adler, A. B. (2019).
Improving the measurement of group-level constructs by optimizing between-group
differentiation. Journal of Applied Psychology, 104, 293-302.
}
\seealso{
 \code{\link{ICC1}}
 \code{\link{sim.icc}}
}
\examples{
\dontrun{
#
# Examine the multilevel properties of work hours and well-being
# in the bh1996 data
#
data(bh1996)
with(bh1996,waba(HRS,WBEING,GRP))
mult.icc(bh1996[,c("HRS","WBEING")],bh1996$GRP)
#
#Estimate true group-mean correlation by adding ICC2 adjusted incremental
#correlation back to within correlation. For value of -.8256
#
-0.110703+(-0.7121729--0.110703)/(sqrt(0.9171286*0.771756))
#
# Simulate data with same properties assuming no measurement error
#
set.seed(578323)
SIM.ML.COR<-sim.mlcor(gsize=75,ngrp=99,gcor=-.8256,wcor=-.1107,
                 icc1x=0.04338,icc1y=0.12924,alphax=1,alphay=1)
#
# Compare Simulated results to results (above) from bh1996
#
with(SIM.ML.COR,waba(X,Y,GRP))
mult.icc(SIM.ML.COR[,c("X","Y")],SIM.ML.COR$GRP)
}
}
\keyword{datagen}