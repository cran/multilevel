\name{sim.icc}
\alias{sim.icc}
\title{Simulate 2-level ICC(1) values with and without level-1 correlation}
\description{ICC(1) values play an important role influencing the form of
relationships among variables in nested data. This simulation
allows one to create data with known ICC(1) values. Multiple variables can be
created both with and without level-1 correlation.
}
\usage{
sim.icc(gsize, ngrp, icc1,nitems=1,item.cor=FALSE)
}
\arguments{
	\item{gsize}{The simulated group size.}
	\item{ngrp}{The simulated number of groups.}
	\item{icc1}{The simulated ICC(1) value.}
	\item{nitems}{The number of items (vectors) to simulate.}
	\item{item.cor}{An option to create level-1 correlation among items. 
                        Provided as a value between 0 and 1. If used, nitems
                        must be larger than 1.}
}
\value{
	\item{GRP}{The grouping designator.}
	\item{VAR1}{The simulated value. Multiple numbered columns if nitems>1}
}

\author{
Paul Bliese
\email{paul.bliese@moore.sc.edu}
}

\references{
Bliese, P. D. (2000).  Within-group agreement, non-independence, and reliability:
Implications for data aggregation and analysis.  In K. J. Klein & S. W. Kozlowski (Eds.),
Multilevel Theory, Research, and Methods in Organizations (pp. 349-381). San
Francisco, CA:  Jossey-Bass, Inc.
}
\seealso{
 \code{\link{ICC1}}
}
\examples{
 \dontrun{
set.seed(1535324)
ICC.SIM<-sim.icc(gsize=10,ngrp=100,icc1=.15)
ICC1(aov(VAR1~as.factor(GRP), ICC.SIM))

# 4 items with no level-1 correlation
set.seed(15324)
ICC.SIM<-sim.icc(gsize=10,ngrp=100,icc1=.15,nitems=4) #items with no level-1 correlation
mult.icc(ICC.SIM[,2:5],ICC.SIM$GRP)
with(ICC.SIM,waba(VAR1,VAR2,GRP))$Cov.Theorem  #Examine CorrW 

# 4 items with a level-1 correlation of .30
set.seed(15324)
ICC.SIM<-sim.icc(gsize=10,ngrp=100,icc1=.15,nitems=4, item.cor=.3) #.30 level-1 item correlations
mult.icc(ICC.SIM[,2:5],ICC.SIM$GRP)
with(ICC.SIM,waba(VAR1,VAR2,GRP))$Cov.Theorem  #Examine CorrW 
}
}
\keyword{datagen}