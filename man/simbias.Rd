\name{simbias}
\alias{simbias}
\title{Simulate Standard Error Bias in Non-Independent Data}
\description{Non-independence due to groups is a common characteristic
of applied data.  In non-independent data, responses from members of the
same group are more similar to each other than would be expected by chance.  
Non-independence is typically measured using the Intraclass Correlation
Coefficient 1 or ICC(1).  When non-independent data is treated as though it is
independent, standard errors will be biased and power can decrease.  This simulation
allows one to estimate the bias and loss of statistical power that occurs when
non-independent data is treated as though it is independent.  The simulation contrasts
a simple Ordinary Least Squares (OLS) model that fails to account for non-independence
with a random coefficient model that accounts for non-independence.  The simulation
assumes that both the outcome (y) and the predictor (x) vary among individuals
in the same group.

}
\usage{
simbias(corr,gsize,ngrp,icc1x,icc1y,nrep)
}
\arguments{
	\item{corr}{The simulated true correlation between x and y.}
	\item{gsize}{The group size from which x and y are drawn.}
	\item{ngrp}{The number of groups.}
	\item{icc1x}{The simulated ICC(1) value for x.}
	\item{icc1y}{The simulated ICC(1) value for y.}
	\item{nrep}{The number of repetitions of simulated data sets.}
}
\value{
	\item{icc1.x}{Observed ICC(1) value for x in the simulation.}
	\item{icc1.y}{Observed ICC(1) value for y in the simulation.}
	\item{lme.coef}{Parameter estimate from the lme model.}
	\item{lme.se}{Standard error estimate from the lme model.}
	\item{lme.tvalue}{t-value from the lme model.}
	\item{lm.coef}{Parameter estimate from the linear model (OLS).}
	\item{lm.se}{Standard error estimate from the linear model (OLS).}
	\item{lm.tvalue}{t-value from the linear model (OLS).}
}

\author{
Paul Bliese
\email{paul.bliese@us.army.mil}
}

\references{
Bliese, P. D. & Hanges, P. J. (2004).  Being both too liberal and too conservative:  The perils of treating grouped data as
though they were independent.  Organizational Research Methods, 7, 400-417.

}
\seealso{
 \code{\link{ICC1}}
}
\examples{

library(nlme)
set.seed(15)
SIM.OUTPUT<-simbias(corr=.15,gsize=10,ngrp=50,icc1x=0.05,
            icc1y=0.35, nrep=100)
apply(SIM.OUTPUT,2,mean)
1-pnorm(1.96-3.39) #Power of the lme model (two-tailed, alpha=.05)
1-pnorm(1.96-2.95) #Power of the OLS model (two-tailed, alpha=.05)

}
\keyword{datagen}