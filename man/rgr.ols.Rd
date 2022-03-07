\name{rgr.ols}
\alias{rgr.ols}
\title{Random Group Resampling OLS Regression}

\description{Uses Random Group Resampling (RGR) within an Ordinary Least Square (OLS)
framework to contrast actual group results with pseudo group results. This specific
function performs an RGR on an OLS hierarchical OLS model with two predictors as
in Bliese & Halverson (2002). To run this analysis on data with more predictors,
the function would have to be modified.}

\usage{
rgr.ols(xdat1,xdat2,ydata,grpid,nreps)
}
\arguments{
   \item{xdat1}{The first predictor.}
   \item{xdat2}{The second predictor.}
   \item{ydata}{The outcome.}
   \item{grpid}{The group identifier.}
   \item{nreps}{The number of pseudo groups to create.}
}
\value{A matrix containing mean squares.  Each row provides mean square values for a single pseudo
group iteration}

\references{Bliese, P. D., & Halverson, R. R. (2002).  Using random group resampling in multilevel research.  Leadership Quarterly, 13, 53-68.}

\author{Paul Bliese
 \email{pdbliese@gmail.com}}


\seealso{\code{\link{mix.data}}}

\examples{
data(lq2002)
RGROUT<-rgr.ols(lq2002$LEAD,lq2002$TSIG,lq2002$HOSTILE,lq2002$COMPID,100)

#Compare values to those reported on p.62 in Bliese & Halverson (2002)
summary(RGROUT)  
}
\keyword{attribute}
