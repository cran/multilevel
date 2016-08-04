\name{rgr.OLS}
\alias{rgr.OLS}
\title{Random Group Resampling OLS Regression}

\description{This function uses Random Group Resampling (RGR) within an Ordinary Least Square (OLS)
framework to allow one to contrast actual group results with pseudo group results.  The number of columns
in the output matrix of the function (OUT) has to correspond to the number of mean squares you want in the output
which in turn is a function of the number of predictors. This specific
function does RGR on an OLS hierarchical OLS model with two predictors as in Bliese & Halverson (2002). To run this
analysis on data with more predictors, the function will have to be modified.}

\usage{
rgr.OLS(xdat1,xdat2,ydata,grpid,nreps)
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

\author{ Paul Bliese
 \email{paul.bliese@moore.sc.edu}}


\seealso{\code{\link{mix.data}}}

\examples{
data(lq2002)
RGROUT<-rgr.OLS(lq2002$LEAD,lq2002$TSIG,lq2002$HOSTILE,lq2002$COMPID,100)

#Compare values to those reported on p.62 in Bliese & Halverson (2002)
summary(RGROUT)  
}
\keyword{attribute}
