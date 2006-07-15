\name{quantile.disagree.sim}
\alias{quantile.disagree.sim}
\title{S3 method for class 'disagree.sim'}
\description{
This function provides a concise quantile summary of objects created using
the function ad.m.sim.  The simulation functions for the average deviation of
the mean (or median) return a limited number of estimated values.  Consequently, the normal
quantile methods are biased.  The quantile methods incorporated in this function
produce unbiased estimates. 
}
\usage{
\method{quantile}{disagree.sim}(x,confint,\dots)
}
\arguments{
	 \item{x}{An object of class 'disagree.sim'.}
	 \item{confint}{The confidence intervals to return.  The values of 0.05 and 0.01
return the approximate 5 percent and 1 percent confidence intervals.  Values equal to or
smaller than these values are significant (p=.05, p=.01).}
	 \item{\dots}{Optional arguments.  None used.}
}

\value{A dataframe with two columns.  The first column contains the quantile value
and the second contains the estimate based on the object.
}

\author{
Paul Bliese
\email{paul.bliese@us.army.mil}
}

\seealso{
 \code{\link{ad.m.sim}}
}
\examples{
#Example from Dunlap et al. (2003), Table 3.  The listed significance value (p=.05)
#for a group of size 5 with a 7-item response format is 0.64 or less 
SIMOUT<-ad.m.sim(gsize=5, nitems=1, nresp=7, itemcors=NULL, type="mean", nrep=1000) 
quantile(SIMOUT, c(.05,.01))
}
\keyword{programming}