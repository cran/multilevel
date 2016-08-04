\name{quantile.agree.sim}
\alias{quantile.agree.sim}
\title{S3 method for class 'agree.sim'}
\description{
This function provides a concise quantile summary of objects created using
the functions rwg.sim and rwg.j.sim.  The simulation functions for rwg and
rwg.j return a limited number of estimated values.  Consequently, the normal
quantile methods are biased.  The quantile methods incorporated in this function
produce unbiased estimates. 
}
\usage{
\method{quantile}{agree.sim}(x,confint,\dots)
}
\arguments{
	 \item{x}{An object of class 'agree.sim'.}
	 \item{confint}{The confidence intervals to return.  The values of 0.95 and 0.99
	   return the approximate 95th and 99th percentile confidence intervals (p=.05 and
p=.01).}
	 \item{\dots}{Optional arguments.  None used.}
}

\value{A dataframe with two columns.  The first column contains the quantile value
and the second contains the estimate based on the object.
}

\author{
Paul Bliese
\email{paul.bliese@moore.sc.edu}
}

\seealso{
 \code{\link{rwg.sim}}
 \code{\link{rwg.j.sim}}
}
\examples{
#An example from Dunlap et al. (2003).  The estimate from Dunlap et al.
#Table 2 is 0.53
RWG.OUT<-rwg.sim(gsize=10,nresp=5,nrep=1000)
quantile(RWG.OUT, c(.95,.99))
}
\keyword{programming}