\name{summary.disagree.sim}
\alias{summary.disagree.sim}
\title{S3 method for class 'disagree.sim'}
\description{
This function provides a concise summary of objects created using
the function ad.m.sim. 
}
\usage{
\method{summary}{disagree.sim}(object,\dots)
}
\arguments{
	 \item{object}{An object of class 'disagree.sim'.}
	 \item{\dots}{Optional additional arguments.  None used.}
}
\value{

	\item{A summary of all the output elements in the disagree.sim class object.}
}

\author{
Paul Bliese
\email{paul.bliese@us.army.mil}
}

\seealso{
 \code{\link{ad.m.sim}}
}
\examples{
#Example from Dunlap et al. (2003), Table 3.  The listed significance value
#for a group of size 5 with a 7-item response format is 0.64 or less 
SIMOUT<-ad.m.sim(gsize=5, nitems=1, nresp=7, itemcors=NULL, type="mean", nrep=1000) 
summary(SIMOUT)
}
\keyword{programming}