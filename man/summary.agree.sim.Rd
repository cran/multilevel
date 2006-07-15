\name{summary.agree.sim}
\alias{summary.agree.sim}
\title{S3 method for class 'agree.sim'}
\description{
This function provides a concise summary of objects created using
the functions rwg.sim and rwg.j.sim. 
}
\usage{
\method{summary}{agree.sim}(object,\dots)
}
\arguments{
	 \item{object}{An object of class 'agree.sim'.}
	 \item{\dots}{Optional additional arguments.  None used.}
}
\value{

	\item{A summary of all the output elements in the agree.sim class object.}
}

\author{
Paul Bliese
\email{paul.bliese@us.army.mil}
}

\seealso{
 \code{\link{rwg.sim}}
 \code{\link{rwg.j.sim}}
}
\examples{
#An example from Dunlap et al. (2003).  The estimate from Dunlap et al.'s (2003)
#Table 2 is 0.53
RWG.OUT<-rwg.sim(gsize=10,nresp=5,nrep=1000)
summary(RWG.OUT)
}
\keyword{programming}