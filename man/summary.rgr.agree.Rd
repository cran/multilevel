\name{summary.rgr.agree}
\alias{summary.rgr.agree}
\title{S3 method for class 'rgr.agree'}
\description{Provides a concise summary of objects created using
the function rgr.agree. 
}
\usage{
\method{summary}{rgr.agree}(object,\dots)
}
\arguments{
	 \item{object}{An object of class 'rgr.agree'.}
	 \item{\dots}{Optional additional arguments.  None used.}
}
\value{

	\item{Summary Statistics for Random and Real Groups}{Number of random groups, Average random group variance,
Standard Deviation of random group variance, Actual group variance, z-value}
	\item{Lower Confidence Intervals (one-tailed)}{Lower confidence intervals based on sorted random group variances.}
	\item{Upper Confidence Intervals (one-Tailed)}{Upper confidence intervals based on sorted random group variances.} 
}

\author{
Paul Bliese
\email{pdbliese@gmail.com}
}

\seealso{
 \code{\link{rgr.agree}}
}
\examples{
# Example with a small number of replications (500). Increase in actual use.
data(bh1996)
RGROUT<-rgr.agree(bh1996$HRS,bh1996$GRP,500)
summary(RGROUT)
}
\keyword{programming}