\name{quantile.rgr.waba}
\alias{quantile.rgr.waba}
\title{S3 method for class 'rgr.waba'}
\description{
This function provides a concise quantile summary of objects created using
the function rgr.waba. 
}
\usage{
\method{quantile}{rgr.waba}(x,confint,\dots)
}
\arguments{
	 \item{x}{An object of class 'rgr.waba'.}
	 \item{confint}{The confidence intervals to return.  The values of 0.025 and 0.975
	   return the approximate 95th percentile confidence intervals
	   (p=.05).}
	 \item{\dots}{Optional arguments.  None used.}
}

\value{A dataframe containing the confidence intervals for each parameter in the rgr.waba model.
}

\author{
Paul Bliese
\email{paul.bliese@us.army.mil}
}

\seealso{
 \code{\link{rgr.waba}}
}
\examples{
data(bh1996)
waba(bh1996$HRS,bh1996$WBEING,bh1996$GRP)        #estimate the actual group model
RWABA<-rgr.waba(bh1996$HRS,bh1996$WBEING,bh1996$GRP,100)  #create 100 pseudo group runs
quantile(RWABA,confint=c(.025,.975))  
}
\keyword{programming}