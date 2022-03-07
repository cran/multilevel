\name{summary.rgr.waba}
\alias{summary.rgr.waba}
\title{S3 method for class 'rgr.waba'}
\description{Provides a concise summary of objects created using
the function rgr.waba. 
}
\usage{
\method{summary}{rgr.waba}(object,\dots)
}
\arguments{
	 \item{object}{An object of class 'rgr.waba'.}
	 \item{\dots}{Optional additional arguments. None used.}
}
\value{A dataframe containing summary statistics in the form of number of repetitions (NRep),
Mean and Standard Deviations (SD) for each parameter in the rgr.waba model.
}

\author{
Paul Bliese
\email{pdbliese@gmail.com}
}

\seealso{
 \code{\link{rgr.waba}}
}
\examples{
data(bh1996)

#estimate the actual group model
waba(bh1996$HRS,bh1996$WBEING,bh1996$GRP)

#create 100 pseudo group runs and summarize results
RWABA<-rgr.waba(bh1996$HRS,bh1996$WBEING,bh1996$GRP,100)
summary(RWABA)  
}
\keyword{programming}