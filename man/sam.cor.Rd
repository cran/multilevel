\name{sam.cor}
\alias{sam.cor}
\title{Generate a Sample that Correlates with a Fixed Set of Observations}
\description{This function will generate a sample of numbers (y) that are correlated with a fixed set of observations (x).
The degree of correlation between x and y is determined by the parameter rho (the population correlation).  Note that
the sample correlation of x and y will vary around the true population correlation, but this variation will decrease
as the size of x increases.}

\usage{
sam.cor(x,rho)
}
\arguments{
  \item{x}{The fixed set of observations.}
  \item{rho}{The targeted population correlation value.}
}
\value{The function prints the sample correlation for the specific set of numbers generated.
	\item{y}{A vector of numbers correlated with x.}
}

\author{ Paul Bliese
 \email{paul.bliese@us.army.mil}}

\seealso{\code{\link{simbias}}}

\examples{
   data(bh1996)
   NEWVAR<-sam.cor(x=bh1996$LEAD,rho=.30)
   cor(bh1996$LEAD,NEWVAR)
}
\keyword{programming}
