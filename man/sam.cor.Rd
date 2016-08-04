\name{sam.cor}
\alias{sam.cor}
\title{Generate a Sample that Correlates with a Fixed Set of Observations}
\description{This function will generate a vector (y) with a known correlation to a given vector (x).
The degree of correlation between x and y is determined by the parameter rho (the population
correlation).  Observed sample correlations between x and y will vary around rho, but this variation
will decrease as the size of x increases.}

\usage{
sam.cor(x,rho)
}
\arguments{
  \item{x}{The given vector.}
  \item{rho}{Population correlation.}
}
\value{The function prints the sample correlation for the specific set of numbers generated.
	\item{y}{A vector of numbers correlated with x.}
}

\author{ Paul Bliese
 \email{paul.bliese@moore.sc.edu}}

\seealso{\code{\link{simbias}}}

\examples{
   data(bh1996)
   NEWVAR<-sam.cor(x=bh1996$LEAD,rho=.30)
   cor(bh1996$LEAD,NEWVAR)
}
\keyword{programming}
