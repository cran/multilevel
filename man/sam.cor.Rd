\name{sam.cor}
\alias{sam.cor}
\title{Generate a Sample that Correlates a Given Vector}
\description{Generate a sample vector (y) with a known population correlation to a given vector (x).
The degree of correlation between x and y is determined by rho (the population
correlation). Observed sample correlations between x and y will vary around rho, but this variation
will decrease as the number of x observations increase.}

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

\author{Paul Bliese
 \email{pdbliese@gmail.com}}

\seealso{\code{\link{sim.mlcor}}}

\examples{
   data(bh1996)
   NEWVAR<-sam.cor(x=bh1996$LEAD,rho=.30)
   cor(bh1996$LEAD,NEWVAR)
}
\keyword{programming}
