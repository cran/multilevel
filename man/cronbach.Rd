\name{cronbach}
\alias{cronbach}
\title{Estimate Cronbach's Alpha}
\description{This function calculates the Cronbach's Alpha estimate of reliability for a multi-item scale.}
\usage{
cronbach(items)
}
\arguments{
  \item{items}{An matrix or data frame where each column represents an item in a multi-item scale.}
}
\value{
	\item{Alpha}{Estimate of Cronbach's Alpha.}
	\item{N}{The number of observations on which the Alpha was estimated.}
}
\references{Cronbach L. J. (1951) Coefficient Alpha and
the internal structure of tests. Psychometrika, 16,297-334}

\author{ Paul Bliese
 \email{paul.bliese@moore.sc.edu}}


\seealso{\code{\link{cronbach}}}

\examples{
   data(bhr2000)
   cronbach(bhr2000[,2:11])
}
\keyword{attribute}
