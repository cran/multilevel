\name{cordif}
\alias{cordif}
\title{Estimate whether two independent correlations differ}
\description{
This function tests for statistical differences between two independent
correlations using the formula provided on page 54 of Cohen & Cohen (1983).
The function returns a z-score estimate.
}
\usage{
cordif(rvalue1,rvalue2,n1,n2)
}
\arguments{
 \item{rvalue1}{Correlation value from first sample.}
 \item{rvalue2}{Correlation value from second sample.}
 \item{n1}{The sample size of the first correlation.}
 \item{n2}{The sample size of the second correlation.}
}
\value{
Produces a single value, the z-score for the differences between the correlations.
}     
\author{
 Paul Bliese
 \email{paul.bliese@moore.sc.edu}
}
\references{
Cohen, J. & Cohen, P. (1983).  Applied multiple regression/correlation analysis
for the behavioral sciences (2nd Ed.).  Hillsdale, NJ:  Lawrence Erlbaum Associates.
}
\seealso{
 \code{\link{rtoz}}
 \code{\link{cordif.dep}}
}
\examples{
cordif(rvalue1=.51,rvalue2=.71,n1=123,n2=305)
}
\keyword{htest}