\name{rtoz}
\alias{rtoz}
\title{Conducts an r to z transformation}
\description{
Transforms a correlation (r) to a z variate using the formula
provided on page 53 of Cohen & Cohen (1983).  The formula is z=.5*((log(1+r))-(log(1-r)))
where r is the correlation.
}
\usage{
rtoz(rvalue)
}
\arguments{
 \item{rvalue}{The correlation to be z transformed.}
}
\value{The z transformation.
}     
\author{
 Paul Bliese
 \email{pdbliese@gmail.com}
}
\references{
Cohen, J. & Cohen, P. (1983).  Applied multiple regression/correlation analysis
for the behavioral sciences (2nd Ed.).  Hillsdale, NJ:  Lawrence Erlbaum Associates.
}
\seealso{
 \code{\link{cordif}}
}
\examples{
rtoz(.84)
}
\keyword{htest}