\name{cordif.dep}
\alias{cordif.dep}
\title{Estimate whether two dependent correlations differ}
\description{
Tests for statistical differences between two dependent
correlations using the formula provided on page 56 of Cohen & Cohen (1983).
The function returns a t-value, the DF, and the p-value.
}
\usage{
cordif.dep(r.x1y,r.x2y,r.x1x2,n)
}
\arguments{
 \item{r.x1y}{The correlation between x1 and y where y is typically the outcome variable.}
 \item{r.x2y}{The correlation between x2 and y where y is typically the outcome variable.}
 \item{r.x1x2}{The correlation between x1 and x2 (the correlation between the two predictors).}
 \item{n}{The sample size.}
}
\value{
Returns three values. A t-value, DF and p-value.
}     
\author{
 Paul Bliese
 \email{pdbliese@gmail.com}
}
\references{
Cohen, J. & Cohen, P. (1983).  Applied multiple regression/correlation analysis
for the behavioral sciences (2nd Ed.).  Hillsdale, nJ:  Lawrence Erlbaum Associates.
}
\seealso{
 \code{\link{cordif}}
}
\examples{
cordif.dep(r.x1y=.30,r.x2y=.60,r.x1x2=.10,n=305)
}
\keyword{htest}