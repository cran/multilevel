\name{item.total}
\alias{item.total}
\title{Item-total correlations}
\description{Calculates item-total correlations in multi-item scales.}
\usage{
item.total(items)
}
\arguments{
  \item{items}{A matrix or dataframe where each column represents an item in a multi-item scale.}
}
\value{
	\item{Variable}{Variable examined in the reliability analyses.}
	\item{Item.Total}{The item correlation with the mean of the other items.}
	\item{Alpha.Without}{The Cronbach Alpha reliability estimate of the scale without the variable.}
	\item{N}{The number of observations on which the analyses were calculated.}
}
\references{Cronbach L. J. (1951) Coefficient Alpha and
the internal structure of tests. Psychometrika, 16,297-334}

\author{ Paul Bliese
 \email{pdbliese@gmail.com}}

\seealso{\code{\link{cronbach}}}

\examples{
   data(bhr2000)
   item.total(bhr2000[,2:11])}

\keyword{reliability}
