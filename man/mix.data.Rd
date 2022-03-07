\name{mix.data}
\alias{mix.data}
\title{Randomly mix grouped data}
\description{Called by graph.ran.mean (and potentially other functions)
to randomly mix data and create new pseudo group ID variables. Pseudo group IDs match 
real group IDs in terms of size.}
\usage{
mix.data(x,grpid)
}
\arguments{
  \item{x}{A matrix or vector containing data to be randomly sorted.}
  \item{grpid}{A vector containing a group identifier.}
}
\value{
	\item{newid}{ A pseudo group ID.}
	\item{grpid}{ The real group ID.}
	\item{x}{ The values in x arranged as belonging to newid.}
}
\references{Bliese, P. D., & Halverson, R. R. (2002).  Using random group resampling in multilevel research.  Leadership Quarterly, 13, 53-68.}

\author{Paul Bliese
 \email{pdbliese@gmail.com}}


\seealso{\code{\link{graph.ran.mean}}}

\examples{
data(bh1996)
mix.data(x=bh1996[c(1:10,200:210,300:310),2:3],
         grpid=bh1996$GRP[c(1:10,200:210,300:310)])
}
\keyword{programming}
