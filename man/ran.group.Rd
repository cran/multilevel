\name{ran.group}
\alias{ran.group}
\title{Randomly mix grouped data and return function results}
\description{This function is called by rgr.agree (and potentially other functions).  The
ran.group function randomly mixes data and applies a function to the pseudo groups. Pseudo group IDs match 
real group IDs in terms of size.}
\usage{
ran.group(x,grpid,fun,...)
}
\arguments{
  \item{x}{A matrix or vector containing data to be randomly sorted.}
  \item{grpid}{A vector containing a group identifier.}
  \item{fun}{A function to be applied to the observations within each random group.}
  \item{...}{Additional arguments to fun.}
}
\value{A vector containing the results of applying the specified function to each random group.}

\references{Bliese, P. D., & Halverson, R. R. (2002).  Using random group resampling in multilevel research.  Leadership Quarterly, 13, 53-68.}

\author{
Paul Bliese
 \email{paul.bliese@us.army.mil}}

\seealso{\code{\link{rgr.agree}}}

\examples{
   data(bh1996)
   ran.group(bh1996$HRS,bh1996$GRP,mean)
}
\keyword{programming}
