\name{rgr.waba}
\alias{rgr.waba}
\title{Random Group Resampling of Covariance Theorem Decomposition}
\description{Performs the covariance theorem decomposition of a raw correlation in
situations where lower-level entities (individuals) are nested in higher-level groups
(see Dansereau, Alutto & Yammarino, 1984; Robinson, 1950). Builds upon previous work by
incorporating Random Group Resampling or RGR. RGR is used
to randomly assign individuals to pseudo groups and create a sampling
distributions of the covariance theorem components. The sampling distribution
provides a way to contrast actual group covariance components to pseudo group covariance
components.

Note that rgr.waba is computationally intensive. 
}
\usage{
rgr.waba(x, y, grpid, nrep)
}
\arguments{
 \item{x}{A vector representing one variable for the correlation.}
 \item{y}{A vector representing the other variable for the correlation.}
 \item{grpid}{A vector identifying the groups from which X and Y originated.}
 \item{nrep}{The number of times that the entire data set is reassigned to
pseudo groups}
}
\value{
Returns an object of class rgr.waba. The object is a list containing
each random run for each component of the covariance theorem.
}
\author{
Paul Bliese
\email{pdbliese@gmail.com}
}
\references{
Bliese, P. D. & Halverson, R. R. (1996).  Individual and nomothetic
models of job stress:  An examination of work hours, cohesion, and well-
being.  Journal of Applied Social Psychology, 26, 1171-1189.

Bliese, P. D., & Halverson, R. R. (2002).  Using random group resampling in multilevel
research.  Leadership Quarterly, 13, 53-68.

Dansereau, F., Alutto, J. A., & Yammarino, F. J. (1984).  Theory
testing in organizational behavior: The varient approach.  Englewood
Cliffs, NJ: Prentice-Hall.

Robinson, W. S. (1950).  Ecological correlations and the behavior of
individuals.  American Sociological Review, 15, 351-357.
}
\seealso{
 \code{\link{waba}}
}
\examples{
# This example is from Bliese & Halverson (1996). Notice that all of the
# values from the RGR analysis differ from the values based on actual
# group membership. Confidence intervals for individual components can
# be estimated using the quantile command. In practice, the nrep option
# should be more than 100

data(bh1996)

#estimate the actual group model
waba(bh1996$HRS,bh1996$WBEING,bh1996$GRP)

#create 100 pseudo group runs and summarize the model        
RWABA<-rgr.waba(bh1996$HRS,bh1996$WBEING,bh1996$GRP,nrep=100)  
summary(RWABA)

#Estimate 95th percentile confidence intervals (p=.05)                     
quantile(RWABA,c(.025,.975))       
}
\keyword{attribute}