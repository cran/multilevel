\name{make.univ}
\alias{make.univ}
\title{Convert data from multivariate to univariate form}
\description{
Longitudinal data is often stored in multivariate or wide form.  In multivariate form, each row contains
data from one subject, and repeated measures variables are indexed by different names (e.g., OUTCOME.T1,
OUTCOME.T2, OUTCOME.T3).  In the case of repeated measures designs and growth modeling, it
is necessary to convert the data to univariate or stacked form where each row represents one of the 
repeated measures indexed by a TIME variable and nested within subject.  In univariate form, each subject
has as many rows of data as there are time points.
}
\usage{
make.univ(x,dvs,tname="TIME", outname="MULTDV")
}
\arguments{
 \item{x}{A dataframe in multivariate form.}
 \item{dvs}{A subset dataframe of x containing the repeated measures columns.  Note that the repeated measures
            must be ordered from Time 1 to Time N for this function to work properly.}
 \item{tname}{An optional name for the new time variable.  Defaults to TIME.}
 \item{outname}{An optional name for the outcome variable name.  Defaults to MULTDV.}
}
\value{
Returns a dataframe in univariate (i.e., stacked) form with a TIME variable representing the repeated observations,
and a variable named MULTDV representing the time-indexed variable.  The TIME variable begins with 0.
}     
\author{
 Paul Bliese
 \email{paul.bliese@us.army.mil}
}
\references{
Bliese, P. D., & Ployhart, R. E. (2002). Growth modeling using random coefficient models:  Model building, testing
and illustrations.  Organizational Research Methods, 5, 362-387.
}
\seealso{
 \code{\link{mult.make.univ}}
}
\examples{
data(univbct)  #a dataframe in univariate form for job satisfaction
TEMP<-univbct[3*1:495,c(22,1:17)]  #convert back to multivariate form

#Transform data to univariate form
TEMP2<-make.univ(x=TEMP,dvs=TEMP[,c(10,13,16)]) 

#Same as above, but renaming repeated variable
TEMP3<-make.univ(x=TEMP,dvs=TEMP[,c(10,13,16)],outname="JOBSAT")  
}
\keyword{manip}