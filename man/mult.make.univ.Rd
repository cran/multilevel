\name{mult.make.univ}
\alias{mult.make.univ}
\title{Convert two or more variables from multivariate to univariate form}
\description{
Longitudinal data is often stored in multivariate or wide form.  In multivariate form, each row contains
data from one subject, and repeated measures variables are indexed by different names (e.g., OUTCOME.T1,
OUTCOME.T2, OUTCOME.T3).  In the case of repeated measures designs and growth modeling, it
is necessary to convert the data to univariate or stacked form where each row represents one of the 
repeated measures indexed by a TIME variable and nested within subject.  In univariate form, each subject
has as many rows of data as there are time points. The make.univ function in the multilevel library will convert a
single item to univariate form while the mult.make.univ function converts two or more variables to 
univariate form.  The mult.make.univ function was developed by Patrick Downes at the University of Iowa, 
and was recommended for inclusion in the multilevel library in January of 2013. 
}
\usage{
mult.make.univ(x,dvlist,tname="TIME", outname="MULTDV")
}
\arguments{
 \item{x}{A dataframe in multivariate form.}
 \item{dvlist}{A list containing the repeated measures.  Note that each element of the
       list must be ordered from Time 1 to Time N for this function to work properly.}
 \item{tname}{An optional name for the new time variable.  Defaults to TIME.}
 \item{outname}{An optional name for the outcome variable name.  Defaults to MULTDV1 to MULTDV(N).}
}
\value{
Returns a dataframe in univariate (i.e., stacked) form with a TIME variable representing the repeated observations,
and new variables representing the time-indexed variables (MULTDV1, MULTDV2, etc.).  The TIME variable begins with 0.
}     
\author{
Patrick Downes
 \email{pat-downes@uiowa.edu} 
Paul Bliese
 \email{paul.bliese@moore.sc.edu}
}
\references{
Bliese, P. D., & Ployhart, R. E. (2002). Growth modeling using random coefficient models:  Model building, testing
and illustrations.  Organizational Research Methods, 5, 362-387.
}
\seealso{
 \code{\link{make.univ}}
}
\examples{
data(univbct) #a dataframe in univariate form for job sat
TEMP<-univbct[3*1:495,c(22,1:17)]  #convert back to multivariate form
names(TEMP) #use the column names to find the column numbers

#Create a list of DV's - each DV should have the same number of obs
dvlist <- list(c(10,13,16),c(11,14,17))
names(dvlist) <- c("JOBSAT","COMMIT") #names for univariate output

#Transform the data into univariate form with multiple level-1 variables
mldata <- mult.make.univ(x=TEMP,dvlist=dvlist)  
}
\keyword{manip}