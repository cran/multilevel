\name{make.univ}
\alias{make.univ}
\title{Convert data from multivariate to univariate form}
\description{
Longitudinal data is typically stored in multivariate form.  In multivariate form, each row contains
data from one subject even when the subject provides data at multiple time periods (eg. OUTCOME1,
OUTCOME2, OUTCOME3).  In the case of repeated measures designs and growth modeling it
is necessary to convert the data to univariate or stacked form where each row represents a single
time indexed by a variable such as TIME.  In the univariate case, a subject
has multiple rows of data.
}
\usage{
make.univ(x,dvs)
}
\arguments{
 \item{x}{A dataframe in typical multivariate form.}
 \item{dvs}{A subset dataframe of x containing the repeated measures columns.  Note that dvs must be time-ordered
from Time 1 to Time N for this function to work properly.}
}
\value{
Returns a new dataframe in univariate (i.e., stacked) form with a new variable named TIME representing the time intervals
and a new variable named MULTDV representing the stacked multilple dependent variables.  The TIME variable begins with 0.
}     
\author{
 Paul Bliese
 \email{paul.bliese@us.army.mil}
}
\references{
Bliese, P. D., & Ployhart, R. E. (2002). Growth modeling using random coefficient models:  Model building, testing
and illustrations.  Organizational Research Methods, 5, 362-387.
}

\examples{
data(univbct)  #a data set already in univariate or stacked form for job satisfaction
TEMP<-univbct[3*1:495,c(22,1:17)]  #converting it back to multivariate form
TEMP2<-make.univ(x=TEMP,dvs=TEMP[,c(10,13,16)])  #transforming it to univariate form again
}
\keyword{manip}