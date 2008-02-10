\name{waba}
\alias{waba}
\title{Covariance Theoreom Decomposition of Bivariate Two-Level Correlation}
\description{
This routine performs the covariance theorem decomposition discussed by
Robinson (1950) and Dansereau, Alutto and Yammarino (1984).  Dansereau et al.
have labeled the variance decomposition Within-And-Between-Analysis II or
WABA II.  The program decomposes a raw correlation from a two-level nested design
into 6 components.  These components are (1) eta-between value for X, (2) eta-between
value for Y, (3) the group-size weighted group-mean correlation, (4) the within-eta value
for X, (5) the within-eta value for Y, and (6) the within-group correlation
between X and Y.  The last value represents the correlation between X and Y
after each variable has been group-mean centered.

The program is designed to automatically perform listwise deletion on missing values; 
consequently, users should pay attention to the diagnostic information (Number of Groups and
Number of Observations) provided as part of the output.

Note that Within-And-Between-Analysis proposed by Dansereau et al.
involves more than covariance theorem decomposition of correlations.
Specifically, WABA involves decision rules based on eta-values.  These
are not replicated in the R multilevel library because the eta based
decision rules have been shown to be highly related to group size (Bliese,
2000; Bliese & Halverson, 1998), a factor not accounted for in the complete
Within-And-Between-Analysis methodology.
}
\usage{
waba(x, y, grpid)
}
\arguments{
 \item{x}{A vector representing one variable in the correlation.}
 \item{y}{A vector representing the other variable in the correlation.}
 \item{grpid}{A vector identifying the groups from which x and y originated.}
}
\value{
Returns a list with three elements.
 \item{Cov.Theorem}{A 1 row dataframe with all of the elements of the covariance theorem.}
 \item{n.obs}{The number of observations used to calculate the covariance theorem.}
 \item{n.grps}{The number of groups in the data set.}
}
\author{
 Paul Bliese
 \email{paul.bliese@us.army.mil}
}
\references{
Bliese, P. D. (2000).  Within-group agreement, non-independence,
and reliability:  Implications for data aggregation and Analysis.  In
K. J. Klein & S. W. Kozlowski (Eds.), Multilevel Theory, Research,
and Methods in Organizations (pp. 349-381). San Francisco, CA:
Jossey-Bass, Inc.

Bliese, P. D., & Halverson, R. R. (1998).  Group size and
measures of group-level properties:  An examination of eta-squared and
ICC values.  Journal of Management, 24, 157-172.

Dansereau, F., Alutto, J. A., & Yammarino, F. J. (1984).  Theory
testing in organizational behavior: The varient approach.  Englewood
Cliffs, NJ: Prentice-Hall.

Robinson, W. S. (1950).  Ecological correlations and the behavior of
individuals.  American Sociological Review, 15, 351-357.
}

\seealso{
 \code{\link{rgr.waba}}
}

\examples{
data(bh1996)
waba(bh1996$HRS,bh1996$WBEING,bh1996$GRP)
}

\keyword{attribute}
