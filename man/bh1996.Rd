\name{bh1996}
\docType{data}
\alias{bh1996}
\title{Data from Bliese and Halverson (1996)}
\description{
This dataset contains the complete data used in Bliese and Halverson
(1996).  The dataset contains 4 variables.  These variables are Cohesion (COHES),
Leadership Climate (LEAD), Well-Being (WBEING) and Work Hours (HRS).  Each of these
variables has two variants -- a group mean version that replicates each group mean
for every individual, and a within-group version where the group mean is subtracted
from each individual response.  The group mean version is designated with a
G. (e.g., G.HRS), and the within-group version is designated with a W. (e.g., W.HRS).
}
\usage{data(bh1996)}
\format{A data frame with 13 columns and 7,382 observations from 99 groups
\tabular{llll}{
[,1] \tab GRP     \tab numeric  \tab Group Identifier\cr
[,2] \tab COHES   \tab numeric  \tab Cohesion\cr
[,3] \tab G.COHES \tab numeric  \tab Average Group Cohesion\cr
[,4] \tab W.COHES \tab numeric  \tab Group-Mean Centered Cohesion\cr
[,5] \tab LEAD    \tab numeric  \tab Leadership \cr
[,6] \tab G.LEAD  \tab numeric  \tab Average Group Leadership\cr
[,7] \tab W.LEAD  \tab numeric  \tab Group-Mean Centered Leadership\cr
[,8] \tab HRS     \tab numeric  \tab Work Hours\cr
[,9] \tab G.HRS   \tab numeric  \tab Average Group Work Hours\cr
[,10] \tab W.HRS  \tab numeric  \tab Group-Mean Centered Work Hours\cr
[,11] \tab WBEING \tab numeric  \tab Well-Being\cr
[,12] \tab G.WBEING \tab numeric  \tab Average Group Well-Being\cr
[,13] \tab W.WBEING \tab numeric  \tab Group-Mean Centered Well-Being
}
}
\references{
Bliese, P. D. & Halverson, R. R. (1996).  Individual and nomothetic
models of job stress:  An examination of work hours, cohesion, and 
well-being.  Journal of Applied Social Psychology, 26, 1171-1189.
}
\keyword{datasets}
