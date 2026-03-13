\name{correct.gcor}
\alias{correct.gcor}
\title{Correct group-level correlations for group-mean reliability}
\description{
The magnitude of the observed correlation between two group means is partially a function
of the group-mean reliabilities of the two means (Bliese, 1998; 2000). To compare
group-mean correlations across studies, it is important to correct group-mean correlations
for the reliability of the group means. As shown by Maltarich, Bliese & Carpenter (2026), 
the correction should only be applied to the incremental difference between the within-group 
and the between-group correlation to avoid over correcting. 

By performing the correction on the incremental difference between the within-group and
Between-group correlations, it is also possible to adjust the group-mean correlation downwards.
A downward correction is important because group-mean correlations can be both larger and
smaller than the within-group correlation.

The correct.gcor function performs a correction on the incremental component of 
the group-mean correlation and returns an estimate of the corrected group-mean 
correlation. The function requires some combination of (a) an estimate of the
ICC(2) for both variables, (b) an estimate of the ICC(1) and group size for 
both variables, or (c) a combination of ICC(1) and average group size for one 
variable and the ICC(2) for the other. The ICC(2) is calculated from the ICC(1) 
and average group size so providing ICC(2) should is preferred if available, but
results should be similar regardless of which values are provided.     
}
\usage{
correct.gcor(w.corr,b.corr,icc1.x, icc2.x=NA, icc1.y, icc2.y=NA, gsize)
}
\arguments{
 \item{w.corr}{The value of the within-group (or raw) correlation.}
 \item{b.corr}{The value of the between-group correlation.}
 \item{icc1.x}{The ICC(1) of the x-variable if the ICC(2) value is unavailable.}
 \item{icc2.x}{The ICC(2) value of the x-variable (preferred if available).}
 \item{icc1.y}{The ICC(1) of the y-variable if the ICC(2) value is unavailable.}
 \item{icc2.y}{The ICC(2) value of the y-variable (preferred if available).}
 \item{gsize}{Average group size required if ICC(1) values are used.}
}
\value{Provides an estimate of the group-mean correlation corrected for
group-mean reliability.}

\author{
Paul Bliese
\email{pdbliese@gmail.com}
}
\references{
Bliese, P. D. (1998). Group size, ICC values, and group-level correlations: A 
simulation.  Organizational Research Methods, 1, 355-373.

Bliese, P. D. (2000). Within-group agreement, non-independence, and reliability: 
Implications for data aggregation and analysis. In K. J. Klein & S. W. Kozlowski (Eds.),
Multilevel Theory, Research, and Methods in Organizations (pp. 349-381). 
San Francisco, CA: Jossey-Bass, Inc.

Li, G., Rubenstein, A. L., Lin, W., Wang, M., & Chen, X. (2018). The curvilinear effect 
of benevolent leadership on team performance: The mediating role of team
action processes and the moderating role of team commitment. Personnel Psychology,
71:369-397.

Maltarich, M., Bliese, P. D., & Carpenter, M. (2026). Relationships between shared 
group properties: Theory, measurement, estimation and adjustment. Personnel Psychology.
DOI: https://doi.org/10.1111/peps.70022 
}

\seealso{
 \code{\link{ICC1}}
 \code{\link{ICC2}}
 \code{\link{sim.icc}}
 \code{\link{sim.mlcor}}
}

\examples{
# The example below is taken from Li et al., (2018). The study reports
# a group-mean correlation of .25 and a within correlation of .35 between 
# benevolent leadership (X) and team commitment (Y). This is a case where the 
# group-mean correlation is smaller than the within-correlation. Li and 
# colleagues also provide ICC(1) and ICC(2) values for both X and Y along with
# average group size. Data are from study 2 in Li et al. (2018). The ICC(1) for
# X (benevolent leadership) is .26 and the ICC(2) is .59. The ICC(1) for Y (team
# commitment) is .21 and the ICC(2) is .52. Average team size was 4.13. The
# corrected value for the group-mean correlation is estimated to be .17 rather
# than .25.

# First example based on ICC(2) values
correct.gcor(w.corr=.35,b.corr=.25,icc2.x=.59,icc2.y=.52)

# Second example based on ICC(1) values and group size
correct.gcor(w.corr=.35,b.corr=.25,icc1.x=.26,icc1.y=.21, gsize=4.13)

# Third example uses ICC(1) values for X and ICC(2) values for Y
correct.gcor(w.corr=.35,b.corr=.25,icc1.x=.26,icc2.y=.52, gsize=4.13)

}
\keyword{attribute}