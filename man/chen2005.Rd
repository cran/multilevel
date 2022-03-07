\name{chen2005}
\docType{data}
\alias{chen2005}
\title{Data from Chen (2005)}
\description{Contains the complete data used in Chen (2005). Chen (2005) examined
newcomer adaptation in 65 project teams. The level of analysis was the team. In the
study, team leaders assessed the initial team performance (TMPRF) at time 1
and then assessed newcomer performance over three additional time points (NCPRF.T1,
NCPRF.T2, NCPRF.T3). Initial team expectations (TMEXP) and initial newcomer
empowerment (NCEMP) were also assessed and modeled, but were not analyzed as repeated
measures. To specify Table 2 model in Chen (2005), these data need to be converted
to univariate or stacked form (see the make.univ function). Using the default values
of make.univ and creating a dataframe called chen2005.univ, the specific lme model for
Table 2 in Chen (2005) is:

lme(MULTDV~NCEMP*TIME+TMEXP*TIME+TMPRF*TIME,random=~TIME|ID,chen2005.univ)
}
\usage{data(chen2005)}
\format{A data frame with 7 columns and 65 team-level observations
\tabular{llll}{
[,1] \tab ID       \tab numeric  \tab Team Identifier\cr
[,2] \tab TMPRF    \tab numeric  \tab Initial Team Performance (time 1 in article)\cr
[,3] \tab TMEXP    \tab numeric  \tab Team Expectations (time 1 in article)\cr
[,4] \tab NCEMP    \tab numeric  \tab Initial Newcomer Empowerment(time 2 in article)\cr
[,5] \tab NCPRF.T1 \tab numeric  \tab Newcomer Performance Time 1 (time 2 in article)\cr
[,6] \tab NCPRF.T2 \tab numeric  \tab Newcomer Performance Time 2 (time 3 in article)\cr
[,7] \tab NCPRF.T3 \tab numeric  \tab Newcomer Performance Time 3 (time 4 in article)\cr
}
}
\references{
Chen, G.(2005). Newcomer adaptation in teams:  Multilevel antecedents and outcomes.
Academy of Management Journal, 48, 101-116.
}
\keyword{datasets}