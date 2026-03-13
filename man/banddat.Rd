\name{banddat}
\docType{data}
\alias{banddat}
\title{Longitudinal Engagement Data from a University Band from Patrick Flynn}
\description{
Contains responses from 131 members of a University band who provided weekly engagement
data over 16 weeks. During the data collection, an unexpected football victory occurred
on the weekend preceding the twelfth data collection so TIME 11 was the first data
collection after the event. The unexpected football win can be used within a discontinuous
growth model (Bliese & Lang, 2016) to see how the trajectories of engagement changed in
response to this event. The TRANS and RECOV variables in the dataset are provided for
estimating a discontinuous growth model. The data also contain person-level measures
of the big five personality; a variable indicating the student's year at the university
(TENURE), and a variable indicating whether the student intended to remain a member of the
band the following year (RET.15).
}

\usage{data(banddat)}
\format{A data frame with 13 columns and 2096 observations
\tabular{llll}{
[,1] \tab ID       \tab numeric  \tab Person ID\cr
[,2] \tab ENGAGE    \tab numeric  \tab Self-Rating of Engagement with Band\cr
[,3] \tab TIME    \tab numeric  \tab Week\cr
[,4] \tab TRANS    \tab numeric  \tab Dummy Code Pre and Post Unexpected Football Win\cr
[,5] \tab RECOV \tab numeric  \tab Time indexed to Event for Post Transition Slope\cr
[,6] \tab FEMALE \tab numeric  \tab Female(1) vs. Male (0)\cr
[,7] \tab EMOSTAB \tab numeric  \tab Big Five Personality: Emotional Stability\cr
[,8] \tab EXTRO \tab numeric  \tab Big Five Personality: Extroversion\cr
[,9] \tab AGEE \tab numeric  \tab Big Five Personality: Agreeableness\cr
[,10] \tab CONSCIEN \tab numeric  \tab Big Five Personality: Conscientiousness\cr
[,11] \tab OPEN \tab numeric  \tab Big Five Personality: Openness\cr
[,12] \tab TENURE \tab numeric  \tab Freshman(1), Sophmore(2), Junior(3), Senior(4)\cr
[,13] \tab RET.15 \tab numeric  \tab Intent to Return: Low (1) to High (5)\cr
}
}
\references{
Bliese, P. D., & Lang, J. W. B. (2016).  Understanding relative and absolute change in
discontinuous growth models: Coding alternatives and implications for hypothesis testing.
Organizational Research Methods, 19, 562-592. https://doi.org/10.1177/1094428116633502
}
\keyword{datasets}