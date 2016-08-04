\name{sherifdat}
\docType{data}
\alias{sherifdat}
\title{Sherif (1935) group data from 3 person teams}
\description{
This data set contains estimates of movement length (in inches) of a light in a
completely dark room.  Eight groups of three individuals provided three estimates
for a total of 72 observations.  In four of the groups, participants first made 
estimates alone prior to providing estimates as a group.  In the other four groups
participants started as groups. Lang and Bliese (forthcoming) used these data to
illustrate how variance functions in mixed-effects models (lme) could be used to
test whether groups displayed consensus emergence. Data were obtained from
https://brocku.ca/MeadProject/Sherif/Sherif_1935a/Sherif_1935a_3.html
}
\usage{data(sherifdat)}
\format{A dataframe with 5 columns and 72 observations
\tabular{llll}{
[,1] \tab person  \tab numeric  \tab Participant ID within a group\cr
[,2] \tab time    \tab numeric  \tab Measurment Occasion\cr
[,3] \tab group   \tab numeric  \tab Group Identifier \cr
[,4] \tab y       \tab numeric  \tab Estimate of movement length in inches \cr
[,4] \tab condition   \tab numeric  \tab Experimental Condition for either starting individually (1) or as a group (0)\cr
}
}
\references{
Sherif, M. (1935). A study of some social factors in perception: Chapter 3. Archives of Psychology, 27, 23- 46.

https://brocku.ca/MeadProject/Sherif/Sherif_1935a/Sherif_1935a_3.html

Lang, J. W. B., & Bliese, P. D. (forthcoming). A Temporal Perspective on Emergence: Using
3-level Mixed Effects Models to Track Consensus Emergence in Groups.
}
\keyword{datasets}