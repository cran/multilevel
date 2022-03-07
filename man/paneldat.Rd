\name{paneldat}
\docType{data}
\alias{paneldat}
\title{Firm-Level Panel Data from Donald J. (DJ) Schepker}
\description{
Contains a random sample of firm-level data from 196 firms for the years 2007 to 2011
based on data from COMPUSTAT and MSCI. The data are balanced in that each firm provided five years of data. The data contains
a time variable and time varying covariates for the discontinuous growth model along
with within-firm predictors related to the composition of boards of directors. The discontinuity was
indexed to 2009 as the first year following the great recession. In addition to
illustrating variants of growth models, the panel data is useful for contrasting
econometric fixed-effect models and mixed-effect models (Bliese et al., 2020).
}

\usage{data(paneldat)}
\format{A data frame with 14 columns and 960 observations
\tabular{llll}{
[,1] \tab companyid       \tab numeric  \tab Firm ID\cr
[,2] \tab companyname    \tab numeric  \tab Firm Name\cr
[,3] \tab year    \tab numeric  \tab Calendar Year\cr
[,4] \tab time    \tab numeric  \tab Calendar Year minus 2007\cr
[,5] \tab trans \tab numeric  \tab Dummy Coded Variable Indexed to Transition\cr
[,6] \tab post \tab numeric  \tab Time Varying Covariate for Post Transition Slope\cr
[,8] \tab roani \tab numeric  \tab Return on Assets\cr
[,9] \tab boardindprop \tab numeric  \tab Yearly Proportion of Independent Members of the Firm's Board\cr
[,10] \tab dirageavg \tab numeric  \tab Yearly Average Age of the Firm's Directors\cr
[,11] \tab femaledirsprop \tab numeric  \tab Yearly Female Board Members Proportion\cr
[,12] \tab debtassets \tab numeric  \tab Yearly Debt to Assets Ratio\cr
[,13] \tab lnemp \tab numeric  \tab Yearly Natural Log of Number of Employees\cr
[,14] \tab ceotenure \tab numeric  \tab Yearly CEO Tenure\cr
}
}
\references{
Bliese, P. D., Schepker, D. J., Essman, S. M., & Ployhart, R. E. (2020).  Bridging
methodological divides between macro- and microresearch: Endogeneity and methods
for panel data. Journal of Management, 46, 70-99.
}
\keyword{datasets}