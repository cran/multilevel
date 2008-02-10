\name{sobel.lme}
\alias{sobel.lme}
\title{Estimate Sobel's (1982) Test for Mediation in Two-Level lme Model}
\description{Estimate Sobel's (1982) indirect test for mediation in a two-level
linear mixed effects (lme) model.  The function provides
an estimate of the magnitude of the indirect effect, Sobel's first-order estimate of the standard error
associated with the indirect effect, and the corresponding z-value.  The estimates are based
upon three models as detailed on page 84 of MacKinnon, Lockwood, Hoffman, West and Sheets (2002). By
estimating the three models in a two-level lme framework, one can account for the effects of 
hierarchical nesting in the data, and thereby obtain more accurate estimates of standard errors (see
Bliese & Hanges, 2004).}

\usage{
sobel.lme(pred, med, out, grpid)
}

\arguments{
  \item{pred}{The predictor or independent variable.}
  \item{med}{The mediating variable.}
  \item{out}{The outcome or dependent variable.}
  \item{grpid}{A variable identifying the groups in which pred, med and out are nested.}
}

\value{
	\item{Model.1}{Summary of coefficients from Model 1 of MacKinnon et al.}
	\item{Model.2}{Summary of coefficients from Model 2 of MacKinnon et al.}
	\item{Model.3}{Summary of coefficients from Model 3 of MacKinnon et al.}
	\item{Indirect.Effect}{Estimate of indirect mediating effect.}
	\item{SE}{Sobel's (1982) Standard Error estimate.}
	\item{z.value}{The estimated z-value.}
	\item{N}{The number of observations used in model estimation.}
	\item{n.grps}{The number of groups used in model estimation.}
}
\references{Bliese, P. D. & Hanges, P. J. (2004).  Being both too liberal and too conservative:  The
perils of treating grouped data as though they were independent.  Organizational Research
Methods, 7, 400-417.

MacKinnon, D. P., Lockwood, C. M., Hoffman, J. M., West, S. G., Sheets, V. (2002). A
comparison of methods to test mediation and other intervening variable effects.  Psychological
Methods, 7, 83-104.

Sobel, M. E., (1982).  Asymptotic confidence intervals for indirect effects in structural equation
models.  In S. Leinhardt (Ed.), Sociological Methodology 1982 (pp. 290-312). Washington, DC:  American
Sociological Association.}

\author{ Paul Bliese
 \email{paul.bliese@us.army.mil}}

\seealso{  \code{\link{sobel}}
}

\examples{
data(bh1996)
library(nlme)

#A small but significant indirect effect indicates leadership mediates
#the relationship between work hours and well-being.

with(bh1996, sobel.lme(pred=HRS,med=LEAD,out=WBEING,grpid=GRP))
}
\keyword{htest}
