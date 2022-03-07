\name{sobel}
\alias{sobel}
\title{Estimate Sobel's (1982) Test for Mediation}
\description{Estimate Sobel's (1982) indirect test for mediation. The function provides
an estimate of the magnitude of the indirect effect, Sobel's first-order estimate of the standard error
associated with the indirect effect, and the corresponding z-value. The estimates are based
upon three models as detailed on page 84 of MacKinnon, Lockwood, Hoffman, West and Sheets (2002).}
\usage{
sobel(pred,med,out)
}
\arguments{
  \item{pred}{The predictor or independent variable (X).}
  \item{med}{The mediating variable (M).}
  \item{out}{The outcome or dependent variable (Y).}
}

\value{
	\item{Mod1: Y~X}{A summary of coefficients from Model 1 of MacKinnon et al., (2002).}
	\item{Mod2: Y~X+M}{A summary of coefficients from Model 2 of MacKinnon et al., (2002).}
	\item{Mod3: M~X}{A summary of coefficients from Model 3 of MacKinnon et al., (2002).}
	\item{Indirect.Effect}{The estimate of the indirect mediating effect.}
	\item{SE}{Sobel's (1982) Standard Error estimate.}
	\item{z.value}{The estimated z-value.}
	\item{N}{The number of observations used in model estimation.}
}
\references{MacKinnon, D. P., Lockwood, C. M., Hoffman, J. M., West, S. G., Sheets, V. (2002). A
comparison of methods to test mediation and other intervening variable effects.  Psychological
Methods, 7, 83-104.

Sobel, M. E., (1982).  Asymptotic confidence intervals for indirect effects in structural equation
models.  In S. Leinhardt (Ed.), Sociological Methodology 1982 (pp. 290-312). Washington, DC:  American
Sociological Association.}

\author{ Paul Bliese
 \email{pdbliese@gmail.com}}

\examples{
data(bh1996)

#A small but significant indirect effect indicates leadership mediates
#the relationship between work hours and well-being.
sobel(pred=bh1996$HRS,med=bh1996$LEAD,out=bh1996$WBEING)
}
\keyword{htest}
