\name{graph.ran.mean}
\alias{graph.ran.mean}
\title{Graph Random Group versus Actual Group distributions}
\description{This function uses random group resampling (RGR) to create a distribution
of pseudo group means.  The pseudo group means are then contrasted with
actual group means to provide a visualization of the group-level properties
of the data. It is, in essense, a way of visualizing an Intraclass Correlation
Coefficient -- ICC(1).}
\usage{
graph.ran.mean(x, grpid, nreps, limits, graph=TRUE, bootci=FALSE)
}

\arguments{
 \item{x}{The vector representing the construct of interest.}
 \item{grpid}{A vector identifying the groups from which x originated.}
 \item{nreps}{A number representing the number of random groups to generate.
Because groups are created with the exact size
characteristics of the actual groups, the total number of pseudo groups created
may be calculated as nreps * Number Actual Groups. The value chosen for nreps only
affects the smoothness of the pseudo group line -- values greater than 25 should
provide sufficiently smooth lines.  Values of 1000 should be used if the bootci
option is TRUE.}
 \item{limits}{Controls the upper and lower limits of the y-axis on the plot.
The default is to set the limits at the 10th and 90th percentiles of the raw data.
This option only affects how the data is plotted.}
 \item{graph}{Controls whether or not a plot is returned.  If graph=FALSE,
the program returns a data frame with two columns.  The first column
contains the sorted means from the actual groups, and the second column contains
the sorted means from the pseudo groups.  This can be useful for plotting results
in other programs.}
 \item{bootci}{Determines whether approximate 95 percent confidence interval estimates
are calculated and plotted.  If bootci is TRUE, the nreps option should be 1000 or more.}
}

\value{Produces either a plot (graph=TRUE) or a data.frame (graph=FALSE)}

\references{Bliese, P. D., & Halverson, R. R. (2002).  Using random group resampling in multilevel research.  Leadership Quarterly, 13, 53-68.}

\author{ Paul Bliese
 \email{paul.bliese@us.army.mil}}


\seealso{
	\code{\link{ICC1}}
	\code{\link{mix.data}}
}

\examples{
   data(bh1996)
   graph.ran.mean(bh1996$HRS,bh1996$GRP,limits=c(8,16),nreps=1000, bootci=TRUE)
   GRAPH.DAT<-graph.ran.mean(bh1996$HRS,bh1996$GRP,limits=c(8,16),nreps=50,graph=FALSE)
}
\keyword{dplot}