\name{rmv.blanks}
\alias{rmv.blanks}
\title{Remove blanks spaces from non-numeric variables imported from SPSS dataframes}
\description{When large SPSS datasets are imported into R, non-numeric fields frequently
have numerous blank spaces prior to the text.  The blank spaces make it difficult to
summarize non-numeric text.  The function is applied to an entire dataframe and removes
the blank spaces.
}
\usage{
rmv.blanks(object)
}
\arguments{
  \item{object}{Typically a dataframe created from an imported SPSS file.}
}
\value{Returns a new dataframe without preceeding }
\author{ Paul Bliese
 \email{paul.bliese@moore.sc.edu}
}

\seealso{  \code{\link[foreign]{read.spss}}
}
\examples{
   \dontrun{library(foreign)
   mydata<-read.spss(file.choose(),to.data.frame=T,use.value.labels=F)
   mydata<-rmv.blanks(mydata)}
}
\keyword{manip}
