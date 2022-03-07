\name{dgm.code}
\alias{dgm.code}
\title{Create Coding Matrix for Discontinuous Growth Model}
\description{Creates time-varying covariates for estimating a discontinuous
growth model (DGM). Creating time-varying covariates requires only simple recoding
of a time vector when data are balanced and the discontinuity event occurs
at the same time for each group. When data are not balanced and one or more
events occurs at different times for each group, coding the time-varying
covariates is complex. For instance, if the event of interest was employee turnover
in a store in a given month, it is likely that stores would differ on how many months 
of data were available and on the specific months when turnover occurred. With these
irregularly-timed events and different time windows for each store, 
it would be challenging to create time-varying covariates for the DGM specific to each
store's circumstances.}
\usage{
dgm.code(grp,time,event,n.events=FALSE,first.obs=FALSE)
}

\arguments{
 \item{grp}{A vector representing the group ID in the panel data. Each group ID
is repeated n times as represented in the panel.}
 \item{time}{A vector from 0:n-1 where n represents the number of observations
within each panel. Within each group, the time vector can vary in length.}
 \item{event}{A vector containing 1 for each time the event of interest for a
group occurs and a 0 in all other cases.}
 \item{n.events}{Controls how many events for which to create time-varying
covariates. In cases where some groups encounter numerous events, this argument
can be used to limit the coding to a small number (e.g., 2 or 3). Default is to create
as many time related covariates as occur in the group with the most events (which
is often not useful).}
 \item{first.obs}{Controls what to do if the first observation is an event. If
TRUE, then the first observation is changed to a zero and treated as a non-event. If FALSE,
the function provides a list of the groups where the event is the first observation.
The reason why the first observation is flagged is that it is not possible to estimate
a discontinuous growth model if the event occurs on the first observation.}
}
\value{Produces a data.frame with the columns
 
[, c("grp","time","event","trans1","post1","time.a","tot.events","event.first")]

If numerous events are coded,
the data.frame will contain more trans and post values (trans2, post2, etc.) corresponding to the maximum
number of events experienced by a group unless n.events is set to limit the number of
events. The data.frame must be merged back to the original data.frame for subsequent
analyses. See examples.}

\references{Bliese, P. D., & Lang, J. W. B. (2016).  Understanding relative and
absolute change in discontinuous growth models: Coding alternatives and implications
for hypothesis testing.  Organizational Research Methods, 19, 562-592.

Bliese, P. D., Kautz, J., & Lang, J. W. (2020). Discontinuous growth models: Illustrations,
recommendations, and an R function for generating the design matrix. In Y. Griep & 
S. D. Hansen (Eds.), Handbook on the Temporal Dynamics of Organizational Behavior
(pp. 319-350). Northampton, MA: Edward Elgar Publishers, Inc.}

\author{Paul Bliese
 \email{pdbliese@gmail.com}}


\examples{
##########
# Example 1: Coding trans, post and time.a in data where
# every event occurs at the same time for each person
##########

# Read data from library
data(tankdat)

# Add a marker (1 or 0) indicating an event
tankdat$taskchange<-0
tankdat$taskchange[tankdat$TIME==6]<-1

# Run function with defaults
OUT<-with(tankdat,dgm.code(ID,TIME,taskchange))
names(OUT)
names(tankdat)

# Merge original data and dgm codes and reorder
tankdat.dgm<-merge(tankdat,OUT,by.x=c("ID","TIME"),by.y=c("grp","time"))
tankdat.dgm<-tankdat.dgm[order(tankdat.dgm$ID,tankdat.dgm$TIME),]

# Examine data
tankdat.dgm[1:12,]

##########
# Example 2: Coding trans, post and time.a in data where every transition
# event occurs at the different times for each person
##########
# Read data from library
data(tankdat)

# Add a marker (1 or 0) indicating an event at random
set.seed(343227)
tankdat$taskchange<-rbinom(nrow(tankdat),1,prob=.1)
tankdat[1:24,] #ID 1 had one event at TIME 10. ID 2 had 3 events

# Run function with defaults
\dontrun{
OUT<-with(tankdat,dgm.code(ID,TIME,taskchange))} 
# returns an error showing the 24 groups that started with an event.
# Either drop these groups or change the first.obs option to TRUE
# which changes these first events to 0 (non-events)

OUT<-with(tankdat,dgm.code(ID,TIME,taskchange,first.obs=TRUE))
OUT[1:24,]
OUT[OUT$grp==9,]
#Notice the event.first value of 1 for group 9 indicating that the
#first value was present and recoded.


# In the default setting, one ID had 4 events. It may be preferable 
# to restrict the number of events to 3 or more and code accordingly

OUT<-with(tankdat,dgm.code(ID,TIME,taskchange,n.events=3,first.obs=TRUE))
OUT[1:24,]
}
\keyword{growth models}