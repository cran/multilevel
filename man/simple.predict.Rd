\name{simple.predict}
\alias{simple.predict}
\title{A simple variant of the generic predict command useful for S4 lmer models}
\description{
Statistical models can be used to predict levels of an outcome variable given specific
values of predictors. R has a number of predict functions linked to specific models
(e.g., predict.lm, predict.glm).  At their core, predict commands multiple
the estimated coefficients and the corresponding predictors, and then sum across the intercept
and predictors.  Given the existing functions, there is little value in using simple.predict
in most cases.  As of 2013, however, models produced with lme4 functions (lmer and
glmer) had no specific predict commands, so the simple.predict function was primarily designed
to provide predicted values from lme4 models.
}
\usage{
simple.predict(orig.mod,formula,newdata,dichot=FALSE)
}
\arguments{
 \item{orig.mod}{The original model providing the coefficient estimates.}
 \item{formula}{A formula that does NOT specify the outcome variable. The formula should mirror
      the formula in the orig.mod minus the variable to the left of the ~ sign (e.g., ~IV1+IV2).}
 \item{newdata}{A dataframe with names corresponding to each predictor in orig.mod and providing
      values upon which to estimate the outcome.}
 \item{dichot}{An option used with binomial models.  If TRUE, returns estimates scaled to a
      response vector between 0 and 1. Equivalent to using type="response" in predict.glm for a
      binomial model.}
}
\value{
Returns a vector of predicted values.
}     
\author{
 Paul Bliese
 \email{paul.bliese@us.army.mil}
}
\seealso{
 \code{\link{predict}}
 \code{\link{predict.lm}}
 \code{\link{predict.glm}}
}
\examples{
data(bh1996)  

# A simple example with lm

tmod<-lm(WBEING~HRS+LEAD,bh1996)
simple.predict(tmod,~HRS+LEAD,data.frame(HRS=c(5,10),LEAD=c(3.0,4.5)))
predict(tmod,data.frame(HRS=c(5,10),LEAD=c(3.0,4.5)))  #Using the predict.lm function

# A simple example with glm

tmod.glm<-glm(ifelse(WBEING>3,1,0)~HRS+LEAD,family="binomial",bh1996)
simple.predict(tmod.glm,~HRS+LEAD,data.frame(HRS=c(5,10),LEAD=c(3.0,4.5)),dichot=TRUE)
predict(tmod.glm,data.frame(HRS=c(5,10),LEAD=c(3.0,4.5)),type="response") #Using the predict.glm function

################################################
# An example including a factor as a predictor #
################################################

#Creating a factor
bh1996$LEAD.F<-ifelse(bh1996$LEAD<2.5,0,ifelse(bh1996$LEAD>3.5,2,1))
bh1996$LEAD.F<-factor(bh1996$LEAD.F,levels=c(0,1,2),labels=c("Lo","Med","Hi"))

#Estimating the original model
tmod<-lm(WBEING~HRS+LEAD.F,bh1996)

#Setting up the predictor dataframe
TDAT<-data.frame(HRS=c(5,10),LEAD.F=c("Lo","Med"))
TDAT$LEAD.F<-factor(TDAT$LEAD.F,levels=c("Lo","Med","Hi"))  #This is key in making this work with factors

#Contrasting the predicted values from approaches
simple.predict(tmod,~HRS+LEAD.F,TDAT)
predict(tmod,data.frame(HRS=c(5,10),LEAD.F=c("Lo","Med")))  #Using the predict.lm function

#################################################
#  An Example using the lme4 package and lmer   #
# Delete the preceding # to execute the example #
#################################################

#library(lme4)
#tmod<-lmer(WBEING~HRS+LEAD+(1|GRP),bh1996)
#simple.predict(tmod,~HRS+LEAD,data.frame(HRS=c(5,10),LEAD=c(3.0,4.5)))

}
\keyword{predict}