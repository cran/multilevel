#ad.m ####
ad.m<-function (x, grpid,type="mean") 
{
NEWDAT <- data.frame(x, grpid = grpid)
NEWDAT <- na.exclude(NEWDAT)
DATSPLIT <- split(NEWDAT[, 1:(ncol(NEWDAT) - 1)], NEWDAT$grpid)
# Code to estimate AD Mean on scales
    if(ncol(as.matrix(x))>1){
         ans1 <- lapply(DATSPLIT, function(Q) {
         if (nrow(Q) > 1) {
             mean(apply(Q,2,function(AD){
               sum(abs(AD-eval(call(paste(type),AD))))/length(AD)}))
         }
         else {
            NA
         }
    })
    ans2 <- lapply(DATSPLIT, nrow)
    ans1 <- unlist(ans1)
    ans2 <- unlist(ans2)
    OUTPUT <- data.frame(grpid = names(DATSPLIT), AD.M = ans1, 
        gsize = ans2)
    return(OUTPUT)  
    stop()
  }
#Code to estimate AD Mean on single items
ans1<-lapply(DATSPLIT,function(AD){
      sum(abs(AD-eval(call(paste(type),AD))))/length(AD)})
ans2 <- lapply(DATSPLIT, length)
ans1 <- unlist(ans1)
ans2 <- unlist(ans2)
ans1[ans2==1]<-NA
OUTPUT <- data.frame(grpid = names(DATSPLIT), AD.M = ans1, 
        gsize = ans2)
    return(OUTPUT)
}

#ad.m.sim####

ad.m.sim<-function (gsize, nitems = 1, nresp, itemcors = NULL, type = "mean", 
    nrep) 
{
    OUT <- rep(NA, nrep)
    if (nitems == 1 & is.null(itemcors)) {
        for (i in 1:nrep) {
            OUT[i] <- ad.m(x = sample(1:nresp, gsize, replace = T), 
                grpid = rep(1, gsize), type)[, 2]
        }
    }
    if (nitems > 1 & is.null(itemcors)) {
        for (i in 1:nrep) {
            OUT[i] <- ad.m(x = matrix(sample(1:nresp, gsize * 
                nitems, replace = T), ncol = nitems), grpid = rep(1, 
                gsize), type)[, 2]
        }
    }
    if (!is.null(itemcors)) {
        nitems <- ncol(itemcors)
        for (i in 1:nrep) {
            SIMDAT <- mvrnorm(n = gsize, mu = rep(0, nitems), 
                itemcors)
            SIMDAT <- apply(SIMDAT, 2, cut, breaks = qnorm(c(0, 
                (1/nresp) * 1:nresp)), include.lowest = T, labels = F)
            OUT[i] <- ad.m(x = SIMDAT, grpid = rep(1, gsize), 
                type)[, 2]
        }
   }
 cumpct <- cumsum(table(OUT)/length(OUT))
 lag1 <- c(NA, cumpct[1:length(cumpct) - 1])
 TDAT <- matrix(c(as.numeric(names(cumpct)),cumpct, lag1,1:length(cumpct)),ncol=4)
 TR <- TDAT[TDAT[,2] > 0.05 & TDAT[,3] <= 0.05,4]
 ad.m.05 <- TDAT[TR - 1, 1]
 estout <- list(ad.m = OUT, gsize = gsize, nresp = nresp, 
           nitems = nitems, ad.m.05 = ad.m.05, pract.sig = nresp/6)
 class(estout) <- "disagree.sim"
 return(estout)
}

#quantile.disagree.sim####

quantile.disagree.sim<-function(x, confint, ...)
{
  out<-data.frame(quantile.values=confint,confint.estimate=rep(NA,length(confint)))
    cumpct<-cumsum(table(x[[1]])/length(x[[1]]))
    lag1<-c(NA,cumpct[1:length(cumpct)-1])
    TDAT<-data.frame(dagree.val=as.numeric(names(cumpct)),cumpct,lag1)
    row.names(TDAT)<-1:nrow(TDAT)
  for(i in 1:length(confint)){
     TR<-as.numeric(row.names(TDAT[TDAT$cumpct>confint[i]&TDAT$lag1<=confint[i],]))
     out[i,2]<-TDAT[TR-1,1]
  }
 return(out)
}


#summary.disagree.sim####

summary.disagree.sim<-function(object, ...)
{
    out<-list(summary(object[[1]]),
              object[[2]],
              object[[3]],
              object[[4]],
                  object[[5]],
              object[[6]])
    names(out)<-names(object)
    return(out)
}

#awg ####
awg<-function (x, grpid, range = c(1, 5)) 
{
    NEWDAT <- data.frame(x, grpid = grpid)
    NEWDAT <- na.exclude(NEWDAT)
    DATSPLIT <- split(NEWDAT[, 1:(ncol(NEWDAT) - 1)], NEWDAT$grpid)
    if (ncol(as.matrix(x)) > 1) {
        ans1 <- lapply(DATSPLIT, function(Q) {
            if (nrow(Q) > 1) {
                mean(apply(Q, 2, function(AW) {
                  H <- range[2]
                  L <- range[1]
                  M <- mean(AW)
                  k <- length(AW)
                  A.WG <- 1 - ((2 * var(AW))/(((H + L) * M - 
                    (M^2) - (H * L)) * (k/(k - 1))))
                  if (M < ((L * (k - 1) + H)/k)) 
                    A.WG <- NA
                  if (M > ((H * (k - 1) + L)/k)) 
                    A.WG <- NA
                  if (M == H | M == L) 
                    A.WG = 1
                  A.WG
                }), na.rm = T)
            }
            else {
                NA
            }
        })
        ans2 <- lapply(DATSPLIT, nrow)
        ans3 <- lapply(lapply(DATSPLIT, var),mean,na.rm=T)
        ans1 <- unlist(ans1)
        ans2 <- unlist(ans2)
        ans3 <- unlist(ans3)
        OUTPUT <- data.frame(grpid = names(DATSPLIT), a.wg = ans1, 
            nitems = ncol(as.matrix(x)), nraters = ans2, avg.grp.var = ans3)
        return(OUTPUT)
        stop()
    }
    ans1 <- lapply(DATSPLIT, function(AW) {
        H <- range[2]
        L <- range[1]
        M <- mean(AW)
        k <- length(AW)
        A.WG <- 1 - ((2 * var(AW))/(((H + L) * M - (M^2) - (H * 
            L)) * (k/(k - 1))))
        if (M < ((L * (k - 1) + H)/k)) 
            A.WG <- NA
        if (M > ((H * (k - 1) + L)/k)) 
            A.WG <- NA
        if (M == H | M == L) 
            A.WG = 1
        A.WG
    })
    ans2 <- lapply(DATSPLIT, length)
    ans3 <- lapply(DATSPLIT, var)
    ans1 <- unlist(ans1)
    ans2 <- unlist(ans2)
    ans3 <- unlist(ans3)
    ans1[ans2 == 1] <- NA
    OUTPUT <- data.frame(grpid = names(DATSPLIT), a.wg = ans1, 
        nraters = ans2, grp.var = ans3)
    return(OUTPUT)
}

#boot.icc####
boot.icc<-function(x, grpid, nboot, aov.est=FALSE){
   if(aov.est){
   data<-data.frame(grpid,x)  #fixes data because code below uses the first column to select level 2
   B.OUT<-rep(NA,nboot)
   ngrp<-length(unique(grpid))
   ugrp<-unique(grpid)
   for(i in 1:nboot) {
       #The code below creates an empty list, selects a sample of level 2
       #units and then goes in and samples level-1 units for each level 2 unit
       ROUT<-list(NA)
       rgrps<-sample(ugrp,ngrp,replace=T)
          for (k in 1:ngrp){
             ROUT[[k]]<-data.frame(newgrp=k,data[is.element(data[,1],rgrps[k]),])
             dindex<-sample(nrow(ROUT[[k]]),nrow(ROUT[[k]]),replace=T)
          ROUT[[k]]<-ROUT[[k]][dindex,]
          }
       ROUT<-(do.call(rbind,ROUT))
       tmod<-aov(x~as.factor(newgrp),data=ROUT)  
       B.OUT[i]<-ICC1(tmod)
       }
    return(B.OUT)
   }
   if(!aov.est){
   data<-data.frame(grpid,x)  #fixes data because code below uses the first column to select level 2
   B.OUT<-rep(NA,nboot)
   ngrp<-length(unique(grpid))
   ugrp<-unique(grpid)
   for(i in 1:nboot) {
       #The code below creates an empty list, selects a sample of level 2
       #units and then goes in and samples level-1 units for each level 2 unit
       ROUT<-list(NA)
       rgrps<-sample(ugrp,ngrp,replace=T)
          for (k in 1:ngrp){
             ROUT[[k]]<-data.frame(newgrp=k,data[is.element(data[,1],rgrps[k]),])
             dindex<-sample(nrow(ROUT[[k]]),nrow(ROUT[[k]]),replace=T)
          ROUT[[k]]<-ROUT[[k]][dindex,]
          }
       ROUT<-(do.call(rbind,ROUT))
       tmod<-lme(x~1, random=~1|newgrp,data=ROUT,control=list(opt="optim"))  
       temp<-VarCorr(tmod)
         Tau<-as.numeric(temp[[1]])
       Sigma.Sq<-(tmod$sigma)^2
       B.OUT[i]<-Tau/(Tau+Sigma.Sq)
       }
   return(B.OUT)
   }
}
#cordif####
cordif<-function(rvalue1,rvalue2,n1,n2){
	zvalue1<-.5*((log(1+rvalue1))-(log(1-rvalue1)))
	zvalue2<-.5*((log(1+rvalue2))-(log(1-rvalue2)))
	zest<-(zvalue1-zvalue2)/sqrt(1/(n1-3)+1/(n2-3))
	out<-list("z value"=zest)
	return(out)
}

#cordif.dep####
cordif.dep<-function(r.x1y, r.x2y, r.x1x2, n)
{
#
# This function tests whether two dependent correlations are significantly
# different from each other.  The formula is taken from Cohen & Cohen (1983)
# p. 56
#
        rbar <- (r.x1y + r.x2y)/2
        barRbar <- 1 - r.x1y^2 - r.x2y^2 - r.x1x2^2 + 2 * r.x1y * r.x2y * r.x1x2
        tvalue.num <- ((r.x1y - r.x2y) * sqrt((n - 1) * (1 + r.x1x2)))
        tvalue.den <- sqrt(((2 * ((n - 1)/(n - 3))) * 
                      barRbar + ((rbar^2)) * (1 - r.x1x2)^3))
        t.value <- tvalue.num/tvalue.den
        DF <- n - 3
        p.value <- (1 - pt(abs(t.value), DF)) * 2
        OUT <- data.frame(t.value, DF, p.value)
        return(OUT)
}
#cronbach####
cronbach<-function(items)
{
	items<-na.exclude(items)	
        N <- ncol(items)
        TOTVAR <- var(apply(items, 1, sum))
        SUMVAR <- sum(apply(items, 2, var))
        ALPHA <- (N/(N - 1)) * (1 - (SUMVAR/TOTVAR))
	OUT<-list(Alpha=ALPHA,N=nrow(items))        
	return(OUT)
}

#dgm.code####
dgm.code<-function(grp,time,event,n.events=FALSE,first.obs=FALSE){
  #ensure data structure is correct
  newdata<-data.frame(grp=grp,time=time,event=event)
  newdata<-na.exclude(newdata)
  newdata<-newdata[order(newdata$grp,newdata$time),]
  
  #If first observation is an event and first.obs is TRUE
  #change the first observation to a non-event
  if(first.obs){
    fobs.grps<-newdata$grp[!duplicated(newdata$grp)&(newdata$event==1)]
    #print(fobs.grps)
    newdata$event[!duplicated(newdata$grp)&(newdata$event==1)]<-0
  }
  
  #Check to see if first observation is an event
  s.event<-nrow(newdata[!duplicated(newdata$grp)&(newdata$event==1),])
  if(s.event>0){
    print("The following groups start with an event")
    print(newdata[!duplicated(newdata$grp)&(newdata$event==1),])
    print("Drop the groups or use the first.obs=TRUE option")
    stop() 
  }
  
  #count the maximum number of events for any group
  max.events<-max(with(newdata,tapply(event,grp,sum)))
  n.grps<-length(unique(newdata$grp))
  
  #adjust the maximum number of events to a specified level
  if(n.events){
    max.events=n.events
  }
  
  # Set up the structure for the output
  ANS<-matrix(0,nrow(newdata),ncol=max.events*2)
  ANS<-data.frame(ANS)
  names(ANS)<-c(paste0("trans",c(1:max.events)),paste0("post",c(1:max.events)))
  ANS<-data.frame(ANS,grp=newdata$grp,time=newdata$time,event=newdata$event)
  g.size<-tapply(newdata$grp,newdata$grp,length)
  #  ANS$cum.event<-unlist(tapply(ANS$event,ANS$grp,cumsum))
  ANS$time.a<-ANS$time
  ANS$cum.event<-do.call(c, tapply(ANS$event, ANS$grp, FUN=cumsum))
  ANS$num.grp<-rep(1:n.grps,times=g.size) #create a numeric group for loops
  
  # Add two check variables for total events and whether an event
  # occured on the first occasion
  ANS$tot.events<-NA
  ANS$event.first<-0
  
  if(first.obs){
    ANS$event.first[ANS$grp%in%fobs.grps]<-1
  }
  
  #collapse the number of events to a specified level
  if(n.events){
    ANS$cum.event<-ifelse(ANS$cum.event>n.events,n.events,ANS$cum.event)
  }
  
  # Set up a factor outside of the loop to get all levels
  ANS$cum.event.f<-factor(ANS$cum.event,levels=c(0:max.events))
  
  # Set up a loop to put values in trans and post variables
  # First skip groups with no events
  for(i in 1:n.grps){
    if(sum(ANS$event[ANS$num.grp==i])==0){
      ANS$tot.events[ANS$num.grp==i]<-0
      next(i)
    }
    # Use model.matrix to set up dummy codes for trans and post
    ANS[ANS$num.grp==i,1:max.events]<-model.matrix(~C(cum.event.f,contr.treatment),data=ANS[ANS$num.grp==i,])[,2:(max.events+1)]
    ANS[ANS$num.grp==i,(max.events+1):(2*max.events)]<-model.matrix(~C(cum.event.f,contr.treatment),data=ANS[ANS$num.grp==i,])[,2:(max.events+1)]
    if(max.events==1){
      ANS[ANS$num.grp==i,2]<-cumsum(ANS[ANS$num.grp==i,2])-1
    }
    if(max.events>1){
      ANS[ANS$num.grp==i,(max.events+1):(2*max.events)]<-apply(ANS[ANS$num.grp==i,(max.events+1):(2*max.events)],2,cumsum)-1
    }
    ANS$time.a[ANS$num.grp==i]<-ifelse(ANS$cum.event[ANS$num.grp==i]==0,ANS$time[ANS$num.grp==i],NA)
    ANS$time.a[is.na(ANS$time.a)&ANS$num.grp==i]<-max(ANS$time.a[!is.na(ANS$time.a)&ANS$num.grp==i])
    ANS$tot.events[ANS$num.grp==i]<-sum(ANS$event[ANS$num.grp==i])
    next(i)
  }
  # Clean up the ANS matrix column by column
  # print(ANS) to see the structure of the previous loop
  for(j in 1:max.events){
    ANS[,max.events+j]<-ifelse(ANS[,j]==0,0,ANS[,max.events+j])
    next(j)
  }
  
  #rearrange the ANS matrix for output
  #print(ANS[1:30,]) to see the first 30 rows of complete data
  ANS[,c((max.events*2)+1:3,1:(max.events*2),(max.events*2)+c(4,7,8))]
}

#GmeanRel####
gmeanrel<-function(object)
{
	OUTFILE<-aggregate(object$group,object$group,length)
	names(OUTFILE)<-c("Group","GrpSize")
	temp<-VarCorr(object)
	Tau<-as.numeric(temp[[1]])
	Sigma.Sq<-(object$sigma)^2
	ICC<-Tau/(Tau+Sigma.Sq)
	OUTFILE$GmeanRel<-(OUTFILE[,2]*ICC)/(1+(OUTFILE[,2]-1)*ICC)
	estout<-list(ICC=ICC,Group=OUTFILE[,1],GrpSize=OUTFILE[,2],MeanRel=OUTFILE[,3])
	class(estout)<-"gmeanrel"
	return(estout)
}
#graph.ran.mean####
graph.ran.mean<-function(x, grpid, nreps, limits, graph=TRUE, bootci=FALSE)
{
 if(bootci){
    if (missing(limits)) 
        limits <- quantile(x[is.na(x) == F], c(0.10, 0.90))
    if (is.factor(grpid)) 
        grpid <- grpid[, drop = TRUE]
    TDAT<-na.exclude(data.frame(x,grpid))
    x<-TDAT[,1]
    grpid<-TDAT[,2]
    NGRPS <- length(tapply(x, grpid, length))
    OUT <- matrix(NA, NGRPS, nreps)
    for (i in 1:nreps) {
        TOUT <- mix.data(x, grpid)
        OUT[, i] <- sort(tapply(TOUT[, 3], TOUT[, 1], mean, na.rm = T))
    }
    REALGRP <- sort(tapply(x, grpid, mean, na.rm = T))
    if (graph) {
        plot(c(REALGRP, max(REALGRP)), type = "h", ylim = limits, 
            ylab = "Group Average")
        lines(c(REALGRP, max(REALGRP)), type = "s")
        PSEUDOMEAN <- apply(OUT, 1, mean)
        lines(PSEUDOMEAN, type = "l")
        PSEUDO.LCI <- apply(OUT, 1, quantile, 0.025)
        lines(PSEUDO.LCI, type = "l",lty=2)
        PSEUDO.HCI <- apply(OUT, 1, quantile, 0.975)
        lines(PSEUDO.HCI, type = "l",lty=2)
    }
    else {
        REALGRP <- sort(tapply(x, grpid, mean, na.rm = T))
        GRPNAMES <- row.names(REALGRP)
        REALGRP <- as.vector(REALGRP)
        PSEUDOMEAN <- apply(OUT, 1, mean)
        PSEUDO.LCI <- apply(OUT, 1, quantile, .025)
        PSEUDO.HCI <- apply(OUT, 1, quantile, .975)
        OUT <- data.frame(GRPNAMES, GRPMEAN = REALGRP, 
            PSEUDOMEAN, PSEUDO.LCI, PSEUDO.HCI)
        return(OUT)
    }
  }
 if(!bootci){
    if (missing(limits)) 
        limits <- quantile(x[is.na(x) == F], c(0.10, 0.90))
    if (is.factor(grpid)) 
        grpid <- grpid[, drop = TRUE]
    TDAT<-na.exclude(data.frame(x,grpid))
    x<-TDAT[,1]
    grpid<-TDAT[,2]
    NGRPS <- length(tapply(x, grpid, length))
    OUT <- matrix(NA, NGRPS, nreps)
    for (i in 1:nreps) {
        TOUT <- mix.data(x, grpid)
        OUT[, i] <- sort(tapply(TOUT[, 3], TOUT[, 1], mean, na.rm = T))
    }
    REALGRP <- sort(tapply(x, grpid, mean, na.rm = T))
    if (graph) {
        plot(c(REALGRP, max(REALGRP)), type = "h", ylim = limits, 
            ylab = "Group Average")
        lines(c(REALGRP, max(REALGRP)), type = "s")
        PSEUDOGRP <- apply(OUT, 1, mean)
        lines(PSEUDOGRP, type = "l")
    }
    else {
        REALGRP <- sort(tapply(x, grpid, mean, na.rm = T))
        GRPNAMES <- row.names(REALGRP)
        REALGRP <- as.vector(REALGRP)
        PSEUDOGRP <- apply(OUT, 1, mean)
        OUT <- data.frame(GRPNAMES = GRPNAMES, GRPMEAN = REALGRP, 
            PSEUDOMEAN = PSEUDOGRP)
        return(OUT)
    }
  }
}

#ICC1####
ICC1<- function(object)
{
	MOD <- summary(object)
	MSB <- MOD[[1]][1, 3]
	MSW <- MOD[[1]][2, 3]
	GSIZE <- (MOD[[1]][2, 1] + (MOD[[1]][1, 1] + 1))/(MOD[[1]][1, 1] + 1)	
	#	print(GSIZE)
	OUT <- (MSB - MSW)/(MSB + ((GSIZE - 1) * MSW))
	return(OUT)
}

#ICC2####
ICC2 <-function(object)
{
	MOD <- summary(object)
	MSB <- MOD[[1]][1, 3]
	MSW <- MOD[[1]][2, 3]
	OUT <- (MSB - MSW)/MSB
	return(OUT)
}
#item.total####
item.total<-function(items)
{
items<-na.exclude(items)
        N <- ncol(items)
        ans <- matrix(0, N, 3)
        ans[, 1] <- labels(items)[[2]]
        for(i in 1:N) {
                ans[i, 2] <- cor(items[, i], apply(items[,  - i], 1, mean))
                ans[i, 3] <- cronbach(items[,  - i])[[1]]
        }
        OUT <- data.frame(Variable=ans[,1],Item.Total=as.numeric(ans[,2]), 
                Alpha.Without=as.numeric(ans[,3]),N=nrow(items))
        return(OUT)
}

#make.univ####
make.univ<-function (x, dvs, tname="TIME", outname="MULTDV") 
{
    NREPOBS <- ncol(dvs)        
    UNIV.DAT <- x[rep(1:nrow(x), rep(NREPOBS, nrow(x))), 1:ncol(x)]
    FINAL.UNIV <- data.frame(timedat = rep(0:(NREPOBS - 1), 
        nrow(x)), outdat = as.vector(t(dvs)))
    names(FINAL.UNIV)<-c(tname,outname)
    FINAL.DAT <- data.frame(UNIV.DAT, FINAL.UNIV)
    return(FINAL.DAT)
}

#mix.data####
mix.data<-function (x, grpid) 
{
    TDAT <- cbind(rnorm(length(grpid)), grpid, x)
    TDAT <- TDAT[is.na(grpid) == F & grpid != "NA", ]
    TDAT <- TDAT[order(TDAT[, 1]),1:ncol(TDAT)]
    TMAT <- tapply(TDAT[, 2], TDAT[, 2], length)
    NGRPS <- length(TMAT)
    newid <- rep(1:NGRPS, TMAT)
    OUT <- cbind(newid, TDAT[, 2:ncol(TDAT)])
    return(OUT)
}
#mult.icc####
mult.icc<-function (x, grpid) 
{
    ans <- data.frame(Variable = names(x[, 1:ncol(x)]), ICC1 = as.numeric(rep(NA, 
        ncol(x))), ICC2 = as.numeric(rep(NA, ncol(x))))
    GSIZE <- mean(aggregate(grpid, list(grpid), length)[,2])
    for (i in 1:ncol(x)) {
        DV <- x[, i]
        tmod <- lme(DV ~ 1, random = ~1 | grpid, na.action = na.omit,
              control=list(opt="optim"))
        TAU <- as.numeric(VarCorr(tmod)[, 1][1])
        SIGMASQ <- tmod$sigma^2
        ICC1 <- TAU/(TAU + SIGMASQ)
        ICC2 <- (GSIZE * ICC1)/(1 + (GSIZE - 1) * ICC1)
        ans[i, 2] <- ICC1
        ans[i, 3] <- ICC2
    }
    return(ans)
}
#mult.make.univ####
mult.make.univ <- function(x,dvlist,tname="TIME",outname="MULTDV")
{
  NREPOBS <- length(dvlist[[1]])
  UNIV.DAT <- x[rep(1:nrow(x), rep(NREPOBS, nrow(x))), 1:ncol(x)]
  FINAL.UNIV <- data.frame(timedat = rep(0:(NREPOBS - 1), nrow(x)), as.data.frame(lapply(dvlist,function(cols) {as.vector(t(x[,cols]))})))
  if (is.null(names(dvlist)))
  {
    names(FINAL.UNIV) <- c(tname,paste(outname,1:(ncol(FINAL.UNIV)-1),sep=''))
  }else{
    names(FINAL.UNIV) <- c(tname,names(dvlist))
  }
  FINAL.DAT <- data.frame(UNIV.DAT,FINAL.UNIV)
  return(FINAL.DAT)
}
#sam.cor####
sam.cor<-function(x,rho)
{
	y <- (rho * (x - mean(x)))/sqrt(var(x)) + sqrt(1 - rho^2) * rnorm(length(x))
	cat("Sample corr = ", cor(x, y), "\n")
	return(y)
}

#rmv.blanks####
rmv.blanks<-function (object) 
{
    OUT <- lapply(object, function(xsub) {
        ANY.BLNK <- grep(" +$", xsub)
        if (length(ANY.BLNK) < length(xsub)) 
            xsub <- xsub
        else xsub <- sub(" +$", "", xsub)
    })
    return(data.frame(OUT))
}
#rgr.agree####
rgr.agree<-function (x, grpid, nrangrps) 
{
    GVARDAT <- tapply(x, grpid, var)
    NGRPS <- length(GVARDAT)
    GSIZE <- tapply(grpid, grpid, length)
    if(min(GSIZE)<2){
       print("One or more groups has only one group member.")
       stop("There must be at least two group members per group to estimate rgr.agree.")
        }
    NREPS <- round((nrangrps/NGRPS), digits = 0)
    ans <- rep(0, (length(GSIZE) * NREPS))
    for (i in 1:NREPS) {
        ans[((i * length(GSIZE)) - (length(GSIZE)) + 1):(i * 
            length(GSIZE))] <- ran.group(x, grpid, var)
    }
    AVGRPVAR <- mean(GVARDAT)
    NGRPS <- length(GVARDAT)
    RGRVAR <- mean(ans)
    RGRSD <- sqrt(var(ans))
    ZVALUE <- (AVGRPVAR - RGRVAR)/(RGRSD/sqrt(NGRPS))
     estout <- list(NRanGrp =length(ans),
                        AvRGRVar = RGRVAR,
                        SDRGRVar = RGRSD,
                        AvGrpVar = AVGRPVAR,
                        zvalue = ZVALUE,
                        RGRVARS =ans)
    class(estout)<-"rgr.agree"
    return(estout)
}

ran.group<-function(x, grpid, fun, ...)
{
        if(!is.null(ncol(x))) {
                GSIZE <- tapply(grpid, grpid, length)
                ans <- rep(0, length(GSIZE))
                if(length(x[, 1]) != sum(GSIZE))
                        stop("The sum of group sizes does not match the number of observations"
                                )
                for(i in 1:length(GSIZE)) {
                        GID2 <- c(1:length(x[, 1]))
                        SAM <- sample(GID2, size = GSIZE[i])
                        ans[i] <- mean(apply(x[SAM,  ], 2, fun))
                        x <- x[ - SAM,  ]
                }
                return(ans)
        }
        GSIZE <- tapply(grpid, grpid, length)
        ans <- rep(0, length(GSIZE))
        if(length(x) != sum(GSIZE,na.rm=T))
                stop("The sum of group sizes does not match the number of observations"
                        )
        for(i in 1:length(GSIZE)) {
                GID2 <- c(1:length(x))
                SAM <- sample(GID2, size = GSIZE[i])
                ans[i] <- fun(x[c(SAM)])
                x <- x[ - SAM]
        }
        ans
}

summary.rgr.agree<-function(object, ...)
{
    Table <- data.frame(object$NRanGrp, object$AvRGRVar, object$SDRGRVar, 
        object$AvGrpVar, object$zvalue)
    names(Table) <- c("N.RanGrps", "Av.RanGrp.Var", "SD.Rangrp.Var", 
        "Av.RealGrp.Var", "Z-value")
    object$Table <- as.matrix(Table)
    object$lowercis<-quantile(object$RGRVARS,c(.005,.01,.025,.05,.10))
    object$uppercis<-quantile(object$RGRVARS,c(.90,.95,.975,.99,.995))
    OUT<-list(object$Table,object$lowercis,object$uppercis)
    names(OUT)<-c("Summary Statistics for Random and Real Groups","Lower Confidence Intervals (one-tailed)",
    "Upper Confidence Intervals (one-Tailed)")
    OUT
}


#rgr.OLS####
rgr.ols<-function(xdat1, xdat2, ydata, grpid, nreps)
{
#
# The number of columns in the output matrix has to correspond to
# the number of mean squares you want in the output.
# This function does RGR on a two IV OLS hierarchical OLS model.
#
OUT <- matrix(0, nreps, 4)
NEWDAT <- cbind(grpid, xdat1, xdat2, ydata)
for(k in 1:nreps) {
	TDAT <- mix.data(NEWDAT, grpid)
	Y <- tapply(TDAT[, 6], TDAT[, 1], mean)
	X1 <- tapply(TDAT[, 4], TDAT[, 1], mean)
	X2 <- tapply(TDAT[, 5], TDAT[, 1], mean)
	MOD <- lm(Y ~ X1 * X2)
	#print(anova(MOD,test="F"))
	TOUT <- anova(MOD, test = "F")[, 3]
	OUT[k,  ] <- TOUT
}
return(OUT)
}
#rgr.waba####
rgr.waba<-function(x, y, grpid, nrep)
{
#
# Create Matrix and sort it by Group ID
#
        SMAT <- cbind(grpid, x, y)
        SMAT <- SMAT[order(SMAT[, 1]), 1:3]
        GID.S <- SMAT[, 1]
        X.S <- SMAT[, 2]
        Y.S <- SMAT[, 3]        #
#
# Create a matrix in which to put the random WABA elements
#
        ans <- matrix(1, nrep, 7)       #
#
# WABA random group loop
#
        for(i in 1:nrep) {
#
# Generate a random number and sort x and y by it
#
                TR <- rnorm(length(X.S))
                T.DAT <- cbind(TR, X.S, Y.S)
                T.DAT <- T.DAT[order(T.DAT[, 1]), 1:3]  #
#
# Create a Matrix in which to put the WABA elements
#
                tmat <- matrix(0, length(X.S), 6)       #
#
# Split up the x observations by the Group ID and make WABA elements
#
                TX <- split(T.DAT[, 2], GID.S)
                TX.M <- unlist(lapply(TX, mean))
                TX.L <- unlist(lapply(TX, length))
                tmat[, 1] <- T.DAT[, 2]
                tmat[, 2] <- rep(TX.M, TX.L)
                tmat[, 3] <- (T.DAT[, 2] - tmat[, 2])   #
#
# Split up the y observations by the Group ID and make WABA elements
#
                TY <- split(T.DAT[, 3], GID.S)
                TY.M <- unlist(lapply(TY, mean))
                TY.L <- unlist(lapply(TY, length))
                tmat[, 4] <- T.DAT[, 3]
                tmat[, 5] <- rep(TY.M, TY.L)
                tmat[, 6] <- T.DAT[, 3] - tmat[, 5]     #
#
# Calculate WABA parameters and put them in a Matrix format
#
                ans[i, 1] <- cor(tmat[, 1], tmat[, 4])
                ans[i, 2] <- cor(tmat[, 1], tmat[, 2])
                ans[i, 3] <- cor(tmat[, 4], tmat[, 5])
                ans[i, 4] <- cor(tmat[, 2], tmat[, 5])
                ans[i, 5] <- cor(tmat[, 1], tmat[, 3])
                ans[i, 6] <- cor(tmat[, 4], tmat[, 6])
                ans[i, 7] <- cor(tmat[, 3], tmat[, 6])
        }
    estout <- data.frame(ans)
    names(estout) <- c("RawCorr", "EtaBx", "EtaBy", "CorrB", "EtaWx", 
        "EtaWy", "CorrW")
    class(estout) <- "rgr.waba"
    return(estout)

}

summary.rgr.waba<-function(object, ...)
{
    T.DAT <- rep(0, 3)
    object2<-data.frame(object$RawCorr,object$EtaBx,object$EtaBy,
             object$CorrB,object$EtaWx,object$EtaWy,object$CorrW)
    ans <- data.frame(RawCorr = T.DAT, EtaBx = T.DAT, EtaBy = T.DAT, 
        CorrB = T.DAT, EtaWx = T.DAT, EtaWy = T.DAT, CorrW = T.DAT, 
        row.names = c("NRep", "Mean", "SD"))
    for (i in 1:7) {
        ans[1, i] <- length(object2[, i])
        ans[2, i] <- mean(object2[, i])
        ans[3, i] <- sqrt(var(object2[, i]))
    }
    return(ans)
}



quantile.rgr.waba<-function (x, confint, ...) 
{
    object2 <- data.frame(x$EtaBx, x$EtaBy, 
        x$CorrB, x$EtaWx, x$EtaWy, x$CorrW)
    names(object2)<-c("EtaBx","EtaBy","CorrB","EtaWx","EtaWy","CorrW")
    ans<-apply(object2,2,quantile,confint)
    return(ans)
}

#rtoz####
rtoz<-function(rvalue){
	zest<-.5*((log(1+rvalue))-(log(1-rvalue)))
	#out<-list("z prime"=zest)
	return(zest)
}

#rwg####
rwg<-function(x, grpid, ranvar=2) 
{

    NEWDAT<-data.frame(x=x,grpid=grpid)
    NEWDAT<-na.exclude(NEWDAT)
    DATSPLIT <- split(NEWDAT$x, NEWDAT$grpid)
    ans1 <- lapply(DATSPLIT, function(Q) {
        if (length(Q) > 1) {
        V <- var(Q)
        if (V > ranvar) 
            V <- ranvar
        out <- 1 - (V/ranvar)
        out}
else {out<-NA
out}
    })
    ans2<-lapply(DATSPLIT,length)
    ans1 <- unlist(ans1)
    ans2<-unlist(ans2)
    OUTPUT <- data.frame(grpid=names(DATSPLIT),rwg = ans1, gsize = ans2)
    return(OUTPUT)
}



#rwg.j####
rwg.j<-function(x, grpid,ranvar=2)
{
    NEWDAT<-data.frame(x,grpid=grpid)
    NEWDAT<-na.exclude(NEWDAT)
    DATSPLIT <- split(NEWDAT[,1:(ncol(NEWDAT)-1)], NEWDAT$grpid)
    ans1 <- lapply(DATSPLIT, function(Q) {
        J <- ncol(Q)
        if (nrow(Q) > 1) {
            S <- mean(apply(Q, 2, var,na.rm=T))
            if (S > ranvar) 
                S <- ranvar
            out <- (J * (1 - (S/ranvar)))/((J * (1 - (S/ranvar))) + 
                (S/ranvar))
            out
        }
        else {out<-NA
out
        }
    })
    ans2<-lapply(DATSPLIT,nrow)
    ans1 <- unlist(ans1)
    ans2 <-unlist(ans2)
    OUTPUT <- data.frame(grpid=names(DATSPLIT),rwg.j = ans1, gsize = ans2)
    return(OUTPUT)
}

#rwg.j.lindell####
rwg.j.lindell<-function (x, grpid, ranvar = 2) 
{
    NEWDAT<-data.frame(x,grpid=grpid)
    NEWDAT<-na.exclude(NEWDAT)
    DATSPLIT <- split(NEWDAT[,1:(ncol(NEWDAT)-1)], NEWDAT$grpid)
    ans1 <- lapply(DATSPLIT, function(Q) {
        if (nrow(Q) > 1) {
            S <- mean(apply(Q, 2, var))
            out <- 1-(S/ranvar)
            out
        }
        else {out<-NA
out
        }
    })
    ans2<-lapply(DATSPLIT,nrow)
    ans1 <- unlist(ans1)
    ans2 <-unlist(ans2)
    OUTPUT <- data.frame(grpid=names(DATSPLIT),rwg.lindell = ans1, gsize = ans2)
    return(OUTPUT)
}



#rwg.sim####
rwg.sim<-function (gsize, nresp, nrep) 
{
    OUT <- rep(NA, nrep)
    for (i in 1:nrep) {
        OUT[i] <- rwg(x = sample(1:nresp, gsize, replace = T), 
            grpid = rep(1, gsize), ranvar = (nresp^2 - 1)/12)[,2]
    }
    cumpct <- cumsum(table(OUT)/length(OUT))
    lag1 <- c(NA, cumpct[1:length(cumpct) - 1])
    lag2 <- c(NA, lag1[1:length(lag1) - 1])
    TDAT <- matrix(c(as.numeric(names(cumpct)),cumpct,lag1,lag2),ncol=4)
    rwg.95 <- TDAT[TDAT[,2] > 0.95 & TDAT[,3] >= 0.95 & TDAT[,4] < 0.95, 1]
    estout <- list(rwg = OUT, gsize = gsize, nresp = nresp, 
        nitems = 1, rwg.95 = rwg.95)
    class(estout) <- "agree.sim"
    return(estout)
}

#rwg.j.sim####

rwg.j.sim<-function(gsize, nitems, nresp, itemcors = NULL, nrep) 
{
    OUT <- rep(NA,nrep)
    if (is.null(itemcors)) {
        for (i in 1:nrep) {
            OUT[i] <- rwg.j(x = matrix(sample(1:nresp, gsize * 
                nitems, replace = T), ncol = nitems), grpid = rep(1, 
                gsize), ranvar = (nresp^2 - 1)/12)[, 2]
        }
    }
    if (!is.null(itemcors)){
        for (i in 1:nrep) {
           nitems <- ncol(itemcors)
           SIMDAT <- mvrnorm(n = gsize, mu = rep(0, nitems), itemcors)
           SIMDAT <- apply(SIMDAT, 2, cut, breaks = qnorm(c(0, (1/nresp) * 
                1:nresp)), include.lowest = T, labels = F)
           OUT[i] <- rwg.j(SIMDAT, grpid = rep(1, gsize), ranvar = ((nresp^2 - 
            1)/12))[, 2]
        }
    }
    cumpct <- cumsum(table(OUT)/length(OUT))
    lag1 <- c(NA, cumpct[1:length(cumpct) - 1])
    lag2 <- c(NA, lag1[1:length(lag1) - 1])
    TDAT <- matrix(c(as.numeric(names(cumpct)),cumpct,lag1,lag2),ncol=4)
    rwg.95 <- TDAT[TDAT[,2] > 0.95 & TDAT[,3] >= 0.95 & TDAT[,4] < 0.95, 1]
    estout <- list(rwg = OUT, gsize = gsize, nresp = nresp, 
        nitems = nitems, rwg.95 = rwg.95)
    class(estout) <- "agree.sim"
    return(estout)
}


#summary.agree.sim####

summary.agree.sim<-function(object, ...)
{
    out<-list(summary(object[[1]]),
              object[[2]],
              object[[3]],
              object[[4]],
              object[[5]])
    names(out)<-names(object)
    return(out)
}

#quantile.agree.sim####

quantile.agree.sim<-function(x, confint, ...)
{
  out<-data.frame(quantile.values=confint,confint.estimate=rep(NA,length(confint)))
    cumpct<-cumsum(table(x[[1]])/length(x[[1]]))
    lag1<-c(NA,cumpct[1:length(cumpct)-1])
    lag2<-c(NA,lag1[1:length(lag1)-1])
    TDAT<-data.frame(agree.val=as.numeric(names(cumpct)),cumpct,lag1,lag2)
  for(i in 1:length(confint)){
    out[i,2]<-TDAT[TDAT$cumpct>confint[i]&TDAT$lag1>=confint[i]&TDAT$lag2<confint[i],1]
  }
 return(out)
}

#simbias####
simbias<-function(corr, gsize, ngrp, icc1x, icc1y, nrep)
{
ANS <- matrix(NA, nrep, 8)
NOBS <- gsize * ngrp
GID <- sort(rep(1:ngrp, gsize))
for(i in 1:nrep) {
X<-rnorm(NOBS)
Y<- scale(sam.cor(X, rho = corr))
GSTD <- sqrt((1/(1 - icc1y)) - 1)
TG.Y <- rnorm(ngrp, 0, GSTD)
TDATY <- sort(rep(TG.Y, gsize))
Y.2 <- X + TDATY
ICC1.Y <- ICC1(aov(Y.2 ~ as.factor(GID)))
GSTDX <- sqrt((1/(1 - icc1x)) - 1)
TG.X <- rnorm(ngrp, 0, GSTDX)
TDATX <- sort(rep(TG.X, gsize))
X.2 <- Y + TDATX
ICC1.X <- ICC1(aov(X.2 ~ as.factor(GID)))
ANS[i, 1] <- ICC1.X
ANS[i, 2] <- ICC1.Y
TEMPFRAME <- data.frame(GID, X.2, Y.2)
TEMPAGG <- aggregate(TEMPFRAME[, 2:3], list(TEMPFRAME[, 1]), mean)
names(TEMPAGG) <- c("GID", "GX.2", "GY.2")
TEMPFRAME <- merge(TEMPFRAME, TEMPAGG, by = "GID", all.x = T)
tmod <- lme(Y.2 ~ X.2 + GX.2, random =  ~ 1 | GID, data = TEMPFRAME)
ANS[i, 3:5] <- summary(tmod)$tTable[2, c(1:2, 4)]
tmod.lm <- lm(Y.2 ~ X.2 + GX.2, data = TEMPFRAME)
ANS[i, 6:8] <- summary(tmod.lm)$coef[2, 1:3]
}
#
#
ANS <- data.frame(ANS)
names(ANS) <- c("icc1.x", "icc1.y", "lme.coef", "lme.se", "lme.tvalue", 
"lm.coef", "lm.se", "lm.tvalue")
return(ANS)
}

#sim.icc####
sim.icc<-function (gsize, ngrp, icc1,nitems=1,item.cor=FALSE) 
{
    if(!item.cor){
    ANS <- matrix(NA, ngrp*gsize, nitems+1)
    GID <- sort(rep(1:ngrp, gsize))
    ANS[,1]<-GID
    NOBS<-gsize*ngrp
    for (i in 1:nitems) {
        X <- rnorm(NOBS)
        GSTD <- sqrt((1/(1 - icc1)) - 1)
        #GSTD <- GSTD/(gsize/(gsize-1))# potential n-1 fix: ngrp or gsize
        TG.X <- rnorm(ngrp, 0, GSTD)
        TDATX <- sort(rep(TG.X, gsize))
        X.2 <- X + TDATX
        ANS[, i+1] <- X.2
    }
    ANS <- data.frame(ANS)
    names(ANS)<-c("GRP",paste0("VAR",c(1:nitems)))
    return(ANS)
    }
    if(item.cor){
    if(nitems==1){
       print("You must generate more than one item")
       stop()
       }
    ANS <- matrix(NA, ngrp*gsize, nitems+1)
    GID <- sort(rep(1:ngrp, gsize))
    ANS[,1]<-GID
    NOBS<-gsize*ngrp
    X <- rnorm(NOBS)
    for (i in 1:nitems) {
        X2<-sqrt(item.cor)*X+sqrt(1-item.cor)*rnorm(NOBS)
        GSTD <- sqrt((1/(1 - icc1)) - 1)
        #GSTD <- GSTD/(gsize/(gsize-1)) #potential n-1 fix: ngrp or gsize
        TG.X <- rnorm(ngrp, 0, GSTD)
        TDATX <- sort(rep(TG.X, gsize))
        VALUE <- X2 + TDATX
        ANS[, i+1] <- VALUE
    }
    ANS <- data.frame(ANS)
    names(ANS)<-c("GRP",paste0("VAR",c(1:nitems)))
    return(ANS)
    }
}
#sim.mlcor ####
sim.mlcor<-function(gsize,ngrp,gcor,wcor,icc1x,icc1y,alphax=1,alphay=1){
  # The group-level correlation is simulated using means from a
  # third variable (GMEANS). Because each variable is correlated with
  # this 3rd variable, we need to "double" the group correlation by taking the 
  # square-root of the correlation. This way the product of the two "paths"
  # give the correct correlation. An if statement handles negative values below.
  gcor.s<-sqrt(abs(gcor))
  #Create two variables: VAR1 for X, VAR2 for Y. The correlation between
  #VAR1 (X) and VAR2 (Y) is the within, but since there are no group-level
  #properties associated with VAR1 (X) it is fine to just run a raw
  #correlation for the within estimate.
  VAR1<-rnorm(gsize*ngrp)
  VAR2<-(wcor * (VAR1 - mean(VAR1)))/sqrt(var(VAR1)) + sqrt(1 - wcor^2) * rnorm(length(VAR1))
  #Adjust the Within for reliability
  ALPHAX<-sqrt(alphax)
  VAR1<-(ALPHAX * (VAR1 - mean(VAR1)))/sqrt(var(VAR1)) + sqrt(1 - ALPHAX^2) *  rnorm(length(VAR1))
  ALPHAY<-sqrt(alphay)
  VAR2<-(ALPHAY * (VAR2 - mean(VAR2)))/sqrt(var(VAR2)) + sqrt(1 - ALPHAY^2) *  rnorm(length(VAR2))
  #The Vector below represents the separate Group Mean vector that G.X and G.Y
  #are built on. Add a GROUP ID
  GMEANS<-rnorm(ngrp)
  GID<-rep(1:ngrp,each=gsize)
  #The Code below calculates a group mean for G.X that is correlated with 
  #the group means (GMEANS). Then it repeats the estimate within each group
  #so the number of values matches the number of level-1 observations, and 
  #it can be combined with the individual variable to create the final X.
  X.G<-(gcor.s * (GMEANS - mean(GMEANS)))/sqrt(var(GMEANS))+sqrt(1-gcor.s^2)*rnorm(ngrp)
  X.G<-rep(X.G,each=gsize)
  #The line below creates the final X variable by combining VAR1 (the Level 1
  #variable) and X.G (the expanded level 2 variable) together with weights
  #determined by the ICC variable for X. The higher the ICC, the more X.G is weighted.
  X<-sqrt(icc1x)*X.G+sqrt(1-icc1x)*VAR1
  # Same set of procedures to create a Y variable
  Y.G<-(gcor.s * (GMEANS - mean(GMEANS)))/sqrt(var(GMEANS))+sqrt(1-gcor.s^2)*rnorm(ngrp)
  # Minor change here to reverse code Y.G to get negative values.
  # Note that this corresponds to the creation of VAR2 above in that VAR2
  # which is Y.G's referent is the negative value created in reference to VAR1
  Y.G<-rep(Y.G,each=gsize)
  if(gcor<0)
    Y.G<-Y.G*-1
  Y<-sqrt(icc1y)*Y.G+sqrt(1-icc1y)*VAR2
  #Return the raw variables
  return(data.frame(GRP=GID,X,Y))
}
#sobel####
sobel<-function(pred,med,out){
	NEWDAT<-data.frame(pred=pred,med=med,out=out)
	NEWDAT<-na.exclude(NEWDAT)
	model1<-lm(out~pred,data=NEWDAT)
	model2<-lm(out~pred+med,data=NEWDAT)
	model3<-lm(med~pred,data=NEWDAT)
	mod1.out<-summary(model1)$coef
	mod2.out<-summary(model2)$coef
	mod3.out<-summary(model3)$coef
	indir<-mod3.out[2,1]*mod2.out[3,1]
	effvar<-(mod3.out[2,1])^2*(mod2.out[3,2])^2+(mod2.out[3,1])^2*(mod3.out[2,2])^2
	serr<-sqrt(effvar)
	zvalue=indir/serr
out<-list('Mod1: Y~X'=mod1.out,'Mod2: Y~X+M'=mod2.out,'Mod3: M~X'=mod3.out,
	Indirect.Effect=indir,SE=serr,z.value=zvalue,N=nrow(NEWDAT))
return(out)
}

#waba####
waba<-function(x, y, grpid)
{
	SMAT <- na.exclude(data.frame(grpid, x, y))
	TEMP <- aggregate(SMAT[, 2:3], list(SMAT$grpid), mean)
	names(TEMP) <- c("grpid", "MEANx", "MEANy")
	SMAT <- merge(SMAT, TEMP, by = "grpid")
	SMAT$WITHx <- SMAT$x - SMAT$MEANx
	SMAT$WITHy <- SMAT$y - SMAT$MEANy
	ans <- matrix(0, 1, 7)
	dimnames(ans) <- list(NULL, c("RawCorr", "EtaBx", "EtaBy", "CorrB",
		"EtaWx", "EtaWy", "CorrW"))
	ans[1, 1] <- cor(SMAT[, 2], SMAT[, 3])
	ans[1, 2] <- cor(SMAT[, 2], SMAT[, 4])
	ans[1, 3] <- cor(SMAT[, 3], SMAT[, 5])
	ans[1, 4] <- cor(SMAT[, 4], SMAT[, 5])
	ans[1, 5] <- cor(SMAT[, 2], SMAT[, 6])
	ans[1, 6] <- cor(SMAT[, 3], SMAT[, 7])
	ans[1, 7] <- cor(SMAT[, 6], SMAT[, 7])
	ngroups <- nrow(TEMP)
	nobs <- nrow(SMAT)
	return(list(Cov.Theorem = data.frame(ans), n.obs = nobs, n.grps = 
		ngroups))
}
