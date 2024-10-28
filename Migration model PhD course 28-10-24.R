#### Preamble ####
rm(list=ls())
library(nls2)
library(rptR)
library(MASS)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(lme4)
library(dplyr)
library(vioplot)
#setwd("M:\\Documents\\Migration")
pal_d <- brewer.pal(8,"Dark2")
pal_v <- viridis(3)

nnzero <- function(x){
  length(which(x!=0))
}

mysample <- function(x,size,replace=T,prob=NULL){
  if(length(x)==1) { rep(x,times=round(size)) }
  else { sample(x,size=round(size),replace=replace,prob=prob)}
}

## is.mature: Quick function to check whether an individual has reached age of maturity.
## Arguments:
#  x  : A vector of ages, can include NAs.
#  a  : Age at maturity.
#  Outputs: A boolean vector of length x.
is.mature <- function(x,a){
  x[is.na(x)] <- -1
  out <- ifelse(x>=a,TRUE,FALSE)
  return(out)
}

## ddfun: Quickly plot survival as a function of density for arbitrarily many m2 parameters
## stars: a matrix of dimensions N x reps of population densities, to plot along the lines.
ddfun <- function(m2,type="exponential",title=NA,stars=matrix(NA,3,3),xmax=2,pal=pal_v,legend=T){
  xs <- seq(0,xmax,by=0.002)
  survs <- numeric(length(xs))
  plot(xs,rep(1,length(xs)),type="n",ylim=c(0,1),col=pal[1],xlab=bquote(italic(n[i]/K)),cex.lab=1.2,ylab="Survival probability",main=title)
  if(length(m2)>1){
    for(i in 1:length(m2)){
      if(type=="exponential"){
        survs <- dexp(xs,m2[i])/m2[i]
        lines(xs,survs,col=pal[i],lwd=2)
        if(sum(stars[i,],na.rm=T)>0){
          for(j in 1:length(stars[1,])){ # For each replicate, plot an open circle
            points(stars[i,j],dexp(stars[i,j],m2[i])/m2[i],pch=1,col=alpha(pal[i],0.5),cex=1.5)
          }
          points(mean(stars[i,]),dexp(mean(stars[i,]),m2[i])/m2[i],pch=20,col=pal[i],cex=1.8)
        }
      } else if(type=="linear"){
        survs <- 1-xs*m2[i]
        lines(xs,survs,col=pal[i],lwd=2)
        if(sum(stars,na.rm=T)>0){
          points(stars[i],1-xs*m2[i],pch=20,col=pal[i])
        }
      }
      #print(survs[101])
    }
  }
  abline(v=1,lty=2)
  if(legend) legend("bottomleft",col=pal,pch=rep(20,length(m2)),legend=as.fractions(m2),title=expression(italic(gamma)[n]),bty="n",cex=1,pt.cex=1.3)
}
#####

#### Functions for wrangling and plotting results ####

### Popstats ###
##  Function to summarize output from simulations into a plottable format (takes some time for very long simulations)
##  Arguments:
#   local    : A simulation result (created locally)
#   saved    : A simulation result (from saved file)
#   max.age  : Number of age classes to bin individuals in for age distribution (individuals of age > max.age are binned into age=max.age)
##  Outputs
#   A list of two elements:
#   [[1]] EVOL: Storage array of format:
#   Row 1:      Mean breeding values
#   Row 2:      Sd breeding values
#   Row 3:(N+2):Frequencies of destination genes
#   Row (N+3):(N+2+max.age):Age distribution
#   per population for each Patch (2nd dim), recorded year (3rd dim) and Replicate (4th dim)
#   [[2]] LIAB: Storage array of format:
#   Row 1:      Mean liability
#   Row 2:      Sd liability
#   Row 3:      Fraction L>0
#   per population for each patch (2nd dim), recorded year (3rd dim) and Replicate (4th dim)
#   [[3]] AGE: Storage array of subpopulation mean age per patch (1st dim), recorded year (2nd dim) and replicate (3rd dim).
#   [[4]] EXT: Storage array of format
#   Row 1: Year of first extinction (pre-breeding population census of 0)
#   Row 2: Most recently recorded mean BV prior to extinction.
#   Row 3: Most recently recorded sd BV prior to extinction.
#   per site (2nd dim) and replicate (3rd dim)
popstats <- function(local=NULL,saved=NULL,max.age=10){
  if(is.null(local)){
    output <- readRDS(saved)
  } else {
    output <- local
  }
  
  popstorage <- output[[1]]
  dyn <- output[[2]]
  N <- dim(dyn)[2] # Number of subpopulations
  K <- dim(popstorage)[1]/N # Max number of individuals per subpopulation (carrying capacity)
  t <- dim(popstorage)[3] # Number of recorded timesteps
  
  reps <- dim(dyn)[3] # Number of replicate simulations
  rec.int <- dim(dyn)[1]/(2*t) # Interval of recording
  years <- (1:t)*rec.int # The years at which recordings took place
  
  ## Create storage matrices ##
  # Mean and sd gene values per pop per recording
  evol <- array(NA,dim=c(N+2+max.age,N,length(years),reps)) # Dims: 1) Stats (see below). 2) Patches. 3) Years. 4) Reps
  liab <- array(NA,dim=c(3,N,length(years),reps)) # Rows: 1) Mean L 2) Sd L 3) Fraction L>0 - per patch, year, rep.
  age <- array(NA,dim=c(N,length(years),reps)) # Subpopulation mean age. Dims: 1) Patch, 2) Years, 3) Rep
  ext <- array(NA,dim=c(3,N,reps)) # Subpop extinction times (1st row), most recent mean (2nd row) and sd (3rd row) BVs pre-extinction, per patch (cols) and rep (3rd dim)
  # Stats format for evol:
  # 1: Mean breeding value
  # 2: Sd breeding value
  # 3,4,5 (or "3:(N+2)"): Frequencies of destination genes
  # 6-16 (or "(N+3):(N+2+max.age)"): Age distribution
  t <- Sys.time()
  for(r in 1:reps){
    for(i in 1:length(years)){
      popstorage[which(popstorage[,3,i,r]>max.age),3,i,r] <- max.age # Final age class is "10 or older".
      for(j in 1:N){
        inds <- (1:K)+(j-1)*K
        evol[1,j,i,r] <- mean(popstorage[inds,1,i,r],na.rm=T) # Mean breeding value
        evol[2,j,i,r] <- sd(popstorage[inds,1,i,r],na.rm=T) # Sd breeding value
        liab[1,j,i,r] <- mean(popstorage[inds,6,i,r],na.rm=T) # Mean liability
        liab[2,j,i,r] <- sd(popstorage[inds,6,i,r],na.rm=T) # Sd liability
        liab[3,j,i,r] <- sum(popstorage[inds,6,i,r]>0,na.rm=T)/sum(!is.na(popstorage[inds,6,i,r])) # Fraction of migrants
        evol[3:(N+2),j,i,r] <- c(tabulate(popstorage[inds,2,i,r]),rep(0,times=N-length(tabulate(popstorage[inds,2,i,r]))))/sum(!is.na(popstorage[inds,1,i,r])) #Or divide by dyn[(i*2)-1,j,r]  # Destination genes distribution
        #print(length(c(tabulate(popstorage[(1:K)+(j-1)*K,3,i,r]),rep(0,times=max.age-length(tabulate(popstorage[(1:K)+(j-1)*K,3,i,r]))))))
        #print(length((N+3):(N+2+max.age))) # These are just for checking that everything runs smoothly.
        evol[(N+3):(N+2+max.age),j,i,r] <- c(tabulate(popstorage[inds,3,i,r])/sum(!is.na(popstorage[inds,1,i,r])),rep(0,times=max.age-length(tabulate(popstorage[inds,3,i,r])))) # Age distribution
        evol[is.infinite(evol)] <- NA # Find a better solution for this!
        liab[is.infinite(liab)] <- NA 
        age[j,i,r] <- mean(popstorage[((j-1)*K+1):(j*K),3,i,r],na.rm=T)
      }
    }
    for(j in 1:N){
    # Use dyn instead  ext[r,j] <- which(is.na(statstmp[[2]][3,3,,1]))[min(which(diff(which(is.na(statstmp[[2]][3,3,,1])))==1))]
      summers <- seq(1,length(dyn[,1,1]),by=2)
      if(any(dyn[summers,j,r]==0)){
        ext[1,j,r] <- min(which(dyn[summers,j,r]==0)) # Year of first extinction
        ext[2,j,r] <- evol[1,j,max(which(years<ext[1,j,r])),r] # Most recently recorded mean breeding value
        ext[3,j,r] <- evol[2,j,max(which(years<ext[1,j,r])),r]
      }
    }
    ##
    
    print(r)
  }
  
  print(Sys.time()-t)
  return(list(evol,liab,age,ext))
}

### Popplots ###
## Function to plot simulation results.
## Arguments:
#  local    : A stats file (generated by popstats) in the local environment
#  saved    : A stats file (generated by popstats) that can be read from disk
#  mainfile : A simulation result (generated by sim), created locally or read from file.
#  nview    : Number of replicates to choose when plotting fewer than all replicates. Default 5.
#  plot     : Shortnames for the different plots wanted: "all" (default),
#             "indBVs","indLs","meanBVs","meanLs","popsize","endpopsize","winteruse","agedist","Dest","meancorr","beforecorr","aftercorr",
#             "EndBVs","EndBVviol","EndLs","EndLiabviol","EndMeans","EndSDs","EndDests","EndDensity",
#             "preshockBVs","preshockBVviol","preshockLs","preshockLviol","preshockMeans","preshockSDs","preshockDests","preshockDensity".
#  pal      : Color palette to use: "vir" (viridis, default), or "dark" (RColorbrewer's Dark2). Both should be colorblind friendly.
#  legend   : Setting FALSE can suppress some more annoying legends: EndDests, (add more here as they are implemented)
#  shocktime: For endpopsize if running a short simulation after a long one - 
#  sitenames: Whether sites have "geographic" names (default) or "alphabetic"
popplots <- function(local=NULL,saved="",mainfile=NULL,nview=5,plot="all",pal="vir",legend=T,shocktime=10000,sitenames="geographic",...){
  st <- Sys.time()
  if(is.null(local)){
    stats <- readRDS(saved)
  } else {
    stats <- local
  }
  evol <- stats[[1]]
  liab <- stats[[2]]
  
  if(is.null(mainfile)){
    tmp <- strsplit(saved," stats")[[1]][1]
    output <- readRDS(paste0(tmp,".R"))
  } else {
    output <- mainfile
  }
  
  popstorage <- output[[1]]
  dyn <- output[[2]]
  N <- dim(dyn)[2] # Number of subpopulations
  K <- dim(popstorage)[1]/N # Max number of individuals per subpopulation (carrying capacity)
  t <- dim(popstorage)[3] # Number of recorded timesteps
  
  reps <- dim(dyn)[3] # Number of replicate simulations
  rec.int <- dim(dyn)[1]/(2*t) # Interval of recording
  years <- (1:t)*rec.int # The years at which recordings took place
  max.age <- dim(evol)[1]-N-2 # Number of age classes individuals are binned into (for age distribution)
  
  pal <- switch(pal,"vir"=viridis(N),"dark"=brewer.pal(N,"Dark2"))
  if(sitenames=="geographic"){
    sites <- switch(as.character(N),
                    "2"=c("North","South"),
                    "3"=c("North","Middle","South"),
                    "4"=c("North","Middle-N","Middle-S","South"),
                    "5"=c("North","Middle-N","Middle","Middle-S","South"),
                    "6"=c("North","","Middle-N","Middle-S","","South"),
                    "7"=c("North","","Middle-N","Middle","Middle-S","","South"),
                    "8"=c("North","","Middle-N","","","Middle-S","","South"),
                    "9"=c("North","","Middle-N","","Middle","","Middle-S","","South"))
  } else if(sitenames=="alphabetic"){
    sites <- LETTERS[1:N]
  }
  if(N>8){
    pal <- c(pal,brewerpal(8,"Set2"))
  }
  
  # Other info that will be useful to generate
  ncorrs <- choose(N,2) # Calculate how many combinations of patches there are
  viewr <- sample(reps,min(nview,reps)) # If there are too many reps to visualize easily, only plot 5 of them.
  shocks <- output[[3]]$shock.size
  print(Sys.time()-st)

  # Mean gene values from random 3 replicates over time. Line type is replicate. Shaded bands are population SD.
  if(any(plot=="all") | any(plot=="indBVs")){
    par(mfrow=c(2,2),mar=c(4,4,2,1))
    
    plot(1:t,rep(0,t),type="n",xlab="Year",xaxt="n",ylab="Breeding values \U00B1 within-pop. sd",ylim=c(-5,5),cex.lab=1.2)
    axis(1,at=seq(0,t,by=1000/rec.int),labels=c(0,years[which(years%%1000==0)]))
    viewr <- sample(reps,size=min(reps,3),replace=F)
    for(r in 1:length(viewr)){
      for(j in 1:N){
        lowerbounds <- evol[1,j,,viewr[r]] - evol[2,j,,viewr[r]]
        upperbounds <- evol[1,j,,viewr[r]] + evol[2,j,,viewr[r]]
        polygon(c(1:t,rev(1:t)),c(lowerbounds,rev(upperbounds)),col=alpha(pal[j],0.1),border=NA)
        lines(1:t,evol[1,j,,viewr[r]],col=pal[j],lty=r)
      }
    }
    if(legend) legend("bottomleft",legend=paste("Site",sites[1:N]),lty=1,col=pal[1:N],bg=alpha("White",alpha=0.5))
    abline(h=0,lty=2)
  }
  
  # Mean liabilities from 3 random replicates over time. Line type is replicate. Shaded bands are population SD.
  if(any(plot=="all") | any(plot=="indLs")){
    plot(1:t,rep(0,t),type="n",xlab="Year",xaxt="n",ylab="Liabilities \U00B1 within-pop. sd",ylim=c(-5,5),cex.lab=1.2)
    axis(1,at=seq(0,t,by=1000/rec.int),labels=c(0,years[which(years%%1000==0)]))
    for(r in 1:length(viewr)){
      for(j in 1:N){
        lowerbounds <- liab[1,j,,viewr[r]] - liab[2,j,,viewr[r]]
        upperbounds <- liab[1,j,,viewr[r]] + liab[2,j,,viewr[r]]
        polygon(c(1:t,rev(1:t)),c(lowerbounds,rev(upperbounds)),col=alpha(pal[j],0.1),border=NA)
        lines(1:t,liab[1,j,,viewr[r]],col=pal[j],lty=r)
      }
    }
    if(legend) legend("bottomleft",legend=paste("Site",sites[1:N]),lty=1,lwd=1,col=pal[1:N],bg=alpha("white",alpha=0.5))
    abline(h=0,lty=2)
  }
  
  #Mean gene values across replicates over time. Shaded bands are standard error of the mean across all reps
  if(any(plot=="all") | any(plot=="meanBVs")){
#    if(any(plot!="all") & any(plot=="meanBVs")){
 #     par(mfrow=c(2,2),mar=c(4,4,2,1))
  #  }
    sigma_e <- output[[3]]$sigmaE_coarse
    plot(1:t,rep(0,t),type="n",xlab="Year",xaxt="n",ylab="Mean breeding values \U00B1 among-pop. SE",ylim=c(-5,5),cex.lab=1.2)
    axis(1,at=seq(0,t,by=1000/rec.int),labels=c(0,years[which(years%%1000==0)]))
    for(j in 1:N){
      means <- apply(evol[1,j,,],FUN=mean,MARGIN=1,na.rm=T)
      lowerbounds <- means - apply(evol[1,j,,],MARGIN=1,FUN=sd,na.rm=T)/sqrt(reps) # Standard error of the mean
      upperbounds <- means + apply(evol[1,j,,],MARGIN=1,FUN=sd,na.rm=T)/sqrt(reps) # Standard of the mean
      lines(1:t,apply(evol[1,j,,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=2)
      polygon(c(1:t,rev(1:t)),c(lowerbounds,rev(upperbounds)),col=alpha(pal[j],alpha=0.1),border=NA)
    }
    ys <- seq(-3*sigma_e,3*sigma_e,by=0.01)
    xs <- dnorm(ys,0,sigma_e)
    lines(xs*t/4,ys)
    abline(h=0,lty=2,lwd=2)
    if(legend) legend("bottomleft",title="Site",legend=sites[1:N],col=pal[1:N],lty=1,lwd=2,bty="n",seg.len=1)
  }
  
  # Mean liabilities across replicates over time. Shaded bands are standard error of the mean across all reps
  if(any(plot=="all") | any(plot=="meanLs")){
    plot(1:t,rep(0,t),type="n",xlab="Year",xaxt="n",ylab="Mean liabilities \U00B1 among pop. SE",ylim=c(-5,5),cex.lab=1.2)
    axis(1,at=seq(0,t,by=1000/rec.int),labels=c(0,years[which(years%%1000==0)]))
    for(j in 1:N){
      means <- apply(liab[1,j,,],FUN=mean,MARGIN=1,na.rm=T)
      #lowerbounds <- apply(liab[1,j,,],MARGIN=1,FUN=quantile,na.rm=T,probs=0.25) # 50 % quantile
      #upperbounds <- apply(liab[1,j,,],MARGIN=1,FUN=quantile,na.rm=T,probs=0.75) # 50 % quantile
      lowerbounds <- means - apply(liab[1,j,,],MARGIN=1,FUN=sd,na.rm=T)/sqrt(reps) # Standard error of the mean
      upperbounds <- means + apply(liab[1,j,,],MARGIN=1,FUN=sd,na.rm=T)/sqrt(reps) # Standard of the mean
      lines(1:t,apply(liab[1,j,,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=2)
      polygon(c(1:t,rev(1:t)),c(lowerbounds,rev(upperbounds)),col=alpha(pal[j],alpha=0.1),border=NA)
    }
    sigma_e <- output[[3]]$sigmaE_coarse
    ys <- seq(-3*sigma_e,3*sigma_e,by=0.01)
    xs <- dnorm(ys,0,sigma_e)
    lines(xs*t/4,ys)
    abline(h=0,lty=2,lwd=2)
    if(legend) legend("bottomleft",title="Site",legend=sites[1:N],col=pal[1:N],lty=1,lwd=2,bty="n",seg.len=1)
  }
  
  # Population sizes over time, means and SDs across all reps.
  if(any(plot=="all") | any(plot=="popsize")){
    par(mfrow=c(2,1),mar=c(4,4,1,0.5))
    plot(1:(rec.int*t),seq(1,K,length.out=rec.int*t),type="n",ylab="Local population size",xlab="Year",ylim=c(0,max(dyn,na.rm=T)),xaxt="n",cex.lab=1.2)
    axis(1,at=seq(0,rec.int*t,by=1000)-1,labels=c(0,years[which(years%%1000==0)]))
    for(j in 1:N){
      lines(1:(rec.int*t),apply(dyn[seq(1,(rec.int*t*2)-1,by=2),j,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=1,lty=2) # Summer
      lines(1:(rec.int*t),apply(dyn[seq(2,(rec.int*t*2),by=2),j,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=1) # Winter
      lowerbounds.b <- apply(dyn[seq(1,(rec.int*t*2)-1,by=2),j,],FUN=mean,MARGIN=1,na.rm=T) - apply(dyn[seq(1,(rec.int*t*2)-1,by=2),j,],FUN=sd,MARGIN=1,na.rm=T) # Summer
      upperbounds.b <- apply(dyn[seq(1,(rec.int*t*2)-1,by=2),j,],FUN=mean,MARGIN=1,na.rm=T) + apply(dyn[seq(1,(rec.int*t*2)-1,by=2),j,],FUN=sd,MARGIN=1,na.rm=T)
      polygon(c(1:(rec.int*t),rev(1:(rec.int*t))),c(lowerbounds.b,rev(upperbounds.b)),col=alpha(pal[j],alpha=0.1),border=NA)
      lowerbounds.w <- apply(dyn[seq(2,rec.int*t*2,by=2),j,],FUN=mean,MARGIN=1,na.rm=T) - apply(dyn[seq(2,rec.int*t*2,by=2),j,],FUN=sd,MARGIN=1,na.rm=T) # Winter
      upperbounds.w <- apply(dyn[seq(2,rec.int*t*2,by=2),j,],FUN=mean,MARGIN=1,na.rm=T) + apply(dyn[seq(2,rec.int*t*2,by=2),j,],FUN=sd,MARGIN=1,na.rm=T)
      polygon(c(1:(rec.int*t),rev(1:(rec.int*t))),c(lowerbounds.w,rev(upperbounds.w)),col=alpha(pal[j],alpha=0.1),border=NA)
    }
    abline(h=K,lty=5) # Breeding season carrying capacity
    if(legend) legend("topleft",legend=c("Breeding","Non-breeding",paste(sites[1:N],"site")),lty=c(2,1,rep(0,N)),pch=c(NA,NA,rep(15,N)),lwd=2,col=c("Black","Black",pal[1:N]),bg=alpha("white",alpha=0.5))
  }
  
  # Population sizes +- sd across reps, zoom in on last 75 time steps. Change hashtags for 175 time steps. 
  if(any(plot=="all") | any(plot=="endpopsize")){
    ts <- seq(rec.int*t*2-350,rec.int*t*2-1,by=2)
    plot(1:length(ts),seq(1,K,length.out=length(ts)),type="n",ylab="Local population size",xlab="Year",ylim=c(0,max(dyn,na.rm=T)),xaxt="n",cex.lab=1.2)
    axis(1,at=c(1,(1:12)*25),labels=rec.int*t-(25*(13:1)))
    #if(rec.int*t>100){
    #  axis(1,at=c(1,25,50,75),labels=((rec.int*t)-25*(3:0)))
    #} else {
    #  axis(1,at=c(1,25,50,75),labels=shocktime+c(25,50,75,100))
    #}
    for(j in 1:N){
      lines(1:length(ts),apply(dyn[ts+1,j,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=2,lty=2) # Summer
      lines(1:length(ts),apply(dyn[ts+2,j,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=2) # Winter
      lowersum <- apply(dyn[ts+1,j,],FUN=mean,MARGIN=1,na.rm=T) - apply(dyn[ts+1,j,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(reps) # Summer
      uppersum <- apply(dyn[ts+1,j,],FUN=mean,MARGIN=1,na.rm=T) + apply(dyn[ts+1,j,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(reps) # Summer
      lowerwin <- apply(dyn[ts+2,j,],FUN=mean,MARGIN=1,na.rm=T) - apply(dyn[ts+2,j,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(reps) # Winter
      upperwin <- apply(dyn[ts+2,j,],FUN=mean,MARGIN=1,na.rm=T) + apply(dyn[ts+2,j,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(reps) # Winter
      polygon(c(1:length(ts),length(ts):1),c(lowersum,rev(uppersum)),col=alpha(pal[j],alpha=0.1),border=NA)
      polygon(c(1:length(ts),length(ts):1),c(lowerwin,rev(upperwin)),col=alpha(pal[j],alpha=0.1),border=NA)
      #for(r in viewr){
      #  lines(1:75,dyn[ts,j,r],col=alpha(pal[j],0.1),lty=2) # Summer
      #  lines(1:75,dyn[ts+1,j,r],col=alpha(pal[j],0.1)) # Winter
      #}
    }
    abline(h=K,lty=5,lwd=2) # Breeding season carrying capacity
    abline(v=25,lty=3,lwd=2) # Shock time
    if(legend){
      legend("topright",legend=c("Breeding", "Non-breeding",paste(sites[1:N],"site")),lty=c(2,1,rep(0,N)),pch=c(NA,NA,rep(15,N)),lwd=2,col=c("Black","Black",pal[1:N]),bg=alpha("white",alpha=0.5))
    }
  }
  
  # Frequencies of individuals overwintering where. Each point is a replicate - max 5 random replicates shown.
  if(any(plot=="all") | any(plot=="winteruse")){
    par(mfrow=c(N,1),mar=c(3.5,4,0.5,0.5))
    means <- array(NA,dim=c(N,N,length(years),reps)) # Dims: From, to, timestep, rep
    for(f in 1:N){
      inds <- ((f-1)*K+1):(f*K)
      for(y in 1:t){
        for(r in 1:length(viewr)){
          means[f,,y,viewr[r]] <- c(tabulate(popstorage[inds,5,y,viewr[r]]),rep(0,times=N-length(tabulate(popstorage[inds,5,y,viewr[r]]))))/sum(!is.na(popstorage[inds,5,y,viewr[r]]))
          if(sum(means[f,,y,viewr[r]],na.rm=T)>1.001 | sum(means[f,,y,viewr[r]],na.rm=T)<0.999){
            print(paste("Overwintering frequencies do not sum to 1, population",f,", year",y*rec.int,", rep",viewr[r]))
            print(sum(means[f,,y,viewr[r]],na.rm=T),digits=10)
          }
        }
      }
      plot(1:t,rep(1,t),type="n",ylim=c(0,1),ylab=paste("Freq of winter locations from",sites[f],"pop"),xlab="",xaxt="n")
      axis(1,at=seq(1,t+1,by=1000/rec.int)-1,labels=c(0,years[which(years%%1000==0)]))
      for(j in 1:N){
        for(r in 1:reps){
          points(jitter(1:t),means[f,j,,r],col=alpha(pal[j],0.5),pty=r)
        }
      }
    }
    if(legend) legend("bottomleft",title="Non-breeding season site",legend=sites[1:N],col=pal,pch=rep(15,N),bg=alpha("white",0.5))
    mtext("Year",side=1,line=2.5)
  }
  
  # Age distribution from up to 5 random replicates.
  if(any(plot=="all") | any(plot=="agedist")){ 
    par(mfrow=c(length(viewr),1),mar=c(3.5,4.2,0.5,0.3),oma=c(0,0,1.5,0))
    for(r in 1:length(viewr)){
      plot(1:max.age,rep(0.1,max.age),type="n",xlab="",ylab="Frequency",cex.lab=1.2,cex.axis=1.2,ylim=c(0,max(evol[(N+3):(N+2+max.age),,,],na.rm=T)))
      for(j in 1:N){
        lines((1:max.age)+0.1*(j-1),evol[(N+3):(N+2+max.age),j,t,viewr[r]],col=pal[j],type="h",lwd=2)
        abline(v=mean(popstorage[((j-1)*K+1):(j*K),3,t,viewr[r]],na.rm=T),lty=2,col=pal[j])
      }
      text(max.age/2,0.25,paste("Replicate",viewr[r]),cex=1.5)
    }
    mtext("Age class",side=1,line=2)
    mtext("   Age distribution, post shock",side=3,line=-0.3,cex=1.2,outer=T)
    if(legend) legend("topright",title="Home site",legend=sites[1:N],col=pal,pch=rep(15,N),bg=alpha("white",0.5))
  }
  
  # Destination gene evolution and migrant fraction over time.
  if(any(plot=="all") | any(plot=="Dest")){
    par(mfrow=c(N,1),mar=c(3.5,4,0.5,0.5))
    for(i in 1:N){
      plot(1:t,rep(1,t),type="n",ylim=c(0,1),ylab=paste("Freq of Dest. gene in",sites[i],"pop"),xlab="",xaxt="n",cex.lab=1+(4-N)/10)
      axis(1,at=seq(1,t+1,by=1000/rec.int)-1,labels=c(0,years[which(years%%1000==0)]))
      for(j in 1:N){
        lines(1:t,apply(evol[3+(j-1),i,1:t,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=2)
        lowerbounds <- apply(evol[3+j-1,i,1:t,],FUN=mean,MARGIN=1,na.rm=T) - apply(evol[3+j-1,i,1:t,],FUN=sd,MARGIN=1,na.rm=T)/reps
        upperbounds <- apply(evol[3+j-1,i,1:t,],FUN=mean,MARGIN=1,na.rm=T) + apply(evol[3+j-1,i,1:t,],FUN=sd,MARGIN=1,na.rm=T)/reps
        polygon(c(1:t,rev(1:t)),c(lowerbounds,rev(upperbounds)),col=alpha(pal[j],alpha=0.1),border=NA)
      }
      lines(1:t,apply(liab[3,i,1:t,],FUN=mean,MARGIN=1,na.rm=T),col="Black",lty=2,lwd=2) # Proportion migrant
      lowerbounds <- apply(liab[3,i,1:t,],FUN=mean,MARGIN=1,na.rm=T) - apply(liab[3,i,1:t,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(reps)
      upperbounds <- apply(liab[3,i,1:t,],FUN=mean,MARGIN=1,na.rm=T) + apply(liab[3,i,1:t,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(reps)
      polygon(c(1:t,rev(1:t)),c(lowerbounds,rev(upperbounds)),col=alpha("Black",alpha=0.1),border=NA)
    }
    mtext("Year",side=1,line=2.5)
    legend("topleft",title="Home site",legend=sites[1:N],col=pal,pch=rep(15,N),bg=alpha("white",0.5),cex=1+(1-N)/20)
    if(legend) legend("topright",legend=c("Allele frequencies","Proportion migrant"),lty=1:2,lwd=2,bg=alpha("white",0.5))
  }
  
  # Correlations between population sizes
  if(any(plot=="all") | any(plot=="meancorr")){
    par(mfrow=c(3,1),mar=c(3,4,1.5,0.5))
    # Same kind of loop as for networks, plotting all combinations
    scorrs <- wcorrs <- matrix(NA,ncorrs,reps)
    names <- character(ncorrs)
    i <- counter <- 1
    summers <- seq(1,rec.int*t-1,by=2)
    winters <- summers+1
    plot(1:ncorrs,rep(0,ncorrs),type="n",xaxt="n",xlim=c(ifelse(N>3,1,0.5),ifelse(N>3,ncorrs,ncorrs+0.5)),ylim=c(-1,1),ylab="Correlation",xlab="",main="All years",bty="L")
    while(i<N){
      for(j in (i+1):N){
        polygon(c(counter-0.3,counter+0.3,counter+0.3,counter-0.3),c(-1.2,-1.2,1.2,1.2),col=alpha("Light grey",alpha=min(1,shocks[i]+shocks[j])),border=NA)
        scorrs[counter,] <- diag(cor(log(dyn[summers,i,]),log(dyn[summers,j,])))
        wcorrs[counter,] <- diag(cor(log(dyn[winters,i,]),log(dyn[winters,j,])))
        names[counter] <- paste0(LETTERS[i],LETTERS[j])
        vioplot(scorrs[counter,],add=T,at=counter-0.1,col=alpha("Dark green",0.2),wex=0.25,border="Dark green",lineCol=NA)
        vioplot(wcorrs[counter,],add=T,at=counter+0.1,col=alpha("Blue",0.2),wex=0.25,border="Blue",lineCol=NA)
        #points(counter+0.1*(1:length(viewr))-0.1*mean(c(1,length(viewr))),scorrs[counter,viewr],pch=2,col="Dark green")
        #points(counter+0.1*(1:length(viewr))-0.1*mean(c(1,length(viewr))),wcorrs[counter,viewr],pch=4,col="Blue")
        counter <- counter+1
      }
      i <- i+1
    }
    axis(1,at=1:ncorrs,labels=names,cex=1.2)
    abline(h=0,lty=2)
    #legend("bottomleft",legend=c("Summer","Winter"),col=c("Dark green", "Blue"),pch=c(2,4))
    if(legend) legend("bottomleft",legend=c("Breeding","Non-breeding"),col=c("Dark green","Blue"),pt.bg=c(alpha("Dark green",0.2),alpha("Blue",0.2)),pch=22)
  }
  
  # Correlations between population sizes before shock
  if(any(plot=="all") | any(plot=="beforecorr")){
    # Same kind of loop as for networks, plotting all combinations
    scorrs <- wcorrs <- matrix(NA,ncorrs,reps)
    i <- counter <- 1
    summers <- seq(rec.int*(t-2)+1,rec.int*(t-1)-1,by=2)
    winters <- summers+1
    plot(1:ncorrs,rep(0,ncorrs),type="n",xaxt="n",xlim=c(ifelse(N>3,1,0.5),ifelse(N>3,ncorrs,ncorrs+0.5)),ylim=c(-1,1),ylab="Correlation",xlab="",main="50 years before shock",bty="L")
    while(i<N){
      for(j in (i+1):N){
        polygon(c(counter-0.3,counter+0.3,counter+0.3,counter-0.3),c(-1.2,-1.2,1.2,1.2),col=alpha("Light grey",alpha=min(1,shocks[i]+shocks[j])),border=NA)
        scorrs[counter,] <- diag(cor(log(dyn[summers,i,]),log(dyn[summers,j,])))
        wcorrs[counter,] <- diag(cor(log(dyn[winters,i,]),log(dyn[winters,j,])))
        names[counter] <- paste0(LETTERS[i],LETTERS[j])
        vioplot(scorrs[counter,],add=T,at=counter-0.1,col=alpha("Dark green",0.2),wex=0.25,border="Dark green",lineCol=NA)
        vioplot(wcorrs[counter,],add=T,at=counter+0.1,col=alpha("Blue",0.2),wex=0.25,border="Blue",lineCol=NA)
        #points(counter+0.1*(1:length(viewr))-0.1*mean(c(1,length(viewr))),scorrs[counter,viewr],pch=2,col="Dark green")
        #points(counter+0.1*(1:length(viewr))-0.1*mean(c(1,length(viewr))),wcorrs[counter,viewr],pch=4,col="Blue")
        counter <- counter+1
      }
      i <- i+1
    }
    axis(1,at=1:ncorrs,labels=names,cex=1.2)
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Breeding","Non-breeding"),col=c("Dark green","Blue"),pt.bg=c(alpha("Dark green",0.2),alpha("Blue",0.2)),pch=22)
  }
  
  # Correlations between population sizes after shock
  if(any(plot=="all") | any(plot=="aftercorr")){
    # Same kind of loop as for networks, plotting all combinations
    scorrs <- wcorrs <- matrix(NA,ncorrs,reps)
    names <- character(choose(N,2))
    i <- counter <- 1
    summers <- seq(rec.int*(t-1)+1,(rec.int*t)-1,by=2)
    winters <- summers+1
    plot(1:ncorrs,rep(0,ncorrs),type="n",xaxt="n",xlim=c(ifelse(N>3,1,0.5),ifelse(N>3,ncorrs,ncorrs+0.5)),ylim=c(-1,1),ylab="Correlation",xlab="",main="50 years after shock",bty="L")
    while(i<N){
      for(j in (i+1):N){
        polygon(c(counter-0.3,counter+0.3,counter+0.3,counter-0.3),c(-1.2,-1.2,1.2,1.2),col=alpha("Light grey",alpha=min(1,shocks[i]+shocks[j])),border=NA)
        scorrs[counter,] <- diag(cor(log(dyn[summers,i,]),log(dyn[summers,j,])))
        wcorrs[counter,] <- diag(cor(log(dyn[winters,i,]),log(dyn[winters,j,])))
        names[counter] <- paste0(LETTERS[i],LETTERS[j])
        vioplot(scorrs[counter,],add=T,at=counter-0.1,col=alpha("Dark green",0.2),wex=0.25,border="Dark green",lineCol=NA)
        vioplot(wcorrs[counter,],add=T,at=counter+0.1,col=alpha("Blue",0.2),wex=0.25,border="Blue",lineCol=NA)
        #points(counter-0.05*(1:length(viewr)),scorrs[counter,viewr],pch=2,col="Dark green")
        #points(counter+0.05*(1:length(viewr)),wcorrs[counter,viewr],pch=4,col="Blue")
        counter <- counter+1
      }
      i <- i+1
    }
    axis(1,at=1:ncorrs,labels=names,cex=1.2)
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Breeding","Non-breeding"),col=c("Dark green","Blue"),pt.bg=c(alpha("Dark green",0.2),alpha("Blue",0.2)),pch=22)
  }
  
  ### Average over all recordings of summary stats for last 1000 years.
  ts <- max(1,(t-1000%/%rec.int)):t # 20 recordings if rec.int=50
  
  if(any(plot=="all") | any(plot=="avgBVs")){
    par(mfrow=c(2,4),mar=c(4,4,2,1),oma=c(0,0,1.5,0))
    
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Breeding values",ylim=range(evol[1:2,,t,],na.rm=T),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),apply(evol[1,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.25),apply(evol[2,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),cex=1.8)
    }
    abline(h=0,lty=2)
    mtext("Average over recordings from last 1000 years",outer=T,cex=1.4,line=-0.5)
    if(legend) legend("bottomleft",legend=c("Mean","Sd"),bty="n",pch=c(16,1),col=alpha("Black",0.5),title="Within-pop BV",pt.cex=1.8)
  }
  
  if(any(plot=="all") | any(plot=="avgBVviol")){
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,1,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Breeding values",bty="L")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,1,t,r]))){
          vioplot(popstorage[inds,1,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,na.rm=T,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
  }
  
  if(any(plot=="all") | any(plot=="avgLs")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Liabilities",ylim=c(range(liab[1:2,,t,],na.rm=T)),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),apply(liab[1,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),pch=16,cex=1.8) # Means: filled points
      points(jitter(1:N,amount=0.25),apply(liab[2,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),cex=1.8) # SDs: open points
    }    
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Mean","Sd"),bty="n",pch=c(16,1),col=alpha("Black",0.5),title="Within-pop L",pt.cex=1.8)
  }
  
  if(any(plot=="all") | any(plot=="avgLiabviol")){ # Plots violins of liabilities in final year for viewr random replicates
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,6,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Liabilities")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,6,t,r]))){
          vioplot(popstorage[inds,6,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
  } 
  
  if(any(plot=="all") | any(plot=="avgMeans")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Within-population means",ylim=c(range(liab[1:2,,t,])),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.2),apply(evol[1,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.2),apply(liab[1,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Breeding value","Liability"),pch=c(16,1),pt.cex=1.8,bty="n")
  }
  
  if(any(plot=="all") | any(plot=="avgSDs")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Within-population SDs",ylim=c(0,max(liab[2,,t,],na.rm=T)),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.2),apply(evol[2,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.2),apply(liab[2,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    if(legend) legend("bottomleft",legend=c("Breeding value","Liability"),pch=c(16,1),pt.cex=1.8,bty="n")
  }
  
  if(any(plot=="all") | any(plot=="avgDests")){
    plot(0:(N+1),0:(N+1),type="n",ylim=c(0,1),ylab="Frequency in subpopulation",xaxt="n",xlab="Subpopulation",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),apply(liab[3,,ts,r],FUN=mean,MARGIN=1,na.rm=T),pch=18,col=alpha("Black",0.5),cex=1.6) # proportion migrants
      for(j in 1:N){ # First plot all the points of 'going to A', then all 'going to B', etc.
        paltemp <- rep(pal[j],N)
        paltemp[j] <- NA
        points(jitter(1:N,amount=0.25),apply(evol[3+(j-1),,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(paltemp[1:N],0.7),cex=1.2)
      }
    }
    if(legend) {
      legend("topleft",legend=c("Migrants","Destination alleles"),pch=c(18,1),pt.cex=c(1.6,1.2),bty="n",bg=alpha("white",0.5),col=alpha("Black",0.5),cex=0.9)
      if(N>3) legend("topright",legend=sites[1:N],col=pal[1:N],bg=alpha("White",0.5),cex=1+(1-N)/20,pch=1,title="Destination")
    }
  }
  
  if(any(plot=="all") | any(plot=="avgDensity")){
    tt <- seq(dim(dyn)[1]-1999,dim(dyn)[1]-1,by=2) # Average pop densities over every year
    plot(0:(N+1),0:(N+1),type="n",xlab="Site",ylab="Local population size",xaxt="n",cex.lab=1.2,ylim=c(0,ifelse(grepl("extreme",saved,fixed=T),3,2)*K))
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),apply(dyn[tt,,r],FUN=mean,MARGIN=2,na.rm=T),pch=16,col=alpha(pal[1:N],0.5),cex=1.8) # Breeding 
      points(jitter(1:N,amount=0.25),apply(dyn[tt+1,,r],FUN=mean,MARGIN=2,na.rm=T),col=alpha(pal[1:N],0.5),cex=1.8) # Non-breeding
    }
    abline(h=K,lty=2)
    if(legend) legend("bottomleft",legend=c("Breeding","Non-breeding"),pch=c(16,1),pt.cex=1.8,col=alpha("Black",0.5),bty="n")
  }
  
  ### Final time step summary stats
  if(any(plot=="all") | any(plot=="EndBVviol")){
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,1,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Breeding values",bty="L")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,1,t,r]))){
          vioplot(popstorage[inds,1,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,na.rm=T,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
    mtext("Post shock",outer=T,cex=1.4,line=-0.5)
  }
  
  if(any(plot=="all") | any(plot=="EndBVbar")){
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,1,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab=expression(paste("Breeding values, ",italic(a))),bty="L")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      inds <- ((j-1)*K+1):(j*K)
      for(r in 1:reps){
        if(any(!is.na(popstorage[inds,1,t,r]))){
          vioplot(popstorage[inds,1,t,r],add=T,at=jitter(j,amount=0.25),col=alpha("white",0),colMed=pal[j],colMed2=alpha(pal[j],0.3),border=NA,na.rm=T,lineCol=NA,rectCol=alpha(pal[j],0.3),cex=1.5,pchMed=1) #pchMed=21 if filled points
        }
      }
    }
    abline(h=0,lty=2)
  } 
  
  if(any(plot=="all") | any(plot=="EndLs")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Liabilities",ylim=c(range(liab[1:2,,t,],na.rm=T)),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),liab[1,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8) # Means: filled points
      points(jitter(1:N,amount=0.25),liab[2,,t,r],col=alpha(pal[1:N],0.5),cex=1.8) # SDs: open points
    }    
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Mean","Sd"),bty="n",pch=c(16,1),col=alpha("Black",0.5),title="Within-pop L",pt.cex=1.8)
  }
  
  if(any(plot=="all") | any(plot=="EndLiabviol")){ # Plots violins of liabilities in final year for viewr random replicates
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,6,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Liabilities")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,6,t,r]))){
          vioplot(popstorage[inds,6,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
  } 
  
  if(any(plot=="all") | any(plot=="EndMeans")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Within-population means",ylim=c(range(liab[1:2,,t,])),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),evol[1,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.25),liab[1,,t,r],col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Breeding value","Liability"),pch=c(16,1),pt.cex=1.8,bty="n")
  }
  
  if(any(plot=="all") | any(plot=="EndSDs")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Within-population SDs",ylim=c(0,max(liab[2,,t,],na.rm=T)),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),evol[2,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.25),liab[2,,t,r],col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    if(legend) legend("bottomleft",legend=c("Breeding value","Liability"),pch=c(16,1),pt.cex=1.8,bty="n")
  }
  
  if(any(plot=="all") | any(plot=="EndDensity")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Site",ylab="Local population size",xaxt="n",cex.lab=1.2,ylim=c(0,ifelse(grepl("extreme",saved,fixed=T),3,2)*K))
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),dyn[(rec.int*t*2)-1,,r],pch=16,col=alpha(pal[1:N],0.5),cex=1.8)
      points(jitter(1:N,amount=0.25),dyn[rec.int*t*2,,r],col=alpha(pal[1:N],0.5),cex=1.8)
    }
    abline(h=K,lty=2)
    
    if(legend) legend("bottomleft",legend=c("Breeding","Non-breeding"),pch=c(16,1),pt.cex=1.8,col=alpha("Black",0.5),bty="n")
  }
  
  if(any(plot=="all") | any(plot=="EndDests")){
    plot(0:(N+1),0:(N+1),type="n",ylim=c(0,1),ylab="Frequency in subpopulation",xaxt="n",xlab="Subpopulation",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),liab[3,,t,r],pch=18,col=alpha("Black",0.5),cex=1.6) # proportion migrants
      for(j in 1:N){ # First plot all the points of 'going to A', then all 'going to B', etc.
        paltemp <- rep(pal[j],N)
        paltemp[j] <- NA
        points(jitter(1:N,amount=0.25),evol[3+(j-1),,t,r],col=alpha(paltemp[1:N],0.7),cex=1.2)
      }
    }
    if(legend) legend("topleft",legend=c("Migrants","Destination allele"),pch=c(18,1),pt.cex=c(1.6,1.2),bty="n",bg=alpha("white",0.5),col=alpha("Black",0.5),cex=0.9)
    #legend("topright",legend=sites[1:N],col=pal[1:N],bg=alpha("White",0.5),cex=1+(1-N)/20,pch=1,title="Destination")
  }
  
  ### End-but-pre-shock summary stats
  t <- max(1,dim(popstorage)[3]-2)
 
  if(any(plot=="all") | any(plot=="preshockBVviol")){
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,1,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Breeding values",bty="L")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,1,t,r]))){
          vioplot(popstorage[inds,1,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,na.rm=T,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
    mtext("Pre shock",outer=T,cex=1.4,line=-0.5)
  }
  
  if(any(plot=="all") | any(plot=="preshockBVbar")){
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,1,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Breeding values",bty="L")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      inds <- ((j-1)*K+1):(j*K)
      for(r in 1:reps){
        if(any(!is.na(popstorage[inds,1,t,r]))){
          vioplot(popstorage[inds,1,t,r],add=T,at=jitter(j,amount=0.25),col=alpha("white",0),border=NA,na.rm=T,lineCol=NA,rectCol=alpha(pal[j],0.5))
        }
      }
    }
    abline(h=0,lty=2)
  } 
  
  if(any(plot=="all") | any(plot=="preshockLs")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Liabilities",ylim=c(range(liab[1:2,,t,],na.rm=T)),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),liab[1,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8) # Means: filled points
      points(jitter(1:N,amount=0.25),liab[2,,t,r],col=alpha(pal[1:N],0.5),cex=1.8) # SDs: open points
    }    
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Mean","Sd"),bty="n",pch=c(16,1),col=alpha("Black",0.5),title="Within-pop L",pt.cex=1.8)
  }
  
  if(any(plot=="all") | any(plot=="preshockLviol")){ # Plots violins of liabilities in final year for three random replicates
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,6,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Liabilities")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,6,t,r]))){
          vioplot(popstorage[inds,6,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
  } 
  
  if(any(plot=="all") | any(plot=="preshockMeans")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Within-population means",ylim=c(range(liab[1:2,,t,])),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),evol[1,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.25),liab[1,,t,r],col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    abline(h=0,lty=2)
    legend("bottomleft",legend=c("Breeding value","Liability"),pch=c(16,1),pt.cex=1.8,bty="n")
  }
  
  if(any(plot=="all") | any(plot=="preshockSDs")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Within-population SDs",ylim=c(0,max(liab[2,,t,],na.rm=T)),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),evol[2,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.25),liab[2,,t,r],col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    legend("bottomleft",legend=c("Breeding value","Liability"),pch=c(16,1),pt.cex=1.8,bty="n")
  }

  if(any(plot=="all") | any(plot=="preshockDensity")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Site",ylab="Local population density",xaxt="n",cex.lab=1.2,ylim=c(0,2*K))
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),dyn[(rec.int*t*2)-1,,r],pch=16,col=alpha(pal[1:N],0.5),cex=1.8) # Breeding
      points(jitter(1:N,amount=0.25),dyn[rec.int*t*2,,r],col=alpha(pal[1:N],0.5),cex=1.8) # Non-breeding
    }
    abline(h=K,lty=2)
    if(legend) legend("bottomleft",legend=c("Breeding","Non-breeding"),pch=c(16,1),pt.cex=1.8,col=alpha("Black",0.5),bty="n")
  }
  
  if(any(plot=="all") | any(plot=="preshockDests")){
    plot(0:(N+1),0:(N+1),type="n",ylim=c(0,1),ylab="Frequency in subpopulation",xaxt="n",xlab="Subpopulation",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),liab[3,,t,r],pch=18,col=alpha("Black",0.5),cex=1.6) # proportion migrants
      for(j in 1:N){ # First plot all the points of 'going to A', then all 'going to B', etc.
        paltemp <- rep(pal[j],N)
        paltemp[j] <- NA
        points(jitter(1:N,amount=0.25),evol[3+(j-1),,t,r],col=alpha(paltemp[1:N],0.7),cex=1.2)
      }
    }
    if(legend){
      legend("topleft",legend=c("Migrants","Destination allele"),pch=c(18,1),pt.cex=c(1.6,1.2),bty="n",bg=alpha("white",0.5),col=alpha("Black",0.5),cex=0.9)
      if(N>3) legend("topright",legend=sites[1:N],col=pal[1:N],bg=alpha("White",0.5),cex=1+(1-N)/20,pch=1,title="Destination")
    }
  }

  # Age distribution from up to 3 random replicates before shocks.
  if(any(plot=="all") | any(plot=="agedist")){ 
    par(mfrow=c(length(viewr),1),mar=c(3.5,4.2,0.5,0.3),oma=c(0,0,1.5,0))
    for(r in 1:length(viewr)){
      plot(1:max.age,rep(0.1,max.age),type="n",xlab="",ylab="Frequency",cex.lab=1.2,cex.axis=1.2,ylim=c(0,max(evol[(N+3):(N+2+max.age),,,],na.rm=T)))
      for(j in 1:N){
        lines((1:max.age)+0.1*(j-1),evol[(N+3):(N+2+max.age),j,t,viewr[r]],col=pal[j],type="h",lwd=2)
        abline(v=mean(popstorage[((j-1)*K+1):(j*K),3,t,viewr[r]],na.rm=T),lty=2,col=pal[j])
      }
      text(max.age/2,0.25,paste("Replicate",viewr[r]),cex=1.5)
    }
    mtext("Age class",side=1,line=2)
    mtext("Age distribution, pre shock",side=3,line=-0.3,cex=1.2,outer=T)
    if(legend) legend("topright",title="Home site",legend=sites[1:N],col=pal,pch=rep(15,N))
  }
  
}

allbutpopsize <- c("indBVs", "indLs", "meanBVs", "meanLs","endpopsize","winteruse","agedist","Dest","meancorr","beforecorr","aftercorr",
                   "EndBVs", "EndBVviol", "EndLs", "EndLiabviol", "EndMeans", "EndSDs", "EndDests", "EndDensity","EndBVbar",
                   "avgBVs", "avgBVviol", "avgLs", "avgLiabviol", "avgMeans", "avgSDs", "avgDests", "avgDensity",
                   "preshockBVs", "preshockBVviol", "preshockLs", "preshockLviol", "preshockMeans", "preshockSDs", "preshockDests","preshockDensity","preshockBVbar")
shorter <- c("meanBVs","meanLs","endpopsize","preshockBVviol","preshockLviol","preshockDests","preshockDensity")

### Evolplots ###
## Function to plot detailed plots of last ~75 years
## Arguments:
#  local    : An alldb file (generated by create.alldb and polish.alldb) in the local environment.
#  saved    : An alldb file (generated by create.alldb and polish.alldb) inthat can be read from disk.
#  mainfile : A simulation result (generated by sim), created locally or read from file. Not needed much.
#  plot     : Shortnames for the different plots wanted: "all" (default),
#             "BVs","PEs","Ls","TEs","Dests"
#  pal      : Color palette to use.
evolplots <- function(local=NULL,saved=NULL,mainfile=NULL,plots="all",pal=viridis(3)){
  if(is.null(local)){
    tmp <- readRDS(file=saved)
  } else {(tmp <- local)}
  
  evoldyn <- tmp[[1]] # Contains means and sds ACROSS INDIVIDUALS of BV and PE per subpop, year and rep, and %migr per subpop, year and rep
  yearmeans <- tmp[[2]] # Contains means and sds ACROSS REPS of subpop-wide meanBV, sdBV, meanPE, sdPE and PropMigr
  dests <- tmp[[3]]
  years <- tmp[[4]]
  N <- max(evoldyn$Home,na.rm=T)
  R <- max(evoldyn$Rep,na.rm=T)
  
  sites <- switch(as.character(N),
                  "2"=c("North","South"),
                  "3"=c("North","Middle","South"),
                  "4"=c("North","Middle-N","Middle-S","South"),
                  "5"=c("North","Middle-N","Middle","Middle-S","South"),
                  "6"=c("North","","Middle-N","Middle-S","","South"),
                  "7"=c("North","","Middle-N","Middle","Middle-S","","South"),
                  "8"=c("North","","Middle-N","","","Middle-S","","South"),
                  "9"=c("North","","Middle-N","","Middle","","Middle-S","","South"))
  if(N>8){
    pal <- c(pal,brewerpal(8,"Set2"))
  }
  
  shockyear <- switch(is.null(mainfile), 
                      max(years,na.rm=T)-51,
                      mainfile[[3]]$shock.time) # Alternative below in case this doesn't work
  #if(is.null(mainfile)){
  #  shockyear <- max(years,na.rm=T)-51
  #} else {shockyear <- mainfile[[3]]$shock.time}
  
  par(mfrow=c(2,2),mar=c(4.2,4.2,1.3,0.4))
  
  if(any(plots=="all")|any(plots=="BVs")){
    plot(1:length(years),rep(0,length(years)),ylim=range(evoldyn$meanBV,na.rm=T),xaxt="n",cex.lab=1.2,
         ylab=expression(paste("Additive genetic effects, ",italic(a))),xlab="Year (+10000)",type="l",lty=2,lwd=2)
    for(n in 1:N){
      lines(1:length(years),filter(yearmeans,Home==n)$yearmeanBV,col=pal_v[n],lwd=2)
      lower <- filter(yearmeans,Home==n)$yearmeanBV - filter(yearmeans,Home==n)$yearsdBV*qt(0.975,R-1)/sqrt(R)
      upper <- filter(yearmeans,Home==n)$yearmeanBV + filter(yearmeans,Home==n)$yearsdBV*qt(0.975,R-1)/sqrt(R)
      polygon(c(1:length(years),length(years):1),c(lower,rev(upper)),col=alpha(pal_v[n],alpha=0.1),border=NA)
    }
    abline(v=25,lty=3,lwd=2)
    axis(1,at=c(1,(1:8)*25),labels=years[c(1,(1:8)*25)])
  }
  
  if(any(plots=="all")|any(plots=="Ls")){
    plot(1:length(years),rep(0,length(years)),ylim=range(evoldyn$meanL,na.rm=T),xaxt="n",ylab=expression(paste("Liabilities, ",italic(L))),
         xlab="Year (+10000)",type="l",lty=2,lwd=2,cex.lab=1.2)
    for(n in 1:N){
      lines(1:length(years),filter(yearmeans,Home==n)$yearmeanL,col=pal_v[n],lwd=2)
      lower <- filter(yearmeans,Home==n)$yearmeanL - filter(yearmeans,Home==n)$yearsdL*qt(0.975,R-1)/sqrt(R)
      upper <- filter(yearmeans,Home==n)$yearmeanL + filter(yearmeans,Home==n)$yearsdL*qt(0.975,R-1)/sqrt(R)
      polygon(c(1:length(years),length(years):1),c(lower,rev(upper)),col=alpha(pal_v[n],alpha=0.1),border=NA)
    }
    abline(v=25,lty=3,lwd=2)
    axis(1,at=c(1,(1:8)*25),labels=years[c(1,(1:8)*25)])
  }
  
  if(any(plots=="all")|any(plots=="PEs")){
    plot(1:length(years),rep(0,length(years)),ylim=c(min(-0.1,min(evoldyn$meanPE,na.rm=T)),max(0.1,max(evoldyn$meanPE,na.rm=T))),cex.lab=1.2,
         xaxt="n",ylab=expression(paste("Permanent effects, ",italic(B),"+",italic(b),"+",italic(C))),xlab="Year (+10000)",type="l",lty=2,lwd=2)
    for(n in 1:N){
      lines(1:length(years),filter(yearmeans,Home==n)$yearmeanPE,col=pal_v[n],lwd=2)
      lower <- filter(yearmeans,Home==n)$yearmeanPE - filter(yearmeans,Home==n)$yearsdPE*qt(0.975,R-1)/sqrt(R)
      upper <- filter(yearmeans,Home==n)$yearmeanPE + filter(yearmeans,Home==n)$yearsdPE*qt(0.975,R-1)/sqrt(R)
      polygon(c(1:length(years),length(years):1),c(lower,rev(upper)),col=alpha(pal_v[n],alpha=0.1),border=NA)
    }
    abline(v=25,lty=3,lwd=2)
    axis(1,at=c(1,(1:8)*25),labels=years[c(1,(1:8)*25)])
  }
  
  if(any(plots=="all")|any(plots=="TEs")){
    plot(1:length(years),rep(0,length(years)),ylim=c(min(-0.1,min(evoldyn$meanE,na.rm=T)),max(0.1,max(evoldyn$meanE,na.rm=T))),cex.lab=1.2,
         xaxt="n",ylab=expression(paste("Temporary effects, ",italic(D),"+",italic(E),"+",italic(e))),xlab="Year (+10000)",type="l",lty=2,lwd=2)
    for(n in 1:N){
      lines(1:length(years),filter(yearmeans,Home==n)$yearmeanE,col=pal_v[n],lwd=2)
      lower <- filter(yearmeans,Home==n)$yearmeanE - filter(yearmeans,Home==n)$yearsdE*qt(0.975,R-1)/sqrt(R)
      upper <- filter(yearmeans,Home==n)$yearmeanE + filter(yearmeans,Home==n)$yearsdE*qt(0.975,R-1)/sqrt(R)
      polygon(c(1:length(years),length(years):1),c(lower,rev(upper)),col=alpha(pal_v[n],alpha=0.1),border=NA)
    }
    abline(v=25,lty=3,lwd=2)
    axis(1,at=c(1,(1:8)*25),labels=years[c(1,(1:8)*25)])
    legend("bottomright",legend=sites,col=pal,bty="n",title="Subpopulation",lty=1,lwd=2,bg=alpha("white",0.3))
  }
  
  if(any(plots=="all")|any(plots=="Dests")){
    for(i in 1:N){ # For each home patch
      plot(1:length(years),rep(0,length(years)),ylim=c(0,1),xaxt="n",ylab="Frequency",xlab="Year (+10000)",type="n",main=paste(sites[i],"subpopulation"),cex.lab=1.2) # 1 plot per subpopulation for this one?
      for(n in 1:N){ # For each going-to patch
        tmp <- as.name(paste0("Dest",n))
        repsums <- filter(dests,Home==i) %>% group_by(Year) %>% summarise(repmeans:=mean({{tmp}}/tot),repsds:=sd({{tmp}}/tot)) # Take avg and SD for each year across reps
        lines(1:length(years),repsums$repmeans,lwd=2,col=pal_v[n])
        upper <- repsums$repmeans + repsums$repsds*qt(0.975,R-1)/sqrt(R)
        lower <- repsums$repmeans - repsums$repsds*qt(0.975,R-1)/sqrt(R)
        polygon(c(1:length(years),length(years):1),c(lower,rev(upper)),col=alpha(pal_v[n],alpha=0.1),border=NA)
        lines(1:length(years),rep(0,length(years)),lwd=2,col="White")
      }
      repstats <- filter(evoldyn,Home==i) %>% group_by(Year) %>% summarise(repmeans=mean(PropMigr),repsds=sd(PropMigr))
      upper <- repstats$repmeans + repstats$repsds*qt(0.975,R-1)/sqrt(R)
      lower <- repstats$repmeans - repstats$repsds*qt(0.975,R-1)/sqrt(R)
      polygon(c(1:length(years),length(years):1),c(lower,rev(upper)),col=alpha("Black",alpha=0.1),border=NA)
      lines(1:length(years),repstats$repmeans,lwd=2,col="Black") # Mean number of migrants
      abline(v=25,lty=3,lwd=2)
      #mtext(paste0("(",letters[i+1],")"),2,line=2.3,las=1,padj=-10.9,cex=1.1)
      axis(1,at=c(1,(1:8)*25),labels=years[c(1,(1:8)*25)])
    }
  }
  
}
#####

#### Functions for individual-level database & analyses ####

### create.alldb(): Function for creating individual-level database, with one row per individual per year. Slow for large files!
## Arguments:
# data     : A file created by sim() - locally or read from file
## Outputs: A data frame with columns: ID, Year, BV, Dest, Age, Home, Winter, L, Mort, Rep, PE, ind, max.age, problem (T/F), migrated (0/1)
create.alldb <- function(data){
  t <- Sys.time()
  colnames(data[[4]]) <- c("ID","Year","BV","Dest","Age","Home","Winter","L","Mort","PE")
  alldf <- data.frame("ID"=numeric(0),"Year"=numeric(0),"BV"=numeric(0),"Dest"=numeric(0),"Age"=numeric(0),"Home"=numeric(0),"Winter"=numeric(0),"L"=numeric(0),"Mort"=numeric(0),"Rep"=numeric(0),"PE"=numeric(0))
  for(i in 1:dim(data[[4]])[3]){
    tmp <- as.data.frame(data[[4]][,,i]) %>% mutate("Rep"=i)
    alldf <- bind_rows(alldf,tmp)
  }
  # Make sure each individual gets its unique identity "ind".
  # First rank the individuals by breeding value.
  alldf <- mutate(alldf,ind=rank(BV,ties.method="min",na.last="keep"))
  
  alldf <- mutate(alldf,new=as.numeric(paste0(ind,ID))) %>% filter(new!="NANA") %>% filter(!is.na(new)) %>% group_by(new)
  alldf <- group_by(alldf,new) %>% mutate(max.age=max(Age,na.rm=T), problem=any(tabulate(Age)>1) | sum(Age==0)>1) %>% ungroup() %>% arrange(new) 
  alldf <- alldf %>% dplyr::select(-ind) %>% rename(ind=new)

  print(Sys.time()-t)
  return(alldf)
}

### polish.alldb(): Function for fixing problems in a database created by create.alldb. Returns a database with same dimensions as create.alldb().
polish.alldb <- function(newdf){
  t <- Sys.time()
  allinds <- 1:10*length(newdf$ind) # or max(newdf$ind,na.rm=T)? Huge when new IDs merged ind and ID. 
  free <- allinds[-na.omit(unique(newdf$ind))] # Find free slots to fill in when we split our individuals
  free <- c(free,tail(allinds,1)+1)
  cycle <- unique(newdf$ind[newdf$problem==T])
  cycle <- cycle[!is.na(cycle)]
  while(length(cycle)*10>length(free)){
    free <- c(free,tail(free,1)+1:length(cycle))
  }
  
  counter <- 1 
  for(i in cycle){ # For each problematic individual
    start <- counter # This is where the counter was when we started splitting this individual.
    n <- which(newdf$ind==i) # First identify the places this individual occupies
    morts1 <- which(diff(as.matrix(newdf[n,5]))!=1) # Where is the age in the next cell not 1 more than the age in the previous cell?
    morts2 <- which(newdf[newdf$ind %in% i,9]>0) # Registered death events
    morts3 <- sort(union(morts1,morts2))
    morts1 <- which(newdf[n,2]==max(newdf$Year,na.rm=T)) # Individual who reach the end of the simulation
    morts <- sort(union(morts1,morts3))
    newdf[n[1]:n[morts[1]],12] <- free[counter]
    counter <- counter+1
    if(length(morts)>1){
      for(j in 2:length(morts)){
        newdf[(n[morts[j-1]]+1):(n[morts[j]]),12] <- free[counter]
        counter <- counter+1
      }
    }
    #    for(j in free[start:(counter-1)]){
    #     newdf[newdf$ind %in% j,12] <- max(newdf[newdf$ind %in% j,5],na.rm=T) #Now assign new max.ages to all individuals split in this process.
    #  }
    print(i)
  }
  newdf <- group_by(newdf,ind) %>% mutate(max.age=max(Age,na.rm=T), problem=any(tabulate(Age)>1) | sum(Age==0)>1 ) %>% ungroup() %>% arrange(ind)
  newdf <- mutate(newdf,migrated=as.numeric(L>0))
  print(Sys.time()-t)
  return(newdf)
}

## Function linkstrength(): For calculating the strength of the link between two patches.
# x is the square matrix of transition probabilities (columns: from summer location; rows: to winter location)
# s1 and s2 are the subpopulation for which link strength is to be calculated. Must be as.numeric(subpopulation name). 0 < s1,s2 =< N.
linkstrength <- function(x,patch1,patch2){
  if(nrow(x)!=ncol(x)){
    print("x not a square matrix")
    break
  }
  if(any(colSums(x)>1.001) | any(colSums(x)<0.999)){
    print("Transition probabilities do not sum to 1")
    break
  }
  N <- dim(x)[1]
  tmp <- x[patch1,patch2]*x[patch2,patch2] + x[patch2,patch1]*x[patch1,patch1]
  for(i in (1:N)[-c(patch1,patch2)]){
    tmp <- tmp + x[patch1,i]*x[patch2,i]
  }
  
  return(tmp)
}

## Function trans.matrix(): For calculating the transition probabilities between patches, using polished data frame alldf.
# Output: A NxN transition matrix from (columns) each patch to (rows) each patch.
trans.matrix <- function(alldf){
  N <- max(alldf$Home)
  tmp <- matrix(NA,N,N) # Columns: From summer location; Rows: To winter location
  for(j in 1:N){
    tmp[,j] <- tabulate(alldf$Winter[alldf$Home==j])
    tmp[,j] <- tmp[,j]/sum(tmp[,j]) # Standardize based on number of observations from each patch?
  }
  return(tmp)
}

## Function for producing network plot and calculating link strengths (calling linkstrength() function). Returns the vector of link strengths and (if plot=T) a plot where edge color and thickness reflects link strength.
network <- function(x,plot=T,widthscale=10,main=NULL){
  N <- dim(x)[1]
  
  i <- counter <-  1
  links <- numeric(choose(N,2))
  names <- character(choose(N,2))
  if(plot){
    pal <- brewer.pal(8,"Dark2")
    origin <- (N+1)/2
    xs <- 1:N
    ys <- sqrt(((N+1)/2)^2-((1:N)-origin)^2)
    plot(xs,ys,xlim=c(0.5,N+0.5),ylim=c(min(ys)-0.3,max(ys)+0.3),xaxt="n",yaxt="n",main=main,xlab="Subpopulation",ylab="")
    axis(1,at=xs,labels=LETTERS[xs])
  }  
  while(i<N){
    for(j in (i+1):N){
      links[counter] <- linkstrength(x,i,j)
      names[counter] <- paste0(LETTERS[i],LETTERS[j]) #Finish this
      if(plot){
        lines(c(xs[i],xs[j]),c(ys[i],ys[j]),lwd=log(widthscale*links[counter]),col=grey((1-links[counter])^2))
      }
      counter <- counter+1
    }
    i <- i+1
  }
  if(plot){
    points(xs,ys,cex=4,pch=21,col=pal[1:N],bg=pal[1:N])
  }
  names(links) <- names
  return(links)
}

### create.evoldyn: Function for gathering detailed data of evolutionary dynamics. Requires an alldb file (paste whole path and alldb.R file name).
## Other arguments:
# years: The years over which you extract data. Default 9926:10000 (tailored for an ECE at t=9950), but 26:100 is also useful.
# adult: Whether to only extract data for the adult part of the population (age >=age.mat) (T), or the whole population (F, default)
## Output:A list of
# [[1]]: evoldyn - data frame grouped by Home patch and Year. pop mean and sd BV, PE, L and E; plus PropMigr; for each rep
# [[2]]: yearmeans - data grame grouped by Home patch. mean and sd across reps of meanBV, PropMigr, meanPE, meanL; plus mean sdBV and mean sdL; for each year.
# [[3]]: dests - data frame grouped by Home patch and Year with destination alleles tabulated, for each rep.
# [[4]]: years - vector of the years used
# [[5]]: popmeans - data frame of seasonal pop sizes, mean and sd across reps in each year and season.
create.evoldyn <- function(alldb,years=9926:10000,adult=F){
  N <- max(alldb$Home,na.rm=T)
  R <- max(alldb$Rep,na.rm=T)
  alldb <- filter(alldb,Year %in% years)
  if(adult){
    alldb <- filter(alldb,Age>=3)
  }
  evoldyn <- alldb %>% group_by(Home,Year,Rep) %>% 
    summarise(meanBV=mean(BV),sdBV=sd(BV),PropMigr=mean(migrated,na.rm=T),meanPE=mean(PE),sdPE=sd(PE),meanL=mean(L,na.rm=T),sdL=sd(L,na.rm=T),
              meanE=mean(L-BV-PE,na.rm=T),sdE=sd(L-BV-PE,na.rm=T)) #ACROSS INDIVIDUALS
  yearmeans <- evoldyn %>% group_by(Home,Year) %>% 
    summarise(yearmeanBV=mean(meanBV),yearsdBV=sd(meanBV),yearmeansdBV=mean(sdBV),yearmeanPropMigr=mean(PropMigr),yearsdPropMigr=sd(PropMigr),
              yearmeanPE=mean(meanPE),yearsdPE=sd(meanPE),yearmeanL=mean(meanL),yearsdL=sd(meanL),yearmeansdL=mean(sdL),yearmeanE=mean(meanE),yearsdE=sd(meanE),yearmeansdE=mean(sdE)) # ACROSS REPS
  # Or would it be better to quantify temporary individual effects rather than liabs?
  dests <- alldb %>% group_by(Home, Year, Rep) # Destination alleles - tabulate
  dests <-  summarise(dests,Dest1=tabulate(Dest)[1],Dest2=tabulate(Dest)[2],Dest3=tabulate(Dest)[3],tot=sum(Dest1,Dest2,Dest3,na.rm=T)) # Edit manually if N!=3
  nb_locs <- alldb %>% filter(!is.na(L)) %>% # Gets rid of individuals without liabilities, i.e. died before migrating, so don't count towards nb season dens
    group_by(Year,Rep) %>% summarise(nb1=tabulate(Winter)[1],nb2=tabulate(Winter)[2],nb3=tabulate(Winter)[3],nbtot=sum(nb1,nb2,nb3,na.rm=T))
  b_locs <- alldb %>% filter(Age>0) %>% # Gets rid of juveniles, who don't count towards breeding season dens (but do to nb-season)
    group_by(Year,Rep) %>% summarise(b1=tabulate(Home)[1],b2=tabulate(Home)[2],b3=tabulate(Home)[3],btot=sum(b1,b2,b3,na.rm=T))
  popdyn <- full_join(ungroup(nb_locs),ungroup(b_locs))
  popmeans <- popdyn %>% group_by(Year) %>% 
    summarize(meannb1=mean(nb1),meannb2=mean(nb2),meannb3=mean(nb3),senb1=sd(nb1)/sqrt(max(Rep)),senb2=sd(nb2)/sqrt(max(Rep)),senb3=sd(nb3)/sqrt(max(Rep)),
              meanb1=mean(b1),meanb2=mean(b2),meanb3=mean(b3),seb1=sd(b1)/sqrt(max(Rep)),seb2=sd(b2)/sqrt(max(Rep)),seb3=sd(b3)/sqrt(max(Rep)))
  
  return(list(evoldyn,yearmeans,dests,years,popmeans))
}
#####

#### SIMULATION ####

## Arguments
#  N           Number of patches (zones)
#  K           Carrying capacity per patch
#  T           End time
#  mu          Per-locus mutation rate. Default 0.01.
#  m1          Breeding season survival. A vector of length N. Default rep(0.9,N), but can use c(0.95,0.9,0.85). But be careful with overall lifespan!
#  m2          Strength of winter density-dependence (higher m2 -> lower survival prob). A vector of length N. Default type is exponential
#              Exponential: Survival=dexp(n_i/K,m2)/m2) Use e.g. c(1.5,1,0.5). But lines cross! see ddfun(). 
#              Linear: Survival=1-(n_i/K)*m2[i]. Use e.g. c(1/2,1/3,1/10)
#  m3          Mortality cost of migration. Default 0.01
#  maxClutch   Maximum (female) clutch size. A vector of length N. Default rep(3,N).
#  age.mat     Age at first reproduction. Default 1.
#  disp        Dispersal rate (between-patch movement). Default 0.01.
#  sigmaPE_fine Size of permanent (micro)environmental effects on liability. Default 0. (Standard deviation of normal distr. from which values are chosen)
#  sigmaPE_coarse Size of permanent (macro)environmental effects on liability. Default 0. (Standard deviation of normal distr. from which values are chosen)
#  lambda_P_dens Size of permanent density effects (at birth) on liability. Default 0.
#  lambda_dens Size of temporary density effects on liability. Default 0.
#  sigmaE_fine Size of temporary (micro)environmental effects on liability. Default 0. (Standard deviation of normal distr. from which values are chosen)
#  sigmaE_coarse Size of temporary (macro)environmental effects on liability. Default 0. (Standard deviation of normal distr. from which values are chosen)
#  ddtype      Whether density dependence is "Exponential" (default) or "Linear".
#  sigma_m1    Stochastic components to breeding season survival. A vector of length N. Default rep(0,N).
#  sigma_m2    Stochastic components to overwinter survival. A vector of length N. Default rep(0,N).
#  rec.int     Time interval for recording full population matrix. Default 50.
#  rec.all.t   Number of final time steps over which to record full population matrix every year. Default 100.
#  rec.all.r   Number of replicates over which to record full population matrices. Default 5.
#  shock.time  When to introduce perturbations (shocks). Can be a vector. Default T-rec.all.t/2
#  shock.size  Size of shocks (fraction of population) for each patch. A vector of length N. Default c(0,0.8,0).
#  distance    Whether all patches are equally easy to reach (FALSE, default) or whether farther away ones are costlier (TRUE).
#  sex         Whether reproduction is sexual (T, default) or asexual (F).
#  v0          Initial additive genetic variance. Default 1.
#  a0          Initial mean breeding value. Default 0.
#  popsaved    Whether to initiate simulation with an already saved evolved population. Default FALSE. If TRUE, load the whole simulation output file, using readRDS(), or local file.
#  immat.return Do immature individuals (age<age.mat) return to breeding grounds (TRUE, default) or not (FALSE)?
#  age.liab    Age component to migration liability - effect = age*age.liab. Default 0.
#  freeze      If TRUE, freezes metapopulation at first element in shock.time and draws new pop each year from that. Default FALSE. Only implemented for sex=T.
#  ...         Additional arguments, such as size of shocks in further patches if using N>3.
## Outputs:
#  A list of 4 elements
#  [[1]] popstorage: Storage matrix of format:
#  Rows:    Individuals
#  Columns: Individual-level information recorded at regular intervals:
#  1: BV    Breeding value for genetic component of liability. Migration occurs if BV>0
#  2: Dest  Gene for determining destination if migrating. Randomly sampled among all patches except home patch
#  3: Age   Current age. Starts at 0
#  4: Home  Each individual's breeding population
#  5: Location: The zone the individual is spending the winter in
#  6: L     Liability value. Individual migrated if L>0.
#  7: PE    Permanent individual effects
#  for each recorded time step (3rd dim) and replicate (4th dim) 
#  [[2]] dyn: Storage matrix of dimensions:
#  1: Time steps (NB! two recordings per year! Odd numbers are summer, even numbers winter)
#  2: Patches
#  3: Replicate simulations.
#  [[3]] pars: list of all function inputs.
#  [[4]] alldata: Storage matrix of format:
#  Rows:    Individuals
#  Rolumns: Individual-level information recorded every year for the last rec.all.t years:
#  1: Individual ID
#  2: Year
#  3: BV gene
#  4: Dest gene
#  5: Age
#  6: Home patch
#  7: Wintering location in year t
#  8: Liability in year t
#  9: Cause of mortality. 0: Still alive. 1: Breeding season mortality. 2: Winter mortality. 3: Migration mortality.
#  10:Permanent individual effect on liability
#  for each replicate simulation (3rd dim)
sim <- function(reps=10,N=3,K=1000,T=5000,m1=rep(0.9,N),m2=c(1/2,1/3,1/5),m3=0.01,maxClutch=rep(3,N),age.mat=1,mu=0.01,m_size=0.1,
                sigmaPE_fine=0,sigmaPE_coarse=0,lambda_P_dens=0,lambda_dens=0,sigmaE_fine=0,sigmaE_coarse=0,ddtype="Exponential",
                disp=0.01,sigma_m1=rep(0,N),sigma_m2=rep(0,N),rec.int=50,rec.all.t=100,rec.all.r=5,shock.time=T-1-rec.all.t/2,
                shock.size=matrix(NA,N,length(shock.time)),info="None",distance=F,sex=TRUE,v0=1,a0=0,popsaved=FALSE,immat.return=TRUE,age.liab=0,freeze=F,...){
  t <- Sys.time()
  if(dim(shock.size)[1]!=N){
    print("Argument shock.size must be a matrix with N rows")
    break
  }
  if(dim(shock.size)[2]!=length(shock.time)){
    print("Wrong number of shocks. Require ncol(shock.size)==length(shock.time)")
  }
  if(length(m2)!=N){
    print(paste("Argument m2 must be a vector of length N"))
    break
  }
  if(length(sigma_m2)!=N){
    print(paste("Argument sigma_m2 must be a vector of length N"))
    break
  }
  if(length(m1)!=N){
    print(paste("Argument m1 must be a vector of length N"))
    break
  }
  if(length(sigma_m1)!=N){
    print(paste("Argument sigma_m1 must be a vector of length N"))
    break
  }
  if(length(maxClutch)!=N){
    print(paste("Argument maxClutch must be a vector of length N"))
    break
  }

  # Create storage vectors for recording data #
  dyn <- array(0,dim=c(T*2,N,reps)) # Dimensions: 1) Time steps (NB! two recordings per year! Odd numbers are summer, even numbers winter), 2) Patches, 3) Replicate simulations.
  colnames(dyn) <- LETTERS[1:N]
  popstorage <- array(NA,c(N*K,7,T%/%rec.int,reps)) # Dims: 1) Individuals, 2) Pop matrix contents (see below), 3) recording intervals, 4) replicates
  pars <- list(reps,N,K,T,m1,m2,m3,maxClutch,age.mat,mu,m_size,sigmaPE_fine,sigmaPE_coarse,lambda_P_dens,lambda_dens,sigmaE_fine,sigmaE_coarse,
               disp,sigma_m1,sigma_m2,rec.int,rec.all.t,rec.all.r,shock.time,shock.size,info,distance,sex,v0,a0,immat.return,age.liab,freeze)
  names(pars) <- c("reps","N","K","maxTime","Breeding-season survival","Winter mort","Migration cost","maxClutch","age@mat","Mut rate","Mut size",
                   "sigmaPE_fine","sigmaPE_coarse","lambda_P_dens","lambda_dens","sigmaE_fine","sigmaE_coarse","Dispersal rate",
                   "sigma_m1","sigma_m2","rec.int","rec.all.t","rec.all.r","shock.time","shock.size","info","distance","sex","v0","Initial mean BV","immat.return","age.liab","freeze")
  allstorage <- array(NA,c(N*K*rec.all.t*1.5,10,min(reps,rec.all.r)))

  for(r in 1:reps){
    #### Initialize population ####
    cue <- dens <- numeric(N*K) # Empty vectors for storing information on migration decision
    
    alldata <- matrix(NA,N*K*rec.all.t*1.5,10) 
    ### alldata matrix: 10 columns, 1 observation per row.
    # 1: Individual ID
    # 2: Year
    # 3: BV gene
    # 4: Dest gene
    # 5: Age
    # 6: Home patch
    # 7: Wintering location in year t
    # 8: Liability in year t
    # 9: Cause of mortality. 0: Still alive. 1: Breeding season mortality. 2: Winter mortality. 3: Migration mortality.
    # 10:Permanent individual effect on liability
    colnames(alldata) <- c("ID","Year","BV","Dest","Age","Home","Winter","L","Mort","PE")
    obs <- 1 # Counter for observations in Alldata
    
    n_off <- integer(1)
    inds <- fill <- rep(NA,K) # Initialize necessary storage vectors
    if(sex){
      mothers <- fathers <- rep(NA,K)
    } else {lucky <- rep(NA,K)}
    
    extinct <- N+1
    
    if(!isFALSE(popsaved)){ # If using a pre-saved population, 
      pop <- popsaved[[1]][,,dim(popsaved[[1]])[3],r]
    } else {
      pop <- matrix(NA,N*K,7)
      colnames(pop) <- c("BV","Dest","Age","Home","Winter","L","PE")
      ### pop matrix: 7 columns, 1 individual per row ###
      # BV:     Breeding value for genetic component of liability. Migration occurs if BV>0
      # Dest:   Gene for determining destination if migrating. Randomly sampled among all patches except home patch
      # Age:    Current age. Starts at 0
      # Home:   Each individual's breeding population
      # Location:The zone the individual is spending the winter in
      # L:      Liability value. Individual migrated if L>0.
      # PE:     Permanent individual component of liability.
      
      pop[,4] <- rep(1:N,each=K) # Assign zones to each slot in pop matrix
      pop[,5] <- pop[,4] # Overwinter at home unless migrating.
      pop[,3] <- rpois(N*K,3) # Assign starting age distribution. This can be made fancier!
      pop[,1] <- rnorm(N*K,a0,sqrt(v0)) # Assign random breeding values, normally distributed around a0 with var=v0.
      pop[,7] <- rep(rnorm(N,0,sigmaPE_coarse),each=K) # Permanent (coarse-grained/macro-) environmental individual effects
      pop[,7] <- rnorm(N*K,pop[,7],sigmaPE_fine) # Permanent (fine-grained/micro-) environmental individual effects
      pop[,7] <- pop[,7] + lambda_P_dens # Permanent individual effect of density at birth - initiated at K so it's 1. If initiating at N!=K, something like a patch-specific *sum(!is.na(pop[,1]))/K is needed
      
      # Assign random destination gene values - can't be the patch you breed in.
      for(i in 1:N){
        inds <- ((i-1)*K+1):(i*K) # IDs of individuals in this patch
        pop[inds,2] <- mysample((1:N)[-i],size=K,replace=T) # Assign random destination gene values 
      }
    }
    #####
    
    #### Start simulation ####
    
    for(i in 1:T){ # For every year
      pop[,3] <- pop[,3]+1 # Age increases for surviving individuals. Conveniently, NA+1=NA.
      
      if(i %% rec.int == 0){
        popstorage[,,i%/%rec.int,r] <- pop # Record individual traits
      }
      
      if(i==shock.time[1]){
        frozen <- pop
      }
      
      ### BREEDING SEASON ###
      
      for(j in (1:N)[-extinct]){ # For each patch
        fill[] <- NA # Reset parent vectors
        if(sex){
          mothers[] <- fathers[] <- NA
        } else {
          lucky[] <- NA
          removed <- 0}
        n_off <- 0  # Reset offspring counter
        
        inds[] <- ((j-1)*K+1):(j*K) # Slots in this patch
        alive <- which(!is.na(pop[inds,1])) + (j-1)*K # IDs of individuals alive 
        dyn[(i*2)-1,j,r] <- length(alive) # Pre-mortality and pre-breeding summer population census
        
        pmort <- rnorm(1,1-m1[j],sigma_m1[j])
        pmort <- min(1,pmort)
        pmort <- max(0,pmort)
        mort <- mysample(alive,size=round(pmort*length(alive)),replace=F) # Deterministic pre-breeding mortality due to seasonal suitability
        
        if(i>(T-rec.all.t) & rec.all.r >=r  & length(mort)>0){
          alldata[obs:(obs+length(mort)-1),1] <- as.numeric(paste0(r,mort)) # ID
          alldata[obs:(obs+length(mort)-1),2] <- i # Year
          alldata[obs:(obs+length(mort)-1),3] <- pop[mort,1] # BV
          alldata[obs:(obs+length(mort)-1),4] <- pop[mort,2] # Dest
          alldata[obs:(obs+length(mort)-1),5] <- pop[mort,3] # Age
          alldata[obs:(obs+length(mort)-1),6] <- j # Home
          alldata[obs:(obs+length(mort)-1),9] <- 1 # Cause of death - breeding season
          alldata[obs:(obs+length(mort)-1),10] <- pop[mort,7] # PE
          obs <- obs+length(mort)
        }
        
        pop[mort,] <- NA 
        
        ## Reproduction ##
        
        n_off <- sum(is.na(pop[inds,1])) # How many empty slots are there for offspring
        if(n_off==K){ # If everyone is dead
          print(paste("Subpopulation",j,"extinct in summer time",i))
          pop[inds,] <- NA # Shouldn't be necessary, but try 15.11 whether this is what causes ages to grow continue growing?
          extinct <- c(j,extinct)
          break
        }
        
        dens[inds] <- sum(!is.na(pop[inds,1]))/K
        
        if(max(pop[inds,3],na.rm=T)>=age.mat & n_off>0){ #Check that there are breeders in the population
          n_off <- min(round(sum(is.mature(pop[inds,3],a=age.mat))*maxClutch[j]),n_off)
          if(n_off==0){
            print(paste("No offspring produced in patch",j,"in time",i))
            print("hi")
            break
          }
          if(freeze==F | i<shock.time[1]){ # Unless freezing
            if(sex){ # Sexual reproduction
              mothers[1:n_off] <- mysample(which(is.mature(pop[inds,3],a=age.mat))+(j-1)*K,size=n_off,replace=T) # Sample mothers
              fathers[1:n_off] <- mysample(which(is.mature(pop[inds,3],a=age.mat))+(j-1)*K,size=n_off,replace=T) # Sample fathers
              
              # Which slots get new offspring?
              fill[1:n_off] <- head(which(is.na(pop[inds,1])),n_off)+(j-1)*K # Empty slots in this patch
              pop[fill[1:n_off],3] <- 0 # Assign ages to all newborns
              pop[fill[1:n_off],4:5] <- j # Assign home and current patch to all newborns
              pop[fill[1:n_off],7] <- rnorm(n_off,rnorm(1,0,sigmaPE_coarse),sigmaPE_fine) # Permanent (developmental plasticity) effect, fine-grained var around coarse-grained mean.
              pop[fill[1:n_off],7] <- pop[fill[1:n_off],7] + lambda_P_dens*dens[j*K] # Permanent effect of adult pop density at birth.
              pop[fill[1:n_off],1] <- apply(cbind(pop[mothers[1:n_off],1],pop[fathers[1:n_off],1]),FUN=mean,MARGIN=1) # Midpoint of parents' BVs
              # Next, newborn's BV is drawn from a Gaussian around midparent BV with sd given by initial genetic variance v0
              pop[fill[1:n_off],1] <- pop[fill[1:n_off],1]+rnorm(n_off,0,sqrt(v0))
              pop[fill[1:n_off],2] <- ifelse(runif(n_off)>0.5,pop[mothers[1:n_off],2],pop[fathers[1:n_off],2]) # Mendelian segregation at the Destination locus
            
              if(any(is.na(pop[fill[1:n_off],1:5]))){
                print("oops")
                print(paste0("N_off=",n_off,", last parent ID=",tail(parents,1)))
              }
            } else { # Asexual reproduction
              if(maxClutch[j]==1){ # Simplify if max clutch size= 1.
                parents <- which(is.mature(pop[inds,3],a=age.mat))+(j-1)*K
                n_off <- min(n_off,parents)
                lucky[1:n_off] <- mysample(parents,size=n_off,replace=F)
              } else {
                lucky[1:n_off] <- mysample(which(is.mature(pop[inds,3],a=age.mat)) + (j-1)*K,size=n_off,replace=T) # IDs of parents who get offspring
                if(any(tabulate(lucky)>maxClutch[j])){ # Make sure no parent gets more than 3 offspring
                  problem <- which(tabulate(lucky)>maxClutch[j]) # IDs of parents who got more than 3 offspring
                  for(p in problem){ # For each problematic parent
                    while(tabulate(lucky)[p]>maxClutch[j]){ # As long as they have more than 3 offspring
                      lucky[1:(n_off-removed-1)] <- lucky[(1:(n_off-removed))[-head(which(lucky==p),1)]] # Remove one of their offspring
                      lucky[(n_off-removed):K] <- NA
                      removed <- removed+1 
                    }
                  }
                }
                n_off <- n_off-removed
              }
              
              # Which slots get new offspring?
              fill[1:n_off] <- head(which(is.na(pop[inds,1])),n_off)+(j-1)*K
              pop[fill[1:n_off],3] <- 0 # Assign ages to all newborns
              pop[fill[1:n_off],4:5] <- j # Assign home and current patch to all newborns
              pop[fill[1:n_off],7] <- rnorm(n_off,rnorm(1,0,sigmaPE_coarse),sigmaPE_fine) # Permanent (developmental plasticity) effect, fine-grained var around coarse-grained mean.
              pop[fill[1:n_off],7] <- pop[fill[1:n_off],7] + lambda_P_dens*dens[j*K] # Permanent effect of adult pop density at birth.
              pop[fill[1:n_off],1] <- pop[lucky[1:n_off],1] # Newborns inherit their parent's genes
              pop[fill[1:n_off],2] <- pop[lucky[1:n_off],2] # Newborns inherit their parent's genes
            }
            
            # Mutation - only for Destination if sexual, BV is taken care of by the standard deviation of the random draw.
            if(N>2){
              mutants <- sample(fill[1:n_off],round(n_off*mu)) # A subset of newborns mutate their Destination
              for(m in mutants){
                pop[m,2] <- mysample((1:N)[-c(j,pop[m,2])],1)
              }
            }
            if(!sex){
              mutants <- sample(fill[1:n_off],round(n_off*mu)) # A subset of newborns mutate their BV. THOUGHT: Since many genes of small effect - mutate more often?
              pop[mutants,1] <- rnorm(length(mutants),pop[mutants,1],m_size)
            }
          } else {
            # Sampling from frozen population
            mothers[1:n_off] <- mysample(which(is.mature(frozen[inds,3],a=age.mat))+(j-1)*K,size=n_off,replace=T) # Sample mothers
            fathers[1:n_off] <- mysample(which(is.mature(frozen[inds,3],a=age.mat))+(j-1)*K,size=n_off,replace=T) # Sample fathers
            
            # Which slots get new offspring?
            fill[1:n_off] <- head(which(is.na(pop[inds,1])),n_off)+(j-1)*K # Empty slots in this patch
            pop[fill[1:n_off],3] <- 0 # Assign ages to all newborns
            pop[fill[1:n_off],4:5] <- j # Assign home and current patch to all newborns
            pop[fill[1:n_off],7] <- rnorm(n_off,rnorm(1,0,sigmaPE_coarse),sigmaPE_fine) # Permanent (developmental plasticity) effect, fine-grained var around coarse-grained mean.
            pop[fill[1:n_off],7] <- pop[fill[1:n_off],7] + lambda_P_dens*dens[j*K] # Permanent effect of adult pop density at birth.
            pop[fill[1:n_off],1] <- apply(cbind(frozen[mothers[1:n_off],1],frozen[fathers[1:n_off],1]),FUN=mean,MARGIN=1) # Midpoint of parents' BVs
            # Next, newborn's BV is drawn from a Gaussian around midparent BV with sd given by distance among parents' BVs
            # pop[fill[1:n_off],1] <- pop[fill[1:n_off],1]+rnorm(n_off,0,v0) # Removed due to no-evol
            pop[fill[1:n_off],2] <- ifelse(runif(n_off)>0.5,frozen[mothers[1:n_off],2],frozen[fathers[1:n_off],2]) # Mendelian segregation at the Destination locus
            
      #      if(any(is.na(pop[fill[1:n_off],1:5]))){
      #        print("oops")
      #        print(paste0("N_off=",n_off," in subpop ",j))
      #      }
            
          }
        } else {
          print(paste("No breeding in patch",j,"in time",i))
        }
      }
      
      if(any(!is.na(pop))==FALSE){
        print(paste("Metapopulation extinct in summer time",i))
        break
      }
      if(length(extinct)==N+1){
        print(paste("Metapopulation extinct in summer time",i))
      }
      ### AUTUMN ###
      
      # Choose stochastic environmental cues
      cue <- rnorm(N*K,rep(rnorm(N,0,sigmaE_coarse),each=K),sigmaE_fine) # Temporary environmental effect, fine-grained around coarse-grained var
      
      # Work out density cue for each patch - hashtagging away because now we only use adults present??
      #for(j in 1:N){
      #  inds <- ((j-1)*K+1):(j*K) # Slots in this patch
      #  dens[((j-1)*K+1):(j*K)] <- sum(!is.na(pop[inds,1]))/K # Higher cue if fewer dead
      #}  
      #dens <- dens*lambda_dens
      
      # Migrate if env val + pop dens cause liability to exceed threshold
      pop[,6] <- pop[,1] + pop[,7] + cue + dens*lambda_dens + pop[,3]*age.liab # LIABILITY EQUATION 
      pop[which(pop[,6]>0),5] <- pop[which(pop[,6]>0),2] # Individuals who migrate go to patch determined by their Dest gene.
      cue[which(pop[,6]>0)] <- 0 # Individuals who migrate are 'released' from the bad weather in their home patch and get a neutral cue.
      
      if(i %% rec.int == 0){
        popstorage[,5:6,i%/%rec.int,r] <- pop[,5:6] # If recording pop matrix, track liabilities and where individuals overwintered.
      }
      
      # Deterministic mortality depending on seasonal suitability and density in post-migration patch
      for(j in 1:N){
        whos_there <- which(pop[,5]==j)
        localN <- length(whos_there)
        dyn[i*2,j,r] <- localN # Pre-mortality population census
        
        if(m2[j]==0){ # Calculate density-dependence term
          ddterm <- 0
        } else {
          ddterm <- switch(ddtype,
                           "Exponential"=1-dexp(localN/K,m2[j])/m2[j],
                           "Linear"=localN/K*m2[j]) 
        # ddterm <- 1-dexp(localN/K,m2[j])/m2[j] # Exponential
        # ddterm <- localN/K*m2[j] # Linear
        }
        pmort <- rnorm(1,ddterm,sigma_m2[j]) # Proportion killed
        pmort <- min(1,pmort)
        pmort <- max(0,pmort)
        if(i %in% shock.time){
          pmort <- max(c(pmort,shock.size[j,which(shock.time==i)]),na.rm=T)
        }

        # Mortality correlated with cue received? Hashtag away the one you're not using
        if(info=="m2"){
          mort <- mysample(whos_there,size=localN*pmort,replace=F,prob=exp(cue[whos_there])) # Cue gives information: Those who didn't migrate get worse weather.
        } else {
          mort <- mysample(whos_there,size=localN*pmort,replace=F) # Cue gives no information
        }
        #mort <- mysample(whos_there,size=min(K,round(localN*localN/K*m2[j])),replace=F)
        
        if(i>(T-rec.all.t) & rec.all.r >=r & length(mort)>0){
          alldata[obs:(obs+length(mort)-1),1] <- as.numeric(paste0(r,mort)) # ID
          alldata[obs:(obs+length(mort)-1),2] <- i  # Year
          alldata[obs:(obs+length(mort)-1),3] <- pop[mort,1] # BV
          alldata[obs:(obs+length(mort)-1),4] <- pop[mort,2] # Dest
          alldata[obs:(obs+length(mort)-1),5] <- pop[mort,3] # Age
          alldata[obs:(obs+length(mort)-1),6] <- pop[mort,4] # Home
          alldata[obs:(obs+length(mort)-1),7] <- j # Current patch
          alldata[obs:(obs+length(mort)-1),8] <- pop[mort,6] # Liab
          alldata[obs:(obs+length(mort)-1),9] <- 2 # Cause of death - winter
          alldata[obs:(obs+length(mort)-1),10] <- pop[mort,7] # PE
          obs <- obs+length(mort)
        }
        
        pop[mort,] <- NA
      }
      ### WINTER ###
      
      if(any(!is.na(pop))==FALSE){
        print(paste("Metapopulation extinct in winter time",i))
        break
      }
      # Choose a stochastic env value for each site
       
      # Individuals still in breeding site migrate if env val + pop dens cause liability to exceed threshold
      
      # Deterministic mortality depending on seasonal suitability in post-migration patch
      
      ### End of winter
      
      # Dispersal: Some individuals from each patch "disperse" to a different patch.
      for(j in 1:N){
        inds <- ((j-1)*K+1):(j*K) # Slots in this patch
        free <- which(is.na(pop[inds,1])) + (j-1)*K # Find free slots
        maxdisp <- length(which(is.na(pop[,1]))[-inds])*disp/(N-1) # Max number of dispersers
        n_imm <- min(length(free),round(maxdisp)) # How many immigrants? Not more than free slots, nor as determined by dispersal parameter.
        if(n_imm>0){
          immigrants <- sample((1:(N*K))[-c(inds,which(is.na(pop[,1])))],size=n_imm,replace=F) # Choose random immigrants from other patches
          pop[free[1:n_imm],] <- pop[immigrants,] # Move their slots to this patch
          pop[free[1:n_imm],4] <- j # Immigrants get new home patch
          pop[free[1:n_imm],2] <- mysample((1:N)[-j],size=n_imm,replace=T) # Immigrants get new Dest gene
          pop[immigrants, ] <- NA # Remove their slots in other patches
        }
        if(j %in% extinct){
          extinct <- extinct[-which(extinct==j)]
        }
      }
      # Cost of migration: A fraction of everyone who migrated still alive are killed.
      if(distance){ # Are adjacent patches easier to reach than farther away ones?
        if(any(pop[,6]>0,na.rm=T)){
          mort <- mysample(which(pop[,6]>0),size=round(m3*sum(pop[,6]>0,na.rm=T)),replace=F,prob=abs(pop[which(pop[,6]>0),4]-pop[which(pop[,6]>0),2]))
        }
      } else{ # All patches are equally easy to reach
        mort <- mysample(which(pop[,6]>0),size=round(m3*sum(pop[,6]>0,na.rm=T)),replace=F)
      }
      
      # If recording individual-level data, track those who died from migration
      if(i>(T-rec.all.t) & rec.all.r >=r  & length(mort)>0){
        alldata[obs:(obs+length(mort)-1),1] <- as.numeric(paste0(r,mort)) # ID
        alldata[obs:(obs+length(mort)-1),2] <- i # Year
        alldata[obs:(obs+length(mort)-1),3] <- pop[mort,1] # BV
        alldata[obs:(obs+length(mort)-1),4] <- pop[mort,2] # Dest
        alldata[obs:(obs+length(mort)-1),5] <- pop[mort,3] # Age
        alldata[obs:(obs+length(mort)-1),6] <- pop[mort,4] # Home
        alldata[obs:(obs+length(mort)-1),7] <- pop[mort,5] # Current (winter) patch
        alldata[obs:(obs+length(mort)-1),8] <- pop[mort,6] # Liab
        alldata[obs:(obs+length(mort)-1),9] <- 3 # Cause of death - migration
        alldata[obs:(obs+length(mort)-1),10] <- pop[mort,7] # PE
        obs <- obs+length(mort)
      }
      
      # Kill individuals who died from migration
      pop[mort,] <- NA
      
      # If recording individual-level data, update all information on survived individuals.
      if(i>(T-rec.all.t) & rec.all.r >= r & length(alive)>0){
        alive <- which(!is.na(pop[,1])) 
        alldata[obs:(obs+length(alive)-1),1] <- as.numeric(paste0(r,alive)) # ID
        alldata[obs:(obs+length(alive)-1),2] <- i # Year
        alldata[obs:(obs+length(alive)-1),3] <- pop[alive,1] # BV
        alldata[obs:(obs+length(alive)-1),4] <- pop[alive,2] # Dest
        alldata[obs:(obs+length(alive)-1),5] <- pop[alive,3] # Age
        alldata[obs:(obs+length(alive)-1),6] <- pop[alive,4] # Home
        alldata[obs:(obs+length(alive)-1),7] <- pop[alive,5] # Current (winter) patch
        alldata[obs:(obs+length(alive)-1),8] <- pop[alive,6] # Liab
        alldata[obs:(obs+length(alive)-1),9] <- 0 # Cause of death - still alive
        alldata[obs:(obs+length(alive)-1),10] <- pop[alive,7] # PE
        obs <- obs+length(alive)
      }
      # All survived individuals return to their home patch # MAKE MATURE ONLY??
      if(immat.return){
        pop[,5] <- pop[,4]
      } else {
        pop[which(pop[,3]>=age.mat),5] <- pop[which(pop[,3]>=age.mat),4]
      }
      #print(i)
    }

    # End simulation
    #####
    
    if(r <= rec.all.r) {
      allstorage[,,r] <- alldata
    }
    print(r)
  }
  print(Sys.time()-t)
  return(list(popstorage,dyn,pars,allstorage))
}
#####

#### Running and saving your first simulation ####

moderate <- c(1/2,1/3,1/5) # A weak environmental gradient
extreme <- c(5,1/2,1/10) # A strong environmental gradient

# Run a full set of evolutionary simulations(takes 4-5 minutes)
test <- sim(T=5000,reps=5,m2=moderate,m3=0.02,disp=0.01,v0=0.5,maxClutch=rep(1,3)) # Simulation. See L1138.
statstmp <- popstats(test) # Create stats
popplots(statstmp,mainfile=tmp,plot="avgBVviol") # or whatever plots you want

# If you want to save
saveRDS(tmp,file="Baseline.R")
saveRDS(statstmp,file="Baseline stats.R")

## Read from file
popplots(saved="Baseline stats.R",plot="avgBVviol") # plotting (uses stats file)
tmp <- readRDS("Baseline.R") # (uses mainfile)
#####


#### Running climate change scenarios (takes <20 secs) ####
# Uses popsaved argument (read local or saved mainfile)

#1: constant poor conditions in south patch.
CC1 <- sim(T=200,N=3,m2=moderate,popsaved=readRDS("Baseline.R"),K=1000,reps=5,v0=0.5,maxClutch=rep(1,3),
            shock.time=50:200,shock.size=matrix(c(rep(0,151*2),rep(0.5,151)),3,151,byrow=T),m3=0.02,rec.all.t=200,rec.all.r=20)
statstmp1 <- popstats(CC1)

#2: gradual worsening in south patch.
CC2 <- sim(T=200,N=3,m2=moderate,popsaved=readRDS("Baseline.R"),K=1000,reps=5,rec.all.t=200,rec.all.r=20,v0=0.5,maxClutch=rep(1,3),
            shock.time=50:200,shock.size=matrix(c(rep(0,151*2),seq(0.3,0.8,by=0.01),rep(0.8,100)),3,151,byrow=T),m3=0.02)
statstmp2 <- popstats(CC2)

#3: regularly staggered shocks in same patch.
CC3 <- sim(T=100,N=3,m2=moderate,popsaved=readRDS("D:\\Migration\\paper2\\m3=0.02\\v0=0.5, f=1\\Control.R"),K=1000,reps=20,rec.all.t=100,rec.all.r=20,
            shock.time=seq(50,90,by=10),shock.size=matrix(c(rep(0,10),rep(0.5,5)),3,5,byrow=T),m3=0.02,v0=0.5,maxClutch=rep(1,3))
statstmp3 <- popstats(CC3)

#4: regularly staggered shocks in different patches.
CC4 <- sim(T=100,N=3,m2=moderate,popsaved=readRDS("Baseline.R"),K=1000,reps=5,rec.all.t=100,rec.all.r=5,v0=0.5,m3=0.02,maxClutch=rep(1,3),
            shock.time=c(50,60,70,80,90),shock.size=matrix(c(c(0,0,0,0.8,0),c(0,0.8,0,0,0.8),c(0.8,0,0.8,0,0)),3,5,byrow=T))
statstmp4 <- popstats(CC4)

#5: irregularly staggered shocks in same patch.
t<-50
ts<-t
counter <- 1
set.seed(1)
tmp <- rgeom(15,0.15) # Assume time between shocks is geometrically distributed
while(t<100){
  t <- t+tmp[counter]+1
  ts <- c(ts,t)
  if(t<100) counter <- counter+1
}

CC5 <- sim(T=100,N=3,m2=moderate,popsaved=baseline,K=1000,reps=5,rec.all.t=100,rec.all.r=5,v0=0.5,m3=0.02,maxClutch=rep(1,3),
            shock.time=ts,shock.size=matrix(c(rep(0,2*length(ts)),rep(0.8,length(ts))),3,length(ts),byrow=T))
statstmp5 <- popstats(CC5)
#####

#### Creating individual-level database (takes ~1 minute) ####
alldb <- create.alldb(test)
alldb <- polish.alldb(alldb)
head(alldb,n=20)
#saveRDS(alldb,file="Baseline alldb.R")

# Extract evolutionary dynamics from alldb - see L1099
evoltmp <- create.evoldyn(alldb,years=4901:5000,adult=F) # Remember to change years if needed!
head(evoltmp[[1]],n=20)
#saveRDS(alldb,file="Baseline evoltmp.R")
##### 

#### Figures of all liability components
pdf(file="Columnfig.R",width=3.5,height=10)
par(mar=c(2,4.2,0.5,0.3),mfrow=c(4,1),oma=c(2,0,0.5,0))
plot(1:175,rep(0,175),ylim=range(c(df$yearmeanBV,df2$yearmeanBV)),ylab="Breeding values",xlab="Year (+10000)",xaxt="n",cex.lab=1.2,type="n")
axis(1,at=c(1,25*(1:7)),labels=25*(1:8))
#plot(1:75,rep(0,75),ylim=range(c(df2$yearmeanBV,df$yearmeanBV)),ylab="Breeding values",xlab="Year (+10000)",xaxt="n",cex.lab=1.2,type="n")
#axis(1,at=c(1,25*(1:3)),labels=25*(1:4))
for(j in 1:3){
  lines(1:175,filter(df,Home==j)$yearmeanBV,col=pal_v[j],lwd=2)
  lines(1:175,filter(df2,Home==j)$yearmeanBV,col=pal_v[j],lwd=2,lty=3)
  upper1 <- filter(df,Home==j)$yearmeanBV+filter(df,Home==j)$yearsdBV/sqrt(20)
  upper2 <- filter(df2,Home==j)$yearmeanBV+filter(df2,Home==j)$yearsdBV/sqrt(20)
  lower1 <- filter(df,Home==j)$yearmeanBV-filter(df,Home==j)$yearsdBV/sqrt(20)
  lower2 <- filter(df2,Home==j)$yearmeanBV-filter(df2,Home==j)$yearsdBV/sqrt(20)
  polygon(c(1:175,175:1),c(lower1,rev(upper1)),col=alpha(pal_v[j],alpha=0.2),border=NA)
  polygon(c(1:175,175:1),c(lower2,rev(upper2)),col=alpha(pal_v[j],alpha=0.2),border=NA)
}
mtext("B",side=3,at=1,line=-2)
abline(v=25,lty=3,lwd=1)
abline(h=0,lty=2)

plot(1:175,rep(0,175),ylim=range(c(df2$yearmeanL,df$yearmeanL)),ylab="Liabilities",xlab="Year(+10000)",xaxt="n",cex.lab=1.2,type="n")
axis(1,at=c(1,25*(1:7)),labels=25*(1:8))
#plot(1:75,rep(0,75),ylim=range(c(df2$yearmeanL,df$yearmeanL)),ylab="Liabilities",xlab="Year(+10000)",xaxt="n",cex.lab=1.2,type="n")
#axis(1,at=c(1,25*(1:3)),labels=25*(1:4))
for(j in 1:3){
  lines(1:175,filter(df,Home==j)$yearmeanL,col=pal_v[j],lwd=2)
  lines(1:175,filter(df2,Home==j)$yearmeanL,col=pal_v[j],lwd=2,lty=3)
  upper1 <- filter(df,Home==j)$yearmeanL+filter(df,Home==j)$yearsdL/sqrt(20)
  upper2 <- filter(df2,Home==j)$yearmeanL+filter(df2,Home==j)$yearsdL/sqrt(20)
  lower1 <- filter(df,Home==j)$yearmeanL-filter(df,Home==j)$yearsdL/sqrt(20)
  lower2 <- filter(df2,Home==j)$yearmeanL-filter(df2,Home==j)$yearsdL/sqrt(20)
  polygon(c(1:175,175:1),c(lower1,rev(upper1)),col=alpha(pal_v[j],alpha=0.2),border=NA)
  polygon(c(1:175,175:1),c(lower2,rev(upper2)),col=alpha(pal_v[j],alpha=0.2),border=NA)
}
mtext("E",side=3,at=1,line=-2)
abline(v=25,lty=3,lwd=1)
abline(h=0,lty=2)

plot(1:175,rep(0,175),ylim=c(min(c(0.4,df2$yearmeanPE,df$yearmeanPE)),max(c(0.02,df$yearmeanPE,df2$yearmeanPE))),ylab="Permanent effects",xlab="",xaxt="n",cex.lab=1.2,type="n")
axis(1,at=c(1,25*(1:7)),labels=25*(1:8))
#plot(1:75,rep(0,75),ylim=c(min(c(0.5,df$yearmeanPE,df2$yearmeanPE)),max(c(0.02,df$yearmeanPE,df2$yearmeanPE))),ylab="Permanent effects",xlab="",xaxt="n",cex.lab=1.2,type="n")
#axis(1,at=c(1,25*(1:3)),labels=25*(1:4))
for(j in 1:3){
  lines(1:175,filter(df,Home==j)$yearmeanPE,col=pal_v[j],lwd=2)
  lines(1:175,filter(df2,Home==j)$yearmeanPE,col=pal_v[j],lwd=2,lty=3)
  upper1 <- filter(df,Home==j)$yearmeanPE+filter(df,Home==j)$yearsdPE/sqrt(20)
  upper2 <- filter(df2,Home==j)$yearmeanPE+filter(df2,Home==j)$yearsdPE/sqrt(20)
  lower1 <- filter(df,Home==j)$yearmeanPE-filter(df,Home==j)$yearsdPE/sqrt(20)
  lower2 <- filter(df2,Home==j)$yearmeanPE-filter(df2,Home==j)$yearsdPE/sqrt(20)
  polygon(c(1:175,175:1),c(lower1,rev(upper1)),col=alpha(pal_v[j],alpha=0.2),border=NA)
  polygon(c(1:175,175:1),c(lower2,rev(upper2)),col=alpha(pal_v[j],alpha=0.2),border=NA)
}
mtext("H",side=3,at=1,line=-2)
abline(v=25,lty=3,lwd=1)
abline(h=0,lty=2)

plot(1:175,rep(0,175),ylim=c(min(c(0.4,df2$yearmeanE,df$yearmeanE)),max(c(0.02,df$yearmeanE,df2$yearmeanE))),ylab="Temporary effects",xlab="Year(+10000)",xaxt="n",cex.lab=1.2,type="n")
axis(1,at=c(1,25*(1:7)),labels=25*(1:8))
#plot(1:75,rep(0,75),ylim=c(min(c(-0.02,df2$yearmeanE,df$yearmeanE)),max(c(0.02,df$yearmeanE,df2$yearmeanE))),ylab="Temporary effects",xlab="",xaxt="n",cex.lab=1.2,type="n")
#axis(1,at=c(1,25*(1:3)),labels=25*(1:4))
for(j in 1:3){
  lines(1:175,filter(df,Home==j)$yearmeanE,col=pal_v[j],lwd=2)
  lines(1:175,filter(df2,Home==j)$yearmeanE,col=pal_v[j],lwd=2,lty=3)
  upper1 <- filter(df,Home==j)$yearmeanE+filter(df,Home==j)$yearsdE/sqrt(20)
  upper2 <- filter(df2,Home==j)$yearmeanE+filter(df2,Home==j)$yearsdE/sqrt(20)
  lower1 <- filter(df,Home==j)$yearmeanE-filter(df,Home==j)$yearsdE/sqrt(20)
  lower2 <- filter(df2,Home==j)$yearmeanE-filter(df2,Home==j)$yearsdE/sqrt(20)
  polygon(c(1:175,175:1),c(lower1,rev(upper1)),col=alpha(pal_v[j],alpha=0.2),border=NA)
  polygon(c(1:175,175:1),c(lower2,rev(upper2)),col=alpha(pal_v[j],alpha=0.2),border=NA)
}
mtext("Year (+10000)",side=1,line=2.3)
mtext("K",side=3,at=1,line=-2)
abline(v=25,lty=3,lwd=1)
abline(h=0,lty=2)
legend("topright",legend=c("Evol","No evol"),lty=c(1,3),lwd=2,col=rep("Black",2),bg=alpha("white",alpha=0.5),box.lty=0)
legend("bottomright",legend=LETTERS[1:3],lty=rep(1,3),lwd=2,col=pal_v[1:3],bg=alpha("white",alpha=0.5),box.lty=0,title="Subpopulation",title.cex=1.1)
dev.off()
#####

####plotting shortsims####
statstmp <- popstats(test)
alldb <- create.alldb(test)
alldb <- polish.alldb(alldb)
#alldb <- filter(alldb,Rep!=4) # If excluding a divergent rep
evoltmp <- create.evoldyn(alldb,years=26:100,adult=F)
pdf(file="D:\\Migration\\paper2\\m3=0.02\\Control - staggered 0.5.pdf",width=7,height=7)
evolplots(local=evoltmp,plots=c("BVs","PEs","Ls","TEs","Dests"))
popplots(statstmp,mainfile=test,plot="endpopsize",legend=T)
dev.off()
#####

# retrieve alldb
alldb <- readRDS("D:\\Migration\\paper2\\m3=0.02\\v0=0.5, f=1\\sigmas_fine=sqrt(0.5) - Constant poor alldb.R")
h2s <- rpts <- BVs <- PEs <- Ls <- matrix(NA,max(alldb$Home),max(alldb$Rep)) # Storage matrices - values for each site and rep

#### Variance components from alldb ####
h2s <- rpts <- BVs <- PEs <- Ls <- matrix(NA,test[[3]]$N,test[[3]]$reps) # Storage matrices - values for each site and rep
t <- Sys.time()
for(i in 1:nrow(Ls)){
  for(j in 1:ncol(Ls)){
    tmp <- filter(alldb,Home==i & Year<50 & Rep==j) # Only pre-CC years. If running longer sims, set e.g. Year<1950
    BVs[i,j] <- var(tmp$BV)
    PEs[i,j] <- var(tmp$PE)
    Ls[i,j] <- var(tmp$L,na.rm=T)
    h2s[i,j] <- BVs[i,j]/Ls[i,j]
    rpts[i,j] <- (BVs[i,j]+PEs[i,j])/Ls[i,j]
  }
}
Sys.time()-t

apply(BVs,FUN=mean,MARGIN=1)
apply(PEs,FUN=mean,MARGIN=1)
apply(Ls,FUN=mean,MARGIN=1)
mean(apply(rpts,FUN=mean,MARGIN=1))
mean(apply(h2s,FUN=mean,MARGIN=1))
apply(BVs,FUN=var,MARGIN=1) # among replicate variance
#####


#For table 2
tmp <- readRDS("D:\\Migration\\paper2\\m3=0.02\\v0=0.5, f=1\\sigmas_fine=sqrt(2) stats.R")
mean(apply(tmp[[2]][3,1,180:200,],FUN=mean,MARGIN=2)) # mean in north
mean(apply(tmp[[2]][3,2,180:200,],FUN=mean,MARGIN=2)) # mean in middle
mean(apply(tmp[[2]][3,3,180:200,],FUN=mean,MARGIN=2)) # mean in south
var(apply(tmp[[2]][3,1,180:200,],FUN=mean,MARGIN=2)) # var in north
var(apply(tmp[[2]][3,2,180:200,],FUN=mean,MARGIN=2)) # var in middle
var(apply(tmp[[2]][3,3,180:200,],FUN=mean,MARGIN=2)) # var in south

#For fig 1
tmp <- readRDS("D:\\Migration\\paper2\\m3=0.02\\v0=0.5, f=1\\sigmas_coarse=sqrt(0.5) stats.R")


#### Compare two scenarios (e.g. with&without evolution) ####
alldb <- readRDS("file 1.R")
alldb2 <- readRDS("file 2.R")
evoltmp <- create.evoldyn(alldb,years=26:100,adult=F) # Remember to change years if needed.
evoltmp2 <- create.evoldyn(alldb2,years=26:100,adult=F)
df <- evoltmp[[2]]
df2 <- evoltmp2[[2]]
pop <- evoltmp[[5]]
pop2 <- evoltmp2[[5]]

pdf(file="Comparison figure.pdf",width=3.5,height=12)
par(mfrow=c(5,1),mar=c(2,4.2,1.5,0.5),oma=c(2,0,0,0))
#plot(1:175,rep(0,175),ylim=range(c(df2$yearmeanBV,df$yearmeanBV)),ylab="Breeding values",xlab="Year (+10000)",xaxt="n",cex.lab=1.2,type="n")
#axis(1,at=c(1,25*(1:7)),labels=25*(1:8))
plot(1:75,rep(0,75),ylim=range(c(df2$yearmeanBV,df$yearmeanBV)),ylab="Breeding values",xlab="Year (+10000)",xaxt="n",cex.lab=1.2,type="n")
axis(1,at=c(1,25*(1:3)),labels=25*(1:4))
for(j in 1:3){
  lines(1:75,filter(df,Home==j)$yearmeanBV,col=pal_v[j],lwd=2)
  lines(1:75,filter(df2,Home==j)$yearmeanBV,col=pal_v[j],lwd=2,lty=3)
  upper1 <- filter(df,Home==j)$yearmeanBV+filter(df,Home==j)$yearsdBV/sqrt(5)
  upper2 <- filter(df2,Home==j)$yearmeanBV+filter(df2,Home==j)$yearsdBV/sqrt(5)
  lower1 <- filter(df,Home==j)$yearmeanBV-filter(df,Home==j)$yearsdBV/sqrt(5)
  lower2 <- filter(df2,Home==j)$yearmeanBV-filter(df2,Home==j)$yearsdBV/sqrt(5)
  polygon(c(1:75,75:1),c(lower1,rev(upper1)),col=alpha(pal_v[j],alpha=0.2),border=NA)
  polygon(c(1:75,75:1),c(lower2,rev(upper2)),col=alpha(pal_v[j],alpha=0.2),border=NA)
}
mtext("A",side=3,at=1,line=-2)
abline(v=25,lty=3,lwd=1)
abline(h=0,lty=2)

#plot(1:175,rep(0,175),ylim=range(c(df2$yearmeanPE,df$yearmeanPE)),ylab="Permanent effects",xlab="Year (+10000)",xaxt="n",cex.lab=1.2,type="n")
#axis(1,at=c(1,25*(1:7)),labels=25*(1:8))
plot(1:75,rep(0,75),ylim=range(c(df2$yearmeanPE,df$yearmeanPE)),ylab="Permanent effects",xlab="Year (+10000)",xaxt="n",cex.lab=1.2,type="n")
axis(1,at=c(1,25*(1:3)),labels=25*(1:4))
for(j in 1:3){
  lines(1:75,filter(df,Home==j)$yearmeanPE,col=pal_v[j],lwd=2)
  lines(1:75,filter(df2,Home==j)$yearmeanPE,col=pal_v[j],lwd=2,lty=3)
  upper1 <- filter(df,Home==j)$yearmeanPE+filter(df,Home==j)$yearsdPE/sqrt(5)
  upper2 <- filter(df2,Home==j)$yearmeanPE+filter(df2,Home==j)$yearsdPE/sqrt(5)
  lower1 <- filter(df,Home==j)$yearmeanPE-filter(df,Home==j)$yearsdPE/sqrt(5)
  lower2 <- filter(df2,Home==j)$yearmeanPE-filter(df2,Home==j)$yearsdPE/sqrt(5)
  polygon(c(1:75,75:1),c(lower1,rev(upper1)),col=alpha(pal_v[j],alpha=0.2),border=NA)
  polygon(c(1:75,75:1),c(lower2,rev(upper2)),col=alpha(pal_v[j],alpha=0.2),border=NA)
}
mtext("D",side=3,at=1,line=-2)
abline(v=25,lty=3,lwd=1)
abline(h=0,lty=2)

#plot(1:175,seq(1,1000,length.out=175),type="n",ylab="Breeding season pop size",xlab="",ylim=c(0,800),xaxt="n",cex.lab=1.2)
#axis(1,at=c(1,25*(1:7)),labels=25*(1:8))
plot(1:75,seq(1,1000,length.out=75),type="n",ylab="Breeding season pop size",xlab="",ylim=c(330,780),xaxt="n",cex.lab=1.2)
axis(1,at=c(1,25*(1:3)),labels=25*(1:4))
lines(1:75,pop$meanb1,col=pal_v[1],lwd=2)
lines(1:75,pop$meanb2,col=pal_v[2],lwd=2)
lines(1:75,pop$meanb3,col=pal_v[3],lwd=2)
lines(1:75,pop2$meanb1,col=pal_v[1],lwd=2,lty=3)
lines(1:75,pop2$meanb2,col=pal_v[2],lwd=2,lty=3)
lines(1:75,pop2$meanb3,col=pal_v[3],lwd=2,lty=3)
polygon(c(1:75,75:1),c(pop$meanb1-pop$seb1,rev(pop$meanb1+pop$seb1)),col=alpha(pal_v[1],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop$meanb2-pop$seb2,rev(pop$meanb2+pop$seb2)),col=alpha(pal_v[2],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop$meanb3-pop$seb3,rev(pop$meanb3+pop$seb3)),col=alpha(pal_v[3],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop2$meanb1-pop2$seb1,rev(pop2$meanb1+pop2$seb1)),col=alpha(pal_v[1],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop2$meanb2-pop2$seb2,rev(pop2$meanb2+pop2$seb2)),col=alpha(pal_v[2],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop2$meanb3-pop2$seb3,rev(pop2$meanb3+pop2$seb3)),col=alpha(pal_v[3],alpha=0.2),border=NA)
mtext("G",side=3,at=1,line=-1.5)
abline(v=25,lty=3,lwd=1)
abline(h=0,lty=2)
legend("bottomleft",legend=c("Evol","No evol","North","Middle","South"),lty=c(1,3,rep(0,3)),pch=c(NA,NA,rep(15,3)),lwd=2,col=c("Black","Black",pal_v[1:3]),bg=alpha("white",alpha=0.5))

#plot(1:175,seq(1,1000,length.out=175),type="n",ylab="Non-breeding season pop size",xlab="",ylim=c(0,2000),xaxt="n",cex.lab=1.2)
#axis(1,at=c(1,25*(1:7)),labels=25*(1:8))
plot(1:75,seq(1,1000,length.out=75),type="n",ylab="Non-breeding season pop size",xlab="",ylim=c(500,1500),xaxt="n",cex.lab=1.2)
axis(1,at=c(1,25*(1:3)),labels=25*(1:4))
lines(1:75,pop$meannb1,col=pal_v[1],lwd=2)
lines(1:75,pop$meannb2,col=pal_v[2],lwd=2)
lines(1:75,pop$meannb3,col=pal_v[3],lwd=2)
lines(1:75,pop2$meannb1,col=pal_v[1],lwd=2,lty=3)
lines(1:75,pop2$meannb2,col=pal_v[2],lwd=2,lty=3)
lines(1:75,pop2$meannb3,col=pal_v[3],lwd=2,lty=3)
polygon(c(1:75,75:1),c(pop$meannb1-pop$senb1,rev(pop$meannb1+pop$senb1)),col=alpha(pal_v[1],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop$meannb2-pop$senb2,rev(pop$meannb2+pop$senb2)),col=alpha(pal_v[2],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop$meannb3-pop$senb3,rev(pop$meannb3+pop$senb3)),col=alpha(pal_v[3],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop2$meannb1-pop2$senb1,rev(pop2$meannb1+pop2$senb1)),col=alpha(pal_v[1],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop2$meannb2-pop2$senb2,rev(pop2$meannb2+pop2$senb2)),col=alpha(pal_v[2],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop2$meannb3-pop2$senb3,rev(pop2$meannb3+pop2$senb3)),col=alpha(pal_v[3],alpha=0.2),border=NA)
mtext("J",side=3,at=1,line=-1.5)
abline(v=25,lty=3,lwd=1)
abline(h=0,lty=2)

#Net deviations- formula for the variance in difference in means is (sdx-sdy)^2/n, so sd = (sdx-sdy)/sqrt(n) = sdx/sqrt(n)-sdy/sqrt(n) = se(x)-se(n)
plot(1:75,seq(1,1000,length.out=75),type="n",ylab="Difference (Evol - No evol)",xlab="",ylim=c(-100,100),xaxt="n",cex.lab=1.2)
axis(1,at=c(1,25*(1:3)),labels=seq(25,100,by=25))
#plot(1:175,seq(1,1000,length.out=175),type="n",ylab="Difference (Evol - No evol)",xlab="",ylim=c(-400,800),xaxt="n",cex.lab=1.2)
#axis(1,at=c(1,seq(25,175,by=25)),labels=seq(25,200,by=25))
lines(1:75,pop$meanb1-pop2$meanb1,col=pal_v[1],lwd=2,lty=5) # Summer
lines(1:75,pop$meannb1-pop2$meannb1,col=pal_v[1],lwd=2,lty=6) # Winter
lines(1:75,pop$meanb2-pop2$meanb2,col=pal_v[2],lwd=2,lty=5) # Summer
lines(1:75,pop$meannb2-pop2$meannb2,col=pal_v[2],lwd=2,lty=6) # Winter
lines(1:75,pop$meanb3-pop2$meanb3,col=pal_v[3],lwd=2,lty=5) # Summer
lines(1:75,pop$meannb3-pop2$meannb3,col=pal_v[3],lwd=2,lty=6) # Winter
polygon(c(1:75,75:1),c(pop$meanb1-pop2$meanb1-pop$seb1-pop2$seb2,rev(pop$meanb1-pop2$meanb1+pop$seb1+pop2$seb2)),col=alpha(pal_v[1],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop$meanb2-pop2$meanb2-pop$seb2-pop2$seb2,rev(pop$meanb2-pop2$meanb2+pop$seb2+pop2$seb2)),col=alpha(pal_v[2],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop$meanb3-pop2$meanb3-pop$seb3-pop2$seb3,rev(pop$meanb3-pop2$meanb3+pop$seb3+pop2$seb3)),col=alpha(pal_v[3],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop$meannb1-pop2$meannb1-pop$senb1-pop2$senb2,rev(pop$meannb1-pop2$meannb1+pop$senb1+pop2$senb2)),col=alpha(pal_v[1],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop$meannb2-pop2$meannb2-pop$senb2-pop2$senb2,rev(pop$meannb2-pop2$meannb2+pop$senb2+pop2$senb2)),col=alpha(pal_v[2],alpha=0.2),border=NA)
polygon(c(1:75,75:1),c(pop$meannb3-pop2$meannb3-pop$senb3-pop2$senb3,rev(pop$meannb3-pop2$meannb3+pop$senb3+pop2$senb3)),col=alpha(pal_v[3],alpha=0.2),border=NA)
abline(h=0,lty=2,lwd=1) # Zero line
abline(v=25,lty=3,lwd=1) # Shock time
mtext("M",side=3,at=1,line=-1.5)
mtext("    Year (+10000)",line=0.5,outer=T,side=1,cex=0.9)
legend("bottomleft",legend=c("A","B","C","Breeding","Non-breeding"),lty=c(rep(0,3),5,6),pch=c(rep(15,3),NA,NA),lwd=2,col=c(pal_v[1:3],"Black","Black"),bg=alpha("white",alpha=0.5))
dev.off()
#####

#### Workflow for large simulations:
# 1: Run the simulation to element tmp and save tmp as file with name "x.R".
# 2: Run function popstats(saved="x") and create file statstmp
# 3: Save statstmp as file with name "x stats.R"
# 4: Run function popplots(saved="x stats.R") and use pdf() save to PDF with name "x plots.pdf"
# 5: Run function create.alldb(readRDS(file="x")) to temporary element alldb
# 6: Run polish.alldb(alldb) to temporary element alldb
# 7: Save alldb as file with name "x alldb.R"