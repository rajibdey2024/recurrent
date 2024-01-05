library(survival);library(rsample);library(dplyr);library(tibble);library(tidyr)

cul.pop=function(data,N,boot,mixed.pois=1,frailty,predict.time,a1,p=.05){
	# data: the dataset
	# N: total number of observations
	#mixed.pois= 1, fitting frailty marginal recurrent model
	#predict.time: prediction window
	#a1: represents "r" in the Estimatnds
  set.seed(123)
  if(mixed.pois==0){
    fun=function(x,J,alpha,beta){
      exp(-alpha*exp(beta*x)+beta*x*J)*(((alpha)^J)/factorial(J))
    }
  } else {
    fun=function(x,J,alpha,beta){
      Con=gamma(J+frailty^{-1})/gamma(frailty^{-1})
      
      prob=(frailty*alpha*exp(beta*x))^{J}/(1+frailty*alpha*exp(beta*x))^{J+frailty^{-1}}
      
      Con*prob*(1/factorial(J))
    }
  }
  
  
  alpha=0.5
  a1=a1
  
  
  FP=seq(0.01,.99,.1)
  sim=boot
  boot.TCP.a.c.t=boot.FCP.a.c.t=boot.PPV.a.c.t=boot.NPV.a.c.t=matrix(NA,boot,length(FP))
  TCP.a.c.t=FCP.a.c.t=PPV.a.c.t=NPV.a.c.t=matrix(NA,sim,length(FP))
  lower.ci.TCP.a.c.t= lower.ci.FCP.a.c.t=lower.ci.PPV.a.c.t=lower.ci.NPV.a.c.t=matrix(NA,sim,length(FP))
  upper.ci.TCP.a.c.t=upper.ci.FCP.a.c.t=upper.ci.PPV.a.c.t=upper.ci.NPV.a.c.t=matrix(NA,sim,length(FP))
  area.est=matrix(NA,boot,1)
  se.TCP.a.c.t=se.FCP.a.c.t=se.PPV.a.c.t=se.NPV.a.c.t=matrix(NA,sim,length(FP))
  
  
 
  
  ps<-NULL
  simdata.full <-dnase.mod #data_final[[k1]]  
  id1=c()
  for (i in 1:dim(simdata.full)[1]){
    if(simdata.full$start[i]==0 && simdata.full$stop[i]>=predict.time){simdata.full$stop[i]=predict.time;simdata.full$status[i]=0}}
  
  for (i in 1:dim(simdata.full)[1]){
    if(simdata.full$start[i]<predict.time && simdata.full$stop[i]>predict.time){simdata.full$stop[i]=predict.time;simdata.full$status[i]=0}}
  
  simdata.full=simdata.full[simdata.full$stop<=predict.time,]
  print(simdata.full[1:20,])
  N_t=c()
  for(i in (unique(simdata.full$id))){
    N_t=c(N_t,sum(simdata.full$status[simdata.full$id==i]==1))
  }
  N_t <- N_t[!is.na(N_t)]
  
  treat=c()
  for(i in (unique(simdata.full$id))){
    treat=c(treat,unique(simdata.full$trt[simdata.full$id==i]))
  }
 
  
  for(b in 1:boot){
    
    print(b)
    D <- simdata.full %>% nest(-id)
    bs <- bootstraps(D, times = N)
    data1= as.tibble(bs$splits[[1]]) %>% unnest(cols = c(data))
    simdata=dnase.mod#data.frame(data1)
    
    if(mixed.pois==0){
      model.2 = coxph(Surv(start,stop,status) ~ x+cluster(id), method="breslow", data = simdata,timefix = FALSE)
    } else
    {
      model.2 = coxph(Surv(start,stop,status) ~ x+frailty(id), method="breslow", data = simdata,timefix = FALSE)
    }
    
    # Obtain the log of the cumulative baseline hazard function
    
    sfit <- survfit(model.2)
    h <- log(-log(sfit$surv))
    h_time <- approxfun(sfit$time, h, yleft=min(h), yright=max(h))(predict.time) # log(\Lambda_0(predict.time))
    h0_time=exp(h_time)
    
    alpha.hat=h0_time # It estimates alpha*time
    
    N_t=c()
    for(i in (unique(simdata$id))){
      N_t=c(N_t,sum(simdata$status[simdata$id==i]==1))
    }
    N_t <- N_t[!is.na(N_t)]
    print((table(N_t)))
    
    J.rec= c(0:max(((N_t))))
    
    
    #################### START: ESTIMATED VALUE CALCULATION #####################
    
    # Calculation of Total number of events up to time t (end of follow-up)
    mark=c()
    for(i in (unique(simdata$id))){
      mark=c(mark,simdata$x[simdata$id==i][1])
    }
    #quantile(mark,FP,na.rm=TRUE) #for ppv calculation
    cutoff.sam=quantile(unique(simdata$x)[N_t<a1],FP,na.rm=TRUE)
    
    Mat=matrix(NA,length(mark),length(J.rec))
    for(l in 1:length(J.rec)){
      Mat[,l]=fun(mark,J.rec[l],alpha.hat,model.2$coefficients[1])
    }
     model = coxph(Surv(start,stop,status) ~ 1, method="breslow", data = simdata,timefix = FALSE)
    sfit1 <- survfit(model)
    h1 <- log(-log(sfit1$surv))
    h1_time <- approxfun(sfit1$time, h1, yleft=min(h1), yright=max(h1))(predict.time) # log(\Lambda_0(predict.time))
    h01_time=exp(h1_time)
    
    alpha1.hat=h01_time # It estimates alpha*time
    
      Mat1=matrix(NA,length(mark),length(J.rec))
    for(l in 1:length(J.rec)){
      Mat1[,l]=fun(mark,J.rec[l],alpha1.hat,beta=0)
    }
    for(i in 1:length(cutoff.sam)){
      
      #print(i)
      mat.c=data.frame(M.c=as.numeric(mark>cutoff.sam[i]),M.nc=as.numeric(mark<=cutoff.sam[i]))
      
      mat.c.new=data.frame(mat.c,Mat)
      num.tp=sum(mat.c.new$M.c*apply(Mat[,(a1+1):max(J.rec)],1,sum))
      boot.TCP.a.c.t[b,i]=num.tp/sum(apply(Mat[,(a1+1):max(J.rec)],1,sum))
      
      boot.FCP.a.c.t[b,i]=sum(mat.c.new$M.c*apply(Mat[,1:a1],1,sum))/sum(apply(Mat[,1:a1],1,sum))
      #print(sum(apply(Mat[,(a1+1):max(J.rec)],1,sum))/N)
      #boot.PPV.a.c.t[b,i]=num.tp/sum(mat.c.new$M.c)#sum((mat.c.new$M.c)*apply(Mat[,(1):max(J.rec)],1,sum))
      p_Nt=sum(apply(Mat1[,(a1+1):max(J.rec)],1,sum)) # pr(N(t)>a1) marginal event for uninformative marker
      print(p_Nt)
      boot.PPV.a.c.t[b,i]=
        (boot.TCP.a.c.t[b,i]*p_Nt)/(sum(mat.c.new$M.c))#sum((mat.c.new$M.c)*apply(Mat[,(1):max(J.rec)],1,sum))
      
      p_Nt_comp=p_Nt=sum(apply(Mat1[,(1):a1],1,sum))
      boot.NPV.a.c.t[b,i]=((1-boot.FCP.a.c.t[b,i])*(p_Nt_comp))/(sum(mat.c.new$M.nc))
        #sum(mat.c.new$M.nc*apply(Mat[,1:a1],1,sum))/sum(mat.c.new$M.nc)
      
    } # End of i loop
    n0 = length(c(1,1-FP,0))
    fPR=c(1,1-FP,0)
    tPR=c(1,boot.TCP.a.c.t[b,],0)                    #****** Look Carefuly *************
    print(cbind(tPR,fPR))
    plot(fPR,tPR,col="blue")
    dx <- fPR[-n0] - fPR[-1]
    mid.y <- (tPR[-n0] + tPR[-1])/2
    area.est[b,] <- sum(dx * mid.y)
    #print(area.est)
  }# end of boot
  TCP.a.c.t=apply(boot.TCP.a.c.t, 2, mean)
  FCP.a.c.t=apply(boot.FCP.a.c.t, 2, mean)
  PPV.a.c.t=round(apply(boot.PPV.a.c.t, 2, mean),3)
  NPV.a.c.t=round(apply(boot.NPV.a.c.t, 2, mean),3)
  auc.a.c.t=apply(area.est,2,mean)
  se.TCP.a.c.t=apply(boot.TCP.a.c.t, 2, sd)
  se.FCP.a.c.t=apply(boot.FCP.a.c.t, 2, sd)
  se.PPV.a.c.t=round(apply(boot.PPV.a.c.t, 2, sd),3)
  se.NPV.a.c.t=round(apply(boot.NPV.a.c.t, 2, sd),3)
  se.auc.a.c.t=apply(area.est,2,sd)
  
  lower.ci.TCP.a.c.t= apply(boot.TCP.a.c.t, 2, function(x){quantile(x,p,na.rm=T)})
  lower.ci.FCP.a.c.t=apply(boot.FCP.a.c.t, 2, function(x){quantile(x,p,na.rm=T)})
  lower.ci.PPV.a.c.t=apply(boot.PPV.a.c.t, 2, function(x){quantile(x,p,na.rm=T)})
  lower.ci.NPV.a.c.t=apply(boot.NPV.a.c.t, 2, function(x){quantile(x,p,na.rm=T)})
  lower.ci.auc.a.c.t=apply(area.est, 2, function(x){quantile(x,p,na.rm=T)})
  
  upper.ci.TCP.a.c.t= apply(boot.TCP.a.c.t, 2, function(x){quantile(x,1-p,na.rm=T)})
  upper.ci.FCP.a.c.t=apply(boot.FCP.a.c.t, 2, function(x){quantile(x,1-p,na.rm=T)})
  upper.ci.PPV.a.c.t=apply(boot.PPV.a.c.t, 2, function(x){quantile(x,1-p,na.rm=T)})
  upper.ci.NPV.a.c.t=apply(boot.NPV.a.c.t, 2, function(x){quantile(x,1-p,na.rm=T)})
  upper.ci.auc.a.c.t=apply(area.est, 2, function(x){quantile(x,1-p,na.rm=T)})
  
  #} # End of Simulation
  ps=cbind(FPR=FP,TP=TCP.a.c.t,TP.se=se.TCP.a.c.t,FP=FCP.a.c.t,
           FP.se=se.FCP.a.c.t,ppv=PPV.a.c.t,ppv.se=se.PPV.a.c.t,
           npv=NPV.a.c.t,Npv.se=se.NPV.a.c.t,AUC=auc.a.c.t,AUC.se=se.auc.a.c.t)
  CI=cbind(l.TP=lower.ci.TCP.a.c.t, l.FP=lower.ci.FCP.a.c.t,l.PPV=lower.ci.PPV.a.c.t,l.NPV=lower.ci.NPV.a.c.t,l.auc=lower.ci.auc.a.c.t,u.TP=upper.ci.TCP.a.c.t,u.FP=upper.ci.FCP.a.c.t,u.PPV=upper.ci.PPV.a.c.t,u.NPV=upper.ci.NPV.a.c.t,u.auc=upper.ci.auc.a.c.t)
  simulated.out=cbind(ps,CI)
  return(simulated.out)
  
}

# Data Pre-processing
library(survival)
first <- subset(rhDNase, !duplicated(id)) #first row for each subject
dnase <- tmerge(first, first, id=id, tstop=as.numeric(end.dt -entry.dt))

# Subjects whose fu ended during the 6 day window are the reason for
#  this next line
temp.end <- with(rhDNase, pmin(ivstop+6, end.dt-entry.dt))
dnase <- tmerge(dnase, rhDNase, id=id,
                infect=event(ivstart),
                end=  event(temp.end))
# toss out the non-at-risk intervals, and extra variables
#  3 subjects had an event on their last day of fu, infect=1 and end=1


id=dnase$id
inst=dnase$inst
trt=dnase$trt
x=-dnase$fev
start=dnase$tstart
stop=dnase$tstop
status=dnase$infect
dnase.mod=data.frame(id,inst,trt,x,start,stop,status)
head(dnase.mod,n=20)

# application of the function

aaa=cul.pop(data=dnase.mod,N=647,a1=2,boot=2,mixed.pois=1,frailty=0.5,predict.time = 196,p=.025)

aaa
