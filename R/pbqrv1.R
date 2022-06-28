#Last modified 5/17/2022
#rm(list=ls(all=TRUE))
library(stats)
library(survey)
library(quantreg)
library(splines)
library(Hmisc)

#data
#from runpbqr_bmi.r
#y=dat$bmi
#X.dis=dat$race
#x=data.frame(age=dat$age,pir=dat$pir,insurance=dat$insurance,sex=dat$sex)
#wt=dat$wt8yr
#strat=dat$strata
#psu= dat$psu


PBQRP=function(y,x, X.dis, M,m,X.ctr,a0,a1, B, K, wt,psu,strat){
  set.seed(123456)
  newdat=data.frame(y=y,x,wt=wt,psu=psu,strat=strat)
  taus=c(seq(.05,.95,.05))
  n.taus=length(taus)
  theta<-seq(1,K,1)/(K+1)
  
  q11.out=matrix(0,nrow=B,ncol=n.taus) #to store white quantile results at .25, .5,.75
  q12.out=matrix(0,nrow=B,ncol=n.taus) 
  q22.out=matrix(0,nrow=B,ncol=n.taus) 
  #q11
  temp.data<-newdat[which(X.dis==M),]
  temp.data=temp.data[temp.data[X.ctr]>=a0&temp.data[X.ctr]<=a1,]
  Q1<-wtd.quantile(temp.data$y, probs=taus, weights=temp.data$wt)
  
  wt.normalized<-temp.data$wt/sum(temp.data$wt) #sum(wt.normalized)
  cumsum.wt.normalized<-cumsum(wt.normalized) #cumsum.wt.normalized
  u.x<-runif(K); 
  temp.index<-NULL
  for (k in 1:K) { temp.index<-c(temp.index, sum(cumsum.wt.normalized<u.x[k])+1) }
  cat1.X<-temp.data[temp.index,]
  
  
  results<-NULL
  temp.data<-newdat[which(X.dis==M),]
  for (k in 1:K) {  
    xx=as.matrix(temp.data[,c(colnames(x))])
  temp<-rq(y~xx, weights=wt, data=temp.data, tau=theta[k])
  results<-rbind(results, temp$coef)
  }
  cat1.results<-results
  for (v in 1:ncol(results)){
  ttemp<-lm(cat1.results[,v]~ns(theta, df=10, intercept=F)); 
  cat1.results[,v]<-ttemp$fitted
}

  
  
 #q12 and q22
  temp.data<-newdat[which(X.dis==m),]
  results<-NULL
  for (k in 1:K) { 
    xx=as.matrix(temp.data[,c(colnames(x))])
    temp<-rq(y~xx, weights=wt, data=temp.data, tau=theta[k])
    results<-rbind(results, temp$coef)
  }
  cat2.results<-results
  for (v in 1:ncol(results)){
    ttemp<-lm(cat2.results[,v]~ns(theta, df=10, intercept=F)); 
    cat2.results[,v]<-ttemp$fitted
  }
  
  temp.data<-newdat[which(X.dis==m),]
  temp.data=temp.data[temp.data[X.ctr]>=a0&temp.data[X.ctr]<=a1,]
  Q2<-wtd.quantile(temp.data$y, probs=taus, weights=temp.data$wt)
  wt.normalized<-temp.data$wt/sum(temp.data$wt) #sum(wt.normalized)
  cumsum.wt.normalized<-cumsum(wt.normalized) #cumsum.wt.normalized
  u.x<-runif(K); 
  temp.index<-NULL
  for (k in 1:K) { temp.index<-c(temp.index, sum(cumsum.wt.normalized<u.x[k])+1) }
  cat2.X<-temp.data[temp.index,]
  
   q11<-apply(as.matrix(cat1.X[,c(colnames(x))])*cat1.results[,-1], 1, sum)+cat1.results[,1]
   q22<-apply(as.matrix(cat2.X[,c(colnames(x))])*cat2.results[,-1], 1, sum)+cat2.results[,1]
   q12<-apply(as.matrix(cat2.X[,c(colnames(x))])*cat1.results[,-1], 1, sum)+cat1.results[,1]
  
  
  unconditional<-rbind(quantile(q11, prob=taus), 
                       quantile(q12, prob=taus),
                       quantile(q22, prob=taus))
  
  rownames(unconditional)=c("q11","q12","q22")
  marginal=rbind(Q1,Q2)
  rownames(marginal)=c("q1","q2")
  ###########################################################
  # Inference of q11-q12: explained disparity
  #compute B of q11
  for(b in 1:B){
    set.seed(b)
    R.hj=rep(0, nrow(newdat))
    S=unique(newdat$strat)
    J=unique(newdat$psu)
    for (h in 1:length(S)){
      for (j in 1:length(J)){
      index=which(newdat$strat==S[h]&newdat$psu==J[j])
      R.hj[index]<-rexp(n=1,rate=1)
      #samples within the same strata & PSU will have the same random R
      }
    }

    newdat$R.hj<-R.hj
    new.wt=newdat$R.hj*newdat$wt #samples within the same strata would have the same R
    newdat$new.wt=new.wt
  
    #q11.b
    temp.data<-newdat[which(X.dis==M),]
    xx=as.matrix(temp.data[,c(colnames(x))])
    results<-NULL
    for (k in 1:K) { 
      temp<-rq(y~xx,weights=new.wt, data=temp.data, tau=theta[k])
    results<-rbind(results, temp$coef)
    }
    cat1.results<-results;#dim(cat1.results)
    for (v in 1:ncol(results)){
      ttemp<-lm(cat1.results[,v]~ns(theta, df=10, intercept=F)); 
      cat1.results[,v]<-ttemp$fitted
    }

    temp.data=temp.data[temp.data[X.ctr]>=a0&temp.data[X.ctr]<=a1,]
    wt.normalized<-temp.data$wt/sum(temp.data$wt)
    cumsum.wt.normalized<-cumsum(wt.normalized);#cumsum.wt.normalized
    u.x<-runif(K); 
    temp.index<-NULL
    for (k in 1:K) { temp.index<-c(temp.index, sum(cumsum.wt.normalized<u.x[k])+1) }
    cat1.X<-temp.data[temp.index,]
   
    q11<-apply(as.matrix(cat1.X[,c(colnames(x))])*cat1.results[,-1], 1, sum)+cat1.results[,1]
    q11.out[b,]<-quantile(q11, prob=taus)
    
    ##compute B of q12 (i.e.,X.black *beta.white)
    temp.data<-newdat[which(X.dis==m),]
    temp.data=temp.data[temp.data[X.ctr]>=a0&temp.data[X.ctr]<=a1,]
    wt.normalized<-temp.data$wt/sum(temp.data$wt)
    cumsum.wt.normalized<-cumsum(wt.normalized)
    u.x<-runif(K); 
    temp.index<-NULL
    for (k in 1:K) { temp.index<-c(temp.index, sum(cumsum.wt.normalized<u.x[k])+1) }
    cat2.X<-temp.data[temp.index,]
  
    q12<-apply(as.matrix(cat2.X[,c(colnames(x))])*cat1.results[,-1], 1, sum)+cat1.results[,1]
    q12.out[b,]<-quantile(q12, prob=taus)
    
    #compute B of q22 
    results<-NULL
    temp.data<-newdat[which(X.dis==m),]
    xx=as.matrix(temp.data[,c(colnames(x))])
    for (k in 1:K) {  
       temp<-rq(y~xx,weights=new.wt, data=temp.data, tau=theta[k])
      results<-rbind(results, temp$coef)
    }
    cat2.results<-results;#dim(cat1.results)
    for (v in 1:ncol(results)){
      ttemp<-lm(cat2.results[,v]~ns(theta, df=10, intercept=F)); 
      cat2.results[,v]<-ttemp$fitted
    }
    
    
    q22<-apply(as.matrix(cat2.X[,c(colnames(x))])*cat2.results[,-1], 1, sum)+cat2.results[,1]
    q22.out[b,]<-quantile(q22, prob=taus)
    print(b)
  }

  out=list(unconditional=unconditional, marginal=marginal, q11.out=q11.out, q12.out=q12.out,q22.out=q22.out)
  out
  }



