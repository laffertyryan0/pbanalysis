##########
#compre inference for telomere example
#05192022
library(magrittr)
library(multipanelfigure)
library(dplyr)

#########
setwd("/Users/hongh9/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/Projects/QR-PB/Analysis")
#data=read.csv("NHANES2011to2018.csv"); dim(data)
data=read.csv("NHANES_Telo.csv")
dat=na.omit(data); dim(dat)
head(dat)
table(dat$sex)
dat.1<-dat[dat$sex==2,]
y=dat.1$telo
X.dis=dat.1$race
x=data.frame(age=dat.1$age,smoke=dat.1$smoke, bmi=dat.1$bmi, pa=dat.1$pa,edu=dat.1$edu)
wt=dat.1$wt4yr
strat=dat.1$strata
psu= ifelse(dat.1$psu>=2,2,1)
xtabs(~strat+psu)



fill1<-c("white","grey70","beige")
color1<-c("pink","black","brown")
race1<-c("Whites","Blacks","Hispanics")

#out=compare.inf(y,x,X.dis=X.dis,M="1",m="2",X.ctr="age",wt=wt,psu=psu,strat=strat,a0=agecut1[1],a1=agecut2[1],B=2,K=500, R=28,f=0.3)

source("compare_inference.r")
plot_list=list()
BB=20
B=40
R=28
f=.3
agecut1<-c(20)
agecut2<-c(90)
taus=seq(5,95,5)/100

result <- vector("list",2)
for(j in 1:2){
  i<-1;mm=c("2","3");tau=seq(5,95,5)/100
out=compare.inf(y,x,X.dis=X.dis,M="1",m=mm[j],X.ctr="age",wt=wt,psu=psu,strat=strat,a0=agecut1[i],a1=agecut2[i],B=B,BB=BB,K=500, R=R,f=0.3)
result[[j]]<-out
}
q11=result[[1]]$unconditional[1,]
q12=result[[1]]$unconditional[2,]
q22=result[[1]]$unconditional[3,]
q11.2=result[[2]]$unconditional[1,]#; should be the same as q11
q13=result[[2]]$unconditional[2,]
q33=result[[2]]$unconditional[3,]
U.1=q12-q22
U.2=q13-q33
#q11;q11.2
#q12;q22;q13;q33

q12.out=result[[1]]$q12.out; q22.out=result[[1]]$q22.out
q13.out=result[[2]]$q12.out; q33.out=result[[2]]$q22.out
U1=q12.out-q22.out; mean.U1=apply(U1,2,mean); sd.U1=apply(U1,2,sd)
U2=q13.out-q33.out; mean.U2=apply(U2,2,mean); sd.U2=apply(U2,2,sd)
UB1=mean.U1+1.96*sd.U1
LB1=mean.U1-1.96*sd.U1
UB2=mean.U2+1.96*sd.U2
LB2=mean.U2-1.96*sd.U2


q12b.out=result[[1]]$q12b.out; q22b.out=result[[1]]$q22b.out
q13b.out=result[[2]]$q12b.out; q33b.out=result[[2]]$q22b.out
U1.b=q12b.out-q22b.out; mean.U1.b=apply(U1.b,2,mean)
var0=sweep(U1.b,2,mean.U1.b)
var1=var0^2
var=(1/R)*(1/((1-f)^2))*apply(var1,2,sum)
sd.U1=sqrt(var)
UB1.b=mean.U1.b+1.96*sd.U1
LB1.b=mean.U1.b-1.96*sd.U1
U2.b=q13b.out-q33b.out; mean.U2.b=apply(U2.b,2,mean)
var0=sweep(U2.b,2,mean.U2.b)
var1=var0^2
var=(1/R)*(1/((1-f)^2))*apply(var1,2,sum)
sd.U2=sqrt(var)
UB2.b=mean.U2.b+1.96*sd.U2
LB2.b=mean.U2.b-1.96*sd.U2


q12c.out=result[[1]]$q12c.out; q22c.out=result[[1]]$q22c.out
q13c.out=result[[2]]$q12c.out; q33c.out=result[[2]]$q22c.out
U1.c=q12c.out-q22c.out; mean.U1.c=apply(U1.c,2,mean)
var0=sweep(U1.c,2,mean.U1.c)
var1=var0^2
var=(1/(BB-1))*apply(var1,2,sum)
sd.U1=sqrt(var)
UB1.c=mean.U1.c+1.96*sd.U1
LB1.c=mean.U1.c-1.96*sd.U1
U2.c=q13c.out-q33c.out; mean.U2.c=apply(U2.c,2,mean)
var0=sweep(U2.c,2,mean.U2.c)
var1=var0^2
var=(1/(BB-1))*apply(var1,2,sum)
sd.U2=sqrt(var)
UB2.c=mean.U2.c+1.96*sd.U2
LB2.c=mean.U2.c-1.96*sd.U2


out1=data.frame(x=taus,U1=U.1, U2=U.2,
                LB1=LB1,UB1=UB1, LB2=LB2, UB2=UB2,
                LB1.b=LB1.b, UB1.b=UB1.b,LB2.b=LB2.b, UB2.b=UB2.b,
                LB1.c=LB1.c, UB1.c=UB1.c,LB2.c=LB2.c, UB2.c=UB2.c)


l1=min(LB1,LB1.b,LB1.c);u1=max(UB1,UB1.b,UB1.c);
l2=min(LB2,LB2.b,LB2.c);u2=max(UB2,UB2.b,UB2.c);


p1.per<-ggplot(out1, aes(x = x, y = U1)) +
  labs(x = "Quantile",y = paste("Unexplained,", race1[2]))+lims(y=c(l1,u1))+
  geom_ribbon(aes(ymin=LB1,ymax=UB1), fill=fill1[2]) +
  geom_line(color=color1[2])

p2.per<-ggplot(out1, aes(x = x, y = U2)) +
  labs(x = "Quantile",y = paste("Unexplained,", race1[3]))+lims(y=c(l2,u2))+
  geom_ribbon(aes(ymin=LB2,ymax=UB2), fill=fill1[3]) +
  geom_line(color=color1[3])

p1.brr<-ggplot(out1, aes(x = x, y = U1)) +
  labs(x = "Quantile",y = paste("Unexplained,", race1[2]))+lims(y=c(l1,u1))+
  geom_ribbon(aes(ymin=LB1.b,ymax=UB1.b), fill=fill1[2]) +
  geom_line(color=color1[2])

p2.brr<-ggplot(out1, aes(x = x, y = U2)) +
  labs(x = "Quantile",y = paste("Unexplained,", race1[3]))+lims(y=c(l2,u2))+
  geom_ribbon(aes(ymin=LB2.b,ymax=UB2.b), fill=fill1[3]) +
  geom_line(color=color1[3])


p1.boot<-ggplot(out1, aes(x = x, y = U1)) +
  labs(x = "Quantile",y = paste("Unexplained,", race1[2]))+lims(y=c(l1,u1))+
  geom_ribbon(aes(ymin=LB1.c,ymax=UB1.c), fill=fill1[2]) +
  geom_line(color=color1[2])

p2.boot<-ggplot(out1, aes(x = x, y = U2)) +
  labs(x = "Quantile",y = paste("Unexplained,", race1[3]))+lims(y=c(l2,u2))+
  geom_ribbon(aes(ymin=LB2.c,ymax=UB2.c), fill=fill1[3]) +
  geom_line(color=color1[3])


gridExtra::grid.arrange(p1.per,p2.per,p1.brr,p2.brr, p1.boot,p2.boot, ncol = 2, nrow = 3)-> panelG
panelG

ggsave(file="inftelo.png", panelG)



