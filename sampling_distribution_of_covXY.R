source('utils_cov.R')

############ 1. NO teleconnection case ############
# True model (NO teleconnection)
p=100; r=0.5
V=genV(p=p,r=0.5,tele=F) # only neighbour effect
S=V%*%t(V)+diag(rep(1,p)) # true cov, generated from factor model S=VV'+I
printM(V)
printM(S)

############ 1.0 sampling dis. of cov of two (univariate) ind. normal ############
# Given X, and Y are two independent normal (sample data with sample size 10^N2),
# what is the distribution of cov(X,Y)? (of course, the true value cov(X,Y)=0).
# This part simply repeat it N times, and find the sampling distribution of cov(X,Y). 
# This is a preparation of the zero(true value) elements in the sample covariance matrix. 

N=1000 # repeat size
N2=5 # sample size: 10 to the power
res=rep(0,N) # each iteration, save the covariance. all repeatants forms a sample
res2=rep(0,N2) # each sample size, save sd

for(j in 1:N2){
  n=10^j # sample size 
  for(i in 1:N){
    x=rnorm(n)
    y=rnorm(n)
    res[i]=cov(x,y)
  }
  res2[j]=sd(res) # standard error of the "sample covariance"
}

# sample size vs sd (of sampling dis.) and 1/sqrt(n)
plot(1:N2,log(res2),main="sd of sampling dis. of cov(true value = 0)", 
     xlab="log(sample size)",ylab="log(standard error)")
lines(1:N2,log(1/sqrt(10^(1:N2))),type="l",col="red")
legend("topright",c("log(sd)","log(1/sqrt(n))"),lty=c(1,1),pch = c(1,NA),col=c("black","red"))

# sampling distribution vs its normal version(same mean & sd)
plot(density(res),main="sampling dis. of cov ind. two normal")
x=seq(min(res),max(res),length=100)
y=dnorm(x,mean=mean(res), sd=sd(res))
lines(x,y, type="l", lwd=1,col="red")
legend("topright",c("sampling dis.","its normal appr."),lty=c(1,1),col=c("black","red"))

############ 1.1 dis of zero/non-0s in sample cov ############
# (for one fixed sample cov) zero elements(off-diagonals) are the sample data
n=1000
data = mvrnorm(n, mu = rep(0,p), Sigma = S)
Shat=cov(data)

res=NULL
for(i in 1:p){
  for(j in 1:p){
    if((i-j)>2){
      res=c(res,Shat[i,j])
    }
  }
}

# sampling distribution vs its normal version(same mean & sd)
plot(density(res),col="blue",main="dis of zeros(true value) in sample cov")
x=seq(min(res),max(res),length=100)
y=dnorm(x,mean=mean(res), sd=sd(res))
lines(x,y, type="l", lwd=1,col="red")
legend("topright",c("dis. of zeros","its normal appr."),lty=c(1,1),col=c("black","red"))



############ 1.2 dis of max zero/non-0 in many sample cov ############
N=1000
res2=rep(0,N)
for(k in 1:N){
  data = mvrnorm(n, mu = rep(0,p), Sigma = S)
  Shat=cov(data)
  # for(i in 1:p){
  #   Shat[i,i]=0
  # }
  # for(i in 1:(p-1)){
  #   Shat[i,i+1]=Shat[i+1,i]=0
  # }
  # for(i in 1:(p-2)){
  #   Shat[i,i+2]=Shat[i+2,i]=0
  # }
  res2[k]=max(abs(Shat*(S==0)))
}

plot(density(res2),lty=2,col="blue",main="dis of max of zeros in many sample cov")
x=seq(min(res2),max(res2),length=100)
y=dnorm(x,mean=mean(res2), sd=sd(res2))
lines(x,y, type="l", lwd=1,col="red")
abline(v=mean(res2),lty=2,col="orange")
legend("topright",c("dis. of max of zeros","its normal appr."),lty=c(2,1),col=c("blue","red"))
y2=dGumbel(x, mu = mean(res2)-0.5772*sd(res2)*sqrt(6)/pi, sigma = sd(res2)*sqrt(6)/pi, log = FALSE) 
lines(x,y2, type="l", lwd=1,col="pink")

summary(res2)
quantile(res2,0.99)



############ 2. WITH teleconnection case ############
# True model (WITH teleconnection)
p=100; n=1000; r=0.9; s=0.5; tele=-1; sigma=1
Ip=diag(rep(1,p))
V=genV(p=p,r=r,s=s,tele=tele)
S=V%*%t(V)+diag(rep(1,p))
M=(S==0)
################### 2.1 one sample cov, distr of zero elements ###################
n=1000
data = mvrnorm(n, mu = rep(0,p), Sigma = S)
Shat=cov(data)
# take zero elements in true cov
res=NULL

for(i in 1:p){
  for(j in 1:p){
    if(M[i,j] & (i-j>2)){
      res=c(res,Shat[i,j])
    }
  }
}

# result: dis. of zeros like normal
plotDen(res) 

################### 2.2 multi sample cov, distr of max elements ###################
# compare the distribution  of:
# (1). max of zeros in true cov (Well, I exclude the margin 2<|i-j|<5, 
#      this part maybe large in the 3 by 3 block version)
# (2). max of teleconnection signal (max off-diagonal non-zero center)
#
# type: 0: sample cov; 1: spectual norm; 2: Frobenius norm; 3: mean; 4: abs(mean); 
#       5: prod; 6: 2nd largest; 8: abs(prod(m0)); 10: sd; 11: skewness; 12: kurtosis;
# r: rho; s: scale of the teleconnection s*(rho,1,rho); 
p=100; n=1000; r=0.9; s=0.5; tele=T; sigma=1

for(i in c(0,4)){ # type = 1,2 seems also work
  maxnoise(p,type=i, n=1000, N=100,r=r,s=s,tele=tele, sigma=sigma)
}
# results: for mean removed M0(3 by 3 block), 1,2,4,6,10 works
# especially: Tony Cai's adaptive threshoulding matrix not work here (type=21)

################### 2.3 covariance estimation ###################
# e.g. s=2, the results usually random #set.seed(1234567)

p=100; n=1000; r=0.9; s=0.4; tele=T; sigma=1
Ip=diag(rep(1,p))
V=genV(p=p,r=r,s=s,tele=tele)
S=V%*%t(V)+diag(rep(1,p))

# Step I: Generate data 
# method 1: generate data directly with sigma
data = mvrnorm(n, mu = rep(0,p), Sigma = S)
# method 2: generate according to factor model
V=genV(p=p,r=r,s=s,tele=tele)
f=mvrnorm(n, mu = rep(0,p), Sigma = Ip)
E=mvrnorm(n, mu = rep(0,p), Sigma = Ip)
data = V%*%t(f)+sigma*t(E) #mvrnorm(n, mu = rep(0,p), Sigma = S)
data=t(data)

# Step II: Calculate Different Covariance Estimator
# NOTE:  each one has a threshould to tune based on different model

# TYPE 1: sample covariance
Shat=cov(data) 
# TYPE 2: (universal) threshoulding sample covariance
ST=Shat*(abs(Shat)>0.4) 
# note: large sample size n=1e5, thre=0.1, frob_thre=0.2, could find all, better than threshoulding (not not very interesting)
# TYPE 3: use (3 by 3 block) variance threshoulding
Snorm=MnormM(Shat,h=2,type=4)
B1=(Snorm>0.35)
B2=matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    if(abs(i-j)<=10){
      B1[i,j]=0
    }
    if(abs(i-j)<=2){
      B2[i,j]=1
    }
  }
}
# SB=Shat*(B1 | B2) # B1: teleconnection; B2: banding, cut nearest-neighbours (too good)
SB=Shat*(Snorm>0.3)
# TYPE 4: Tony Cai's Adaptive threshould
SAT=MnormMvar(data) # Calculate variance of covariance elements
SA=regAT(Shat,sqrt(SAT),n,const=1.5)

#printM(Shat)
printM(ST)
printM(SB)
printM(SA)
#printM(S)
base::norm(S-Shat,"F")
base::norm(S-ST,"F")
base::norm(S-SB,"F")
base::norm(S-SA,"F")



round(S[4:8,89:93],2)
round(Shat[4:8,89:93],2)
round(Snorm[4:8,89:93],2)


################### 2.4 others ###################
# example of ROC (not very useful, currently.)
df=NULL

for(i in 1:p){
  for(j in 1:p){
    if((i-j)>=0){
      #t1=c(i,j,S[i,j]!=0,abs(ST[i,j]),abs(SB[i,j]),abs(SA[i,j]))
      t1=c(i,j,S[i,j]!=0,abs(Shat[i,j]),abs(Snorm[i,j]))
      df=rbind(df,t1)
    }
  }
}

colnames(df)=c("i","j","truth","ST","SB")

roc_ST=roc(response=df[,"truth"], predictor=df[,"ST"], positive = 1,smooth = F)
roc_SB=roc(response=df[,"truth"], predictor=df[,"SB"], positive = 1,smooth = F)
#roc_SA=roc(response=df[,"truth"], predictor=df[,"SA"], positive = 1,smooth = T)
plot(roc_ST,col="blue")
plot(roc_SB,col="red",add=TRUE)
plot(roc_SA,col="green",add=TRUE)

data.frame(roc_ST$sensitivities, roc_ST$specificities, roc_ST$thresholds)

# example of extreme distribution
library(evd) #fgev(...)$estimate to fit model parameters
library(QRM) # xGumbel like function

rGumbelSim = rGumbel(1000,mu=0,sigma=1)
plot(density(rGumbelSim))

fgev(rGumbelSim)$estimate
fgev(rGumbelSim,shape=0)$estimate








