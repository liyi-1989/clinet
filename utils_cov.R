library(MASS) # use mvnorm to generate multivariate normal
library(Matrix) # use "dgCMatrix" matrix class, so that to make plot of a matrix like hand written order
library(moments) # calculate higher order moments, skewness, kurtosis. not useful
library(pROC) # make ROC curve, not very useful currently

# Generate V in the factor model
genV=function(p=100,r=0.5,s=0.2,sp,tele=F){
  # p: number of features, assume equal number of latent factors. So V is p by p
  # r: neighbour kernel bump
  # s: tele connection bump scale
  # sp positions of teleconnection, like c(5,90)
  # tele: teleconnection signal. 1: s*(r,1,r) -1: s*(-r,1,-r)
  if(missing(sp)){
    sp=c(0.05*p,0.9*p)
  }
  V=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if((i==j) | (i-j)==2){
        V[i,j]=r
      }else if(i-j==1){
        V[i,j]=1
      }
    }
  }
  if(tele==-1){
    V[sp[2]:(sp[2]+2),sp[1]]=c(-r,1,-r)*s
  }else if(tele==1){
    V[sp[2]:(sp[2]+2),sp[1]]=c(r,1,r)*s
  }
  return(V)
}

# print a matrix with heatmap, in a hand written order, left->right, up->down (no need to rotate)
printM=function(M){
  image(as(M, "dgCMatrix"),main=deparse(substitute(M)))
}

# plot the density of a sequence and its normal approx with same mean & sd
plotDen=function(res){
  plot(density(res),col="blue",main=paste("Density of",deparse(substitute(res))))
  x=seq(min(res),max(res),length=100)
  y=dnorm(x,mean=mean(res), sd=sd(res))
  lines(x,y, type="l", lwd=1,col="red")
  abline(v=mean(res),lty=2,col="orange")
  text(mean(res)+2*sd(res),range(density(res)$y)[2]/2,paste("mean:",round(mean(res),4)),col="green")
  # fit Gumbel distribution
  #y2=dGumbel(x, mu = mean(res)-0.5772*sd(res)*sqrt(6)/pi, sigma = sd(res)*sqrt(6)/pi, log = FALSE) 
  #lines(x,y2, type="l", lwd=1,col="pink")
}

# universal threshoulding
regT=function(M,tres){
  return(M*(M>tres))
}

# adaptive threshoulding "Shat" according to "Snorm"
# (times const*sqrt(log(p)/n))
regAT=function(Shat,Snorm,n,const=2){
  p=dim(Shat)[1]
  Shat=Shat*(Shat>const*(Snorm*sqrt(log(p)/n)))
}

# Calculate "different norms" of the h-nearest-neighbour small block of a matrix M
MnormM=function(M,h=2,type=1,RM=F){
  # h is the bandwidth
  p=dim(M)[1]
  normM=M
  
  # for(i in h:(p-h+1)){
  #   for(j in h:(p-h+1)){
  for(i in 1:(p)){
    for(j in 1:(p)){    
      il=max(i-h+1,1)
      iu=min(i+h-1,p)
      jl=max(j-h+1,1)
      ju=min(j+h-1,p)
      n0=(iu-il)*(ju-jl)
      M0=M[il:iu,jl:ju]
      if(RM){
        M0=M0-mean(M0) ########
      }
      
      if(type==1){
        normM[i,j]=base::norm(M0,"2") 
      }else if(type==2){
        normM[i,j]=base::norm(M0,"F")
      }else if(type==3){
        normM[i,j]=mean(M0)
      }else if(type==4){
        normM[i,j]=mean(abs(M0))
      }else if(type==5){
        normM[i,j]=log((prod(abs(M0)))^{1/n0})
      }else if(type==6){
        normM[i,j]=sort(M0,decreasing = TRUE)[2]
      }else if(type==7){
        m0=c(M0[2,2],M0[2,1],M0[2,3],M0[1,2],M0[3,2])
        normM[i,j]=mean(m0)
      }else if(type==8){
        m0=c(M0[2,2],M0[2,1],M0[2,3],M0[1,2],M0[3,2])
        normM[i,j]=abs(prod(m0))
      }else if(type==9){
        normM[i,j]=exp(-abs(mean(M0)))
      }else if(type==10){
        normM[i,j]=sd(c(M0))
      }else if(type==11){
        normM[i,j]=skewness(c(M0))
      }else if(type==12){
        normM[i,j]=kurtosis(c(M0))
      }
      
    }
  }
  
  return(normM)
}

# variance of a covariance elements (used in Tony Cai's Adaptive threshouding) 
MnormMvar=function(X){
  n=dim(X)[1]
  Sigma_hat=cov(X)
  X1=apply(X,2,function(x){return (x-mean(x))}) # remove column mean version
  Sigma_var=((t(X1^2)%*%(X1^2)))/n-2*Sigma_hat*(t(X1)%*%X1)/n+(Sigma_hat)^2
  return(Sigma_var)
  # n=dim(X)[1]; p=dim(X)[2]
  # Xbar=as.matrix(apply(X,2,mean))
  # S1=matrix(0,p,p)
  # for(i in 1:p){
  #   for(j in 1:p){
  #     for(k in 1:n){
  #       S1[i,j]=S1[i,j]+((X[k,i]-Xbar[i,1])*(X[k,j]-Xbar[j,1])-Sigma_hat[i,j])^2
  #     }
  #     
  #   }
  # }
  # S1=S1/n
  # return(list(Sigma_var=Sigma_var,S1=S1))
}

# main simulation function
# See after "block norm", if the "teleconnection signal" could be stand out from the "random noise" (zeros). 
maxnoise=function(p=p,type=0, n=2000, N=100, r,s,tele, sigma,RM=F){
  # S: ground truth
  # type: Estimator of S. 0 is sample covariance
  # n: sample size to generate multivariate normal to get Shat (or other estimator)
  # N: repeat time
  # Goal: find sampling distribution for 
  # (1) random error for zeros in S; 
  # (2) (max) off-diagonal signal
  #p=dim(S)[1]
  Ip=diag(rep(1,p))
  #r=0.95; s=0.5; tele=T
  
  resZ=resOffDiag=res691=rep(0,N) # max for 1. zero elements, 2. off diagonal elements, 3. the teleconnection center
  V=genV(p=p,r=r,s=s,tele=tele)
  S=V%*%t(V)+diag(rep(1,p))*sigma^2
  M=(S==0)
  
  for(i in 1:p){
    for(j in 1:p){
      if(abs(i-j)<5){
        M[i,j]=0
      }
    }
  }
  
  for(k in 1:N){
    
    f=mvrnorm(n, mu = rep(0,p), Sigma = Ip)
    E=mvrnorm(n, mu = rep(0,p), Sigma = Ip)
    data = V%*%t(f)+sigma*t(E) #mvrnorm(n, mu = rep(0,p), Sigma = S)
    Shat=cov(t(data))
    if(0<type & type<=20){
      Snorm=MnormM(Shat,h=2,type=type,RM=RM)
      Shat=Snorm
      #Shat=Shat*(Shat>2*(Snorm*log(p)/n))
    }
    if(type==21){
      Shat=MnormMvar(t(data))
    }
    
    resZ[k]=max(abs(Shat*M))
    res691[k]=Shat[round(0.05*p)+1,round(0.9*p)+1] # Shat[6,91]
    # for(i in 1:p){
    #   Shat[i,i]=0
    # }
    # for(i in 1:(p-1)){
    #   Shat[i,i+1]=Shat[i+1,i]=0
    # }
    # for(i in 1:(p-2)){
    #   Shat[i,i+2]=Shat[i+2,i]=0
    # }
    # 
    # resOffDiag[k]=max(abs(Shat))
  }
  
  #x=seq(min(resZ),max(resZ),length=100)
  #y=dnorm(x,mean=mean(resZ), sd=sd(resZ))
  #y2=dGumbel(x, mu = mean(resZ)-0.5772156649*sd(resZ)*sqrt(6)/pi, sigma = sd(resZ)*sqrt(6)/pi, log = FALSE) 
  
  
  plot(density(resZ),xlim=c(min(c(resZ,res691)),max(c(resZ,res691))),col="blue",main=type)
  lines(density(res691),lty=2,type="l",col="red")
  abline(v=(mean(resZ)+mean(res691))/2,lty=3,col="orange")
  legend("topright",c("random error max","max off-diag signal","(mean1+mean2)/2"),col=c("blue","red","orange"),lty=c(1,2,3),cex=0.5)
  text(mean(resZ)+2*sd(resZ),range(density(resZ)$y)[2]/2,paste("mean:",round(mean(resZ),4)),col="blue")
  text(mean(res691)+sd(res691),range(density(res691)$y)[2]/2,paste("mean:",round(mean(res691),4)),col="red")
  
  
  #lines(x,y, type="l", lwd=1,col="green")
  #lines(x,y2, type="l", lwd=1,col="pink")
  #plotDen(resZ)
}


