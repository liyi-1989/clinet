library(MASS)
library(Matrix)

p=500
Ip=diag(rep(1,p))
t1=1
v1=c(rep(0,30),1,rep(0,40),1,rep(0,p-72))
S=Ip+t1*v1%*%t(v1)

n=100
data = mvrnorm(n, mu = rep(0,p), Sigma = S)

Shat=cov(data)
Stap=Shat*(Shat>0.5)

image(S)
image(Shat)
image(Stap)

base::norm(S-Shat*(Shat>0.5),"F")
base::norm(S-Shat,"F")

base::norm(S-Shat*(Shat>0.5),"2")
base::norm(S-Shat,"2")

#=================================================
p=7
Ip=diag(rep(1,p))
r=0.5
v1=c(r,1,r,0,r,1,r)
v2=c(0,0,0,1,0,0,0)
V=cbind(v1,v2)
V=bdiag(c(r,1,r),Ip,c(r,1,r))
V[(p+4):(p+6),1]=c(r,1,r)
V=V[,1:(1+p)]
V=bdiag(Ip,V,Ip)

p=30
Ip=diag(rep(1,p))
V=matrix(0,p,3)
r=0.5
s=c(r,1,r)
p1=5
V[(p1-1):(p1+1),1]=s
p1=27
V[(p1-1):(p1+1),1]=s
p1=7
V[(p1-1):(p1+1),3]=s
p1=9
V[(p1-1):(p1+1),2]=s
p1=20
V[(p1-1):(p1+1),2]=s


p1=5
V[(p1-1):(p1+1),2]=s
p1=10
V[(p1-1):(p1+1),1]=s
p1=15
V[(p1-1):(p1+1),3]=s

S=V%*%t(V)
image(as(V, "dgCMatrix"))
image(as(S, "dgCMatrix"))
image(as(S+Ip, "dgCMatrix"))
image(as(solve(S+Ip), "dgCMatrix"))


image(V)
image(S)

n=100
data = mvrnorm(n, mu = rep(0,dim(S)[1]), Sigma = S)

Shat=cov(data)
Stap=Shat*(abs(Shat)>0.25)

image(S)
image(as(Shat, "dgCMatrix"))
image(as(Stap, "dgCMatrix"))
image(Shat)
image(Stap)


#=================================================
V=c(r,0,0,0,0,0,0,
    r,r,r,0,1,0,0,
    0,r,r,r,0,0,0,
    0,0,r,r,r,0,0,
    0,0,0,r,r,r,0,
    0,0,0,0,r,r,r,
    0,0,0,0,0,0,r)
V=matrix(V,nrow=7,byrow=T)

H=V%*%t(V)
image(V)
image(H)
