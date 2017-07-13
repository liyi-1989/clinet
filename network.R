library(rcd)
library(igraph)
library("maps")
library("geosphere")
source('utils_network.R')
library(scatterplot3d)
library(ggplot2)
library(glmnet)
library(e1071)
library(kernlab)

############ 1. Load GCM Data ##############
# 3D wide to 2D long
load("./data/air_mon_mean_mon_mean_removed_sub.RData")

dfv=NULL
count=0
for(i in 1:NLON){
  for(j in 1:NLAT){
    count=count+1
    dfv=rbind(dfv,c(count,i,j,LON[i],LAT[j],X[i,j,372+360+1]))
  }
}
colnames(dfv)=c("vertex","idlon","idlat","lon","lat","x")
p=dim(dfv)[1]
X1=NULL # Long Data (2D)
for(i in 1:NLON){
  for(j in 1:NLAT){
    X1=cbind(X1,X[i,j,])
  }
} # X1 is the final data matrix to work on. Vertex data frame is in dfv. Edge data frame need analysis with correlation

############ 2. Correlation Matrix ##############

S1=cor(X1)
S2=S1

for(i in 1:p){
  print(i)
  for(j in i:p){
    S2[i,j]=S2[j,i]=rcd(X1[,i],X1[,j],method="kde")
  }
}

save(S1,S2,file="cor_rcd_matrix.RData")

S1=cor(X1[id1,])
S1hat=doubeltaper(S1,NLON,NLAT,k=NLON/2,l=NLAT/2)


############ 3. Network Analysis of the features ##############
# thres=5
# Magic number: (abs(S1)>0.55 ~ S2>0.25373 ~133points
# Magic number: (abs(S1)>0.5565 ~ S2>0.25795 ~100points
# Magic number: (abs(S1)>0.597 ~ S2>0.284 ~10points
# Magic number: (abs(S1)>0.4219 ~ S2>0.19838 ~1000points
load("./data/cor_rcd_matrix.RData")

###### 3.1 Define the local nearest neighbour ######
thres=5
S0=nn_ind(NLON,NLAT,round(thres*NLON/NLAT),thres) # nearest neighbour indicator matrix

###### 3.2 Create threshoulded network for analysis ######
# threshoulded network
net1=graph_from_adjacency_matrix((abs(S1)>0.5565)*(!S0),mode = "undirected")
net2=graph_from_adjacency_matrix((S2>0.25795)*(!S0),mode = "undirected")
# threshoulded edge
dfe1=as_edgelist(net1)
dfe2=as_edgelist(net2)

###### 3.3 Visualization ######
# plot as network
par(mfrow=c(1,2))
plot_circle(net1)
plot_circle(net2)
# plot on map
plot_arc(dfv,dfe1)
plot_arc(dfv,dfe2)

png("cor.png")
par(mfrow=c(1,1))
plot_arc(dfv,dfe1)
dev.off()

png("rcd.png")
plot_arc(dfv,dfe2)
dev.off()

#----------------------
plot(X1[,232],X1[,778])
# https://www.darrinward.com/lat-long/?id=3104911

plot(X1[,805],X1[,2467])

plot(X1[,859],X1[,1518])
plot(X[24,31,],X[43,6,])
plot(X1[,859],X1[,1554])


plot(X1[,805],X1[,2431])
plot(X[23,13,],X[68,19,])

###### 3.4 Statistics of the Medium Sparse(1000) Network ######
Net1=graph_from_adjacency_matrix((abs(S1)>0.4219)*(!S0),mode = "undirected")
Net2=graph_from_adjacency_matrix((S2>0.19838)*(!S0),mode = "undirected")
Dfe1=as_edgelist(Net1)
Dfe2=as_edgelist(Net2)
par(mfrow=c(1,2))
plot_circle(Net1)
plot_circle(Net2)

# Degree Distribution 
plot(density(degree(Net1)))
lines(density(degree(Net2)),col="red")

# Betweeness 
B1=betweenness(Net1)
EB1=edge.betweenness(Net1)
B2=betweenness(Net2)
EB2=edge.betweenness(Net2)

plot(density(EB1))
lines(density(EB2),col="red")

# Community
EBC1 = edge.betweenness.community(Net1)

#member <- community.to.membership(g, eb$merges, step=nrow(eb$merges)-10L+1L)

member=membership(EBC1)

plot(Net1,
     vertex.color= rainbow(10, .8, .8, alpha=.8)[member+2],
     vertex.size=0.1, layout=layout,  vertex.label=NA,
     edge.arrow.size=0)

ec <- evcent(Net1)$vector
plot(Net1, layout=layout, vertex.size=map(ec, c(1,20)), vertex.label=NA, edge.arrow.size=.2)


############ 4. Clustering Analysis of the features ##############

D0=1-abs(S1hat)
D1=1-abs(S1)
D2=1-S2

C0=hclust(as.dist(D0))
C1=hclust(as.dist(D1))
C2=hclust(as.dist(D2))

CM0=cutree(C0, k=500) # clustering member
CM1=cutree(C1, k=500) # clustering member
CM2=cutree(C2, k=500) # clustering member

length(unique(CM1))
length(unique(CM2))

dfv0=cbind(dfv,CM0)
dfv1=cbind(dfv,CM1)
dfv2=cbind(dfv,CM2)

par(mfrow=c(2,1))
map("world",col="skyblue",border="gray10",fill=T,bg="gray30")
points(x=dfv1[,4],y=dfv1[,5],col=dfv1[,7],pch=19,cex=0.5)

map("world",col="skyblue",border="gray10",fill=T,bg="gray30")
points(x=dfv2[,4],y=dfv2[,5],col=dfv2[,7],pch=19,cex=0.5)

# given features X1, we have the cluster information of the features in dfvi
# when doing prediction, we can use representation of the clusters to do dimensionality reduction

############ 5. Summary of each cluster - dimension reduction ##############

fscluster=function(X1,dfv,C1,k=500){
  CM1=cutree(C1, k=k)
  dfv1=cbind(dfv,CM1)
  X1S=NULL
  for(i in unique(CM1)){
    idx=(dfv1[,6]==i)
    X1sub=X1[,idx]
    if(sum(idx)==1){
      tmp=X1sub
    }else{
      tmp=rowMeans(X1sub)
      #pca=prcomp(X1sub, retx=TRUE, center=TRUE, scale=TRUE)
      #tmp=pca$x[,1]
      #kpc=kpca(X1sub)
      #tmp=kpc@pcv[,1]
    }
    X1S=cbind(X1S, tmp)
  }
  return(X1S)
}

K=100
X0S=fscluster(X1,dfv,C0,K)
X1S=fscluster(X1,dfv,C1,K)
X2S=fscluster(X1,dfv,C2,K)

#save(X0S,X1S,file="reanalysis_features_clustered.RData")
############ 6. Downscaling ##############
load("../obs_data_mon_avg/monavg_tmmx_1979_2008.RData")
load("../obs_data_mon_avg/monavg_tmmx_2009_2016.RData")

# GCM has more time points: choose overlap with observation (match X, Y's n)
id1=373:(372+360) # training index
id2=(372+360+1):828 # testing index
# select (subset of) observations to do downscaling individually 
iLon=1:nlon
iLat=1:nlat
#Y1=y_train[iLon,iLat,]
#Y2=y_test[iLon,iLat,]

############ 6.1 Prediction with Lasso (or SVM) ##############
pred=function(x1,x2,Y1,Y2,iLon,iLat){
  Y1_pred=Y1
  Y2_pred=Y2
  
  set.seed(100)
  for(i in 1:length(iLon)){
    cat("lon",i,"...\n")
    for(j in 1:length(iLat)){
      y1=Y1[i,j,]
      #y2=Y2[i,j,]
      if(any(is.na(y1))==F){
        glmmod1=cv.glmnet(x1, y1, nfolds=10, alpha=1, standardize = T,parallel=F)
        Y1_pred[i,j,]=predict(glmmod1, x1, s = "lambda.min")
        Y2_pred[i,j,]=predict(glmmod1, x2, s = "lambda.min")
        #svm_model = svm(x1,y1)
        #Y1_pred[i,j,]=predict(svm_model, x1)
        #Y2_pred[i,j,]=predict(svm_model, x2)
      }
    }
  }
 
  return(list(Y1_pred=Y1_pred,Y2_pred=Y2_pred))
}

## sample code for running prediction 
## slow to run all, use cluster (in sub.R) parallel
## results for comparing lasso for S1 (cov) and S0 (tapered cov) is saved in y_test_pred_lasso_[01].RData
# RR=pred(X1[id1,],X1[id2,],Y1,Y2,iLon,iLat)
# R0=pred(X0S[id1,],X0S[id2,],Y1,Y2,iLon,iLat)
# R1=pred(X1S[id1,],X1S[id2,],Y1,Y2,iLon,iLat)
# R2=pred(X2S[id1,],X2S[id2,],Y1,Y2,iLon,iLat)

# DF=(Y2-R0$Y2_pred)^2-(Y2-R1$Y2_pred)^2
# DF=apply((Y1-R0$Y1_pred)^2,c(1,2),mean)-apply((Y1-R1$Y1_pred)^2,c(1,2),mean)
# DF=apply((Y2-R0$Y2_pred)^2,c(1,2),mean)-apply((Y2-R1$Y2_pred)^2,c(1,2),mean)
# DF1=apply((Y1-R0$Y1_pred)^2,c(1,2),mean)
# DF2=apply((Y1-R1$Y1_pred)^2,c(1,2),mean)

############ 6.2 Load Lasso Prediction Results (or svm) ##############
load("../obs_data_mon_avg/y_test_pred_lasso_0.RData")
load("../obs_data_mon_avg/y_test_pred_lasso_1.RData")
# ============= 6.2.1 comparing times =============
# each year
DF1=apply((y_test-y_test_pred_lasso_0)^2,3,mean,na.rm=T)
DF2=apply((y_test-y_test_pred_lasso_1)^2,3,mean,na.rm=T)
# each month
df1=df2=rep(0,12)
for(i in 1:12){
  df1[i]=mean(DF1[(0:3)*12+i])
  df2[i]=mean(DF2[(0:3)*12+i])
}
# each season
df10=c(mean(df1[c(1,2,12)]),mean(df1[c(3,4,5)]),mean(df1[c(6,7,8)]),mean(df1[c(9,10,11)]))
df20=c(mean(df2[c(1,2,12)]),mean(df2[c(3,4,5)]),mean(df2[c(6,7,8)]),mean(df2[c(9,10,11)]))

plot(1:96,DF1,col="red",type="l") # plot year (avg mse)
lines(1:96,DF2)
plot(1:12,df1,col="red",type="l") # plot month
lines(1:12,df2)
plot(1:4,df10,col="red",type="l") # plot season
lines(1:4,df20)

# bar plot of season mse average
df=data.frame(x=c("DJF","MAM","JJA","SON","DJF","MAM","JJA","SON"),y=c(df10,df20),type=c(rep("Double Taper",4),rep("Sample",4)))

ggplot(df, aes(x=x, y=y, fill=type)) + 
  geom_bar(stat="identity", position=position_dodge())+ 
  xlab("Season")+ylab("Testing MSE")+ggtitle("Prediction Results")+theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", x=3, y=250, label= "Overall MSE (Test)",fontface =2)+
  annotate("text", x=3, y=230, label= "Double Taper: 82.46\n Sample Cov: 117.06",fontface =1)

# #-----------------------------------------------------------------------------------------
# load("../obs_data_mon_avg/y_test_pred_svm_0.RData")
# load("../obs_data_mon_avg/y_test_pred_svm_1.RData")
# # cross-section
# DF1=apply((y_test-y_test_pred_svm_0)^2,c(1,2),mean)
# DF2=apply((y_test-y_test_pred_svm_1)^2,c(1,2),mean)
# # each year
# DF1=apply((y_test-y_test_pred_svm_0)^2,3,mean,na.rm=T)
# DF2=apply((y_test-y_test_pred_svm_1)^2,3,mean,na.rm=T)
# # each month
# df1=df2=rep(0,12)
# for(i in 1:12){
#   df1[i]=mean(DF1[(0:7)*12+i])
#   df2[i]=mean(DF2[(0:7)*12+i])
# }
# # each season
# df10=c(mean(df1[c(1,2,12)]),mean(df1[c(3,4,5)]),mean(df1[c(6,7,8)]),mean(df1[c(9,10,11)]))
# df20=c(mean(df2[c(1,2,12)]),mean(df2[c(3,4,5)]),mean(df2[c(6,7,8)]),mean(df2[c(9,10,11)]))
# 
# plot(1:4,df10,col="red",type="l")
# lines(1:4,df20)

# ============= 6.2.2 comparing cross-sections =============
# cross-section
#DF1=apply((y_test-y_test_pred_lasso_0)^2,c(1,2),mean)
#DF2=apply((y_test-y_test_pred_lasso_1)^2,c(1,2),mean)
DF1=apply((y_test[,,1:48]-y_test_pred_lasso_0[,,1:48])^2,c(1,2),mean)
DF2=apply((y_test[,,1:48]-y_test_pred_lasso_1[,,1:48])^2,c(1,2),mean)


dfDF1=NULL
count=0
for(i in 1:dim(DF1)[1]){
  for(j in 1:dim(DF1)[2]){
    count=count+1
    if((i%%4==0)&(j%%4==0)){
      dfDF1=rbind(dfDF1,c(count,lon[iLon[i]],lat[iLat[j]],y_test_pred_lasso_0[i,j,1]))
    }
    
  }
}

dfDF2=NULL
count=0
for(i in 1:dim(DF2)[1]){
  for(j in 1:dim(DF2)[2]){
    count=count+1
    if((i%%4==0)&(j%%4==0)){
      dfDF2=rbind(dfDF2,c(count,lon[iLon[i]],lat[iLat[j]],y_test_pred_lasso_1[i,j,1]))
    }
  }
}

dfDF0=NULL
count=0
for(i in 1:dim(DF2)[1]){
  for(j in 1:dim(DF2)[2]){
    count=count+1
    if((i%%4==0)&(j%%4==0)){
      dfDF0=rbind(dfDF0,c(count,lon[iLon[i]],lat[iLat[j]],y_test[i,j,1]))
    }
  }
}

colnames(dfDF0)=c("vertex","lon","lat","value")
plot_lonlat_df(dfDF0,region="usa")
#scatterplot3d(dfDF1)

dfDF1=dfDF1[!is.na(dfDF1[,4]),]
dfDF2=dfDF2[!is.na(dfDF2[,4]),]
dfDF0=dfDF0[!is.na(dfDF0[,4]),]

ft=function(x){
  ifelse(x<0,0,ifelse(x>1,1,x))
}

par(mfrow=c(2,2))

map("world",col="skyblue",border="gray10",fill=T,bg="gray30")
points(x=dfv[,4],y=dfv[,5],col=rgb((dfv[,6]-min(dfv[,6]))/(max(dfv[,6])-min(dfv[,6])),0,0),pch=19,cex=0.3,main="Reanalysis")
title("Renalysis")

map("usa",col="skyblue",border="gray10",fill=T,bg="gray30")
points(x=dfDF0[,2],y=dfDF0[,3],col=rgb((dfDF0[,4]-min(dfDF0[,4]))/(max(dfDF0[,4])-min(dfDF0[,4])),0,0),pch=19,cex=0.01)
title("Observation")

map("usa",col="skyblue",border="gray10",fill=T,bg="gray30")
points(x=dfDF1[,2],y=dfDF1[,3],col=rgb(ft((dfDF1[,4]-min(dfDF1[,4]))/(max(dfDF1[,4])-min(dfDF1[,4]))),0,0),pch=19,cex=0.01)
title("Prediction (Sample Cov.)")

map("usa",col="skyblue",border="gray10",fill=T,bg="gray30")
points(x=dfDF2[,2],y=dfDF2[,3],col=rgb(ft((dfDF2[,4]-min(dfDF2[,4]))/(max(dfDF2[,4])-min(dfDF2[,4]))),0,0),pch=19,cex=0.01)
title("Prediction (Double Tapered Cov.)")

image(DF1,axes = 0)
image(DF2)

png("test.png")
image(y_test[,,1],axes = 0)
dev.off()

################### 99. Useless code #######################

source("http://michael.hahsler.net/SMU/ScientificCompR/code/map.R")

g <- barabasi.game(1000, power=1)
layout <- layout_nicely(Net1)
layout <- layout_with_mds(Net1)
layout <- layout.fruchterman.reingold(net)
layout <- layout.sphere(Net1)
layout <- layout.circle(g)
plot(g, layout=layout, vertex.size=2,
     vertex.label=NA, edge.arrow.size=.2)


degree(g)

B=betweenness(g)
EB=edge.betweenness(g)
plot(density(EB))

plot(g, layout=layout,
     vertex.size=map(betweenness(g),c(1,15)),
     edge.width=map(edge.betweenness(g), c(1,10)),vertex.label=NA,)

eb <- edge.betweenness.community(g)

#member <- community.to.membership(g, eb$merges, step=nrow(eb$merges)-10L+1L)

member=membership(eb)

plot(g,
     vertex.color= rainbow(10, .8, .8, alpha=.8)[member+2],
     vertex.size=5, layout=layout,  vertex.label=NA,
     edge.arrow.size=.2)


ec <- evcent(g)$vector
plot(g, layout=layout, vertex.size=map(ec, c(1,20)), vertex.label=NA, edge.arrow.size=.2)



