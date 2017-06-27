############ 1. Load GCM Data ##############
# 3D wide to 2D long
load("./data/air_mon_mean_mon_mean_removed_sub.RData")

dfv=NULL
count=0
for(i in 1:NLON){
  for(j in 1:NLAT){
    count=count+1
    dfv=rbind(dfv,c(count,i,j,LON[i],LAT[j]))
  }
}
colnames(dfv)=c("vertex","idlon","idlat","lon","lat")
p=dim(dfv)[1]
X1=NULL # Long Data (2D)
for(i in 1:NLON){
  for(j in 1:NLAT){
    X1=cbind(X1,X[i,j,])
  }
} # X1 is the final data matrix to work on. Vertex data frame is in dfv. Edge data frame need analysis with correlation

############ 2. Correlation Matrix ##############
library(rcd)
S1=cor(X1)
S2=S1

for(i in 1:p){
  print(i)
  for(j in i:p){
    S2[i,j]=S2[j,i]=rcd(X1[,i],X1[,j],method="kde")
  }
}

save(S1,S2,file="cor_rcd_matrix.RData")

############ 3. Network Analysis of the features ##############
# thres=5
# Magic number: (abs(S1)>0.55 ~ S2>0.25373 ~133points
# Magic number: (abs(S1)>0.5565 ~ S2>0.25795 ~100points
# Magic number: (abs(S1)>0.597 ~ S2>0.284 ~10points
# Magic number: (abs(S1)>0.4219 ~ S2>0.19838 ~1000points
load("./data/cor_rcd_matrix.RData")
library(igraph)
library("maps")
library("geosphere")
source('utils_network.R')
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

D1=1-abs(S1)
D2=1-S2

C1=hclust(as.dist(D1))
C2=hclust(as.dist(D2))

CM1=cutree(C1, k=500) # clustering member
CM2=cutree(C2, k=500) # clustering member

length(unique(CM1))
length(unique(CM2))

dfv1=cbind(dfv,CM1)
dfv2=cbind(dfv,CM2)

par(mfrow=c(2,1))
map("world",col="skyblue",border="gray10",fill=T,bg="gray30")
points(x=dfv1[,4],y=dfv1[,5],col=dfv1[,6],pch=19,cex=0.5)

map("world",col="skyblue",border="gray10",fill=T,bg="gray30")
points(x=dfv2[,4],y=dfv2[,5],col=dfv2[,6],pch=19,cex=0.5)

# given features X1, we have the cluster information of the features in dfvi
# when doing prediction, we can use representation of the clusters to do dimensionality reduction

############ 5. Summary of each cluster - dimension reduction ##############
X1S=NULL
for(i in unique(CM1)){
  idx=(dfv1[,6]==i)
  X1sub=X1[,idx]
  if(sum(idx)==1){
    tmp=X1sub
  }else{
    #tmp=rowMeans(X1sub)
    pca=prcomp(X1sub, retx=TRUE, center=TRUE, scale=TRUE)
    tmp=pca$x[,1]
  }
  X1S=cbind(X1S, tmp)
}

X2S=NULL
for(i in unique(CM2)){
  idx=(dfv2[,6]==i)
  X1sub=X1[,idx]
  if(sum(idx)==1){
    tmp=X1sub
  }else{
    #tmp=rowMeans(X1sub)
    pca=prcomp(X1sub, retx=TRUE, center=TRUE, scale=TRUE)
    tmp=pca$x[,1]
  }
  X2S=cbind(X2S, tmp)
}

############ 6. Downscaling ##############
library(glmnet)
load("../obs_data_mon_avg/monavg_tmmx_1979_2008.RData")
load("../obs_data_mon_avg/monavg_tmmx_2009_2016.RData")

x1=X1[373:(372+360),]
x2=X1[(372+360+1):828,]

x1=X1S[373:(372+360),]
x2=X1S[(372+360+1):828,]
ilon=500
ilat=500
y1=y_train[ilon,ilat,]
y2=y_test[ilon,ilat,]


set.seed(100)
glmmod1=cv.glmnet(as.matrix(x1), y1, nfolds=10, alpha=1, standardize = FALSE,parallel=F)
predTest <- predict(glmmod1, as.matrix(x2), s = "lambda.min")
predTest_train <- predict(glmmod1, as.matrix(x1), s = "lambda.min")

mean((predTest_train-y1)^2)
mean((predTest-y2)^2)


#######################################################################

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



