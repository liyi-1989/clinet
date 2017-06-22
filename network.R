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

nn_ind=function(NLON,NLAT,t1,t2){
  S0=matrix(0,NLON*NLAT,NLON*NLAT)
  l1=0
  for(i1 in 1:NLON){
    for(j1 in 1:NLAT){
      l1=l1+1
      l2=0
      for(i2 in 1:NLON){
        for(j2 in 1:NLAT){
          l2=l2+1
          d1=abs(i1-i2)
          if((ifelse(d1<(NLON/2),d1,NLON-d1)<=t1)|(abs(j1-j2)<=t2)){
            S0[l1,l2]=1
          }
        }
      }
    }
  }
  
  return(S0)
}



############ 3. Network Analysis of the features ##############
# Magic number: (abs(S1)>0.55 ~ S2>0.25373 ~133points
# Magic number: (abs(S1)>0.5565 ~ S2>0.25795 ~133points
# Magic number: (abs(S1)>0.597 ~ S2>0.284 ~10points
load("./data/cor_rcd_matrix.RData")
library(igraph)
#A=matrix(c(2,0,3,1),2,2)
thres=5
S0=nn_ind(NLON,NLAT,round(thres*NLON/NLAT),thres)
net1=graph_from_adjacency_matrix((abs(S1)>0.5565)*(!S0),mode = "undirected")
layout=layout.circle(net1)
plot(net1, layout=layout, vertex.size=0.1,vertex.label=NA, vertex.frame.color="blue",vertex.color="blue",
     edge.arrow.size=0,edge.color=rgb(0.5,0.5,0.5,alpha = 0.1))
plot(X1[,232],X1[,778])
# https://www.darrinward.com/lat-long/?id=3104911
thres=5
#S0=nn_ind(NLON,NLAT,round(thres*NLON/NLAT),thres)
net2=graph_from_adjacency_matrix((S2>0.25795)*(!S0),mode = "undirected")
layout=layout.circle(net2)
plot(net2, layout=layout, vertex.size=0.01,vertex.label=NA, 
     vertex.frame.color=rgb(0,0,1,alpha = 0.1),vertex.color=rgb(0,0,1,alpha = 0.1),
     edge.arrow.size=0,edge.color=rgb(0.5,0.5,0.5,alpha = 0.1))


el1=as_edgelist(net1)
el2=as_edgelist(net2)
dim(el1)
dim(el2)
el1
el2

cbind(el1,el2)

n1=graph_from_edgelist(el1)
n2=graph_from_edgelist(el2)
layout=layout.circle(n1)
layout=layout.circle(n2)
plot(n1, layout=layout,vertex.size=0.01,edge.arrow.size=0)
plot(n2, layout=layout,vertex.size=0.01,edge.arrow.size=0)


par(mfrow=c(1,2))
plot(net1, layout=layout, vertex.size=0.1,vertex.label=NA, vertex.frame.color="blue",vertex.color="blue",
     edge.arrow.size=0,edge.color=rgb(0.5,0.5,0.5,alpha = 0.1))
plot(net2, layout=layout, vertex.size=0.01,vertex.label=NA, 
     vertex.frame.color=rgb(0,0,1,alpha = 0.1),vertex.color=rgb(0,0,1,alpha = 0.1),
     edge.arrow.size=0,edge.color=rgb(0.5,0.5,0.5,alpha = 0.1))

plot(X1[,805],X1[,2467])

plot(X1[,859],X1[,1518])
plot(X[24,31,],X[43,6,])
plot(X1[,859],X1[,1554])


plot(X1[,805],X1[,2431])
plot(X[23,13,],X[68,19,])

library("maps")
library("geosphere")

png("cor.png")
par(mfrow=c(1,1))

map("world",col="skyblue",border="gray10",fill=T,bg="gray30")
points(x=dfv[unique(c(el1[,1],el1[,2])),4],y=dfv[unique(c(el1[,1],el1[,2])),5],col="orange",pch=19)
for(i in 1:nrow(el1)){
  node1=dfv[dfv[,1]==el1[i,1],]
  node2=dfv[dfv[,1]==el1[i,2],]
  arc=gcIntermediate(c(node1[4],node1[5]),c(node2[4],node2[5]),n=1000,addStartEnd = T)
  lines(arc,col="yellow")
}
dev.off()

png("rcd.png")
map("world",col="skyblue",border="gray10",fill=T,bg="gray30")
points(x=dfv[unique(c(el2[,1],el2[,2])),4],y=dfv[unique(c(el2[,1],el2[,2])),5],col="orange",pch=19)
for(i in 1:nrow(el1)){
  node1=dfv[dfv[,1]==el2[i,1],]
  node2=dfv[dfv[,1]==el2[i,2],]
  arc=gcIntermediate(c(node1[4],node1[5]),c(node2[4],node2[5]),n=1000,addStartEnd = F,breakAtDateLine=F)
  lines(arc,col="yellow")
}
dev.off()
#######################################################################

source("http://michael.hahsler.net/SMU/ScientificCompR/code/map.R")

g <- barabasi.game(1000, power=1)
layout <- layout_nicely(net)
layout <- layout_with_mds(net)
layout <- layout.fruchterman.reingold(net)
layout <- layout.sphere(g)
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



