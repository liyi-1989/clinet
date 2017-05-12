library(igraph)

id2ij=function(id,ny){
  # id=(i-1)*ny+j  
  j=id%%ny
  i=id%/%ny+1
  if(i==0){
    j=ny
    i=i-1
  }
  return(c(i,j))
}

allfiles=list.files(path="../data/results_cor",pattern = "*.csv")
nfiles=length(allfiles)

df=NULL
score=NULL
for(i in 2:nfiles){
  X0=read.csv(paste0("../data/results_cor/",allfiles[i]),header = F)
  score=c(score,X0[,7])
  NODE_FROM=(X0[,1]-1)*73+X0[,2]
  NODE_TO=(X0[,3]-1)*73+X0[,4]
  df=rbind(df,cbind(NODE_FROM,NODE_TO,X0[,7]))
}

cutoff=quantile(score,0.99)
df0=df[df[,3]>cutoff,]

net <- graph_from_data_frame(d=df0, vertices=sort(unique(c(df0[,1],df0[,2]))), directed=F) 
deg <- degree(net, mode="all")
V(net)$size <- deg
V(net)$size[vcount(net)]=5
l=layout_in_circle(net)
l=layout_on_sphere(net)
l=layout_with_fr(net)
plot(net,vertex.size=1, vertex.label.cex=0.75, vertex.color="blue",vertex.frame.color=NA,
     edge.width=0.5,layout=l)
 


netm <- as_adjacency_matrix(net, sparse=F)
palf <- colorRampPalette(c("gold", "dark orange")) 
heatmap(netm[,17:1], Rowv = NA, Colv = NA, col = palf(20), 
        scale="none", margins=c(10,10) )

# Degree distribution
deg.dist <- degree_distribution(net, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.4, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")


sort(closeness(net), decreasing = TRUE)
sort(betweenness(net),decreasing = TRUE)
sort(evcent(net)[[1]], decreasing = TRUE)
