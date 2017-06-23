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

plot_circle=function(net1){
  layout=layout.circle(net1)
  plot(net1, layout=layout, vertex.size=0.1,vertex.label=NA, vertex.frame.color="orange",vertex.color="orange",
       edge.arrow.size=0,edge.color=rgb(0,0.5,0,alpha = 0.1))
}

plot_arc=function(dfv,dfe){
  map("world",col="skyblue",border="gray10",fill=T,bg="gray30")
  points(x=dfv[unique(c(dfe[,1],dfe[,2])),4],y=dfv[unique(c(dfe[,1],dfe[,2])),5],col="orange",pch=19)
  for(i in 1:nrow(dfe)){
    node1=dfv[dfv[,1]==dfe[i,1],]
    node2=dfv[dfv[,1]==dfe[i,2],]
    arc=gcIntermediate(c(node1[4],node1[5]),c(node2[4],node2[5]),n=1000,addStartEnd = T)
    lines(arc,col="yellow")
  }
}


