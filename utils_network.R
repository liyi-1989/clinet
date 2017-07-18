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

# plot the map data (lon*lat, lon, lat)
plot_lonlat=function(X,lon,lat,vcol="value",region="world",cap=""){
  # Input:
  # X: a map lon*lat, 2d array
  # lon, lat: 1d array corresponds to X
  # Output:
  # dfv: vertex dataframe:
  #       vertex: vertex id
  #       lon,lat: lon,lat for each point
  #       value: value for each point
  nlon=length(lon)
  nlat=length(lat)
  
  dfv=NULL
  count=0
  for(i in 1:nlon){
    for(j in 1:nlat){
      count=count+1
      dfv=rbind(dfv,c(count,i,j,lon[i],lat[j],X[i,j]))
    }
  }
  colnames(dfv)=c("vertex","idlon","idlat","lon","lat",vcol)
  plot_lonlat_df(dfv,vcol=vcol,region=region,cap=cap)
  return(dfv)
}

# plot vertex dataframe (vertex,lon,lat,value)
plot_lonlat_df=function(dfv,vcol="value",region="world",cap="",CEX=0.3){
  # Input:
  # dfv: vertex dataframe, with lon,lat,vcol
  # vcol: column name for the vertex value column
  # region: world/usa, etc
  # cap: title of the plot
  # Output: 
  # just plot points on map
  
  x = dfv[,vcol]
  vcolor = (x-min(x))/(max(x)-min(x))
  vred=(vcolor>0.5)*2*(vcolor-0.4)/1.2
  vblue=(vcolor<=0.5)*2*(0.6-vcolor)/1.2
  map(region,col="skyblue",border="gray10",fill=T,bg="gray30")
  points(x=dfv[,"lon"],y=dfv[,"lat"],col=rgb(vcolor,0,0),pch=15,cex=CEX,main="Reanalysis")
  if(region=="world"){
    plot(wrld_simpl,add=T)
  }
  title(cap)
}





