load("monavg_tmmx_1979_2008.RData")
load("monavg_tmmx_2009_2016.RData")
load("reanalysis_features_clustered.RData")
library(glmnet)

id1=373:(372+360) # training index
id2=(372+360+1):828 # testing index

iLon=1:nlon
iLat=1:nlat
Y1=y_train[iLon,iLat,]
Y2=y_test[iLon,iLat,]

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


R0=pred(X0S[id1,],X0S[id2,],Y1,Y2,iLon,iLat)
save(R0,file = "R0.RData")

R1=pred(X1S[id1,],X1S[id2,],Y1,Y2,iLon,iLat)
save(R1,file = "R1.RData")

