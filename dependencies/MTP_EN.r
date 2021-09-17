###########################################
##R script for MTP-EN #####################
## Jie Liu, USC       #####################
###########################################

# The script below provides a fully self-contained classification analysis using the multiple penalty elastic net regression (MTP EN).

require(glmnet); require(mvtnorm); require(ROCR)


#Function to simulate variance-covariance matrix for the independent variables
vcm.fun=function(pho_gene,pho_meth,pho_inter,c1, c2, p1, p2){
  p=p1+p2
  vcm=array(0,c(p,p))
  for(i in 1:c1){
    vcm[i,c(i:p)]=c(1,rep(pho_gene,c(c1-i)),rep(0,c(p-c1-c2)),rep(pho_inter,c2))
  }
  
  for(i in c(p-c2+1):p){
    vcm[i,i:p]=c(1,rep(pho_meth,p-i))
  }
  for(i in c(c1+1):c(p-c2)){
    vcm[i,i]=1
  }
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      vcm[j,i]=vcm[i,j]
    }
  }
  return(vcm)
}

#Functions to simulate independent variables (multivariate normal distribution) and dependent variables (binomial dist)
simu.data=function(vcm,n,p1, p2, r1,r2r1ratio,coef1,coef2,beta0){
  p=p1+p2
  simu.x=mvtnorm::rmvnorm(n,mean=rep(0,p),sigma=vcm,method=c("eigen"))
  ##Simulate outcome variables
  r2=r1*r2r1ratio
  r0=p-r1-r2
  coeff=rep(c(coef1,0,coef2),times=c(r1,r0,r2)) 
  z=simu.x%*%coeff+beta0
  pr=1/(1+exp(-z))
  y = rbinom(n,1,pr) 
  return(list(simu.x=simu.x,y=y,pr=pr))}

#Function to calculate AUC for a penalty ratio parameter in MTP-EN model
predict.fun=function(model,x,y){
  predict.standard.high.new=predict(model,x,s="lambda.1se",type="response")
  pred.ori.high.new=ROCR::prediction(predict.standard.high.new,y)
  perf2.ori.high.new=ROCR::performance(pred.ori.high.new,"auc")
  auc.high=unlist(perf2.ori.high.new@y.values)
  return(auc.high)
}

#Function to calculate AUCs for a sequence of penalty ratio parameters in MTP-EN model
mtp.fun=function(n,p1,p2,r1,r2r1ratio,coef1,coef2,beta0,pho_gene,pho_meth,pho_inter,c1,c2, weight.seq){
  auc.record=rep(NA,length(weight.seq))
  vcm=vcm.fun(pho_gene,pho_meth,pho_inter,c1,c2, p1, p2)
  training=simu.data(vcm,n,p1, p2, r1,r2r1ratio,coef1,coef2,beta0)
  testing=simu.data(vcm,n, p1, p2, r1,r2r1ratio,coef1,coef2,beta0)
  for(i in 1:length(weight.seq)){
    lasso=glmnet::cv.glmnet(training[[1]],training[[2]], family=c("binomial"),type.measure='auc', 
                            intercept=T,penalty.factor=rep(c(1,weight.seq[i]), times=c(p1,p2)),alpha=0.5)
    auc.record[i]=predict.fun(lasso,testing[[1]],testing[[2]])
  }
  return(auc.record)
}

#Function to calculate 1- mis-classification error, sensitivity and specificity for a penalty ratio parameter
metrics.fun=function(model,x, y, r1,r2r1ratio,p1,p2){
  predict.standard.high.new=predict(model,x,s="lambda.1se",type="response")
  mypred=ifelse(predict.standard.high.new<=0.5, 0, 1)
  mis.error=sum(mypred==y)/length(y)
  
  mycoef=as.numeric(coef(model))[-1]
  sens=(sum(mycoef[1:r1]!=0)+sum(mycoef[(p1+p2-(r1*r2r1ratio)+1):(p1+p2)]!=0))/(r1+r1*r2r1ratio)
  spec=sum(mycoef[(r1+1):(p1+p2-r1*r2r1ratio)]==0)/(p1+p2-r1-r1*r2r1ratio)
  return(c(mis.error,sens,spec))
}

#Function to calculate 1- mis-classification error, sensitivity and specificity for a sequence of penalty ratio parameters 
metrics.mtp.fun=function(n,p1,p2,r1,r2r1ratio,coef1,coef2,beta0,pho_gene,pho_meth,pho_inter,c1,c2, weight.seq){
  sens.record=spec.record=mis.record=rep(NA,length(weight.seq))
  vcm=vcm.fun(pho_gene,pho_meth,pho_inter,c1,c2, p1, p2)
  training=simu.data(vcm,n,p1, p2, r1,r2r1ratio,coef1,coef2,beta0)
  testing=simu.data(vcm,n, p1, p2, r1,r2r1ratio,coef1,coef2,beta0)
  
  for(i in 1:length(weight.seq)){
    lasso=cv.glmnet(training[[1]],training[[2]], family=c("binomial"),type.measure='auc', 
                    intercept=T,penalty.factor=rep(c(1,weight.seq[i]), times=c(p1,p2)),alpha=0.5)
    pred.result=metrics.fun(lasso,testing[[1]], testing[[2]], r1,r2r1ratio, p1,p2)
    #mis.class.record[i]=pred.result[1]
    mis.record[i]=pred.result[1]
    sens.record[i]=pred.result[2]
    spec.record[i]=pred.result[3]
  }
  
  return(list(mis.record=mis.record, sensitivity=sens.record,specificity=spec.record))
}

####################################################################################################
# This is the full analysis example corresponding to Scenario 1 in the paper.


#Get AUCs for different penalty ratio parameters

weight.seq=seq(0.2, 1.8, 0.1) #In paper: weight.seq=seq(0.2, 1.8, 0.05)

#Parameters in mtp.fun:
  # sample size: 200
  # number of total features in both platforms: 250
  # number of informative features in platform 1: 5
  # number of informative features in platform 2: 5*4=20
  # coefficient for informative features in platform 1: 0.6
  # coefficient for informative features in platform 2: 0.8
  # intercept in the model: 0  
  # no correlations among features: pho_gene= pho_meth= pho_inter= 0
  # number of informative features that are correlated: here we can input any number>1 
auc=mtp.fun(200,250,250,5,4,0.6,0.8,0,0,0,0,3,3, weight.seq=weight.seq)
names(auc)=paste("Weight=", weight.seq, sep="")
auc

#Caculate 1-misclassification error, sensitivity and specificity for Scenario #1
#Penalty ratio parameters for Standard EN and MTP-EN are 1 and 0.55, respectively. 
weight.seq=c(1,0.55)

metrics.mtp.fun(200,250,250,5,4,0.6,0.8,0,0,0,0,1,1,weight.seq)