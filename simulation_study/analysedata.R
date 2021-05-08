##### Set environment #####
rm(list=ls());graphics.off()
setwd("simulation_study/")
source("utilities.R"); load("design.rda")


##### Get estimates and compute statistics #####
B=5000
Rho = array(NA,c(NROW(design),8,3),dimnames = list(NULL,c("mean","var","bias","rmse","meanSE","NAs","ratioBias","propPosBias"),c("fEM","maxdef","meandef")))
TauX1 = TauX2 = array(NA,c(NROW(design),5,3),dimnames = list(NULL,c("var","bias","rmse","ratioBias","propPosBias"),c("fEM","maxdef","meandef")))
Tau_fEM = Tau_maxdef = Tau_meandef = array(NA,c(NROW(design),2,10),dimnames = list(NULL,c("mean","var"),NULL))

iid=1:NROW(design)
for(k in iid){
  print(paste0("Processing design: ",k))
  
  load(paste0("/home/antonio/mount_server/estimdata_",k,".rda"))
  rho0=design$rho[k]
  
  Rho[k,"NAs",] = apply(Y[,"rho",],1,function(x)sum(is.na(x)))
  Rho[k,"mean",] = apply(Y[,"rho",],1,function(x)mean(x,na.rm=TRUE))
  Rho[k,"var",] = apply(Y[,"rho",],1,function(x)var(x,na.rm=TRUE))
  Rho[k,"bias",] = apply(Y[,"rho",],1,function(x)mean((x-rho0),na.rm=TRUE))
  Rho[k,"rmse",] = apply(Y[,"rho",],1,function(x)sqrt(mean((x-rho0)^2,na.rm=TRUE)))
  Rho[k,"meanSE",] = apply(Y[,"rhoSE",],1,function(x)mean(x,na.rm=TRUE))
  Rho[k,"ratioBias",] = apply(Y[,"rho",],1,function(x)sum(x>rho0,na.rm=TRUE))/apply(Y[,"rho",],1,function(x)sum(x<rho0,na.rm=TRUE))
  Rho[k,"propPosBias",] = apply(Y[,"rho",],1,function(x)sum(x>rho0,na.rm=TRUE))/apply(Y[,"rho",],1,function(x)sum(!is.na(x)))
  
  tau0 = rep(seq(from=-2,to=2,length.out = design$M[k]-1),1)
  jjd=grep(x = colnames(Y),pattern = "tauX")
  for(b in 1:B){X=Y[,jjd,b]; X[is.infinite(X)]=NA; Y[,jjd,b]=X}
  TauX1[k,"bias",] = apply(abs(apply(Y[,jjd,],c(1,2),function(x)mean(x,na.rm=TRUE)) - pracma::repmat(tau0,3,1)),1,mean)
  TauX1[k,"var",] = apply(apply(Y[,jjd,],c(1,2),function(x)var(x,na.rm=TRUE)),1,mean)
  TauX1[k,"rmse",] = apply(t(mapply(function(m)sqrt(apply((Y[m,jjd,]-tau0)^2,1,function(x)mean(x,na.rm=TRUE))),1:3)),1,mean)
  TauX1[k,"ratioBias",] = apply(t(mapply(function(m)apply(Y[m,jjd,]-tau0,1,function(x)sum(x>0,na.rm=TRUE)/sum(x<0,na.rm=TRUE)),1:3)),1,mean)
  TauX1[k,"propPosBias",] = apply(t(mapply(function(m)apply(Y[m,jjd,]-tau0,1,function(x)sum(x>0,na.rm=TRUE)/sum(!is.na(x))),1:3)),1,mean)
  
  tau0 = rep(seq(from=-2,to=2,length.out = design$M[k]-1),1)
  jjd=grep(x = colnames(Y),pattern = "tauY")
  for(b in 1:B){X=Y[,jjd,b]; X[is.infinite(X)]=NA; Y[,jjd,b]=X}
  TauX2[k,"bias",] = apply(abs(apply(Y[,jjd,],c(1,2),function(x)mean(x,na.rm=TRUE)) - pracma::repmat(tau0,3,1)),1,mean)
  TauX2[k,"var",] = apply(apply(Y[,jjd,],c(1,2),function(x)var(x,na.rm=TRUE)),1,mean)
  TauX2[k,"rmse",] = apply(t(mapply(function(m)sqrt(apply((Y[m,jjd,]-tau0)^2,1,function(x)mean(x,na.rm=TRUE))),1:3)),1,mean)
  TauX2[k,"ratioBias",] = apply(t(mapply(function(m)apply(Y[m,jjd,]-tau0,1,function(x)sum(x>0,na.rm=TRUE)/sum(x<0,na.rm=TRUE)),1:3)),1,mean)
  TauX2[k,"propPosBias",] = apply(t(mapply(function(m)apply(Y[m,jjd,]-tau0,1,function(x)sum(x>0,na.rm=TRUE)/sum(!is.na(x))),1:3)),1,mean)
  
  Z=apply(Y[,grep(x = colnames(Y),pattern = "tau"),],c(1,2),function(x)mean(x,na.rm=TRUE))
  Tau_fEM[k,"mean",1:((design$M[k]-1)*2)] = Z[1,];Tau_maxdef[k,"mean",1:((design$M[k]-1)*2)] = Z[2,];Tau_meandef[k,"mean",1:((design$M[k]-1)*2)] = Z[3,]

  Z=apply(Y[,grep(x = colnames(Y),pattern = "tau"),],c(1,2),function(x)var(x,na.rm=TRUE))
  Tau_fEM[k,"var",1:((design$M[k]-1)*2)] = Z[1,];Tau_maxdef[k,"var",1:((design$M[k]-1)*2)] = Z[2,];Tau_meandef[k,"var",1:((design$M[k]-1)*2)] = Z[3,]
}  
  
Y1=data.frame(cbind(pracma::repmat(as.matrix(design[iid,]),3,1),rbind(Rho[iid,,1],Rho[iid,,2],Rho[iid,,3])))
Y1$method=rep(c("fEM","maxdef","meandef"),each=length(iid));Y1$method=as.factor(Y1$method); names(Y1)[1:4]=colnames(design)
Y1$mse=Y1$bias^2+Y1$var

Y2=data.frame(cbind(pracma::repmat(as.matrix(design[iid,]),3,1),rbind(TauX1[iid,,1],TauX1[iid,,2],TauX1[iid,,3])))
Y2$method=rep(c("fEM","maxdef","meandef"),each=length(iid));Y2$method=as.factor(Y2$method); names(Y2)[1:4]=colnames(design)
Y2$mse=Y2$bias^2+Y2$var

Y3=data.frame(cbind(pracma::repmat(as.matrix(design[iid,]),3,1),rbind(TauX2[iid,,1],TauX2[iid,,2],TauX2[iid,,3])))
Y3$method=rep(c("fEM","maxdef","meandef"),each=length(iid));Y3$method=as.factor(Y3$method); names(Y3)[1:4]=colnames(design)
Y3$mse=Y3$bias^2+Y3$var

#save(Rho,TauX1,TauX2,Y1,Y2,Y3,Tau_fEM,Tau_maxdef,Tau_meandef,file = "sim_results.rda")


