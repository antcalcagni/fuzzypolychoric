##### Set environment #####
rm(list=ls());graphics.off()
setwd("simulation_study/")
source("utilities.R"); require(doParallel)
load("design.rda")
ncores=16
cl = parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl,cores=ncores)
combAbind = abind::abind


estimparsfun = function(N_tilde,Tau0,rho0,M){
  out1=out2=out3=NULL
  EstimPars=matrix(NA,3,(2+(M-1)+(M-1))) #rows: methods; cols: rho|SE_rho|tau1|tau2
  
  tryCatch(expr={out1=estimPars_EM(R=M,C=M,FFtab=N_tilde,use.optim=TRUE,K=5000,eps=1e-6,Tau0=Tau0,rho0=rho0,verbose=FALSE,starting.method="random")},error = function(e){print(e)}) #EM
  tryCatch(expr={out2=estimPars_maxdefuzz(FFtab=N_tilde,R=M,C=M)},error = function(e){print(e)}) #max-defuzz
  tryCatch(expr={out3=estimPars_meandefuzz(FFtab=N_tilde,R=M,C=M)},error = function(e){print(e)}) #max-defuzz
  
  if(!is.null(out1)){if(out1$conv==1){EstimPars[1,] = as.numeric(c(unlist(out1$rho_EM),as.vector(t(out1$Tau_EM))))}}
  if(!is.null(out2)){EstimPars[2,] = as.numeric(c(out2$rho,NA,as.vector(t(out2$Tau))))}
  if(!is.null(out3)){EstimPars[3,] = as.numeric(c(out3$rho,NA,as.vector(t(out3$Tau))))}
  #rownames(EstimPars)=c("fEM","maxdef","meandef")
  #colnames(EstimPars)=c("rho","rhoSE",paste0("tauX",1:(M-1)),paste0("tauY",1:(M-1)))
  
  return(EstimPars)
}


##### Run algorithms on simulated data #####
B = 5000
chunk.size = floor(B/ncores)
n_resid = rep(0,ncores)
if((ncores*chunk.size)<B){n_resid[ncores] = B - ncores*chunk.size}

for(k in 1:NROW(design)){
  load(paste0("gendata_",k,".rda"))
  
  cat(paste("\n Design cell n.: ",k),sep="")
  cat("\n > Running algorithms")
  
  Tau0 = matrix(rep(seq(from=-2,to=2,length.out = Ycurr$M-1),2),nrow = 2,byrow = TRUE) + matrix(rnorm(2*(Ycurr$M-1),sd=0.1),2) 
  rho0 = Ycurr$rho + rnorm(1,0,0.05)
  
  funToExport=dataToExport=NULL
  Y <- foreach(h=1:ncores, .combine = "combAbind", .export = c(funToExport,dataToExport)) %dopar% { 
    res = array(NA, dim = c(3,(2+(Ycurr$M-1)*2),chunk.size+n_resid[h]),dimnames = list(c("fEM","maxdef","meandef"),c("rho","rhoSE",paste0("tauX",1:(Ycurr$M-1)),paste0("tauY",1:(Ycurr$M-1))),NULL)) #..+2 because nsup=0:nk
    iid_core = ((h-1)*chunk.size+1):((h*chunk.size)+n_resid[h])
    for(u in iid_core){
      res[,,u-(h-1)*chunk.size] = mapply(function(u)estimparsfun(N_tilde = Ycurr$N_fuzzy[,,u],Tau0 = Tau0,rho0 = rho0,M = Ycurr$M),u)
    }
    res #it is needed to populate Y along the third-dimension
  }
  cat("\n > Saving current data")
  save(Y,file=paste("estimdata_",k,".rda",sep=""))
  }

doParallel::stopImplicitCluster(); parallel::stopCluster(cl)
cat("\n Done.")





