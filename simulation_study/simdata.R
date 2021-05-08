##### Set environment #####
rm(list=ls());graphics.off()
setwd("simulation_study/")
source("utilities.R"); require(doParallel)
ncores=4
cl = parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl,cores=ncores)
combAbind = abind::abind
set.seed(28032020)

gendatafun = function(n,rho,Tau,M,sx){
  #generate crisp freq table in fuzzy freq form
  tab = approx_lambda_rc(n,Tau,rho)
  N = matrix(0,n+1,M*M); colnames(N)=apply(expand.grid(1:M,1:M),1,function(x)paste(x,collapse="")); rownames(N)=0:n;
  fftab = round(as.vector(tab)); if(sum(fftab)<n){fftab[which.max(fftab)]=fftab[which.max(fftab)]+(n-sum(fftab))}
  for(j in 1:NCOL(N)){N[fftab[j]+1,j]=1}
  
  #fuzzification
  N_tilde = apply(N,2,function(x)bigammaFuzzify_nrc(x,sx=sx)); rownames(N_tilde)=0:n;
  
  #arrange output
  Y = rbind(N_tilde,N)
  return(Y)
}


##### Simulation design #####
n = c(150,250,500,1000) #sample sizes
rho = c(0.15,0.5,0.85) #true rho
m = c(4,6) #number of categories R/C
sx = 0.25 #degree of fuzzification
design = expand.grid(n,m,rho,sx); colnames(design)=c("n","M","rho","sx")
save(design,file="design.rda")

##### Generate data for each cell of design #####
B=5000
gendata = list()

chunk.size = floor(B/ncores)
n_resid = rep(0,ncores)
if((ncores*chunk.size)<B){n_resid[ncores] = B - ncores*chunk.size}

for(k in 1:NROW(design)){
  print(paste0("Processing design no.: ",k))
  
  #get current design 
  nk=design$n[k]; Mk=design$M[k]; sxk=design$sx[k]; rhok=design$rho[k]
  
  #simulate B samples for the k-th design
  Tau = matrix(rep(seq(from=-2,to=2,length.out = Mk-1),2),nrow = 2,byrow = TRUE)
  
  funToExport=dataToExport=NULL
  Y <- foreach(h=1:ncores, .combine = "combAbind", .export = c(funToExport,dataToExport)) %dopar% { 
    res = array(NA, dim = c((nk*2+2),Mk*Mk,chunk.size+n_resid[h]),dimnames = list(NULL,NULL,NULL)) #..+2 because nsup=0:nk
    iid_core = ((h-1)*chunk.size+1):((h*chunk.size)+n_resid[h])
    for(u in iid_core){
      res[,,u-(h-1)*chunk.size] = mapply(function(u)gendatafun(nk,rhok,Tau,Mk,sxk),u)
    }
    res #it is needed to populate Y along the third-dimension
  }
  
  #save generated data
  gendata[[k]] = list(n=nk,M=Mk,rho=rhok,sx=sxk,N_fuzzy=Y[1:(nk+1),,],N_crisp=Y[(nk+2):(nk*2+2),,])
  Ycurr=gendata[[k]];save(Ycurr,file=paste0("gendata_",k,".rda"))
  
}

#stop clusters
doParallel::stopImplicitCluster() 
parallel::stopCluster(cl)
