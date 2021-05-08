trapezoidal_fn = function(x=NULL,lb=0,ub=1,m1=0.5,m2=0.5){
  #It generalizes: (i) Triangular case (m1=m2), (ii) Rectangular case (lb=m1,ub=m2)
  if(is.null(x)){x=seq(from=lb,to=ub,length.out=101)}
  y=rep(1e-99,length(x))
  y[x>=lb&x<m1] = (x[x>=lb&x<m1]-lb)/(m1-lb)
  y[x>=m1&x<=m2] = 1
  y[x>m2&x<=ub] = (ub-x[x>m2&x<=ub])/(ub-m2)
  return(y)
}


inclusion_degree = function(fn_y=NULL,fn_c=NULL,yinf=NULL,ysup=NULL){ 
  #it works for a trapezoidal/triangular category and a trapezoidal/triangular observation
  if(is.null(yinf)){yinf=min(min(fn_y),min(fn_c))}; if(is.null(ysup)){ysup=max(max(fn_y),max(fn_c))}
  if(fn_y[4]-fn_y[1]<1e-02){
    #non-fuzzy observation
    eps_yc = trapezoidal_fn(x=fn_y[1],lb=fn_c[1],m1=fn_c[2],m2=fn_c[3],ub=fn_c[4])
  }else{
    #observation
    card_y = integrate(function(x)trapezoidal_fn(x=x,lb=fn_y[1],m1=fn_y[2],m2=fn_y[3],ub=fn_y[4]),yinf,ysup,subdivisions=2000)$value
    #category
    card_c = integrate(function(x){apply(cbind(trapezoidal_fn(x=x,lb=fn_y[1],m1=fn_y[2],m2=fn_y[3],ub=fn_y[4]),
                                               trapezoidal_fn(x=x,lb=fn_c[1],m1=fn_c[2],m2=fn_c[3],ub=fn_c[4]))
                                         ,1,min)},yinf,ysup,subdivisions=2000)$value
    eps_yc = card_c/card_y
  }
  return(eps_yc)
}


inclusion_degree_vectorize = function(Ydata=NULL,C=NULL){
  I=NROW(Ydata); M=NROW(C)
  return(mapply(function(m)mapply(function(i)inclusion_degree(fn_y = Ydata[i,],fn_c = C[m,]),1:I),1:M))
}


FCounts = function(e=NULL,normalize=TRUE){
  n=length(e); K=0:n; fgc=flc=fec=matrix(0,n+1,1)
  
  x=mapply(function(a)sum(e>=a),e)
  fgc = mapply(function(k)max(0,e[x>=k],na.rm = TRUE),K)
  flc = 1-mapply(function(k)max(0,e[x>=k],na.rm=TRUE),K+1)
  fec = apply(cbind(fgc,flc),1,min)
  
  if(sum(fec<1e-04)>=(n+1)){fec[1]=1}; if(normalize){fec=fec/max(fec)}
  return(list(K=K,FGCount=fgc,FLCount=flc,FECount=fec))
}


FECount_joint = function(Y1=NULL,Y2=NULL,c1=NULL,c2=NULL,normalize=TRUE,verbose){ #it works given two categories only (i.e., a given cell of Fuzzy Freq Matrix)
  c1 = matrix(c1,nrow=1); c2 = matrix(c2,nrow=1)
  e1 = inclusion_degree_vectorize(Y1,c1)
  e2 = inclusion_degree_vectorize(Y2,c2)
  e = apply(cbind(e1,e2),1,min)
  fec_joint = FCounts(e,normalize)$FECount
  return(fec_joint)
}

plot_FuzzyFrequency_table = function(FFtab=NULL,M,plotx_xlim=TRUE,new.dev=TRUE,dir=1){
  xnames=colnames(FFtab); xsup=0:(NROW(FFtab)-1); H=NCOL(FFtab);
  if(new.dev){x11();if(dir==1){par(mfcol=c(M,M))}else{par(mfrow=c(M,M))}}else{if(dir==1){par(mfcol=c(M,M))}else{par(mfrow=c(M,M))}}
  for(h in 1:H){
    iid=FFtab[,h]>1e-03
    if(plotx_xlim){plot(xsup[iid],FFtab[iid,h],type="p",bty="n",xlab="",ylab="",main=xnames[h],lwd=1.25,xlim=c(min(xsup[iid])-5,max(xsup[iid])+5),ylim=c(0,1.05))}
    else{plot(xsup[iid],FFtab[iid,h],type="p",bty="n",xlab="",ylab="",main=xnames[h],lwd=1.25,xlim=c(-1,max(xsup)),ylim=c(0,1))}
    segments(x0=xsup[iid],x1=xsup[iid],y0=0,y1=FFtab[iid,h],lty=2,lwd=1.25) 
  }
}

FuzzyFrequency_table = function(Y1=NULL,Y2=NULL,C1=NULL,C2=NULL,normalize=TRUE,verbose=TRUE,plotx=FALSE,plotx_xlim=TRUE){
  n = NROW(Y1)
  M = NROW(C1); #both categories C1 and C2 should have the same length!
  Iid = expand.grid(1:M,1:M) #overall indices for the categories (pairs)
  H = NROW(Iid); FF = matrix(NA,n+1,H)
  
  if(verbose){pracma::tic()}
  for(h in 1:H){
    j = Iid[h,1]; k = Iid[h,2]
    FF[,h] = FECount_joint(Y1,Y2,C1[j,],C2[k,],normalize)
  }
  rownames(FF)=0:n; colnames(FF)=apply(Iid,1,function(x)paste(x,collapse=""))
  if(verbose){pracma::toc()}
  if(plotx){plot_FuzzyFrequency_table(FF,M,plotx_xlim = plotx_xlim)}
  
  return(round(FF,4))
}

FuzzyFrequency_list = function(Y=NULL,Cats=NULL,verbose=TRUE,normalize=TRUE){
  J=length(Y); D1=expand.grid(1:J,1:J); Iid = D1[D1[,1]<D1[,2],]
  FFlist = list()
  for(u in 1:NROW(Iid)){
    j=Iid[u,1]; k=Iid[u,2]
    if(verbose){cat(paste0("@ Computing pair: ",paste(Iid[u,],collapse=",")),"\n")}
    FFlist[[u]] = FuzzyFrequency_table(Y[[j]],Y[[k]],Cats[[j]],Cats[[k]],verbose=verbose,normalize=normalize)
  }
  names(FFlist) = as.character(apply(Iid,1,function(x)paste(x,collapse=",")))
  if(verbose){cat("@ Done.","\n")}
  return(list(FFtables=FFlist,Pairs=Iid))
}

loglik = function(rho,rc,cc,tab){
  P = matrix(0, (length(rc)-1),(length(cc)-1))
  R = matrix(c(1,rho,rho,1),2,2)
  for (i in 1:(length(rc)-1)){for(j in 1:(length(cc)-1)){P[i,j]=mvtnorm::pmvnorm(lower = c(rc[i], cc[j]), upper = c(rc[i+1], cc[j+1]), mean=rep(0,2),corr = R)}}
  P[P<=0] = NA
  #print(round(P,3))
  lP = log(P)
  lP[lP == -Inf] = NA; lP[lP == Inf] = NA
  loglikx = -sum(tab * lP,na.rm=TRUE) #loglikelhood of bivariate normal model
  
  return(loglikx)
}

loglik_stirling = function(rho,rc,cc,tab){
  P = matrix(0, (length(rc)-1),(length(cc)-1)); R = matrix(c(1,rho,rho,1),2,2)
  for (i in 1:(length(rc)-1)){for(j in 1:(length(cc)-1)){P[i,j]=mvtnorm::pmvnorm(lower = c(rc[i], cc[j]), upper = c(rc[i+1], cc[j+1]), mean=rep(0,2),corr = R)}}
  P[P<=0] = NA; lP = log(P); lP[lP == -Inf] = NA; lP[lP == Inf] = NA
  tP = log(tab); tP[tP == -Inf] = NA; tP[tP == Inf] = NA
  loglikx = sum(tab,na.rm=TRUE) - 0.5*sum(tab*tP ,na.rm=TRUE) - sum(tab * lP,na.rm=TRUE) 
  return(loglikx)
}


approx_lambda_rc = function(n=NULL,Tau=NULL,rho=NULL){
  ## Using results from Shiina, K., Ueda, T., & Kubo, S. (2017, July). Polychoric correlations for ordered categories using the EM algorithm. In The Annual Meeting of the Psychometric Society (pp. 247-259). Springer
  if(is.null(dim(Tau))){Tau=matrix(Tau,2)}
  rc = c(-Inf,Tau[1,],Inf); cc = c(-Inf,Tau[2,],Inf)
  Lambda=matrix(NA,(length(rc)-1),(length(cc)-1))
  R = matrix(c(1,rho,rho,1),2,2)
  if(!pracma::isposdef(R)){R = make.positive.definite(R)}
  for(r in 1:(length(rc)-1)){
    for(c in 1:(length(cc)-1)){
      Lambda[r,c]=tmvtnorm::ptmvnorm(lowerx = c(rc[r], cc[c]), upperx = c(rc[r+1], cc[c+1]), mean=rep(0,2),sigma = R)*n
      #Lambda[r,c]=mvtnorm::pmvnorm(lower = c(rc[r], cc[c]), upper = c(rc[r+1], cc[c+1]), corr = R,mean = rep(0,2))*n
    }
  }
  Lambda[Lambda<0]=0
  return(Lambda)
}

compute_nrc = function(FFtab=NULL,tab=NULL){
  nsup=0:(NROW(FFtab)-1)
  tabx=as.vector(tab);
  #Px=sapply(tabx,function(x)dpois(nsup,x))
  Px=sapply(tabx,function(x)dbinom(nsup,sum(tabx),x/sum(tabx)))
  Pn_tilde = matrix(mapply(function(j)sum(FFtab[,j]*Px[,j]),1:length(tabx)),NROW(tab),NCOL(tab)); pn_tilde = as.vector(Pn_tilde)
  N_est = matrix(mapply(function(j)sum(nsup*((FFtab[,j]*Px[,j])/pn_tilde[j])),1:length(tabx)),NROW(tab),NCOL(tab)); N_est[is.na(N_est)]=0
  return(N_est)
}

compute_nrc_var = function(FFtab=NULL,tab=NULL){
  nsup=0:(NROW(FFtab)-1)
  tabx=as.vector(tab);
  #Px=sapply(tabx,function(x)dpois(nsup,x))
  Px=sapply(tabx,function(x)dbinom(nsup,sum(tabx),x/sum(tabx)))
  Pn_tilde = matrix(mapply(function(j)sum(FFtab[,j]*Px[,j]),1:length(tabx)),NROW(tab),NCOL(tab)); pn_tilde = as.vector(Pn_tilde) #denominatore conditional density
  N_est = matrix(mapply(function(j)sum(nsup*((FFtab[,j]*Px[,j])/pn_tilde[j])),1:length(tabx)),NROW(tab),NCOL(tab)); N_est[is.na(N_est)]=0; n_est = as.vector(N_est) #expectations
  V_est = matrix(mapply(function(j)sum((nsup-n_est[j])^2*((FFtab[,j]*Px[,j])/pn_tilde[j])),1:length(tabx)),NROW(tab),NCOL(tab)); V_est[is.na(V_est)]=0; #variances
  return(V_est)
}


make.positive.definite = function(m, tol){
  if (!is.matrix(m)) m = as.matrix(m)
  
  d = dim(m)[1] 
  if ( dim(m)[2] != d ) stop("Input matrix is not square!")
  
  es = eigen(m, symmetric=TRUE)
  esv = es$values
  
  if (missing(tol))
    tol = d*max(abs(esv))*.Machine$double.eps 
  delta =  2*tol # factor to is just to make sure the resulting
  # matrix passes all numerical tests of positive definiteness
  
  tau = pmax(0, delta - esv)
  dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)    
  
  #print(max(DA))
  #print(esv[1]/delta)
  
  return( m +  dm )
}

gen_counts = function(n=NULL,Tau=NULL,rho,seed=NULL){ #generate counts according to a gaussian latent model (bivariate)
  if (!is.null(seed)){set.seed(seed)}
  tab = matrix(0,NCOL(Tau)+1,NCOL(Tau)+1)
  R=diag(2);R[1,2]=R[2,1]=rho
  X = MASS::mvrnorm(n = n,mu = rep(0,NROW(Tau)),Sigma = R); Y=X
  Y = mapply(function(j){as.numeric(cut(x = X[,j],breaks = c(-Inf,Tau[j,],Inf),right=TRUE,labels = seq(1,NCOL(Tau)+1)))},1:NCOL(X))
  tab = tableFast(Y[,1],Y[,2],min(Y[,1]),max(Y[,1]),min(Y[,2]),max(Y[,2]))
  return(list(X=X,Y=Y,tab=tab))
}

bigammaFuzzify_nrc = function(Ff=NULL,mx=1,sx=1.5,plotx=FALSE,only.spread=FALSE,normalize=TRUE){ #fuzzification of a crisp count according to a Gamma-based model (not centered around crisp nrc)
  m0=as.numeric(which.max(Ff)-1)
  xsup=0:(NROW(Ff)-1)
  
  # modeling alpha (model for the spread of nrc)
  b = (mx + sqrt(mx^2 + 4*sx^2)) / (2*sx^2)
  a = 1 + mx*b
  s0=max(1,extraDistr::rdgamma(1,a,b)) #via discrete gamma dist
  
  # modeling beta (model for the shift of mode for nrc)
  if(only.spread==FALSE){
    b = (m0 + sqrt(m0^2 + 4*s0^2)) / (2*s0^2)
    a = 1 + m0*b
    my=extraDistr::rdgamma(1,a,b) #via discrete gamma dist
  }else{my=m0}
  
  # modeling n_rc fuzzy
  b = (my + sqrt(my^2 + 4*s0^2)) / (2*s0^2)
  a = 1 + my*b
  mux=extraDistr::ddgamma(xsup,a,b); if(normalize){mux=mux/max(mux)}
  
  if(plotx){
    x11();
    iid=mux>1e-03
    plot(xsup[iid],Ff[iid],type="p",bty="n",xlab="",ylab="",main="",lwd=1.25,xlim=c(min(xsup[iid])-5,max(xsup[iid])+5),ylim=c(0,1))
    segments(x0=xsup[iid],x1=xsup[iid],y0=0,y1=Ff[iid],lty=2,lwd=1.25) 
    points(xsup,mux,type="p",lwd=1,col=2); segments(x0=xsup,x1=xsup,y0=0,y1=mux,lty=2,lwd=1,col=2) 
  }
  return(mux)
}

diff_loglik = function(tab=NULL,rho=NULL,Tau=NULL){
  cc=c(-Inf,Tau[1,],Inf); rc=c(-Inf,Tau[2,],Inf)
  len1 = length(rc)-1; len2 = length(cc)-1
  P = dP = matrix(0,len1,len2)
  
  # Compute P
  R = matrix(c(1,rho,rho,1),2,2)
  for (i in 1:len1){for(j in 1:len2){P[i,j]=mvtnorm::pmvnorm(lower = c(rc[i], cc[j]), upper = c(rc[i+1], cc[j+1]), corr = R)}}
  P[P<=0] = NA
  
  # Compute dP
  for(i in 2:(len1+1)){
    for(j in 2:(len2+1)){
      p1=p2=p3=p4=0
      if(i<len1&j<len2){p1 = mvtnorm::dmvnorm(c(rc[i],cc[j]),sigma=R)}
      if(i>2&j<len2){p2 = mvtnorm::dmvnorm(c(rc[i-1],cc[j]),sigma=R)}
      if(i<len1&j>2){p3 = mvtnorm::dmvnorm(c(rc[i],cc[j-1]),sigma=R)}
      if(i>2&j>2){p4 = mvtnorm::dmvnorm(c(rc[i-1],cc[j-1]),sigma=R)}
      dP[i-1,j-1]=p1-p2-p3+p4
    }
  }
  
  # Compute d_Loglik/d_rho
  dloglik = sum(tab/P*dP,na.rm=TRUE)
  
  return(dloglik)
}

compute_SEs = function(tab=NULL,rho=NULL,Tau=NULL){
  R = NROW(tab); C = NCOL(tab); I_empirical = 0
  for(r in 1:R){
    for(c in 1:C){
      sc = matrix(diff_loglik(matrix(tab[r,c],R,C),rho,Tau))
      I_empirical = I_empirical + sc%*%sc
    }
  }
  ses = sqrt(1/I_empirical)
  
  return(ses)
}

startingValues_EM = function(method=c("max.defuzz","mean.defuzz","random"),Tau0=NULL,rho0=NULL,FFtab=NULL,R=NULL,C=NULL){
  methodx = match.arg(method)
  out = switch(methodx,
         max.defuzz = estimPars_maxdefuzz(FFtab,R,C),
         mean.defuzz = estimPars_meandefuzz(FFtab,R,C),
         random = list(rho=rnorm(1,0,0.15),
                       Tau=matrix(rep(seq(from=-2,to=2,length.out = R-1),2),nrow = 2,byrow = TRUE) +
                         matrix(rnorm((R-1)*2,sd=0.1),2)) #R=C for Tau
                      
           )
  
  if(is.null(Tau0)&is.null(rho0)){ #both Tau and rho need to be estimated
    Tau0 = out$Tau
    rho0 = out$rho
  }
  if(is.null(rho0)&!is.null(Tau0)){ #only rho needs to be computed
    rho0 = out$rho
  }
  if(!is.null(rho0)&is.null(Tau0)){ #only rho needs to be computed
    Tau0 = out$Tau
  }
  
  return(list(rho0=rho0,Tau0=Tau0))
}


loglik_EM = function(rho,rc,cc,tab,tab_var){
  P = matrix(0, (length(rc)-1),(length(cc)-1)); R = matrix(c(1,rho,rho,1),2,2)
  for (i in 1:(length(rc)-1)){for(j in 1:(length(cc)-1)){P[i,j]=mvtnorm::pmvnorm(lower = c(rc[i], cc[j]), upper = c(rc[i+1], cc[j+1]), mean=rep(0,2),corr = R)}}
  P[P<=0] = NA; lP = log(P); lP[lP == -Inf] = NA; lP[lP == Inf] = NA
  tP = log(tab); tP[tP == -Inf] = NA; tP[tP == Inf] = NA
  loglikx = -sum(tab * lP,na.rm=TRUE) - 0.5*sum( tab*tP + (tab_var/(2*tab)) ,na.rm=TRUE) + sum(tab) - sum(tab)*0.5*log(2*pi)
  return(loglikx)
}


#it works for a given pair (j,k)
estimPars_EM = function(K=500,R=NULL,C=NULL,FFtab=NULL,Tau0=NULL,rho0=NULL,only.rho=FALSE,use.optim=FALSE,eps=1e-7,verbose=TRUE,starting.method=c("max.defuzz","mean.defuzz","random")){
  if(is.null(C)){C=R}
  n=NROW(FFtab)-1
  
  ## Starting values (using max-defuzzification)
  methodx = match.arg(starting.method)
  out = startingValues_EM(method = methodx,Tau0 = Tau0,rho0 = rho0,FFtab = FFtab,R = R,C = C)
  
  ## EM routine
  Taux = array(NA,c(2,R-1,K)); Rhox = matrix(NA,K,1); loglikx = matrix(NA,K,1)
  conv=1; k=2; 
  Taux[,,1] = out$Tau; Rhox[1,1] = out$rho; loglikx[1] = -100000
  
  while(conv>eps){
    if(verbose){message(paste0("@ Iter no.: ",k))}
    ## Compute expectations
    N_star = approx_lambda_rc(n,Taux[,,k-1],Rhox[k-1,1])
    N_est = compute_nrc(FFtab,N_star)
    V_est = compute_nrc_var(FFtab,N_star)
    
    #nn_est = apply(rmultinom(n=10000,size = n,prob = as.vector(round(N_est/sum(N_est),3))),1,mean)
    #N_est = matrix(nn_est,R,C)
    
    ## Maximization 
    tab = N_est/sum(N_est)
    if(!only.rho){
      csum = apply(tab,2,sum); rsum = apply(tab,1,sum)
      Tau_k = rbind(qnorm(cumsum(csum)[-length(csum)]),qnorm(cumsum(rsum)[-length(rsum)]))
    }else{Tau_k = Tau0}
    cc =  c(-Inf,Tau_k[1,],Inf); rc =  c(-Inf,Tau_k[2,],Inf)
    #if(use.optim==FALSE){rho_k = uniroot(f = function(x)diff_loglik(tab,x,Tau_k),lower = -0.99,upper = 0.99)$root
    #}else{rho_k = optim(par = 0.1,fn = function(x)loglik_EM(x,rc,cc,N_est,V_est),method = "L-BFGS-B",lower = -0.99,upper = 0.99)$par}
    rho_k = optim(par = 0.1,fn = function(x)loglik_EM(x,rc,cc,N_est,V_est),method = "L-BFGS-B",lower = -0.99,upper = 0.99)$par
    Taux[,,k] = Tau_k; Rhox[k,1] = rho_k
    if(verbose){message(paste0("Current params: ",paste(paste(round(Taux[1,,k],4),collapse = ", "),paste(round(Taux[2,,k],4),collapse = ", "),round(Rhox[k,1],4),sep = " | ")))}
    ## Compute loglik and check for convergence
    #loglikx[k] = loglik(rho_k,rc,cc,tab)
    loglikx[k] = loglik_EM(rho_k,rc,cc,N_est,V_est)  
    conv = abs(loglikx[k]-loglikx[k-1])/abs(loglikx[k])
    if(verbose){message(paste0(paste0("Current loglik: ", round(loglikx[k],4))," (",round(conv,19),") ","\n"))}
    if(k>=K){break};
    k=k+1
  }
  
  ## Prepare output
  N_star = approx_lambda_rc(n,Taux[,,k-1],Rhox[k-1,1]); N_est = compute_nrc(FFtab,N_star)
  ses = compute_SEs(tab = N_est/sum(N_est), rho = Rhox[k-1,1], Tau = Taux[,,k-1])
  outlist = list(init_par=list(method=methodx,rho=out$rho,Tau=out$Tau),Theta=list(Tau=Taux[,,1:(k-1)],rho=Rhox[1:(k-1)]),loglikel=loglikx[1:(k-1)],rho_EM=list(rho=Rhox[k-1,1],se=ses),Tau_EM=Taux[,,k-1],conv=ifelse(conv<=eps,1,-1),N=list(N_est=N_est,N_star=N_star))
  
  return(invisible(outlist))
}

estimPars_maxdefuzz = function(FFtab=NULL,R=NULL,C=NULL){
  if(is.null(C)){C=R}
  tab_max = matrix(apply(FFtab,2,function(x)which.max(x)-1),R,C)
  tab_max=tab_max/sum(tab_max)
  csum = apply(tab_max,2,sum); rsum = apply(tab_max,1,sum)
  Taux = rbind(qnorm(cumsum(csum)[-length(csum)]),qnorm(cumsum(rsum)[-length(rsum)]))
  cc =  c(-Inf,Taux[1,],Inf);rc =  c(-Inf,Taux[2,],Inf)
  rhox = optim(par = 0.1,fn = function(x)loglik(x,rc,cc,tab_max),method = "L-BFGS-B",lower = -0.99,upper = 0.99)$par
  return(list(rho=rhox,Tau=Taux))
}

estimPars_meandefuzz = function(FFtab=NULL,R=NULL,C=NULL){
  if(is.null(C)){C=R}
  tab_mean = matrix(apply(FFtab,2,function(x)sum(x*seq(from=0,to=length(x)-1))/sum(x)),R,C)
  tab_mean=tab_mean/sum(tab_mean)
  csum = apply(tab_mean,2,sum); rsum = apply(tab_mean,1,sum)
  Taux = rbind(qnorm(cumsum(csum)[-length(csum)]),qnorm(cumsum(rsum)[-length(rsum)]))
  cc =  c(-Inf,Taux[1,],Inf);rc =  c(-Inf,Taux[2,],Inf)
  rhox = optim(par = 0.1,fn = function(x)loglik(x,rc,cc,tab_mean),method = "L-BFGS-B",lower = -0.99,upper = 0.99)$par
  return(list(rho=rhox,Tau=Taux))
}

add_legend = function(...) {
  #From: https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}



