rm(list=ls());graphics.off()
source("utilities.R")

##### Figure 1 #####
set.seed(20210408)
I=25

# Trapezoidal categories (R=C=2)
C1 = C2 = matrix(rbind(c(1,1,3,5),
                       c(3.01,5.01,6,7)),ncol = 4)

# Fuzzy observations (two samples)
X1 = round(matrix(runif(I*3,1,7),ncol = 3),3); X1 = t(apply(X1,1,sort))
X1 = matrix(cbind(X1[,1],X1[,2],X1[,2],X1[,3]),nrow = I) #trapezoidal form

#tikzDevice::tikz(file='draft_paper/fig1.tex',width=8.5,height=4.5)
layout(mat = matrix(c(1,3,4,2,5,6), 2, 3, byrow = TRUE))
glbls=c("G1","G2");xsup=seq(from=1,to=7,length.out=1001)

# fuzzy obs (panel A)
C=C1;j=1;fx=trapezoidal_fn(x=xsup,lb = C[j,1],m1 = C[j,2],m2 = C[j,3],ub = C[j,4]); #fx=fx/(sum(fx)*2)
plot(xsup,fx,axes=FALSE,lwd=1.5,bty="n",ylab="",xlab="",type="l",cex.main=2,xlim=c(0.5,7.5),adj=0,ylim=c(0,1.15),col="deepskyblue4",lty=1,main="(A)"); text(mean(C[j,2:3])-0.5,1.1,glbls[j],col = "deepskyblue4",cex=1.25,cex.axis=1.25)
axis(side = 1, at=c(1,4,7),labels = c(1,4,7));axis(side = 2, at=c(0,1),labels = c(0,1));
for(j in 2:NROW(C)){fx=trapezoidal_fn(x=xsup,lb = C[j,1],m1 = C[j,2],m2 = C[j,3],ub = C[j,4]); points(xsup,fx,col="darkorange4",lwd=1.5,lty=1,type="l");text(mean(C[j,2:3])-0.5,1.1,glbls[j],col = "darkorange4",cex=1.25,cex.axis=1.25)}
D=X1; i=1; points(xsup,trapezoidal_fn(x=xsup,lb = D[i,1],m1 = D[i,2],m2 = D[i,3],ub = D[i,4]),col="gray30",lwd=1.25,lty=2,type="l")
for(i in 2:I){points(xsup,trapezoidal_fn(x=xsup,lb = D[i,1],m1 = D[i,2],m2 = D[i,3],ub = D[i,4]),col="gray30",lwd=1.25,lty=2,type="l")}

# crisp obs (panel B)
C=C1;j=1;fx=trapezoidal_fn(x=xsup,lb = C[j,1],m1 = C[j,2],m2 = C[j,3],ub = C[j,4]); #fx=fx/(sum(fx)*2)
plot(xsup,fx,axes=FALSE,lwd=1.5,bty="n",ylab="",xlab="",type="l",cex.main=2,xlim=c(0.5,7.5),ylim=c(0,1.15),adj=0,col="deepskyblue4",lty=1,main="(B)"); text(mean(C[j,2:3])-0.5,1.1,glbls[j],col = "deepskyblue4",cex=1.25,cex.axis=1.25)
axis(side = 1, at=c(1,4,7),labels = c(1,4,7));axis(side = 2, at=c(0,1),labels = c(0,1));
for(j in 2:NROW(C)){fx=trapezoidal_fn(x=xsup,lb = C[j,1],m1 = C[j,2],m2 = C[j,3],ub = C[j,4]); points(xsup,fx,col="darkorange4",lwd=1.5,lty=1,type="l");text(mean(C[j,2:3])-0.5,1.1,glbls[j],col = "darkorange4",cex=1.25,cex.axis=1.25)}
D=X1; for(i in 1:I){segments(x0=D[i,1],y0=0,x1=D[i,1],y1=1,col="gray30",lwd=1.25,lty=2)}

# fuzzy counts (panel C)
FFtab = FuzzyFrequency_table(X1,X1,C1,C2)
xsup=0:(NROW(FFtab)-1);
iid=FFtab[,1]>1e-03; plot(xsup[iid],FFtab[iid,1],type="p",bty="n",xlab="",ylab="",main=glbls[1],col.main="deepskyblue4",cex.main=1.15,lwd=1.25,xlim=c(-1,max(xsup)),ylim=c(0,1),yaxt = "n",col="deepskyblue4");segments(x0=xsup[iid],x1=xsup[iid],y0=0,y1=FFtab[iid,1],lty=2,lwd=1.25,col="deepskyblue4"); axis(side = 2, at=c(0,1),labels = c(0,1)) 
iid=FFtab[,4]>1e-03; plot(xsup[iid],FFtab[iid,4],type="p",bty="n",xlab="",ylab="",main=glbls[2],col.main="darkorange4",cex.main=1.15,lwd=1.25,xlim=c(-1,max(xsup)),ylim=c(0,1),yaxt = "n",col="darkorange4");segments(x0=xsup[iid],x1=xsup[iid],y0=0,y1=FFtab[iid,4],lty=2,lwd=1.25,col="darkorange4"); axis(side = 2, at=c(0,1),labels = c(0,1))

# fuzzy counts (panel D)
FFtab = FuzzyFrequency_table(pracma::repmat(matrix(X1[,1]),1,4),pracma::repmat(matrix(X1[,1]),1,4),C1,C2)
xsup=0:(NROW(FFtab)-1);
iid=FFtab[,1]>1e-05; plot(xsup[iid],FFtab[iid,1],type="p",bty="n",xlab="",ylab="",xlim=c(14,27),main=paste0("",glbls[1]),col.main="deepskyblue4",cex.main=1.15,lwd=1.25,ylim=c(0,1),yaxt="n",col="deepskyblue4");segments(x0=xsup[iid],x1=xsup[iid],y0=0,y1=FFtab[iid,1],lty=2,lwd=1.25,col="deepskyblue4");axis(side = 2, at=c(0,1),labels = c(0,1));
iid=FFtab[,4]>1e-03; plot(xsup[iid],FFtab[iid,4],type="p",bty="n",xlab="",ylab="",xlim=c(0,10),main=paste0("",glbls[2]),col.main="darkorange4",cex.main=1.15,lwd=1.25,ylim=c(0,1),yaxt="n",col="darkorange4");segments(x0=xsup[iid],x1=xsup[iid],y0=0,y1=FFtab[iid,4],lty=2,lwd=1.25,col="darkorange4");axis(side = 2, at=c(0,1),labels = c(0,1))
#dev.off()


##### Simulation study #####
load("simulation_study/sim_results.rda")
load("simulation_study/design.rda")

## Tables 2a,2b
X1 = cbind(Y1$n[Y1$M==4 & Y1$method=="fEM"],
           Y1$bias[Y1$M==4 & Y1$method=="fEM"],Y1$rmse[Y1$M==4 & Y1$method=="fEM"],
           Y1$bias[Y1$M==4 & Y1$method=="maxdef"],Y1$rmse[Y1$M==4 & Y1$method=="maxdef"],
           Y1$bias[Y1$M==4 & Y1$method=="meandef"],Y1$rmse[Y1$M==4 & Y1$method=="meandef"])

X2 = cbind(Y1$n[Y1$M==4 & Y1$method=="fEM"],
           Y1$bias[Y1$M==6 & Y1$method=="fEM"],Y1$rmse[Y1$M==6 & Y1$method=="fEM"],
           Y1$bias[Y1$M==6 & Y1$method=="maxdef"],Y1$rmse[Y1$M==6 & Y1$method=="maxdef"],
           Y1$bias[Y1$M==6 & Y1$method=="meandef"],Y1$rmse[Y1$M==6 & Y1$method=="meandef"])

Xtab = X1 #2a
Xtab_tex = xtable::xtable(Xtab,digits = c(0,0,rep(5,6)))
attributes(Xtab_tex)$caption = "Simulation study: "
attributes(Xtab_tex)$label = "tab2a"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
Xtab_tex$`1` = paste0("$I=",Xtab_tex$`1`,"$")
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x},include.rownames = FALSE)

Xtab = X2 #2b
Xtab_tex = xtable::xtable(Xtab,digits = c(0,0,rep(5,6)))
attributes(Xtab_tex)$caption = "Simulation study: "
attributes(Xtab_tex)$label = "tab2b"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
Xtab_tex$`1` = paste0("$I=",Xtab_tex$`1`,"$")
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x},include.rownames = FALSE)

## Tables 2c,2d
X3 = cbind(Y2$n[ Y2$M==4 & Y2$method=="fEM" ],
           Y2$bias[ Y2$M==4 & Y2$method=="fEM" ],Y2$rmse[ Y2$M==4 & Y2$method=="fEM" ],
           Y2$bias[ Y2$M==4 & Y2$method=="maxdef" ],Y2$rmse[ Y2$M==4 & Y2$method=="maxdef" ],
           Y2$bias[ Y2$M==4 & Y2$method=="meandef" ],Y2$rmse[ Y2$M==4 & Y2$method=="meandef" ])

Xtab = X3 #2c
Xtab_tex = xtable::xtable(Xtab,digits = c(0,0,rep(5,6)))
attributes(Xtab_tex)$caption = "Simulation study: "
attributes(Xtab_tex)$label = "tab2c"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
Xtab_tex$`1` = paste0("$I=",Xtab_tex$`1`,"$")
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x},include.rownames = FALSE)


X4 = cbind(Y2$n[Y2$n>50 & Y2$M==6 & Y2$method=="fEM" & Y2$sx==0.25],
           Y2$bias[Y2$n>50 & Y2$M==6 & Y2$method=="fEM" & Y2$sx==0.25],Y2$rmse[Y2$n>50 & Y2$M==6 & Y2$method=="fEM" & Y2$sx==0.25],
           Y2$bias[Y2$n>50 & Y2$M==6 & Y2$method=="maxdef" & Y2$sx==0.25],Y2$rmse[Y2$n>50 & Y2$M==6 & Y2$method=="maxdef" & Y2$sx==0.25],
           Y2$bias[Y2$n>50 & Y2$M==6 & Y2$method=="meandef" & Y2$sx==0.25],Y2$rmse[Y2$n>50 & Y2$M==6 & Y2$method=="meandef" & Y2$sx==0.25])

Xtab = X4 #2d
Xtab_tex = xtable::xtable(Xtab,digits = c(0,0,rep(5,6)))
attributes(Xtab_tex)$caption = "Simulation study: "
attributes(Xtab_tex)$label = "tab2d"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
Xtab_tex$`1` = paste0("$I=",Xtab_tex$`1`,"$")
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x},include.rownames = FALSE)

## Tables S3
iid=which(design$M==4)
X1=cbind(design$n[iid],Tau_fEM[iid,1,1:3])
iid=which(design$M==6)
X1=cbind(X1,Tau_fEM[iid,1,1:5])

Xtab =  X1#S3
Xtab_tex = xtable::xtable(Xtab,digits = c(0,0,rep(5,8)))
attributes(Xtab_tex)$caption = "Simulation study: "
attributes(Xtab_tex)$label = "tabS3"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
Xtab_tex$`1` = paste0("$I=",Xtab_tex$`1`,"$")
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x},include.rownames = FALSE)

## Tables S4
iid=which(design$M==4)
X1=cbind(design$n[iid],Tau_fEM[iid,2,1:3])
iid=which(design$M==6)
X1=cbind(X1,Tau_fEM[iid,2,1:5])

Xtab =  X1#S4
Xtab_tex = xtable::xtable(Xtab,digits = c(0,0,rep(5,8)))
attributes(Xtab_tex)$caption = "Simulation study: "
attributes(Xtab_tex)$label = "tabS4"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
Xtab_tex$`1` = paste0("$I=",Xtab_tex$`1`,"$")
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x},include.rownames = FALSE)

## Figure S1
iid=which(design$M==4)
cols=c("chocolate4","darkorange3","darkolivegreen")
xsup=c(150,250,500,1000); rhox=c(0.15,0.50,0.85); pm=0.30
#tikzDevice::tikz(file='draft_paper/figS1.tex',width=7.5,height=6.5)
x11();
par(mfcol=c(3,2),mai=c(0.7, 0.8, 0.3, 0.15))
for(h in 1:2){ #over tau
  for(j in 1:3){ #over rho
    if(h==1&j==1){xmain=latex2exp::TeX("$\\tau_1")}else{if(h==1&j>1){xmain=""}}
    if(h==2&j==1){xmain=latex2exp::TeX("$\\tau_3")}else{if(h==2&j>1){xmain=""}}
    if(h==1){ymain=paste0("p = ",rhox[j])}else{ymain=""}
    
    # drawing fEM
    ym=t(mapply(function(k)Tau_fEM[k,c("mean"),1:((design$M[k]-1))][c(1,(design$M[k]-1))],iid))-matrix(c(-2,2),12,2,byrow = TRUE)
    plot(xsup,ym[design$rho[iid]==rhox[j],h],bty="n",ylab=ymain,xlab="",ylim=c(-pm,+pm),xlim=c(50,1050),type="b",xaxt="none",col=cols[3],lwd=1.5,main=xmain,cex.main=2,cex.lab=1.5,cex.axis=1.25) 
    abline(h = 0,lty=2,lwd=1.25); axis(side = 1,at = xsup,labels = xsup,cex.axis=1.25)
    points(xsup,ym[design$rho[iid]==rhox[j],h],col=cols[3],pch=20)
    
    # drawing dML-max
    ym=t(mapply(function(k)Tau_maxdef[k,c("mean"),1:((design$M[k]-1))][c(1,(design$M[k]-1))],iid))-matrix(c(-2,2),12,2,byrow = TRUE)
    points(xsup,ym[design$rho[iid]==rhox[j],h],type="b",col=cols[1],lwd=1.5);
    points(xsup,ym[design$rho[iid]==rhox[j],h],col=cols[1],pch=20)
    
    # drawing dML-mean
    ym=t(mapply(function(k)Tau_meandef[k,c("mean"),1:((design$M[k]-1))][c(1,(design$M[k]-1))],iid))-matrix(c(-2,2),12,2,byrow = TRUE)
    points(xsup,ym[design$rho[iid]==rhox[j],h],type="b",col=cols[2],lwd=1.5);
    points(xsup,ym[design$rho[iid]==rhox[j],h],col=cols[2],pch=20)
  }
}
add_legend("bottom",fill = cols,legend = c("dML-max","dML-mean","fEM"),border = FALSE,bty = "n",ncol = 3,cex=1.5)
#dev.off()


## Figure S2 
iid=which(design$M==6)
cols=c("chocolate4","darkorange3","darkolivegreen")
xsup=c(150,250,500,1000); rhox=c(0.15,0.50,0.85); pm=0.35
#tikzDevice::tikz(file='draft_paper/figS2.tex',width=5,height=7.5)
x11();
par(mfcol=c(3,2),mai=c(0.7, 0.8, 0.3, 0.15))
for(h in 1:4){ #over tau
  for(j in 1:3){ #over rho
    xmain=""
    if(h==1&j==1){xmain=latex2exp::TeX("$\\tau_1")};if(h==2&j==1){xmain=latex2exp::TeX("$\\tau_2")}
    if(h==3&j==1){xmain=latex2exp::TeX("$\\tau_4")};if(h==4&j==1){xmain=latex2exp::TeX("$\\tau_5")}
    if(h==1){ymain=paste0("p = ",rhox[j])}else{ymain=""}
    
    # drawing fEM
    ym=t(mapply(function(k)Tau_fEM[k,c("mean"),1:((design$M[k]-1))][c(1:2,4:5)],iid))-matrix(c(-2,-1,1,2),12,4,byrow = TRUE)
    plot(xsup,ym[design$rho[iid]==rhox[j],h],bty="n",ylab=ymain,xlab="",ylim=c(-pm,+pm),xlim=c(50,1050),type="b",xaxt="none",col=cols[3],lwd=1.5,main=xmain,cex.main=2,cex.lab=1.5,cex.axis=1.25) 
    abline(h = 0,lty=2,lwd=1.25); axis(side = 1,at = xsup,labels = xsup,cex.axis=1.25)
    points(xsup,ym[design$rho[iid]==rhox[j],h],col=cols[3],pch=20)
    
    # drawing dML-max
    ym=t(mapply(function(k)Tau_maxdef[k,c("mean"),1:((design$M[k]-1))][c(1:2,4:5)],iid))-matrix(c(-2,-1,1,2),12,4,byrow = TRUE)
    points(xsup,ym[design$rho[iid]==rhox[j],h],type="b",col=cols[1],lwd=1.5);
    points(xsup,ym[design$rho[iid]==rhox[j],h],col=cols[1],pch=20)
    
    # drawing dML-mean
    ym=t(mapply(function(k)Tau_meandef[k,c("mean"),1:((design$M[k]-1))][c(1:2,4:5)],iid))-matrix(c(-2,-1,1,2),12,4,byrow = TRUE)
    points(xsup,ym[design$rho[iid]==rhox[j],h],type="b",col=cols[2],lwd=1.5);
    points(xsup,ym[design$rho[iid]==rhox[j],h],col=cols[2],pch=20)
  }
}
add_legend("bottom",fill = cols,legend = c("dML-max","dML-mean","fEM"),border = FALSE,bty = "n",ncol = 3,cex=1.5)
#dev.off()

##### Application 1 #####
load("casestudy1.rda")

# Table 3a
X = cbind(1:5,Cats$V1,Cats$V2,Cats$V3)
Xtab = X
Xtab_tex = xtable::xtable(Xtab,digits = c(0,rep(2,NCOL(X))))
attributes(Xtab_tex)$caption = "Application 1: "
attributes(Xtab_tex)$label = "tab3a"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
Xtab_tex$`1` = paste0("$r=",Xtab_tex$`1`,"$")
Xtab_tex
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x},include.rownames = FALSE)

# Figure 2
cols=c("deepskyblue4","darkolivegreen4","darkorange2","darkorange4","deepskyblue3")
ax=c(0.05,0.25,0.50,0.75,0.99); pchs=c(0,1,5,2,0)
xsup=seq(from=0.5,to=10,length.out=1001)
glbls=c("G0","G1","G2","G3","G4")

#tikzDevice::tikz(file='draft_paper/fig2.tex',width=7.5,height=2.5)
x11()
par(mfrow=c(1,3))
for(j in 1:length(Ydata)){
  C=Cats[[j]];r=1;mux=trapezoidal_fn(x=xsup,lb = C[r,1],m1 = C[r,2],m2 = C[r,3],ub = C[r,4]); plot(xsup,mux,lwd=1.5,bty="n",ylab="",xlab="",type="l",main=paste0("X",j),cex.main=1.25,xlim=c(0,10),ylim=c(0,1.15),col=cols[r],lty=1); text(mean(C[r,2:3])-0.5,1.1,glbls[r],col = cols[r],cex=1.15,cex.axis=1.25)
  #for(k in 1:length(ax)){points(c(min(xsup[mux>ax[k]]),max(xsup[mux>ax[k]])),c(min(mux[mux>ax[k]]),min(mux[mux>ax[k]])),col=cols[r],pch=pchs[r],lty=2,cex=1.5)}
  for(r in 2:NROW(C)){
    mux=trapezoidal_fn(x=xsup,lb = C[r,1],m1 = C[r,2],m2 = C[r,3],ub = C[r,4]); points(xsup,mux,lwd=1.5,type="l",col=cols[r],lty=1); text(mean(C[r,2:3]),1.1,glbls[r],col = cols[r],cex=1.15)
    #for(k in 1:length(ax)){points(c(min(xsup[mux>ax[k]]),max(xsup[mux>ax[k]])),c(min(mux[mux>ax[k]]),min(mux[mux>ax[k]])),col=cols[r],pch=pchs[r],cex=1.5)}
  }
  D=Ydata[[j]]; for(i in 1:NROW(D)){segments(x0=D[i,1],y0=0,x1=D[i,1],y1=1,lwd=1.25,lty=2,col="azure4")}
}
#dev.off()

# Computing fuzzy frequency tables for each pair of (X1,X2,X3)
FFtab_obj = FuzzyFrequency_list(Ydata,Cats) #it returns a list

# Figure 3 (only pair X2,X3)
#tikzDevice::tikz(file='draft_paper/fig3.tex',width=8.5,height=7)
plot_FuzzyFrequency_table(FFtab_obj$FFtables$`2,3`,NROW(Cats[[1]]),new.dev = TRUE) #new.dev = FALSE) deactivate x11() inside the function
#dev.off()

# Estimating R via EM and defuzz-methods
J=length(Ydata)
M=NROW(Cats[[1]])
R_EM=R_maxd=R_meand=matrix(0,J,J);
Tau_EM=list() 
for(j in 1:NROW(FFtab_obj$Pairs)){
  jjd = as.numeric(unlist(strsplit(x = (names(FFtab_obj$FFtables)[j]),split = ",")))
  print(jjd)
  out = estimPars_EM(R=M,C=M,FFtab=FFtab_obj$FFtables[[j]],use.optim=TRUE,K=1000,eps=1e-7,starting.method = "max",verbose=FALSE)
  R_EM[jjd[1],jjd[2]] = out$rho_EM$rho
  R_maxd[jjd[1],jjd[2]] = estimPars_maxdefuzz(FFtab =FFtab_obj$FFtables[[j]],M,M)$rho
  R_meand[jjd[1],jjd[2]] = estimPars_meandefuzz(FFtab =FFtab_obj$FFtables[[j]],M,M)$rho
  Tau_EM[[j]] = out$Tau_EM
}
R_EM = t(R_EM)+R_EM
R_maxd = t(R_maxd)+R_maxd
R_meand = t(R_meand)+R_meand
diag(R_EM)=diag(R_maxd)=diag(R_meand)=1
colnames(R_EM)=rownames(R_EM)=paste0("X",1:3)

# Table 3b
Xtab = R_EM
Xtab_tex = xtable::xtable(Xtab,digits = rep(5,NCOL(Xtab)+1))
attributes(Xtab_tex)$caption = "Application 1: "
attributes(Xtab_tex)$label = "tab3b"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
#Xtab_tex$`1` = paste0("$r=",Xtab_tex$`1`,"$")
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x},include.rownames = FALSE)


##### Application 2 #####
# Note: Here we start from degrees of inclusion instead of crisp/fuzzy data and fuzzy categories
load("casestudy2.rda") 
nms = c("Sun","Hum","Pre","Alt","Max")
M = 3; I = NROW(datax); J = 5

# All pairs of variables (indices)
Jjd = expand.grid(1:J,1:J); Jjd = Jjd[Jjd[,1]<Jjd[,2],]
Iid = expand.grid(1:M,1:M) #overall indices for the categories (pairs)

# Computing (pairwise) fuzzy frequencies matrix for all variables
FFlist = list()
for(u in 1:NROW(Jjd)){
  v=as.numeric(Jjd[u,])
  E1=datax[,grep(x=colnames(datax),pattern = nms[v[1]])] #Sun
  E2=datax[,grep(x=colnames(datax),pattern = nms[v[2]])] #Hum
  
  H = NROW(Iid); FF = matrix(NA,I+1,H)
  for(h in 1:H){
    j = Iid[h,1]; k = Iid[h,2]
    e = apply(cbind(E1[,j],E2[,k]),1,min)
    FF[,h] = FCounts(e,TRUE)$FECount #matrix of fuzzy frequencies
  }
  rownames(FF)=0:I; colnames(FF)=apply(Iid,1,function(x)paste(x,collapse=""))
  FFlist[[u]]=FF
}
names(FFlist) = as.character(apply(Jjd,1,function(x)paste(x,collapse=",")))

# Figure 4 (only pair X2,X3: Humidity vs. Rainfall)
#tikzDevice::tikz(file='draft_paper/fig4.tex',width=8.5,height=7)
plot_FuzzyFrequency_table(FFlist$`2,3`,M,new.dev = TRUE)#new.dev = FALSE)
#dev.off()

# Estimating R via EM and defuzz-methods
R_EM=R_maxd=R_meand=matrix(0,J,J);
Tau_EM=out=list() 
for(j in 1:NROW(FFlist)){
  jjd = as.numeric(unlist(strsplit(x = (names(FFlist)[j]),split = ",")))
  print(jjd)
  out[[j]] = estimPars_EM(R=M,C=M,FFtab=FFlist[[j]],use.optim=TRUE,K=1000,eps=1e-7,starting.method = "max",verbose=FALSE)
  R_EM[jjd[1],jjd[2]] = out[[j]]$rho_EM$rho
  R_maxd[jjd[1],jjd[2]] = estimPars_maxdefuzz(FFtab =FFlist[[j]],M,M)$rho
  R_meand[jjd[1],jjd[2]] = estimPars_meandefuzz(FFtab =FFlist[[j]],M,M)$rho
  Tau_EM[[j]] = out$Tau_EM
}
R_EM = t(R_EM)+R_EM
R_maxd = t(R_maxd)+R_maxd
R_meand = t(R_meand)+R_meand
diag(R_EM)=diag(R_maxd)=diag(R_meand)=1
colnames(R_EM)=rownames(R_EM)=nms

# Table 4a
Xtab = R_EM
Xtab_tex = xtable::xtable(Xtab,digits = rep(5,NCOL(Xtab)+1))
attributes(Xtab_tex)$caption = "Application 2: "
attributes(Xtab_tex)$label = "tab4a"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
#Xtab_tex$`1` = paste0("$r=",Xtab_tex$`1`,"$")
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x},include.rownames = FALSE)

# Path analysis on R_EM
isSymmetric(R_EM)
pracma::isposdef(R_EM)
R_EM_pd = as.matrix(Matrix::nearPD(R_EM,corr = TRUE)$mat)

library(lavaan)
m1 = "Pre~Hum+Sun \n Hum~Max; Max~Alt \n Sun~~0*Alt"
sem.out = lavaan::sem(m1,sample.cov=R_EM_pd,sample.nobs=NROW(datax))
lavaan::summary(sem.out,standardized=T,rsquare=T)

#norm(lavInspect(sem.out, "cov.ov"))^2/norm(R_EM_pd)^2

# Figure 5
#tikzDevice::tikz(file='draft_paper/fig5.tex',width=8.5,height=4)
semPlot::semPaths(object = sem.out, what = "eq", nCharNodes=6, sizeMan=10,edge.label.cex=1,style = "lisrel",
                  layout = "spring",residuals = FALSE,color = c("darkorange2",rep("darkolivegreen3",4)),
                  edge.color = "gray30",edge.width=1,esize=4.1,asize=4.1,minimum=0.1)
#dev.off()

# Table 4b
A=summary(sem.out)
Xtab = cbind(c(A$PE$rhs[1:4],NA),c(A$PE$lhs[1:4],NA),c(round(A$PE$est[1:4],4),NA),c(round(A$PE$se[1:4],4),NA),round(A$PE$est[6:10],4),round(A$PE$se[6:10],4))
#Xtab = cbind(A$PE$lhs[-1][1:4],round(cbind(A$PE$est[1:4],A$PE$se[1:4])[-5],4))
Xtab_tex = xtable::xtable(Xtab,digits = c(0,rep(3,NCOL(Xtab))))
attributes(Xtab_tex)$caption = "Application 2: "
attributes(Xtab_tex)$label = "tab4b"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x},include.rownames = FALSE)




