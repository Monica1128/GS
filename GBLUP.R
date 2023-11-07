setwd("  ") #输入路径
gen<-as.matrix(read.csv(file="gen.csv",header=F)) #导入基因型
kin<-as.matrix(read.csv(file ="K.csv", header=F)) #导入kinship关系矩阵
phe<-as.matrix(read.csv(file="phe.csv",header=T)) #导入表型
foldidID<-as.matrix(read.csv(file="foldidID.csv",header=F)) #用于10-fold交叉验证
phe<-apply(phe[,-1],2,as.numeric)
x<-apply(gen,2,as.numeric)
nfold<-10

y<-as.matrix(phe[,1]) #第一个性状

kk<-kin
id<-NULL
fold<-NULL
GBLUP_r2<-NULL
GBLUP_yp<-NULL
GBLUP_yt<-NULL
for (i in 1:5){
  cat(i)
  foldid<-foldidID[,i]
  load(file="mixed.RData")
  n<-length(y)
  x.<-matrix(1,n,1)
  par<-mixed(x=x.,y=y,kk=kk,method="REML",eigen=T)
  cv.mixed<-function(x.,y,kk,nfold=nfold,foldid=NULL){
    n<-length(y)
    yp_Gblup<-NULL
    yt_Gblup<-NULL
    for(k in 1:nfold){
      i1<-which(foldid!=k)
      i2<-which(foldid==k)
      x1<-x.[i1,,drop=F]
      y1<-y[i1,,drop=F]
      k11<-kk[i1,i1]
      parm<-mixed(x=x1,y=y1,kk=k11,method="REML",eigen=TRUE)
      qq<-parm[[2]]
      delta<-qq[[1]]
      u<-qq[[2]]
      beta<-as.matrix(parm[[1]]$beta,ncol(x),1)
      va<-parm[[1]]$va[1]
      ve<-parm[[1]]$ve[1]
      x2<-x.[i2,,drop=F]
      y2<-y[i2,,drop=F]
      k41<-kk[i2,i1]
      h<-1/(delta*va+ve)
      y3<-x2%*%beta+va*k41%*%u%*%diag(h)%*%t(u)%*%(y1-x1%*%beta)
      #y3<-x2%*%beta+va*k41%*%solve(k11*va+diag(length(y1))*ve)%*%(y1-x1%*%beta)
      fold<-c(fold,rep(k,length(y2)))
      yt_Gblup<-c(yt_Gblup,y2)
      yp_Gblup<-c(yp_Gblup,y3)
      id<-c(id,i2)
    }
    pred<-data.frame(fold,id,yt_Gblup,yp_Gblup)
    R2<-cor(pred$yt_Gblup,pred$yp_Gblup)^2
    return(list(data.frame(r2=R2),pred))
  }
  Gblup<-cv.mixed(x=x.,y=y,kk=kk,nfold=nfold,foldid=foldid)
  Gblup_r2<-as.numeric(Gblup[[1]])
  GBLUP_r2<-c(GBLUP_r2,Gblup_r2)
  gblup<-Gblup[[2]]
  yp_Gblup<-as.numeric(gblup[,4])
  yt_Gblup<-as.numeric(gblup[,3])
  GBLUP_yp<-cbind(GBLUP_yp,yp_Gblup)
  GBLUP_yt<-cbind(GBLUP_yt,yt_Gblup)
}
R2_Gblup<-as.matrix(GBLUP_r2)
GBLUP_R2<-colMeans(R2_Gblup)
write.csv(x=R2_Gblup,file="./SCI/SCI_GBLUP_R2.csv",row.names=T)
write.csv(x=GBLUP_R2,file="./SCI/SCI_GBLUP_mean_R2.csv",row.names=T)
write.csv(x=GBLUP_yt,file="./SCI/SCI_GBLUP_yt.csv",row.names=T)
write.csv(x=GBLUP_yp,file="./SCI/SCI_GBLUP_yp.csv",row.names=T)


y<-as.matrix(phe[,2]) #第二个性状

kk<-kin
id<-NULL
fold<-NULL
GBLUP_r2<-NULL
GBLUP_yp<-NULL
GBLUP_yt<-NULL
for (i in 1:5){
  cat(i)
  foldid<-foldidID[,i]
  load(file="mixed.RData")
  n<-length(y)
  x.<-matrix(1,n,1)
  par<-mixed(x=x.,y=y,kk=kk,method="REML",eigen=T)
  cv.mixed<-function(x.,y,kk,nfold=nfold,foldid=NULL){
    n<-length(y)
    yp_Gblup<-NULL
    yt_Gblup<-NULL
    for(k in 1:nfold){
      i1<-which(foldid!=k)
      i2<-which(foldid==k)
      x1<-x.[i1,,drop=F]
      y1<-y[i1,,drop=F]
      k11<-kk[i1,i1]
      parm<-mixed(x=x1,y=y1,kk=k11,method="REML",eigen=TRUE)
      qq<-parm[[2]]
      delta<-qq[[1]]
      u<-qq[[2]]
      beta<-as.matrix(parm[[1]]$beta,ncol(x),1)
      va<-parm[[1]]$va[1]
      ve<-parm[[1]]$ve[1]
      x2<-x.[i2,,drop=F]
      y2<-y[i2,,drop=F]
      k41<-kk[i2,i1]
      h<-1/(delta*va+ve)
      y3<-x2%*%beta+va*k41%*%u%*%diag(h)%*%t(u)%*%(y1-x1%*%beta)
      #y3<-x2%*%beta+va*k41%*%solve(k11*va+diag(length(y1))*ve)%*%(y1-x1%*%beta)
      fold<-c(fold,rep(k,length(y2)))
      yt_Gblup<-c(yt_Gblup,y2)
      yp_Gblup<-c(yp_Gblup,y3)
      id<-c(id,i2)
    }
    pred<-data.frame(fold,id,yt_Gblup,yp_Gblup)
    R2<-cor(pred$yt_Gblup,pred$yp_Gblup)^2
    return(list(data.frame(r2=R2),pred))
  }
  Gblup<-cv.mixed(x=x.,y=y,kk=kk,nfold=nfold,foldid=foldid)
  Gblup_r2<-as.numeric(Gblup[[1]])
  GBLUP_r2<-c(GBLUP_r2,Gblup_r2)
  gblup<-Gblup[[2]]
  yp_Gblup<-as.numeric(gblup[,4])
  yt_Gblup<-as.numeric(gblup[,3])
  GBLUP_yp<-cbind(GBLUP_yp,yp_Gblup)
  GBLUP_yt<-cbind(GBLUP_yt,yt_Gblup)
}
R2_Gblup<-as.matrix(GBLUP_r2)
GBLUP_R2<-colMeans(R2_Gblup)
write.csv(x=R2_Gblup,file="./MAT/MAT_GBLUP_R2.csv",row.names=T)
write.csv(x=GBLUP_R2,file="./MAT/MAT_GBLUP_mean_R2.csv",row.names=T)
write.csv(x=GBLUP_yt,file="./MAT/MAT_GBLUP_yt.csv",row.names=T)
write.csv(x=GBLUP_yp,file="./MAT/MAT_GBLUP_yp.csv",row.names=T)
