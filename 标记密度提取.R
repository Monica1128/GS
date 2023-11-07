setwd("   ") #输入文件的路径
gen<-read.csv("gen.csv",header=F) #基因型数据  #怎么重复！！！！！！
number<-seq(1,157989,by=200) #790
gen1<-gen[number,]
write.csv(gen1,"790.csv",row.names=F)
number<-seq(1,157989,by=100) #1580
gen2<-gen[number,]
write.csv(gen2,"1580.csv",row.names=F)
number<-seq(1,157989,by=50) #3160
gen3<-gen[number,]
write.csv(gen3,"3160.csv",row.names=F)
number<-seq(1,157989,by=25) #6320
gen4<-gen[number,]
write.csv(gen4,"6320.csv",row.names=F)
number<-seq(1,157989,by=15) #10533
gen5<-gen[number,]
write.csv(gen5,"10533.csv",row.names=F)
number<-seq(1,157989,by=8) #19749
gen6<-gen[number,]
write.csv(gen6,"19749.csv",row.names=F)
number<-seq(1,157989,by=5) #31598
gen7<-gen[number,]
write.csv(gen7,"31598.csv",row.names=F)
number<-seq(1,157989,by=3) #52663
gen8<-gen[number,]
write.csv(gen8,"52663.csv",row.names=F)
number<-seq(1,157989,by=2) #78995
gen8<-gen[number,]
write.csv(gen8,"78995.csv",row.names=F)


#生成kinship用于GBLUP
G<-read.csv("78995.csv",header=F) 
G<-apply(G,2,as.numeric)
Gt<-t(G)
Gt<-apply(Gt,2,as.numeric)
k<-Gt%*%G
k<-k/mean(diag(k))
write.csv(k,"K78995.csv",row.names=F)

G<-read.csv("52663.csv",header=F) 
G<-apply(G,2,as.numeric)
Gt<-t(G)
Gt<-apply(Gt,2,as.numeric)
k<-Gt%*%G
k<-k/mean(diag(k))
write.csv(k,"K52663.csv",row.names=F)

G<-read.csv("31598.csv",header=F) 
G<-apply(G,2,as.numeric)
Gt<-t(G)
Gt<-apply(Gt,2,as.numeric)
k<-Gt%*%G
k<-k/mean(diag(k))
write.csv(k,"K31598.csv",row.names=F)

G<-read.csv("19749.csv",header=F) 
G<-apply(G,2,as.numeric)
Gt<-t(G)
Gt<-apply(Gt,2,as.numeric)
k<-Gt%*%G
k<-k/mean(diag(k))
write.csv(k,"K19749.csv",row.names=F)
  
G<-read.csv("10533.csv",header=F) 
G<-apply(G,2,as.numeric)
Gt<-t(G)
Gt<-apply(Gt,2,as.numeric)
k<-Gt%*%G
k<-k/mean(diag(k))
write.csv(k,"K10533.csv",row.names=F)

G<-read.csv("6320.csv",header=F) 
G<-apply(G,2,as.numeric)
Gt<-t(G)
Gt<-apply(Gt,2,as.numeric)
k<-Gt%*%G
k<-k/mean(diag(k))
write.csv(k,"K6320.csv",row.names=F)

G<-read.csv("3160.csv",header=F) 
G<-apply(G,2,as.numeric)
Gt<-t(G)
Gt<-apply(Gt,2,as.numeric)
k<-Gt%*%G
k<-k/mean(diag(k))
write.csv(k,"K3160.csv",row.names=F)

G<-read.csv("1580.csv",header=F) 
G<-apply(G,2,as.numeric)
Gt<-t(G)
Gt<-apply(Gt,2,as.numeric)
k<-Gt%*%G
k<-k/mean(diag(k))
write.csv(k,"K1580.csv",row.names=F)

G<-read.csv("790.csv",header=F) 
G<-apply(G,2,as.numeric)
Gt<-t(G)
Gt<-apply(Gt,2,as.numeric)
k<-Gt%*%G
k<-k/mean(diag(k))
write.csv(k,"K790.csv",row.names=F)

