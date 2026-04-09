setwd("pathto/rawData/LungAdeno/")
library(limma)
library(readxl)
############ Differentially Expression Analysis RNA data
tableS= read_excel("./1-s2.0-S0092867420307443-mmc2.xlsx",sheet = "Table S2D")
subtype=tableS[70,]
matTableS=as.matrix(tableS)
ids=matTableS[71:nrow(matTableS),1]
subSolid=which(subtype=="solid")
subPapillary=which(subtype=="papillary")

tableS2solid=tableS[,subSolid]
tableS2papillary=tableS[,subPapillary]
tableS2=cbind(tableS2solid,tableS2papillary)

expressionData=tableS2[71:nrow(tableS2),]


covar_df=tableS2[c(70,10,12,43,50,51,49,65),]
covar_df=t(covar_df)
naInEGFR=which(covar_df[,4]=="NA")
covar_df[naInEGFR,4]="0"
covar_df[-naInEGFR,4]="1"


data_f=data.frame(rep(0,18))
for(i in 1:8){
  data_f=cbind(data_f,as.factor(covar_df[,i]))
}

data_f=data_f[,-1]


data_f=data_f[,-8]

expressionData=as.matrix(expressionData)
exprmat=matrix(0, nrow = nrow(expressionData), ncol = ncol(expressionData))
for(i in 1:ncol(expressionData)){
  exprmat[,i]=as.numeric(expressionData[,i])
}



row.names(exprmat)=ids
dim(exprmat)


colnames(data_f)=c("histologic","smoking","R_O","EGFR","KRAS","STK11","TP53")

design <- model.matrix(~ 0+histologic+smoking+R_O+EGFR+KRAS+STK11+TP53, data=data_f) 





fit <- lmFit(exprmat, design)
dim(fit$coefficients)

cont.matrix=makeContrasts(histologicpapillary-histologicsolid,levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2,0.05)


topTrans=topTable(fit2,adjust.method = "fdr",sort.by = "p",number = Inf)

idduplicated1=which(topTrans[,1]=="43892")
idduplicated2=which(topTrans[,1]=="43891")
topTrans=topTrans[-c(idduplicated1,idduplicated2),]
write.table(topTrans,"./DegTrans.txt",col.names = T,row.names = T)

########### Differentially Expression Analysis protein data
tableS= read_excel("./1-s2.0-S0092867420307443-mmc3.xlsx",sheet = "Table S3A")
subtype=tableS[70,]
matTableS=as.matrix(tableS)
ids=matTableS[71:nrow(matTableS),3]
subSolid=which(subtype=="solid")
subPapillary=which(subtype=="papillary")
tableS2solid=tableS[,subSolid]
tableS2papillary=tableS[,subPapillary]

tableS2=cbind(tableS2solid,tableS2papillary)

expressionData=tableS2[71:nrow(tableS2),]




covar_df=tableS2[c(70,10,12,43,50,51,49,65),]

covar_df=t(covar_df)

naInEGFR=which(covar_df[,4]=="NA")
covar_df[naInEGFR,4]="0"
covar_df[-naInEGFR,4]="1"

#covar_df=covar_df[,-4]

data_f=data.frame(rep(0,nrow(covar_df)))
for(i in 1:8){
  data_f=cbind(data_f,as.factor(covar_df[,i]))
}

data_f=data_f[,-1]

data_f=data_f[,-8]
#data_f=data_f[,-7]

expressionData=as.matrix(expressionData)
exprmat=matrix(0, nrow = nrow(expressionData), ncol = ncol(expressionData))
for(i in 1:ncol(expressionData)){
  exprmat[,i]=as.numeric(expressionData[,i])
}

row.names(exprmat)=ids
dim(exprmat)

colnames(data_f)=c("histologic","smoking","R_O","EGFR","KRAS","STK11","TP53")
design <- model.matrix(~ 0+histologic+smoking+R_O+EGFR+KRAS+STK11+TP53, data=data_f) 
fit <- lmFit(exprmat, design)
cont.matrix=makeContrasts(histologicpapillary-histologicsolid,levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2,0.01)
topProt=topTable(fit2,adjust.method = "fdr",sort.by = "p",number = Inf)
which(topProt[,1]=="PTK7")
which(topProt[,1]=="IDH3A")
topProt[which(topProt[,1]=="PTK7"),]
topProt[which(topProt[,1]=="IDH3A"),]


simpleRemoveDups<-function(sortedTable) {
  indup=which(duplicated(sortedTable[,1]))
  sortedTable=sortedTable[-indup,]
  return(sortedTable)
}
topProtFiltered=simpleRemoveDups(topProt)
write.table(topProtFiltered,file = "./DepProt.txt",row.names = T,col.names = T)

############## Merging Deg and Dep

DegTrans=read.table("./DegTrans.txt")
rownames(DegTrans)=DegTrans$ID
DegProt=read.table("./DepProt.txt")
rownames(DegProt)=DegProt$ID


DegDep=merge(x=DegTrans[,c("ID","P.Value","adj.P.Val","logFC")],y=DegProt[,c("ID","P.Value","adj.P.Val","logFC")],by="ID",all.x=TRUE)
indna=which(is.na(DegDep$P.Value.x))

DegDep=DegDep[-indna,]

which(is.na(DegDep$adj.P.Val.x))

colnames(DegDep)=c("ID","pvalRNA","qvalRNA","foldRNA","pvalProt","qvalProt","foldProt")

rna=DegDep[,c("ID","pvalRNA")]
prot=DegDep[,c("ID","pvalProt")]
write.table(rna,"./omics1.txt",col.names = F,row.names = F)
write.table(prot,"./omics2.txt",col.names = F,row.names = F)
