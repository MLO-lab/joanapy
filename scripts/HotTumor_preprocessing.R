setwd("pathto/rawData/HotTumor/")

library(readxl)

RNA=read_excel("./NIHMS1694599-supplement-Supplementary_table_2.xlsx",sheet = "DEG (RNA)")
colnames(RNA)=RNA[3,]
RNA=RNA[-c(1,2,3),]
RNA=data.frame(RNA)
rownames(RNA)=RNA[,1]
RNA=RNA[,-1]
RNA2=apply(RNA, 2, as.numeric)
rownames(RNA2)=rownames(RNA)



prot=read_excel("./NIHMS1694599-supplement-Supplementary_table_2.xlsx",sheet = "DEG (global)")
colnames(prot)=prot[3,]
prot=prot[-c(1,2,3),]
prot=data.frame(prot)
rownames(prot)=prot[,1]
prot=prot[,-c(1,2)]
prot2=apply(prot, 2, as.numeric)
rownames(prot2)=rownames(prot)


x="Hot.Tumor"    
    
logFCRNA=RNA2[,x]
logFCRNA=logFCRNA*abs(1/logFCRNA)

logFCprot=prot2[,x]
logFCprot=logFCprot*abs(1/logFCprot)

pvalsRNA=abs(RNA2[,x])
pvalsRNA=10^(-pvalsRNA)
pvalsRNA=data.frame("ID"=names(pvalsRNA),pvalsRNA)

pvalsProt=abs(prot2[,x])
pvalsProt=10^(-pvalsProt)
pvalsProt=data.frame("ID"=names(pvalsProt),pvalsProt)

tr_pr=merge(pvalsRNA,pvalsProt,by="ID",all=TRUE)

rna=tr_pr[,c("ID","pvalsRNA")]
prot=tr_pr[,c("ID","pvalsProt")]
write.table(rna,"./omics1.txt",col.names = F,row.names = F)
write.table(prot,"./omics2.txt",col.names = F,row.names = F)








