setwd("pathto/rawdata/LDCmouse/")

rna <- read.table("./RNA_seq_mouse_with_entrez.tsv",
                  header = T,
                  sep = "\t",
                  stringsAsFactors = FALSE)
chip <- read.table("./chip_seq_mouse_with_entrez.tsv",
                  header = T,
                  sep = "\t",
                  stringsAsFactors = FALSE)




rna=na.omit(rna)
chip=na.omit(chip)
dim(rna)
dim(chip)

colnames(rna)
colnames(chip)

rna_joana=rna[,c("entrez","qval")]
chip_joana=chip[,c("entrez","qval")]
rna_joana=rna_joana[-which(duplicated(rna_joana$entrez)),]
Deg_Dep=merge(rna_joana,chip_joana,all = T,by="entrez")
colnames(Deg_Dep)=c("entrez","rna_qval","RNA_logFC","chip_qval")
Deg_Dep=Deg_Dep[-which(is.na(Deg_Dep$rna_qval)),]
dim(Deg_Dep)
dim(rna_joana)

#write.table(Deg_Dep,"./Deg_Dep.txt",row.names = F,col.names = T)
write.table(rna_joana,"./RNA-1.txt",row.names = F,col.names = F)
write.table(chip_joana,"./chip-2.txt",row.names = F,col.names = F)
