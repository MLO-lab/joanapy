setwd("pathto/rawData/single_cell/SMM_MM/")

BMAsym_norm=read.table("./SMM_MMDegBM.txt")



indZero=which(BMAsym_norm[,"p_val_adj"]==0)
tiny_value <- 1e-10
near_one <- 1-1e-10
BMAsym_norm[BMAsym_norm==0]=tiny_value

BMAsym_norm[BMAsym_norm==1]=near_one



BMpval=data.frame(rownames(BMAsym_norm),BMAsym_norm[,"p_val"])

BMqval=data.frame(rownames(BMAsym_norm),BMAsym_norm[,"p_val_adj"])

write.table(BMqval,"./BMqval.txt",col.names = F,row.names = F)

CircAsym_norm=read.table("./SMM_MMDegCirc.txt")

indZero=which(CircAsym_norm[,"p_val_adj"]==0)


CircAsym_norm[CircAsym_norm==0]=tiny_value

CircAsym_norm[CircAsym_norm==1]=near_one



Circpval=data.frame(rownames(CircAsym_norm),CircAsym_norm[,"p_val"])

Circqval=data.frame(rownames(CircAsym_norm),CircAsym_norm[,"p_val_adj"])

write.table(Circqval,"./CIRCqval.txt",col.names = F,row.names = F)


