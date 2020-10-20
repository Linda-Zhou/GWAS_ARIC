
pheno1=read.table("invreslogNBL1_SeqId_2944_66.visit3.W.FAST.runfile.iid.txt")

pheno2=read.table("SeqId_2944_66.pheno", header=T)

merged=merge(pheno1, pheno2, by.x="V2", by.y="IID", all.x=TRUE, sort=FALSE)

merged$V6=merged$SeqId_2944_66

merged=merged[,!names(merged) %in% c("FID","SeqId_2944_66")]
merged=merged[,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20")]
write.table(merged,"invreslogNBL1_SeqId_2944_66.visit3.W.FAST.runfile.iid.txt.jn",sep="\t", col.names=F, quote=F,row.names=F)

system("cat invreslogNBL1_SeqId_2944_66.visit3.W.FAST.runfile.iid.txt.head invreslogNBL1_SeqId_2944_66.visit3.W.FAST.runfile.iid.txt.jn >invreslogNBL1_SeqId_2944_66.visit3.W.FAST.runfile.iid.txt.jn.wheader")

