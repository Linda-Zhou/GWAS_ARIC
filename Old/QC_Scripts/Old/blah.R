
hits<-read.csv("/users/asurapan/04_ARIC.GWAS.with.FAST/invZHX3_Visit3_W/invZHX3_Visit3_hits.txt")
head(hits)

SNPlookup=read.csv("/dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/rsid/GRCh38_dbSNP151_rsid_final_USE.txt")
head(SNPlookup)
names(SNPlookup)<-c("ChrId", "rsID_dbSNP151")
hits=merge(hits, SNPlookup)
names(hits)
hits=hits[ChrId,SNPId, other stuff]



