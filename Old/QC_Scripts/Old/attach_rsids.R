
hits<-read.csv("/users/asurapan/04_ARIC.GWAS.with.FAST/invZHX3_Visit3_W/invZHX3_Visit3_hits.txt", sep="\t")
hits<-head(hits,50)
#summary(hits)

#HitSNPs=as.vector(hits$SNPID)
#HitSNPs=HitSNPs[1:2]

rsid_from_snp <- function(SNP) {
	snp_rsid=system(paste0("grep ",SNP," /dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/rsid/GRCh38_dbSNP151_rsid_final_USE.txt "), intern=TRUE)
        #cat(snp_rsid[1],"\n")
        snp_rsid=strsplit(snp_rsid, "\t")[1]
        rsid<-as.vector(snp_rsid)[[1]][2]
	return(rsid)
}
   
#for (SNIP in HitSNPs){
#	cat(rsid_from_snp(SNIP))
#}

hits$rsid <- mapply(rsid_from_snp,hits$SNP)


write.table(hits, "/users/asurapan/04_ARIC.GWAS.with.FAST/invZHX3_Visit3_W/invZHX3_Visit3_hits.txt",sep="\t",quote=FALSE,row.names=FALSE)

