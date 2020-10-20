###############################
#Generating indiv file for FAST#
################################

#Key is located here: dcs01/arking/ARIC_static/ARIC_Data/gwa_official_idlist_rev121023.csv

#In bash: Get a list of the samples from the TOPMed Imputation files
#cd /dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/matchingids
#bcftools query -l /dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/chr10.dose.vcf.gz > W_ARIC_sampIDs.from.topmed

##Run in R
# Adi setwd("/dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/matchingids")

original_phenofile <- "/dcs01/arking/arkinglab/static/aric/phenotypes/f3v2/aric.allphenotypes.trim.final.v4.f3v2.txt"
ids_order_file <- "/dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/matchingids/W_ARIC_sampIDs.from.topmed"
key <- "/dcs01/arking/ARIC_static/ARIC_Data/gwa_official_idlist_rev121023.csv"

ophe <- read.table(original_phenofile,head=T,strings=F)
help.key <- read.csv(key, head = T)
rownames(ophe) <- ophe$IID

ids.topmed <- scan(ids_order_file,what=character())
ids.tomatch = gsub(".*_","",ids.topmed)

#Make key with all samples that matched
key.matched = help.key[help.key$gwasid %in% ids.tomatch,]

#Make one table with all IDs. Need to preseve row order
ids.topmed.table = as.data.frame(matrix(ncol = 2, nrow = 9345))
ids.topmed.table[,1] = ids.topmed
ids.topmed.table[,2] = gsub(".*_","",ids.topmed.table[,1])
names(ids.topmed.table) = c("full", "gwasid")

library(plyr) #join command allows for retention of row order
key.mod = join(ids.topmed.table, help.key, type = "inner")
#Order of listing the data frame in the above command matters. The order retained is that of the FIRST data frame.

#Make individual ID file, both with full id and IID
# Adi setwd("/dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation")
write.table(key.mod$pid, "../datasets/W_ARIC_topmed.imputation.iid", quote=F,row.names=F,col.names=F)
#This is the file that will be used in future analyses
