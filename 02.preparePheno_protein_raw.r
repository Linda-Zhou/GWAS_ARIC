################################################
#Generating phenotype file: Spatial QRS-T Angle#
################################################

#
# File containing all phenotypes
# Variables should be separated by space or tab, missing values coded as NA
# Header line should contain variable names
# Compulsory key variable should be named "IID'
# All people listed in ids_order_file should be present in original_phenofile
#

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Sytax is Rscript 02.preparePheno_protein_raw.r transformation(log/inv) protein visit(3/5)", call.=FALSE)
 }

transformation<-args[1]
protien<-args[2]
protein<-paste0(transformation,protien)
visit<-args[3]


original_phenofile <- "/dcs01/arking/arkinglab/static/aric/phenotypes/f3v2/aric.allphenotypes.trim.final.v4.f3v2.txt"
pcs_file <-"../datasets/pc.txt"

# Files contating the SOMA protein and ARIC visit 3 covariates. Made in $ARIC_SOMA/users/asurapa2/04
if (visit==3){
aric_covariatesfile <- 	"../datasets/ARICv3_covariates.csv"
aric_proteinsfile   <-	"../datasets/SOMAv3_proteins.csv"
}
if (visit==5){
aric_covariatesfile <- "../datasets/ARICv5_covariates.csv"
aric_proteinsfile   <- "../datasets/SOMAv5_proteins.csv"
}

#
# IDs present in TOPMed Imputation files
# IDs are in order as presented in .dose.vcf.gz files
#
ids_order_blacks <- "/dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/B_ARIC_topmed.imputation.iid"
ids_order_whites <- "/dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/W_ARIC_topmed.imputation.iid"

#
# Output file names (the one used by linear/logistic regression)
#
output_phenofile.W.iid <- paste0(protein,".visit",visit,".W.FAST.runfile.iid.txt")
output_phenofile.B.iid <- paste0(protein,".visit",visit,".B.FAST.runfile.iid.txt")

#
# Trait to be analyzed; if running CoxPH, format is (trait, time to event)
#

trait <- c(protein)

#
# Covariates to be included into model
# Adjusting for center??
covars <- c("AGE", "PC1", "PC2", "PC3", "PC4","PC5", "PC6", "PC7", "PC8", "PC9", "PC10","Forsyth","Minneapolis","Washington","gPC1", "gPC2", "gPC3", "gPC4","gPC5","gPC6","gPC7","gPC8","gPC9","gPC10","sex")

#
# Directory to save file
#
#directory = "/dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/GEH/spatial/raw/"
directory =  "/users/asurapan/04_ARIC.GWAS.with.FAST/"
directory_W = paste0("../results/" , protein, "_Visit", visit, "_W/")
directory_B = paste0("../results/" , protein, "_Visit", visit, "_B/")

#
# Code
#

##################
#Reading in files#
##################

aric_covars<-read.csv(aric_covariatesfile)

library(data.table)
setnames(aric_covars, "id", "IID")
setnames(aric_covars, "age", "AGE")


aric_proteins<-read.csv(aric_proteinsfile)
setnames(aric_proteins, "SampleId", "IID")

merged = merge(aric_proteins, aric_covars, by = "IID")
rownames(merged) <- merged$IID


#####################
#Phenotype exclusion#
#####################

# What is this exclusion for? There was something to do with bad genetic data that needed to be excluded..
#merged$QRS_TmV1<-ifelse(merged$exclude==1, NA, merged$QRS_TmV1) 

####################################
#Add paternal and maternal columns.#
####################################
#Phenotype file does not have information about this variable so just code it all as 0

ophe <- read.table(original_phenofile,head=T,strings=F)
ophe <- ophe[,c('FID','IID','GENDER','CENTER',"PC1", "PC2", "PC3", "PC4","PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
names(ophe) <- c('FID','IID','GENDER','CENTER',"gPC1", "gPC2", "gPC3", "gPC4","gPC5", "gPC6", "gPC7", "gPC8", "gPC9", "gPC10")


pcs <- read.table(pcs_file, head=T, strings=F)
setnames(pcs,"SampleId","IID")
ophe <-merge(ophe, pcs, by="IID", all.x=TRUE)


ophe$Forsyth     <-ifelse(ophe$CENTER=="F",1,0)
ophe$Minneapolis <-ifelse(ophe$CENTER=="M",1,0)
ophe$Washington  <-ifelse(ophe$CENTER=="W",1,0)


rownames(ophe) <- ophe$IID

merged<-merge(ophe,merged,by="IID", all.x=TRUE)
rownames(merged) <- merged$IID

merged$FAT_ID<-0
merged$MAT_ID<-0

#################
#Fix the centers#
#################

#There are three centers in the whites. Add a new CENTER variables to account for that
#merged$CENTERM<-ifelse(merged$CENTER=="M",1,0)

###################################
#NA samples with no BMI and Height#
###################################
#merged$QRS_TmV1<-ifelse(is.na(merged$height01), NA, merged$QRS_TmV1)

################################
#Write out a new phenotype file#
################################

ids <- scan(ids_order_whites,what=character())
print(ids)
head(merged)
merged$sex<-merged$GENDER

newphe.iid <- merged[ids, c("FID", "IID", "FAT_ID", "MAT_ID", "GENDER", trait, covars)]
names(newphe.iid)[1]<-"#FID"
system(paste0("mkdir -p ",directory_W))
#write.table(newphe.iid, file=paste0(directory_W, output_phenofile.W.iid), quote=F, row.names=F, col.names=T, sep="\t")
summary(newphe.iid)

#ids <- scan(ids_order_blacks,what=character())
#newphe.iid <- merged[ids, c("FID", "IID", "FAT_ID", "MAT_ID", "GENDER", trait, covars)]
#names(newphe.iid)[1]<-"#FID"
#system(paste0("mkdir ",directory_B))
#write.table(newphe.iid, file=paste0(directory_B, output_phenofile.B.iid), quote=F, row.names=F, col.names=T, sep="\t")
#summary(newphe.iid)


#tail(newphe.iid)
#tail(ids)
#cat("\n\n\n")
#summary(newphe.iid)
names(newphe.iid)
model= paste(c(trait,paste(covars,collapse="+")), collapse=" ~ ")
fit=lm(as.formula(model), data=newphe.iid, na.action="na.exclude")
newphe.iid$res=residuals(fit)

newphe.iid[[protein]]=newphe.iid$res
#newphe.iid=newphe.iid[,!names(newphe.iid) %in% c("res")]
newphe.iid <- newphe.iid[, c(c("#FID", "IID", "FAT_ID", "MAT_ID", "GENDER"),protein)]
summary(newphe.iid)

write.table(newphe.iid, file=paste0(directory_W, output_phenofile.W.iid), quote=F, row.names=F, col.names=T, sep="\t")



#newphe.iid$invres=qnorm((rank(newphe.iid$res)-0.5)/sum(!is.na(newphe.iid$res)))
#library(RNOmni)
#newphe.iid$invres=rankNorm(newphe.iid$res)

#tmp=newphe.iid[!is.na(newphe.iid$res),]

#tmp$invres=rankNorm(tmp$res)
#tmp$invres2=qnorm((rank(tmp$res)-0.5)/sum(!is.na(tmp$res)))
#tmp=tmp[c("#FID","IID","res")]
#write.table(tmp, file=paste0(directory_W, "invres",output_phenofile.W.iid,"tmp"), quote=F, row.names=F, col.names=T, sep="\t")

#for (var in covars){
	#cat(var)
#	newphe.iid[[var]]<-0
#}

#newphe.iid[[protein]]=newphe.iid$invres
#newphe.iid=newphe.iid[,!names(newphe.iid) %in% c("res","invres")]
#summary(newphe.iid)
#write.table(newphe.iid, file=paste0(directory_W, "invres",output_phenofile.W.iid), quote=F, row.names=F, col.names=T, sep="\t")


