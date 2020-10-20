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
original_phenofile <- "/dcs01/arking/arkinglab/static/aric/phenotypes/f3v2/aric.allphenotypes.trim.final.v4.f3v2.txt"

#
# Most recent QRS duration data provided by Alvaro on 11.29.2017 combined with LBBB and RBBB exclusion criteria (generated 09.24.2018 by ThuyVy)
# Individuals are flagged based on exclusion criteria listed in QT/JT/QRS 2018 Analysis Plan
#
updated_qrs <- "/dcs01/arking/arkinglab/static/aric/phenotypes/traits/qt/exclusion_qt.jt.qrs_20180924.txt"

#
# IDs present in TOPMed Imputation files
# IDs are in order as presented in .dose.vcf.gz files
#
ids_order_file <- "/dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation"

#
# Output file names (the one used by linear/logistic regression)
#
output_phenofile.iid <- "FAST.runfile.iid.txt"

#
# Trait to be analyzed; if running CoxPH, format is (trait, time to event)
#
trait <- c("QRS_TmV1")

#
# Covariates to be included into model
#
covars <- c("V1AGE01", "GENDER", "RRInterval", "QT_Dur", "height01", "BMI01", "CENTER", "PC1", "PC2", "PC3", "PC4")

#
# Directory to save file
#
directory = "/dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/GEH/spatial/raw/"

#
# Code
#

##################
#Reading in files#
##################

#Extract R-R Interval and QT Interval
ecg.data = "/dcs01/arking/arkinglab/static/aric/phenotypes/raw.data/visit1/v1ecg.sas7bdat"
library(haven)
ecg = read_sas(ecg.data)
HR = as.data.frame(ecg[,c("ID", "v1ecg6", "v1ecg9")])
names(HR)[c(2,3)] = c("HeartRate", "QT_Dur")
HR$RRInterval = (60/ecg$v1ecg6)*1000

ophe <- read.table(original_phenofile,head=T,strings=F)
rownames(ophe) <- ophe$IID
names(ophe)[c(67,72)] = c("qrs_dur.old", "qrs_axis.old")
ids <- scan(ids_order_file,what=character())
new.qrsphe <- read.table(updated_qrs, head = T, strings = F, sep = "\t")
new.qrsphe <- merge(new.qrsphe, HR, by = "ID", all.x = T, sort = F)
names(new.qrsphe)[1:3] = c("IID", "qrs_dur", "qrs_axis")
rownames(new.qrsphe) <- new.qrsphe$IID

merged = merge(ophe, new.qrsphe, by = "IID")
rownames(merged) <- merged$IID

#Extract spatial QRS-T angle phenotype
geh.pheno<-read.table("/dcs01/arking/arkinglab/static/aric/analysis.plans/GEH/ARIC GEH visit1.txt",header=T)
geh.key<-read.table("/dcs01/arking/arkinglab/static/aric/analysis.plans/GEH/ARIC_IDs_visit1_fromCC.txt",header=T)
geh.key$id_n<-as.factor(geh.key$id_n)
geh.pheno$id_n<-as.factor(geh.pheno$id_n)
geh<-merge(geh.key,geh.pheno,by="id_n")
names(geh)[2]<-"IID"

merged = merge(merged, geh, by = "IID", all.x = TRUE)
rownames(merged) <- merged$IID

#####################
#Phenotype exclusion#
#####################
merged[, 94:105][is.na(merged[, 94:105])] <- 0
merged$QRS_TmV1<-ifelse(merged$afib_v1==1, NA, merged$QRS_TmV1)
merged$QRS_TmV1<-ifelse(merged$long_qrs==1, NA, merged$QRS_TmV1)
merged$QRS_TmV1<-ifelse(merged$comp_lbbb==1, NA, merged$QRS_TmV1)
merged$QRS_TmV1<-ifelse(merged$comp_rbbb==1, NA, merged$QRS_TmV1)
merged$QRS_TmV1<-ifelse(merged$pacemaker==1, NA, merged$QRS_TmV1)
merged$QRS_TmV1<-ifelse(merged$wpw.y==1, NA, merged$QRS_TmV1)
#merged$QRS_TmV1<-ifelse(merged$avblock==1, NA, merged$QRS_TmV1)
merged$QRS_TmV1<-ifelse(merged$PREVHF01.y==1, NA, merged$QRS_TmV1)
merged$QRS_TmV1<-ifelse(merged$PREVMI05.y==1, NA, merged$QRS_TmV1)
merged$QRS_TmV1<-ifelse(merged$antiarrh==1, NA, merged$QRS_TmV1)
merged$QRS_TmV1<-ifelse(merged$meds==1, NA, merged$QRS_TmV1)
merged$QRS_TmV1<-ifelse(merged$exclude==1, NA, merged$QRS_TmV1)

####################################
#Add paternal and maternal columns.#
####################################
#Phenotype file does not have information about this variable so just code it all as 0
merged$FAT_ID<-0
merged$MAT_ID<-0

#################
#Fix the centers#
#################
#There are three centers in the whites. Add a new CENTER variables to account for that
merged$CENTERM<-ifelse(merged$CENTER=="M",1,0)
merged$CENTER<-ifelse(merged$CENTER=="F",1,0)

###################################
#NA samples with no BMI and Height#
###################################
merged$QRS_TmV1<-ifelse(is.na(merged$BMI01), NA, merged$QRS_TmV1)
merged$QRS_TmV1<-ifelse(is.na(merged$height01), NA, merged$QRS_TmV1)

################################
#Write out a new phenotype file#
################################
newphe.iid <- merged[ids, c("FID", "IID", "FAT_ID", "MAT_ID", "GENDER", trait, covars, "CENTERM")]
names(newphe.iid)[1]<-"#FID"
write.table(newphe.iid, file=paste0(directory, output_phenofile.iid), quote=F, row.names=F, col.names=T, sep="\t")
