#module load conda_R
#if (!require("devtools")) { install.packages("devtools") } else {}
#devtools::install_github("rondolab/MR-PRESSO")
#install.packages("MRCIEU/TwoSampleMR") 
#devtools::install_github("MRCIEU/TwoSampleMR") #to update the package
#devtools::install_github("MRCIEU/MRInstruments")
#install_github("WSpiller/RadialMR")
#devtools::install_github("gqi/MRMix")
library(dplyr)
#library(ggrepel)
#library(ggthemes)
library(devtools)
library(TwoSampleMR)
#library(MRInstruments)
#library(MRPRESSO)
#library(MRMix)
#library(ggplot2)
#library(png)
#library(psych)
#library(knitr)
#library(ggrepel)
set.seed(629)
library(data.table)
library(readr)
library(haven)
library(stringr)

#[compute-062 /dcl02/leased/kidney/CRIC/projects/04_GWAS_Metabolites/results/black/logKDComp100_pos_asgpcmpc_black]$ head allchr.FAST.hits.txt
#SNPID   Chr     Position        non_coded_allele        coded_allele    Beta    SE      CAF     Qual    ESampleSize     P

transformation <-"log"
metabolite <-"KDComp100"
adj        <-"asgpcmpc"
race       <-"white"
dataset    <-"c8pos"




#annot <- read_stata("/dcl02/leased/kidney/CRIC/static/Metabolomics/Original_Data/CKD_HILIC-pos_plasma_Metabolite_FullInfor.dta")
#annot <- read_stata("/dcl02/leased/kidney/CRIC/static/Metabolomics/Original_Data/CKD_C8-pos_plasma_Metabolite_FullInfor.dta")
morgan_request <- read.table('morgan_request.sh')

#KN=annot[annot$mettype==1| annot$mettype==2,]$metid

proteins=morgan_request$V4
proteins=as.vector(proteins)

proteins<-c("ME0_v5","ME1_v5","ME2_v5","ME3_v5")

#KN <- str_replace(KN,"_","")
#KN <-sort(KN)

#KN_1 <-c("KNComp213")
#KN_2 <-tail(KN,90)

#print(KN)
print(proteins)
#print(KN_2)

#stop("")
for (protein in proteins) {
 print(protein)
 #system("sleep 1s")

try({

#filename <- paste0("/dcl02/leased/kidney/CRIC/projects/04_GWAS_Metabolites/results/",race,"/",transformation,metabolite,"_",dataset,"_",adj,"_",race,"/",transformation,metabolite,"_",adj,"_",race,"_hits")
filename <- paste0("/dcl02/leased/kidney/ARIC/projects/01_GWAS_Proteins/results/",transformation,protein,"_Visit5_W/",transformation,protein,"_Visit5_hits")
print(filename)
dt <- read.table(file=paste0(filename,".txt"), header = TRUE, sep='\t', stringsAsFactors = FALSE)
print(dim(dt))
print(names(dt))

#dt<-dt[dt$ESampleSize>25,]

print(dim(dt))

#stop("")
exposure <- dt %>% 
            dplyr::select(c("rsid","Chr","Position","Beta","SE","CAF","coded_allele","non_coded_allele","Pvalue")) %>%
            dplyr::rename("SNP"=rsid,"eaf"=CAF,"effect_allele"=coded_allele,"other_allele"=non_coded_allele, "pvalue"=Pvalue, "beta"=Beta, "se"=SE) 

exp_data <- format_data(exposure, type="exposure")
exp_data$exposure <- "SNP"
exp_data$chr_name=exp_data$chr.exposure
exp_data$chrom_start=exp_data$pos.exposure

                
# LD pruning
exp_data_pruned <- clump_data(exp_data, pop="EUR",clump_p1=0.00000005, clump_p2=0.01)
exp_data_pruned$rsid <-exp_data_pruned[['SNP']]

exp_data_pruned <-exp_data_pruned[c('rsid')]

dt<-left_join(exp_data_pruned,dt)        
print(dim(dt))
write.table(dt,paste0(filename,"_clumped.txt"), quote=FALSE, row.names=FALSE, sep="\t")

})
}
