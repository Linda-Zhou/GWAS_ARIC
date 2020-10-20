#args = commandArgs(trailingOnly=TRUE)
#filename<-args[1]

#QUAL filter to decrease # of SNPs read into R
x="awk -v OFS='\\t\' '{if ($9>0.3) {print $0,$2\":\"$3} }'\ allchr.FAST.txt > allchr.FAST.tmp\n"
system(x)

library(data.table)
#filename<-"/dcs01/arking/arkinglab/active/projects/aric/analyses/f3v2.1000G.phase3/qrs.dur.blacks/raw/allchr.FAST.txt"
data<-fread("allchr.FAST.tmp",header=T)
data<-na.omit(data)
names(data)[ncol(data)] = "chrpos"
#data$chrpos<-paste0(data$Chr, data$Position)

#Remove problematic SNPs
snps<-fread("/dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/scripts/topmed.blacks.snps.to.remove.chrpos.txt",header=F)
trim<-subset(data, !(chrpos %in% snps$V1))
trim<-subset(trim, CAF>=0.1 & CAF <=0.99)
#trim<-subset(trim, Qual>0.3)
write.table(trim,"allchr.FAST.cleanSNPs.txt",sep="\t",row.names=FALSE,quote=FALSE)

system("rm allchr.FAST.tmp")
