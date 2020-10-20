#TO RUN SCRIPT: Rscript TOPMed.GWAS.qc.R trait file > TOPMed.GWAS.qc.Rout#
#SCRIPT MUST BE RUN INTERACTIVELY#
install.packages("sqldf",repos="http://cran.r-project.org")
args = commandArgs(trailingOnly=TRUE)
print(length(args))
if (length(args)!=3) {
  stop("Syntax: Rscript script transformation protein visit", call.=FALSE)
 }
transformation<-args[1]
trait<-args[2]
visit<-args[3]
directory<-paste0("../results/",transformation,trait,"_Visit",visit,"_W/")
full_filename<-paste0("../results/",transformation,trait,"_Visit",visit,"_W/allchr.FAST.txt")
full_filename_trim<-paste0(full_filename,".trim")

if (1==1){

#BEGIN QQ PLOT Code


##trim results file by using every 20th SNP##
#system("sed 's/ \\+/\\t/g' allchr.FAST.txt >allchr.FAST.txt.fix")
#system("mv allchr.FAST.txt.fix allchr.FAST.txt")
#system("sed -n '1~20p' allchr.FAST.txt > allchr.FAST.txt.trim")

system(paste0("sed -n '1~200p' ", full_filename, " > ",full_filename_trim)) 

#to speed up reading gwas input, use sqldf
library(sqldf)
f <- file(full_filename_trim)
data <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = T, row.names = F, sep="\t"))
data<-na.omit(data)
z.sq<-(data$Beta / data$SE)^2
lambda<-median(z.sq,na.rm=TRUE) / 0.456
cat(paste("Lambda =",lambda,"\n"))
##working from FAST output###
data$Pvalue<-data$P
data$imp_quality<-data$Qual

##qq maf plot##
png(paste(directory,transformation,trait,"Visit_",visit,".qqplot.maf.png",sep=""),res=400,units="cm",width = 16, height = 16, type='cairo')
s<-paste("lambda=",lambda,sep="")
pvals <- data$Pvalue
observed <- sort(pvals)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))

plot(c(0,9), c(0,9), col="red", lwd=4, type="l", xlab="Expected (-logP)",
	ylab="Observed (-logP)", xlim=c(0,9), ylim=c(0,9), las=1, xaxs="i", yaxs="i", bty="l",
	sub=s, main=trait)
points(lexp, lobs, pch=23, cex=.5, col="black", bg="black")

##maf005##
if (min(data$CAF) <0.005 ) {
maf<-subset(data,CAF<0.005 | CAF >0.995)
pvals<-maf$Pvalue
observed <- sort(pvals)
nmaf005<-length(observed)
nmaf005<-paste("MAF <0.005 [",nmaf005,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="brown", bg="brown")
}
##maf01##
if (min(data$CAF) <0.01 ) {
maf<-subset(data,CAF<0.01 & CAF >=0.005 | CAF >0.99 & CAF <=0.995)
pvals<-maf$Pvalue
observed <- sort(pvals)
nmaf01<-length(observed)
nmaf01<-paste("0.005>= MAF <0.01 [",nmaf01,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="purple", bg="purple")
}
##maf05##
maf<-subset(data,(CAF>=0.01 & CAF <0.05) | (CAF<=0.99 & CAF >0.95))
pvals<-maf$Pvalue
observed <- sort(pvals)
nmaf05<-length(observed)
nmaf05<-paste("0.01>= MAF <0.05 [",nmaf05,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="blue", bg="blue")
##maf05-20##
maf<-subset(data,(CAF>=0.05 & CAF <0.2) | (CAF<=0.95 & CAF >0.8))
pvals<-maf$Pvalue
observed <- sort(pvals)
nmaf10<-length(observed)
nmaf10<-paste("0.05>= MAF <0.20 [",nmaf10,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="orange", bg="orange")
##maf20##
maf<-subset(data,CAF>=0.2 & CAF <=0.8)
pvals<-maf$Pvalue
observed <- sort(pvals)
nmaf20<-length(observed)
nmaf20<-paste("MAF >=0.20 [",nmaf20,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="green", bg="green")
legend("bottomright", c("Expected","Observed",nmaf20,nmaf10,nmaf05), cex=0.8, col=c("red","black","green","orange","blue","purple","brown"), pch=18)
dev.off()


##qq imp plot##
png(paste(directory,transformation,trait,"Visit_",visit,".qqplot.imp.png",sep=""),res=400,units="cm",width = 16, height = 16, type='cairo')
s<-paste("lambda=",lambda,sep="")
pvals <- data$Pvalue
observed <- sort(pvals)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))

plot(c(0,9), c(0,9), col="red", lwd=4, type="l", xlab="Expected (-logP)",
	ylab="Observed (-logP)", xlim=c(0,9), ylim=c(0,9), las=1, xaxs="i", yaxs="i", bty="l",
	sub=s, main=trait)
points(lexp, lobs, pch=23, cex=.5, col="black", bg="black")
##impqual0.25##
if (min(data$imp_quality, na.rm=T)<0.25) {
maf<-subset(data,imp_quality<0.25)
pvals<-maf$Pvalue
observed <- sort(pvals)
nimp25<-length(observed)
nimp25<-paste("ImpQual <0.25 [",nimp25,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="purple", bg="purple")
}
##impqual0.50##
if (min(data$imp_quality, na.rm=T)<0.5) {
maf<-subset(data,imp_quality>=0.25 & imp_quality <0.5)
pvals<-maf$Pvalue
observed <- sort(pvals)
nimp50<-length(observed)
nimp50<-paste("0.25>= ImpQual <0.50 [",nimp50,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="blue", bg="blue")
}
##impqual0.75##
maf<-subset(data,imp_quality>=0.5 & imp_quality <0.75)
pvals<-maf$Pvalue
observed <- sort(pvals)
nimp75<-length(observed)
nimp75<-paste("0.50>= ImpQual <0.75 [",nimp75,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="orange", bg="orange")
##impqual0.75+##
maf<-subset(data,imp_quality>=0.75)
pvals<-maf$Pvalue
observed <- sort(pvals)
nimp90<-length(observed)
nimp90<-paste("ImpQual >0.75 [",nimp90,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="green", bg="green")
legend("bottomright", c("Expected","Observed",nimp90,nimp75,nimp50,nimp25), cex=0.8, col=c("red","black","green","orange","blue","purple"), pch=18)
dev.off()

##qq eSampleSize plot##
png(paste(directory,transformation,trait,"Visit_",visit,".qqplot.effectiveSampleSize.png",sep=""),res=400,units="cm",width = 16, height = 16, type='cairo')
s<-paste("lambda=",lambda,sep="")
pvals <- data$Pvalue
observed <- sort(pvals)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))

plot(c(0,9), c(0,9), col="red", lwd=4, type="l", xlab="Expected (-logP)",
	ylab="Observed (-logP)", xlim=c(0,9), ylim=c(0,9), las=1, xaxs="i", yaxs="i", bty="l",
	sub=s, main=trait)
points(lexp, lobs, pch=23, cex=.5, col="black", bg="black")
##ESS5##
if (min(data$ESampleSize) <5) {
maf<-subset(data,ESampleSize<5)
pvals<-maf$Pvalue
observed <- sort(pvals)
nmaf005<-length(observed)
nmaf005<-paste("ESampleSize <5 [",nmaf005,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="brown", bg="brown")
}
##ESS10##
if (min(data$ESampleSize) <10 ) {
maf<-subset(data,ESampleSize<10 & ESampleSize >5)
pvals<-maf$Pvalue
observed <- sort(pvals)
nmaf01<-length(observed)
nmaf01<-paste("5>= ESampleSize <10 [",nmaf01,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="purple", bg="purple")
}
#ESS25##
maf<-subset(data,(ESampleSize>=10 & ESampleSize <25))
pvals<-maf$Pvalue
observed <- sort(pvals)
nmaf05<-length(observed)
nmaf05<-paste("10>= ESampleSize <25 [",nmaf05,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="blue", bg="blue")
##ESS50##
maf<-subset(data,(ESampleSize>=25 & ESampleSize <50))
pvals<-maf$Pvalue
observed <- sort(pvals)
nmaf10<-length(observed)
nmaf10<-paste("25>= ESampleSize <50 [",nmaf10,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="orange", bg="orange")
##ESS100##
maf<-subset(data,ESampleSize>=50)
pvals<-maf$Pvalue
observed <- sort(pvals)
nmaf20<-length(observed)
nmaf20<-paste("ESampleSize >=50 [",nmaf20,"]",sep="")
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
points(lexp, lobs, pch=23, cex=.5, col="green", bg="green")
legend("bottomright",c("Expected","Observed",nmaf20,nmaf10,nmaf05,nmaf01,nmaf005), cex=0.8, col=c("red","black","green","orange","blue","purple","brown"), pch=18)
dev.off()



}

# END QQ Plot code





library(sqldf)
library(data.table)

####manhattan plot for all data####

# CAF between 0.01 and 0.99. Imputation quality >0.1
system(paste0("cat ",full_filename, " |awk -F'\t' '$8 > 0.01 && $8 <0.99 && $9 >0.01 && $11<0.01 { print }' > " ,full_filename_trim))
#system(paste0("cat ",full_filename, " |awk -F'\t' '$11<0.0001 { print }' > " ,full_filename_trim))

system(paste0("head -1 ",full_filename, " > ",full_filename,".head"))
system(paste0("cat ",full_filename,".head ",full_filename,".trim > ",full_filename,".trim.wheader"))



f <- file((paste0(full_filename,".trim.wheader")))

if (1==1){
data <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = T, row.names = F, sep="\t"))

#data<-fread(paste0(full_filename,".trim.wheader"))
data<-na.omit(data)


setnames(data, "P", "Pvalue")
setnames(data, "Qual", "imp_quality")

yminl<- -log(min(data$Pvalue),base=10)+1


###output file with best hits###
hits<-data
hits<-subset(data,Pvalue<=0.00000001)
#hits<-hits[ , !(names(hits) %in% c("P","Qual"))]
hits<-hits[order(hits$Pvalue),] 
#hits<-head(hits,50)

rsid_from_snp <- function(SNP) {
	rsid=""
	try(snp_rsid<-system(paste0("grep ",SNP," /dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/rsid/GRCh38_dbSNP151_rsid_final_USE.txt "), intern=TRUE))
        #cat(snp_rsid[1],"\n")
        try(snp_rsid<-strsplit(snp_rsid, "\t")[1])
	try(rsid<-as.vector(snp_rsid)[[1]][2])
        return(rsid)
}

try(hits$rsid <- mapply(rsid_from_snp,hits$SNP))

hits$len<- mapply(length,hits$rsid)
hits<-hits[hits$len>0,]
hits$rsid<- as.character(hits$rsid)

write.table(hits, paste0(directory,transformation,trait,"_Visit",visit,"_hits.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)



#chr <- c(1:22)
chr<-c(1:max(data$Chr))

#Summary statistics
data$position<-round(data$Position/1000000,digits=0)
#print(summary(data))
print(table(data$Chr))
par(mar=c(5,5,4,2))
phy.max<-tapply(data$position, data$Chr,max,na.rm=T)
cumlen=0
for(i in chr){
	cat(paste("Now working on chromosome ",i,"\r"))
	data[data$Chr==i,"loc"]<-data[data$Chr==i,"position"]+cumlen
	cumlen<-cumlen+phy.max[i]
}

phy.med<-tapply(data$loc,data$Chr,median,na.rm=T)[chr]
print(phy.med)
data$mlgpval<--log(data[,"Pvalue"], base=10)
png(paste(directory,transformation,trait,"Visit_",visit,".manhattanplot.png",sep=""),res=400,width = 32, height = 16,units="cm",type='cairo')
plot(data[,"loc"],data[,"mlgpval"],type="n",yaxt="n",xaxt="n",xlab="Chromosome",
	ylab=expression(-log[10]* (P)),main=trait,
	xlim=c(0,max(data$loc,na.rm=T)),ylim=c(0,yminl))
col=rep(c("blue","red"),13)
axis(side=2, at=seq(from=0,to=yminl,by=1), labels=seq(from=0,to=yminl,by=1),
	tick=T,cex.axis=0.5,las=1)
axis(side=1, at=phy.med[c(1:max(data$Chr))], labels=chr[c(1:max(data$Chr))],
	tick=T,cex.axis=0.5,las=1)
#axis(side=1, at=phy.med[c(20:23)], labels=chr[c(20:23)],
#	tick=T,cex.axis=1,las=1)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
for(i in chr){
	cat(paste("Now working on chromosome ",i,"\r"))
	#if(which(chr==i) %in% seq(2,30,2)) col="blue" else col="red"
	points(data[data$Chr==i & data$ESampleSize<10,"loc"],data[data$Chr==i & data$ESampleSize<10,"mlgpval"],col="gray70",pch=20,cex=0.8)
	points(data[data$Chr==i & data$ESampleSize<25 & data$ESampleSize>=10,"loc"],data[data$Chr==i & data$ESampleSize<25 & data$ESampleSize>=10,"mlgpval"],col="gray40",pch=20,cex=0.8)
	points(data[data$Chr==i & data$ESampleSize>=25,"loc"],data[data$Chr==i & data$ESampleSize>=25,"mlgpval"],col=col[i],pch=20,cex=0.8)
}
legend(0.3,0.3,c("ESampleSize <10","10>= ESampleSize <25","ESampleSize >=25","ESampleSize >=25"), cex=0.8, col=c("gray70","gray40","blue","red"),horiz=T, pch=18)
#adding broken lines
#for (i in 0:9){
	abline(h=7.3,lty="dotted",col="dark grey")
#}
dev.off()

}

system(paste0("rm ",full_filename,".trim "))
system(paste0("rm ",full_filename,".head "))
system(paste0("mv ",full_filename,".trim.wheader ",full_filename,".trim.wheader.txt"))
#system(paste0("rm ",full_filename))

system(paste0("rm Rplots.pdf"))


