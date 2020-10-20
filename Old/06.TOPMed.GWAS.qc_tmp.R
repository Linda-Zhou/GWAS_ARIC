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

library(data.table)

library(sqldf)

f <- file((paste0(full_filename,".trim.wheader.txt")))

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
png(paste(directory,transformation,trait,"Visit_",visit,".manhattanplot.png",sep=""),res=100,width = 8, height = 4,units="cm",type='cairo')
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


system(paste0("rm Rplots.pdf"))


