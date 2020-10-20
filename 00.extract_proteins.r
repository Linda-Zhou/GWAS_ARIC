rm(list=ls())

# 1. Read in Visit 3 and 5 protein data
f.cohort_v3 = paste("/dcl02/leased/kidney/ARIC/static/Proteomics/Data/", "soma_visit_3_log2_SMP.txt", sep="")
d.cohort_v3 = read.table(f.cohort_v3, header=T, as.is=T, sep="\t")

f.cohort_v5 = paste( "/dcl02/leased/kidney/ARIC/static/Proteomics/Data/", "soma_visit_5_log2_SMP.txt", sep="")
d.cohort_v5 = read.table(f.cohort_v5, header=T, as.is=T, sep="\t")


# 2. Read in the annotation data
#install.packages("readxl") # CRAN version
library(readxl)
f.annot = paste("../datasets/" , "Abbreviated annotation visits 3 and 5.xlsx", sep="")
d.annot=read_excel(f.annot, sheet = "Annotation")


# 3. Select proteins of interest from the annotation data.
d.annot_proteins=subset(d.annot,grepl("seqid_10036_201", tolower(seqid_in_sample)) |
                                grepl("seqid_10440_26",  tolower(seqid_in_sample)) |
                                grepl("seqid_10514_5",   tolower(seqid_in_sample)) |
                                grepl("seqid_11388_75",  tolower(seqid_in_sample)) | 
                                grepl("seqid_13126_52",  tolower(seqid_in_sample)) |
                                grepl("seqid_2654_19",   tolower(seqid_in_sample)) |
                                grepl("seqid_2944_66",   tolower(seqid_in_sample)) |
                                grepl("seqid_3152_57",   tolower(seqid_in_sample)) |
                                grepl("seqid_3438_10",   tolower(seqid_in_sample)) |
                                grepl("seqid_4721_54",   tolower(seqid_in_sample)) |
                                grepl("seqid_8323_163",  tolower(seqid_in_sample)) |
                                grepl("seqid_8368_102",  tolower(seqid_in_sample)) |
                                grepl("seqid_9359_9",    tolower(seqid_in_sample)) |
                                grepl("seqid_5879_51",    tolower(seqid_in_sample)) |
                                grepl("seqid_9468_8",    tolower(seqid_in_sample)) )



#d.annot_proteins=d.annot_proteins[c('seqid_in_sample','target','targetfullname','entrezgenesymbol')]
#View(d.annot_proteins)
seq_ids=as.vector(d.annot_proteins$seqid_in_sample)
write.csv(d.annot_proteins,"Selected_proteins.csv")

protein_ids=as.vector(d.annot_proteins$entrezgenesymbol)

rm(list=ls())

# 1. Read in Visit 3 and 5 protein data
f.cohort_v3 = paste("O:/Data/Visit_3/", "soma_visit_3_log2_SMP.txt", sep="")
d.cohort_v3 = read.table(f.cohort_v3, header=T, as.is=T, sep="\t")

f.cohort_v5 = paste( "O:/Data/Visit_3/", "soma_visit_5_log2_SMP.txt", sep="")
d.cohort_v5 = read.table(f.cohort_v5, header=T, as.is=T, sep="\t")


# 2. Read in the annotation data
#install.packages("readxl") # CRAN version
library(readxl)
#f.annot = paste("O:/Annotation/" , "Abbreviated annotation visits 3 and 5.xlsx", sep="")
f.annot = "/dcl02/leased/kidney/ARIC/static/Proteomics/Annotation/Complete annotation visit 3.txt"
d.annot = read.table(f.annot, header=T, as.is=T, sep="\t")


# 3. Select proteins of interest from the annotation data.
d.annot_proteins=subset(d.annot,grepl("seqid_10036_201", tolower(seqid_in_sample)) |
                                grepl("seqid_10440_26",  tolower(seqid_in_sample)) |
                                grepl("seqid_10514_5",   tolower(seqid_in_sample)) |
                                grepl("seqid_11388_75",  tolower(seqid_in_sample)) | 
                                grepl("seqid_13126_52",  tolower(seqid_in_sample)) |
                                grepl("seqid_2654_19",   tolower(seqid_in_sample)) |
                                grepl("seqid_2944_66",   tolower(seqid_in_sample)) |
                                grepl("seqid_3152_57",   tolower(seqid_in_sample)) |
                                grepl("seqid_3438_10",   tolower(seqid_in_sample)) |
                                grepl("seqid_4721_54",   tolower(seqid_in_sample)) |
                                grepl("seqid_8323_163",  tolower(seqid_in_sample)) |
                                grepl("seqid_8368_102",  tolower(seqid_in_sample)) |
                                grepl("seqid_9359_9",    tolower(seqid_in_sample)) |
                                grepl("seqid_5879_51",    tolower(seqid_in_sample))|
                                grepl("seqid_9468_8",    tolower(seqid_in_sample)) )



d.annot_proteins=d.annot_proteins[c('seqid_in_sample','target','targetfullname','entrezgenesymbol')]
View(d.annot_proteins)
seq_ids=as.vector(d.annot_proteins$seqid_in_sample)
write.csv(d.annot_proteins,"Selected_proteins.csv")
protein_ids=as.vector(d.annot_proteins$entrezgenesymbol)

n=length(protein_ids)
for (i in 1:n){
  protein_ids[i]=paste0(protein_ids[i],"_",seq_ids[i])
  
}

# 4. Select proteins of interest from the SOMA data. Inverse normal transform.
# (optionally rename proteins from seq_ids)
d.cohort_v5_proteins=d.cohort_v5[c("SampleId",seq_ids)]
names(d.cohort_v5_proteins)=c("SampleId",protein_ids)

for (name in protein_ids) {
  d.cohort_v5_proteins$inv=qnorm((rank(d.cohort_v5_proteins[name],na.last="keep")-0.5)/sum(!is.na(d.cohort_v5_proteins[name])))
  colnames(d.cohort_v5_proteins)[colnames(d.cohort_v5_proteins) == name] <-paste0('log',name)
  colnames(d.cohort_v5_proteins)[colnames(d.cohort_v5_proteins) == 'inv'] <-paste0('inv',name)
}

d.cohort_v3_proteins=d.cohort_v3[c("SampleId",seq_ids)]
names(d.cohort_v3_proteins)=c("SampleId",protein_ids)
for (name in protein_ids) {
  d.cohort_v3_proteins$inv=qnorm((rank(d.cohort_v3_proteins[name],na.last="keep")-0.5)/sum(!is.na(d.cohort_v3_proteins[name])))
  colnames(d.cohort_v3_proteins)[colnames(d.cohort_v3_proteins) == name] <-paste0('log',name)
  colnames(d.cohort_v3_proteins)[colnames(d.cohort_v3_proteins) == 'inv'] <-paste0('inv',name)
}


# 5. Export data as csv
write.csv(d.cohort_v5_proteins,"SOMAv5_proteins.csv")
write.csv(d.cohort_v3_proteins,"SOMAv3_proteins.csv")

# 4. Select proteins of interest from the SOMA data. Inverse normal transform.
# (optionally rename proteins from seq_ids)
d.cohort_v5_proteins=d.cohort_v5[c("SampleId",seq_ids)]
names(d.cohort_v5_proteins)=c("SampleId",protein_ids)
for (name in protein_ids) {
  d.cohort_v5_proteins$inv=qnorm((rank(d.cohort_v5_proteins[name],na.last="keep")-0.5)/sum(!is.na(d.cohort_v5_proteins[name])))
  colnames(d.cohort_v5_proteins)[colnames(d.cohort_v5_proteins) == name] <-paste0('log',name)
  colnames(d.cohort_v5_proteins)[colnames(d.cohort_v5_proteins) == 'inv'] <-paste0('inv',name)
}

d.cohort_v3_proteins=d.cohort_v3[c("SampleId",seq_ids)]
names(d.cohort_v3_proteins)=c("SampleId",protein_ids)
for (name in protein_ids) {
  d.cohort_v3_proteins$inv=qnorm((rank(d.cohort_v3_proteins[name],na.last="keep")-0.5)/sum(!is.na(d.cohort_v3_proteins[name])))
  colnames(d.cohort_v3_proteins)[colnames(d.cohort_v3_proteins) == name] <-paste0('log',name)
  colnames(d.cohort_v3_proteins)[colnames(d.cohort_v3_proteins) == 'inv'] <-paste0('inv',name)
}


# 5. Export data as csv
write.csv(d.cohort_v5_proteins,"../datasets/SOMAv5_proteins.csv")
write.csv(d.cohort_v3_proteins,"../datasets/SOMAv3_proteins.csv")
