#!/bin/bash
#$ -N FAST
#$ -cwd
#$ -l mem_free=600M,h_vmem=600M
#$ -e ../log
#$ -o ../log
##$ -m e
##$ -M asurapa2@jhu.edu

module load perl
module unload conda_R


dir_data=('/dcs01/arking/ARIC_static/ARIC_Data/GWAS/TOPMed/EA/by.chunks/dosage')
dir_analysis=.
dir_output=.
CHR=$(printf $SGE_TASK_ID)
covar=0
# age, center1,2,3, PC1,-10


time /users/asurapan/Software/FAST.2.4.bitbucket/FAST \
--mode genotype \
--tped-file $dir_data/$part.tped.gz \
--snpinfo-file $dir_data/$part.snpinfo \
--mlinfo-file $dir_data/$part.mlinfo \
--indiv-file /dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/W_ARIC_topmed.imputation.iid \
--trait-file logNBL1_SeqId_2944_66.visit3.W.FAST.runfile.iid.txt \
--out-file $dir_output/$part \
--chr $CHR \
--imputation-quality 0.1 \
--verbose \
--num-covariates $covar \
--linear-snp


