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
dir_analysis=../results/${1}${2}_Visit${3}_W
dir_output=../results/${1}${2}_Visit${3}_W
#CHR=$(printf $SGE_TASK_ID)
covar=0
# age, center1,2,3, PC1,-10,genetic Pcs 1-10


time /dcl02/leased/kidney/software/FAST.2.4.bitbucket/FAST  \
--mode genotype \
--tped-file $dir_data/$part.tped.gz \
--snpinfo-file $dir_data/$part.snpinfo \
--mlinfo-file $dir_data/$part.mlinfo \
--indiv-file /dcs01/arking/arkinglab/active/projects/aric/analyses/topmed.imputation/W_ARIC_topmed.imputation.iid \
--trait-file $dir_analysis/${1}${2}.visit${3}.W.FAST.runfile.iid.txt \
--out-file $dir_output/$part \
--chr $CHR \
--imputation-quality 0.1 \
--verbose \
--num-covariates $covar \
--linear-snp


