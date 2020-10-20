#!/bin/bash
#$ -N TOPMedQC
#$ -cwd
#$ -l mem_free=6G,h_vmem=6G
#$ -e /users/asurapan/04_ARIC.GWAS.with.FAST/QC_Scripts/Log
#$ -o /users/asurapan/04_ARIC.GWAS.with.FAST/QC_Scripts/Log
##$ -m e
##$ -M asurapa2@jhu.edu


module load R

Rscript TOPMed.GWAS.qc_adi.R $1 $2 $3 

