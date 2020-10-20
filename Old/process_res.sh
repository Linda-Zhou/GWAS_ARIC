#!/bin/bash
#$ -N FAST_GWAS
#$ -cwd
#$ -l mem_free=80G,h_vmem=80G
#$ -e /users/asurapan/04_ARIC.GWAS.with.FAST/FAST_Scripts/Log
#$ -o /users/asurapan/04_ARIC.GWAS.with.FAST/FAST_Scripts/Log
##$ -m e
##$ -M asurapa2@jhu.edu

module load perl
module unload conda_R
module load R

read -p "$1 transformed protien $2 at visit $3. Press enter if correct"

# 4. Make QQPlots and Manhattan plots
Rscript 06.TOPMed.GWAS.qc_adi.R $1 $2 $3

# 5. Clean up log files and temp scripts
rm /users/asurapan/04_ARIC.GWAS.with.FAST/${1}${2}_Visit${3}_W/format.FAST.results.pl
rm /users/asurapan/04_ARIC.GWAS.with.FAST/${1}${2}_Visit${3}_W/launch_whites_topmed.sh
rm /users/asurapan/04_ARIC.GWAS.with.FAST/FAST_Scripts/Log/*

# 6.
stata -b  07_merge_CKDGen.do $1 $2 $3


# qmem check how much memory did a job use??
#  try running without the sleeps
# do if [-f filename] 
# that'll check if a file exists
# fi

