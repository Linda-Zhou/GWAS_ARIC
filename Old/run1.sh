#!/bin/bash
#$ -N FAST_GWAS
#$ -cwd
#$ -l mem_free=80G,h_vmem=80G
#$ -e ../log
#$ -o ../log
##$ -m e
##$ -M asurapa2@jhu.edu

module load perl
module unload conda_R

read -p "$1 transformed protien $2 at visit $3. Press enter if correct"


# 4. Make QQPlots and Manhattan plots
module load R
Rscript 06.TOPMed.GWAS.qc_adi.R $1 $2 $3

# 5. Clean up log files and temp scripts
rm ../results/${1}${2}_Visit${3}_W/format.FAST.results.pl
rm ../results/${1}${2}_Visit${3}_W/launch_whites_topmed.sh
#rm ../log/*

# 6.
stata -b  07.merge_CKDGen.do $1 $2 $3


# qmem check how much memory did a job use??
#  try running without the sleeps
# do if [-f filename] 
# that'll check if a file exists
# fi

