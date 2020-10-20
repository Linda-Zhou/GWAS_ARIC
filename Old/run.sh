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

#read -p "$1 transformed protien $2 at visit $3. Press enter if correct"

# 1. Preparing the phenotype file
module load R
Rscript 02.preparePheno_protein_raw.r $1 $2 $3


# 2. Make the script to launch the GWAS scripts 
cat 04.launch_whites_topmed.sh | sed "s/transformation/$1/" | sed "s/protein/$2/" | sed "s/visit/$3/" > /users/asurapan/04_ARIC.GWAS.with.FAST/${1}${2}_Visit${3}_W/launch_whites_topmed.sh
bash /users/asurapan/04_ARIC.GWAS.with.FAST/${1}${2}_Visit${3}_W/launch_whites_topmed.sh
sleep 2h  


# 3. Process the GWAS results
cd /users/asurapan/04_ARIC.GWAS.with.FAST/${1}${2}_Visit${3}_W
cp /users/asurapan/04_ARIC.GWAS.with.FAST/FAST_Scripts/05.format.FAST.results.pl /users/asurapan/04_ARIC.GWAS.with.FAST/${1}${2}_Visit${3}_W/format.FAST.results.pl
perl /users/asurapan/04_ARIC.GWAS.with.FAST/${1}${2}_Visit${3}_W/format.FAST.results.pl
sleep 2h

# 4. Make QQPlots and Manhattan plots
Rscript 06.TOPMed.GWAS.qc_adi.R $1 $2 $3

# 5. Clean up log files and temp scripts
rm /users/asurapan/04_ARIC.GWAS.with.FAST/${1}${2}_Visit${3}_W/format.FAST.results.pl
rm /users/asurapan/04_ARIC.GWAS.with.FAST/${1}${2}_Visit${3}_W/launch_whites_topmed.sh
rm /users/asurapan/04_ARIC.GWAS.with.FAST/FAST_Scripts/Log/*


# qmem check how much memory did a job use??
#  try running without the sleeps
# do if [-f filename] 
# that'll check if a file exists
# fi

