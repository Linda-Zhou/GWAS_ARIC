#$ -N FAST_GWAS
#$ -cwd
#$ -l mem_free=60G,h_vmem=60G
#$ -e ../log
#$ -o ../log
##$ -m e

#module load perl
module unload conda_R

#read -p "$1 transformed protien $2 at visit $3. Press enter if correct"
#
## 1. Preparing the phenotype file
module load R
Rscript 02.preparePheno_protein_raw.r $1 $2 $3
Rscript 02.preparePheno_eigenprotein_raw.r $1 $2 $3

#
#
## 2. Make the script to launch the GWAS scripts 
cat 04.launch_whites_topmed.sh | sed "s/transformation/$1/" | sed "s/protein/$2/" | sed "s/visit/$3/" > ../results/${1}${2}_Visit${3}_W/launch_whites_topmed.sh
bash ../results/${1}${2}_Visit${3}_W/launch_whites_topmed.sh
sleep 120m  
#
#
## 3. Process the GWAS results
cp 05.format.FAST.results.pl ../results/${1}${2}_Visit${3}_W/format.FAST.results.pl
pushd ../results/${1}${2}_Visit${3}_W
perl format.FAST.results.pl
#
#
## 4. Make QQPlots and Manhattan plots
popd
Rscript 06.TOPMed.GWAS.qc_adi.R $1 $2 $3
#
## 5. Clean up log files and temp scripts
rm ../results/${1}${2}_Visit${3}_W/format.FAST.results.pl
rm ../results/${1}${2}_Visit${3}_W/launch_whites_topmed.sh
rm ../log/FAST.o*
rm ../log/FAST.e*
rm Rplots.pdf
rmdir ../results/${1}${2}_Visit${3}_W/rawdata
#
## 6.
stata -b  07_merge_CKDGen.do $1 $2 $3
mv 07_merge_CKDGen.log ../results/${1}${2}_Visit${3}_W/07_merge_CKDGen.log
#
##
##
### qmem check how much memory did a job use??
###  try running without the sleeps
### do if [-f filename] 
### that'll check if a file exists
### fi
#
