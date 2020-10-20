
# Run the result of this from the local side
cd ../results
for f in * 
do 
echo sftp -r  asurapan@jhpce01.jhsph.edu:/dcl02/leased/kidney/ARIC/projects/01_GWAS_Proteins/results/$f \""O:/Users/asurapa2/04\ Protein\ GWAS/02\ GWAS\ Results"\"
done


#sftp
# put /drives/o/Users/asurapa2/04\ Protein\ GWAS/01\ Extracting\ proteins\ and\ ARIC\ covariates/SOMAv5_proteins.csv /users/asurapan/04_ARIC.GWAS.with.FAST
