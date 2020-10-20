

/*
insheet using "/dcl02/leased/kidney/ARIC/static/Genetics/CKDGen/20171017_MW_eGFR_overall_EA_nstud42.dbgap.sig01.txt", clear delim(" ")

keep  rsid effect stderr pvalue allele1  allele2
rename ( rsid effect stderr pvalue) ( rsid effect_eGFR stderr_eGFR pvalue_eGFR)
rename (allele1  allele2) (allele1_eGFR  allele2_eGFR)
save "/dcl02/leased/kidney/ARIC/static/Genetics/CKDGen/20171017_MW_eGFR_overall_EA_nstud42.dbgap.sig01.dta", replace

insheet using "/dcl02/leased/kidney/ARIC/static/Genetics/CKDGen/CKD_overall_EA_JW_20180223_nstud23.dbgap.sig01.txt", clear delim(" ")
keep  rsid effect stderr pvalue allele1  allele2
rename (allele1  allele2) (allele1_CKD  allele2_CKD)
rename (rsid effect stderr pvalue) ( rsid effect_CKD stderr_CKD pvalue_CKD)
save "/dcl02/leased/kidney/ARIC/static/Genetics/CKDGen/CKD_overall_EA_JW_20180223_nstud23.dbgap.sig01.dta", replace




insheet using "/dcl02/leased/kidney/ARIC/static/Genetics/CKDGen/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt", clear delim(" ")

keep  rsid effect stderr pvalue allele1  allele2 freq1
rename ( rsid effect stderr pvalue) ( rsid effect_eGFR stderr_eGFR pvalue_eGFR)
rename (allele1  allele2) (allele1_eGFR  allele2_eGFR)

save "/dcl02/leased/kidney/ARIC/static/Genetics/CKDGen/20171017_MW_eGFR_overall_EA_nstud42.dbgap.dta", replace

insheet using "/dcl02/leased/kidney/ARIC/static/Genetics/CKDGen/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt", clear delim(" ")

keep  rsid effect stderr pvalue allele1  allele2 freq1
rename (allele1  allele2) (allele1_CKD  allele2_CKD)
rename (rsid effect stderr pvalue) ( rsid effect_CKD stderr_CKD pvalue_CKD)

save "/dcl02/leased/kidney/ARIC/static/Genetics/CKDGen/CKD_overall_EA_JW_20180223_nstud23.dbgap.dta", replace
*/



insheet using " /dcl02/leased/kidney/ARIC/projects/01_GWAS_Proteins/results/`1'`2'_Visit`3'_W/`1'`2'_Visit`3'_hits.txt", clear
rename pos pos
merge 1:1  rsid using "/dcl02/leased/kidney/ARIC/static/Genetics/CKDGen/20171017_MW_eGFR_overall_EA_nstud42.dbgap.dta", nogen keep(master match)
merge 1:1  rsid using "/dcl02/leased/kidney/ARIC/static/Genetics/CKDGen/CKD_overall_EA_JW_20180223_nstud23.dbgap.dta", nogen keep(master match)
outsheet using  "../results/`1'`2'_Visit`3'_W/`1'`2'_Visit`3'_hits_CKDGen.txt" , replace

