# GWAS_ARIC
Genome-wide association study (or GWAS), also known as whole genome association study, is an observational study of a genome-wide set of genetic variants in 
different individuals to see if any variant is associated with a trait. GWAS studies typically focus on associations between single-nucleotide polymorphisms (SNPs)
and traits like major human diseases, but can equally be applied to any other genetic variants and any other organisms. This project aims at investigating the association 
between SNPs and kidney function traits in ARIC (Atherosclerosis Risk in Communities) study. This GWAS project uses FAST.2.4 package. 

## Steps to implement GWAS
 - Prepare Input Genotype Data in FAST format, individual id file (option --indiv-file), Phenotype + Covariate File (option --trait-file). Note that For both single SNP and gene-based Cox methods, the phenotype 
 file requires 7 mandatory columns: Family ID, Individual ID, Paternal ID, Maternal ID, Sex (1=male; 2=female; other=unknown), Status, Time to Event. Details refer to: 
 https://bitbucket.org/baderlab/fast/wiki/InputFileFormats (script 00 - 03)
 - Run FAST in linux enviornment (script 04 - 05)
 - perform QC (quality control) to filter MAF > 0.05 and produce Manhattan plot. （script 06 - 07）
 - merge the GWAS summary statistics with CKDgen kidney database (script 08)
 
 ## Results 
 Grams ME, Surapaneni A, Chen J, et al. Proteins associated with risk of kidney function decline in the general population. J Am Soc Nephrol. 2021;32(9):2291. http://jasn.asnjournals.org/content/32/9/2291.abstract. doi: 10.1681/ASN.2020111607.

### ABSTRACT
#### Background 
Proteomic profiling may allow identification of plasma proteins that associate with subsequent
changesin kidney function, elucidating biologic processes underlying the development and progression of CKD.
#### Methods 
We quantified the association between 4877 plasma proteins and a composite outcome of ESKD or
decline in eGFR by 50% among 9406 participants in the Atherosclerosis Risk in Communities (ARIC) Study
(visit 3; mean age, 60 years) who were followed for a median of 14.4 years. We performed separate analyses
for these proteins in a subset of 4378 participants(visit 5), whowere followed at a later time point, for amedian
of 4.4 years. For validation, we evaluated proteins with significant associations (false discovery rate <5%) in
both time periods in 3249 participants in the Chronic Renal Insufficiency Cohort (CRIC) and 703 participants
in the African American Study of Kidney Disease and Hypertension (AASK). We also compared the genetic
determinants of protein levels with those from a meta-analysis genome-wide association study of eGFR.
#### Results 
In models adjusted for multiple covariates, including baseline eGFR and albuminuria, we identified 13
distinct proteins that were significantly associated with the composite end point in both time periods, 
including TNF receptor superfamily members 1A and 1B, trefoil factor 3, and b-trace protein. Of these proteins, 12
were also significantly associated in CRIC, and nine were significantly associated in AASK. Higher levels of
each protein associated with higher risk of 50% eGFR decline or ESKD.We found genetic evidence for a causal
role for one protein, lectin mannose-binding 2 protein (LMAN2).
#### Conclusions 
Large-scale proteomic analysis identified both known and novel proteomic risk factors for eGFR
decline.
