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
 https://pubmed.ncbi.nlm.nih.gov/33838163/
 https://jasn.asnjournals.org/content/32/9/2291
 
