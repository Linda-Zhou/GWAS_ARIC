
capture rm "../results/ARIC_Morgan_indexSNPs.dta"
foreach prot in NBL1_SeqId_2944_66 SVEP1_SeqId_11109_56 SVEP1_SeqId_11178_21 FBLN5_SeqId_15585_304 PLA2G2A_SeqId_2692_74 IGFBP4_SeqId_2950_57 ASGR1_SeqId_5452_71 HAVCR1_SeqId_9021_1 MMP7_SeqId_2789_26 GM2A_SeqId_15441_6 LY86_SeqId_3623_84 HPX_SeqId_15347_12 DOK2_SeqId_19578_19 SAA2_SeqId_18832_65 ZHX3_SeqId_10036_201 CLMP_SeqId_10440_26 PTGDS_SeqId_10514_5 WFDC2_SeqId_11388_75 DSC2_SeqId_13126_52 TNFRSF1A_SeqId_2654_19 TNFRSF1B_SeqId_3152_57 FSTL3_SeqId_3438_10 TFF3_SeqId_4721_54 TFF3_SeqId_8323_163 TNFRSF1B_SeqId_8368_102 DLK2_SeqId_9359_9 LMAN2_SeqId_9468_8 DCTN2_SeqId_5879_51 {
	
        di "/dcl02/leased/kidney/ARIC/projects/01_GWAS_Proteins/results/log`prot'_Visit3_W/log`prot'_Visit3_hits.txt"
	import delimited "/dcl02/leased/kidney/ARIC/projects/01_GWAS_Proteins/results/log`prot'_Visit3_W/log`prot'_Visit3_hits.txt", clear asdouble
	sort pvalue
	//	di "`prot'"
	capture  list in 1
	if _rc>0{
	capture import delimited "/dcl02/leased/kidney/ARIC/projects/01_GWAS_Proteins/results/log`prot'_Visit3_W/allchr.FAST.txt.trim.wheader.txt", clear asdouble
	if _rc>0 {
        capture import delimited "/dcl02/leased/kidney/ARIC/projects/01_GWAS_Proteins/results/log`prot'_Visit3_W/allchr.FAST.txt.trim", clear asdouble
	
	}
	}	
	capture gen rsid=""
	capture gen protein="`prot'"
	capture rename p pvalue
	keep snpid rsid protein beta se caf pvalue
	order protein snpid rsid pvalue
	noisily	list in 1
 	keep in 1	
	capture append using "../results/ARIC_Morgan_indexSNPs"
	save "../results/ARIC_Morgan_indexSNPs", replace
		
}
use "../results/ARIC_Morgan_indexSNPs", clear
sort protein
export delimited using "../results/ARIC_Morgan_indexSNPs.txt", replace
