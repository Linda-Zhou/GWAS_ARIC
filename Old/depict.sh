

trans=$1
seqid=$2
visit=$3


awk '{print $12}' /dcl02/leased/kidney/ARIC/projects/01_GWAS_Proteins/results/${trans}${seqid}_Visit${visit}_W/${trans}${seqid}_Visit${visit}_hits.txt > /dcl02/leased/kidney/software/depict/${trans}${seqid}_Visit${visit}_hits_rsids.txt


cd /dcl02/leased/kidney/software/depict
echo "python depict.py ${trans}${seqid}_Visit${visit}_hits_rsids.txt  ${trans}${seqid}_Visit${visit}_W"
python depict.py ${trans}${seqid}_Visit${visit}_hits_rsids.txt  ${trans}${seqid}_Visit${visit}_W

mv ${trans}${seqid}_Visit${visit}_hits_rsids.txt  /dcl02/leased/kidney/ARIC/projects/01_GWAS_Proteins/results/${trans}${seqid}_Visit${visit}_W
mv results/${trans}${seqid}_Visit${visit}_W* /dcl02/leased/kidney/ARIC/projects/01_GWAS_Proteins/results/${trans}${seqid}_Visit${visit}_W








