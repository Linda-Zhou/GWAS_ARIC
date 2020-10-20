#Script is modified from Dan's format.FAST.results.pl script.

#!/usr/local/bin/perl -w



for ($i=1;$i<23;$i++) {
	system "rm chr$i.FAST.out; for f in \$(ls chr$i.*.Linear.txt | sort -t \".\" -nk 2,3); do grep -v SNP.id \$f | cut -f1-7,10-12,15 >>chr$i.FAST.out; done";
}

system 'echo -e \'SNPID\\tChr\\tPosition\\tnon_coded_allele\\tcoded_allele\\tBeta\\tSE\\tCAF\\tQual\\tESampleSize\\tP\' >allchr.FAST.txt';
for ($i=1;$i<23;$i++) {
	system "cat chr$i.FAST.out >>allchr.FAST.txt; rm chr$i.FAST.out;";
}

#Note, there is no X chr data for TOPMed Imputation

#Move all immediate FAST output file into rawdata dir
system "mkdir -p rawdata";
system "rm chr* ";
system "rm FAST.o* ";
system "rm FAST.e* ";

#REMINDER TO REMOVE PROBLEMATIC SNPS FROM AA ANALYSES#
system 'echo -e \'For blacks, run remove.snps.topmed.R\\n\'';


exit;

