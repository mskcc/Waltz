#! /bin/bash

# aggregate genotypes from sample mafs in the current folder into one genotypes.maf file



cat *-genotypes.maf | grep Hugo_Symbol | head -1 > genotypes.maf

sampleField=`cat genotypes.maf | awk 'BEGIN{FS=OFS="\t"}{for(i=1; i<=NF; i++){if($i=="Tumor_Sample_Barcode"){print i; break}}}'`

# combine sample mafs
for f in `ls *-genotypes.maf`
do
	sample=`echo $f | awk '{if(index($1, "_")!=0) split($1, a, "_"); else split($1, a, "-genotypes.maf"); print a[1]}'`
  
  	# add genotypes to genotypes.txt
	awk -v sampleField=$sampleField -v sample=$sample 'BEGIN{FS=OFS="\t";}{if(NR==1) next; $sampleField=sample; print}' $f >> genotypes.maf

done




#
