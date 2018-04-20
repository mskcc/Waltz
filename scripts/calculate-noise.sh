#!/bin/bash

# calculate substituion rate ("noise") for all pileup and pileup-without-duplicates files in the current folder. Run it in a folder with Waltz output


echo -e "Sample\tGenotypeCount\tAltCount\tAltPercent\tMethod" > noise.txt

for f in `ls *-pileup.txt`
do
  sampleName=`basename $f`
  sampleName=${sampleName/-IGO*/}

  awk -v sample=$sampleName 'BEGIN{FS=OFS="\t"; base[5]="A"; base[6]="C"; base[7]="G"; base[8]="T"}{max=-1; garbage=0; genotype=-1; for(i=5; i<=8; i++){if($i>max){max=$i; genotype=i}} for(i=5; i<=8; i++){if(i==genotype) continue; if($i>max*.1) next}   for(i=5; i<=8; i++){if(i==genotype) continue; substitution=base[genotype]">"base[i]; genotypeCount[substitution]+=max; altCount[substitution]+=$i;}}END{for(s in genotypeCount){g+=genotypeCount[s]; a+=altCount[s]} print sample, g, a, 100*a/(g+a), "Total" }' $f >> noise.txt

  g=${f/-pileup.txt/-pileup-without-duplicates.txt}

  awk -v sample=$sampleName 'BEGIN{FS=OFS="\t"; base[5]="A"; base[6]="C"; base[7]="G"; base[8]="T"}{max=-1; garbage=0; genotype=-1; for(i=5; i<=8; i++){if($i>max){max=$i; genotype=i}} for(i=5; i<=8; i++){if(i==genotype) continue; if($i>max*.1) next}   for(i=5; i<=8; i++){if(i==genotype) continue; substitution=base[genotype]">"base[i]; genotypeCount[substitution]+=max; altCount[substitution]+=$i;}}END{for(s in genotypeCount){g+=genotypeCount[s]; a+=altCount[s]} print sample, g, a, 100*a/(g+a), "Picard" }' $g >> noise.txt


done







#
