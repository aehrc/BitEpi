#!/bin/bash

plink --bfile $1 --recode vcf bgz --out $1
awk '{if($6==1)print($1"_"$2"\tFalse"); else print($1"_"$2"\tTrue");}' $1.fam | sort -k1b,1 > $1.an
bcftools query -l $1.vcf.gz | awk '{print($1"\t"NR)}' | sort -k1b,1 > $1.so
join $1.an $1.so | tr ' ' \\t | sort -k3g,3 | cut -f 1-2 > $1.anso
echo "Sample-isCase" | tr - \\t | cat - $1.anso > $1.label
bash VCF_2_EPI.sh $1.vcf.gz $1.label  2

rm $1.an
rm $1.anso
rm $1.label
rm $1.so
rm $1.vcf.gz

echo "the output is $1.vcf.gz.epi.csv"