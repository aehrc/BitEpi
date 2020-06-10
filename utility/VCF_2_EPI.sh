#!/bin/bash

VcfFile=$1          # Pass to the VCF file
SampleAnnot=$2      # Pass to the sample annotation.csv file
isCase=$3           # the column number in sample annotation.csv file that has the isCase (True/False)

echo "=============================================================="

echo " Argument 1: path to vcf or vcf.bgz file"
echo " Argument 2: path tab seprated (tsv) sample annotation file"
echo " Argument 3: isCase column number in sample annotation file"

echo "=============================================================="

echo " Argument 1: $VcfFile"
echo " Argument 2: $SampleAnnot"
echo " Argument 3: $isCase"

echo "=============================================================="

echo " ***** This scripts uses bcftools."
echo " ***** This scripts uses vcf file rsid column as SNP id. It should be uniq id and should not be empty."
echo " ***** This scripts does not check if rsid exist and if it is uniqu for each snp"
echo " ***** The order of samples in VCF file and sample annotation file should be the same."
echo " ***** The first line of the sample annotation file should be header line."
echo " ***** The sample annotation file should tab separated (tsv) file."
echo " ***** isCase column should only include 'True' and 'False'"
echo " ***** The output is vcffile.epi.csv"
echo " ***** This scripts does not check input arguments."


echo " >>>>> Test if sampels are in the same order in vcf and sample annotation file"
bcftools query -l $VcfFile > $VcfFile.s.tsv
cut -f 1 $SampleAnnot | tail -n +2 > $SampleAnnot.s.tsv
if cmp -s  $VcfFile.s.tsv $SampleAnnot.s.tsv; then
    echo "Sample list are matched"
else
   echo "Sample list are not matched"
   exit
fi

echo " >>>>> Converting...."

cut -f $isCase $SampleAnnot | tail -n +2 | sed 's/False/0/g' | sed 's/True/1/g' | tr \\n \\t | awk '{print("VarSam\t"$0)}' > $SampleAnnot.l.tsv
bcftools query -f '%ID\n' $VcfFile > $VcfFile.id.tsv
bcftools query -f '[%GT\t]\n' $VcfFile | tr -d '|' | tr -d '/' | tr '\.' '0' | sed 's/00/0/g' | sed 's/01/1/g' | sed 's/10/1/g' | sed 's/11/2/g' | paste $VcfFile.id.tsv - | cat $SampleAnnot.l.tsv - | tr \\t ',' | sed 's/.$//g' > $VcfFile.epi.csv

echo " >>>>> Done"

echo " >>>>> Remove temporary files"
rm $VcfFile.s.tsv $VcfFile.id.tsv $SampleAnnot.s.tsv $SampleAnnot.l.tsv

echo " >>>>> Output: $VcfFile.epi.csv"
