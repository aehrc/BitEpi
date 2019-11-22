#!/bin/bash

#set -x

# check if there is only one argument 
if [ "$#" -ne 1 ]; then
    echo "The program accepts only one argument that is the GAMETES file name without .txt"
	exit
fi

# check if the input file exist
FILE=$1.txt
if test -f "$FILE"; then
    echo "$FILE exists"
else
	echo "Error: $FILE does not exist"
	exit;
fi

# check if datamesh installed
pkgs='datamash'
if ! dpkg -s $pkgs >/dev/null 2>&1; then
	echo "Error: Install $pkgs first using"
	echo ">>> sudo apt-get install $pkgs"
	exit
else
	echo "$pkgs is installed"
fi

#check if plink is installed
pkgs='plink'
if hash $pkgs 2>/dev/null; then
	echo "$pkgs is installed"
else
	echo "Error: Install $pkgs 1.9 first from"
	echo ">>> https://www.cog-genomics.org/plink/1.9/"
	exit
fi

SNP=$(head -n 1 $1.txt | tr \\t \\n | tail -n +2 | wc -l)
SAM=$(tail -n +2 $1.txt | wc -l)

echo ">>>>> There are $SNP SNPs and $SAM samples in the file"

tail -n +2 $1.txt | cut -f 1-$SNP | tr 0 x | tr 1 y | tr 2 z | sed 's/x/0\/0/g' | sed 's/y/0\/1/g' | sed 's/z/1\/1/g' | awk '{print("S"NR"\t"$0)}' | datamash transpose | awk 'BEGIN{print("##fileformat=VCFv4.0");print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");}{if(NR==1) print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"$0); else print("1\t"NR-1"\t.\tA\tC\t.\t.\t.\tGT\t"$0)}' > $1.vcf  
cut -f $(( $SNP+1 )) $1.txt | tail -n +2 | awk '{print("S"NR" S"NR" "$1+1)}' > $1.pheno

cut -d ' ' -f 3 $1.pheno | sed 's/1/0/g' | sed 's/2/1/g' | tr \\n \\t | sed 's/.$//g' | awk '{print("VarSam\t"$0)}' > $1.l.tsv

cut -f 1-$SNP $1.txt | datamash transpose | cat $1.l.tsv - | tr \\t ','  > $1.epi.csv

plink --vcf $1.vcf --pheno $1.pheno --recode transpose -out $1
cat $1.tped | tr ' ' \\t | cut -f 2-1000000000 > $1.x
cut -f 1 -d ',' $1.epi.csv | tail -n +2 | paste - $1.x | awk '{$2=$1;$1="1";print($0)}' | tr \\t ' '> $1.tped

rm $1.log $1.l.tsv $1.nosex $1.pheno $1.vcf $1.x
