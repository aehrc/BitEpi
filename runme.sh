#!/bin/bash
set -x

#Examples on how to run BitEpi
./BitEpi.o  -i sampleData/data.csv -o sampleData/out1 -t 1 -b1 10 -a1 0.001 -b2 0.502 -a2 35 -a3 0.003
./BitEpi.o  -i sampleData/data.csv -o sampleData/out2 -t 2 -b3 0.505 -b4 40 -a4 30 -sort
./BitEpi.o  -i sampleData/data.csv -o sampleData/out3 -t 1 -best

#An example how to convert GAMETES simulated data to the format accepted by BitEpi and transposed plink file (for MPI3SNP and BOOST)
bash GAMETES_2_EPI.sh sampleData/GAMETES_Example_Output

#An example how to convert VCF file to BitEpi input file.
bash VCF_2_EPI.sh sampleData/Hipster.vcf.bgz sampleData/Hipster.tsv 2

#Examles on how to run cluster mode (not implemented yet)
#./BitEpi.o  -i sampleData/data.csv -o sampleData/out4 -c -j 6 -f 3 -t 2 -a4 100

#Examples to run Pvalue module
BitEpi -i sampleData/bfile -bfile -o sampleData/bfile -sort -a2 100
python3 Pvalue.py rnd sampleData/bfile 5 100 sampleData/pval-rnd.tsv 2
python3 Pvalue.py epi sampleData/bfile 5 100 sampleData/pval-epi.tsv sampleData/bfile.Alpha.2.csv

set +x

echo ">>>>>>> Output files"
ls -l sampleData/out*
echo "<<<<<<<<<<<<<<<<<<<<"

exit
