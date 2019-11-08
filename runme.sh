#!/bin/bash
set -x

#Examples on how to run BitEpi
./BitEpi.o  -i sampleData/data.csv -o sampleData/out -t 1 -b1 10 -a1 0.001 -b2 0.502 -a2 35 -a3 0.003
./BitEpi.o  -i sampleData/data.csv -o sampleData/out -t 2 -b3 0.505 -b4 40 -a4 30 -sort
./BitEpi.o  -i sampleData/data.csv -o sampleData/out -t 1 -best

#An example how to convert GAMETES simulated data to the format accepted by BitEpi and transposed plink file (for MPI3SNP and BOOST)
bash GAMETES_2_EPI.sh sampleData/GAMETES_Example_Output

set +x

echo ">>>>>>> Output files"
ls -l sampleData/out*
echo "<<<<<<<<<<<<<<<<<<<<"

exit
