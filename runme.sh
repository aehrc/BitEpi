#!bin/bash
set -x

./BitEpi.o  -i sampleData/data.csv -o sampleData/out -t 2 -b1 0 -b2 0 -b3 0 -b4 0 -a1 0 -a2 0 -a3 0 -a4 0 -sort
./BitEpi.o  -i sampleData/data.csv -o sampleData/out -t 2 -best

set +x

echo ">>>>>>> Output files"
ls -l sampleData/out*
echo "<<<<<<<<<<<<<<<<<<<<"

exit
