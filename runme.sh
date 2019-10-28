#!bin/bash
set -x

./BitEpi.o  -i sampleData/data.csv -o sampleData/out -t 2 -b4 0.3 -a4 0.1 -sort
./BitEpi.o  -i sampleData/data.csv -o sampleData/out -t 1 -best

set +x

echo ">>>>>>> Output files"
ls -l sampleData/out*
echo "<<<<<<<<<<<<<<<<<<<<"

exit
