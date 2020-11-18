BitEpi.o -i DS -bfile -sort -a1 50 -a2 50 -a3 50 -a4 50 -o DS

python3 ../../Pvalue.py rnd DS 10 1000 p-val-1r.tsv 1 &
python3 ../../Pvalue.py rnd DS 10 1000 p-val-2r.tsv 2 &
wait
python3 ../../Pvalue.py rnd DS 10 1000 p-val-3r.tsv 3 &
python3 ../../Pvalue.py rnd DS 10 1000 p-val-4r.tsv 4 &
wait
exit
python3 ../../Pvalue.py epi DS 10 1000 p-val-1.tsv DS.Alpha.1.csv &
python3 ../../Pvalue.py epi DS 10 1000 p-val-2.tsv DS.Alpha.2.csv &
wait
python3 ../../Pvalue.py epi DS 10 1000 p-val-3.tsv DS.Alpha.3.csv &
python3 ../../Pvalue.py epi DS 10 1000 p-val-4.tsv DS.Alpha.4.csv &
wait
