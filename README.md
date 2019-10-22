# BitEpi

BitEpi is a software to identify epistasis interactions between genomic variants (SNPs). The method is described in the paper (link will be available once published). BitEpi uses fast bitwise operations to speedup epistasis analysis. That is why we call it BitEpi. It uses an entropy-based metric to evaluate association power and the interaction effect size. BitEpi performs an exhaustive search for 2-SNP, 3-SNP and 4-SNP interactions.

**BitEpi has two interfaces: Command-Line and Python**

# Python Interface:
BitEpi Python Interface is implemented in a separate GitHub page. You can find examples and descriptions there.
https://github.com/bhosking/bitepi-python

# Command-Line Interface (CLI) Arguments:
| Argument       | Description                                                                                     |
|----------------|-------------------------------------------------------------------------------------------------|
| -i [str]       | Path to input CSV file                                                                          |
| -o [str]       | Output prefix                                                                                   |
| -sort          | If present sorts the output of alpha and beta test by  alpha and beta value in descending order |
| -best          | Sort output files by Beta and Information-Gained                                                |
| -t [int]       | Number of threads                                                                               |
| -b1 [thr]      | Perform 1-SNP beta test with the given threshold [thr]                                          |
| -b2 [thr]      | Perform 2-SNP beta test with the given threshold [thr]                                          |
| -b3 [thr]      | Perform 3-SNP beta test with the given threshold [thr]                                          |
| -b4 [thr]      | Perform 4-SNP beta test with the given threshold [thr]                                          |
| -a1 [thr]      | Perform 1-SNP alpha test with the given threshold [thr]                                         |
| -a2 [thr]      | Perform 2-SNP alpha test with the given threshold [thr]                                         |
| -a3 [thr]      | Perform 3-SNP alpha test with the given threshold [thr]                                         |
| -a4 [thr]      | Perform 4-SNP alpha test with the given threshold [thr]                                         |

# Command-Line Interface (CLI) Arguments for Cluster/Cloud mode:
| Argument       | Description                                                                                     |
|----------------|-------------------------------------------------------------------------------------------------|
| -c             | Cluster/Cloud Mode                                                                              |
| -j [int]       | Total number of jobs (best value: total number of threads in the cluster)                       |
| -f [int]       | First job index to be processed on this computer (starting from 0)                              |
| -t [int]       | Number of jobs in parallel on this computer (best value: number of threads on this computer)    |

**notes**
- [thr] is the threshold between 0 and 1 (on alpha and beta). All SNPs or interactive SNPs that exceed the threshold [thr] will be listed in the output file. 
- [thr] is optional and if you don't pass a [thr] the program computes the metric but it does not report anything (does not create output file). For performance testing only.
- If you want all interactions set thr to 0.
- If you set thr to 1 the program creates an empty output file.
- The best mode does not work yet with the cluster mode.
- if t>j-f then the program runs the remaining jobs in parallel. For example, if there are 100 jobs and the first job index is 95, only 5 jobs left to be done [95,96,97,98 and 99]. In this case, if you set the -t to be larger than 5 (i.e 20), the program only runs the remaining 5 jobs in parallel (use 5 parallel threads)

**Input Format**
- The first row includes labels: 1 and 0 for cases and controls respectively
- The first column includes SNP unique ids. BitEpi does not check the uniqueness.
- The first entry (first col and first row) is ignored.
- All other entry represents genotype and can be 0, 1 or 2 for 0/0, 0/1 and 1/1 genotype (phased genotype are not supported). 

**Output File Naming**
- output_prefix.best.csv
- output_prefix.Alpha.Order.csv
- output_prefix.Beta.Order.csv
- output_prefix.Alpha.Order.FirstJobIndex.LastJobIndex.csv
- output_prefix.Beta.Order.FirstJobIndex.LastJobIndex.csv

# To Compile:
Use the g++ commands in compile.sh or just run it
```sh
$ cd BitEpi
$ bash compile.sh
```

#To Run Test:
See two examples in runme.sh or just run it
```sh
$ cd BitEpi
$ bash runme.sh
```
