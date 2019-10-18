# BitEpi

BitEpi is a software to identify epistasis interactions between genomic variants (SNPs). The method is described in the paper (link will be available once published). BitEpi uses fast bitwise operations to speedup epistasis analysis. That is why we call it BitEpi. It uses an entropy-based metric to evaluate association power and the interaction effect size. BitEpi performs an exhaustive search for 2-SNP, 3-SNP and 4-SNP interactions.

**BitEpi has two interfaces: Command-Line and Python**

# Python Interface:
BitEpi Python Interface is implemented in a separate GitHub page. You can find examples and descriptions there.
https://github.com/bhosking/bitepi-python

# Command-Line Interface (CLI) Arguments:
| Argument       | Description                                                                                     |
|----------------|-------------------------------------------------------------------------------------------------|
| -t             | Input CSV file                                                                                  |
| -o             | Output prefix                                                                                   |
| -t             | Number of threads (parallel Processing)                                                         |
| -sort          | If present sorts the output of alpha and beta test by  alpha and beta value in descending order |
| -best          | Sort output files by Beta and Information-Gained                                                |
| -b1 [thr]      | Perform 1-SNP beta test with the given threshold [thr]                                          |
| -b2 [thr]      | Perform 2-SNP beta test with the given threshold [thr]                                          |
| -b3 [thr]      | Perform 3-SNP beta test with the given threshold [thr]                                          |
| -b4 [thr]      | Perform 4-SNP beta test with the given threshold [thr]                                          |
| -a1 [thr]      | Perform 1-SNP alpha test with the given threshold [thr]                                         |
| -a2 [thr]      | Perform 2-SNP alpha test with the given threshold [thr]                                         |
| -a3 [thr]      | Perform 3-SNP alpha test with the given threshold [thr]                                         |
| -a4 [thr]      | Perform 4-SNP alpha test with the given threshold [thr]                                         |

**notes**
- [thr] is the threshold between 0 and 1 (on alpha and beta). All SNPs or interactive SNPs that exceed the threshold [thr] will be listed in the output file. 
- [thr] is optional and if you don't pass a [thr] the program computes the metric but it does not report anything (does not create output file). For performance testing only.
- If you want all interactions set thr to 0.
- If you set thr to 1 the program creates an empty output file.  

**Input Format**
- The first row includes labels: 1 and 0 for cases and controls respectively
- The first column includes SNP unique ids. BitEpi does not check the uniqueness.
- The first entry (first col and first row) is ignored.
- All other entry represents genotype and can be 0, 1 or 2 for 0/0, 0/1 and 1/1 genotype (phased genotype are not supported). 

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
