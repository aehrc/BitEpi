# BitEpi: An exhaustive search of higher-order epistatic interactions

BitEpi is software for an exhaustive search of higher-order epistatic interactions between genomic variants (SNPs) given a binary phenotype (Case/Control).

It can search for 2-SNP (pairwise) interactions as well as 3-SNP and 4-SNP (**Higher-Order**) interactions. It also processes SNPs individually (1-SNP test)

BitEpi performs an **exhaustive** search that means it tests all possible combinations of SNPs. For example, the 4-SNP search of 50 SNPs (sampleData/data.csv) requires to test all "4 out of 50" 4-SNP combinations (230300 4-SNP combinations).

BitEpi implements efficient **multi-threading** (parallelization) such that each thread tests almost the same number of combinations.

**BitEpi could perform 2 different  association analysis: Beta and Alpha**

**Beta** analysis measures the association power of a SNP (1-SNP), a pair of SNPs (2-SNP), a combination of three SNPs (3-SNP) and a combination of four SNPs (4-SNP). Beta represents the combined association power.  Beta is computed as the weighted average purity (Gini-Index) of each row of the contingency table where the weight is the fraction of samples in that row. The below example shows how we compute Beta for an individual SNPs (1-SNP Beta).

| Genotype | # Cases | # Controls |                Purity              |    Weight   | Purity * Weight |
|:--------:|:-------:|:----------:|:----------------------------------:|:-----------:|:---------------:|
|    0/0   |    5    |     10     | ((5\*5)+(10\*10))/(15\*15) = 0.555 | 15/30 = 0.5 |      0.228      |
|    0/1   |    2    |      4     | ((2\*2)+(4\*4))/(6\*6)     = 0.555 | 6/30  = 0.2 |      0.111      |
|    1/1   |    8    |      1     | ((8\*8)+(1\*1))/(9\*9)     = 0.802 | 9/30  = 0.3 |      0.241      |
|          |         |            |                                    | Beta (sum): |    **0.580**    |

2-SNP, 3-SNP, and 4-SNP Beta are computed in the same way. the only difference is that there are more rows in the contingency table (9, 27 and 81 rows respectively).

**Alpha** analysis measure the interactions effect size between SNPs. In other words, Alpha is the gain in the association power due to interaction between SNPs. Alpha computed as below where B_ indicates Beta value. 

- (4-SNP) Alpha WXYZ = B_WXYZ - Max(B_WXY, B_WXZ, B_WYZ, B_XYZ)
- (3-SNP) Alpha XYZ = B_XYZ - Max(B_XY, B_XZ, B_YZ)
- (2-SNP) Alpha XY = B_XY - Max(B_X, B_Y)
- (1-SNP) Alpha X = B_X - B_0

B_0 is the purity (Gini-Index) of the samples in the dataset (CONSTANT). Given r case and q controls, B_0 is computed as (r^2 + q^2) / (r+q)^2. 1-SNP Alpha is similar to 1-SNP Beta with an offset (B_0). While the minimum value for Beta depends on the dataset characteristics (B_0), the minimum value for Alpha is always 0.

**High value for Beta does not necessarily indicate strong interactions**.
For example, let's say

- Beta A  = 0.90 (association power of SNP A)
- Beta B  = 0.03 (association power of SNP B)
- Beta AB = 0.92 (association power of pair of A and B)

In this example, the Beta AB  is high. However, it is driven mainly by Beta A but not the interaction between A and B. Another example is

- Beta X  = 0.51
- Beta Y  = 0.53 
- Beta XY = 0.75

In this example, Beta XY is not as high as Beta AB. However, it is not driven by individual SNPs but by the strong interaction between X and Y.

The Alpha analysis reveals this:

- Alpha AB = Beta AB - Max( Beta A, Beta B) = 0.92 - Max(0.90, 0.03) = **0.02** (Low Alpha -> no interaction)
- Alpha XY = Beta XY - Max( Beta X, Beta Y) = 0.75 - Max(0.51, 0.53) = **0.22** (High Alpha -> strong interaction)

**best Mode**
If you run BitEpi in the best mode, BitEpi finds the best 2-SNP, 3-SNP and 4-SNP interactions for each SNP (maximize Alpha) and then report Beta and Alpha for each SNPs and its best 2-SNP, 3-SNP and 4-SNP interactions. 

**Notes**
- To computer (N)-SNP Alpha the program needs (N-1)-SNP Beta. To avoid computational redundancy, BitEpi computes all (N-1)-SNP Beta and store them (n-1) dimension array. For 4-SNP search, the size of this array is (v^3)*8 bytes where v is the number of SNPs in the dataset (i.e for 4000 SNPs 256 GB of Memory). This is a bottleneck of the program. Using 3 dimension array to store 3-SNP combination results in memory redundancy. However, these data are accessed so frequently and the array leads to the fastest access. 
- You can combine different analyses, for example, you can do 3-SNP Beta and 2-SNP alpha in the same run.
- You can choose to report all SNP combinations that exceed specific Alpha and Beta threshold. However, there are no guidelines on how to choose such a threshold. It would be easier if you ask for the top N SNP combinations (with the highest Alpha or Beta) to be reported. This option is implemented in BitEpi with minor computational cost (see command line parameter Notes)
- If you parallelize the program on so many threads such that each thread process a very small number of combinations to test (less than MIN_COMB_IN_JOB in the code). The program exit without performing the job. In this case, you should use a smaller number of threads.

**BitEpi is Fast**
BitEpi uses fast bitwise operations to count samples with a different genotype. That is why we call it BitEpi.
BOOST (pairwise search only)  and MPI3SNP (3-SNP search only) also use bitwise operations for this purpose. However, the bitwise method used in BitEpi is different and better suits 3-SNP and 4-SNP (higher-order) analysis. As a consequence, BitEpi is slower than BOOST for 2-SNP search but faster than anything else for 3-SNP and 4-SNP search (to the best of our knowledge). 

# BitEpi has two interfaces: Command-Line and Python

**Python Interface**
BitEpi Python Interface is implemented in a separate GitHub page. You can find examples and descriptions there.
It may take time for changes on this repository to be reflected on the other below repository. 
https://github.com/bhosking/bitepi-python

**Command-Line Interface (CLI) Arguments**
| Argument       | Description                               |
|----------------|-------------------------------------------|
| -i [str]       | Path to input CSV file                    |
| -o [str]       | Output prefix                             |
| -sort          | If present sorts the output               |
| -best          | find the best interactions for each SNP   |
| -t [int]       | Number of threads                         |
| -b1 [thr]      | Perform 1-SNP Beta  analysis              |
| -b2 [thr]      | Perform 2-SNP Beta  analysis              |
| -b3 [thr]      | Perform 3-SNP Beta  analysis              |
| -b4 [thr]      | Perform 4-SNP Beta  analysis              |
| -a1 [thr]      | Perform 1-SNP Alpha analysis              |
| -a2 [thr]      | Perform 2-SNP Alpha analysis              |
| -a3 [thr]      | Perform 3-SNP Alpha analysis              |
| -a4 [thr]      | Perform 4-SNP Alpha analysis              |

**notes**
- If thr<1 then thr is the minimum threshold on alpha and beta to be reported in the output.
- If thr>=1 then thr is the number of top hits in alpha and beta to be reported in the output.
- If thr>=1 and more than 1 thread are used, then each thread reports thr top hits and the output files are merged. So (t*thr) top hits will be reported. You can use -sort option and only consider the top thr record.
- thr is optional. If you don't pass a thr the program computes the metric but it does not report anything (performance testing).
- If all interactions should be reported set thr to 0.

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

# Cite BitEpi
The paper is not published yet (will be available soon). You may cite our GitHub page for now.
