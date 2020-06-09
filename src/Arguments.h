#include "def.h"

class Arguments
{
public:
    char plinkPrefix[1024];
    char outputPrefix[1024];
    char setFileName[1024];
    uint32 numJob;
    uint32 firstJob;
    uint32 lastJob;
    uint32 jobsToDo;  // {(lastJob - firstJob) + 1} number of jobs to do on this computer (could be less than number of threads)
    bool beta;        // if true only compute Beta but not Alpha
    uint32 order;     // the order of analysis
    bool keepTop;     // if (-x argument) >= 1
    uint32 numTop;    // driven by -x argument
    double threshold; // driven by -x argument

    double bufRatio;
    double numComb;
    double avgJobNumCombinations;

    Arguments()
    {
        memset(this, 0, sizeof(Arguments));

        numJob = 1;
        firstJob = 1;
        lastJob = numJob;
        jobsToDo = (lastJob - firstJob) + 1;
        beta = false;
        order = 2;
        keepTop = true;
        numTop = 1000;
    }

    void PrintHelp()
    {
        printf("\n========================================================\n");
        printf(" -i [path]  Input prefix (PLINK 1.9 bfile)\n");
        printf(" -o [path]  Output prefix\n");
        printf(" -s [path]  Variants on each line (separated by tab) are considered as a sets\n");
        printf("            If not presented a set of all variants are considered\n");
        printf(" -t [int]   Number of parallel threads (jobs) (default 1)\n");
        printf(" -f [int]   First job index (default 1)\n");
        printf(" -l [path]  Last job index (default -t value)\n");
        printf(" -b         If presents perform Beta analysis otherwise perform Alpha analysis\n");
        printf(" -h [int]   Order of analysis\n");
        printf(" -x [int]   (default 1000) \n");
        printf("            If -s is less than 1 then it acts like minimum value of Alpha or Beta\n");
        printf("            (depending on the analysis) to be reported.\n");
        printf("            If -s is grater or equal to 1 then it identify the number of top most\n");
        printf("            significant hits to be reported by each thread\n");
        printf("\n==========================================================\n");
        return;
    }

    void Print() // print the given option for debugging
    {
        printf("\n=========================================");
        printf("\n Given Arguments:");
        printf("\n Input Prefix                %s", plinkPrefix);
        printf("\n Output Prefix               %s", outputPrefix);
        printf("\n Set File Name               %s", setFileName);
        printf("\n Number of Jobs              %u", numJob);
        printf("\n First Job                   %u", firstJob);
        printf("\n Last Job                    %u", lastJob);
        printf("\n Order                       %u", order);
        if (keepTop)
            printf("\n Num Top                     %u", numTop);
        else
            printf("\n Threshold                   %f", threshold);
        if (beta)
            printf("\n Compute Beta");
        else
            printf("\n Compute Alpha");

        printf("=========================================\n");
    }

    // void Parse(int argc, char *argv[])
    // {
    //     double d = -1;
    //     uint32 o = 0;
    //     char str[100];
    //     bool next = false;

    //     // Read arguments

    //     for (int i = 1; i < argc; i++)
    //     {
    //         next = false;

    //         // check if beta flag is passed for any order
    //         for (uint32 o = 0; o < MAX_ORDER; o++)
    //         {
    //             sprintf(str, "-b%u", o + 1);
    //             if (!strcmp(argv[i], str))
    //             {
    //                 // set the computep flag
    //                 computeBeta[o] = true;
    //                 betaGiven[o] = true;
    //                 sprintf(clusterCmd, "%s -b%u", clusterCmd, o + 1);
    //                 // check if there is any threashold argument to this option
    //                 if ((i + 1) != argc)
    //                 {
    //                     if (argv[i + 1][0] != '-')
    //                     {
    //                         d = atof(argv[i + 1]);
    //                         if ((d == 0) && (argv[i + 1][0] != '0'))
    //                         {
    //                             printf("\n Cannot parse -b%u [%s]", o + 1, argv[i + 1]);
    //                             PrintHelp(argv[0]);
    //                         }
    //                         // set the threshold and print flag for beta
    //                         sprintf(clusterCmd, "-b%u %f", o + 1, d);
    //                         beta[o] = d;
    //                         printBeta[o] = true;
    //                         if (beta[o] >= 1)
    //                             topNbeta[o] = true;
    //                         i++;
    //                     }
    //                 }
    //                 next = true;
    //                 break;
    //             }
    //         }
    //         // if this option is processed move to the next option
    //         if (next)
    //             continue;

    //         // check if  flag is passed for any order
    //         for (uint32 o = 0; o < MAX_ORDER; o++)
    //         {
    //             sprintf(str, "-a%u", o + 1);
    //             if (!strcmp(argv[i], str))
    //             {
    //                 // set the computep flag and beta flag of previous order
    //                 computeBeta[o] = computeAlpha[o] = true;
    //                 alphaGiven[o] = true;
    //                 sprintf(clusterCmd, "%s -a%u", clusterCmd, o + 1);
    //                 if (o > 0)
    //                     computeBeta[o - 1] = saveBeta[o - 1] = true;

    //                 // check if there is any threashold argument to this option
    //                 if ((i + 1) != argc)
    //                 {
    //                     if (argv[i + 1][0] != '-')
    //                     {
    //                         d = atof(argv[i + 1]);
    //                         if ((d == 0) && (argv[i + 1][0] != '0'))
    //                         {
    //                             printf("\n Cannot parse -a%u [%s]", o + 1, argv[i + 1]);
    //                             PrintHelp(argv[0]);
    //                         }
    //                         // set the threshold and print flag for a
    //                         sprintf(clusterCmd, "-a%u %f", o + 1, d);
    //                         alpha[o] = d;
    //                         printAlpha[o] = true;
    //                         if (alpha[o] >= 1)
    //                             topNalpha[o] = true;
    //                         i++;
    //                     }
    //                 }
    //                 next = true;
    //                 break;
    //             }
    //         }
    //         if (next)
    //             continue;

    //         // read input file name
    //         if (!strcmp(argv[i], "-i"))
    //         {
    //             if ((i + 1) == argc)
    //             {
    //                 printf("\n Please enter path to the input file");
    //                 PrintHelp(argv[0]);
    //             }

    //             if (argv[i + 1][0] != '-')
    //             {
    //                 strcpy(input, argv[i + 1]);
    //                 inputGiven = true;
    //             }
    //             else
    //             {
    //                 printf("\n Please enter path to the input file");
    //                 PrintHelp(argv[0]);
    //             }
    //             i++;
    //             continue;
    //         }

    //         // read output file prefix
    //         if (!strcmp(argv[i], "-o"))
    //         {
    //             if ((i + 1) == argc)
    //             {
    //                 printf("\n Please enter an output prefix");
    //                 PrintHelp(argv[0]);
    //             }

    //             if (argv[i + 1][0] != '-')
    //                 strcpy(output, argv[i + 1]);
    //             else
    //             {
    //                 printf("\n Please enter an output prefix");
    //                 PrintHelp(argv[0]);
    //             }
    //             i++;
    //             continue;
    //         }

    //         // read number of threads
    //         if (!strcmp(argv[i], "-t"))
    //         {
    //             if ((i + 1) == argc)
    //             {
    //                 printf("\n Please enter the number of threads");
    //                 PrintHelp(argv[0]);
    //             }

    //             if (argv[i + 1][0] != '-')
    //             {
    //                 numThreads = atoi(argv[i + 1]);
    //                 if (numThreads == 0)
    //                 {
    //                     printf("\n Please enter a valid number of threads greater than 0 but not [%s]", argv[i + 1]);
    //                     PrintHelp(argv[0]);
    //                 }
    //             }
    //             else
    //             {
    //                 printf("\n Please enter the number of threads");
    //                 PrintHelp(argv[0]);
    //             }
    //             i++;
    //             continue;
    //         }

    //         // read bfile flag
    //         if (!strcmp(argv[i], "-bfile"))
    //         {
    //             readBfile = true;
    //             continue;
    //         }

    //         // read best flag
    //         if (!strcmp(argv[i], "-best"))
    //         {
    //             best = true;
    //             continue;
    //         }

    //         // read sort flag
    //         if (!strcmp(argv[i], "-sort"))
    //         {
    //             sprintf(clusterCmd, "%s -sort", clusterCmd);
    //             sort = true;
    //             continue;
    //         }

    //         // read bufRatio
    //         if (!strcmp(argv[i], "-r"))
    //         {
    //             if ((i + 1) == argc)
    //             {
    //                 printf("\n Please enter the buffer ratio");
    //                 PrintHelp(argv[0]);
    //             }

    //             if (argv[i + 1][0] != '-')
    //             {
    //                 bufRatio = atof(argv[i + 1]);
    //                 if (bufRatio < 2)
    //                 {
    //                     printf("\n Please enter a valid ratio for buffer greater than 2 but not [%s]", argv[i + 1]);
    //                     PrintHelp(argv[0]);
    //                 }
    //             }
    //             else
    //             {
    //                 printf("\n Please enter the buffer ratio");
    //                 PrintHelp(argv[0]);
    //             }
    //             i++;
    //             continue;
    //         }
    //         printf("\n Invalid option %s", argv[i]);
    //         PrintHelp(argv[0]);
    //     }

    //     // check arguments

    //     // There should be an input file in all cases
    //     // if the output prefix does not pass we use default output prefix that is "OUTPUT_BitEpi"
    //     if (!inputGiven)
    //     {
    //         printf("\n You should specify the path to the input file (-i)");
    //         PrintHelp(argv[0]);
    //     }

    //     if (master && best)
    //     {
    //         printf("\n best mode does not work in master mode yet");
    //         PrintHelp(argv[0]);
    //     }

    //     if (!clusterMode)
    //     {
    //         if (numJobs != -1 || firstJobIdx != -1)
    //         {
    //             printf("\n -j and -f are only used in cluster/cloud mode (-c presented)");
    //             PrintHelp(argv[0]);
    //         }
    //         firstJobIdx = 0;
    //         numJobs = numThreads;
    //     }
    //     else
    //     {
    //         if (master)
    //         {
    //             printf("\n You run the program in master mode. You cannot use -c in master mode");
    //             PrintHelp(argv[0]);
    //         }
    //         if (best)
    //         {
    //             printf("\n best mode does not work in cluster/cloud mode yet");
    //             PrintHelp(argv[0]);
    //         }

    //         uint anCnt = 0;
    //         for (uint32 i = 0; i < 4; i++)
    //         {
    //             if (betaGiven[i])
    //             {
    //                 anCnt++;
    //                 clusterAlpha = false;
    //                 clusterOrder = i;
    //             }
    //             if (alphaGiven[i])
    //             {
    //                 anCnt++;
    //                 clusterAlpha = true;
    //                 clusterOrder = i;
    //             }
    //         }
    //         if (anCnt > 1)
    //         {
    //             printf("\n In cluster mode only one analysis (alpha or beta) can be executed");
    //             PrintHelp(argv[0]);
    //         }
    //         if (clusterOrder == 0)
    //         {
    //             printf("\n cluster mode works only for 2-SNP, 3-SNP and 4-SNP tests but not 1-SNP test");
    //             PrintHelp(argv[0]);
    //         }

    //         if (numJobs == -1 || firstJobIdx == -1)
    //         {
    //             printf("\n Specify number of jobs (-j) and first job index (-f) in cluster/cloud mode (-c presented)");
    //             PrintHelp(argv[0]);
    //         }
    //     }
    //     jobsToDo = ((numJobs - firstJobIdx) < numThreads) ? (numJobs - firstJobIdx) : numThreads;
    //     lastJobIdx = firstJobIdx + jobsToDo - 1;

    //     if (firstJobIdx >= numJobs)
    //     {
    //         printf("\n First job index out of range [%u>=%u]", firstJobIdx, numJobs);
    //         PrintHelp(argv[0]);
    //     }

    //     // apply best
    //     if (best)
    //         for (uint32 o = 0; o < MAX_ORDER; o++)
    //         {
    //             computeBeta[o] = saveBeta[o] = computeAlpha[o] = true;
    //         }

    //     if (computeBeta[0])
    //         maxOrder = 1;
    //     if (computeBeta[1])
    //         maxOrder = 2;
    //     if (computeBeta[2])
    //         maxOrder = 3;
    //     if (computeBeta[3])
    //         maxOrder = 4;

    //     if (!maxOrder)
    //     {
    //         printf("\n No test to be performed.");
    //         PrintHelp(argv[0]);
    //     }
    // }
};
