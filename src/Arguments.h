#include "def.h"

class Arguments
{
public:
    char inputFilePrefix[MAX_FILE_NAME_LEN + 100];  // [-input]   Input file prefix including path and without extension (csv or plink bfile)
    bool isPlink;                                   // [-isPlink] Is the input Plink bfile (default: false)
    char outputFilePrefix[MAX_FILE_NAME_LEN + 100]; // [-output]  Output file prefix including path and without extension
    char setFileName[MAX_FILE_NAME_LEN + 100];      // [-set]     Set file name
    char orderFileName[MAX_FILE_NAME_LEN + 100];    // [-order]   Order file name
    uint32 numJob;                                  // [-jobs]
    uint32 firstJob;                                // [-first]
    uint32 lastJob;                                 // [-last]
    uint32 jobsToDo;                                // [(lastJob  - firstJob) + 1] number of jobs to run in parallel on this computer
    bool isBeta;                                    // [-isBeta]  should the program find top hits by alpha or beta (default: false)
    bool keepTop;                                   // [-top]     (default true) if (-top argument) >= 1
    uint32 numTop;                                  // [-top]     (default 0)
    double threshold;                               // [-top]     (default 1000)
    double bufRatio;                                // [-br]      (default 2)

    Arguments()
    {
        memset(this, 0, sizeof(Arguments));

        isPlink = false;
        numJob = 1;
        firstJob = 1;
        lastJob = numJob;
        jobsToDo = (lastJob - firstJob) + 1;
        isBeta = false;
        keepTop = true;
        numTop = 1000;
        bufRatio = 2;
    }

    void PrintHelp()
    {
        printf("\n========================================================");
        printf("\n -input   [path] Input file prefix without extension (csv or plink bfile)");
        printf("\n -isPlink        If presents, input is plink bfile otherwise CSV");
        printf("\n -output  [path] Output file prefix");
        printf("\n -set     [path] (Optional) Set file name (see Manual)");
        printf("\n -order   [path] (Optional) Order file name (see Manual)");
        printf("\n -jobs    [int]  Number of jobs (default 1)");
        printf("\n -first   [int]  First job index on this computer (default 1)");
        printf("\n -last    [path] Last job index on this computer (default -t value)");
        printf("\n -isBeta         If presents, the program perform Beta analysis otherwise Alpha analysis");
        printf("\n -top     [num]  see Manual (default 1000) ");
        printf("\n -br      [num]  Buffer Ratio (see Manual) (default 2)");
        printf("\n==========================================================\n");
        return;
    }

    void Print()
    {
        printf("\n=========================================");
        printf("\n Given arguments:");
        printf("\n Input file prefix                         %s", inputFilePrefix);
        if (isPlink)
            printf("\n Input file is Plink bfile");
        else
            printf("\n Input file is in CSV format");
        printf("\n Output file prefix                        %s", outputFilePrefix);
        if (setFileName[0])
            printf("\n Set file name                             %s", setFileName);
        if (orderFileName[0])
            printf("\n Order file name                           %s", orderFileName);
        printf("\n Number of jobs                            %u", numJob);
        printf("\n First job                                 %u", firstJob);
        printf("\n Last job                                  %u", lastJob);
        printf("\n Number of parallel jobs on this computer  %u", jobsToDo);
        if (keepTop)
            printf("\n Number of top hits to report              %u", numTop);
        else
            printf("\n Threshold to report top hits              %f", threshold);
        if (isBeta)
            printf("\n Perform Beta analysis");
        else
            printf("\n Perform Alpha analysis");
        printf("\n Buffer ratio                              %f", bufRatio);
        printf("\n=========================================\n");
    }

    void Parse(int argc, char *argv[])
    {
        if (argc == 1)
        {
            PrintHelp();
            ERROR("Please Enter Arguments.");
        }

        int i = 1;
        while (true)
        {
            if (!strcmp(argv[i], "-input"))
            {
                i++;
                if (i == argc)
                    ERROR("Please enter input file prefix");
                if (strlen(argv[i]) > MAX_FILE_NAME_LEN)
                    ERROR("Input file prefix is too long");
                strcpy(inputFilePrefix, argv[i]);
            }
            else if (!strcmp(argv[i], "-output"))
            {
                i++;
                if (i == argc)
                    ERROR("Please enter output file prefix");
                if (strlen(argv[i]) > MAX_FILE_NAME_LEN)
                    ERROR("Output file prefix is too long");
                strcpy(outputFilePrefix, argv[i]);
            }
            else if (!strcmp(argv[i], "-set"))
            {
                i++;
                if (i == argc)
                    ERROR("Please enter set file name");
                if (strlen(argv[i]) > MAX_FILE_NAME_LEN)
                    ERROR("Set file name is too long");
                strcpy(setFileName, argv[i]);
            }
            else if (!strcmp(argv[i], "-order"))
            {
                i++;
                if (i == argc)
                    ERROR("Please enter order file name");
                if (strlen(argv[i]) > MAX_FILE_NAME_LEN)
                    ERROR("Order file name is too long");
                strcpy(orderFileName, argv[i]);
            }
            else if (!strcmp(argv[i], "-isPlink"))
            {
                isPlink = true;
            }
            else if (!strcmp(argv[i], "-isBeta"))
            {
                isBeta = true;
            }
            else if (!strcmp(argv[i], "-jobs"))
            {
                i++;
                if (i == argc)
                    ERROR("Please enter number of jobs");
                numJob = atoi(argv[i]);
                if (!numJob)
                    ERROR("Please enter valid number of jobs");
            }
            else if (!strcmp(argv[i], "-first"))
            {
                i++;
                if (i == argc)
                    ERROR("Please enter first job index");
                firstJob = atoi(argv[i]);
                if (!firstJob)
                    ERROR("Please enter valid first job index");
            }
            else if (!strcmp(argv[i], "-last"))
            {
                i++;
                if (i == argc)
                    ERROR("Please enter last job index");
                lastJob = atoi(argv[i]);
                if (!lastJob)
                    ERROR("Please enter valid last job index");
            }
            else if (!strcmp(argv[i], "-br"))
            {
                i++;
                if (i == argc)
                    ERROR("Please enter buffer ratio");
                bufRatio = atof(argv[i]);
                if (!bufRatio)
                    ERROR("Please enter valid buffer ratio");
            }
            else if (!strcmp(argv[i], "-top"))
            {
                i++;
                if (i == argc)
                    ERROR("Please enter top value");
                threshold = atof(argv[i]);
                if (!threshold)
                    ERROR("Please enter valid threshold");
                if (threshold < 1)
                {
                    keepTop = false;
                    numTop = 0;
                }
                else
                {
                    keepTop = true;
                    numTop = (uint32)threshold;
                    threshold = 0;
                }
            }
            else
            {
                printf("%s is not a valid arguments", argv[i]);
                ERROR("Please enter valid arguments");
            }
            i++;
            if (i == argc)
                break;
        }
        jobsToDo = (lastJob - firstJob) + 1;
    };
};