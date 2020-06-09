#include "def.h"

class EpiStat
{
public:
    Dataset *dataset;

    ARGS args;
    uint32 threadIdx;
    uint32 jobIdx;

    void *(*threadFunction[4])(void *);

    FILE **topBetaFile;
    FILE **topAlphaFile;

    // below items must be allocated by each thread separately
    word *epiCaseWord[3];
    word *epiCtrlWord[3];

    sampleIdx *contingencyCase;
    sampleIdx *contingencyCtrl;

    STORAGE *topAlpha;
    STORAGE *topBeta;

    void OpenFiles(uint32 order)
    {
        topBetaFile = new FILE *[args.jobsToDo];
        NULL_CHECK(topBetaFile);

        topAlphaFile = new FILE *[args.jobsToDo];
        NULL_CHECK(topAlphaFile);

        char *fn = new char[strlen(args.output) + 20];
        NULL_CHECK(fn);
        for (uint32 j = 0; j < args.jobsToDo; j++)
        {
            char postfix[100];
            sprintf(postfix, "%u.%u", order + 1, args.firstJobIdx + j);

            if (args.printBeta[order])
            {
                sprintf(fn, "%s.Beta.%s.csv", args.output, postfix);
                topBetaFile[j] = fopen(fn, "w");
                NULL_CHECK(topBetaFile)
            }
            if (args.printAlpha[order])
            {
                sprintf(fn, "%s.Alpha.%s.csv", args.output, postfix);
                topAlphaFile[j] = fopen(fn, "w");
                NULL_CHECK(topAlphaFile)
            }
        }
        delete[] fn;
    }

    void CloseFiles(uint32 order)
    {
        for (uint32 t = 0; t < args.jobsToDo; t++)
        {
            if (args.printBeta[order])
            {
                fclose(topBetaFile[t]);
            }
            if (args.printAlpha[order])
            {
                fclose(topAlphaFile[t]);
            }
        }
    }

    EpiStat()
    {
        return;
    }

    ~EpiStat()
    {
    }

    EpiStat(EpiStat *ref)
    {
        memcpy(this, ref, sizeof(EpiStat));
    }

    void Init(Dataset *d, ARGS a, void *(*tf1)(void *), void *(*tf2)(void *), void *(*tf3)(void *), void *(*tf4)(void *))
    {
        dataset = d;
        args = a;
        threadFunction[0] = tf1;
        threadFunction[1] = tf2;
        threadFunction[2] = tf3;
        threadFunction[3] = tf4;
    }

    void AllocateThreadMemory(uint32 o)
    {
        for (uint32 i = 0; i < MAX_ORDER - 1; i++)
        {
            epiCaseWord[i] = new word[dataset->numWordCase];
            epiCtrlWord[i] = new word[dataset->numWordCtrl];

            NULL_CHECK(epiCaseWord[i]);
            NULL_CHECK(epiCtrlWord[i]);
        }

        contingencyCase = new sampleIdx[(uint32)pow(2, MAX_ORDER * 2)];
        contingencyCtrl = new sampleIdx[(uint32)pow(2, MAX_ORDER * 2)];

        NULL_CHECK(contingencyCase);
        NULL_CHECK(contingencyCtrl);

        if (args.topNalpha[o])
        {
            topAlpha = new STORAGE((uint32)args.alpha[o], args.bufRatio);
            NULL_CHECK(topAlpha);
        }
        if (args.topNbeta[o])
        {
            topBeta = new STORAGE((uint32)args.beta[o], args.bufRatio);
            NULL_CHECK(topBeta);
        }
    }

    void FreeThreadMemory(uint32 o)
    {
        for (uint32 i = 0; i < MAX_ORDER - 1; i++)
        {
            delete[] epiCaseWord[i];
            delete[] epiCtrlWord[i];
        }

        delete[] contingencyCase;
        delete[] contingencyCtrl;

        if (args.topNalpha[o])
        {
            delete topAlpha;
        }
        if (args.topNbeta[o])
        {
            delete topBeta;
        }
    }

    void countVariantCombinations(sampleIdx *contingencyTable, uint64_t variantsIn8Samples)
    {
        //each of the 8 bytes contains the variants of a single sample
        for (int byte = 0; byte < 8; byte++)
        {
            //increment count of given byte representing a variant combination
            contingencyTable[(variantsIn8Samples >> (byte * 8)) & 0xFF]++;
        }
    }

    template <uint32_t N, bool CountCombinations = false>
    void OR(varIdx idx)
    {
        const uint32 OIDX = N - 1; // SNPs
        word *caseData = dataset->GetVarCase(OIDX, idx);
        word *ctrlData = dataset->GetVarCtrl(OIDX, idx);

        for (uint32 i = 0; i < dataset->numWordCase; i++)
        {
            word variants = caseData[i];

            if (N > 1)
                variants |= epiCaseWord[OIDX - 1][i];

            if (CountCombinations)
            {
                countVariantCombinations(contingencyCase, variants);
            }
            else
            {
                epiCaseWord[OIDX][i] = variants;
            }
        }

        for (uint32 i = 0; i < dataset->numWordCtrl; i++)
        {
            word variants = ctrlData[i];

            if (N > 1)
                variants |= epiCtrlWord[OIDX - 1][i];

            if (CountCombinations)
            {
                countVariantCombinations(contingencyCtrl, variants);
            }
            else
            {
                epiCtrlWord[OIDX][i] = variants;
            }
        }
    }

    template <uint32_t N>
    void resetContigencyTable()
    {
        size_t arraySize = (1 << (2 * N)) * sizeof(sampleIdx);
        memset(contingencyCtrl, 0, arraySize);
        memset(contingencyCase, 0, arraySize);
    }

    template <uint32_t N>
    double Gini()
    {
        const uint32 entry = integerPow<N>(3);
        double beta = 0;

        for (uint32 i = 0; i < entry; i++)
        {
            uint32 index = cti[i];
            double nCase = (double)contingencyCase[index];
            double nCtrl = (double)contingencyCtrl[index];
            double sum = nCase + nCtrl;
            if (sum)
                beta += (P2(nCase) + P2(nCtrl)) / sum;
        }
        return beta / dataset->numSamples;
    }

    void Epi_1(ThreadData *td)
    {
        const uint32 OIDX = 0; // SNP
        threadIdx = td->threadId;
        jobIdx = td->jobId;

        AllocateThreadMemory(OIDX);

        COMBIN c;
        varIdx *idx = c.idx;

        uint64 cnt = 0;
        for (idx[0] = args.jobs[jobIdx].s[0]; idx[0] <= args.jobs[jobIdx].e[0]; idx[0]++)
        {
            cnt++;
            resetContigencyTable<1>();
            OR<1, true>(idx[0]);
            // compute beta
            double b = Gini<1>();
            c.power = b;
            // report SNP combination if beta meet threshold
            if (args.printBeta[OIDX])
                if (args.topNbeta[OIDX])
                {
                    topBeta->Add(c);
                }
                else if (c.power >= args.beta[OIDX])
                    c.Print(topBetaFile[threadIdx], dataset->nameVariable, OIDX + 1);

            // Save Beta to compute Alpha of next order
            if (args.saveBeta[OIDX])
                SnpBeta[idx[0]] = b;

            // compute Information Gained
            if (args.computeAlpha[OIDX])
            {
                double max_p = dataset->setBeta;

                double a = b - max_p;
                c.power = a;

                // report SNP combination if Alpha meet threshold
                if (args.printAlpha[OIDX])
                    if (args.topNalpha[OIDX])
                    {
                        topAlpha->Add(c);
                    }
                    else if (c.power >= args.alpha[OIDX])
                        c.Print(topAlphaFile[threadIdx], dataset->nameVariable, OIDX + 1);

                // compute the best
                if (args.best)
                    dataset->results[threadIdx].Max_1(a, b, idx);
            }
        }
        args.jobs[jobIdx].counted = cnt;
        if (args.topNalpha[OIDX])
        {
            topAlpha->Print(topAlphaFile[threadIdx], dataset->nameVariable, OIDX + 1);
        }
        if (args.topNbeta[OIDX])
        {
            topBeta->Print(topBetaFile[threadIdx], dataset->nameVariable, OIDX + 1);
        }
        FreeThreadMemory(OIDX);
    }

    void Epi_2(ThreadData *td)
    {
        const uint32 OIDX = 1; // Pair
        threadIdx = td->threadId;
        jobIdx = td->jobId;

        AllocateThreadMemory(OIDX);

        COMBIN c;
        varIdx *idx = c.idx;

        uint64 cnt = 0;
        for (idx[0] = args.jobs[jobIdx].s[0]; idx[0] <= args.jobs[jobIdx].e[0]; idx[0]++)
        {
            OR<1>(idx[0]);
            for (idx[1] = idx[0] + 1; idx[1] < (dataset->numVariables - (OIDX - 1)); idx[1]++)
            {
                cnt++;
                resetContigencyTable<2>();
                OR<2, true>(idx[1]);
                // compute beta
                double b = Gini<2>();
                c.power = b;
                // report SNP combination if beta meet threshold
                if (args.printBeta[OIDX])
                    if (args.topNbeta[OIDX])
                    {
                        topBeta->Add(c);
                    }
                    else if (c.power >= args.beta[OIDX])
                        c.Print(topBetaFile[threadIdx], dataset->nameVariable, OIDX + 1);

                // Save Beta to compute Alpha of next order
                if (args.saveBeta[OIDX])
                    PairBeta[idx[0]][idx[1]] = b;

                // compute Information Gained
                if (args.computeAlpha[OIDX])
                {
                    double max_p = (SnpBeta[idx[1]] > SnpBeta[idx[0]]) ? SnpBeta[idx[1]] : SnpBeta[idx[0]];

                    double a = b - max_p;
                    c.power = a;

                    // report SNP combination if Alpha meet threshold
                    if (args.printAlpha[OIDX])
                        if (args.topNalpha[OIDX])
                        {
                            topAlpha->Add(c);
                        }
                        else if (c.power >= args.alpha[OIDX])
                            c.Print(topAlphaFile[threadIdx], dataset->nameVariable, OIDX + 1);

                    // compute the best
                    if (args.best)
                        dataset->results[threadIdx].Max_2(a, b, idx);
                }
            }
        }
        args.jobs[jobIdx].counted = cnt;
        if (args.topNalpha[OIDX])
        {
            topAlpha->Print(topAlphaFile[threadIdx], dataset->nameVariable, OIDX + 1);
        }
        if (args.topNbeta[OIDX])
        {
            topBeta->Print(topBetaFile[threadIdx], dataset->nameVariable, OIDX + 1);
        }
        FreeThreadMemory(OIDX);
    }

    void Epi_3(ThreadData *td)
    {
        const uint32 OIDX = 2; // Triplet
        threadIdx = td->threadId;
        jobIdx = td->jobId;

        AllocateThreadMemory(OIDX);

        COMBIN c;
        varIdx *idx = c.idx;

        uint64 cnt = 0;
        for (idx[0] = args.jobs[jobIdx].s[0]; idx[0] <= args.jobs[jobIdx].e[0]; idx[0]++)
        {
            varIdx s, e;
            if (idx[0] == args.jobs[jobIdx].s[0])
                s = args.jobs[jobIdx].s[1];
            else
                s = idx[0] + 1;
            if (idx[0] == args.jobs[jobIdx].e[0])
                e = args.jobs[jobIdx].e[1];
            else
                e = (dataset->numVariables - (OIDX - 1)) - 1;

            OR<1>(idx[0]);
            for (idx[1] = s; idx[1] <= e; idx[1]++)
            {
                OR<2>(idx[1]);
                for (idx[2] = idx[1] + 1; idx[2] < dataset->numVariables; idx[2]++)
                {
                    cnt++;
                    resetContigencyTable<3>();
                    OR<3, true>(idx[2]);
                    // compute beta
                    double b = Gini<3>();
                    c.power = b;
                    // report SNP combination if beta meet threshold
                    if (args.printBeta[OIDX])
                        if (args.topNbeta[OIDX])
                        {
                            topBeta->Add(c);
                        }
                        else if (c.power >= args.beta[OIDX])
                            c.Print(topBetaFile[threadIdx], dataset->nameVariable, OIDX + 1);

                    // Save Beta to compute Alpha of next order
                    if (args.saveBeta[OIDX])
                        tripletBeta[idx[0]][idx[1]][idx[2]] = b;

                    // compute Information Gained
                    if (args.computeAlpha[OIDX])
                    {
                        double max_p = (PairBeta[idx[0]][idx[1]] > PairBeta[idx[0]][idx[2]]) ? PairBeta[idx[0]][idx[1]] : PairBeta[idx[0]][idx[2]];
                        max_p = (PairBeta[idx[1]][idx[2]] > max_p) ? PairBeta[idx[1]][idx[2]] : max_p;

                        double a = b - max_p;
                        c.power = a;

                        // report SNP combination if Alpha meet threshold
                        if (args.printAlpha[OIDX])
                            if (args.topNalpha[OIDX])
                            {
                                topAlpha->Add(c);
                            }
                            else if (c.power >= args.alpha[OIDX])
                                c.Print(topAlphaFile[threadIdx], dataset->nameVariable, OIDX + 1);

                        // compute the best
                        if (args.best)
                            dataset->results[threadIdx].Max_3(a, b, idx);
                    }
                }
            }
        }
        args.jobs[jobIdx].counted = cnt;
        if (args.topNalpha[OIDX])
        {
            topAlpha->Print(topAlphaFile[threadIdx], dataset->nameVariable, OIDX + 1);
        }
        if (args.topNbeta[OIDX])
        {
            topBeta->Print(topBetaFile[threadIdx], dataset->nameVariable, OIDX + 1);
        }
        FreeThreadMemory(OIDX);
    }

    void Epi_4(ThreadData *td)
    {
        const uint32 OIDX = 3; // Quadlet
        threadIdx = td->threadId;
        jobIdx = td->jobId;

        AllocateThreadMemory(OIDX);

        COMBIN c;
        varIdx *idx = c.idx;

        uint64 cnt = 0;
        for (idx[0] = args.jobs[jobIdx].s[0]; idx[0] <= args.jobs[jobIdx].e[0]; idx[0]++)
        {
            varIdx s, e;
            if (idx[0] == args.jobs[jobIdx].s[0])
                s = args.jobs[jobIdx].s[1];
            else
                s = idx[0] + 1;
            if (idx[0] == args.jobs[jobIdx].e[0])
                e = args.jobs[jobIdx].e[1];
            else
                e = (dataset->numVariables - (OIDX - 1)) - 1;

            OR<1>(idx[0]);
            for (idx[1] = s; idx[1] <= e; idx[1]++)
            {
                OR<2>(idx[1]);
                for (idx[2] = idx[1] + 1; idx[2] < (dataset->numVariables - (OIDX - 2)); idx[2]++)
                {
                    OR<3>(idx[2]);
                    for (idx[3] = idx[2] + 1; idx[3] < dataset->numVariables; idx[3]++)
                    {
                        cnt++;
                        resetContigencyTable<4>();
                        OR<4, true>(idx[3]);
                        // compute beta
                        double b = Gini<4>();
                        c.power = b;

                        // report SNP combination if beta meet threshold
                        if (args.printBeta[OIDX])
                            if (args.topNbeta[OIDX])
                            {
                                topBeta->Add(c);
                            }
                            else if (c.power >= args.beta[OIDX])
                                c.Print(topBetaFile[threadIdx], dataset->nameVariable, OIDX + 1);

                        // compute Information Gained
                        if (args.computeAlpha[OIDX])
                        {
                            double max_p = (tripletBeta[idx[0]][idx[1]][idx[2]] > tripletBeta[idx[0]][idx[1]][idx[3]]) ? tripletBeta[idx[0]][idx[1]][idx[2]] : tripletBeta[idx[0]][idx[1]][idx[3]];
                            max_p = (tripletBeta[idx[0]][idx[2]][idx[3]] > max_p) ? tripletBeta[idx[0]][idx[2]][idx[3]] : max_p;
                            max_p = (tripletBeta[idx[1]][idx[2]][idx[3]] > max_p) ? tripletBeta[idx[1]][idx[2]][idx[3]] : max_p;

                            double a = b - max_p;
                            c.power = a;

                            // report SNP combination if Alpha meet threshold
                            if (args.printAlpha[OIDX])
                                if (args.topNalpha[OIDX])
                                {
                                    topAlpha->Add(c);
                                }
                                else if (c.power >= args.alpha[OIDX])
                                    c.Print(topAlphaFile[threadIdx], dataset->nameVariable, OIDX + 1);

                            // compute the best
                            if (args.best)
                                dataset->results[threadIdx].Max_4(a, b, idx);
                        }
                    }
                }
            }
        }
        args.jobs[jobIdx].counted = cnt;
        if (args.topNalpha[OIDX])
        {
            topAlpha->Print(topAlphaFile[threadIdx], dataset->nameVariable, OIDX + 1);
        }
        if (args.topNbeta[OIDX])
        {
            topBeta->Print(topBetaFile[threadIdx], dataset->nameVariable, OIDX + 1);
        }
        FreeThreadMemory(OIDX);
    }

    void MultiThread(uint32 o, uint32 numThread, uint32 firstJobIdx)
    {
        ThreadData *td = new ThreadData[numThread];
        NULL_CHECK(td);

        pthread_t *threads = new pthread_t[numThread];
        NULL_CHECK(threads);

        for (uint32 i = 0; i < numThread; i++)
        {
            td[i].epiStat = (void *)new EpiStat(this);
            td[i].threadId = i;
            td[i].jobId = i + firstJobIdx;
        }
        for (uint32 i = 0; i < numThread; i++)
        {
            pthread_create(&threads[i], NULL, threadFunction[o], &td[i]);
        }
        for (uint32 i = 0; i < numThread; i++)
        {
            pthread_join(threads[i], NULL);
            if (args.jobs[td[i].jobId].counted != args.jobs[td[i].jobId].comb)
            {
                printf("\n >>> counted: %llu expected %15.0f", args.jobs[td[i].jobId].counted, args.jobs[td[i].jobId].comb);
                ERROR("Problem in parallelisation please report on GitHub issue page");
            }
        }

        delete[] threads;
        delete[] td;
    }

    void Run()
    {
        AllocateBeta(dataset->numVariables, args);

        for (int o = 0; o < MAX_ORDER; o++)
        {
            if (args.computeBeta[o])
            {

                printf("\n >>>>>> %u-SNP exhaustive search", o + 1);
                printf("\n");
                printf("\n Processing %u jobs [%u..%u] in parallel", args.jobsToDo, args.firstJobIdx, args.lastJobIdx);

                auto begin = std::chrono::high_resolution_clock::now();

                OpenFiles(o);
                args.WorkloadDivider(o + 1, dataset->numVariables, args.numJobs);
                MultiThread(o, args.jobsToDo, args.firstJobIdx);
                CloseFiles(o);

                auto end = std::chrono::high_resolution_clock::now();
                double time_spent = difftime(end, begin);
                if (time_spent == 0)
                    time_spent = 1;

                double c = 0;
                for (uint32 i = 0; i < args.jobsToDo; i++)
                {
                    c += args.jobs[i + args.firstJobIdx].comb;
                }
                printf("\n All jobs are compeleted in %10.3f seconds (%10.0f tests per second)", time_spent, c / time_spent);
                printf("\n");
            }
        }

        printf("\n Merge output of multiple threads (stored in separate files). In Linux it uses command line operation (also echo commands in stdout). In Windows it only merge the best output file.");

#ifndef _MSC_VER
        {
            char *cmd = new char[1024];
            NULL_CHECK(cmd);
            for (uint32 order = 0; order < MAX_ORDER; order++)
            {
                char header[1024];
                switch (order)
                {
                case 0:
                    sprintf(header, "SNP_A");
                    break;
                case 1:
                    sprintf(header, "SNP_A,SNP_B");
                    break;
                case 2:
                    sprintf(header, "SNP_A,SNP_B,SNP_C");
                    break;
                case 3:
                    sprintf(header, "SNP_A,SNP_B,SNP_C,SNP_D");
                    break;
                }

                char sortCmd[100];
                if (args.sort)
                    sprintf(sortCmd, "%s", "| sort -g -r -k1,1 -t ','");
                else
                    sprintf(sortCmd, " ");

                char postfix[100];
                sprintf(postfix, "%u.%u.%u", order + 1, args.firstJobIdx, args.lastJobIdx);

                if (args.printBeta[order])
                {
                    printf("\n Merge Beta%u files...", order + 1);
                    // create a merged output file
                    sprintf(cmd, "cat %s.Beta.%u.*.csv %s | awk 'BEGIN{print(\"Beta,%s\")}{print}' > %s.Beta.%u.csv", args.output, order + 1, sortCmd, header, args.output, order + 1);

                    printf("\n >>> %s", cmd);
                    if (system(cmd) == -1)
                        ERROR("Cannot merge output files");

                    sprintf(cmd, "rm %s.Beta.%u.*.csv", args.output, order + 1);

                    printf("\n >>> %s", cmd);
                    if (system(cmd) == -1)
                        ERROR("Cannot delete temp files");
                }
                if (args.printAlpha[order])
                {
                    printf("\n Merge Alpha%u files...", order + 1);
                    // create a merged output file

                    sprintf(cmd, "cat %s.Alpha.%u.*.csv %s | awk 'BEGIN{print(\"Alpha,%s\")}{print}' > %s.Alpha.%u.csv", args.output, order + 1, sortCmd, header, args.output, order + 1);

                    printf("\n >>> %s", cmd);
                    if (system(cmd) == -1)
                        ERROR("Cannot merge output files");

                    sprintf(cmd, "rm %s.Alpha.%u.*.csv", args.output, order + 1);

                    printf("\n >>> %s", cmd);
                    if (system(cmd) == -1)
                        ERROR("Cannot delete temp files");
                }
            }
            delete[] cmd;
        }
#endif
        if (args.best)
        {
            for (uint32 i = 1; i < args.numThreads; i++)
                dataset->results[0].Max(dataset->results[i]);

            char *fn = new char[strlen(args.output) + 20];
            NULL_CHECK(fn);

            sprintf(fn, "%s.best.csv", args.output);

            dataset->results->toCSV(fn, dataset->nameVariable);
            delete[] fn;
        }

        FreeBeta(dataset->numVariables, &args);
        dataset->FreeMemory(&args);
    }

    void RunCluster() // to be implemented
    {
        AllocateBeta(dataset->numVariables, args);

        if (args.clusterAlpha)
        {
            printf("\n >>>>> Compute and store lower order beta on each node independently (number of jobs is set to number of threads).");

            uint32 o = args.clusterOrder - 1;

            printf("\n >>>>>> %u-SNP exhaustive search", o + 1);
            printf("\n");
            printf("\n Processing %u jobs [%u..%u] in parallel", args.jobsToDo, args.firstJobIdx, args.lastJobIdx);

            time_t begin = time(NULL);

            args.WorkloadDivider(o + 1, dataset->numVariables, args.numThreads);
            MultiThread(o, args.numThreads, 0);

            time_t end = time(NULL);
            double time_spent = difftime(end, begin);
            if (time_spent == 0)
                time_spent = 1;

            double c = 0;
            for (uint32 i = 0; i < args.numThreads; i++)
            {
                c += args.jobs[i].comb;
            }
            printf("\n All jobs are compeleted in %10.0f seconds (%10.0f tests per second)", time_spent, c / time_spent);
            printf("\n");
        }

        if (args.clusterAlpha)
            printf("\n >>>>> Compute alpha on the cluster.");
        else
            printf("\n >>>>> Compute beta on the cluster.");

        uint32 o = args.clusterOrder;

        printf("\n >>>>>> %u-SNP exhaustive search", o + 1);
        printf("\n");
        printf("\n Processing %u jobs [%u..%u] in parallel", args.jobsToDo, args.firstJobIdx, args.lastJobIdx);

        time_t begin = time(NULL);

        OpenFiles(o);
        args.WorkloadDivider(o + 1, dataset->numVariables, args.numJobs);
        MultiThread(o, args.jobsToDo, args.firstJobIdx);
        CloseFiles(o);

        time_t end = time(NULL);
        double time_spent = difftime(end, begin);
        if (time_spent == 0)
            time_spent = 1;

        double c = 0;
        for (uint32 i = 0; i < args.jobsToDo; i++)
        {
            c += args.jobs[i + args.firstJobIdx].comb;
        }
        printf("\n All jobs are compeleted in %10.0f seconds (%10.0f tests per second)", time_spent, c / time_spent);
        printf("\n");

        FreeBeta(dataset->numVariables, &args);
        dataset->FreeMemory(&args);
    }
};
