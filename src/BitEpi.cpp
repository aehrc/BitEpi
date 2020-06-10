
#include "def.h"
#include "Arguments.h"

double ***tripletBeta;
double **PairBeta;
double *SnpBeta;

template <class D1, class D2>
double difftime(D1 end, D2 begin)
{
	return std::chrono::duration<double, std::ratio<1, 1>>(end - begin).count();
}

double Combination(uint32 v, uint32 o)
{
	double nc = (double)v;

	if (o == 1)
		return nc;

	if (o == v)
		return 1;

	if (o > v)
		return 0;

	switch (o)
	{
	case 2:
		return (nc * (nc - 1)) / 2;
	case 3:
		return (nc * (nc - 1) * (nc - 2)) / 6;
	case 4:
		return (nc * (nc - 1) * (nc - 2) * (nc - 3)) / 24;
	default:
		ERROR("Does not support combination above 4");
	}
}

struct JOB
{
	uint32 id;
	uint32 s[2];
	uint32 e[2];
	double comb;
	double diff;  // Difference to average.
	double aDiff; // Accumulative difference to average.
	uint64 counted;
	void Print()
	{
		printf("\n %6u (%6u,%6u) (%6u,%6u) %15.0f %15.0f %15.0f", id + 1, s[0] + 1, s[1] + 1, e[0] + 1, e[1] + 1, comb, diff, aDiff);
	}
	void PrintHead()
	{
		printf("\n Outer  loop iterates from S1 to E1");
		printf("\n Second loop iterates from S2 to E2");
		printf("\n Combinations   : SNP combinations to be tested in each job");
		printf("\n diffToAvg      : Difference to average number of combination to be tested in each job");
		printf("\n AccumDiffToAvg : Accumulative diffToAvg");
		printf("\n Job ID (    S1,    S2) (    E1,    E2)    Combinations       diffToAvg  AccumDiffToAvg");
	}
	void Print1SNP()
	{
		printf("\nJob:%6u process %15.0f SNPs (%10u ... %10u)", id + 1, comb, s[0] + 1, e[0] + 1);
	}
};

struct ThreadData
{
	void *epiStat;	 // epi class
	uint32 threadId; // thread id
	uint32 jobId;	 // job Id to be executed on that threads
};

void *EpiThread_1(void *t)
{
	ThreadData *td = (ThreadData *)t;
	EpiStat *epiStat = (EpiStat *)td->epiStat;

	printf("\n Thread %5u processing Job %5u ...", td->threadId + 1, td->jobId + 1);
	printf("\n");

	auto begin = std::chrono::high_resolution_clock::now();

	epiStat->Epi_1(td);

	auto end = std::chrono::high_resolution_clock::now();
	double time_spent = difftime(end, begin);
	if (time_spent == 0)
		time_spent = 1;
	printf("\n Thread %5u processed Job %5u in %10.3f seconds (%10.0f tests per second)", td->threadId + 1, td->jobId + 1, time_spent, epiStat->args.jobs[td->jobId].comb / time_spent);
	printf("\n");
	return NULL;
}

void *EpiThread_2(void *t)
{
	ThreadData *td = (ThreadData *)t;
	EpiStat *epiStat = (EpiStat *)td->epiStat;

	printf("\n Thread %5u processing Job %5u ...", td->threadId + 1, td->jobId + 1);
	printf("\n");

	auto begin = std::chrono::high_resolution_clock::now();

	epiStat->Epi_2(td);

	auto end = std::chrono::high_resolution_clock::now();
	double time_spent = difftime(end, begin);
	if (time_spent == 0)
		time_spent = 1;
	printf("\n Thread %5u processed Job %5u in %10.3f seconds (%10.0f tests per second)", td->threadId + 1, td->jobId + 1, time_spent, epiStat->args.jobs[td->jobId].comb / time_spent);
	printf("\n");
	return NULL;
}

void *EpiThread_3(void *t)
{
	ThreadData *td = (ThreadData *)t;
	EpiStat *epiStat = (EpiStat *)td->epiStat;

	printf("\n Thread %5u processing Job %5u ...", td->threadId + 1, td->jobId + 1);
	printf("\n");

	auto begin = std::chrono::high_resolution_clock::now();

	epiStat->Epi_3(td);

	auto end = std::chrono::high_resolution_clock::now();
	double time_spent = difftime(end, begin);
	if (time_spent == 0)
		time_spent = 1;
	printf("\n Thread %5u processed Job %5u in %10.3f seconds (%10.0f tests per second)", td->threadId + 1, td->jobId + 1, time_spent, epiStat->args.jobs[td->jobId].comb / time_spent);
	printf("\n");
	return NULL;
}

void *EpiThread_4(void *t)
{
	ThreadData *td = (ThreadData *)t;
	EpiStat *epiStat = (EpiStat *)td->epiStat;

	printf("\n Thread %5u processing Job %5u ...", td->threadId + 1, td->jobId + 1);
	printf("\n");

	auto begin = std::chrono::high_resolution_clock::now();

	epiStat->Epi_4(td);

	auto end = std::chrono::high_resolution_clock::now();
	double time_spent = difftime(end, begin);
	if (time_spent == 0)
		time_spent = 1;
	printf("\n Thread %5u processed Job %5u in %10.3f seconds (%10.0f tests per second)", td->threadId + 1, td->jobId + 1, time_spent, epiStat->args.jobs[td->jobId].comb / time_spent);
	printf("\n");
	return NULL;
}

#ifndef TEST
int main1(int argc, char *argv[])
{
	printf("\n=============Start=============");
	printf("\n===============================\n");

	Arguments args;
	args.PrintHelp();
	//args.Parse(argc, argv);
	//args.Print();

	// dataset.ReadDatasetBfile(args.input);

	// dataset.Init(args);

	// EpiStat epiStat;
	// epiStat.Init(&dataset, args, EpiThread_1, EpiThread_2, EpiThread_3, EpiThread_4);

	// epiStat.Run();

	printf("\n=============Finish=============\n\n\n");
	return 0;
}
#endif