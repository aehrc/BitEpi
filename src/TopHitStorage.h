#include "def.h"
#define MIN_TOP 1000

// A combination of variables and its score
struct Combination
{
    varIdx ids[MAX_ORDER]; // Ids of variables in a combination
    double score;          // Score of the combination (Alpha or Beta)

    // print: score,var1,var2,...,varN into a file
    void Print(FILE *f, char **varName, uint32 order)
    {
        fprintf(f, "%f", score);
        for (uint32 i = 0; i < order; i++)
            fprintf(f, ",%s", varName[ids[i]]);
        fprintf(f, "\n");
    }
};

// This function is used to sort list of Combinations by their score
// using qsort implemented in stdlib
int CompareCombination(const void *a, const void *b)
{
    Combination *x = (Combination *)a;
    Combination *y = (Combination *)b;
    // Descending order for alpha and Beta
    if (y->score > x->score)
        return 1;
    else
        return -1;
}

// This storage is used to keep top combinations
class TopHitStorage
{
public:
    uint32 numTop;     // Number of top items to be reported
    uint32 numKeep;    // {numTop>MIN_TOP ? numTop : MIN_TOP}
    uint32 numBuf;     // {numKeep*bufRatio} Number of extra item to keep in buffer (to minimise sort operation).
    uint32 numAll;     // {numKeep + numBuf}
    double bufRatio;   // ratio of extra item to the item kept in storage
    Combination *keep; // List of numAll items to keep
    Combination *buf;  // {&keep[numKeep]}
    uint32 bufIdx;     // number of item in the list;
    double minKeep;    // minumum power in numKeep element

    TopHitStorage(uint32 t, double br)
    {
        numTop = t;
        numKeep = (numTop > MIN_TOP) ? numTop : MIN_TOP;
        bufRatio = br;
        numBuf = (uint32)(numKeep * bufRatio);
        numAll = numKeep + numBuf;
        minKeep = 0;
        bufIdx = 0;
        keep = new Combination[numAll];
        NULL_CHECK(keep);
        memset(keep, 0, sizeof(Combination) * numAll);
        buf = &keep[numKeep];
    }

    ~TopHitStorage()
    {
        delete[] keep;
    }

    void Sort()
    {
        qsort(keep, bufIdx, sizeof(Combination), CompareCombination); // to be optimised only the top numKeep are needed
        minKeep = keep[numKeep - 1].score;
        bufIdx = (bufIdx < numKeep) ? bufIdx : numKeep;
    }

    void Add(Combination &c)
    {
        if (c.score > minKeep)
        {
            buf[bufIdx] = c;
            bufIdx++;
            if (bufIdx == numAll)
                Sort();
        }
    }

    void Print(FILE *f, char **varName, uint32 order)
    {
        Sort();
        uint32 numPrint = (numTop > bufIdx) ? bufIdx : numTop;
        for (uint32 i = 0; i < numPrint; i++)
            keep[i].Print(f, varName, order);
    }
};
