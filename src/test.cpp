#include "def.h"
// #include "Arguments.h"
// #include "TopHitStorage.h"
#include "WorkLoad.h"

void TestLoopTime(uint32 n)
{
    uint32 num = 0;
    for (uint32 i = 0; i < n - 3; i++)
        for (uint32 j = i + 1; j < n - 2; j++)
            for (uint32 k = j + 1; k < n - 1; k++)
                for (uint32 l = k + 1; l < n; l++)
                    num++;
    printf("\n>>>> %u <<<<\n", num);
}

#ifdef TEST
int main(int argc, char *argv[])
{
    printf("\n=============Start Test=============");
    printf("\n====================================\n");

    // Arguments args;
    // args.PrintHelp();
    // args.Parse(argc, argv);
    // args.Print();

    WorkLoad wl;
    wl.Test();
    printf("\n=======Multi Threaded==========");
    RunAllJobs(wl);

    // TestLoopTime(5000); // 24 second
    printf("\n =========== End Test\n");
    return 0;
}
#endif