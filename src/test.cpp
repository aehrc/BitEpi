#include "def.h"
#include "Arguments.h"
#include "TopHitStorage.h"

#ifdef TEST
int main(int argc, char *argv[])
{
    printf("\n=============Start Test=============");
    printf("\n====================================\n");

    Arguments args;
    args.PrintHelp();
    args.Parse(argc, argv);
    args.Print();
}
#endif