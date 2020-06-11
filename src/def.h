#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pthread.h"
#include "unistd.h"
#include "math.h"
#include <vector>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>

#include "csvparser.h"

#define TEST

#define MAX_FILE_NAME_LEN 1000

typedef unsigned char uint8;
typedef unsigned short int uint16;
typedef unsigned int uint32;
typedef unsigned long long int uint64;
typedef int int32;

typedef unsigned int varIdx;
typedef unsigned short int sampleIdx; // This type used in contingency table. This table should be kept in the cache so choose the smallest possible type here. Note that short int is 16 bit and can deal with up to 2^16 (~65,000) samples.
typedef unsigned long long int word;  // for parallel processing
#define byte_in_word 8

#define MAX_ORDER 4

#define P2(X) (X * X)
#define P3(X) (X * X * X)
#define P4(X) (X * X * X * X)

#define ERROR(X)                                                                 \
    {                                                                            \
        printf("\n *** ERROR: %s (line:%u - File %s)\n", X, __LINE__, __FILE__); \
        exit(0);                                                                 \
    }

#define NULL_CHECK(X)                                                                         \
    {                                                                                         \
        if (!X)                                                                               \
        {                                                                                     \
            printf("\n *** ERROR: %s is null (line:%u - File %s)\n", #X, __LINE__, __FILE__); \
            exit(0);                                                                          \
        }                                                                                     \
    }

template <uint32 EXP>
uint32 integerPow(uint32 x);

template <>
uint32 integerPow<1>(uint32 x);

template <>
uint32 integerPow<0>(uint32);

// uint64 NchoosK(uint32 N, uint32 K);