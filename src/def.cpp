#include "def.h"

// we use 2 bits (4 states) to represent a genotype. However a genotype has only 3 states.
// for 4-SNP, blow table is used to translate (4 states)^(4 SNPs) state to (3 states)^(4 SNPs)
const uint8 cti[81] = {0, 1, 2, 4, 5, 6, 8, 9, 10, 16, 17, 18, 20, 21, 22, 24, 25, 26, 32, 33, 34, 36, 37, 38, 40, 41, 42, 64, 65, 66, 68, 69, 70, 72, 73, 74, 80, 81, 82, 84, 85, 86, 88, 89, 90, 96, 97, 98, 100, 101, 102, 104, 105, 106, 128, 129, 130, 132, 133, 134, 136, 137, 138, 144, 145, 146, 148, 149, 150, 152, 153, 154, 160, 161, 162, 164, 165, 166, 168, 169, 170};
const uint32 byte_in_word = sizeof(word);

template <uint32 EXP>
uint32 integerPow(uint32 x)
{
    return integerPow<EXP - 1>(x) * x;
}

template <>
uint32 integerPow<1>(uint32 x) { return x; }

template <>
uint32 integerPow<0>(uint32) { return 1; }