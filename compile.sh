#!/bin/bash
set -x

g++ -o BitEpi.o -O3 BitEpi.cpp csvparser.c -pthread
chmod 777 BitEpi.o

exit
