#!/bin/bash
# change variables below to your evio library and include file path
EVLIBDIR=/site/coda/3.02/Linux-x86_64/lib
EVINCDIR=/site/coda/3.02/common/include

gcc -o evSplitter main.c -I$EVINCDIR -L$EVLIBDIR -levio -lpthread
