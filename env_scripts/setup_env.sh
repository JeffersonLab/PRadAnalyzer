#!/bin/bash
# You will need root and evio to build the program
# Setup the path for these dependencies first

#source /home/chao/root-6.06.00/bin/thisroot.sh

export PRAD_PATH=~/PRadAnalyzer

if [ `uname -m` == 'x86_64' ]; then
    export THIRD_LIB=$PRAD_PATH/thirdparty/lib64
else
    export THIRD_LIB=$PRAD_PATH/thirdparty/lib
fi

export ET_LIB=/home/chao/PRad/coda/Linux-x86_64/lib
export ET_INC=/home/chao/PRad/coda/Linux-x86_64/include

#export ET_LIB=$THIRD_LIB
#export ET_INC=$PWD/thirdparty/include

export PRAD_LIB=$PRAD_PATH/lib
export PRAD_INC=$PRAD_PATH/include

export LD_LIBRARY_PATH=$PRAD_LIB:$THIRD_LIB:$ET_LIB:$LD_LIBRARY_PATH
