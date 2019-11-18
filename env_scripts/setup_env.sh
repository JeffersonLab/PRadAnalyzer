#!/bin/bash
# You will need root and evio to build the program
# Setup the path for these dependencies first

#source /home/chao/root-6.06.00/bin/thisroot.sh

# get directory of this script
SOURCE="${BASH_SOURCE[0]}"
# resolve $SOURCE until the file is no longer a symlink
while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

export PRAD_PATH=$(dirname "$DIR")

if [ `uname -m` == 'x86_64' ]; then
    export THIRD_LIB=$PRAD_PATH/thirdparty/lib64
else
    export THIRD_LIB=$PRAD_PATH/thirdparty/lib
fi

export ET_LIB=/home/chao/PRad/coda/Linux-x86_64/lib
export ET_INC=/home/chao/PRad/coda/Linux-x86_64/include

#export ET_LIB=$THIRD_LIB
#export ET_INC=$PWD/thirdparty/include

PRAD_LIB=$PRAD_PATH/lib
PRAD_INC=$PRAD_PATH/include

export LD_LIBRARY_PATH=$PRAD_LIB:$THIRD_LIB:$LD_LIBRARY_PATH
