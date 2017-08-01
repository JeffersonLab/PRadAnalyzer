#!/bin/bash
# usage: event_sel.sh <old_dir> <new_dir>
# auto search the replayed dst files and event selection list to cut off bad events

# define program position
prog=${PRAD_PATH}"/examples/bin/updateDST"

# pre-test
if [ ! -f "$prog" ]; then
    echo "$prog does not exist!"
    exit
fi

# get input dir
if [ -z $1 ]; then
    echo "usage: event_sel.sh <old_dir> <new_dir>"
    exit
else
    in_dir=$1
fi
LEN=${#in_dir}-1
if [ "${in_dir:LEN}" != "/" ]; then
    in_dir=$in_dir"/"
fi

# get output dir
if [ -z $2 ]; then
    out_dir="usage: event_sel.sh <old_dir> <new_dir>"
    exit
else
    out_dir=$2
fi
LEN=${#out_dir}-1
if [ "${out_dir:LEN}" != "/" ]; then
    out_dir=$out_dir"/"
fi

# check dirs
if [ ! -d "$in_dir" ] || [ ! -d "$out_dir" ]; then
    echo "$in_dir or $out_dir does not exist!"
    exit
fi

# execute program
for file in $in_dir*.dst; do
    fname=${file//$in_dir/}
    if [ ! -f "$out_dir$fname" ]; then
        echo $prog $file $out_dir$fname
        $prog $file $out_dir$fname
    else
        echo "$out_dir$fname exists, abort overwriting."
    fi
done
