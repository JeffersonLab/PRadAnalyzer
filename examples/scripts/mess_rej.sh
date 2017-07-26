#!/bin/bash
# usage: cosmic_rej.sh <begin_run> <end_run>
# auto search the replayed dst files and reject cosmic backgrounds by using a specified program

# define program position
prog="/work/hallb/prad/PRadAnalyzer/examples/bin/messReject"
# define replay files directory
dir="/lustre/expphy/work/hallb/prad/replay/event_sel/"
# define event file format
file_format="prad_[run]_sel.dst"

# pre-test
if [ ! -f "$prog" ]; then
    echo "$prog does not exist!"
    exit
fi
if [ ! -d "$dir" ]; then
    echo "$dir does not exist!"
    exit
fi

# check run number range
# begin
if [ -z $1 ]; then
    run_begin=0
else
    run_begin=$1
fi

# end
if [ -z $2 ]; then
    run_end=999999
else
    run_end=$2
fi

# execute program
search_pattern=${file_format//"[run]"/"[0-9]+"}
for file in $dir*; do
    fname=${file//$dir/}
    if [[ ! $fname =~ $search_pattern ]]; then
        continue
    fi
    run=`echo $fname | egrep -o "[0-9]+"`
    if [ "$run" -ge "$run_begin" ] && [ "$run" -le "$run_end" ]; then
        echo $prog $file
        $prog $file
    fi
done
