#!/bin/bash
# usage: event_sel.sh <begin_run> <end_run>
# auto search the replayed dst files and event selection list to cut off bad events

# define program position
prog="/work/hallb/prad/PRadAnalyzer/examples/bin/eventSelect"
# define replay files directory
dir="/lustre/expphy/work/hallb/prad/replay/"
# define file format
file_format="prad_[run].dst"
# define event list
bad_base="/lustre/expphy/work/hallb/prad/replay_EnS/EventSelect_[run].txt"

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
if [[ $1 -eq 0 ]]; then
    run_begin=0
else
    run_begin=$1
fi

# end
if [[ $2 -eq 0 ]]; then
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
    bad_file=${bad_base//\[run\]/$run}
    if [ "$run" -ge "$run_begin" ] && [ "$run" -le "$run_end" ]; then
        echo $prog $file $bad_file
        $prog $file $bad_file
    fi
done
