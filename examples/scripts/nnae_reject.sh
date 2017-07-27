#!/bin/bash
# usage: event_sel.sh <old_dir> <new_dir>
# auto search the replayed dst files and event selection list to cut off bad events

# define program position
prog=${PRAD_PATH}"/examples/bin/NNAE_Reject"
# root file position
root_dir="/lustre/expphy/work/hallb/prad/NNAE_Files/output/files_with_kin_cuts/"
root_format="PRad_ID_[run]_NNAE_[cut].root"
# dst file position
dst_dir="/lustre/expphy/work/hallb/prad/replay/event_sel/"
dst_format="prad_[run]_sel.dst"
run_pattern="([[:digit:]]{6})"

# get output dir
if [ -z $1 ]; then
    echo "usage: nnae_reject.sh <cut_type>"
    exit
else
    cut_type=$1
fi

if [ "$cut_type" = "cut1" ]; then
    cut_str="--cut-angle-min=5.0 --cut-energy-min=1800 --cut-energy-max=2500"
elif [ "$cut_type" = "cut2" ]; then
    cut_str="--cut-angle-min=2.0 --cut-angle-max=5.0 --cut-energy-min=1800 --cut-energy-max=2500"
else
    echo "unsupported cut type = $cut_type"
    exit
fi

root_format2=${root_format//"[cut]"/$cut_type}
search_pattern=${root_format2//"[run]"/$run_pattern}

# execute program
for file in $root_dir*; do
    fname=${file//$root_dir/}
    if [[ ! $fname =~ $search_pattern ]]; then
        continue
    fi
    run=${BASH_REMATCH[1]}
    dst_file=${dst_format//"[run]"/$run}
    echo $prog $dst_dir$dst_file $file $cut_type $cut_str
    $prog $dst_dir$dst_file $file $cut_type $cut_str
done
