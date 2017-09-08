#!/usr/bin/perl -w
# usage: replay.pl
# reading the listed run numbers from a text file (replaylist.txt)
# form the file name as (BASE + RUN + SUFF)
# get the raw data files from mss and then replay it with a specified program (RE_FUNC)

$FUNC = "/work/hallb/prad/PRadAnalyzer/examples/bin/NNAE_Prepare2";
$IN_DIR = "/work/hallb/prad/replay/event_sel";
$OUT_DIR = "/work/hallb/prad/NNAE_Files";

sub get_input {
    return $IN_DIR."/prad_".$_[0]."_sel.dst";
}

sub get_output {
    return $OUT_DIR."/prad_".$_[0]."_NNAE2.root";
}

open (INPUT, "<./runlist.txt");

foreach $i (<INPUT>) {
    chomp($i);
    $file = get_input($i);
    $outf = get_output($i);

    if(-f $outf)
    {
        print "run $i is already replayed, file exists at $outf.\n";
        next;
    }

    $cmd = $FUNC." ".$file." ".$outf;
    print "$cmd \n";
    system("$cmd");
}
close(INPUT);
