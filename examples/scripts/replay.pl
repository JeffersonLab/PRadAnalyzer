#!/usr/bin/perl -w
# usage: replay.pl
# reading the listed run numbers from a text file (replaylist.txt)
# form the file name as (BASE + RUN + SUFF)
# get the raw data files from mss and then replay it with a specified program (RE_FUNC)

$MSS = "/mss/hallb/prad/data/";
$PATH = "/work/hallb/prad/data/";
$RE_FUNC = "/work/hallb/prad/PRadAnalyzer/examples/bin/replay";
$RE_PATH = "/work/hallb/prad/replay/";
$BASE = "prad_";
$SUFF = ".evio";

open (INPUT, "<./replaylist.txt");

foreach $i (<INPUT>) {
    chomp($i);
    $file = $BASE. $i;
    $re_file = $RE_PATH. $file. ".dst";

    if(-f $re_file)
    {
        print "run $i is already replayed, file exists at $re_file.\n";
        next;
    }
    $jget = "jget ". $MSS. $file. $SUFF. ".*". " ". $PATH;
    print "$jget \n";
    $job = system("$jget");
    print "$job \n";

    $replay = $RE_FUNC." ". $PATH. $file. $SUFF." ". $re_file." -s 1500";
    print "$replay \n";
    $job = system("$replay");
    print "$job \n";

    if(-f $re_file)
    {
        $del = "rm -rf ". $PATH. $file. $SUFF. ".*";
        print "$del \n";
        $job = system("$del");
        print "$job \n";
    }
}

close(INPUT);
