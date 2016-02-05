#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use ScoreSpliceSite;

my $seq = "tgttttctcaacagGTTATTATTGCAGCAGAA" ;#17 nt , 12 upstream, splice site and 3 downstream
$seq = shift if @ARGV; 
foreach my $type ( qw 
/
    AT_AC_U12  
    GC_AG_U2  
    GT_AG_U12  
    GT_AG_U2
/){
    my $score = ScoreSpliceSite::score
    (
        seq  => $seq,
        type => $type,
        site => 'A',
        species => 9606,
    );

    print "$type = $score\n";
}

