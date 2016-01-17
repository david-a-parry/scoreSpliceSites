#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use ScoreSpliceSite;

my $seq = 'AATgtatccttgtggttgtcaga';#'acaatatccttttta';
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
        site => 'D',
        species => 9606,
    );

    print "$type = $score\n";
}

