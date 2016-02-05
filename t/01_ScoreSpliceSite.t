use strict;
use warnings;
use Test::More tests => 9;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";

BEGIN 
{ 
    use_ok("ScoreSpliceSite");
}

my $dseq = 'AATgtatccttgtggttgtcaga';#'acaatatccttttta';

my %donors = 
(
    AT_AC_U12 => 49.9565371851927,
    GC_AG_U2 => 30.9109949942965,
    GT_AG_U12 => 92.0290805124072,
    GT_AG_U2 => 54.1714920340292,
);

foreach my $type (keys %donors){
    my $score = ScoreSpliceSite::score
    (
        seq  => $dseq,
        type => $type,
        site => 'D',
        species => 9606,
    );
    is
    (
        $score,
        $donors{$type},
        "score $type for $dseq",
    );
}

my $aseq = "tgttttctcaacagGTTATTATTGCAGCAGAA" ;#17 nt , 12 upstream, splice site and 3 downstream
my %acceptors = 
(
    AT_AC_U12 => 39.9764369292882,
    GC_AG_U2 => 83.8341905752822,
    GT_AG_U12 => 74.5219436575304,
    GT_AG_U2 => 82.6019538530399,
);
foreach my $type (keys %acceptors){
    my $score = ScoreSpliceSite::score
    (
        seq  => $aseq,
        type => $type,
        site => 'A',
        species => 9606,
    );
    is
    (
        $score,
        $acceptors{$type},
        "score $type for $aseq...",
    );
}
