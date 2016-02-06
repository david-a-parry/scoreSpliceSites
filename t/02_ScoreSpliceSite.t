use strict;
use warnings;
use Test::More tests => 14;

BEGIN 
{ 
    use_ok("ScoreSpliceSite");
}
#precalculated donor scores for sequence below
my $dseq = 'AATgtatccttgtggttgtcaga';#'acaatatccttttta';
my %donors = 
(
    AT_AC_U12 => 49.9565371851927,
    GC_AG_U2 => 30.9109949942965,
    GT_AG_U12 => 92.0290805124072,
    GT_AG_U2 => 54.1714920340292,
);

#check score for each type
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

#precalculated acceptor scores for sequence below
my $aseq = "tgttttctcaacagGTTATTATTGCAGCAGAA" ;#17 nt , 12 upstream, splice site and 3 downstream
my %acceptors = 
(
    AT_AC_U12 => 39.9764369292882,
    GC_AG_U2 => 83.8341905752822,
    GT_AG_U12 => 74.5219436575304,
    GT_AG_U2 => 82.6019538530399,
);
#check score for each type
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
        "score $type for $aseq",
    );
}

#dummy intronic sequence containing 'perfect' U12 branch site
my $u12_branch = 'aaaaaaaacacagataaaaaaaaccaatttttccttaactaaaaaaaaaagagcgtgttttctcaacagGTTATTATTGCAGCAGAA';

my ($score, $bseq) = ScoreSpliceSite::scanForBranchPoint
(
    type => 'AT_AC_U12',
    seq => $u12_branch,
    species => 9606,
);
is (lc($bseq), 'ttttccttaact', "find perfect U12 branch point");
is ($score, 100, "score perfect U12 branch point 100");

my $worst = ScoreSpliceSite::getWorstSeq("AT_AC_U12", 'B', 9606);
($score, $bseq) = ScoreSpliceSite::scanForBranchPoint
(
    type => 'AT_AC_U12',
    seq => $worst,
    species => 9606,
);
is($score, 0, "score worst U12 branch site 0");

#dummy intronic sequence containing 'perfect' U2 branch site
my $u2_branch = 'AAAAAAAAAAAAAAAAAAAAAATTCTCATTCAAAAAAAAAAAAAAAAAAAAAAAAAA';
($score, $bseq) = ScoreSpliceSite::scanForBranchPoint
(
    type => 'GT_AG_U2',
    seq => $u2_branch,
    species => 9606,
);
is (lc($bseq), 'ttctcattc', "find perfect U2 branch point");
is ($score, 100, "score perfect U2 branch point 100");

