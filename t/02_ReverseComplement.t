use strict;
use warnings;
use Test::More tests => 4;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";

BEGIN 
{ 
    use_ok("ReverseComplement", qw /reverse_complement complement/);
}

my $seq     = 'AATgtatccttgtggttgtcaga';
my $rev     = 'agactgttggtgttcctatgTAA';
my $comp    = 'TTAcataggaacaccaacagtct';
my $revcomp = 'tctgacaaccacaaggatacATT';

is(
    reverse($seq),
    $rev,
    "reverse sequence"
);

is(
    complement($seq),
    $comp,
    "complement sequence"
);

is(
    reverse_complement($seq),
    $revcomp,
    "reverse complement sequence"
);

