use strict;
use warnings;
use Test::More tests => 6;

require_ok("Bio::Tools::GFF");
require_ok("Bio::DB::Sam");
require_ok("Bio::SeqFeature::Generic");
require_ok("Sort::External");
require_ok("ScoreSpliceSite");
require_ok("ReverseComplement");
