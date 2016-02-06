#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw( min max );
use POSIX qw/strftime/;
use Bio::Tools::GFF;
use FindBin qw($RealBin);
use Bio::DB::Sam;
use Bio::SeqFeature::Generic;
use lib "$RealBin/lib";
use ScoreSpliceSite;
use ReverseComplement qw/ reverse_complement /;

my %opts = ();
GetOptions(
    \%opts,
    'f|fasta=s',
    'g|gff=s',
    'o|output=s',
    'h|?|help',
) or usage("Error getting options!");
usage() if $opts{h};
usage("-f/--fasta argument is required.") if not $opts{f};
usage("-g/--gff argument is required.") if not $opts{g};

my $IN;
if ($opts{g} =~ /\.gz$/){
    open ($IN, "gzip -dc $opts{g} |") or die "Can't open $opts{g} via gzip: $!\n";
}else{
    open ($IN, "<", $opts{g}) or die "Can't open $opts{g} for reading: $!\n";
}

my $fai = Bio::DB::Sam::Fai->load($opts{f});#should create index if it doesn't exist
my $gff = Bio::Tools::GFF->new
(
    -gff_version => 3,
    -fh          => $IN,
);

my $OUT = \*STDOUT;
if ($opts{o}){
    open ($OUT, ">", $opts{o}) or die "Can't open $opts{o} for writing: $!\n";    
}
print $OUT join
(   
    "\t",
    qw(
        INTRON_TYPE
        SUBTYPE
        EXON_ID
        PERCENT_GC
        TOTAL_REPEAT_LENGTH
        LONGEST_REPEAT
        TOTAL_HOMOPOLYMER_LENGTH
        LONGEST_HOMOPOLYMER
    )
) . "\n";
my %exon_seqs = ();
my @introns = (); 
my %names = (); 
while (my $feat = $gff->next_feature() ) {
    if($feat->has_tag('gene_id')){
        #parse previous introns
        parseIntrons();
        #clear collected exon seqs and introns
        %exon_seqs = (); 
        @introns = ();
        #get name and gene id
        my ($id) = $feat->get_tag_values('ID'); 
        $id =~ s/^gene://;
        my $name = '.';
        if ($feat->has_tag('Name')){
            ($name) = $feat->get_tag_values('Name'); 
        }
        $names{$id} = $name; 
    }elsif ($feat->has_tag('transcript_id')){
        #parse introns in case we missed a gene_id tag
        parseIntrons();
        #...and clear collected exons in case we missed a gene_id tag
        %exon_seqs = (); 
        @introns = ();
        #get transcript name and associate with gene id
        my ($tr) = $feat->get_tag_values('transcript_id'); 
        my ($parent) = $feat->get_tag_values('Parent');
        my $biotype = join(",", $feat->get_tag_values('biotype'));
    }elsif ($feat->primary_tag eq 'exon'){
        #collect exons 
        getExonSequence($feat);
    }elsif ($feat->primary_tag eq 'intron'){
        #collect intron type related to each exon
        push @introns, $feat;
    }
}
parseIntrons();

$gff->close();
close $OUT;

#################################################
sub parseIntrons{
    #if an exon has a neighbouring U12 intron classify as U12
    #if an exon has no neighbouring U12 introns and its neighbouring
    # introns are confirmed U2 classify as U2
    #if an exon has a neighbouring 'unknown' intron and does not have 
    # a neighbouring U12 intron classify as UNKNOWN
    my %u12 = ();
    my %u2  = ();
    my %unknown = ();
    my %all = ();
    foreach my $intr (@introns){
        my ($next)     = $intr->get_tag_values('next_exon_id');
        my ($previous) = $intr->get_tag_values('previous_exon_id');
        my ($type)     = $intr->get_tag_values('intron_type'); 
        if ($type =~ /U12$/){
            $u12{$next}     = $type;
            $u12{$previous} = $type;
        }elsif($type =~ /U2$/){
            $u2{$next}     = $type;
            $u2{$previous} = $type;
        }else{
            $unknown{$next}     = $type;
            $unknown{$previous} = $type;
        }
        $all{$next}     = undef;
        $all{$previous} = undef;
    }
    foreach my $k (%all){
        my ($class, $subclass);
        if (exists $u12{$k}){
            $class = 'U12';
            $subclass = $u12{$k};
        }elsif (exists $unknown{$k}){
            $class = 'UNKNOWN';
            $subclass = $unknown{$k};
        }else{
            $class = 'U2';
            $subclass = $u2{$k};
        }
        writeExonStats($class, $subclass, $k);
    }
}

#################################################
sub writeExonStats{
    my ($class, $subclass, $exon) = @_;
    my $gc = getGcPercentage($exon_seqs{$exon});
    my @r = getRepeats($exon_seqs{$exon});
    my ($longest, $longest_homo, $homo, $total) = (0, 0, 0, 0);
    foreach my $rep (@r){
        my $l = length($rep);
        $longest = $l > $longest ? $l : $longest;
        $total += $l;
        if ($rep =~ /^(\w)(\1+)+$/){
            $homo += $l;
            $longest_homo = $l > $longest_homo ? $l : $longest_homo;
        }
    }
    print $OUT join
    (
        "\t",
        $class,
        $subclass,
        $exon,
        sprintf("%.3f", $gc),
        $total,
        $longest,
        $homo,
        $longest_homo,
    ) . "\n";
    #intron type (U2 or U12 or UNKNOWN), subtype, exon ID, % GC,
    # total length of repeats, longest repeat, total homopolymer length, 
    # longest homopolymer
}

#################################################
sub getRepeats{
#return array of 1-4 nucleotide repeats at least 3 nt long found in $seq
    my $seq = shift;
    my @h = ();
    while ($seq =~ /(\w{1,4})(\1)(\1+)*/g){
        my $rep = "$1$2";
        $rep .= $3 if $3;
        push @h, $rep if length($rep) > 2;
    }
    return @h;
}

#################################################
sub getGcPercentage{
    my $seq = shift;
    my $n = 0;
    my $gc = 0;
    while ($seq =~ /(.)/g){
        if ($1 =~ /[cgCG]/){
            $gc++;
        }
        $n++;
    }
    return 0 if $n == 0; 
    return 100 * $gc / $n;
}

#################################################
sub getExonSequence{
    my $exon = shift;
    my $chrom = $exon->seq_id; 
    my $strand = $exon->strand;
    my $start = $exon->start;
    my $end = $exon->end;
    my $seq = $fai->fetch("$chrom:$start-$end");
    if ($strand < 0){
        $seq = reverse_complement($seq);
    }
    my ($id) = $exon->get_tag_values('exon_id');
    $exon_seqs{$id} = $seq;
}

#################################################
sub usage{
    my $msg = shift;
    print "\n$msg\n" if $msg;

    print <<EOT

TODO 

USAGE: $0 -f genome_fasta.fa -g introns.gff3

OPTIONS:
    
    -f,--fasta FILE
        Genome fasta file for retrieving DNA sequences for intron-exon boundaries

    -g,--gff FILE
        GFF3 intron file created by spliceScorer.pl

    -o,--output FILE
        Optional output file. Default = STDOUT.

    -h,--help 
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}




