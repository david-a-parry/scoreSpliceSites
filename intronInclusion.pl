#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use POSIX qw/strftime/;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;

my %opts = ();
GetOptions(
    \%opts,
    'g|gff=s',
    'o|output=s',
    'h|?|help',
) or usage("Error getting options!");
usage() if $opts{h};

usage("-g/--gff argument is required.") if not $opts{g};
my $IN;
if ($opts{g} =~ /\.gz$/){
    open ($IN, "gzip -dc $opts{g} |") or die "Can't open $opts{g} via gzip: $!\n";
}else{
    open ($IN, "<", $opts{g}) or die "Can't open $opts{g} for reading: $!\n";
}
my $gff = Bio::Tools::GFF->new
(
    -gff_version => 3,
    -fh          => $IN,
);


my $OUT = \*STDOUT;
if ($opts{o}){
    open ($OUT, ">", $opts{o}) or die "Can't open $opts{o} for writing: $!\n";    
}

my %gene = ();
print $OUT join("\t", 
qw/
    Gene
    Symbol
    Intron
    Intron_Type
    Transcripts
    U12_Transcripts
    Coding_Transcripts
    U12_Coding_Transcripts
    Fraction_U12_Transcripts
    Fraction_U12_Coding_Transcripts
/) . "\n";
 
while (my $feat = $gff->next_feature() ) {
   if($feat->has_tag('gene_id')){
        parseGene();
        %gene = ();
        my ($id) = $feat->get_tag_values('ID'); 
        $id =~ s/^gene://;
        my $name = '.';
        if ($feat->has_tag('Name')){
            ($name) = $feat->get_tag_values('Name'); 
        }
        $gene{name} = $name; 
        $gene{gene} = $id;
    }else{
        parseFeature($feat);
    }
}

#################################################
sub parseGene{
    return if not exists $gene{U12};
    my $n_trans  = scalar(keys %{$gene{transcripts}});
    my $n_coding = 0; 
    foreach my $tr (keys %{$gene{transcripts}}){
        $n_coding++ if $gene{transcripts}->{$tr}->{biotype} =~ /protein_coding/;
    }
    foreach my $type (keys %{$gene{U12}}){
        foreach my $u12intron (keys %{$gene{U12}->{$type}}){
            my $u12_transcripts = 0;
            my $u12_coding_transcripts = 0;
            foreach my $tr (keys %{$gene{transcripts}}){
                if (exists $gene{transcripts}->{$tr}->{introns}->{$u12intron}){
                    $u12_transcripts++;
                    if ($gene{transcripts}->{$tr}->{biotype} =~ /protein_coding/){
                        $u12_coding_transcripts++ ;
                    }
                }
            }
            my $u12_frac;
            my $u12_coding_frac;
            eval
            {
                $u12_frac = $u12_transcripts / $n_trans;
            };
            eval
            {
                $u12_coding_frac = $u12_coding_transcripts / $n_coding;
            };
            $u12_frac ||= 0;
            $u12_coding_frac ||= 0;
            print $OUT join("\t",
                $gene{gene},
                $gene{name},
                $u12intron,
                $type,
                $n_trans,
                $u12_transcripts,
                $n_coding,
                $u12_coding_transcripts,
                $u12_frac,
                $u12_coding_frac,
            ) . "\n"; 
        }
    }
}

#################################################
sub parseFeature{
    my $feat = shift;
    
    die "Attempt to parse non-gene feature without any gene entries! "
      if not keys %gene;
    
    if ($feat->has_tag('transcript_id')){
        my ($tr) = $feat->get_tag_values('transcript_id'); 
        my ($parent) = $feat->get_tag_values('Parent');
        $parent =~ s/^gene://;
        die "Transcript parent ($parent) does not match current gene ($gene{gene}) " 
          if $parent ne $gene{gene};
        $gene{transcripts}->{$tr}->{biotype} = join(",", $feat->get_tag_values('biotype'));
    }elsif ($feat->primary_tag eq 'exon'){
        #DO WE ACTUALLY NEED ANY EXON INFO?
#        foreach my $tr ($feat->get_tag_values('Parent') ){
#            $tr =~ s/transcript://;
#            my $name = $feat->get_tag_values('Name');
#            die "Encountered exon ($name) for unencountered transcript: $tr "
#              if not exists  $gene{transcripts}->{$tr};
#            foreach my $ex ($feat->get_tag_values('rank') ){
#                $gene{transcripts}->{$tr}->{exons}->{$name} = $ex; 
#            }
#        }
    }elsif ($feat->primary_tag eq 'intron'){
        #my ($prev) = $feat->get_tag_values('previous_exon_id');
        #my ($next) = $feat->get_tag_values('next_exon_id');
        my ($i)    = $feat->get_tag_values('previous_rank');
        #my $intron = "$prev-$next";
        my $start = $feat->start;
        my $end   = $feat->end;
        my $chrom = $feat->seq_id;
        my $intron = "$chrom:$start-$end";
        my ($type) =  $feat->get_tag_values('intron_type');
        if ($type =~ /\_U12$/){
            $gene{U12}->{$type}->{$intron} = undef;
        }
        foreach my $tr ($feat->get_tag_values('Parent') ){
            $tr =~ s/transcript://;
            $gene{transcripts}->{$tr}->{introns}->{$intron} = $i;
        }
    }
}

#################################################
sub usage{
    my $msg = shift;
    print "\n$msg\n" if $msg;

    print <<EOT

Count the number of U12 introns per transcript for U12 intron containing genes

USAGE: $0 -g genes_and_introns.gff3

OPTIONS:
    
    -g,--gff FILE
        GFF3 file generated using spliceScorer.pl

    -o,--output FILE
        Optional output file. Default = STDOUT.

    -h,--help 
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}



