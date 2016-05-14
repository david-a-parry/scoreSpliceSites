#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use POSIX qw/strftime/;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;

my %opts = ();
my $intron_mode = 'U12';
GetOptions(
    \%opts,
    'g|gff=s',
    't|tsl=i',
    'o|output=s',
    'u|u2mode',
    'h|?|help',
) or usage("Error getting options!");
usage() if $opts{h};
usage("-g/--gff argument is required.") if not $opts{g};
$intron_mode = 'U2' if $opts{u};
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
my $header = join("\t", 
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
if ($opts{u}){
    $header =~ s/U12/U2/g;
}
print $OUT $header;
 
while (my $feat = $gff->next_feature() ) {
    if ($feat->primary_tag eq 'gene'){
    #if($feat->has_tag('gene_id')){
        parseGene();
        %gene = ();
        my ($id) = $feat->get_tag_values('ID'); 
        $id =~ s/^gene://;
        my $name = '.';
        if ($feat->has_tag('Name')){
            ($name) = $feat->get_tag_values('Name'); 
        }elsif ($feat->has_tag('gene_name')){
            ($name) = $feat->get_tag_values('gene_name'); 
        }
        $gene{name} = $name; 
        $gene{gene} = $id;
    }else{
        parseFeature($feat);
    }
}

#################################################
sub parseGene{
    return if not exists $gene{$intron_mode};
    my @valid_tr = ();
    my @valid_coding_tr = (); 
    foreach my $tr (keys %{$gene{transcripts}}){
        if (defined $opts{t}){
            next if $gene{transcripts}->{$tr}->{tsl} eq 'NA';
            next if  $gene{transcripts}->{$tr}->{tsl} > $opts{t};
        }
        push @valid_tr, $tr;
        push @valid_coding_tr, $tr if $gene{transcripts}->{$tr}->{gene_type} =~ /protein_coding/;
    }
    return if not @valid_tr;
    foreach my $type (keys %{$gene{$intron_mode}}){
        foreach my $u12intron (keys %{$gene{$intron_mode}->{$type}}){
            my $u12_transcripts = 0;
            my $u12_coding_transcripts = 0;
            foreach my $tr (keys %{$gene{transcripts}}){
                #filter on transcript support level if specified
                if (defined $opts{t}){
                    next if $gene{transcripts}->{$tr}->{tsl} eq 'NA';
                    next if  $gene{transcripts}->{$tr}->{tsl} > $opts{t};
                }
                if (exists $gene{transcripts}->{$tr}->{introns}->{$u12intron}){
                    $u12_transcripts++;
                    if ($gene{transcripts}->{$tr}->{gene_type} =~ /protein_coding/){
                        $u12_coding_transcripts++ ;
                    }
                }
            }
            next if not $u12_transcripts; #no transcript meeting our required transcript support level 
            my $u12_frac;
            my $u12_coding_frac;
            eval
            {
                $u12_frac = $u12_transcripts / @valid_tr;
            };
            eval
            {
                $u12_coding_frac = $u12_coding_transcripts / @valid_coding_tr;
            };
            $u12_frac ||= 0;
            $u12_coding_frac ||= 0;
            print $OUT join("\t",
                $gene{gene},
                $gene{name},
                $u12intron,
                $type,
                scalar(@valid_tr),
                $u12_transcripts,
                scalar(@valid_coding_tr),
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
    
    if ($feat->primary_tag eq 'transcript'){
    #if ($feat->has_tag('transcript_id')){
        my ($tr) = $feat->get_tag_values('ID'); 
        my ($parent) = $feat->get_tag_values('Parent');
        $parent =~ s/^gene://;
        die "Transcript parent ($parent) does not match current gene ($gene{gene}) " 
          if $parent ne $gene{gene};
        $gene{transcripts}->{$tr}->{gene_type} = join(",", $feat->get_tag_values('gene_type'));
        eval
        {
            ($gene{transcripts}->{$tr}->{tsl}) = $feat->get_tag_values('transcript_support_level');
        };#not all transcripts have tsl tag(?)
        $gene{transcripts}->{$tr}->{tsl} ||= 'NA';
        #clear cruft from tsl
        $gene{transcripts}->{$tr}->{tsl} =~ s/\s+.*//;
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
        my ($i)    = $feat->get_tag_values('previous_exon_number');
        #my $intron = "$prev-$next";
        my $start = $feat->start;
        my $end   = $feat->end;
        my $chrom = $feat->seq_id;
        my $intron = "$chrom:$start-$end";
        my ($type) =  $feat->get_tag_values('intron_type');
        if ($type =~ /_$intron_mode$/){
            $gene{$intron_mode}->{$type}->{$intron} = undef;
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

    -t,--tsl INT
        Maximum transcript support level for inclusion in analysis. By Default
        all transcripts will be analyzed. Specify the maximum value for 
        inclusion here (transcripts with a tsl of 1 have the most support, 
        those with a tsl of 5 have the least. See Ensembl's help:
        http://www.ensembl.org/Help/Glossary?id=492).

    -u,--u2mode
        Use this flag to output counts for U2 introns instead of U12

    -o,--output FILE
        Optional output file. Default = STDOUT.

    -h,--help 
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}



