#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use POSIX qw/strftime/;
use Bio::Tools::GFF;
use FindBin qw($RealBin);
use Bio::DB::Sam;
use Bio::SeqFeature::Generic;
use lib "$RealBin/lib";
use ScoreSpliceSite;

my %opts = (s => 9606);
GetOptions(
    \%opts,
    'f|fasta=s',
    'g|gff=s',
    's|species=s',
    'o|output=s',
    't|transcripts=s',
    'h|?|help',
) or usage("Error getting options!");
usage() if $opts{h};
usage("-f/--fasta argument is required.") if not $opts{f};
usage("-g/--gff argument is required.") if not $opts{g};

if (not grep {$opts{s} eq $_} ScoreSpliceSite::getSpecies){
    die "Splice consensus for species '$opts{s}' not available.\n";
}

my $fai = Bio::DB::Sam::Fai->load($opts{f});#should create index if it doesn't exist
my $gff = Bio::Tools::GFF->new
(
    -gff_version => 3,
    -file        => $opts{g},
);

my $OUT = \*STDOUT;
if ($opts{o}){
    open ($OUT, ">", $opts{o}) or die "Can't open $opts{o} for writing: $!\n";    
}

my $gffwriter = Bio::Tools::GFF->new(
    -gff_version => 3,
    -fh        => $OUT,
);

my $TRANS;
if ($opts{t}){
    open ($TRANS, ">", "$opts{t}") or die "Could not open $opts{t} for writing: $!\n";
    print $TRANS "#NAME\tGENE_ID\tTRANSCRIPT_ID\tTRANSCRIPT_BIOTYPE\tU12_introns\tU2_introns\tUNKNOWN_introns\n";
}

my %exons = ();
my %transcripts = ();
my %names = ();
my $i = 0;
while (my $feat = $gff->next_feature() ) {
    #if ($feat->primary_tag =~ /^(gene|processed_transcript|miRNA_gene|RNA)$/){
    if($feat->has_tag('gene_id')){
        parseExons(\%exons);
        %exons = ();
        my ($id) = $feat->get_tag_values('ID'); 
        $id =~ s/^gene://;
        my ($name) = $feat->get_tag_values('Name'); 
        $names{$id} = $name; 
        $gffwriter->write_feature($feat);
    }elsif ($feat->has_tag('transcript_id')){
        my ($tr) = $feat->get_tag_values('transcript_id'); 
        ($transcripts{$tr}->{parent}) = $feat->get_tag_values('Parent');
        $transcripts{$tr}->{biotype} = join(",", $feat->get_tag_values('biotype'));
        parseExons(\%exons);
        %exons = ();
        #foreach my $tr ($feat->get_tag_values('ID') ){
        #    $transcripts{$tr} = $feat;
        #}
        $gffwriter->write_feature($feat);
    }elsif ($feat->primary_tag eq 'exon'){
        foreach my $tr ($feat->get_tag_values('Parent') ){
            $tr =~ s/transcript://;
            foreach my $ex ($feat->get_tag_values('rank') ){
                $exons{$tr}->{$ex} = $feat; 
            }
        }
#debug    }elsif ($feat->primary_tag ne 'CDS' and $feat->primary_tag !~ /UTR/){
#debug        print Dumper $feat;
    }
}
parseExons(\%exons);
my $time = strftime( "%H:%M:%S", localtime );
print STDERR "[$time] Finished - processed $i introns\n";
if ($TRANS){
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] writing transcript intron counts to $opts{t}...\n";
    writeTranscriptCounts();
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished - processing " . scalar(keys %transcripts) . " transcripts\n";
}

#################################################
sub writeTranscriptCounts{
    foreach my $k (sort keys %transcripts){
        my ($unknown, $u2, $u12, $gene, $biotype) = (0, 0, 0, '.', '.'); 
        foreach my $type (keys %{$transcripts{$k}}){
            if ($type eq 'parent'){
                ($gene = $transcripts{$k}->{$type}) =~ s/^gene://;
            }elsif($type eq 'biotype'){
                $biotype = $transcripts{$k}->{$type};
            }elsif($type eq '0'){
                $unknown += $transcripts{$k}->{$type};
            }elsif($type =~ /U12$/){
                $u12 += $transcripts{$k}->{$type};
            }elsif($type =~ /U2$/){
                $u2 += $transcripts{$k}->{$type};
            }else{
                warn "Don't recognise hash key '$type' for $k in transcripts hash!\n";
            }
        }
        print $TRANS join("\t", 
            $names{$gene},
            $gene,
            $k,
            $biotype,
            $u12,
            $u2,
            $unknown
        ) . "\n";
    }
}
#################################################
sub parseExons{
    my $exons = shift;
    my $n = 0;
    foreach my $tr (keys %$exons){
        foreach my $ex (sort {$a <=> $b} keys %{$exons->{$tr}}){
            if (exists $exons{$tr}->{$ex-1}){
                writeIntron
                (
                    $exons{$tr}->{$ex-1}, 
                    $exons{$tr}->{$ex}, 
                    $tr,
                );
                $i++;
                reportProgress($i);
            }
            $gffwriter->write_feature($exons{$tr}->{$ex});
        }
    }
}

#################################################
sub reportProgress{
    my $n = shift;
    return if not $n;
    return if ($n % 10000); 
    my $time = strftime( "%H:%M:%S", localtime );
    $n =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g; #add commas for readability
    print STDERR "[$time] processed $n introns\n";
}

#################################################
sub writeIntron{
    my ($exon1, $exon2, $tr) = @_;

    my $chrom = $exon1->seq_id; 
    my $strand = $exon1->strand;
    my ($intron_start, $intron_stop); 
    if ($strand > 0){
        $intron_start = $exon1->end + 1; 
        $intron_stop = $exon2->start - 1;
    }else{
        ($exon2, $exon1) = ($exon1, $exon2);
        $intron_start = $exon1->start - 1; 
        $intron_stop = $exon2->end + 1;
    }
    my $intron = new Bio::SeqFeature::Generic
    (
        -seq_id      => $chrom,
        -start       => $intron_start,
        -end         => $intron_stop,
        -strand      => $strand,
        -source_tag  => 'spliceScorer',
        -primary_tag => 'intron',
    
    );
    
    $intron->add_tag_value
    (
        'Parent',
        $exon1->get_tag_values('Parent'),
    );
    foreach my $t 
    (qw /
        exon_id
        rank
        version
        ensembl_phase
        constitutive
        ensembl_end_phase
    /){
        $intron->add_tag_value
        (
            "previous_$t",
            $exon1->get_tag_values($t),
        ); 
        $intron->add_tag_value
        (
            "next_$t",
            $exon2->get_tag_values($t),
        ); 
    }
    my ($d_start, $d_end) = sort {$a <=> $b} 
    (
        ($intron_start - (3 * $strand) ) ,
        ($intron_start + (10 * $strand) ) ,
    );
    my ($a_start, $a_end) = sort {$a <=> $b} 
    (
        ($intron_stop - (13 * $strand)) ,
        ($intron_stop + (3 * $strand) ) ,
    );
    my $donor    = $fai->fetch("$chrom:$d_start-$d_end");
    my $acceptor = $fai->fetch("$chrom:$a_start-$a_end");

    if ($strand < 0){
        $donor = revcomp($donor);
        $acceptor = revcomp($acceptor);
    }
    $intron->add_tag_value
    (
        'donor_seq',
        uc ( substr($donor, 0, 3) ) . 
        lc (substr($donor, 3, ) ),
    );
    $intron->add_tag_value
    (
        'acceptor_seq',
        lc ( substr($acceptor, 0, 14) ) . 
        uc (substr($acceptor, 14, ) ),
    );
    
    my %scores = ();
    foreach my $type (ScoreSpliceSite::getIntronTypes){
        $scores{'D'}->{$type} = ScoreSpliceSite::score
        (
            seq  => $donor,
            type => $type,
            site => 'D',
            species => $opts{s},
        );
        $scores{'A'}->{$type} = ScoreSpliceSite::score
        (
            seq  => $acceptor,
            type => $type,
            site => 'A',
            species => $opts{s},
        );
    }
    my $intron_type = 0;
    foreach my $type (ScoreSpliceSite::getIntronTypes){
        $intron->add_tag_value
        (
            "donor_score_$type",
            $scores{'D'}->{$type},
        );
        $intron->add_tag_value
        (
            "acceptor_score_$type",
            $scores{'A'}->{$type},
        );
        if ($scores{'A'}->{$type} > 50 and 
            $scores{'D'}->{$type} > 50){
            if ($intron_type){
                if ( $scores{'D'}->{$type} >= $scores{'D'}->{$intron_type}){
                #this is a better score than our previous designation
                    if ($type =~ /U12$/ and $intron_type =~ /U2$/){
                    #we only annotate as U12 if the donor site score is 
                    #at least 25 more than a U2 site
                        $intron_type = $type if $scores{'D'}->{$type} - $scores{'D'}->{$intron_type} >= 25;
                    }else{
                        $intron_type = $type;
                    }
                    #TODO - consider branch site
                }
            }else{
                $intron_type = $type;
            }
        }
    }
    $intron->add_tag_value("intron_type", $intron_type); 
    $gffwriter->write_feature($intron);
    $transcripts{$tr}->{$intron_type}++;
}

#################################################
sub revcomp{
    my $seq = shift;
    $seq =~ tr/acgtACGT/tgcaTGCA/;
    return reverse($seq);
}

#################################################
sub usage{
    my $msg = shift;
    print "\n$msg\n" if $msg;

    print <<EOT

Create a GTF3 file of introns scored for U2 and U12 splice sites

USAGE: $0 -f genome_fasta.fa -g genes_and_exons.gff3

OPTIONS:
    
    -f,--fasta 
        Genome fasta file for retrieving DNA sequences for intron-exon boundaries
    -g,--gff
        GTF3 file containing information on the genes and exons to use for intron file creation
    -s,--species
        Taxonomic code for species to use for splice prediction. Default is 9606 (human). 
        Available species are 10090 (mouse), 3702 (A. thaliana), 6239 (C. elegans), 7227 (D. melanogaster) and 9606 (human).
    -o,--output
        Optional output file. Default = STDOUT.
    -h,--help 
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}




