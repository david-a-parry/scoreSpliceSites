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

my $OUT;
if ($opts{o}){
    open ($OUT, ">", $opts{o}) or die "Can't open $opts{o} for writing: $!\n";    
}else{
    $OUT = \*STDOUT;
}

my $gffwriter = Bio::Tools::GFF->new(
    -gff_version => 3,
    -fh        => $OUT,
);

my $TRANS;
if ($opts{t}){
    open ($TRANS, ">", "$opts{t}") or die "Could not open $opts{t} for writing: $!\n";
    print $TRANS "NAME\tGENE_ID\tTRANSCRIPT_ID\tTRANSCRIPT_BIOTYPE\tTSL\tU12_introns\tU2_introns\tUNKNOWN_introns\tLENGTH\tSPLICED_LENGTH\n";
}

my %exons = ();
my %intron_counts = (U12 => 0, U2 => 0, unknown => 0);
my %transcripts = ();
my %names = ();
my $i = 0;
while (my $feat = $gff->next_feature() ) {
    #print STDERR "TAG: " . $feat->primary_tag . "\n";
    #if ($feat->primary_tag =~ /^(gene|processed_transcript|miRNA_gene|RNA)$/){
    if ($feat->primary_tag eq 'gene' or #works for all gencode
       ( $feat->primary_tag =~ /gene|RNA|processed_transcript/ 
        and $feat->has_tag('gene_id') ) #for ensembl GFFs
    ){
        parseExons(\%exons);
        %exons = ();
        my ($id) = $feat->get_tag_values('ID'); 
        $id =~ s/^gene://;
        my $name = '.';
        if ($feat->has_tag('Name')){
            ($name) = $feat->get_tag_values('Name'); 
        }elsif ($feat->has_tag('gene_name')){
            ($name) = $feat->get_tag_values('gene_name'); 
        }
        $names{$id} = $name; 
        $gffwriter->write_feature($feat);
    }elsif ($feat->primary_tag eq 'transcript' or #works for gencode
           ( $feat->primary_tag =~ /gene|transcript|RNA/ 
             and $feat->has_tag('transcript_id') ) #for ensembl GFFs
    ){
        my ($tr) = $feat->get_tag_values('ID'); 
        $tr =~ s/transcript://;
        ($transcripts{$tr}->{parent}) = $feat->get_tag_values('Parent');
        $transcripts{$tr}->{length} = 1 + $feat->end - $feat->start;
        if ($feat->has_tag('gene_type')){
            $transcripts{$tr}->{gene_type} = join(",", $feat->get_tag_values('gene_type'));
        }elsif ($feat->has_tag('biotype')){
            $transcripts{$tr}->{gene_type} = join(",", $feat->get_tag_values('biotype'));
        }else{
            die "Could not determine biotype/gene_type for $tr!\n";
        }

        eval{
            ($transcripts{$tr}->{tsl}) = $feat->get_tag_values('transcript_support_level')
        };#not all transcripts have tsl tag(?)
        $transcripts{$tr}->{tsl} ||= 'NA';
        #clear cruft from tsl
        $transcripts{$tr}->{tsl} =~ s/\s+.*//;

        parseExons(\%exons);
        %exons = ();
        #foreach my $tr ($feat->get_tag_values('ID') ){
        #    $transcripts{$tr} = $feat;
        #}
        $gffwriter->write_feature($feat);
    }elsif ($feat->primary_tag eq 'exon'){
        foreach my $tr ($feat->get_tag_values('Parent') ){
            $tr =~ s/transcript://;
            my $ex_tag = '';
            if ($feat->has_tag('exon_number') ){
                $ex_tag = 'exon_number';#gencode GFFs
            }elsif ($feat->has_tag('rank') ){
                $ex_tag = 'rank';#ensembl GFFs
            }else{
                die "Could not determine exon number for ". $feat->{ID} . "\n";
            }
            foreach my $ex ($feat->get_tag_values($ex_tag) ){
                $exons{$tr}->{$ex} = $feat; 
            }
        }
#debug    }elsif ($feat->primary_tag ne 'CDS' and $feat->primary_tag !~ /UTR/){
#debug        print Dumper $feat;
    }
}
$gff->close();
parseExons(\%exons);
$gffwriter->close();
my $time = strftime( "%H:%M:%S", localtime );
$i =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g; #add commas for readability
printf STDERR
(
    "[$time] processed $i introns - %d U12 introns, %d U2 introns,"
    . " %d unknown intron types\nFinished\n",
    $intron_counts{U12},
    $intron_counts{U2},
    $intron_counts{unknown},
);
if ($TRANS){
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] writing transcript intron counts to $opts{t}...\n";
    writeTranscriptCounts();
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished - processing " . scalar(keys %transcripts) . " transcripts\n";
    close $TRANS;
}

#################################################
sub writeTranscriptCounts{
    foreach my $k (sort keys %transcripts){
        my ($unknown, $u2, $u12, $gene, $gene_type, $tsl, $length, $spliced_length) 
        =  (0, 0, 0, '.', '.', '.', 0, 0); 
        foreach my $type (keys %{$transcripts{$k}}){
            if ($type eq 'parent'){
                ($gene = $transcripts{$k}->{$type}) =~ s/^gene://;
            }elsif($type eq 'gene_type'){
                $gene_type = $transcripts{$k}->{$type};
            }elsif($type eq 'tsl'){
                $tsl = $transcripts{$k}->{$type};
            }elsif($type eq 'length'){
                $length = $transcripts{$k}->{$type};
            }elsif($type eq 'spliced_length'){
                $spliced_length = $transcripts{$k}->{$type};
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
            $gene_type,
            $tsl,
            $u12,
            $u2,
            $unknown,   
            $length,
            $spliced_length,
        ) . "\n";
    }
}
#################################################
sub parseExons{
    my $exons = shift;
    my $n = 0;
    foreach my $tr (keys %$exons){
        foreach my $ex (sort {$a <=> $b} keys %{$exons->{$tr}}){
            $transcripts{$tr}->{spliced_length} += 1 + $exons{$tr}->{$ex}->end - $exons{$tr}->{$ex}->start;
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
    printf STDERR
    (
        "[$time] processed $n introns - %d U12 introns, %d U2 introns,"
        . " %d unknown intron types\n",
        $intron_counts{U12},
        $intron_counts{U2},
        $intron_counts{unknown},
    );
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
        $intron_stop = $exon1->end + 1; 
        $intron_start = $exon2->start - 1;
    }
    my $intron = new Bio::SeqFeature::Generic
    (
        -seq_id      => $chrom,
        -start       => min($intron_start, $intron_stop),
        -end         => max($intron_start, $intron_stop),
        -strand      => $strand,
        -source_tag  => 'spliceScorer',
        -primary_tag => 'intron',
    
    );
    $intron->add_tag_value
    (
        'previous_exon_start',
        $exon1->start,
    );

    $intron->add_tag_value
    (
        'previous_exon_end',
        $exon1->end,
    );

    $intron->add_tag_value
    (
        'next_exon_start',
        $exon2->start,
    );

    $intron->add_tag_value
    (
        'next_exon_end',
        $exon2->end,
    );

    $intron->add_tag_value
    (
        'Parent',
        $exon1->get_tag_values('Parent'),
    );
    foreach my $t 
    (qw /
        exon_id
        exon_number
        rank
    /){
        if ($exon1->has_tag($t)){
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
    }
    my ($d_start, $d_end) = sort {$a <=> $b} 
    (
        ($intron_start - (3 * $strand) ) ,
        ($intron_start + (10 * $strand) ) ,
    );
    my ($a_start, $a_end) = sort {$a <=> $b} 
    (
        ($intron_stop - (100 * $strand)) ,
        ($intron_stop + (3 * $strand) ) ,
    );
    my $donor    = $fai->fetch("$chrom:$d_start-$d_end");
    my $acc_and_branch = $fai->fetch("$chrom:$a_start-$a_end");
    if ($strand < 0){
        $donor = reverse_complement($donor);
        $acc_and_branch = reverse_complement($acc_and_branch);
    }
    my $acceptor = substr($acc_and_branch, 100 - 13,); 
    my $branch = substr($acc_and_branch, 0, 92); 
    $intron->add_tag_value
    (
        'donor_seq',
        lc ( substr($donor, 0, 3) ) . 
        uc (substr($donor, 3, ) ),
    );
    $intron->add_tag_value
    (
        'acceptor_seq',
        uc ( substr($acceptor, 0, 14) ) . 
        lc (substr($acceptor, 14, ) ),
    );
    
    my %scores = ();
    my %branch_seqs = (); 
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
        ($scores{'B'}->{$type}, $branch_seqs{$type}) = 
        ScoreSpliceSite::scanForBranchPoint
        (
            seq  => $branch,
            type => $type,
            species => $opts{s},
        );
    }
    my $u12_b_score;
    my $u12_b_seq;
    if ($scores{'B'}->{AT_AC_U12} > $scores{'B'}->{GT_AG_U12}){
        $u12_b_score = $scores{'B'}->{AT_AC_U12};
        $u12_b_seq = uc ($branch_seqs{AT_AC_U12} );
    }else{
        $u12_b_score = $scores{'B'}->{GT_AG_U12};
        $u12_b_seq = uc ($branch_seqs{GT_AG_U12} );
    }
    $intron->add_tag_value
    (
        "U12_branch_score",
        $u12_b_score,
    );
    $intron->add_tag_value
    (
        "U12_branch_best_seq",
        $u12_b_seq,
    );
    
    $intron->add_tag_value
    (
        "U2_branch_score",
        sprintf("%.2f", $scores{'B'}->{GT_AG_U2}),
    );
    $intron->add_tag_value
    (
        "U2_branch_best_seq",
        $branch_seqs{GT_AG_U2},
    );
    my $best_u2;
    my $best_u12;
    foreach my $type (ScoreSpliceSite::getIntronTypes){
        $intron->add_tag_value
        (
            "donor_score_$type",
            sprintf("%.2f", $scores{'D'}->{$type}),
        );
        $intron->add_tag_value
        (
            "acceptor_score_$type",
            sprintf("%.2f", $scores{'A'}->{$type}),
        );

    }
    

    my $intron_type = ScoreSpliceSite::determineIntronType
    (
            donor    => $donor,
            acceptor => $acceptor,
            branch   => $branch,
            species  => $opts{s},
    ); 

    $intron->add_tag_value("intron_type", $intron_type); 
    $gffwriter->write_feature($intron);
    $transcripts{$tr}->{$intron_type}++;
    if ($intron_type){
        if ($intron_type =~ /U12$/){
            $intron_counts{U12}++;
        }else{
            $intron_counts{U2}++;
        }
    }else{
        $intron_counts{unknown}++;
    }
}

#################################################
sub usage{
    my $msg = shift;
    print "\n$msg\n" if $msg;

    print <<EOT

Create a GFF3 file of introns scored for U2 and U12 splice sites

USAGE: $0 -f genome_fasta.fa -g genes_and_exons.gff3

OPTIONS:
    
    -f,--fasta FILE
        Genome fasta file for retrieving DNA sequences for intron-exon boundaries

    -g,--gff FILE
        GFF3 file containing information on the genes and exons to use for intron file creation

    -s,--species INT
        Taxonomic code for species to use for splice prediction. Default is 9606 (human). 
        Available species are 10090 (mouse), 3702 (A. thaliana), 6239 (C. elegans), 7227 (D. melanogaster) and 9606 (human).

    -o,--output FILE
        Optional output file. Default = STDOUT.

    -t,--transcripts FILE
        Optional output file giving counts of intron types per transcript.

    -h,--help 
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}




