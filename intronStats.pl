#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw( min max sum );
use POSIX qw/strftime/;
use Bio::Tools::GFF;
use FindBin qw($RealBin);
use File::Temp qw/ tempfile /;
use Sort::External;
use Bio::DB::Sam;
use Bio::SeqFeature::Generic;
use lib "$RealBin/lib";
use ScoreSpliceSite;
use ReverseComplement qw/ reverse_complement /;

my %opts = 
(
    m => 3, 
    x => 999999999999999999999999,
    u => 1,
    y => 6,
);
GetOptions(
    \%opts,
    'f|fasta=s',
    'g|gff=s',
    'b|biotype=s',
    'u|min_repeat_unit_length=i',
    'y|max_repeat_unit_length=i',
    'm|min_repeat_length=i',
    'x|max_repeat_length=i',
    'o|output=s',
    'e|exon_seqs=s',
    'h|?|help',
) or usage("Error getting options!");
usage() if $opts{h};
checkOptions();

my $IN;
if ($opts{g} =~ /\.gz$/){
    open ($IN, "gzip -dc $opts{g} |") or die "Can't open $opts{g} via gzip: $!\n";
}else{
    open ($IN, "<", $opts{g}) or die "Can't open $opts{g} for reading: $!\n";
}

my $OUT = setupOutput();
my ($TMPEX, $ex_seq_file) = setupExonSeqFile();

my %exon_seqs = ();#all values undef, records ids of exons we've already retrieved seq for
my %u12 = ();#records exon IDs for exons with a neighbouring U12 intron
my %u2  = ();#as above for U2
my %unknown = ();#as above for unknown intron types
#my %names = ();

processIntronGff();

sortExonSeqs();

outputExonStats();


#################################################
sub outputExonStats{
    open (my $EX, '<', $opts{e}) or die 
     "Can't read exon sequence file '$opts{e}': $!\n";
    my $n = 0; 
    while (my $line = <$EX>){
        my ($id, $seq) = split("\t", $line); 
        my ($class, $subclass);
        if (exists $u12{$id}){
            $class = 'U12';
            $subclass = $u12{$id};
        }elsif (exists $unknown{$id}){
            $class = 'UNKNOWN';
            $subclass = $unknown{$id};
        }elsif (exists $u2{$id}){
            $class = 'U2';
            $subclass = $u2{$id};
        }else{
            #some exons (e.g. one-gene-exon exons) will not have an associated
            # intron, so no need to warn
            next;
        }
        writeExonStats($class, $subclass, $id, $seq);
        $n++;
        if (not $n % 10000){
            my $time = strftime( "%H:%M:%S", localtime );
            print STDERR "[$time] Wrote stats for $n exons...\n" ;
        }
    }
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished writing stats for $n exons.\n";
    close $OUT;
}
#################################################
sub sortExonSeqs{
    return if not $TMPEX;
    close $TMPEX;
    open (my $EXIN, '<', $ex_seq_file) or die 
     "Can't open temporary exon sequence file '$ex_seq_file' for reading: $!\n";
    my $sortex = Sort::External->new(mem_threshold => 1024**2 * 16);
    #exon IDs should all be same length so no need for special sort sub
    my @feeds = ();
    my $n = 0;
    my $SORTOUT;
    if ($opts{e}){
        open ($SORTOUT, '>', $opts{e}) or die 
         "Can't open exon sequence output file '$opts{e}' for writing: $!\n";
    }else{
        ($SORTOUT, $opts{e}) = tempfile("ex_seqs_XXXX", UNLINK => 1);
    }
    while (my $line = <$EXIN>){
        next if $line =~ /^#/;
        push @feeds, $line;
        $n++;
        if (@feeds > 9999){
            $sortex->feed(@feeds);
            @feeds = ();
            my $time = strftime( "%H:%M:%S", localtime );
            print STDERR "[$time] Fed $n exons to sort...\n";
        }
    }
    $sortex->feed(@feeds) if @feeds;
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished feeding $n exons to sort...\n";
    @feeds = ();
    $sortex->finish; 
    $n = 0;
    print $SORTOUT <<EOT
#sorted exon sequence file
#EXON_ID\tSEQ
EOT
;
    while ( defined( $_ = $sortex->fetch ) ) {
        print $SORTOUT $_;
    }
    close $SORTOUT;
}

#################################################
sub setupExonSeqFile{
    my $ex;
    my $TMP;
    if ($opts{e}){
    #if user supplied --exon_seqs file see it it already exists 
    # - will read from it instead of retrieving seqs if it does
    # we write to a tmp file even if we will be creating a new --exon_seqs
    # file cos we need to sort the output once we've got all the seqs
        $ex = $opts{e};
        if (not -e $ex){
            ($TMP, $ex) = tempfile("ex_seqs_XXXX", UNLINK => 1);
        }else{
            return (undef, undef);
        }
    }else{
        ($TMP, $ex) = tempfile("ex_seqs_XXXX", UNLINK => 1);
    }
    print $TMP <<EOT
#unsorted temporary exon sequence file
#EXON_ID\tSEQ
EOT
;
    return ($TMP, $ex);
}

#################################################
sub setupOutput{
    my $FH;
    if ($opts{o}){
        open ($FH, ">", $opts{o}) or die "Can't open $opts{o} for writing: $!\n";    
    }else{
        $FH = \*STDOUT;
    }
    print $FH join
    (   
        "\t",
        qw(
            INTRON_TYPE
            SUBTYPE
            EXON_ID
            EXON_LENGTH
            PERCENT_GC
            REPEATS
            TOTAL_REPEAT_LENGTH
            LONGEST_REPEAT
            MEAN_REPEAT_LENGTH
            HOMOPOLYMERS
            TOTAL_HOMOPOLYMER_LENGTH
            LONGEST_HOMOPOLYMER
            MEAN_HOMOPOLYMER_LENGTH
        )
    ) . "\n";
    return $FH;
}

#################################################
sub processIntronGff{
    my $fai = Bio::DB::Sam::Fai->load($opts{f});#should create index if it doesn't exist
    my $gff = Bio::Tools::GFF->new
    (
        -gff_version => 3,
        -fh          => $IN,
    );
    my $intron_count = 0;
    my $biotype = ''; 
    my @introns = (); #collect all intron features per transcript for processing
    while (my $feat = $gff->next_feature() ) {
        if($feat->has_tag('gene_id')){
            #parse previous introns
            parseIntrons(\@introns, \$intron_count);
            #get name and gene id
            my ($id) = $feat->get_tag_values('ID'); 
            $id =~ s/^gene://;
            my $name = '.';
    #        if ($feat->has_tag('Name')){
    #            ($name) = $feat->get_tag_values('Name'); 
    #        }
    #        $names{$id} = $name; 
        }elsif ($feat->has_tag('transcript_id')){
            #parse introns in case we missed a gene_id tag
            parseIntrons(\@introns, \$intron_count);
            #get transcript name and associate with gene id
    #        my ($tr) = $feat->get_tag_values('transcript_id'); 
    #        my ($parent) = $feat->get_tag_values('Parent');
            $biotype = join(",", $feat->get_tag_values('biotype'));
        }elsif ($feat->primary_tag eq 'exon'){
            #collect exons 
            if ($opts{b}){
                if (grep {$_ eq $opts{b}} split(",", $biotype)){
                    writeTempExonSequence($feat, $fai);
                }
            }else{
                writeTempExonSequence($feat, $fai);
            }
        }elsif ($feat->primary_tag eq 'intron'){
            #collect intron type related to each exon
            if ($opts{b}){
                if (grep {$_ eq $opts{b}} split(",", $biotype)){
                    push @introns, $feat;
                }
            }else{
                push @introns, $feat;
            }
        }
    }
    parseIntrons(\@introns, \$intron_count);
    $gff->close();
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Done reading input - processed $intron_count introns.\n";
}

#################################################
sub parseIntrons{
    #if an exon has a neighbouring U12 intron classify as U12
    #if an exon has no neighbouring U12 introns and its neighbouring
    # introns are confirmed U2 classify as U2
    #if an exon has a neighbouring 'unknown' intron and does not have 
    # a neighbouring U12 intron classify as UNKNOWN
    my $in = shift;
    my $count = shift;
    my %all = ();
    foreach my $intr (@$in){
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
        $$count++;
        reportProgress($$count);
    }
    @$in = (); 
}


#################################################
sub writeExonStats{
    my ($class, $subclass, $exon, $seq) = @_;
    my $gc = getGcPercentage($seq);
    my @r = getRepeats($seq);
    my @r_lengths;
    my @h_lengths;
    my ($longest, $mean, $total) = (0, 0, 0);
    my ($h_longest, $h_mean, $h_total) = (0, 0, 0);
    foreach my $rep (@r){
        my $l = length($rep);
        push @r_lengths, $l;
        if ($rep =~ /^(\w)(\1+)+$/){
            push @h_lengths, $l;
        }
    }
    if (@r_lengths){
        $longest = max(@r_lengths);
        $total = sum(@r_lengths);
        $mean = $total/@r_lengths;
    }
    if (@h_lengths){
        $h_longest = max(@h_lengths);
        $h_total = sum(@h_lengths);
        $h_mean = $total/@h_lengths;
    }
    print $OUT join
    (
        "\t",
        $class,
        $subclass,
        $exon,
        length($seq),
        sprintf("%.3f", $gc),
        scalar(@r),
        $total,
        $longest,
        $mean,
        scalar(@h_lengths),
        $h_total,
        $h_longest,
        $h_mean,
    ) . "\n";
    #intron type (U2 or U12 or UNKNOWN), subtype, exon ID, exon length, % GC,
    # no. repeats, total length of repeats, longest repeat, mean repeat length,
    # no.homopolymers, total length of homopolymers, longest homopolymer, mean homopolymer length
}

#################################################
sub getRepeats{
#return array of 1-6 (or rather $opts{u}-$opts{y}) nucleotide repeats at least 3 nt long found in $seq
    my $seq = shift;
    my @h = ();
    while ($seq =~ /(\w{$opts{u},$opts{y}})(\1)(\1+)*/g){
        my $rep = "$1$2";
        $rep .= $3 if $3;
        push @h, $rep if length($rep) >= $opts{m} and length($rep) <= $opts{x};
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
sub writeTempExonSequence{
    return if not $TMPEX;
    my $exon = shift;
    my $fai = shift;
    my ($id) = $exon->get_tag_values('exon_id');
    return if exists $exon_seqs{$id};
    my $chrom = $exon->seq_id; 
    my $strand = $exon->strand;
    my $start = $exon->start;
    my $end = $exon->end;
    my $seq = $fai->fetch("$chrom:$start-$end");
    if ($strand < 0){
        $seq = reverse_complement($seq);
    }
    print $TMPEX "$id\t$seq\n";
    $exon_seqs{$id} = undef;
}

#################################################

sub getExonSequence{
    my $exon = shift;
    my $fai = shift;
    my ($id) = $exon->get_tag_values('exon_id');
    return if exists $exon_seqs{$id};
    my $chrom = $exon->seq_id; 
    my $strand = $exon->strand;
    my $start = $exon->start;
    my $end = $exon->end;
    my $seq = $fai->fetch("$chrom:$start-$end");
    if ($strand < 0){
        $seq = reverse_complement($seq);
    }
    $exon_seqs{$id} = $seq;
}

#################################################
sub checkOptions{
    usage("-f/--fasta argument is required.") if not $opts{f};
    usage("-g/--gff argument is required.") if not $opts{g};
    my %lengths = 
    (
        u => "min_repeat_unit_length",
        y => "max_repeat_unit_length",
        m => "min_repeat_length",
        x => "max_repeat_length",
    );
    foreach my $k (keys %lengths){
        usage("-$k/--$lengths{$k} argument must be greater than 0.") if $opts{$k} < 1;
    }
    if ($opts{m} > $opts{x}){
        usage("-m/--min_repeat_length argument ($opts{m}) is greater than "
              ."-x/--max_repeat_length argument ($opts{x}).") ;
    }
    if ($opts{u} > $opts{m}){
        usage("-u/--min_repeat_unit_length ($opts{u}) is greater than "
              ."-m/--min_repeat_length argument argument ($opts{m}).") ;
    }
    if ($opts{u} > $opts{y}){
        usage("-u/--min_repeat_unit_length argument ($opts{u}) is greater than "
              ."-y/--max_repeat_unit_length argument ($opts{y}).") ;
    }
}

#################################################
sub reportProgress{
    my $n = shift;
    return if not $n;
    return if ($n % 10000); 
    my $time = strftime( "%H:%M:%S", localtime );
    $n =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g; #add commas for readability
    print STDERR "[$time] processed $n introns...\n";
}

#################################################
sub usage{
    my $msg = shift;
    print STDERR "\n$msg\n" if $msg;

    print STDERR <<EOT

TODO 

USAGE: $0 -f genome_fasta.fa -g introns.gff3

OPTIONS:
    
    -f,--fasta FILE
        Genome fasta file for retrieving DNA sequences for intron-exon boundaries

    -g,--gff FILE
        GFF3 intron file created by spliceScorer.pl

    -b,--biotypes STRING
        Only include exons from transcripts of this biotype (e.g. protein_coding)

    -m,--min_repeat_length INT
        Only count repeats at least this long (default = 3)
    
    -x,--max_repeat_length INT
        Only count repeats of this value or shorter (default = 999999999999999999999999)
    
    -u,--min_repeat_unit_length INT
        Only count repetetive sequences where the repeated unit is at least this long (default = 1)
    
    -y,--max_repeat_unit_length INT
        Only count repetetive sequences where the repeated unit is this long or shorter (default = 6)
    
    -o,--output FILE
        Optional output file. Default = STDOUT.
    
    -e,--exon_seqs FILE
        Optional file for reading/writing DNA sequences of exons retrieved. 
        If this file exists it will be read for assessing exon stats; if not 
        it will be created for future use.

    -h,--help 
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}


