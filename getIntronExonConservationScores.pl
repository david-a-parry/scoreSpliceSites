#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use POSIX qw/strftime/;
use File::Temp qw/ tempfile /;
use Sort::External;
use Bio::Tools::GFF;
use FindBin qw($RealBin);
use Bio::DB::BigWig 'binMean';

my %opts = ();
GetOptions(
    \%opts,
    'g|gff=s',
    'b|biotype=s',
    's|tsl=i',
    'j|merge',
    'o|output=s',
    'w|bigwig=s',
    'h|?|help',
) or usage("Error getting options!");
usage() if $opts{h};
usage("-g/--gff option is required!\n") if not $opts{g};
usage("-w/--bigwig option is required!\n") if not $opts{w};
usage("-o/--output argument is required.") if not $opts{o};

my $IN;
if ($opts{g} =~ /\.gz$/){
    open ($IN, "gzip -dc $opts{g} |") or die "Can't open $opts{g} via gzip: $!\n";
}else{
    open ($IN, "<", $opts{g}) or die "Can't open $opts{g} for reading: $!\n";
}
my $wig  = Bio::DB::BigWig->new
(
    -bigwig => $opts{w},
);

my %exons = ();
my %introns = ();
my %u12_introns = ();#records U12 intron IDs 
my %u2_introns = ();#records U2 intron IDs 
my %unknown_introns = ();#records IDs of introns of unknown type
my %u12_exons = ();#records exon IDs for exons with a neighbouring U12 intron
my %u2_exons  = ();#as above for U2
my %unknown_exons = ();#as above for unknown intron types

my ($I_OUT, $E_OUT) = setupOutput();
my ($TMPIN, $tmp_intron) = tempfile( UNLINK => 1);
my ($TMPEX, $tmp_exon)   = tempfile( UNLINK => 1);

processGff();

close $TMPIN;
close $TMPEX;

$tmp_exon   = sortRegions($tmp_exon, "exon");
outputExonStats();
$tmp_intron = sortRegions($tmp_intron, "intron",);
outputIntronStats();

#################################################
sub setupOutput{
    my $ex = "$opts{o}_exon_cons_score.tsv";
    my $in = "$opts{o}_intron_cons_score.tsv";
    open (my $EX, ">", "$ex") or die "Can't open $ex for writing: $!\n";    
    open (my $IN, ">", "$in") or die "Can't open $in for writing: $!\n";    
    print $EX join
    (
        "\t", 
        qw /
            Region
            Strand
            ID
            Class
            Subclass
            MeanScore
            WithoutFlanking8Score
            Length
        /
    ) . "\n";
    print $IN join
    (
        "\t", 
        qw /
            Region
            Strand
            ID
            Class
            Subclass
            MeanScore
            DonorConsensusScore
            AcceptorConsensusScore
            DonorFlanking50Score
            AcceptorFlanking50Score
            WithoutFlanking50Score
            Length
        /
    ) . "\n";
    return ($IN, $EX);
}

#################################################
sub outputIntronStats{
    open (my $IN, '<', $tmp_intron) or die 
     "Can't read intron sequence file '$tmp_intron': $!\n";
    my $n = 0;
    my %prev_intron;
    while (my $line = <$IN>){
        chomp $line;
        my ($chrom, $start, $end, $id, undef, $strand) = split("\t", $line); 
        my $coord = "$chrom:$start-$end";
        my ($class, $subclass);
        if (exists $u12_introns{$coord}){
            $class = 'U12';
            $subclass = $u12_introns{$coord};
        }elsif (exists $unknown_introns{$coord}){
            $class = 'UNKNOWN';
            $subclass = $unknown_introns{$coord};
        }elsif (exists $u2_introns{$coord}){
            $class = 'U2';
            $subclass = $u2_introns{$coord};
        }else{
            #some exons (e.g. one-gene-exon exons) will not have an associated
            # intron, so no need to warn
            #TODO - check whether we should warn for some exons? (i.e. if not a single exon gene)
            next;
        }
        my %this_intron = 
        (
            class    => $class,
            subclass => $subclass,
            chrom    => $chrom,
            id       => $id,
            start    => $start,
            end      => $end,
            strand   => $strand,
        );
        if (not %prev_intron){
            %prev_intron = %this_intron;
        }
        my $written = checkAndWrite(\%prev_intron, \%this_intron, $I_OUT);
        if($written){
            $n++;
            if (not $n % 10000){
                my $time = strftime( "%H:%M:%S", localtime );
                print STDERR "[$time] Processed scores for $n introns...\n" ;
            }
        }
    }
    writeScores
    (
        class    => $prev_intron{class},
        subclass => $prev_intron{subclass},
        chrom    => $prev_intron{chrom},
        start    => $prev_intron{start},
        end      => $prev_intron{end},
        strand   => $prev_intron{strand},
        id       => $prev_intron{id},
        fh       => $I_OUT,
    ) if (%prev_intron);
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished writing scores for $n introns.\n";
    close $IN;
    close $I_OUT;
}

#################################################
sub checkAndWrite{
#read our hashes of intron features
#if they overlap and are of same type, merge
#otherwise print previous and store current 
    my $prev    = shift;
    my $current = shift;
    my $FH = shift;
    if ($opts{j} and 
        $prev->{chrom} eq $current->{chrom} and 
        $prev->{end} >= $current->{start} and #should be sorted in coordinate order
        $prev->{strand} eq $current->{strand} and #should be sorted in coordinate order
        $prev->{subclass} eq $current->{subclass}
    ){#overlaps
        if ($current->{end} > $prev->{end}){
            $prev->{id} = $prev->{id} . "/" . $current->{id};
            $prev->{end} = $current->{end};
        }
        return 0;
    }else{#doesn't overlap
        writeScores
        (
            class    => $prev->{class},
            subclass => $prev->{subclass},
            chrom    => $prev->{chrom},
            start    => $prev->{start},
            end      => $prev->{end},
            strand   => $prev->{strand},
            id       => $prev->{id},
            fh       => $FH,
        );
        %{$prev} = %{$current};
        return 1;
    }
}

#################################################
sub outputExonStats{
    open (my $EX, '<', $tmp_exon) or die 
     "Can't read exon sequence file '$tmp_exon': $!\n";
    my $n = 0; 
    my %prev_exon;
    while (my $line = <$EX>){
        chomp $line;
        my ($chrom, $start, $end, $id, undef, $strand) = split("\t", $line); 
        my $coord = "$chrom:$start-$end";
        my ($class, $subclass);
        if (exists $u12_exons{$coord}){
            $class = 'U12';
            $subclass = $u12_exons{$coord};
        }elsif (exists $unknown_exons{$coord}){
            $class = 'UNKNOWN';
            $subclass = $unknown_exons{$coord};
        }elsif (exists $u2_exons{$coord}){
            $class = 'U2';
            $subclass = $u2_exons{$coord};
        }else{
            #some exons (e.g. one-gene-exon exons) will not have an associated
            # intron, so no need to warn
            #TODO - check whether we should warn for some exons? (i.e. if not a single exon gene)
            next;
        }
        my %this_exon = 
        (
            class    => $class,
            subclass => $subclass,
            chrom    => $chrom,
            id       => $id,
            start    => $start,
            end      => $end,
            strand   => $strand,
        );
        if (not %prev_exon){
            %prev_exon = %this_exon;
            next;
        }
        my $written = checkAndWrite(\%prev_exon, \%this_exon, $E_OUT);
        if($written){
            $n++;
            if (not $n % 10000){
                my $time = strftime( "%H:%M:%S", localtime );
                print STDERR "[$time] Wrote scores for $n exons...\n" ;
            }
        }
    }
    writeScores
    (
        class    => $prev_exon{class},
        subclass => $prev_exon{subclass},
        chrom    => $prev_exon{chrom},
        start    => $prev_exon{start},
        end      => $prev_exon{end},
        strand   => $prev_exon{strand},
        id       => $prev_exon{id},
        fh       => $E_OUT,
    ) if %prev_exon;
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished writing scores for $n exons.\n";
    close $EX;
    close $E_OUT;
}


#################################################
sub sortRegions{
    my ($tmp_seq_file, $type) = @_;
    open (my $REGIN, '<', $tmp_seq_file) or die 
     "Can't open temporary $type regions file '$tmp_seq_file' for reading: $!\n";
    my $sortex = Sort::External->new(mem_threshold => 1024**2 * 16);
    #exon IDs should all be same length so no need for special sort sub
    my @feeds = ();
    my $n = 0;
    my $time;
    my ($SORTOUT, $sortfile) = tempfile(UNLINK => 1);
    while (my $line = <$REGIN>){
        next if $line =~ /^#/;
        my ($chrom, $start, $end, $id, undef, $strand) = split("\t", $line); 
        my $s_chrom = sprintf("%-25s", $chrom); 
        my $packstart = pack("N", $start); 
        my $packend = pack("N", $end); 
        push @feeds, "$s_chrom,$packstart,$packend,$line";
        $n++;
        if (@feeds > 9999){
            $sortex->feed(@feeds);
            @feeds = ();
            $time = strftime( "%H:%M:%S", localtime );
            print STDERR "[$time] Fed $n $type"."s to sort...\n";
        }
    }
    close $REGIN;
    $sortex->feed(@feeds) if @feeds;
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished feeding $n $type"."s to sort...\n";
    @feeds = ();
    $sortex->finish; 
    $n = 0;
    while ( defined( $_ = $sortex->fetch ) ) {
        print $SORTOUT substr($_, 36,);
        $n++;
        if (not $n % 100000){
            $time = strftime( "%H:%M:%S", localtime );
            print STDERR "[$time] Wrote $n $type"."s to regions file...\n";
        }
    }
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished writing $n $type regions to regions file...\n";
    close $SORTOUT;
    return $sortfile;
}

#################################################
sub processGff{
    my $gff = Bio::Tools::GFF->new
    (
        -gff_version => 3,
        -fh          => $IN,
    );
    my $intron_count = 0;
    my $gene_type = ''; 
    my $tsl = '';
    my @introns = (); #collect all intron features per transcript for processing
    while (my $feat = $gff->next_feature() ) {
        if ($feat->primary_tag eq 'gene'){
            #parse previous introns
            collectIntrons(\@introns, \$intron_count);
            #get name and gene id
            my ($id) = $feat->get_tag_values('ID'); 
            $id =~ s/^gene://;
            my $name = '.';
    #        if ($feat->has_tag('Name')){
    #            ($name) = $feat->get_tag_values('Name'); 
    #        }
    #        $names{$id} = $name; 
        }elsif ($feat->primary_tag eq 'transcript'){
            #parse introns in case we missed a gene_id tag
            collectIntrons(\@introns, \$intron_count);
            #get transcript name and associate with gene id
    #        my ($tr) = $feat->get_tag_values('transcript_id'); 
    #        my ($parent) = $feat->get_tag_values('Parent');
            $gene_type = join(",", $feat->get_tag_values('gene_type'));
            eval
            {
                ($tsl) = $feat->get_tag_values('transcript_support_level'); 
            };#not all transcripts have tsl tag - e.g. pseudogenes
            $tsl ||= 'NA';
            #clear cruft from tsl
            $tsl =~ s/\s+.*//;
        }elsif ($feat->primary_tag eq 'exon'){
            #collect exons 
            my $do_write = 1;
            if ($opts{b}){
                if (not grep {$_ eq $opts{b}} split(",", $gene_type)){
                    $do_write = 0;
                }
            }
            if ($opts{s}){
                if ($tsl eq 'NA'){
                    $do_write = 0;
                }elsif($tsl > $opts{s}){
                    $do_write = 0;
                }
            }
            if ($do_write){
                writeTempExonRegion($feat);
            }
        }elsif ($feat->primary_tag eq 'intron'){
            #collect intron type related to each exon
            my $do_write = 1;
            if ($opts{b}){
                if (not grep {$_ eq $opts{b}} split(",", $gene_type)){
                    $do_write = 0;
                }
            }
            if ($opts{s}){
                if ($tsl eq 'NA'){
                    $do_write = 0;
                }elsif($tsl > $opts{s}){
                    $do_write = 0;
                }
            }
            if ($do_write){
                push @introns, $feat;
            }
        }
    }
    collectIntrons(\@introns, \$intron_count);
    $gff->close();
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Done reading input - found $intron_count introns.\n";
}

#################################################
sub writeTempExonRegion{
    my $exon   = shift;

    my ($id)   = $exon->get_tag_values('exon_id');
    my $chrom  = $exon->seq_id; 
    my $strand = $exon->strand;
    my $start  = $exon->start;
    my $end    = $exon->end;
    my $coord  = "$chrom:$start-$end";
    return if exists $exons{$coord};
    print $TMPEX join
    (
        "\t",
        $chrom,
        $start,
        $end,
        $id,
        0,
        $strand,
    ) . "\n";
    $exons{$coord} = undef;
}

#################################################
sub writeTempIntronRegion{
    my $intron  = shift;
    my $next_ex = shift;
    my $prev_ex = shift;

    my $id      = "$prev_ex-$next_ex";
    my $chrom   = $intron->seq_id; 
    my $strand  = $intron->strand;
    my $start   = $intron->start;
    my $end     = $intron->end;
    my $coord   = "$chrom:$start-$end";
    return if exists $introns{$coord};
    print $TMPIN join
    (
        "\t",
        $chrom,
        $start,
        $end,
        $id,
        0,
        $strand,
    ) . "\n";
    $introns{$coord} = undef;
}   

#################################################
sub writeScores{
    my %args = @_;
#            class    => $prev->{class},
#            subclass => $prev->{subclass},
#            chrom    => $prev->{chrom},
#            start    => $prev->{start}
#            end      => $prev->{end},
#            strand   => $prev->{strand},
#            id       => $prev->{id},
#            fh       => $FH,
    my $coord = $args{chrom} . ":" . $args{start} . "-" . $args{end};
    my $gmean = getScoreMean
    (
        -seq_id => $args{chrom},
        -start  => $args{start},
        -end    => $args{end},
    );
    if ($args{fh} eq $I_OUT){#intron 
        #get scores for splice consensus sequences
        my (@donor, @acceptor, @acceptor50, @donor50); 
        push @donor, $args{strand} > 0 ? $args{start} - 3  : $args{end} + 3;
        push @donor, $args{strand} > 0 ? $args{start} + 10 : $args{end} - 10;
        push @acceptor, $args{strand} > 0 ? $args{end} - 13 : $args{start} + 13;
        push @acceptor, $args{strand} > 0 ? $args{end} + 3  : $args{start} - 3;
        push @acceptor50, $args{strand} > 0 ? $args{end} - 50 : $args{start} + 50;
        push @acceptor50, $args{strand} > 0 ? $args{end} : $args{start} ;
        push @donor50, $args{strand} > 0 ? $args{start} + 50 : $args{end} - 50;
        push @donor50, $args{strand} > 0 ? $args{start} : $args{end} - 50;
        @donor      = sort {$a <=> $b } @donor;
        @acceptor   = sort {$a <=> $b } @acceptor;
        @donor50    = sort {$a <=> $b } @donor50;
        @acceptor50 = sort {$a <=> $b } @acceptor50;
        my $donor_mean = getScoreMean
        (
            -seq_id => $args{chrom},
            -start  => $donor[0],
            -end    => $donor[1],
        );
        
        my $acceptor_mean = getScoreMean
        (
            -seq_id => $args{chrom},
            -start  => $acceptor[0],
            -end    => $acceptor[1],
        );
        
        my $donor50_mean = getScoreMean
        (
            -seq_id => $args{chrom},
            -start  => $donor50[0],
            -end    => $donor50[1],
        );
        
        my $acceptor50_mean = getScoreMean
        (
            -seq_id => $args{chrom},
            -start  => $acceptor50[0],
            -end    => $acceptor50[1],
        );

        my $non_flanking_mean = getScoreMean
        (
            -seq_id => $args{chrom},
            -start  => $args{start} + 50,
            -end    => $args{end} - 50,
        );

        $args{fh}->print
        (
            join
            (
                "\t", 
                $coord,
                $args{strand},
                $args{id},
                $args{class},
                $args{subclass},
                $gmean,
                $donor_mean,
                $acceptor_mean,
                $donor50_mean,
                $acceptor50_mean,
                $non_flanking_mean,
                $args{end} - $args{start} + 1,
            ) . "\n"
        );
    }elsif($args{fh} eq $E_OUT){#exon
        my $non_flanking_mean = getScoreMean
        (
            -seq_id => $chr,
            -start  => $args{start} + 8,
            -end    => $args{end} - 8,
        );

        $args{fh}->print 
        (join
            (
                "\t", 
                $coord,
                $args{strand},
                $args{id},
                $args{class},
                $args{subclass},
                $gmean,
                $non_flanking_mean,
                $args{end} - $args{start} + 1,
            ) . "\n"
        );
    }else{
        die "Error determing filehandle - exiting ";
    }
}


#################################################
sub collectIntrons{
    #if an exon has a neighbouring U12 intron classify as U12
    #if an exon has no neighbouring U12 introns and its neighbouring
    # introns are confirmed U2 classify as U2
    #if an exon has a neighbouring 'unknown' intron and does not have 
    # a neighbouring U12 intron classify as UNKNOWN
    my $in = shift;
    my $count = shift;
    my %all = ();
    foreach my $intr (@$in){
        my ($previous_start) = $intr->get_tag_values('previous_exon_start');
        my ($previous_end)   = $intr->get_tag_values('previous_exon_end');
        my ($next_start)     = $intr->get_tag_values('next_exon_start');
        my ($next_end)       = $intr->get_tag_values('next_exon_end');
        my $previous         = $intr->seq_id . ":" . $previous_start . "-" . $previous_end;
        my $next             = $intr->seq_id . ":" . $next_start . "-" . $next_end;
        my ($previous_id)    = $intr->get_tag_values('previous_exon_id');
        my ($next_id)        = $intr->get_tag_values('next_exon_id');
        my ($type)           = $intr->get_tag_values('intron_type'); 
        my $intron_coords    = $intr->seq_id . ":" . $intr->start . "-" . $intr->end;
        if ($type =~ /U12$/){
            $u12_exons{$next}     = $type;
            $u12_exons{$previous} = $type;
            $u12_introns{$intron_coords} = $type;
        }elsif($type =~ /U2$/){
            $u2_exons{$next}     = $type;
            $u2_exons{$previous} = $type;
            $u2_introns{$intron_coords} = $type;
        }else{
            $unknown_exons{$next}     = $type;
            $unknown_exons{$previous} = $type;
            $unknown_introns{$intron_coords} = $type;
        }
        writeTempIntronRegion($intr, $next_id, $previous_id);
        $$count++;
        reportProgress($$count);
    }
    @$in = (); 
}
#################################################
sub getScoreMean{
    my %args = @_;
    my @score = $wig->features(%args, -type => 'summary');
    my $mean = 'NA';
    eval{
        $mean =  binMean( $score[0]->statistical_summary->[0] );
    };
    return $mean;
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


DESCRIPTION: Prints scores for introns and exons from a GFF file generated by spliceScorer.pl

USAGE: $0 -g introns.gff3 -w genomeBigwigFile.bw -o output_prefix

OPTIONS:
    
    -g,--gff FILE
        GFF3 intron file created by spliceScorer.pl

    -w,--bigwig
        Genome bigwig file containing GERP, phylop or phastcons scores

    -o,--output FILE
        Output file prefix. Two output files will be produced - one with the 
        suffix 'intron_stats.tsv' the other with the suffix 'exon_stats.tsv'.
    
    -b,--biotypes STRING
        Only include exons/introns from transcripts of this biotype (e.g. 
        protein_coding)
   
    -j,--merge
        Use this flag to merge overlapping introns/exons before calculating
        stats.

    -s,--tsl INT
        Only include exons/introns from transcripts with this transcript 
        support level or higher. Transcripts with a tsl of 1 have the most 
        support, those with a tsl of 5 have the least. See Ensembl's help:
        http://www.ensembl.org/Help/Glossary?id=492.

    -t,--trim INT
        Trim the 5' and 3' ends of exons by this amount. The default value of 3
        is used to remove bias from splice site consensus sequences used to 
        determine the intron types.

    -h,--help 
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}
