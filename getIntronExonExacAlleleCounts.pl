#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use POSIX qw/strftime/;
use File::Temp qw/ tempfile /;
use List::Util qw /sum /;
use Sort::External;
use Bio::Tools::GFF;
use FindBin qw($RealBin);
use Bio::DB::HTS::Tabix;
use lib "$RealBin/lib/dapPerlGenomicLib";
use VcfReader;

my %opts = ();
GetOptions(
    \%opts,
    'g|gff=s',
    'e|exac_vcf=s',
    'c|coverage_dir=s',
    'b|biotype=s',
    's|tsl=i',
    'j|merge',
    'o|output=s',
    'h|?|help',
) or usage("Error getting options!");
usage() if $opts{h};
usage("-g/--gff option is required!\n") if not $opts{g};
usage("-e/--exac_vcf option is required!\n") if not $opts{e};
usage("-c/--coverage_dir option is required!\n") if not $opts{c};
usage("-o/--output argument is required.") if not $opts{o};

my $IN;
if ($opts{g} =~ /\.gz$/){
    open ($IN, "gzip -dc $opts{g} |") or die "Can't open $opts{g} via gzip: $!\n";
}else{
    open ($IN, "<", $opts{g}) or die "Can't open $opts{g} for reading: $!\n";
}

my %search_args = VcfReader::getSearchArguments( $opts{e});
my %cov_iters = ();

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
    my $ex = "$opts{o}_exon_exac_freqs.tsv";
    my $in = "$opts{o}_intron_exac_freqs.tsv";
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
            AN
            AC
            AF
            MeanCoverage
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
            DonorConsensusAN
            DonorConsensusAC
            DonorConsensusAF
            DonorConsensusMeanCoverage
            AcceptorConsensusAN
            AcceptorConsensusAC
            AcceptorConsensusAF
            AcceptorConsensusMeanCoverage
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
    if ($args{fh} eq $I_OUT){#intron 
        #get scores for splice consensus sequences
        my (@donor, @acceptor);
        push @donor, $args{strand} > 0 ? $args{start} - 3  : $args{end} + 3;
        push @donor, $args{strand} > 0 ? $args{start} + 10 : $args{end} - 10;
        push @acceptor, $args{strand} > 0 ? $args{end} - 13 : $args{start} + 13;
        push @acceptor, $args{strand} > 0 ? $args{end} + 3  : $args{start} - 3;
        @donor      = sort {$a <=> $b } @donor;
        @acceptor   = sort {$a <=> $b } @acceptor;
        my %donor_af = getAfs
        (
            chrom => $args{chrom},
            start => $donor[0],
            end   => $donor[1],
        );
        
        my %acceptor_af = getAfs
        (
            chrom => $args{chrom},
            start => $acceptor[0],
            end   => $acceptor[1],
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
                $donor_af{an},
                $donor_af{ac},
                $donor_af{af},
                $donor_af{coverage},
                $acceptor_af{an},
                $acceptor_af{ac},
                $acceptor_af{af},
                $acceptor_af{coverage},
                $args{end} - $args{start} + 1,
            ) . "\n"
        );
    }elsif($args{fh} eq $E_OUT){#exon
        my %exon_af = getAfs
        (
            chrom => $args{chrom},
            start => $args{start},
            end   => $args{end},
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
                $exon_af{an},
                $exon_af{ac},
                $exon_af{af},
                $exon_af{coverage},
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
sub getAfs{
    my %args = @_;
    my %results = ();
    $args{chrom} =~ s/^chr//;#ExAC on v37 without chr 
    my @hits = VcfReader::searchByRegion
    (
        %search_args,
        %args,
    );
    foreach my $h (@hits){ 
        my @splt = split("\t", $h);
        my $an = VcfReader::getVariantInfoField(\@splt, "AN");
        my $ac = VcfReader::getVariantInfoField(\@splt, "AC");
        my @acs = split(",", $ac);
        $results{an} += $an;
        $results{ac} += sum(@acs);
    }
    if ($results{an}){
        $results{af} = $results{ac}/$results{an};
    }
    foreach my $k (qw /an ac af /){
        if (not defined $results{$k}){
            $results{$k} = 0;
        }
    }
    $results{coverage} = getCoverage(%args);
    return %results;
}

#################################################
sub getCoverage{
    my %args = @_;
    $args{chrom} =~ s/^chr//;#ExAC on v37 without chr 
    my $cov_file = "$opts{c}/Panel.chr$args{chrom}.coverage.txt.gz";
    if (not -e $cov_file){
        die "Could not find coverage file for $args{chrom} - '$cov_file'\n";
    }
    if (not exists $cov_iters{$cov_file}){
        $cov_iters{$cov_file} = Bio::DB::HTS::Tabix->new(filename =>  $cov_file);
    }
    my $it = $cov_iters{$cov_file}->query
    (
        "$args{chrom}:$args{start}-" . (1 + $args{end})
    );
    my $total_mean = 0;
    my $n = 0;
    while (my $result = $it->next){
        $n++;
        chomp($result);
        my (undef, undef, $mean) = split("\t", $result);
        $total_mean += $mean;
    }
    return $n ? $total_mean/$n : 0;
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

USAGE: $0 -g introns.gff3 -e exac_vcf -c exac_coverage/ -o output_prefix

OPTIONS:
    
    -g,--gff FILE
        GFF3 intron file created by spliceScorer.pl

    -e,--exac_vcf
        VCF of ExAC variants.

    -c,--coverage_dir
        Directory of coverage data from ExAC (as downloaded from 
        ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/coverage).

    -o,--output FILE
        Output file prefix. Two output files will be produced - one with the 
        suffix 'intron_exac_freqs.tsv' the other with the suffix 'exon_exac_freqs.tsv'.
    
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
