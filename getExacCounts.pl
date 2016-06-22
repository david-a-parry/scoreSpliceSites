#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use POSIX qw(strftime);
use Sort::External;
use FindBin qw($RealBin);
use List::Util qw(sum);
use Term::ProgressBar;
use Bio::DB::HTS;
use lib "$RealBin/lib/dapPerlGenomicLib/";
use VcfReader;

my %cov_iters = ();
my @an_fields = ();
my @ac_fields = ();
my @min_an    = ();#min allele number for a call to be counted
my %opts = 
(
    t  => 20,
    p  => 0.9,
    ac => \@ac_fields,
    an => \@an_fields,
    m  => \@min_an,
);
GetOptions(
    \%opts,
    'i|intron_bed=s',
    'e|exon_bed=s',
    'o|output=s',
    'v|exac_vcf=s',
    'c|exac_coverage_dir=s',
    't|coverage_threshold=i',
    'p|proportion_at_threshold=f',
    'ac=s{,}',
    'an=s{,}',
    'm|min_an=i{,}',
    'h|?|help',
) or usage("Error getting options!");
usage() if $opts{h};
usage("ERROR: -i/--intron_bed or -e/--exon_bed option is required!\n") 
      if not $opts{i} and not $opts{e};
usage("ERROR: -v/--exac_vcf option is required!\n") if not $opts{v};
usage("ERROR: -o/--output option is required!\n") if not $opts{o};
my %cov_columns = checkCoverageArgs();

if (not @ac_fields){
    @ac_fields = ("AC");
}
if (not @an_fields){
    @an_fields = ("AN");
}

if (@min_an){
    die "ERROR: --ac and --an options must be lists of the same length when used ".
        "with -m/--min_an option!\n"
        if @ac_fields != @an_fields;
    if (@min_an > 1 and @min_an != @ac_fields){
        die "ERROR: Please provide either only one --min_an threshold or provide one".
            " for each --ac field.\n";
    }
}

checkVcf();
my %search_args = VcfReader::getSearchArguments( $opts{v} );


processIntronBed() if $opts{i};
processExonBed() if $opts{e};


#################################################
sub processExonBed{
    print STDERR "Processing exon bed file $opts{e}...\n";
    my $EXON_BED;
    if ($opts{e} =~ /\.gz$/){
        open ($EXON_BED, "gzip -dc $opts{e} |") or die "Can't open $opts{e} via gzip: $!\n";
    }else{
        open ($EXON_BED, "<", $opts{e}) or die "Can't open $opts{e} for reading: $!\n";
    }
    my $line_count = 0;
    $line_count += tr/\n/\n/ while sysread($EXON_BED, $_, 2 ** 20);
    print STDERR "$opts{e} has $line_count lines...\n";
    my $progressbar = Term::ProgressBar->new
    (
        { 
          name => "Counting",
          count => ($line_count), 
          ETA => "linear" 
        } 
    );
    my $next_update      = 0;
    if ($opts{e} =~ /\.gz$/){
        open ($EXON_BED, "gzip -dc $opts{e} |") or die "Can't open $opts{e} via gzip: $!\n";
    }else{
        open ($EXON_BED, "<", $opts{e}) or die "Can't open $opts{e} for reading: $!\n";
    }
    my $output = "$opts{o}_exon_allele_counts.tsv";
    open (my $OUT, ">", $output) or die "Can't open $opts{o} for writing: $!\n";
    print $OUT join
    (
        "\t",
        "exon",
        "strand",
        "type",
        map { "exon_$_"    } @ac_fields,
    ) . "\n";
    my $n = 0;
    while (my $line = <$EXON_BED>){
        $n++;
        next if $line =~ /^#/;
        chomp $line;
        my ($chr, $start, $end, $type, $score, $strand) = split("\t", $line); 
        $start++; #BED format is 0-based
        my @counts = getVariantCounts 
        (
            $chr,
            $start,
            $end,
        );
        print $OUT join
        (
            "\t",
            "$chr:$start-$end",
            $strand,
            $type,
            @counts
        ) . "\n";
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }
    $next_update = $progressbar->update($n) if $n >= $next_update;
    close $EXON_BED;
}

#################################################
sub processIntronBed{
    print STDERR "Processing intron bed file $opts{i}...\n";
    my $INTRON_BED;
    if ($opts{i} =~ /\.gz$/){
        open ($INTRON_BED, "gzip -dc $opts{i} |") or die "Can't open $opts{i} via gzip: $!\n";
    }else{
        open ($INTRON_BED, "<", $opts{i}) or die "Can't open $opts{i} for reading: $!\n";
    }
    my $line_count = 0;
    $line_count += tr/\n/\n/ while sysread($INTRON_BED, $_, 2 ** 20);
    print STDERR "$opts{i} has $line_count lines...\n";
    my $progressbar = Term::ProgressBar->new
    (
        { 
          name => "Counting",
          count => ($line_count), 
          ETA => "linear" 
        } 
    );
    my $next_update      = 0;
    if ($opts{i} =~ /\.gz$/){
        open ($INTRON_BED, "gzip -dc $opts{i} |") or die "Can't open $opts{i} via gzip: $!\n";
    }else{
        open ($INTRON_BED, "<", $opts{i}) or die "Can't open $opts{i} for reading: $!\n";
    }
    my $output = "$opts{o}_intron_allele_counts.tsv";
    open (my $OUT, ">", $output) or die "Can't open $opts{o} for writing: $!\n";
    print $OUT join
    (
        "\t",
        "intron",
        "strand",
        "type",
        (map { "donor_invariant_$_"    } @ac_fields),
        (map { "acceptor_invariant_$_" } @ac_fields),
        (map { "donor_consensus_$_"    } @ac_fields),
        (map { "acceptor_consensus_$_" } @ac_fields),
    ) . "\n";
    my $n = 0;
    while (my $line = <$INTRON_BED>){
        $n++; 
        next if $line =~ /^#/;
        chomp $line;
        my ($chr, $start, $end, $type, $score, $strand) = split("\t", $line); 
        $start++; #BED format is 0-based
        if ($strand ne '+' and $strand ne '-'){
            die "Could not determine strand for region:\n$line\n";
        }
        my (@d_inv, @a_inv, @d_cons, @a_cons);
        push @d_inv, $strand eq '+' ? $start     : $end ;
        push @d_inv, $strand eq '+' ? $start + 1 : $end - 1;
        push @a_inv, $strand eq '+' ? $end - 1 : $start + 1;
        push @a_inv, $strand eq '+' ? $end     : $start;

        push @d_cons, $strand eq '+' ? $start - 3  : $end + 3;
        push @d_cons, $strand eq '+' ? $start + 10 : $end - 10;
        push @a_cons, $strand eq '+' ? $end - 13 : $start + 13;
        push @a_cons, $strand eq '+' ? $end + 3  : $start - 3;
        
        @d_inv  = sort { $a <=> $b } @d_inv;
        @a_inv  = sort { $a <=> $b } @a_inv;
        @d_cons = sort { $a <=> $b } @d_cons;
        @a_cons = sort { $a <=> $b } @a_cons;
        
        my @counts = ();         
        foreach my $coords 
        (
            \@d_inv,
            \@a_inv,
            \@d_cons,
            \@a_cons,
        ){
            push @counts,getVariantCounts 
            (
                $chr,
                $coords->[0],
                $coords->[1],
            );
        }
        print $OUT join
        (
            "\t",
            "$chr:$start-$end",
            $strand,
            $type,
            @counts
        ) . "\n";
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }
    $next_update = $progressbar->update($n) if $n >= $next_update;
    close $INTRON_BED;    
}

#################################################
sub getVariantCounts{
    my ($chr, $start, $end) = @_;
    $chr =~ s/^chr//;#ExAC is on b37 style refererence
    if (coverageAboveThreshold($chr, $start, $end)){
       #TODO get allele counts form VCF 
        my @hits = VcfReader::searchByRegion(
            %search_args,
            chrom => $chr,
            start => $start,
            end   => $end,
        );
        my %counts = map {$_ => 0} @ac_fields; 
        foreach my $h (@hits){
            my @split = split("\t", $h);
            for (my $i = 0; $i < @ac_fields; $i++){
                my $ac = VcfReader::getVariantInfoField
                ( 
                    \@split, 
                    $ac_fields[$i] 
                );
                my @ac_per_allele = split(",", $ac); 
                if (@min_an){
                    my $an = VcfReader::getVariantInfoField
                    ( 
                        \@split, 
                        $an_fields[$i] 
                    );
                    my $min = @min_an > 1 ? $min_an[$i] : $min_an[0];
                    if ($an >= $min){
                        $counts{$ac_fields[$i]} += sum(@ac_per_allele);
                    }
                }else{
                    $counts{$ac_fields[$i]} += sum(@ac_per_allele);
                }
            }
        }
        return map { $counts{$_} } @ac_fields;
    }else{
        return ("NA") x @ac_fields;
    }    
}

#################################################
sub coverageAboveThreshold{
    my ($chr, $start, $end) = @_;
    return 1 if not $opts{p};
    return 1 if not $opts{t};
    if (not exists $cov_iters{$chr}){
        my $cov_file = "$opts{c}/Panel.chr$chr.coverage.txt.gz";
        if (not -e $cov_file){
            warn "Could not find coverage file ($cov_file) for chromosome $chr ".
                "in $opts{c}!\n";
            return 0;
        }
        $cov_iters{$chr} = VcfReader::getTabixIterator($cov_file);
    }
    my $iter = $cov_iters{$chr}->query("$chr:$start-" . ($end + 1) );
    my $col = $cov_columns{$opts{t}};
    my $base_count = 0;
    while (my $m = $iter->next() ){ 
        my $cov = (split("\t", $m))[$col];
        return 0 if $cov < $opts{p};#require all coordinates to be > $opts{p}
        $base_count++;
    }
    return 0 if $base_count < (1 + $start - $end); #require all bases to be present in cov file
    return 1;
}
    
#################################################
sub checkVcf{
    my @head = VcfReader::getHeader($opts{v});
    die "Header not ok for VCF input ($opts{v}) "
        if not VcfReader::checkHeader( header => \@head );
    my %info = VcfReader::getInfoFields( header => \@head);
    foreach my $c (@ac_fields){
        if (not exists $info{$c}){
            die "Allele count field '$c' does not exist in INFO header for $opts{v}!\n";
        }
    }
    foreach my $c (@an_fields){
        if (not exists $info{$c}){
            die "Allele number field '$c' does not exist in INFO header for $opts{v}!\n";
        }
    }
}

#################################################
sub checkCoverageArgs{
    my @cov = qw /1 5 10 15 20 25 30 50 100/;
    return if not $opts{t};
    if (not $opts{c}){
        usage("ERROR: -c/--exac_coverage_dir option is required for use with ".
              "-t/--coverage_threshold argument!\n");
    }
    my $n = 0;
    my %valid = map {$_ => $n++} @cov;
    if (not exists $valid{$opts{t}}){
        die "Coverage threshold (-t/--coverage_threshold) must be one of the following values: " 
            .join(", ", @cov) . "\n";
    }
    if ($opts{p} < 0 or $opts{p} > 1){
        die "-p/--proportion_at_threshold argument must be between 0.00 and 1.00\n"
    }
    return %valid;
}

#################################################
sub usage{
    my $msg = shift;
    print STDERR "\n$msg\n" if $msg;

    print STDERR <<EOT


DESCRIPTION: Prints no. variants found in ExAC at a given coverage level for intron consensus sequences

USAGE: $0 -b introns.bed -v exac.vcf -c /path/to/exac/coverage/files/ -o output_prefix

OPTIONS:
    
    -i,--intron_bed FILE
        intron BED file created by intronsToBed.pl
    
    -e,--exon_bed FILE
        exon BED file created by intronsToBed.pl
    
    -v,--exac_vcf FILE
        VCF file of ExAC variants

    -c,--exac_coverage_dir DIR
        Directory containing ExAC coverage data

    -o,--output FILE
        Output prefix. Output files will be named 
        <prefix>_intron_allele_counts.tsv and 
        <prefix>_exon_allele_counts.tsv. 
    
    -t,--coverage_threshold INT
        Minimum coverage threshold for including a region. 
        Allowable values are 1, 5, 10, 15, 20, 25, 30, 50 or 100.
        Default = 20.

    -p,--proportion_at_threshold FLOAT
        Minimum proportion of samples (between 0.00 and 1.00) at coverage 
        threshold for including a reigon. Default = 0.9.

    -h,--help 
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}

