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
use Bio::DB::Sam;
use lib "$RealBin/lib";
use ScoreSpliceSite;
use ReverseComplement qw/ reverse_complement /;
use lib "$RealBin/lib/dapPerlGenomicLib/";
use VcfReader;

my %opts = (s => 9606);
my %info_fields = (); 
GetOptions(
    \%opts,
    'i|intron_bed=s',
    'f|fasta=s',
    's|species=s',
    'l|lenient',
    'o|output=s',
    'v|vcf=s',
    'a|all',
    'c|chr',
    'h|?|help',
) or usage("Error getting options!");

usage() if $opts{h};
checkArgs();
my @head = VcfReader::getHeader($opts{v});
die "Header not ok for VCF input ($opts{v}) "
    if not VcfReader::checkHeader( header => \@head );
my %search_args = VcfReader::getSearchArguments( $opts{v} );
 
my $fai = Bio::DB::Sam::Fai->load($opts{f});#should create index if it doesn't exist

my $OUT;
if ($opts{o}){
    open ($OUT, ">", $opts{o}) or die "Can't open $opts{o} for writing: $!\n";    
}else{
    $OUT = \*STDOUT;
}

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

printHeader();

my $n = 0;
while (my $line = <$INTRON_BED>){
    $n++; 
    next if $line =~ /^#/;
    chomp $line;
    my ($chr, $start, $end, $type, $score, $strand) = split("\t", $line); 
    $type =~ s/\w+\///;
    $start++; #BED format is 0-based
    #get donor and acceptor flank regions
    my @regions = getSpliceRegions
    (
        chrom  => $chr,
        start  => $start,
        end    => $end,
        strand => $strand
    );
 
    #get variants in these regions
    my @variants = ();
    if ($opts{c}){
        $chr =~ s/^chr//;
    }
    foreach my $reg (@regions){
        push @variants, VcfReader::searchByRegion
        (
            %search_args,
            chrom => $chr,
            start => $reg->{start},
            end   => $reg->{end},
        );
    
    }
    
    #get DNA for region and flanks - we'll only deal with variants altering
    # 100 bp span of intron
    my $dna = $fai->fetch("$chr:" . ($start - 100) . "-" . ($end + 100) );
    #score original intron
    my ($original_type, $original_scores) = scoreIntron
    (
        seq    => $dna, 
        offset => $start - 100,
        chrom  => $chr,
        start  => $start,
        end    => $end,
        strand => $strand,
    ); 
    #for each variant mutate sequence and rescore splice site 
VAR: foreach my $var (@variants){
        my @var_split = split("\t", $var); 
        my %var_scores = (); 
        my @intron_classes = ($original_type); 
        foreach my $site (keys %info_fields){
            foreach my $class (keys %{$info_fields{$site}}){
                push @{$var_scores{$site}->{$class}}, 
                  $original_scores->{$class}->{$site};
            }
        }
        my ($mutated, $start_offset, $end_offset) = mutateSeq
        (
            seq     => $dna, 
            variant => \@var_split,
            offset  => $start - 100,
            start   => $start,
            end     => $end,
        ); 
        for (my $i = 0; $i <@$mutated; $i++){#mutated is in order of allele number
            if (not $mutated->[$i] #allele not scored
                or ($end + $end_offset->[$i]) - ($start + $start_offset->[$i]) 
                < 20 #bulk of intron is deleted 
            ){
                push @intron_classes, '.';
                foreach my $site (keys %info_fields){
                    foreach my $class (keys %{$info_fields{$site}}){
                        push @{$var_scores{$site}->{$class}}, '.';
                    }
                }
            }else{
                my ($mut_type, $mut_scores) = scoreIntron
                (
                    seq    => $mutated->[$i], 
                    offset => $start - 100,
                    chrom  => $chr,
                    start  => $start + $start_offset->[$i],
                    end    => $end + $end_offset->[$i],
                    strand => $strand,
                ); 
                push @intron_classes, $mut_type;
                foreach my $site (keys %info_fields){
                    foreach my $class (keys %{$info_fields{$site}}){
                        push @{$var_scores{$site}->{$class}}, 
                          $mut_scores->{$class}->{$site};
                    }
                }
            }
        }
        unless ($opts{a}){
            my %all_classes =  map {$_ => undef} grep { $_ ne '0' } grep { $_ ne '.' } @intron_classes;
            next VAR if (keys %all_classes < 2);#no change to intron class
        }
        $var = VcfReader::addVariantInfoField
        (
            line  => \@var_split, 
            id    => 'IntronClass',
            value => join(",", @intron_classes),
        );
        $var = VcfReader::addVariantInfoField
        (
            line  => $var, 
            id    => 'IntronStrand',
            value => $strand,
        );
        $var = VcfReader::addVariantInfoField
        (
            line  => $var, 
            id    => 'IntronStart',
            value => $start,
        );
        $var = VcfReader::addVariantInfoField
        (
            line  => $var, 
            id    => 'IntronEnd',
            value => $end,
        );
        foreach my $site (keys %info_fields){
            foreach my $class (keys %{$info_fields{$site}}){
                $var = VcfReader::addVariantInfoField
                (
                    line  => $var, 
                    id    => $info_fields{$site}->{$class},
                    value => join(",", @{$var_scores{$site}->{$class}}),
                );
            }
        } 
        print $OUT join("\t", @$var) . "\n";
    }
    $next_update = $progressbar->update($n) if $n >= $next_update;
}
$next_update = $progressbar->update($n) if $n >= $next_update;
close $INTRON_BED;


#################################################
sub mutateSeq{
    my %args = @_;
    my @new_seqs = ();
    my @start_offsets  = (); #indels will shift the relative position of the splice site
    my @end_offsets = (); #indels will shift the relative position of the splice site
    #get minimal representation of each allele
    my %min = VcfReader::minimizeAlleles($args{variant});
    foreach my $k (sort {$a <=> $b} keys %min){
        #make sure we have a valid allele (not <DEL> or *)
        if ($min{$k}->{ALT} !~ /^[ACTG]$/ 
            #make sure start position of allele is within range of our seq
          or $min{$k}->{POS} - $args{offset} < 0 
            #make sure end position of allele is within range of our seq
          or $min{$k}->{POS} > $args{offset} + length($args{seq})
        ){
            push @new_seqs, undef;
            push @start_offsets, undef;
            push @end_offsets, undef;
            next;
        }
        #get sequence leading up to mutation
        my $mut = substr
        (
            $args{seq}, 
            0,
            $min{$k}->{POS} - $args{offset},
        );
        #append mutant allele
        $mut .= $min{$k}->{ALT};
        #add sequence after mutation
        $mut .= substr
        (
            $args{seq},
            $min{$k}->{POS} - $args{offset} + length($min{$k}->{REF}),
        );    
        push @new_seqs, $mut;
        my $diff = (length($min{$k}->{ALT})) - (length($min{$k}->{REF}));
        if ($min{$k}->{POS} < $args{start}){
            my $s_diff = $diff; 
            if ( $diff < 0 ){ #deletion 
            # - make sure we only increment by the amount deleted prior to start
                if ($min{$k}->{POS} + length($min{$k}->{REF}) > $args{start}){
                    $s_diff = $min{$k}->{POS} - $args{start}; 
                }
            }
            push @start_offsets, $s_diff;
        }else{
            push @start_offsets, 0;
        }
        if ($min{$k}->{POS} < $args{end}){
            my $s_diff = $diff; 
            if ( $diff < 0 ){ #deletion 
            # - make sure we only increment by the amount deleted prior to end 
                if ($min{$k}->{POS} + length($min{$k}->{REF}) > $args{end}){
                    $s_diff = $min{$k}->{POS} - $args{end}; 
                }
            }
            push @end_offsets, $s_diff;
        }else{
            push @end_offsets, 0;
        }
    }
    return (\@new_seqs, \@start_offsets, \@end_offsets);
}

#################################################
sub scoreIntron{
    my %args = @_;
    my ($strand, $i_start, $i_end) = getIntronStartAndEnd(%args);
    my %seqs = ();
    my %scores = ();
    my ($d_start, $d_end) = ScoreSpliceSite::getDonorConsensusCoords
    (
        $i_start,
        $strand,
    );
    $seqs{D} = substr
    (
        $args{seq}, 
        $d_start - $args{offset}, 
        $d_end - $d_start + 1
    );
    my ($a_start, $a_end) = ScoreSpliceSite::getAcceptorConsensusCoords
    (
        $i_end,
        $strand,
    );
    $seqs{A} = substr
    (
        $args{seq}, 
        $a_start - $args{offset}, 
        $a_end - $a_start + 1
    );
    my ($b_start, $b_end) = ScoreSpliceSite::getBranchRegionCoords
    (
        $i_end,
        $strand,
    ); 
    $seqs{B} = substr
    (
        $args{seq}, 
        $b_start - $args{offset}, 
        $b_end - $b_start + 1
    );
    if ($strand < 0 ){
        foreach my $k (keys %seqs){
            $seqs{$k} = reverse_complement($seqs{$k});
        }
    }
    my $seq_type = ScoreSpliceSite::determineIntronType
    (
            donor    => $seqs{D}, 
            branch   => $seqs{B},
            acceptor => $seqs{A},
            species  => $opts{s},
    ); 
    foreach my $type (ScoreSpliceSite::getIntronTypes()){
        foreach my $s (qw /A D/){
            $scores{$type}->{$s} = ScoreSpliceSite::score
            (
                seq     => $seqs{$s},
                type    => $type,
                site    => $s,
                species => $opts{s},
            );
        }
        ($scores{$type}->{B}) = ScoreSpliceSite::scanForBranchPoint
        (
            seq     => $seqs{B},
            type    => $type,
            species => $opts{s},
        );
    }
    return ($seq_type, \%scores); 
}

#################################################
sub printHeader{
    #print original meta header lines
    print $OUT join("\n", @head[0..$#head-1]) . "\n";
    
    #create new INFO entries for changed splice scores
    foreach my $site (qw /donor acceptor branch/){
        my $initial = uc(substr($site, 0, 1));
        foreach my $type (ScoreSpliceSite::getIntronTypes()){
            my $id = join("_", $site, "score", $type);
            print $OUT "##INFO=<ID=$id,Number=R,"
              . "Type=Float,Description=\"$type $site splice "
              . "consensus score for each allele\">\n";
            $info_fields{$initial}->{$type} = $id;
        }
    }
    print $OUT "##INFO=<ID=IntronClass,Number=R,Type=STRING,"
      . "Description=\"Class of intron being classified for each allele\">\n";
    print $OUT "##INFO=<ID=IntronStart,Number=1,Type=Integer,"
      . "Description=\"Start coordinate of intron being classified\">\n";
    print $OUT "##INFO=<ID=IntronEnd,Number=1,Type=Integer,"
      . "Description=\"End coordinate of intron being classified\">\n";
    print $OUT "##INFO=<ID=IntronStrand,Number=1,Type=String,"
      . "Description=\"Strand of intron being classified\">\n";
    print $OUT "$head[-1]\n";
}



#################################################
sub getIntronStartAndEnd{
    my %args = @_; 
    my ($strand, $i_start, $i_end);
    if ($args{strand} eq '+'){
        $strand = 1;
        ($i_start, $i_end) = ($args{start}, $args{end});
    }elsif($args{strand} eq '-'){
        $strand = -1;
        ($i_start, $i_end) = ($args{end}, $args{start});
    }else{
        die "Do not understand strand value '$strand' for region "
          . "$args{chrom}:$args{start}-$args{end}\n";
    }
    return ($strand, $i_start, $i_end);
}

#################################################
sub getSpliceRegions{
    my %args = @_;
    my ($strand, $i_start, $i_end) = getIntronStartAndEnd(%args);
    
    my ($d_start, $d_end) = ScoreSpliceSite::getDonorConsensusCoords
    (
        $i_start,
        $strand,
    );
    my %donor = 
    (
        start => $d_start,
        end   => $d_end,
    );
    my ($a_start, $a_end) = ScoreSpliceSite::getAcceptorConsensusCoords
    (
        $i_end,
        $strand,
    );
    my ($b_start, $b_end) = ScoreSpliceSite::getBranchRegionCoords
    (
        $i_end,
        $strand,
    );
    my $f_start = $strand > 0 ? $b_start : $a_start;
    my $f_end   = $strand > 0 ? $a_end   : $b_end;

    my %acc_branch = 
    (
        start => $f_start,
        end   => $f_end,
    );
    return (\%donor, \%acc_branch);
}

#################################################
sub checkArgs{
    usage("ERROR: -i/--intron_bed option is required!\n") if not $opts{i};
    usage("ERROR: -f/--fasta option is required!\n") if not $opts{f};
    usage("ERROR: -v/--vcf option is required!\n") if not $opts{v};

    if (not grep {$opts{s} eq $_} ScoreSpliceSite::getSpecies){
        die "ERROR: Splice consensus for species '$opts{s}' not available.\n";
    }
}

#################################################
sub usage{
    my $msg = shift;
    print STDERR "\n$msg\n" if $msg;

    print STDERR <<EOT


DESCRIPTION: TODO

USAGE: $0 -i introns.bed -f genome.fasta -v variants.vcf[options]

OPTIONS:
    
    -i,--intron_bed FILE
        intron BED file created by intronsToBed.pl
    
    -f,--fasta FILE
        Genome fasta file for retrieving DNA sequences for donor and acceptor 
        splice sequences.

    -v,--vcf FILE
        VCF file of variants to check for those altering splice class.

    -o,--output FILE
        Output file name. Default = STDOUT.
    
    -s,--species INT
        Taxonomic code for species to use for splice prediction. 
        Default is 9606 (human). 
        Available species are 10090 (mouse), 3702 (A. thaliana), 
        6239 (C. elegans), 7227 (D. melanogaster) and 9606 (human).
   
    -c,--chr
        Use this option to trim the 'chr' from the beginning of contig names 
        when searching VCFs and retrieving DNA sequences from FASTA files (e.g.
        if your bed file uses the 'chr' naming convention but your VCF and 
        FASTA files BOTH do not). 

    -a,--all
        Use this option to output all variants even if they do not alter
        an intron type. 

    -h,--help 
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}
