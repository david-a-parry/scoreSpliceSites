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

my @mutate_from = ();
my @mutate_to   = ();
    
my %opts = 
(
    m => \@mutate_from,
    t => \@mutate_to,
    s => 9606,
);

GetOptions(
    \%opts,
    'i|intron_bed=s',
    'f|fasta=s',
    's|species=s',
    'l|lenient',
    'o|output=s',
    'm|mutate=s{,}', 
    't|to=s{,}',
    'h|?|help',
) or usage("Error getting options!");

usage() if $opts{h};
checkArgs();

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
my $n = 0;
while (my $line = <$INTRON_BED>){
    $n++; 
    next if $line =~ /^#/;
    chomp $line;
    my ($chr, $start, $end, $type, $score, $strand) = split("\t", $line); 
    $type =~ s/\w+\///;
    next if not grep {$_ eq $type} @mutate_from;#skip types we're not mutating
    my ($donor_seq, $acceptor_seq, $branch_seq, $acc_flank) = getSpliceSeqs
    (
        $chr,
        $start,
        $end,
        $strand,
        $type,
    );
    #check score for donor sequence matches score in BED
    #this is really just a sanity check and can be skipped
    #using --lenient flag
    checkScore($donor_seq, $score, $type, $line);
    
    #mutate sequence
    my ($mutations, $m_seqs, $new_type) = mutateSpliceSeqs
    (
        $donor_seq,
        $acceptor_seq,
        $acc_flank ,
        $type,
    ); 
    
    #output results
    if (not defined $mutations){
        print $OUT "$line" . ("\t-" x 7) . "\n";
    }else{
        print $OUT "$line";
        print $OUT "\tD=$donor_seq|A=$acceptor_seq|B=$branch_seq";
        print $OUT "\t" . join("\|", map { "$_=$m_seqs->{$_}" } qw /D A B/) ;
        my @mut_strings = ();
        my $n_mutations = 0; 
        foreach my $s (qw /D A B/){
            my @site_mut_strings = ();
            foreach my $pos (
              sort { $a <=> $b } keys %{$mutations->{$s}}
            ){
                push @site_mut_strings, $pos . $mutations->{$s}->{$pos};
                $n_mutations++;
            }
            push @mut_strings, "$s=" . join(",", @site_mut_strings) ;
        }
        my $new_score = ScoreSpliceSite::score
        (
            seq  => $m_seqs->{D},
            type => $new_type,
            site => 'D',
            species => $opts{s},
        );
        my $old_score = ScoreSpliceSite::score
        (
            seq  => $m_seqs->{D},
            type => $type,
            site => 'D',
            species => $opts{s},
        );
        print $OUT "\t" . join("\|", @mut_strings);
        print $OUT "\t" . join("\t", $n_mutations, $new_type, $old_score, $new_score);
        print $OUT "\n";
    }
    $next_update = $progressbar->update($n) if $n >= $next_update;
}
$next_update = $progressbar->update($n) if $n >= $next_update;
close $INTRON_BED;    

#################################################
sub mutateSpliceSeqs{
    my ($donor, $acceptor, $flank, $type) = @_;
    my %mutations = (); 
    my %mut_seqs  = ();
    my @don = split(//, uc($donor));
    my @ideal_donor = split 
    (   //,
        ScoreSpliceSite::getBestSeq
        (
            $type,
            'D',
            $opts{s},
        )
    );
    my $pos_by_weight_d = ScoreSpliceSite::getPosByWeight
    (
        $type,
        'D',
        $opts{s},
    );
    my @ideal_branch = split 
    (   //,
        ScoreSpliceSite::getBestSeq
        (
            $type,
            'B',
            $opts{s},
        )
    );
    my $pos_by_weight_b = ScoreSpliceSite::getPosByWeight
    (
        $type,
        'B',
        $opts{s},
    );
    my @ideal_acceptor = split
    (   
        //,
        ScoreSpliceSite::getBestSeq
        (
            $type,
            'A',
            $opts{s},
        )
    );
    my $pos_by_weight_a = ScoreSpliceSite::getPosByWeight
    (
        $type,
        'A',
        $opts{s},
    );
    foreach my $m (@mutate_to){
        my ($new_acceptor, $new_flank) = ($acceptor, $flank);
        $mutations{$m} = 
        {
            D => {},
            B => {},
            A => {},
        };
        #check we have a branch site
        my ($b_score, $branch) = ScoreSpliceSite::scanForBranchPoint
        (
            seq  => $flank,
            type => $m,
            species => $opts{s},
        );
        my $new_branch = $branch;
        if ($b_score < 60){
            ($mutations{$m}->{B}, $new_branch) = mutateToType
            (
                site      => 'B',
                to_type   => $m,
                seq       => $branch,
                score     => 60,
            );
            #in case of multiple matching branch sequences
            # replace last occurence of branch sequence - 
            # i.e. the one closest to splice junction
            $new_flank =~ s/(.*)$branch/$1$new_branch/; 
        }
        my $a_score = ScoreSpliceSite::score
        (
            seq  => $acceptor,
            type => $m,
            site => 'A',
            species => $opts{s},
        );
        #check we have acceptor site
        if ($a_score < 60){
            ($mutations{$m}->{A}, $new_acceptor) = mutateToType
            (
                site      => 'A',
                to_type   => $m,
                seq       => $acceptor,
                score     => 60,
            );
        }
        my @ideal_mutant_d = split
        (
            //,
            ScoreSpliceSite::getBestSeq
            (
                $m,
                'D',
                $opts{s},
            )
        );
        my $mutant_pbw_d = ScoreSpliceSite::getPosByWeight
        (
            $m,
            'D',
            $opts{s},
        );
DPOS:   foreach my $pos (@$mutant_pbw_d){
            #skip for now if same nt at both original and target consensus
            next if $ideal_mutant_d[$pos] eq $ideal_donor[$pos];
            if ($don[$pos] ne $ideal_mutant_d[$pos]){
                $mutations{$m}->{D}->{$pos} = "$don[$pos]>$ideal_mutant_d[$pos]";
                $don[$pos] = $ideal_mutant_d[$pos];
                my $new_d_score = ScoreSpliceSite::score
                (
                    seq  => join("", @don),
                    type => $type,
                    site => 'D',
                    species => $opts{s},
                );
                my $new_d_target_score = ScoreSpliceSite::score
                (
                    seq  => join("", @don),
                    type => $m,
                    site => 'D',
                    species => $opts{s},
                );
                if ($new_d_target_score > 70){
#if new sequence is favourable towards target consensus (value of 70 is arbitrary)
                    if ($new_d_target_score >= $new_d_score + 25){
#and is 25 better than its original designation... 
                        last DPOS;
                    }elsif($new_d_score < 50){
#or score for original type is less than 50... 
                        last DPOS;
                    }
                }
            }
        }
        
        my $success = intronIsOfType
        (
            donor    => join("", @don), 
            branch   => $new_flank,
            type     => $m,
            species  => $opts{s},
        );
        if ($success){
            $mut_seqs{$m} = 
            {
                D => join("", @don), 
                B => $new_branch,
                A => $new_acceptor,
            };
        }else{
            delete $mutations{$m};
        }
    }
    return (undef, undef, undef) if not keys %mutations;
    my $mut_type = conversionWithFewestMutations(%mutations);
    return ($mutations{$mut_type}, $mut_seqs{$mut_type}, $mut_type);
}

#################################################
sub conversionWithFewestMutations{
    my %mut = @_;
    my %changes = ();
    if (keys %mut == 1){
        return (keys %mut)[0];
    }
    foreach my $m (keys %mut){
        my $changes = 0;
        foreach my $site (keys %{$mut{$m}}){
            $changes{$m} += keys(%{$mut{$m}->{$site}});
        }
    }
    return (sort { $changes{$a} <=> $changes{$b} } keys %changes)[0];
}
#################################################
sub mutateToType{
    my %args = @_;
    my %mut = (); 
    my @seq = split(//, $args{seq});
    my @ideal_mutant = split
    (
        //,
        ScoreSpliceSite::getBestSeq
        (
            $args{to_type},
            $args{site},
            $opts{s},
        )
    );
    my $mutant_pbw = ScoreSpliceSite::getPosByWeight
    (
        $args{to_type},
        $args{site},
        $opts{s},
    );
    foreach my $pos (@$mutant_pbw){
        #skip for now if same nt at both original and target consensus
        if ($seq[$pos] ne $ideal_mutant[$pos]){
            $mut{$pos} = "$seq[$pos]>$ideal_mutant[$pos]";
            $seq[$pos] = $ideal_mutant[$pos];
            my $new_score = ScoreSpliceSite::score
            (
                seq  => join("", @seq),
                type => $args{to_type},
                site => $args{site},
                species => $opts{s},
            );
            last if $new_score >= $args{score};
        }
    }
    return (\%mut, join("", @seq));
}


#################################################
sub intronIsOfType{
    my %args = @_;
    my $seq_type = ScoreSpliceSite::determineIntronType(%args);
    return 1 if $args{type} eq $seq_type;
    return 0;
}

#################################################
sub getSpliceSeqs{
    my ($chr, $start, $end, $strand, $type) = @_;
    $start++; #BED format is 0-based
    my ($i_start, $i_end);
    if ($strand eq '+'){
        $strand = 1;
        ($i_start, $i_end) = ($start, $end);
    }elsif($strand eq '-'){
        $strand = -1;
        ($i_start, $i_end) = ($end, $start);
    }else{
        die "Do not understand strand value '$strand' for region "
          . "$chr:$start-$end\n";
    }
    
    #get sequence for donor consensus
    my ($d_start, $d_end) = ScoreSpliceSite::getDonorConsensusCoords
    (
        $i_start,
        $strand,
    );
    my $donor = $fai->fetch("$chr:$d_start-$d_end");
    #get sequence for acceptor consensus
    my ($a_start, $a_end) = ScoreSpliceSite::getAcceptorConsensusCoords
    (
        $i_end,
        $strand,
    );
    my $acceptor = $fai->fetch("$chr:$a_start-$a_end");
    my ($b_start, $b_end) = ScoreSpliceSite::getBranchRegionCoords
    (
        $i_end,
        $strand,
    );
    my $flank = $fai->fetch("$chr:$b_start-$b_end");
    my ($b_score, $branch) = ScoreSpliceSite::scanForBranchPoint
    (
        seq  => $flank,
        type => $type,
        species => $opts{s},
    );
    
    if ($strand < 0){
        foreach my $s ($donor, $acceptor, $branch, $flank){
            $s = reverse_complement($s);
        }
    } 
    return ($donor, $acceptor, $branch, $flank);
}

#################################################
sub checkScore{
    return if $opts{l};
    my ($d_seq, $score, $type, $line) = @_;
    
    my $d_score = ScoreSpliceSite::score
    (
        seq  => $d_seq,
        type => $type,
        site => 'D',
        species => $opts{s},
    );
    if (abs($d_score - $score) > 0.01){#comparing floats with differing precision
        die "Score '$score' does not match computed score ($d_score) for donor "
          . "sequence '$d_seq' from line:\n$line\n";
    }
}

#################################################
sub checkArgs{
    usage("ERROR: -i/--intron_bed option is required!\n") if not $opts{i};
    usage("ERROR: -f/--fasta option is required!\n") if not $opts{f};

    if (not grep {$opts{s} eq $_} ScoreSpliceSite::getSpecies){
        die "ERROR: Splice consensus for species '$opts{s}' not available.\n";
    }
    if (not @mutate_from){
        @mutate_from = qw / AT_AC_U12 GT_AG_U12 / ;
    }
    if (not @mutate_to){
        @mutate_to   = qw / GT_AG_U2 GC_AG_U2 / ;
    }
    @mutate_from = map {uc($_)} @mutate_from;
    @mutate_to   = map {uc($_)} @mutate_to;
    
    my @valid = ScoreSpliceSite::getIntronTypes();
    foreach my $m (@mutate_to){
        if (not grep {$_ eq $m} @valid){
            usage("ERROR: Intron type '$m' specified with --to argument is "
              . " not a valid intron type.\n");
        }
    }
    
    foreach my $m (@mutate_from){
        if (not grep {$_ eq $m} @valid){
            usage("ERROR: Intron type '$m' specified with --mutate argument is "
              . " not a valid intron type.\n");
        }
        if (grep {$_ eq $m} @mutate_to){
            usage("ERROR: Intron type '$m' is in both --mutate and --to lists\n");
        }
    }
}

#################################################
sub usage{
    my $msg = shift;
    print STDERR "\n$msg\n" if $msg;

    print STDERR <<EOT


DESCRIPTION: TODO

USAGE: $0 -i introns.bed -f genome.fasta [options]

OPTIONS:
    
    -i,--intron_bed FILE
        intron BED file created by intronsToBed.pl
    
    -o,--output FILE
        Output file name. Default = STDOUT.
    
    -f,--fasta FILE
        Genome fasta file for retrieving DNA sequences for donor and acceptor 
        splice sequences.

    -s,--species INT
        Taxonomic code for species to use for splice prediction. 
        Default is 9606 (human). 
        Available species are 10090 (mouse), 3702 (A. thaliana), 
        6239 (C. elegans), 7227 (D. melanogaster) and 9606 (human).
   
    -m,--mutate STRING [STRING2]
        Mutate introns of these types. Others will be ignored.
        Default = AT_AC_U12 and GT_AG_U12

    -t.--to STRING [STRING2]
        Mutate introns to one of these types (if more than one is specified the
        type that requires fewest mutations will be chosen). 
        Default = GT_AG_U2 and GC_AC_U2.

    -l.--lenient
        Use this option to prevent strict checking that score in BED file 
        matches computed donor site consensus score.

    -h,--help 
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}
