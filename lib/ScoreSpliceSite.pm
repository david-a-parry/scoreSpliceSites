package ScoreSpliceSite;

use strict; 
use warnings; 
use Carp;
use FindBin qw($RealBin);

my $data_dir = "$RealBin/data/species";

my %matrices = ();
my %min_scores = ();#min scores for each matrix
my %max_scores = ();#max scores for each matrix
my %best_seqs = ();
my %worst_seqs = ();
my %pos_by_weight = (); 
my @species = qw
(
    10090
    3702 
    6239 
    7227
    9606
);
my @introns = qw 
/
    AT_AC_U12 
    GC_AG_U2
    GT_AG_U12
    GT_AG_U2
/;

foreach my $s (@species){
    foreach my $int (@introns){
        my $s_dir = "$data_dir/$s/$int";
        if (not -d $s_dir){
            croak "Could not find splice data directory '$s_dir' ";
        }
        foreach my $pwm (qw / A D B / ){ 
            #my $f = "$s_dir/$pwm"."_logo";
            my $f = "$s_dir/$pwm.matrix";
            next if $int =~ /U12$/ and $s == 6239;
            next if $pwm eq 'B' and not -e $f;#branch site file not available for every species
            $matrices{$s}->{$int}->{$pwm} = _readLogoFile($f);
            (
                $min_scores{$s}->{$int}->{$pwm}, 
                $max_scores{$s}->{$int}->{$pwm}, 
                $worst_seqs{$s}->{$int}->{$pwm}, 
                $best_seqs{$s}->{$int}->{$pwm}, 
            ) = matrixMinMax($matrices{$s}->{$int}->{$pwm});
 
            @{$pos_by_weight{$s}->{$int}->{$pwm}} = _sortPosByWeight
            (
                $matrices{$s}->{$int}->{$pwm},
                $best_seqs{$s}->{$int}->{$pwm},
            );
        }
=cut
        if ($int =~ /U12$/){
            my $f = "$s_dir/B_logo";
            $matrices{$s}->{$int}->{B} = _readLogoFile($f);
        }
=cut
    }
}

sub _readLogoFile{
    my $f = shift;
    open (my $PWM, $f) or croak "Can't open splice data file '$f': $! "; 
    my $in_matrix = 0;
    my @nt = qw /A C G T/; 
    my @weights = ();
    while (my $line = <$PWM>){
        chomp $line;
        if ($in_matrix){
            my @s = split(/\s+/, $line); 
            my %h = map { $nt[$_] => $s[$_] } 0 .. $#s; 
            push @weights, \%h; 
        }elsif ($line =~ /A\s+C\s+G\s+T/){
            $in_matrix++;
        }
    }
    unless ( @weights ){
        croak "No splice weight data found in $f ";
    }
    return \@weights;
}

sub score{
#acceptor splites require 12 upstream, the splice site (2nt) and 3 downstream
#donor splites require 3 upstream (exonic), the splice site (2nt) and 8 downstream
    my %args = @_;
    my $score = 0;
    foreach my $r (qw / seq species type site / ){ 
        croak "$r argument is required for scoreSpliceSite method " 
          if not exists $args{$r};
    }
    my $m = $matrices{$args{species}}->{$args{type}}->{$args{site}};
    if (length($args{seq})  <  @$m){
        my $l = @$m;
        carp "Sequence '$args{seq}' is too short (< $l) to score ";
        return ;
    }
    my $min = $min_scores{$args{species}}->{$args{type}}->{$args{site}};
    my $max = $max_scores{$args{species}}->{$args{type}}->{$args{site}};
    return _scoreSequence($args{seq}, $m, $min, $max);
}

sub scanForBranchPoint{
    my %args = @_;
    my $score = 0;
    foreach my $r (qw / seq species type / ){ 
        croak "$r argument is required for scanForBranchPoint method " 
          if not exists $args{$r};
    }
    my $m = $matrices{$args{species}}->{$args{type}}->{B};
    my $l = @$m;
    if (length($args{seq})  <  $l){
        carp "Sequence '$args{seq}' is too short (< $l) to score ";
        return ;
    }
    my $min = $min_scores{$args{species}}->{$args{type}}->{'B'};
    my $max = $max_scores{$args{species}}->{$args{type}}->{'B'};
    my $best_score = $min;
    my $best_seq = '';
    for (my $i = 0; $i <= length($args{seq}) - $l; $i++){
        my $subseq = substr($args{seq}, $i, $l); 
        my $score = _scoreSequence($subseq, $m, $min, $max);
        if ($score >= $best_score){
            $best_score = $score;
            $best_seq = $subseq;
        }
    }
    return $best_score, $best_seq;
}

sub _scoreSequence{
    my ($seq, $m, $min, $max) = @_;
    my $score = 0;
    for (my $i = 0; $i < @$m; $i++){
        my $n = uc (substr($seq, $i, 1) );
        my $p = exists $m->[$i]->{$n} ? $m->[$i]->{$n}/0.25 : 0; 
        $p ||= 0.0001; 
        $score += log($p);
    }
#'A score close to 50 means the background model is favored, a score close 
#to 100 favors the splice-site motif in question, and a score close to 0 
#implies that the splice site is different from both the background and 
#the splice-site distribution' doi: 10.1093/nar/gkl556
    if ($score < 0){
        return 50 - (50 * $score / $min) ;
    }else{
        return 50 * ($score/$max) + 50;
    }
}


sub matrixMinMax{
    my $m = shift;
    my ($minscore, $maxscore) = 0, 0;
    my $best;
    my $worst;
    for (my $i = 0; $i < @$m; $i++){
        my $max = 0;
        my $min = 1;
        my $best_nt = 'A'; 
        my $worst_nt = 'A'; 
        foreach my $k (keys %{$m->[$i]}){   
            if ($max < $m->[$i]->{$k}){
                $best_nt = $k;
                $max = $m->[$i]->{$k};
            }
            if ($min > $m->[$i]->{$k}){
                $worst_nt = $k;
                $min = $m->[$i]->{$k};
            }
            #$max = $max >= $m->[$i]->{$k} ? $max : $m->[$i]->{$k};
            #$min = $min <= $m->[$i]->{$k} ? $min : $m->[$i]->{$k};
        }
        my $minp = $min/0.25; 
        $minp ||= 0.0001; 
        $minscore += log($minp);
        my $maxp = $max/0.25; 
        $maxp ||= 0.0001; 
        $maxscore += log($maxp);
        $best .= $best_nt;
        $worst .= $worst_nt;
    }
    return $minscore, $maxscore, $worst, $best;
}
 
sub _sortPosByWeight{
    my $m = shift;
    my $best_seq = shift;
    my @seq = split("", $best_seq);
    my %pos_to_score = ();
    for (my $i = 0; $i < @$m; $i++){
        $pos_to_score{$i} = $m->[$i]->{$seq[$i]};
    }
    return sort { $pos_to_score{$b} <=> $pos_to_score{$a} } keys %pos_to_score;
}
       
sub getIntronTypes{
    return @introns;
}

sub getSpecies{
    return @species;
}

sub getWorstSeq{
    my ($type, $site,  $species) = @_;
    return $worst_seqs{$species}->{$type}->{$site};
} 
sub getBestSeq{
    my ($type, $site,  $species) = @_;
    return $best_seqs{$species}->{$type}->{$site};
} 
 
sub getDonorConsensusCoords{
#for a given intron start site coordinate
#return the coordinates of the consensus seq start and end
#strand can be specified as 1 for + strand or -1 for - strand
#coordinates are returned in coordinate order, not relative to strand
    my ($intron_start, $strand) = @_;
    $strand ||= 1;
    return sort {$a <=> $b} 
    (
        ($intron_start - (3 * $strand) ) ,
        ($intron_start + (10 * $strand) ) ,
    ); 
}

sub getAcceptorConsensusCoords{
#for a given intron end site coordinate
#return the coordinates of the consensus seq start and end
#strand can be specified as 1 for + strand or -1 for - strand
#coordinates are returned in coordinate order, not relative to strand
    my ($intron_stop, $strand) = @_;
    $strand ||= 1;
    return sort {$a <=> $b} 
    (
        ($intron_stop - (13 * $strand)) ,
        ($intron_stop + (3 * $strand) ) ,
    );
}

sub getBranchRegionCoords{
#for a given intron end site coordinate
#return the coordinates of the consensus seq start and end
#strand can be specified as 1 for + strand or -1 for - strand
#coordinates are returned in coordinate order, not relative to strand
    my ($intron_stop, $strand) = @_;
    $strand ||= 1;
    return sort {$a <=> $b} 
    (
        ($intron_stop - (100 * $strand)) ,
        ($intron_stop - (8 * $strand) ) ,
    );
}


sub getPosByWeight{
#returns array ref of consensus seq positions (0-based) in order 
#of the positions with the greatest weight according to the PWM
    my ($type, $site,  $species) = @_;
    return $pos_by_weight{$species}->{$type}->{$site};
} 

sub determineIntronType{
#requires: 
#          donor => last 3 bp of exon plus first 10 bp of intron
#          branch => around the last 100 bp of intron, not including last 8 bp
#          species => taxonomic code (e.g. 9606 for human)
#          acceptor => last 14bp of intron + 3 of exon
    my %args = @_;
    my %scores = ();
    my %branch_seqs = (); 
    my $d_term = substr($args{donor}, 3, 2);
    my $a_term = substr($args{acceptor}, 12, 2);
    my $term_type = $d_term . "_" . $a_term; #e.g. GT_AG
    foreach my $type (@introns){
        $scores{'D'}->{$type} = score
        (
            seq  => $args{donor},
            type => $type,
            site => 'D',
            species => $args{species},
        );
        $scores{'A'}->{$type} = score
        (
            seq  => $args{acceptor},
            type => $type,
            site => 'A',
            species => $args{species},
        );
        ($scores{'B'}->{$type}, $branch_seqs{$type}) = 
        scanForBranchPoint
        (
            seq  => $args{branch},
            type => $type,
            species => $args{species},
        );
    }
    my $u12_b_score;
    if ($scores{'B'}->{AT_AC_U12} > $scores{'B'}->{GT_AG_U12}){
        $u12_b_score = $scores{'B'}->{AT_AC_U12};
    }else{
        $u12_b_score = $scores{'B'}->{GT_AG_U12};
    }
    my $best_u2;
    my $best_u12;
    foreach my $type (@introns){
        if ($type =~ /U12$/){
            if ($best_u12){
                 if ($scores{'D'}->{$type} > $scores{'D'}->{$best_u12} and 
                     $scores{'A'}->{$type} > 50
                 ){
                    $best_u12 = $type;
                 }
            }else{
                $best_u12 = $type;
            }
        }elsif($type =~ /U2$/){
            if ($best_u2){
                 if ($scores{'D'}->{$type} > $scores{'D'}->{$best_u2} and 
                     $scores{'A'}->{$type} > 50
                 ){
                    $best_u2 = $type;
                 }
            }else{
                $best_u2 = $type;
            }
        }
    }
    
    my $pick = pickU12orU2 
    (
        scores    => \%scores,
        u12branch => $u12_b_score,
        U12       => $best_u12, 
        U2        => $best_u2, 
    );
    return 0 if not $pick;
    if ($pick =~ /^$term_type(_U12|_U2)$/){
        return $pick;
    }else{
        if ($scores{'D'}->{$pick} > 60){
            if ($pick =~ /(_U1*2)$/){
                return "NonCanon_$term_type$1";
            }else{
                carp "Could not parse intron type '$pick' ";
            }
        }
        return 0;
    }
}

sub pickU12orU2{
    my %args = @_;
    if ($args{scores}->{'D'}->{$args{U12}} < 50 
        and $args{scores}->{'D'}->{$args{U2}} < 50){
        #classify anything with both scores below 50 as
        # UNKNOWN
        return 0;
    }
    if ($args{scores}->{'D'}->{$args{U12}} - 
        $args{scores}->{'D'}->{$args{U2}} >= 25 and 
        $args{scores}->{'A'}->{$args{U12}} > 60
    ){
        #we annotate as U12 if the donor site score is 
        #at least 25 more than a U2 site and acceptor site is above 60
        return $args{U12};
    }elsif ($args{scores}->{'D'}->{$args{U12}} - 
            $args{scores}->{'D'}->{$args{U2}} >= 10 and
            $args{scores}->{'A'}->{$args{U12}} > 60
    ){
        #otherwise if U12 score is at least 10 better we annotate
        # as U12 if there's a 'good' (score >= 65) branch point
        #and acceptor site is above 60
        if ($args{u12branch} >= 65){
            return $args{U12};
        }
    }

    if ($args{scores}->{'A'}->{$args{U2}} > 50){
        return $args{U2};
    }else{
        return 0;
    }
}
1;


