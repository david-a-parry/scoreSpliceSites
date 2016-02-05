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
    
1;


