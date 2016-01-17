package ScoreSpliceSite;

use strict; 
use warnings; 
use Carp;
use FindBin qw($RealBin);

my $data_dir = "$RealBin/data/species";

my %matrices = ();
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
        foreach my $pwm (qw / A D / ){ 
            my $f = "$s_dir/$pwm"."_logo";
            $matrices{$s}->{$int}->{$pwm} = _readLogoFile($f);
        }
        if ($int =~ /U12$/){
            my $f = "$s_dir/B_logo";
            $matrices{$s}->{$int}->{B} = _readLogoFile($f);
        }
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
    my ($min, $max) = matrixMinMax($m); 
    for (my $i = 0; $i < @$m; $i++){
        my $n = uc (substr($args{seq}, $i, 1) );
        my $p = $m->[$i]->{$n}/0.25; 
        $p ||= 0.0001; 
        $score += log($p)/log(2);
    }
    if ($score < 0){
        return 50 - (50 * $score / $min) ;
    }else{
        return 50 * ($score/$max) + 50;
    }
    #return 100 * ($score - $min)/($max - $min);
}

sub matrixMinMax{
    my $m = shift;
    my ($minscore, $maxscore) = 0, 0;
    for (my $i = 0; $i < @$m; $i++){
        my $max = 0;
        my $min = 0;
        foreach my $k (keys %{$m->[$i]}){
            $max = $max >= $m->[$i]->{$k} ? $max : $m->[$i]->{$k};
            $min = $min <= $m->[$i]->{$k} ? $min : $m->[$i]->{$k};
        }
        my $minp = $min/0.25; 
        $minp ||= 0.0001; 
        $minscore += log($minp)/log(2);
        my $maxp = $max/0.25; 
        $maxp ||= 0.0001; 
        $maxscore += log($maxp)/log(2);
    }
    return $minscore, $maxscore;
}
        

1;


