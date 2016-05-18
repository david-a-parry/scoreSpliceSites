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
use File::Basename;
use Sort::External;
use Bio::DB::Sam;
use Bio::SeqFeature::Generic;
#use re 'debugcolor';
use lib "$RealBin/lib";
use ReverseComplement qw/ reverse_complement /;

my @repeat_units = ();  #if we only want to match specific repeats put
                        #the repeat units we want to match in this array
my $MAX_REP_SCAN = 100; #in the first instance only look for a 
                        # maximum of this many repeats in regex to prevent
                        # deep recursion
my %opts = 
(
    m => 3, 
    x => 9999999999,
    u => 1,
    y => 6,
    n => 2,
    t => 3,
    r => \@repeat_units, 
);
GetOptions(
    \%opts,
    'f|fasta=s',
    'g|gff=s',
    'b|biotype=s',
    's|tsl=i',
    'merge|j',
    'u|min_repeat_unit_length=i',
    'y|max_repeat_unit_length=i',
    'm|min_repeat_length=i',
    'x|max_repeat_length=i',
    'n|min_number_of_repeats=i',
    't|trim_exons=i',
    'r|repeat_units=s{,}',
    'o|output=s',
    'e|exon_seqs=s',
    'i|intron_seqs=s',
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
my $fai = Bio::DB::Sam::Fai->load($opts{f});#should create index if it doesn't exist

my ($I_OUT, $E_OUT) = setupOutput();
my ($TMPEX, $ex_seq_file) = setupExonSeqFile();
my ($TMPIN, $in_seq_file) = setupIntronSeqFile();

my %exon_seqs = ();#all values undef, records ids of exons we've already retrieved seq for
my %intron_seqs = ();#all values undef, records ids of exons we've already retrieved seq for
my %u12_introns = ();#records U12 intron IDs 
my %u2_introns = ();#records U2 intron IDs 
my %unknown_introns = ();#records IDs of introns of unknown type
my %u12_exons = ();#records exon IDs for exons with a neighbouring U12 intron
my %u2_exons  = ();#as above for U2
my %unknown_exons = ();#as above for unknown intron types
#my %names = ();

processIntronGff();

$opts{e} = sortSeqs($TMPEX, $ex_seq_file, $opts{e}, "exon");
$opts{i} = sortSeqs($TMPIN, $in_seq_file, $opts{i}, "intron");
writeHeaders();
outputExonStats();
outputIntronStats();


#################################################
sub outputIntronStats{
    open (my $IN, '<', $opts{i}) or die 
     "Can't read intron sequence file '$opts{i}': $!\n";
    my $n = 0;
    my %prev_intron;
    while (my $line = <$IN>){
        my ($coord, $id, $seq) = split("\t", $line); 
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
        my ($chrom, $start, $end) = split(/[:-]/, $coord);
        my %this_intron = 
        (
            class    => $class,
            subclass => $subclass,
            chrom    => $chrom,
            id       => $id,
            start    => $start,
            end      => $end,
            seq      => $seq,
        );
        if (not %prev_intron){
            %prev_intron = %this_intron;
        }
        my $written = checkAndWrite(\%prev_intron, \%this_intron, $I_OUT);
        if($written){
            $n++;
            if (not $n % 10000){
                my $time = strftime( "%H:%M:%S", localtime );
                print STDERR "[$time] Processed stats for $n introns...\n" ;
            }
        }
    }
    writeExonStats
    (
        $prev_intron{class},
        $prev_intron{subclass},
        "$prev_intron{chrom}:$prev_intron{start}-$prev_intron{end}",
        $prev_intron{id},
        $prev_intron{seq},
        $I_OUT,
    );
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished writing stats for $n introns.\n";
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
    if ($opts{merge} and 
        $prev->{chrom} eq $current->{chrom} and 
        $prev->{end} >= $current->{start} and #should be sorted in coordinate order
        $prev->{subclass} eq $current->{subclass}
    ){#overlaps
        if ($current->{end} > $prev->{end}){
            $prev->{id} = $prev->{id} . "/" . $current->{id};
            $prev->{end} = $current->{end};
            if ($current->{start} == $prev->{start}){
                $prev->{seq} = $current->{seq};
            }else{
                $prev->{seq} = $fai->fetch("$prev->{chrom}:$prev->{start}-$prev->{end}");
            }
        }
        return 0;
    }else{#doesn't overlap
        writeExonStats
        (
            $prev->{class},
            $prev->{subclass},
            "$prev->{chrom}:$prev->{start}-$prev->{end}",
            $prev->{id},
            $prev->{seq},
            $FH,
        );
        %{$prev} = %{$current};
        return 1;
    }
}

#################################################
sub outputExonStats{
    open (my $EX, '<', $opts{e}) or die 
     "Can't read exon sequence file '$opts{e}': $!\n";
    my $n = 0; 
    my %prev_exon;
    while (my $line = <$EX>){
        my ($coord, $id, $seq) = split("\t", $line); 
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
        my ($chrom, $start, $end) = split(/[:-]/, $coord);
        my %this_exon = 
        (
            class    => $class,
            subclass => $subclass,
            chrom    => $chrom,
            id       => $id,
            start    => $start,
            end      => $end,
            seq      => $seq,
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
                print STDERR "[$time] Wrote stats for $n exons...\n" ;
            }
        }
    }
    writeExonStats
    (
        $prev_exon{class},
        $prev_exon{subclass},
        "$prev_exon{chrom}:$prev_exon{start}-$prev_exon{end}",
        $prev_exon{id},
        $prev_exon{seq},
        $E_OUT,
    );
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished writing stats for $n exons.\n";
    close $EX;
    close $E_OUT;
}

#################################################
sub sortSeqs{
    my ($TMP_FH, $tmp_seq_file, $outfile, $type) = @_;
    return $outfile if not $TMP_FH;
    close $TMP_FH;
    open (my $EXIN, '<', $tmp_seq_file) or die 
     "Can't open temporary $type sequence file '$tmp_seq_file' for reading: $!\n";
    my $sortex = Sort::External->new(mem_threshold => 1024**2 * 16);
    #exon IDs should all be same length so no need for special sort sub
    my @feeds = ();
    my $n = 0;
    my $SORTOUT;
    my $time;
    if ($outfile){
        open ($SORTOUT, '>', $outfile) or die 
         "Can't open sequence output file '$outfile' for writing: $!\n";
    }else{
        ($SORTOUT, $outfile) = tempfile("ex_seqs_XXXX", UNLINK => 1);
    }
    while (my $line = <$EXIN>){
        next if $line =~ /^#/;
        my @split = split("\t", $line); 
        my ($chrom, $start, $end) = split(/[:-]/, $split[0]);
        my $packstart = pack("N", $start); 
        my $packend = pack("N", $end); 
        push @feeds, "$chrom,$packstart,$packend\t$line";
        $n++;
        if (@feeds > 9999){
            $sortex->feed(@feeds);
            @feeds = ();
            $time = strftime( "%H:%M:%S", localtime );
            print STDERR "[$time] Fed $n $type"."s to sort...\n";
        }
    }
    $sortex->feed(@feeds) if @feeds;
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished feeding $n $type"."s to sort...\n";
    @feeds = ();
    $sortex->finish; 
    $n = 0;
    print $SORTOUT <<EOT
#sorted $type sequence file
#COORDS\tID\tSEQ
EOT
;
    $n = 0;
    while ( defined( $_ = $sortex->fetch ) ) {
        my @split = split("\t", $_);
        print $SORTOUT join("\t", @split[1..$#split]);
        $n++;
        if (not $n % 100000){
            $time = strftime( "%H:%M:%S", localtime );
            print STDERR "[$time] Wrote $n $type"."s to sequence file...\n";
        }
    }
    $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished writing $n $type"."s to sequence file...\n";
    close $SORTOUT;
    return $outfile;
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
            ($TMP, $ex) = tempfile(UNLINK => 1);
        }else{
            return (undef, undef);
        }
    }else{
        ($TMP, $ex) = tempfile(UNLINK => 1);
    }
    print $TMP <<EOT
#unsorted temporary exon sequence file
#COORDs\tEXON_ID\tSEQ
EOT
;
    return ($TMP, $ex);
}

#################################################
sub setupIntronSeqFile{
    my $in;
    my $TMP;
    if ($opts{i}){
    #if user supplied --exon_seqs file see it it already exists 
    # - will read from it instead of retrieving seqs if it does
    # we write to a tmp file even if we will be creating a new --exon_seqs
    # file cos we need to sort the output once we've got all the seqs
        $in = $opts{i};
        if (not -e $in){
            ($TMP, $in) = tempfile( UNLINK => 1);
        }else{
            return (undef, undef);
        }
    }else{
        ($TMP, $in) = tempfile( UNLINK => 1);
    }
    print $TMP <<EOT
#unsorted temporary intron sequence file
#COORDS\tINTRON_ID\tSEQ
EOT
;
    return ($TMP, $in);
}

#################################################
sub writeHeaders{
    my @opt_string = ();
    foreach my $k ( sort keys %opts ) {
        if ( not ref $opts{$k} ) {
            push @opt_string, "$k=$opts{$k}";
        }elsif ( ref $opts{$k} eq 'SCALAR' ) {
            if ( defined ${ $opts{$k} } ) {
                push @opt_string, "$k=${$opts{$k}}";
            }else {
                push @opt_string, "$k=undef";
            }
        }elsif ( ref $opts{$k} eq 'ARRAY' ) {
            if ( @{ $opts{$k} } ) {
                push @opt_string, "$k=" . join( ",", @{ $opts{$k} } );
            }else {
                push @opt_string, "$k=undef";
            }
        }
    }
    my $caller = fileparse($0);
    foreach my $FH ($E_OUT, $I_OUT){
        print $FH "##$caller\"" . join( " ", @opt_string ) . "\"\n";
        my @cols =  qw(
            INTRON_TYPE
            SUBTYPE
            COORDS
            ID
            LENGTH
            PERCENT_GC
            REPEATS
            TOTAL_REPEAT_LENGTH
            LONGEST_REPEAT
            MEAN_REPEAT_LENGTH
        );
        for my $l ($opts{u}..$opts{y}){
            for my $s (qw/
                _NT_REPEATS
                _NT_REPEATs_TOTAL_LENGTH
                _NT_REPEATS_LONGEST
                _NT_REPEATS_MEAN_LENGTH
            /){
                push @cols, "$l$s";
            }
        }
        print $FH join
        (   
            "\t",
            @cols,
        ) . "\n";
    }
}

#################################################
sub setupOutput{
    my $ex = "$opts{o}_exon_stats.tsv";
    my $in = "$opts{o}_intron_stats.tsv";
    open (my $EX, ">", "$ex") or die "Can't open $ex for writing: $!\n";    
    open (my $IN, ">", "$in") or die "Can't open $in for writing: $!\n";    
    return ($IN, $EX);
}

#################################################
sub processIntronGff{
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
            parseIntrons(\@introns, \$intron_count, $fai);
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
            parseIntrons(\@introns, \$intron_count, $fai);
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
                writeTempExonSequence($feat, $fai);
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
    parseIntrons(\@introns, \$intron_count, $fai);
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
    my $fai = shift;
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
        writeTempIntronSequence($intr, $fai, $next_id, $previous_id);
        $$count++;
        reportProgress($$count);
    }
    @$in = (); 
}


#################################################
sub sortRepArray{
    my $hashes = shift;
    @$hashes = sort 
    { 
        $a->{start} <=> $b->{start} ||
        $a->{end} <=> $b->{end}
    } @$hashes;
}
#################################################
sub getRepeatStats{
#get length of each repeat in each category
#and calculate longest, cumulative lengths and mean lengths
    my $r = shift; 
    my @stats = (); 
    foreach my $k ( "all", $opts{u}..$opts{y},){
        if (exists $r->{$k}){
            if ($k eq 'all'){
                sortRepArray($r->{$k});#all other array should be sorted already
            }
            my @l = ();
            my $overlaps = 0; 
            my $prev_hash; #for overlaps
            foreach my $rep_hash (@{$r->{$k}}){
                push @l, length($rep_hash->{seq});
                #calculate overlaps (if any)
                if ($prev_hash){
                    if ($rep_hash->{start} <= $prev_hash->{end}){
                        $overlaps += $prev_hash->{end} - $rep_hash->{start} + 1;
                    }
                }
                $prev_hash = $rep_hash;
                
            }
            my $longest = max(@l);
            my $total = sum(@l);
            #if $r->{$k} exists we must have at least one array entry, no need 
            # to worry about divide by 0 error
            my $mean = sprintf("%.3f", $total/@l);
            #some repeats may overlap so we need to subtract overlaps from total
            # but only after calculating mean repeat length
            $total -= $overlaps;
            push @stats, scalar(@l), $total, $longest, $mean;
        }else{
            push @stats, 0, 0, 0, 0;
        }
    }
    return @stats;
}
        
#################################################
sub writeExonStats{
    my ($class, $subclass, $coords, $id, $seq, $FH) = @_;
    my $trimmed;
    if ($opts{t} > 0){
        my $trimmed_length = length($seq) - $opts{t} - $opts{t};
        return if $trimmed_length < 1;
        $trimmed = substr($seq, $opts{t}, $trimmed_length - 1 );
    }else{
        $trimmed = $seq;
    }
    my $gc = getGcPercentage($seq);
    my %reps = getRepeats($seq);
    my @repeat_stats = getRepeatStats(\%reps); 
    print $FH join
    (
        "\t",
        $class,
        $subclass,
        $coords,
        $id,
        length($seq),
        sprintf("%.3f", $gc),
        @repeat_stats,
    ) . "\n";
    #intron type (U2 or U12 or UNKNOWN), subtype, exon ID, exon length, % GC,
    # no. repeats, total length of repeats, longest repeat, mean repeat length,
    # no.homopolymers, total length of homopolymers, longest homopolymer, mean homopolymer length
}

#################################################
sub getRepeats{
# return array of 1-6 (or rather $opts{u}-$opts{y}) nucleotide repeats 
# between $opts{m} and $opts{x} nt long found in $seq
# repeats must consist of at least $opts{n} repeat units (e.g. A * $opts{n} or CAG * $opts{n})
    my $seq = shift;
    my $min_n_rep = $opts{n} - 1;
    if ($MAX_REP_SCAN < $min_n_rep){
        $MAX_REP_SCAN = $min_n_rep + 1;
    }
=cut
    while ($seq =~ /((\w{$opts{u},$opts{y}})(\2){$min_n_rep,})/g){#does not always get longest full match cos it will find the longest possible repeat in $2
        my $rep = $1;
        my $l   = length($rep);
        if ($l >= $opts{m} and $l <= $opts{x}){
            push @{$r{all}}, $rep;
            push @{$r{$l}}, $rep;
        }
    }
=cut

    #test all repeats between $opts{u} and $opts{y} nt 
    # (ensuring we reduce to the simplest repeat unit)
    # record all repeats even if overlapping, but calculate total repeat length based on 
    # non-overlapping length of repeats
    my %r = (); #hash of arrays of hashes for each repeat found
                   #key is repeat unit length
                   #value is array of hashes with keys:   
                   #    start => start coordinate of repeat 
                   #    end   => end coordinate of repeat
                   #    seq   => repeat sequence
    #for each possible repeat unit length
    for my $n ($opts{u} .. $opts{y}){
        my $longest = ''; 
        my $rep_length = 0;
        for (my $i = 0; $i < length($seq); ){
            my $s = substr($seq, $i, $n * $MAX_REP_SCAN); 
            if ($s =~ /^((\w{$n})(\2){$min_n_rep,$MAX_REP_SCAN})/){
            #$MAX_REP_SCAN is there to prevent deep recursion in complex regex
                my $rep = $1;
                #make sure we don't mistake a 2nt repeat for a 4nt repeat etc.
                my $repeat_length = getRepeatLength($rep);
                my $r = substr($rep, 0, $repeat_length);#actual repeat unit
                if ($repeat_length != $n){
                    $i++;
                    next;
                }elsif(@repeat_units){#if we've specified the repeat units we're looking at, check if it matches
                    if (not grep {$r eq $_} @repeat_units){
                        $i++;
                        next;
                    }
                }
                #we can get full repeat with a simpler regex if we have fewer 
                # than $MAX_REP_SCAN repeats 
                $s = substr($seq, $i,);
                if ($s =~ /^(($r){$min_n_rep,})/){
                    $rep = $1;
                }else{
                    warn "INTERNAL ERROR FINDING FULL LENGTH OF REPEAT for $rep!";
                }
                #check total length of repeat is within our min/max limits
                my $l = length($rep); 
                if ($l >= $opts{m} and $l <= $opts{x}){
                    #add to arrays of hashes for this repeat length
                    my $h = 
                    {
                        start => $i, 
                        end   => $i + $l, 
                        seq   => $rep,
                    };
                    push @{$r{$n}}, $h;
                    push @{$r{all}}, $h;
                }
                #we've found a repeat, next different repeat could overlap
                # by $n - 1 nucleotides with end of this repeat
                # - increment $i accordingly
                $i += ($l - ($n - 1));
            }else{
                $i++;#no repeat found, increment by 1
            }
        }
    }
    return %r;
}

###########################################################
sub getRepeatLength{
#finds shortest repeating unit of string
#if string is not repetitive returns length of string
    my $string = shift;
    my $l = length($string);
    #get first half of divisors (probable needless optimization, could just do this in one pass)
    my @div = grep{ $l % $_ == 0 } 1 .. sqrt($l);
    #get upper half of divisors 
    push @div, map {$l == $_*$_ ? () : $l/$_} reverse @div;
    for (my $i = 0; $i < @div; $i++){
        my $s = substr($string, 0, $div[$i]);
        my $j = $l/$div[$i];
        if ( $string eq ($s x $j)){
            return $div[$i];
        }
    }
    return $l;
}

###########################################################
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
sub writeTempIntronSequence{
    return if not $TMPIN;
    my $intron = shift;
    my $fai = shift;
    my $next_ex = shift;
    my $prev_ex = shift;
    my $id = "$prev_ex-$next_ex";
    my $chrom = $intron->seq_id; 
    my $strand = $intron->strand;
    my $start = $intron->start;
    my $end = $intron->end;
    my $coord = "$chrom:$start-$end";
    return if exists $intron_seqs{$coord};
    my $seq = $fai->fetch($coord);
    if ($strand < 0){
        $seq = reverse_complement($seq);
    }
    print $TMPIN "$coord\t$id\t$seq\n";
    $intron_seqs{$coord} = undef;
}

#################################################
sub writeTempExonSequence{
    return if not $TMPEX;
    my $exon = shift;
    my $fai = shift;
    my ($id) = $exon->get_tag_values('exon_id');
    my $chrom = $exon->seq_id; 
    my $strand = $exon->strand;
    my $start = $exon->start;
    my $end = $exon->end;
    my $coord = "$chrom:$start-$end";
    return if exists $exon_seqs{$coord};
    my $seq = $fai->fetch($coord);
    if ($strand < 0){
        $seq = reverse_complement($seq);
    }
    print $TMPEX "$coord\t$id\t$seq\n";
    $exon_seqs{$coord} = undef;
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
    usage("-o/--output argument is required.") if not $opts{o};
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
    if ($opts{n} < 2){
        usage("-n/--min_number_of_repeats argument must be greater than 1.");
    }
    my @simp_rep = (); 
    foreach my $r (@repeat_units){
        my $l = getRepeatLength($r);
        if (length($l) < length($r) ){
            my $s = substr($r, 0, $l); 
            warn "Simplifying '$r' to '$s' repeat.\n";
            push @simp_rep, $s;
        }else{
            push @simp_rep, $r;
        }
    }
    @repeat_units = @simp_rep;
#TODO - 
#check --min_number_of_repeats * repeat_unit_lengths not greater than --max_repeat_length?
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


DESCRIPTION: Prints stats (length, %GC and repeat information) for introns and 
exons from a GFF file generated by spliceScorer.pl

USAGE: $0 -f genome_fasta.fa -g introns.gff3

OPTIONS:
    
    -f,--fasta FILE
        Genome fasta file for retrieving DNA sequences for intron-exon 
        boundaries. Must not be repeat masked (soft or hard). 

    -g,--gff FILE
        GFF3 intron file created by spliceScorer.pl

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

    -m,--min_repeat_length INT
        Only count repeats at least this long (default = 3)
    
    -x,--max_repeat_length INT
        Only count repeats of this value or shorter 
        (default = 9999999999) 
    
    -u,--min_repeat_unit_length INT
        Only count repetetive sequences where the repeated unit is at least 
        this long (default = 1)
    
    -y,--max_repeat_unit_length INT
        Only count repetetive sequences where the repeated unit is this long or
        shorter (default = 6)
    
    -n,--min_number_of_repeats INT 
        Minimum number of repeat units to consider (default = 2).
    
    -r,--repeat_units STRING [STRING2 STRING3 ... STRINGn]
        One or more specific repeats to look for

    -o,--output FILE
        Output file prefix. Two output files will be produced - one with the 
        suffix 'intron_stats.tsv' the other with the suffix 'exon_stats.tsv'.
    
    -e,--exon_seqs FILE
        Optional file for reading/writing DNA sequences of exons retrieved. 
        If this file exists it will be read for assessing exon stats; if not 
        it will be created for future use.

    -i,--intron_seqs FILE
        Optional file for reading/writing DNA sequences of introns retrieved. 
        If this file exists it will be read for assessing intron stats; if not 
        it will be created for future use.

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


