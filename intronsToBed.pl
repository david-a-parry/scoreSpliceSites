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

my %opts = ();
GetOptions(
    \%opts,
    'g|gff=s',
    'b|biotype=s',
    't|tsl=i',
    'j|merge',
    'o|output=s',
    'h|?|help',
) or usage("Error getting options!");
usage() if $opts{h};
usage("-g/--gff option is required!\n") if not $opts{g};
usage("-o/--output argument is required.") if not $opts{o};

my $IN;
if ($opts{g} =~ /\.gz$/){
    open ($IN, "gzip -dc $opts{g} |") or die "Can't open $opts{g} via gzip: $!\n";
}else{
    open ($IN, "<", $opts{g}) or die "Can't open $opts{g} for reading: $!\n";
}

my %exons = ();
my %introns = ();
my %u12_introns = ();#records U12 intron IDs 
my %u2_introns = ();#records U2 intron IDs 
my %unknown_introns = ();#records IDs of introns of unknown type
my %u12_exons = ();#records exon IDs for exons with a neighbouring U12 intron
my %u2_exons  = ();#as above for U2
my %unknown_exons = ();#as above for unknown intron types

my ($U12_INTR, $U2_INTR, $UNKNOWN_INTR) = setupOutput("introns");
my ($U12_EX, $U2_EX, $UNKNOWN_EX) = setupOutput("exons");
my ($TMPIN, $tmp_intron) = tempfile( UNLINK => 1);
my ($TMPEX, $tmp_exon) = tempfile( UNLINK => 1);

processGff();

close $TMPIN;
close $TMPEX;

$tmp_exon   = sortRegions($tmp_exon, "exon");
outputExons();

$tmp_intron = sortRegions($tmp_intron, "intron");
outputIntrons();

close $U2_INTR;
close $U12_INTR;
close $UNKNOWN_INTR;

#################################################
sub setupOutput{
    my $type = shift;
    my $u12 = "$opts{o}_u12-$type.bed";
    my $u2  = "$opts{o}_u2-$type.bed";
    my $un  = "$opts{o}_unknown-$type.bed";
    open (my $U12, ">", "$u12") or die "Can't open $u12 for writing: $!\n";    
    open (my $U2,  ">", "$u2")  or die "Can't open $u2 for writing: $!\n";    
    open (my $UN,  ">", "$un")  or die "Can't open $un for writing: $!\n";    
    return ($U12, $U2, $UN);
}

#################################################
sub outputExons{
    open (my $EX, '<', $tmp_exon) or die 
     "Can't read exon sequence file '$tmp_exon': $!\n";
    my $n = 0; 
    my %prev_exon;
    while (my $line = <$EX>){
        chomp $line;
        my $FH;
        my ($chrom, $start, $end, $id, undef, $strand) = split("\t", $line); 
        my $coord = "$chrom:$start-$end";
        my ($class, $subclass, $score);
        if (exists $u12_exons{$coord}){
            $class = 'U12';
            ($subclass, $score) = split(/\|/,$u12_exons{$coord});
            $FH = $U12_EX;
        }elsif (exists $unknown_exons{$coord}){
            $class = 'UNKNOWN';
            ($subclass, $score) = split(/\|/,$unknown_exons{$coord});
            $FH = $UNKNOWN_EX;
        }elsif (exists $u2_exons{$coord}){
            $class = 'U2';
            ($subclass, $score) = split(/\|/,$u2_exons{$coord});
            $FH = $U2_EX;
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
            score    => $score,
            fh       => $FH,
        );
        if (not %prev_exon){
            %prev_exon = %this_exon;
            next;
        }
        my $written = checkAndWrite(\%prev_exon, \%this_exon, );
        if($written){
            $n++;
            if (not $n % 10000){
                my $time = strftime( "%H:%M:%S", localtime );
                print STDERR "[$time] Wrote scores for $n exons...\n" ;
            }
        }
    }
    writeBedLine
    (
        class    => $prev_exon{class},
        subclass => $prev_exon{subclass},
        chrom    => $prev_exon{chrom},
        start    => $prev_exon{start},
        end      => $prev_exon{end},
        strand   => $prev_exon{strand},
        id       => $prev_exon{id},
        score    => $prev_exon{score},
        fh       => $prev_exon{fh},
    ) if %prev_exon;
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished writing scores for $n exons.\n";
    close $EX;
}


#################################################
sub outputIntrons{
    open (my $IN, '<', $tmp_intron) or die 
     "Can't read intron sequence file '$tmp_intron': $!\n";
    my $n = 0;
    my %prev_intron;
    while (my $line = <$IN>){
        chomp $line;
        my $FH;
        my ($chrom, $start, $end, $id, $score, $strand) = split("\t", $line); 
        my $coord = "$chrom:$start-$end";
        my ($class, $subclass);
        if (exists $u12_introns{$coord}){
            $class = 'U12';
            $subclass = $u12_introns{$coord};
            $FH = $U12_INTR;
        }elsif (exists $unknown_introns{$coord}){
            $class = 'UNKNOWN';
            $subclass = $unknown_introns{$coord};
            $FH = $UNKNOWN_INTR;
        }elsif (exists $u2_introns{$coord}){
            $class = 'U2';
            $subclass = $u2_introns{$coord};
            $FH = $U2_INTR;
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
            score    => $score,
            fh       => $FH,
        );
        if (not %prev_intron){
            %prev_intron = %this_intron;
        }
        my $written = checkAndWrite(\%prev_intron, \%this_intron);
        if($written){
            $n++;
            if (not $n % 10000){
                my $time = strftime( "%H:%M:%S", localtime );
                print STDERR "[$time] Processed scores for $n introns...\n" ;
            }
        }
    }
    writeBedLine
    (
        class    => $prev_intron{class},
        subclass => $prev_intron{subclass},
        chrom    => $prev_intron{chrom},
        start    => $prev_intron{start},
        end      => $prev_intron{end},
        strand   => $prev_intron{strand},
        id       => $prev_intron{id},
        score    => $prev_intron{score},
        fh       => $prev_intron{fh},
    ) if (%prev_intron);
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] Finished writing scores for $n introns.\n";
    close $IN;
}

#################################################
sub checkAndWrite{
#read our hashes of intron features
#if they overlap and are of same type, merge
#otherwise print previous and store current 
    my $prev    = shift;
    my $current = shift;

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
        writeBedLine
        (
            class    => $prev->{class},
            subclass => $prev->{subclass},
            chrom    => $prev->{chrom},
            start    => $prev->{start},
            end      => $prev->{end},
            strand   => $prev->{strand},
            id       => $prev->{id},
            score    => $prev->{score},
            fh       => $prev->{fh},
        );
        %{$prev} = %{$current};
        return 1;
    }
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
        my ($chrom, $start, $end, $id, $score, $strand) = split("\t", $line); 
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
    my @introns = (); #collect all intron features per transcript for processing
    my %transcripts = (); 
    my %genes = (); 
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
            my ($tr) = $feat->get_tag_values('transcript_id'); 
    #        my ($parent) = $feat->get_tag_values('Parent');
            my $gene_type = join(",", $feat->get_tag_values('gene_type'));
            my $tsl;
            eval
            {
                ($tsl) = $feat->get_tag_values('transcript_support_level'); 
            };#not all transcripts have tsl tag - e.g. pseudogenes
            $tsl ||= 'NA';
            #clear cruft from tsl
            $tsl =~ s/\s+.*//;
            $transcripts{$tr}->{gene_type} = $gene_type;
            $transcripts{$tr}->{tsl} = $tsl;
        }elsif ($feat->primary_tag eq 'exon'){
            #collect exons
            my $do_write = 1;
            my ($parent) = $feat->get_tag_values('Parent');
            if (not exists $transcripts{$parent}){
                warn "Encountered exon without associated transcript! ";
                next;
            }
            if ($opts{b}){
                if (not grep {$_ eq $opts{b}} split
                (
                        ",",
                     $transcripts{$parent}->{gene_type}
                )
            ){
                    $do_write = 0;
                }
            }
            if ($opts{s}){
                if ($transcripts{$parent}->{tsl} eq 'NA'){
                    $do_write = 0;
                }elsif($transcripts{$parent}->{tsl} > $opts{s}){
                    $do_write = 0;
                }
            }
            if ($do_write){
                writeTempExonRegion($feat, 0);
            }
        }elsif ($feat->primary_tag eq 'intron'){
            #collect intron type related to each exon
            my $do_write = 1;
            my ($parent) = $feat->get_tag_values('Parent');
            if (not exists $transcripts{$parent}){
                warn "Encountered intron without associated transcript! ";
                next;
            }
            if ($opts{b}){
                if (not grep {$_ eq $opts{b}} split
                (
                        ",",
                     $transcripts{$parent}->{gene_type}
                )
            ){
                    $do_write = 0;
                }
            }
            if ($opts{s}){
                if ($transcripts{$parent}->{tsl} eq 'NA'){
                    $do_write = 0;
                }elsif($transcripts{$parent}->{tsl} > $opts{s}){
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
    my $score   = shift;

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
        $score,
        $strand,
    ) . "\n";
    $exons{$coord} = undef;
}

#            scores   => $prev->{score},
#################################################
sub writeTempIntronRegion{
    my $intron  = shift;
    my $next_ex = shift;
    my $prev_ex = shift;
    my $score   = shift;

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
        $score,
        $strand,
    ) . "\n";
    $introns{$coord} = undef;
}   

#################################################
sub writeBedLine{
    my %args = @_;
    my $chr = $args{chrom};
    my $strand = $args{strand};
    if ($strand ne '+' and $strand ne '-'){
        $strand = $strand > 0 ? '+' : '-';
    }
    $args{fh}->print
    (
        (join
            (
                "\t", 
                $args{chrom},
                $args{start} - 1,#bed format is 0-based
                $args{end},
                "$args{class}/$args{subclass}",
                "$args{score}",
                $strand,
            ) . "\n"
        )
    );
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
        my $score = 0;
        if ($type =~ /U12$/){
            $u12_exons{$next}     = "$type|$score";
            $u12_exons{$previous} = "$type|$score";
            $u12_introns{$intron_coords} = $type;
            ($score) = $intr->get_tag_values("donor_score_$type");
        }elsif($type =~ /U2$/){
            $u2_exons{$next}     = "$type|$score";
            $u2_exons{$previous} = "$type|$score";
            $u2_introns{$intron_coords} = $type;
            ($score) = $intr->get_tag_values("donor_score_$type");
        }else{
            $unknown_exons{$next}     = "$type|$score";
            $unknown_exons{$previous} = "$type|$score";
            $unknown_introns{$intron_coords} = $type;
        }
        writeTempIntronRegion($intr, $next_id, $previous_id, $score);
        $$count++;
        reportProgress($$count);
    }
    @$in = (); 
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


DESCRIPTION: Prints bed regions for introns and exons from a GFF file generated by spliceScorer.pl

USAGE: $0 -g introns.gff3 -o output_prefix

OPTIONS:
    
    -g,--gff FILE
        GFF3 intron file created by spliceScorer.pl

    -o,--output FILE
        Output file prefix. 
    
    -b,--biotypes STRING
        Only include exons/introns from transcripts of this biotype (e.g. 
        protein_coding)
   
    -j,--merge
        Use this flag to merge overlapping introns/exons before calculating
        stats.

    -t,--tsl INT
        Only include exons/introns from transcripts with this transcript 
        support level or higher. Transcripts with a tsl of 1 have the most 
        support, those with a tsl of 5 have the least. See Ensembl's help:
        http://www.ensembl.org/Help/Glossary?id=492.

    -h,--help 
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}
