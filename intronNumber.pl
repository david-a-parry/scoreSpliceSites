#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw( min max );
use POSIX qw/strftime/;
use Bio::Tools::GFF;
use FindBin qw($RealBin);
use Bio::DB::Sam;
use Bio::SeqFeature::Generic;

my %opts = (); 
GetOptions(
    \%opts,
    'g|gff=s',
    'o|output=s',
    'h|?|help',
) or usage("Error getting options!");
usage() if $opts{h};
usage("-g/--gff argument is required.") if not $opts{g};

my $IN;
if ($opts{g} =~ /\.gz$/){
    open ($IN, "gzip -dc $opts{g} |") or die "Can't open $opts{g} via gzip: $!\n";
}else{
    open ($IN, "<", $opts{g}) or die "Can't open $opts{g} for reading: $!\n";
}

my $gff = Bio::Tools::GFF->new
(
    -gff_version => 3,
    -fh          => $IN,
);
my %intron_counts = ();
my $OUT;
if ($opts{o}){
    open ($OUT, ">", $opts{o}) or die "Can't open $opts{o} for writing: $!\n";    
}else{
    $OUT = \*STDOUT;
}

my %introns = ();
my %transcripts = ();
my %names = ();
my $i = 0;
while (my $feat = $gff->next_feature() ) {
    #print STDERR "TAG: " . $feat->primary_tag . "\n";
    #if ($feat->primary_tag =~ /^(gene|processed_transcript|miRNA_gene|RNA)$/){
    if ($feat->primary_tag eq 'gene' or #works for all gencode
       ( $feat->primary_tag =~ /gene|RNA|processed_transcript/ 
        and $feat->has_tag('gene_id') ) #for ensembl GFFs
    ){
        parseIntrons(\%introns);
        %introns = ();
        my ($id) = $feat->get_tag_values('ID'); 
        $id =~ s/^gene://;
        my $name = '.';
        if ($feat->has_tag('Name')){
            ($name) = $feat->get_tag_values('Name'); 
        }elsif ($feat->has_tag('gene_name')){
            ($name) = $feat->get_tag_values('gene_name'); 
        }
        $names{$id} = $name; 
    }elsif ($feat->primary_tag eq 'transcript' or #works for gencode
           ( $feat->primary_tag =~ /gene|transcript|RNA/ 
             and $feat->has_tag('transcript_id') ) #for ensembl GFFs
    ){
        my ($tr) = $feat->get_tag_values('ID'); 
        $tr =~ s/transcript://;
        ($transcripts{$tr}->{parent}) = $feat->get_tag_values('Parent');
        $transcripts{$tr}->{length} = 1 + $feat->end - $feat->start;
        if ($feat->has_tag('gene_type')){
            $transcripts{$tr}->{gene_type} = join(",", $feat->get_tag_values('gene_type'));
        }elsif ($feat->has_tag('biotype')){
            $transcripts{$tr}->{gene_type} = join(",", $feat->get_tag_values('biotype'));
        }else{
            die "Could not determine biotype/gene_type for $tr!\n";
        }

        eval{
            ($transcripts{$tr}->{tsl}) = $feat->get_tag_values('transcript_support_level')
        };#not all transcripts have tsl tag(?)
        $transcripts{$tr}->{tsl} ||= 'NA';
        #clear cruft from tsl
        $transcripts{$tr}->{tsl} =~ s/\s+.*//;

        parseIntrons(\%introns);
        %introns = ();
        #foreach my $tr ($feat->get_tag_values('ID') ){
        #    $transcripts{$tr} = $feat;
        #}
    }elsif ($feat->primary_tag eq 'intron'){
        foreach my $tr ($feat->get_tag_values('Parent') ){
            $tr =~ s/transcript://;
            my $ex_tag = '';
            if ($feat->strand > 0){
                $ex_tag = 'previous_';
            }else{
                $ex_tag = 'next_';
            }
            if ($feat->has_tag($ex_tag . 'exon_number') ){
                $ex_tag = $ex_tag . 'exon_number';#gencode GFFs
            }elsif ($feat->has_tag($ex_tag . 'rank') ){
                $ex_tag = $ex_tag . 'rank';#ensembl GFFs
            }else{
                die "Could not determine exon number for ". $feat->{ID} . "\n";
            }
            foreach my $ex ($feat->get_tag_values($ex_tag) ){
                $introns{$tr}->{$ex} = $feat;
                foreach my $t ($feat->get_tag_values('intron_type') ){
                    ($transcripts{$tr}->{intron_type}->{$ex}) = $t;
                }
            }
        }
#debug    }elsif ($feat->primary_tag ne 'CDS' and $feat->primary_tag !~ /UTR/){
#debug        print Dumper $feat;
    }
}
$gff->close();
parseIntrons(\%introns);
my $time = strftime( "%H:%M:%S", localtime );
$i =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g; #add commas for readability
print STDERR "[$time] processed $i introns\nFinished\n";

#################################################
sub parseIntrons{
    my $introns = shift;
    my $n = 0;
    foreach my $tr (keys %$introns){
        my %u12 = (); 
        my $last_intr;
        foreach my $intr (sort {$a <=> $b} keys %{$introns->{$tr}}){
            $i++;
            (my $ex_start)  = $introns{$tr}->{$intr}->get_tag_values('previous_exon_start');
            (my $ex_end)    = $introns{$tr}->{$intr}->get_tag_values('previous_exon_end');
            ($ex_start, $ex_end) = sort { $a <=> $b } ($ex_start, $ex_end);
            $transcripts{$tr}->{spliced_length} += 1 + $ex_end - $ex_start;
            $transcripts{$tr}->{introns}++;  
            if ($transcripts{$tr}->{intron_type}->{$intr} =~ /U12$/){
                $u12{$intr}->{length} = $transcripts{$tr}->{spliced_length};
                $u12{$intr}->{strand} = $introns{$tr}->{$intr}->strand;
            }
            reportProgress($i);
            $last_intr = $intr;
        }
        (my $ex_start)  = $introns{$tr}->{$last_intr}->get_tag_values('next_exon_start');
        (my $ex_end)    = $introns{$tr}->{$last_intr}->get_tag_values('next_exon_end');
        ($ex_start, $ex_end) = sort { $a <=> $b } ($ex_start, $ex_end);
        $transcripts{$tr}->{spliced_length} += 1 + $ex_end - $ex_start;
        foreach my $intr(keys %u12){
            my $c_pos = $u12{$intr}->{length};
            my $in = $intr;
            if ($u12{$intr}->{strand} < 1 ){
                $c_pos = 1 + $transcripts{$tr}->{spliced_length} - $u12{$intr}->{length};
                $n = $transcripts{$tr}->{introns} - $intr + 1;
            }
            my $name = $names{$transcripts{$tr}->{parent}};
            $name ||= $transcripts{$tr}->{parent};
            print $OUT join("\t", 
                 $tr,
                 $transcripts{$tr}->{parent},
                 $name,
                 $transcripts{$tr}->{gene_type},
                 $u12{$intr}->{strand},
                 $transcripts{$tr}->{tsl},
                 $u12{$intr}->{length},
                 $transcripts{$tr}->{spliced_length},
                 $in,
                 $transcripts{$tr}->{introns},
                 $transcripts{$tr}->{introns} - $in + 1,
            ) . "\n";
        } 
    }
}

#################################################
sub reportProgress{
    my $n = shift;
    return if not $n;
    return if ($n % 10000); 
    my $time = strftime( "%H:%M:%S", localtime );
    $n =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g; #add commas for readability
    printf STDERR
    (
        "[$time] processed $n introns\n"
    );
}

#################################################
sub usage{
    my $msg = shift;
    print "\n$msg\n" if $msg;

    print <<EOT

Print intron number and CDS position of U12 introns

USAGE: $0 -g genes_and_introns.gff3

OPTIONS:
    
    -g,--gff FILE
        GFF3 file containing information on the genes and introns to use for intron file creation

    -o,--output FILE
        Optional output file. Default = STDOUT.

    -h,--help 
        Show this message and exit

EOT
;
    exit 1 if $msg;
    exit;
}




