#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2011-06-16 23:02:40 -0700 (Thu, 16 Jun 2011) $';

our ($verbose, $help, $man, $output, $numwindow, $backgroundgc);

GetOptions ('verbose'=>\$verbose, 'help'=>\$help, 'man|m'=>\$man, 'output=s'=>\$output, 'numwindow=i'=>\$numwindow, 'backgroundgc=f'=>\$backgroundgc) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV>=2 or pod2usage ("Syntax error: <inputfile> missing");
our ($gcfile, $signalfile) = @ARGV;

$numwindow ||= 100;
$backgroundgc ||= 0.42;

if ($output) {
	open (STDOUT, ">$output") or confess "Error: cannot write to output file $output: $!";	#!
}

my ($header, $header_seg, $lrr_index, $pos_index, $chr_index, $name_index);

open (SIG, $signalfile) or confess "Error: cannot read from file $signalfile: $!";
$header = <SIG>;
$header =~ m/(.+)Pos/ or confess "error: the header file of signalfile $signalfile does not contain Pos annotation";
$header_seg = $1;
$pos_index = ($header_seg =~ tr/\t/\t/+0);
$header =~ m/(.+)Chr/ or confess "error: the header file of signalfile $signalfile does not contain Chr annotation";
$header_seg = $1;
$chr_index = ($header_seg =~ tr/\t/\t/+0);
$header =~ m/(.*)Name/ or confess "error: the header file of signalfile $signalfile does not contain Chr annotation";
$header_seg = $1;
$name_index = ($header_seg =~ tr/\t/\t/+0);

my %alldata;
my $count = 0;
my (%snp_chr, %snp_pos);
while (<SIG>) {
	chomp;
	my @record = split (/\t/, $_);
	my ($curname, $curchr, $curpos) = @record[$name_index, $chr_index, $pos_index];
	$curchr or next;				#sometimes chr is 0
	$curchr =~ m/^\w+$/ or next;			#sometimes chr is ---
	push @{$alldata{$curchr}}, [$curpos, $curname];
	$snp_chr{$curname} = $curchr;
	$snp_pos{$curname} = $curpos;
	$count++;
}
for my $key (keys %alldata) {
	@{$alldata{$key}} = sort {$a->[0] <=> $b->[0]} @{$alldata{$key}};
}
print STDERR "NOTICE: Finished reading chr and position information for $count markers in ${\(scalar keys %alldata)} chromosomes\n";


my (%snp_count, %snp_sum);
my ($prechr, $prestart, $preend, $windowcount, $windowsum) = qw/-1 -1 0 0 0/;

$count=0;
my (%chr_index, %seen_chr);
open (GC, $gcfile) or confess "Error: cannot read from file $gcfile: $!";
while (<GC>) {
	my @record = split (/\t/, $_);
	my ($curchr, $curstart, $curend, $curcount, $cursum) = @record[1..3, 11..12];
	
	#skip the random chromosomes!
	$curchr =~ m/random/ and next;
	#only handle autosomes
	$curchr =~ s/^chr//;
	
	#skip these annotations for 1Mb sliding windows (since we only care about 5kb window)
	m/1K.wib/ and next;
	my $curdata = $alldata{$curchr};
	$curdata or next;
		
	if ($curchr eq $prechr) {
		$curstart > $prestart or confess "Error: a record in chr$curchr has position $curstart, less then the previous position $prestart";
	} else {
		$seen_chr{$curchr} and confess "Error: chr$curchr occur multiple times in non-continuous segment of the <gcfile>: at chr$curchr:$curstart";
		$seen_chr{$curchr}++;
	}
	
	$chr_index{$curchr} ||= 0;
	for my $i ($chr_index{$curchr} .. @$curdata-1) {
		if ($curdata->[$i][0] < $curstart - $numwindow*5120) {
			$chr_index{$curchr} = $i;
			next;
		}
		if ($curdata->[$i][0] > $curend + $numwindow*5120) {
			last;
		}
		$snp_count{$curdata->[$i][1]} += $curcount;
		$snp_sum{$curdata->[$i][1]} += $cursum;
	}
	$count++;
	
	$prestart = $curstart;
	$prechr = $curchr;
}

print STDERR "NOTICE: Finish processing $count lines in GC file\n";

print "Name\tChr\tPosition\tGC\n";
for my $name (keys %snp_count) {
	my $gc = $backgroundgc;
	$snp_count{$name} and $gc = $snp_sum{$name}/$snp_count{$name};
	print "$name\t$snp_chr{$name}\t$snp_pos{$name}\t", sprintf("%.3f", $gc), "\n";
}



=head1 SYNOPSIS

 cal_gc_snp.pl [arguments] <gcfile> <snpfile>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --numwindow <int>		number of sliding window (default=100, or 500kb on each side)
 	    --backgroundgc <float>	backgroud GC frequency (default=0.42)
 	    --output <file>		write output to this file
 	
 Function: calculate GC content surrounding each marker within specified sliding 
 window, using the UCSC GC annotation file (for example, 
 http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/gc5Base.txt.gz for 
 human NCBI36 genome assembly) that is also sorted

 Example: cal_gc_snp.pl gc5Base.txt.sorted signalfile -output file.gcmodel

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--numwindow>

the number of non-overlapping sliding window on each side of the SNP.

=item B<--backgroundgc>

background GC level (for genomic regions without base information). By default 
it is 0.42 for human genome.

=item B<--output>

specify the output file name

=back

=head1 DESCRIPTION

This program is used to calculate GC content for genomic region surrounding SNPs. 

The <gcfile> is UCSC Genome Browser annotation file 
http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/gc5Base.txt.gz. Despite 
the file name, it actually contains GC content per 5120bp. (The 5-bp GC content 
file name and path is annotated in this file, though). If you use other genome 
build (like hg19, mm9, etc), you should download from the corresponding 
directory names.

After downloading the gcfile from UCSC genome browser and unzip the file, next 
use

	sort -k 2,2 -k 3,3n <input > output

to sort this file such that chromosome and positions are sorted. When using an 
unsorted <gcfile> with this progrma, a warning message will be issued.

A sample <gcfile> is below (first 10 lines):

	585     chr1    0       5120    chr1.0  5       1024    0       /gbdb/hg18/wib/gc5Base.wib      0       100     1024    59840   3942400
	585     chr1    5120    10240   chr1.1  5       1024    1024    /gbdb/hg18/wib/gc5Base.wib      0       100     1024    59900   3904400
	585     chr1    10240   15360   chr1.2  5       1024    2048    /gbdb/hg18/wib/gc5Base.wib      0       100     1024    55120   3411200
	585     chr1    15360   20480   chr1.3  5       1024    3072    /gbdb/hg18/wib/gc5Base.wib      0       100     1024    49900   3078800
	585     chr1    20480   25600   chr1.4  5       1024    4096    /gbdb/hg18/wib/gc5Base.wib      0       100     1024    47600   2682400

The <snpfile> is a tab-delimited file, and should contain at least three 
columns: name, chr, position. The first line of the file (so-called header line) 
should tell which columns corresond to these three fields. A sample <snpfile> is 
given below (but it has an extra column):

	Name    Chr     Position        PFB
	rs1000113       5       150220269       0.564615751221256
	rs1000115       9       112834321       0.565931333264192
	rs10001190      4       6335534 0.5668604380025
	rs10002186      4       38517993        0.57141752993563
	rs10002743      4       6327482 0.567557695424774

The output will be printed to STDOUT, but can be redirected to a file. 
Alternatively, the --output argument can be used to write the output to a file. 
A sample output file look like this:

	Name    Chr     Position        GC
	rs4961  4       2876505 48.6531211131841
	rs3025091       11      102219838       38.1080923507463
	rs3092963       3       46371942        44.4687694340796
	rs3825776       15      56534122        40.7894123134328
	rs17548390      6       3030512 45.3604050062189
	rs2276302       11      113355350       43.8200598569652

The GC column contains the GC percentage for a region surrounding the SNP. By 
default a 200kb region will be used (100kb each side of the SNP), but the 
--numwindow argument can be used to adjust the window size.

For questions, bugs or comments, please email kai@openbioinformatics.org.

=cut      
