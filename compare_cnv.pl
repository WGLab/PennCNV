#!/usr/bin/env perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2011-05-03 11:31:25 -0700 (Tue, 03 May 2011) $';

our ($verbose, $help, $man);
our ($operation, $cnvfile);
our ($logfile, $outfile, $listfile, $minoverlap, $reciprocal, $minmarkeroverlap, $pfbfile, $cnv2file);
my ($snp_chr, $snp_index);

GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'logfile=s'=>\$logfile, 'outfile=s'=>\$outfile, 'listfile=s'=>\$listfile,
	'minoverlap=f'=>\$minoverlap, 'reciprocal'=>\$reciprocal, 'minmarkeroverlap=f'=>\$minmarkeroverlap, 'pfbfile=s'=>\$pfbfile, 'cnv2file=s'=>\$cnv2file) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error: <inputfile> missing");

($operation, $cnvfile) = @ARGV;
$minoverlap ||= 0;
$minoverlap >=0 and $minoverlap <= 1 or pod2usage ("Error in argument: the --minoverlap argument should be between 0 and 1");

if ($minmarkeroverlap) {
	$minmarkeroverlap >=0 and $minmarkeroverlap <= 1 or pod2usage ("Error in argument: the --minmarkeroverlap argument should be between 0 and 1");
	defined $pfbfile or pod2usage ("Error in argument: please specify --pfbfile when --minmarkeroverlap is specified");
	($snp_chr, $snp_index) = readPFB ($pfbfile);
}

if ($logfile) {
	print STDERR "NOTICE: The program log message is written to $logfile\n";
	open (STDERR, ">$logfile") or confess "Error: cannot write to log file $logfile: $!\n";
}
if ($outfile) {
	print STDERR "NOTICE: The program output is written to $outfile\n";
	open (STDOUT, ">$outfile") or confess "Error: cannot write to outfile $outfile: $!\n";
}

if ($operation eq 'dup' or $operation eq 'compdup') {
	$listfile or pod2usage ("Error: please specify the --listfile argument for the $operation operation");
	print STDERR "NOTICE: <operation>='dup/compdup': find concordant CNV calls between duplicated samples (same sample but different file names)\n";
	compareDuplicatedSample ($cnvfile, $listfile, $snp_chr, $snp_index);
} elsif ($operation eq 'compcall') {
	$listfile or pod2usage ("Error in argument: please specify the --listfile arguments for the $operation operation");
	$cnv2file or pod2usage ("Error in argument: please specify the --cnv2file arguments for the $operation operation");
	print STDERR "NOTICE: <operation>='compcall': compare CNV calls from two files on specified pairs of samples\n";
	compareCNVCall ($cnvfile, $cnv2file, $listfile, $snp_chr, $snp_index);
} else {
	pod2usage ("Syntax Error: unknown <operation> is specified: $operation (valid operations include dup/compdup, compcall");
}

sub compareDuplicatedSample {
	my ($cnvfile, $listfile, $snp_chr, $snp_index) = @_;
	
	open (LIST, $listfile) or confess "Error: cannot read from list file $listfile: $!";

	print "subject1\tsubject2\tnumcnv1\tnumcnv2\toverlap1\toverlap2\n";

	my (@pair, %ind, %cnv);
	while (<LIST>) {
		s/[\r\n]+$//;
		my @ind = split (/\t/, $_);
		@ind == 2 or confess "Error: invalid record found in listfile $listfile (2 tab-delimited fields expected): <$_>\n";
		push @pair, [@ind];
		$ind{$ind[0]}++;
		$ind{$ind[1]}++;
	}
	close (LIST);
	
	if ($cnvfile eq 'stdin') {
		*CNV = *STDIN;
	} else {
		open (CNV, $cnvfile) or confess "Error: cannot read from cnvfile $cnvfile: $!";
	}
	while (<CNV>) {
		s/[\r\n]+$//;
		m/^chr((\w+):(\d+)-(\d+))\s+numsnp=(\d+)\s+length=([\d\,]+)\s+state\d+,cn=(\d+)\s+(.+?)\s+startsnp=(\S+)\s+endsnp=(\S+)/ or confess "Error: invalid record found in $cnvfile: <$_>";
		if ($ind{$8}) {
			push @{$cnv{$8}}, $_;
		}
	}
	close (CNV);
	
	for my $i (0 .. @pair-1) {
		my ($ind1, $ind2) = @{$pair[$i]};
		$cnv{$ind1} ||= [];
		$cnv{$ind2} ||= [];
		print $ind1, "\t", $ind2, "\t", scalar (@{$cnv{$ind1}}), "\t", scalar (@{$cnv{$ind2}}), "\t";
		print STDOUT scalar (cnvOverlap ($cnv{$ind1}, $cnv{$ind2}, $snp_chr, $snp_index)), "\t", scalar (cnvOverlap ($cnv{$ind2}, $cnv{$ind1}, $snp_chr, $snp_index)), "\n";
	}
}


sub compareCNVCall {
	my ($cnvfile, $cnv2file, $listfile, $snp_chr, $snp_index) = @_;

	my (@pair, %ind1, %ind2, %cnv1, %cnv2);
	open (LIST, $listfile) or confess "Error: cannot read from list file $listfile: $!";
	while (<LIST>) {
		s/[\r\n]+$//;
		my @ind = split (/\t/, $_);
		@ind == 2 or confess "Error: invalid record found in listfile $listfile (2 tab-delimited fields expected): <$_>\n";
		push @pair, [@ind];
		$ind1{$ind[0]}++;	#indicate that this ID corresponds to a subject in the first CNV call file
		$ind2{$ind[1]}++;	#indicate that this ID corresponds to a subject in the second CNV call file
	}
	close (LIST);
	
	if ($cnvfile eq 'stdin') {
		*CNV = *STDIN;
	} else {
		open (CNV, $cnvfile) or confess "Error: cannot read from cnvfile $cnvfile: $!";
	}
	while (<CNV>) {
		s/[\r\n]+$//;
		m/^chr((\w+):(\d+)-(\d+))\s+numsnp=(\d+)\s+length=([\d\,]+)\s+state\d+,cn=(\d+)\s+(.+?)\s+startsnp=(\S+)\s+endsnp=(\S+)/ or confess "Error: invalid record found in $cnvfile: <$_>";
		if ($ind1{$8}) {
			push @{$cnv1{$8}}, $_;
		}
	}
	close (CNV);
	if ($cnv2file eq 'stdin') {
		*CNV = *STDIN;
	} else {
		open (CNV, $cnv2file) or confess "Error: cannot read from cnv2file $cnv2file: $!";
	}
	while (<CNV>) {
		s/[\r\n]+$//;
		m/^chr((\w+):([\d]+)-([\d]+))\s+numsnp=(\d+)\s+length=([\d\,]+)\s+state\d+,cn=(\d+)\s+(.+?)\s+startsnp=(\S+)\s+endsnp=(\S+)/ or confess "Error: invalid record found in $cnvfile: <$_>";
		if ($ind2{$8}) {
			push @{$cnv2{$8}}, $_;
		}
	}
	close (CNV);

	print "subject1\tsubject2\tnumcnv1\tnumcnv2\toverlap1\toverlap2\n";
	
	for my $i (0 .. @pair-1) {
		my ($ind1, $ind2) = @{$pair[$i]};
		exists $cnv1{$ind1} or print STDERR "WARNING: No CNV call found for '$ind1' in $cnvfile\n";
		exists $cnv2{$ind2} or print STDERR "WARNING: No CNV call found for '$ind2' in $cnv2file\n";
		$cnv1{$ind1} ||= [];
		$cnv2{$ind2} ||= [];
		print $ind1, "\t", $ind2, "\t", scalar (@{$cnv1{$ind1}}), "\t", scalar (@{$cnv2{$ind2}}), "\t";
		print STDOUT scalar (cnvOverlap ($cnv1{$ind1}, $cnv2{$ind2}, $snp_chr, $snp_index)), "\t", scalar (cnvOverlap ($cnv2{$ind2}, $cnv1{$ind1}, $snp_chr, $snp_index)), "\n";
	}
}

sub cnvOverlap {
	my ($cnv1, $cnv2, $snp_chr, $snp_index) = @_;
	my @cnv1 = @$cnv1;
	my @cnv2 = @$cnv2;
	my %region2;
	my @return;
	
	for my $i (0 .. @cnv2-1) {
		$cnv2[$i] =~ m/^chr(\w+):(\d+)-(\d+).+startsnp=(\S+)\s+endsnp=(\S+)/ or confess "Error: invalid CNV call format found: <$cnv2[$i]>\n";
		push @{$region2{$1}}, [$2, $3, $4, $5];
	}
	for my $i (0 .. @cnv1-1) {
		$cnv1[$i] =~ m/^chr(\w+):(\d+)-(\d+).+startsnp=(\S+)\s+endsnp=(\S+)/ or confess "Error: invalid CNV call format found: <$cnv1[$i]>\n";
		my ($qchr, $qs, $qe, $qms, $qme) = ($1, $2, $3, $4, $5);
		my ($foundmatch, $over);
		if ($region2{$qchr}) {
			for my $j (0 .. @{$region2{$qchr}}-1) {
				my ($ts, $te, $tms, $tme) = @{$region2{$qchr}->[$j]};
				if ($qs <= $ts and $qe >= $ts) {
					if ($qe <= $te) {
						if ($minmarkeroverlap) {
							$over = $snp_index->{$qme}-$snp_index->{$tms}+1;
							if ($over/($snp_index->{$qme}-$snp_index->{$qms}+1)>=$minmarkeroverlap or $over/($snp_index->{$tme}-$snp_index->{$tms}+1)>=$minmarkeroverlap) {
								$foundmatch++;
							}
						} else {
							$over = $qe-$ts+1;
							if ($reciprocal) {
								if ($over/($qe-$qs+1)>=$minoverlap and $over/($te-$ts+1)>=$minoverlap) {
									$foundmatch++;
								}
							} else {
								if ($over/($qe-$qs+1)>=$minoverlap or $over/($te-$ts+1)>=$minoverlap) {
									$foundmatch++;
								}
							}
						}
					} else {
						if ($reciprocal) {
							if (($te-$ts+1) / ($qe-$qs+1)>=$minoverlap) {
								$foundmatch++;
							}
						} else {
							$foundmatch++;
						}
					}
				} elsif ($qs >= $ts and $qs <= $te) {
					if ($qe >= $te) {
						if ($minmarkeroverlap) {
							$over = $snp_index->{$tme}-$snp_index->{$qms}+1;
							if ($over/($snp_index->{$qme}-$snp_index->{$qms}+1)>=$minmarkeroverlap or $over/($snp_index->{$tme}-$snp_index->{$tms}+1)>=$minmarkeroverlap) {
								$foundmatch++;
							}
						} else {
							$over = $te-$qs+1;
							if ($reciprocal) {
								if ($over/($qe-$qs+1)>=$minoverlap and $over/($te-$ts+1)>=$minoverlap) {
									$foundmatch++;
								}
							} else {
								if ($over/($qe-$qs+1)>=$minoverlap or $over/($te-$ts+1)>=$minoverlap) {
									$foundmatch++;
								}
							}
						}
					} else {
						if ($reciprocal) {
							if (($qe-$qs+1) / ($te-$ts+1) >=$minoverlap) {
								$foundmatch++;
							}
						} else {
							$foundmatch++;
						}
					}
				}
			}
		}
		if ($foundmatch) {
			push @return, $cnv1[$i];
		}
	}
	$verbose and print STDERR "\n---------\nOVERLAP CNVs ARE\n", join ("\n", @return), "\n---------\n";
	return (@return);
}

sub readPFB {
	my ($pfbfile) = @_;
	
	my (%pfb, %snp_chr, %snp_index);
	my ($numrecord) = (0);
	open (PFB, $pfbfile) or confess "\nERROR: cannot read from pfb file $pfbfile: $!\n";
	print STDERR "NOTICE: Reading marker coordinates and population frequency of B allele (PFB) from $pfbfile ...";	
	while (<PFB>) {
		s/[\r\n]+$//;							#delete line feed and return characters
		m/^(\S+)\t(\S+)\t(\S+)\t(\S+)$/ or confess "\nERROR: invalid record found in PFB file $pfbfile (4 tab-delimited records expected): <$_>\n";
		my ($name, $chr, $pos, $pfb) = ($1, $2, $3, $4);
		$chr =~ m/^chr/i and next;					#this is the header line in a regular PFB file
		$numrecord++;
		push @{$pfb{$chr}}, [$pos, $name];
		$snp_chr{$name} = $chr;

	}
	close (PFB);
	print STDERR " Done with $numrecord records\n";

	for my $chr (keys %pfb) {
		@{$pfb{$chr}} = sort {$a->[0] <=> $b->[0]} @{$pfb{$chr}};
		for my $i (0 .. @{$pfb{$chr}}-1) {
			$snp_index{$pfb{$chr}->[$i][1]} = $i;
		}
	}
	return (\%snp_chr, \%snp_index);
}


sub generateRandFile {
	my ($num) = @_;
	$num < 100 or confess "Error: the maximum number of files are limited to 100";
	my (@file, %file);
	while (1) {
		my $nextfile = "rand" . substr (rand()*1000000, 0, 6);
		-f $nextfile and next;
		$file{$nextfile} and next;
		push @file, $nextfile;
		$file{$nextfile}++;
		@file==$num and last;
	}
	return @file;		
}


=head1 SYNOPSIS

 compare_cnv.pl [arguments] <operation> <cnvfile>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --listfile <file>		a tab-delimited 2-column file with two identifiers
 	    --minoverlap <float>	minimum fraction of overlapped bases in either call (default:0)
 	    --reciprocal		minoverlap refers to reciprocal overlap (both call)
 	    --minmarkeroverlap <float>	mimimum fraction of overlapped markers
 	    --pfbfile <file>		specify the PFB file with marker locations
 	    --logfile <file>		write log message to this file
 	    --outfile <file>		write output to this file

 Function: perform comparative analysis of CNV calls on duplicated samples in 
 the same file ('dup' operation), or on same sample called by different 
 algorithms in differnet files ('compcall' operation).

 Example: compare_cnv.pl compdup callfile -list idfile
          compare_cnv.pl compdup callfile -list idfile -minmarkeroverlap 0.8 -pfb hh550.pfb
          compare_cnv.pl compcall callfile1 -list idfile -cnv2 callfile2

 Version: '$LastChangedDate: 2011-05-03 11:31:25 -0700 (Tue, 03 May 2011) $';

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--listfile>

specify a file that contains a pair of identifiers in each line. This pair of 
identifier could be duplicated samples, or could represent the same individual but 
different file names were used in different calling algorithm.

=item B<--minoverlap>

the mimimum fraction of overlapped base pairs between two CNV calls to declare a 
match between two calls. By default this value is zero, meaning that two 
overlapped CNV calls are considered a match.

=item B<--reciprocal>

the minimum overlap refers to reciprocal overlap, that is, the fraction must be 
higher than the threshold for both CNV calls (as opposed to either CNV calls). 
For example, if two CNV calls have size of 300bp and 100bp, respectively and if 
the former encompass the latter, they would NOT be considered as a match if a -
minoverlap of 0.5 is specified with -reciprocal argument.

=item B<--minmarkeroverlap>

the mimimum fraction of overlapped markers between two CNV calls to declare a 
match between two calls. This argument requires the -pfbfile argument, which 
tells the marker locations so that the program can figure out what markers are 
contained within a CNV call.

=item B<--pfbfile>

the population frequency of B allele file, which is used here to supply marker 
chromosome coordinates for the -minmarkeroverlap argument.IT IS VERY IMPORTANT 
THAT THE PFBFILE CONTAINS ONLY THE SET OF MARKERS THAT ARE USED IN THE ACTUAL 
CNV CALLING, SO THAT THE FRACTION OF MARKER OVERLAP CAN BE ACCURATELY 
CALCULATED: for example, one cannot use the hhall.hg18.pfb file for a 
HumanCNV370 array, since there are extra markers annotated in the PFB file; 
instead, the user has to manually collect the subset of markers in the 370 array 
from the hhall.hg18.pfb file and then use the new PFB file for the 
compare_cnv.pl analysis.

=item B<--logfile>

write log information to this file

=item B<--outfile>

write output to this file

=back

=head1 DESCRIPTION

This program is used to perform comparative analysis of CNV calls on duplicated 
samples by the same algorithm, or on same sample called by different algorithms.

=over 8

=item * B<CNV call file>

The CNV call file must be in PennCNV call format. An example is given below:

	chr3:37957465-37961253        numsnp=3      length=3,789       state2,cn=1 father.txt startsnp=rs9837352 endsnp=rs9844203 conf=15.133
	chr3:75511365-75650909        numsnp=7      length=139,545     state2,cn=1 father.txt startsnp=rs4677005 endsnp=rs2004089 conf=26.862
	chr11:55127597-55204003       numsnp=11     length=76,407      state2,cn=1 mother.txt startsnp=rs2456022 endsnp=rs7934845 conf=23.106
	chr11:539119-548884           numsnp=4      length=9,766       state5,cn=3 mother.txt startsnp=rs4963136 endsnp=rs2061586 conf=7.342

Each line in the file corresponds to one CNV call. The various fields represent 
chromosome position, the number of markers (numsnp), the CNV length, the copy 
number (cn), the signal file name, the starting marker, the ending marker and 
the confidence score. The above file contains CNV calls for two files, including 
father.txt and mother.txt. The father.txt contains two deletions on chr3 (cn=1), 
while the mother.txt file contains one deletion and one duplication on chr11 
(cn=1 and cn=3).

=item * B<Comparing CNV calls on duplicated samples>

When comparing the CNV calls on duplicated samples by the same algorithm, an 
list file is necessary that lists one pair of sample identifiers per line 
separated by tab character. For example, a sample list file is given below:

	signalfile1	signalfile1_dup
	signalfile2	signalfile2_dup
	signalfile3	signalfile3_dup

In the file, each line lists one pair of samples, and the sample name is 
separated by a tab character.

When running the program, the command line should supply the operation (dup), 
the input CNV call file (ex1.rawcnv), and use the -list argument to specify the 
list file that contains pairs of samples. For example, the following command may 
be issued:

	[kaiwang@cc ~/project/penncnv/example]$ compare_cnv.pl dup ex1.rawcnv -list complist 
	NOTICE: <operation>='dup/compdup': find concordant CNV calls between duplicated samples (same sample but different file names)
	subject1        subject2        numcnv1 numcnv2 overlap1        overlap2
	signalfile1     signalfile1_dup   4       4       2       2
	signalfile2	signalfile2_dup   2       4       1       1

The output is tab-delimited to facilitate post-processing in Excel or other 
software. For each line, the identifier of two subjects (subject1 and subject2) 
are first printed, then the number of CNV calls for each subject is printed, and 
then the overlapping CNV calls for each subject is printed. For example, the 
signalfile2 and signalfile2_dup contains 2 and 4 CNV calls, respectively. Among 
the calls in signalfile2, 1 of them is detected in signalfile2_dup; among the 4 
CNV calls in signalfile2_dup, 1 of htem is detected in signalfile2.

=item * B<Comparing CNV calls on the same sample by different algorithms>

Sometimes it would be interesting to compare the CNV calls on the same sample by 
different algorithms, for example, calls by Birdeye and PennCNV on Affymetrix 
arrays. Suppose the two CNV call files by two algorithms are cnvfile1 and 
cnvfile2, and both of them should be in PennCNV formats (the convert_cnv.pl 
program can be used for converting CNV calls from one format to another). A list 
file should also be given that contains a pair of identifiers per line, 
representing the signal file names in the first call file and second call file, 
respectively. Note that different algorithms may have used different data 
processing models, so the signal file names given by different CNV calling 
algorithms may differ. For example, an example list file is given below:

	79585_PERDU_g_Plate_No._30_GenomeWideSNP_6_C09_51218.CEL	gw6.GenomeWideSNP_6_C09_51218.CEL
	79585_PERDU_g_Plate_No._30_GenomeWideSNP_6_C09_51218.CEL	gw6.GenomeWideSNP_6_C09_51218.CEL

The first file in each line is the file name used in one CNV calling algorithm, 
and was shown as the fifth field in the cnvfile1; the second file name 
corresponds to the name specified in cnvfile2.

To compare the CNV calls generated on the same subject, the following command 
may be issued:

	[kaiwang@cc ~/project/penncnv/example]$ compare_cnv.pl compcall ex1.rawcnv -cnv2 ex2.triocnv -list complist 
	NOTICE: <operation>='compcall': compare CNV calls from two files on specified pairs of samples
	subject1        subject2        numcnv1 numcnv2 overlap1        overlap2
	bird.file1      gw6.file1   4       6      3       6
	bird.file2      gw6.file2   2       4      2       4

The output is tab-delimited to facilitate post-processing in Excel or other 
software.  For each line, the identifier of two subjects (subject1 and subject2) 
are first printed, then the number of CNV calls for each algorithm is printed, 
and then the overlapping CNV calls is printed. For example, the bird.file1 and 
gw6.file1 contains 2 and 4 CNV calls by two CNV calling algorithms, 
respectively. Among the calls given by the first algorithm, 3 of them are 
detected by the second algorithm; among the 6 CNV calls given by the second 
algorithm, 7 of htem are detected by the first algorithm.

For questions, comments or bug reports, please contact me at 
kai@openbioinformatics.org.

=back

