#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;

our ($help, $verbose, $man);
our ($path_detect_cnv, $path_visualize_cnv, $path_convert_cnv, $path_filter_cnv, $path_compare_cnv, $path_infer_allele, $command);

GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'path_detect_cnv=s'=>\$path_detect_cnv, 'path_visualize_cnv=s'=>\$path_visualize_cnv,
	'path_convert_cnv=s'=>\$path_convert_cnv, 'path_filter_cnv=s'=>\$path_filter_cnv, 'path_compare_cnv=s'=>\$path_compare_cnv, 'path_infer_allele=s'=>\$path_infer_allele) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error");


$path_detect_cnv ||= "detect_cnv.pl";
$path_visualize_cnv ||= "visualize_cnv.pl";
$path_convert_cnv ||= "convert_cnv.pl";
$path_filter_cnv ||= "filter_cnv.pl";
$path_compare_cnv ||= "compare_cnv.pl";
$path_infer_allele ||= "infer_snp_allele.pl";

if ($ARGV[0] == 1) {
	$command = "$path_detect_cnv -test -hmm example.hmm -pfb example.pfb -conf -log ex1.log -out ex1.rawcnv father.txt mother.txt offspring.txt";
	print STDERR "Exercise 1: individual-based calling and write the output to ex1.rawcnv\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 2) {
	$command = "$path_detect_cnv -trio -hmm example.hmm -pfb example.pfb -log ex2.log -out ex2.triocnv -cnv ex1.rawcnv father.txt mother.txt offspring.txt";
	print STDERR "Exercise 2: tiro-based calling and write the output to ex2.triocnv\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 3) {	
	$command = "$path_detect_cnv -test -hmm example.hmm -pfb example.pfb -log ex3.log -out ex3.rawcnv -gcmodel ../lib/hh550.hg18.gcmodel -conf -list inputlist";
	print STDERR "Exercise 3: individual-based calling with GCmodel signal adjustment, write to ex3.rawcnv\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 4) {
	$command = "$path_detect_cnv -validate -hmm example.hmm -pfb example.pfb -log ex4.log -out ex4.rawcnv -startsnp rs8114269 -endsnp rs682562 -delfreq 0.005 -list inputlist";
	print STDERR "Exercise 4: validation-based calling for the region between rs8114969 and rs682562\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 5) {
	$command = "$path_detect_cnv -validate -hmm example.hmm -pfb example.pfb -log ex5.log -out ex5.rawcnv -candlist ccnvr -list inputlist";
	print STDERR "Exercise 5: validation-based calling for all CNV regions annotated in 'ccnvr' file\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 6) {
	$command = "$path_detect_cnv -joint -hmm example.hmm -pfb example.pfb -log ex6.log -out ex6.jointcnv father.txt mother.txt offspring.txt";
	print STDERR "Exercise 6: joint CNV calling for a trio and write to ex6.jointcnv\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 7) {
	$command = "$path_visualize_cnv -intype cnv -format bed -track 'example track' -out ex1.bed ex1.rawcnv";
	print STDERR "Exercise 7: convert CNV call in ex1.rawcnv to BED format and write to ex1.bed\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 8) {
	$command = "$path_convert_cnv -intype penncnv -outtype tab -output ex1.tabcnv ex1.rawcnv";
	print STDERR "Exercise 8: convert CNV call in ex3.rawcnv to tab-delimited format and write to ex1.tabcnv\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 9) {
	$command = "$path_convert_cnv -intype tab -outtype penncnv ex1.tabcnv";
	print STDERR "Exercise 9: convert CNV call in ex1.newcnv from tab-delimited format to PennCNV format\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 10) {
	$command = "$path_filter_cnv -numsnp 10 -length 50k -type del ex1.rawcnv";
	print STDERR "Excisise 10: filter CNV calls in ex1.rawcnv and print out only deletions with >=10 SNPs and >=50kb\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 11) {
	$command = "$path_compare_cnv compdup ex1.rawcnv -list list.compdup";
	print STDERR "Exercise 11: compare CNV calls between pairs of samples\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 12) {
	$command = "$path_compare_cnv compcall ex1.rawcnv -list list.compcall -cnv2 ex2.triocnv";
	print STDERR "Exercise 12: compare CNV calls between pairs of samples\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 13) {
	$command = "$path_infer_allele -pfb example.pfb -hmm example.hmm -allcn 221 -start rs11716390 -end rs17039742 -out ex13.geno father.txt mother.txt offspring.txt";
	print STDERR "Exercise 13: generate CNV-based genotype calls)\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 14) {
	$command = "$path_infer_allele -pfb example.pfb -hmm example.hmm -denovocn 1 -start rs11716390 -end rs17039742 -out ex14.geno father.txt mother.txt offspring.txt";
	print STDERR "Exercise 14: validate de novo CNVs and assign P-values\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 15) {
	$command = "convert_cnv.pl -intype canary -outtype penncnv -canarydef GenomeWideSNP_6.hg18.cnv_defs -output ex15.rawcnv example.canary_calls";
	print STDERR "Exercise 15: convert Canary CNV calls to PennCNV format\n\tRunning command <$command>\n";
} elsif ($ARGV[0] == 16) {
	$command = "visualize_cnv.pl -format plot -signal offspring.txt ex1.rawcnv";
	print STDERR "Exercise 16: generate plots of LRR and BAF values for CNV calls in JPG format (R must be installed in the system for plotting)\n\tRunning command <$command>\n";
}




print STDERR "\n******************************************************************************\n";
system ($command);
print STDERR "\n******************************************************************************\n";






=head1 SYNOPSIS

 runex.pl <1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --path_detect_cnv <string>		path to detect_cnv.pl
 	    --path_visualize_cnv <string>	path to visualize_cnv.pl
 	    --path_convert_cnv <string>		path to convert_cnv.pl
 	    --path_filter_cnv <string>		path to filter_cnv.pl
 	    --path_compare_cnv <string>		path to compare_cnv.pl

 Function: test-drive PennCNV and related scripts in your system

 Example: runex.pl 1 (run PennCNV to call CNV on three signal files)
 	  runex.pl 2 (run posterior CNV calling algorithm on a trio)
 	  runex.pl 3 (run PennCNV with GCmodel adjustment of signal intensity)
 	  runex.pl 4 (validation-based CNV calling on a candidate region)
 	  runex.pl 5 (validation-based CNV calling on all candidate regions in a file)
 	  runex.pl 6 (run joint CNV calling on a trio)
 	  runex.pl 7 (convert CNV call to BED format)
 	  runex.pl 8 (convert CNV call to tab-delimited format)
 	  runex.pl 9 (convert tab-delimited calls to PennCNV format)
 	  runex.pl 10 (filter CNV and print out a subset of calls)
 	  runex.pl 11 (compare CNV calls on duplicated samples)
 	  runex.pl 12 (compare CNV calls on same sample called by different algorithms)
 	  runex.pl 13 (generate CNV-based genotype calls)
 	  runex.pl 14 (validate de novo CNV calls and assign a P-value)
 	  runex.pl 15 (convert Canary CNV calls to PennCNV format)
 	  runex.pl 16 (plot signal intensity for CNV calls in JPG formats)

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=back

=head1 DESCRIPTION

This program is used to test-drive PennCNV and related scripts using several 
simplified example commands.

=cut
