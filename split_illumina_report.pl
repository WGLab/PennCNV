#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2010-11-28 12:47:21 -0800 (Sun, 28 Nov 2010) $';

our ($verbose, $help, $man);
our ($reportfile, $prefix, $suffix, $numeric_name, $comma, $tolerate, $revised_file);

GetOptions ('verbose'=>\$verbose, 'help'=>\$help, 'man'=>\$man, 'prefix=s'=>\$prefix, 'suffix=s'=>\$suffix, 'numeric_name'=>\$numeric_name, 'comma'=>\$comma, 'tolerate'=>\$tolerate, 
	'revised_file=s'=>\$revised_file) or pod2usage ();

# If used, %revised is populated with content from $revised_file.
our %revised = (); # Key = sample ID in Illumina report; value = sample ID to use.

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error");

($reportfile) = @ARGV;
$prefix ||= '';
$suffix ||= '';


#
# Populate %revised if necessary:

if ( $revised_file ) {
	open (REV, $revised_file) or confess "Error: cannot read file $revised_file :$!\n";
	while ( my $line = <REV> ) {
		chomp($line);
		my ($orig_ID, $revised_ID) = split /\s+/, $line;
		if ( exists $revised{$orig_ID} ) {
			warn "\nOriginal ID $orig_ID already encountered in revised_file; keeping old value.\n";
		}
		else {
			$revised{$orig_ID} = $revised_ID;
		}
	}
}



splitIlluminaReport ($reportfile, $prefix, $suffix, $numeric_name, $comma);



sub splitIlluminaReport {
	my ($reportfile, $prefix, $suffix, $numeric_name, $comma) = @_;
	my ($count_line, $count_file, $pre, $name_index, $sample_index, $lrr_index, $baf_index) = (0, 0, 'NA');
	my (@field);
	
	open (REPORT, $reportfile) or confess "Error: cannot read from input report file $reportfile: $!\n";
	while (<REPORT>) {
		$count_line++;
		m/^\[Data\]/ and last;
		$count_line > 1000 and confess "Error: after reading 1000 lines in $reportfile, still cannot find [Data] section. The $reportfile file may not be in Illumina report format.\n";
	}
	
	$_ = <REPORT>;
	s/[\r\n]+$//;
	$count_line++;
	@field = $comma ? (split (/,/, $_)) : (split (/\t/, $_));
	@field >= 4 or confess "Error: invalid header line (at least 4 tab-delimited fields, including 'SNP Name', 'Sample ID', 'B Allele Freq', 'Log R Ratio' expected) in report file $reportfile: <$_>\n";
	
	for my $i (0 .. @field-1) {
		$field[$i] eq 'SNP Name' and $name_index = $i;
		$field[$i] eq 'Sample ID' and $sample_index = $i;
		$field[$i] eq 'B Allele Freq' and $baf_index = $i;
		$field[$i] eq 'Log R Ratio' and $lrr_index = $i;
	}
	defined $name_index or confess "Error: the 'SNP Name' field is not found in header line in report file $reportfile: <$_>\n";
	defined $sample_index or confess "Error: the 'Sample ID' field is not found in header line in report file $reportfile: <$_>\n";
	defined $baf_index or confess "Error: the 'B Allele Freq' field is not found in header line in report file $reportfile: <$_>\n";
	defined $lrr_index or confess "Error: the 'Log R Ratio' field is not found in header line in report file $reportfile: <$_>\n";
	
	my $got_sample_ID = 0; # Use this so we don't report on multiple adjacent lines

	while (<REPORT>) {
		s/[\r\n]+$//;
		$count_line++;
		@field = $comma ? (split (/,/, $_)) : (split (/\t/, $_));
		@field >= 4 or confess "Error: invalid data line (at least 4 tab- or comma-delimited fields expected) in report file $reportfile: <$_>\n";
		defined $field[$name_index] or confess "Error: the 'SNP Name' field is not found in data line in report file $reportfile: <$_>\n";

		#defined $field[$sample_index] or confess "Error: the 'Sample ID' field is not found in data line in report file $reportfile line $.; line is:\n<$_>\n";
	
		# Sometimes a "blank" sampleID can be produced in BeadStuidio files
		# (i,.e. just whitespace). Be resilient to this.

		if ( $field[$sample_index] !~ /\S+/ ) {
		  if ( $got_sample_ID == 0 ) {
			# already reported this; keep quiet.
		  }
		  else { # Just found this bad'un.
			$got_sample_ID = 0;
			print "\n***************\nWARNING: No sample ID found at line $.; line is:\n<$_>\nSkipping...\n***************\n\n";
		  }
		  next;
		}
		else {
		  $got_sample_ID = 1; # OK!
		}

		if (not defined $field[$baf_index]) {
			if ($tolerate) {
				print STDERR "WARNING: Skipping marker $field[$name_index] for $field[$sample_index] due to lack of BAF information\n";
				next;
			} else {
				confess "Error: the 'B Allele Freq' field is not found in data line in report file $reportfile: <$_>\n";
			}
		}
		if (not defined $field[$lrr_index]) {
			if ($tolerate) {
				print STDERR "WARNING: Skipping marker $field[$name_index] for $field[$sample_index] due to lack of LRR information\n";
				next;
			} else {
				confess "Error: the 'Log R Ratio' field is not found in data line in report file $reportfile: <$_>\n";
			}
		}
		
		if ($field[$sample_index] eq $pre) {
			print OUT join ("\t", @field[$name_index, $lrr_index, $baf_index]), "\n";
   		} else {
			$count_file++;
			my $outname;
			my $revised_sample;
			if ($numeric_name) {
				$outname = "$prefix"."split$count_file$suffix";
			} else {
				$outname = "$prefix$field[$sample_index]$suffix";
				if ($revised_file) {
					my $key = $field[$sample_index];
					$revised_sample = $revised{$key};
					$outname = "$prefix$revised_sample$suffix";
					# print "\noutname is $outname.";
				}
			}
			print STDERR "NOTICE: Writing to output signal file $outname\n";
			open (OUT, ">$outname") or confess "Error: cannot write to output file $outname:$!\n";
			if ( $revised_file ) {
				print OUT "Name\t$revised_sample.Log R Ratio\t$revised_sample.B Allele Freq\n";
			}
			else {
				print OUT "Name\t$field[$sample_index].Log R Ratio\t$field[$sample_index].B Allele Freq\n";
			}
			print OUT join("\t",@field[$name_index, $lrr_index, $baf_index]), "\n";
		}
		$pre = $field[$sample_index];
	}
	print STDERR "NOTICE: Finished processing $count_line lines in report file $reportfile, and generated $count_file output signal intensity files\n";
}


=head1 SYNOPSIS

 split_illumina_report.pl [arguments] <reportfile>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	-p, --prefix <string>		prefix of output file name
 	-s, --suffix <string>		suffix of output file name
 	-n, --numeric_name		use numeric file name (default: Sample ID is file name)
 	-c, --comma			fields are comma-delimited (default: fields are tab-delimited)
 	-t, --tolerate			tolerate records without LRR/BAF information
 	-r, --revised_file <file>	path to "revised" file of alternate sample IDs
 	    --tolerate			tolerate when LRR/BAF do not exist in input file


 Function: split the Illumina report file to individual signal intensity files

 Example: split_illumina_report.pl -prefix signal/ -suffix .txt HapMap.report
          split_illumina_report.pl -num HapMap.report

 Version: $LastChangedDate: 2010-11-28 12:47:21 -0800 (Sun, 28 Nov 2010) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--prefix>

specify the prefix for output file name. By default the file name is the "Sample 
ID" field in the report file.

=item B<--suffix>

specify the suffix for output file name. By default the file name is the "Sample 
ID" field in the report file.

=item B<--numeric_name>

specify that the output file name should be arranged in a numeric manner (such 
as split1, split2, split3, etc.).

=item B<--comma>

specify that the data fields in report file is comma-delimited. By default, the 
tab-delimited format is assumed.

=item B<--tolerate>

tolerate the problem when Log R Ratio (LRR) or B Allele Frequency (BAF) values 
do not exist in input file. In this case, the program will still run, but the 
output will not be immediately useful for CNV inference. For example, sometimes 
one may get a Report file with X and Y values, and can generate splitted 
individual files with X/Y values, and then use a simple script to convert the 
X/Y to LRR/BAF for each individual.

=item B<--revised_file>

specify path to file listing sampleIDs to use instead of those appearing in
the Illumina report file. The revised file must have 2 whitespace-separated
columns, without a header; the first column is the original IDs (i.e. those
in the Illumina report file), and the second the corresponding (revised) ID
to use instead. You must provide one such mapping for each sampleID in the
Illumina file.


=back

=head1 DESCRIPTION

This program is used to convert the "Report file" from Illumina BeadStudio 
software to a format that can be used by the detect_cnv.pl program for CNV 
detection.

The report file should contain at least four tab-delimited fields, including 
'SNP Name', 'Sample ID', 'B Allele Freq', 'Log R Ratio'. If not, try to generate 
the report file again in BeadStudio, and make sure that these fields are 
selected to be exported.

It is okay for the report files to contain additional columns, but they will be 
ignored and not used in the output files. For example, the first a few lines of 
a typical Illumina report file is presented below:

	[Header]
	GSGT Version    1.1.4
	Processing Date 3/17/2009 11:12 AM
	Content         HumanCNV100W_B.bpm
	Num SNPs        95484
	Total SNPs      95484
	Num Samples     89
	Total Samples   285
	[Data]
	SNP Name        Sample ID       Allele1 - Forward       Allele2 - Forward       GC Score        X       Y       X Raw   Y Raw   B Allele Freq   Log R Ratio
	cnvi0048890     NA06985 -       -       0.0000  0.211   1.793   3410    23376   0.9847  0.0003
	cnvi0048891     NA06985 -       -       0.0000  0.832   0.031   10373   966     0.0000  0.3365
	cnvi0048892     NA06985 -       -       0.0000  1.297   0.058   15961   1400    0.0027  0.1453
	cnvi0048894     NA06985 -       -       0.0000  0.451   0.193   5841    2974    0.0000  0.1799
	cnvi0048895     NA06985 -       -       0.0000  0.677   0.028   8514    904     0.0000  0.4093
	cnvi0048896     NA06985 -       -       0.0000  2.708   0.638   33072   9072    0.0816  0.1205
	cnvi0048897     NA06985 -       -       0.0000  2.757   0.718   33681   10101   0.0579  0.0254
	cnvi0048898     NA06985 -       -       0.0000  2.409   0.506   29441   7323    0.0119  0.0685
	cnvi0048899     NA06985 -       -       0.0000  2.217   0.590   27166   8369    0.0494  -0.0141
	cnvi0048900     NA06985 -       -       0.0000  2.515   0.177   30625   3139    0.0198  -0.0611
	cnvi0048901     NA06985 -       -       0.0000  0.146   2.119   2715    27531   0.9895  0.0646
	cnvi0048902     NA06985 -       -       0.0000  1.474   0.130   18110   2353    0.0106  -0.2371
	cnvi0048903     NA06985 -       -       0.0000  0.107   1.963   2204    25532   0.9970  -0.1069

Each line in the file include measurements for one marker for one sample. The 
'SNP Name', 'Sample ID', 'B Allele Freq', 'Log R Ratio' fields will be examined 
and be written to output files.


=over 8

=item * B<Explanation of the -revised_file argument>

Sometimes users may want to use a different sample ID in the splitted output 
file, and the --revised_file argument can be useful. The revised file must have 
2 whitespace-separated columns, without a header; the first column is the 
original IDs (i.e. those in the Illumina report file), and the second the 
corresponding (revised) ID to use instead. For example, given an Illumina file 
with lines such as:

	...
	200003  4605306064_R01C01       A       A       0.9175  0.028   1.042   0.998   0.044   14503   1933    0.0000  -0.0296
	200006  4605306064_R01C01       A       G       0.7567  0.441   2.002   1.093   0.908   16015   16156   0.4824  -0.0901
	...

by default the original program will create a file called "4605306064_R01C01", with header line:

	Name    4605306064_R01C01.Log R Ratio   4605306064_R01C01.B Allele Freq

However, alternative ID can be used in the output. Suppose this is 
internal_ID_123. Then, given a file with entries such as:

	...
	4605306064_R01C01   internal_ID_123
	...

One can run the script with the -r option and get a file called "internal_ID_123", with header:

	Name    internal_ID_123.Log R Ratio   internal_ID_123.B Allele Freq

=item * B<File format issues>

It is possible for Illumina final call report files to be output such that they 
are grouped by SNP and not sample. This means that the results from different 
samples would be interleaved with the results from other samples. Currently the 
script cannot cope with this; all it will do is to output the data from the last 
occurrence of each sample.

=back

This program was modified with functional enhancement by Matthew Gillman @ 
Wellcome Trust Sanger.

For questions, comments or bug reports, please contact me at 
kai@openbioinformatics.org.

=cut
