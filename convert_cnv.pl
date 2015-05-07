#!/usr/bin/env perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2011-06-16 22:13:39 -0700 (Thu, 16 Jun 2011) $';

our ($verbose, $help, $man);
our ($intype, $outtype, $output, $snplocfile, $canarydeffile, $oncosnpqcfile);
our ($cnvfile);

GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'intype=s'=>\$intype, 'outtype=s'=>\$outtype, 'output=s'=>\$output, 
	'snplocfile=s'=>\$snplocfile, 'canarydeffile=s'=>\$canarydeffile, 'oncosnpqcfile=s'=>\$oncosnpqcfile) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error: <inputfile> missing");

($cnvfile) = @ARGV;

$intype or print STDERR "NOTICE: the --intype argument is set as 'penncnv' by default\n" and $intype = 'penncnv';
$outtype or print STDERR "NOTICE: the --outtype argument is set as 'tab' by default\n" and $outtype = 'tab';
$intype eq $outtype and pod2usage ("Error in argument: the --intype ($intype) cannot be identical to --outtype ($outtype)");

$intype =~ m/^(penncnv|tab|xml|birdseye|birdsmerge|canary|oncosnp)$/ or pod2usage ("Error in argument: the --intype argument can be only 'penncnv' or 'tab' or 'birdseye' or 'birdsmerge' or 'canary' or 'oncosnp' at this time");
$outtype =~ m/(penncnv|tab)$/ or pod2usage ("Error in argument: the --outtype argument can be only 'penncnv' or 'tab' at this time");

if ($output) {
	open (STDOUT, ">$output") or confess "Error: cannot write to output file $output: $!\n";
}

our $infh;
if ($cnvfile eq 'stdin') {
	$infh = *STDIN;
} else {
	open ($infh, $cnvfile) or confess "Error: cannot read from input CNV file $cnvfile: $!\n";
}

if ($snplocfile) {
	$intype eq 'xml' or $intype eq 'birdseye' or $intype eq 'birdsmerge' or $intype eq 'canary' or pod2usage ("Error in argument: the --snplocfile is useful only for --intype of xml or birdseye or birdsmerge or canary");
}

if ($intype eq 'penncnv') {
	if ($outtype eq 'tab') {
		convertPennCNVToTab ($infh);
	} else {
		pod2usage ("Error: the --outtype $outtype is not supported for the --intype $intype");
	}
} elsif ($intype eq 'tab') {
	if ($outtype eq 'penncnv') {
		convertTabToPennCNV ($infh);
	} else {
		pod2usage ("Error: the --outtype $outtype is not supported for the --intype $intype");
	}
} elsif ($intype eq 'birdseye') {
	if ($outtype eq 'penncnv') {
		convertBirdseyeToPennCNV ($infh, $snplocfile);
	} else {
		pod2usage ("Error: the --outtype $outtype is not supported for the --intype $intype");
	}
} elsif ($intype eq 'birdsmerge') {
	if ($outtype eq 'penncnv') {
		convertBirdsmergeToPennCNV ($infh, $snplocfile);
	} else {
		pod2usage ("Error: the --outtype $outtype is not supported for the --intype $intype");
	}
} elsif ($intype eq 'canary') {
	$canarydeffile or pos2usage ("Error: please specify the --canarydeffile argument");
	if ($outtype eq 'penncnv') {
		convertCanaryToPennCNV ($infh, $canarydeffile);
	} else {
		pod2usage ("Error: the --outtype $outtype is not supported for the --intype $intype");
	}
} elsif ($intype eq 'xml') {
	if ($outtype eq 'penncnv') {
		convertXMLToPennCNV ($infh, $snplocfile);
	} else {
		pod2usage ("Error: the --outtype $outtype is not supported for the --intype $intype");
	}
} elsif ($intype eq 'oncosnp') {
	my $qcresult;
	$oncosnpqcfile and $qcresult = readOncosnpQC ($oncosnpqcfile);
	if ($outtype eq 'penncnv') {
		convertOncosnpToPennCNV ($infh, $qcresult);
	} elsif ($outtype eq 'tab') {
		convertOncosnpToTab ($infh, $qcresult);
	} else {
		pod2usage ("Error: the --outtype $outtype is not supported for the --intype $intype");
	}
}

sub readOncosnpQC {
	my ($oncosnpqcfile) = @_;
	my (%qcresult);
	my (@alpha, @lrr_shift, @ll, $index);
	open (FH, $oncosnpqcfile) or die "Error: cannot read from oncosnpqcfile $oncosnpqcfile: $!\n";
	while (<FH>) {
		if (m/Outlier Rate:\s+(\S+)$/) {
			$qcresult{outlier} = $1;
		} elsif (m/^Std\. Dev\. LRR:\s+(\S+)$/) {
			$qcresult{lrr_sd} = $1;
		} elsif (m/^Std\. Dev\. BAF:\s+(\S+)$/) {
			$qcresult{baf_sd} = $1;
		} elsif (m/^Stromal Contamination:\s+(\S+)\s+(\S+)\s*(\S+?)/) {
			@alpha = ($1, $2, $3);
		} elsif (m/^Log R Ratio Baseline Shift:\s+(\S+)\s+(\S+)\s*(\S+?)/) {
			@lrr_shift = ($1, $2, $3);
		} elsif (m/^Log-likelihood:\s+(\S+)\s+(\S+)\s*(\S+?)/) {
			@ll = ($1, $2, $3);
		}
	}
	if ($ll[0]<$ll[1]) {
		$index = 1;
	} else {
		$index = 0;
	}
	if (@ll==3 and $ll[$index] < $ll[2]) {
		$index = 2;
	}
	$qcresult{stromal} = $alpha[$index];
	$qcresult{lrr_shift} = $lrr_shift[$index];
	$qcresult{ll} = $ll[$index];
	defined $qcresult{stromal} or die "Error: stromal contamination not found in oncosnpqc file\n";
	$verbose and print STDERR "NOTICE: LRR shift values are $lrr_shift[$index] for $oncosnpqcfile\n";
	return (\%qcresult);
}
		
		
sub printPennCNVFormat {
	my ($chr, $start, $end, $cn, $sample, $startsnp, $endsnp, $conf, $numsnp, $state) = @_;
	my $cnvregion = "chr$chr:$start-$end";
	if (length ($cnvregion) < 29) {
		$cnvregion = substr ("$cnvregion                              ", 0, 29);
	}

	if (length($numsnp) < 6) {
		$numsnp = substr ("$numsnp      ", 0, 6);
	}
	
	my $cnvlength = join ('', reverse split (//, $end-$start+1)); $cnvlength =~ s/(...)/$1,/g; $cnvlength =~ s/,$//; $cnvlength = join ('', reverse split (//, $cnvlength));
	if (length($cnvlength) < 11) {
		$cnvlength = substr ("$cnvlength            ", 0, 11);
	}
	
	if (not defined $state) {
		$state=$cn+1; $cn>=2 and $state++;
	}
	
	print "$cnvregion numsnp=$numsnp length=$cnvlength state$state,cn=$cn $sample startsnp=$startsnp endsnp=$endsnp";
	$conf ne 'NA' and print " conf=$conf";
	print "\n";
}
	
sub convertBirdseyeToPennCNV {
	my ($fh, $snplocfile) = @_;
	
	$verbose and print STDERR "NOTICE: Converting Birdseye CNV calls to PennCNV format\n";
	
	my ($pos_name);
	if ($snplocfile) {
		$pos_name = readSNPLocFile ($snplocfile);
	}

	
	$_ = <$fh>;
	s/[\r\n]+$//;
	$_ =~ m/^sample\tsample_index\tcopy_number\tchr\tstart\tend\tper_probe_score\tsize\tnum_probes\tlod_score/ or confess "Error: invalid first line in the Birdeye call file: <$_>\n";
	
	my $count_nonautosome = 0;
	while (<$fh>) {
		s/[\r\n]+$//;
		my @field = split (/\t/, $_);
		@field == 10 or confess "Error: invalid record found in Birdeye call file (10 tab-delimited fields expected): <$_>\n";
		my ($sample, $index, $cn, $chr, $start, $end, $perprobescore, $length, $numsnp, $conf) = @field;
		my ($startsnp, $endsnp) = ("NA", "NA");
		
		unless ($chr =~ m/^\d+$/ and $chr>=1 and $chr <= 22) {
			$count_nonautosome++;
			next;
		}
		
		if ($snplocfile) {
			$startsnp = $pos_name->{$chr, $start} || 'NA';
			$endsnp = $pos_name->{$chr, $end} || 'NA';
		}
		
		if ($chr =~ m/^\d+$/) {
			if ($chr >= 1 and $chr <= 22) {
				$cn eq '2' and next;
			}
		}
		
		printPennCNVFormat ($chr, $start, $end, $cn, $sample, $startsnp, $endsnp, $conf, $numsnp);
	}
	$count_nonautosome and print STDERR "WARNING: Total of $count_nonautosome non-autosome CNV calls are skipped for conversion\n";
}

sub convertOncosnpToPennCNV {
	my ($fh, $qcresult) = @_;
	
	$verbose and print STDERR "NOTICE: Converting Birdseye CNV calls to PennCNV format\n";
	
	$_ = <$fh>;
	s/[\r\n]+$//;
	$_ =~ m#^Chromosome\tStart Position \(bp\)\tEnd Position \(bp\)\tLength / Mb# or confess "Error: invalid first line in the Oncosnp call file: <$_>\n";
	my @field = split (/\t/, $_);
	@field >= 32 or confess "Error: invalid record found in Oncosnp call file (at least 8 tab-delimited fields expected): <$_>\n";
	my $version;
	if ($field[10] eq 'Max. Log BF') {
		$version = 1.0;
	} elsif ($field[10] eq '% Normal') {
		$version = 1.1;
	} else {
		die "Error: unable to identify Oncosnp version from inputfile\n";
	}
	
	my $count_nonautosome = 0;
	while (<$fh>) {
		s/[\r\n]+$//;
		my @field = split (/\t/, $_);
		@field >= 32 or die "Error: invalid record found in Oncosnp call file (at least 8 tab-delimited fields expected): <$_>\n";
		my ($chr, $start, $end, $length, $startsnp, $endsnp, $numsnp, $cn, $loh, $state, $percent_norm, $conf) = @field;
		$version == 1.0 and $conf = $percent_norm;
		$start += 0;		#convert string to number
		$end += 0;		#convert string to number
		$length = $end-$start+1;		#the length in the inputfile is calculated in Mb and is not precise
		
		my $sample = $cnvfile;
		
		$conf >= 10 or next;
		printPennCNVFormat ($chr, $start, $end, $cn, $sample, $startsnp, $endsnp, $conf, $numsnp, $state);
	}
	$count_nonautosome and print STDERR "WARNING: Total of $count_nonautosome non-autosome CNV calls are skipped for conversion\n";
}

sub convertOncosnpToTab {
	my ($fh, $qcresult) = @_;
	
	$verbose and print STDERR "NOTICE: Converting Birdseye CNV calls to PennCNV format\n";
	
	$_ = <$fh>;
	s/[\r\n]+$//;
	$_ =~ m#^Chromosome\tStart Position \(bp\)\tEnd Position \(bp\)\tLength / Mb# or die "Error: invalid first line in the Oncosnp call file: <$_>\n";
	my @field = split (/\t/, $_);
	@field >= 32 or confess "Error: invalid record found in Oncosnp call file (at least 8 tab-delimited fields expected): <$_>\n";
	my $version;
	if ($field[10] eq 'Max. Log BF') {
		$version = 1.0;
		$qcresult->{stromal} or die "Error: oncosnp 1.0 output file detected. Please specify --oncosnpqc argument\n";
	} elsif ($field[10] eq '% Normal') {
		$version = 1.1;
	} else {
		die "Error: unable to identify Oncosnp version from inputfile\n";
	}
	
	print STDOUT "chr	start	end	state	cn	sample	snp1	snp2	alpha\n";
	my $count_nonautosome = 0;
	while (<$fh>) {
		s/[\r\n]+$//;
		my @field = split (/\t/, $_);
		@field >= 32 or confess "Error: invalid record found in Oncosnp call file (at least 8 tab-delimited fields expected): <$_>\n";
		my ($chr, $start, $end, $length, $startsnp, $endsnp, $numsnp, $cn, $loh, $state, $percent_norm, $conf) = @field;
		$version == 1.0 and $conf = $percent_norm;
		$start += 0;		#convert string to number
		$end += 0;		#convert string to number
		$length = $end-$start+1;		#the length in the inputfile is calculated in Mb and is not precise
		
		my $sample = $cnvfile;
		
		$conf >= 10 or next;
		if ($version == 1.0) {
			print STDOUT join ("\t", $chr, $start, $end, $state, $cn, $sample, $startsnp, $endsnp, $qcresult->{stromal}), "\n";
		} elsif ($version == 1.1) {
			print STDOUT join ("\t", $chr, $start, $end, $state, $cn, $sample, $startsnp, $endsnp, $percent_norm), "\n";
		}
	}
	$count_nonautosome and print STDERR "WARNING: Total of $count_nonautosome non-autosome CNV calls are skipped for conversion\n";
}

sub convertBirdsmergeToPennCNV {
	my ($fh, $snplocfile) = @_;
	
	$verbose and print STDERR "NOTICE: Converting Birdseye CNV calls to PennCNV format\n";
	
	my ($pos_name);
	if ($snplocfile) {
		$pos_name = readSNPLocFile ($snplocfile);
	}

	
	$_ = <$fh>;
	$_ =~ m/^sample\tsample_index\tcopy_number\tchr\tstart\tend\tconfidence/ or confess "Error: invalid first line in the Birdmerge call file: <$_>\n";
	
	my $count_nonautosome = 0;
	while (<$fh>) {
		s/[\r\n]+$//;
		my @field = split (/\t/, $_);
		@field == 7 or confess "Error: invalid record found in Birdsmerge call file (7 tab-delimited fields expected): <$_>\n";
		my ($sample, $index, $cn, $chr, $start, $end, $conf) = @field;
		my ($startsnp, $endsnp) = ("NA", "NA");
		
		unless ($chr =~ m/^\d+$/ and $chr>=1 and $chr <= 22) {
			$count_nonautosome++;
			next;
		}
		
		if ($snplocfile) {
			$startsnp = $pos_name->{$chr, $start} || 'NA';
			$endsnp = $pos_name->{$chr, $end} || 'NA';
		}
		
		if ($chr =~ m/^\d+$/) {
			if ($chr >= 1 and $chr <= 22) {
				$cn eq '2' and next;
			}
		}
		
		printPennCNVFormat ($chr, $start, $end, $cn, $sample, $startsnp, $endsnp, $conf, 0);
	}
	$count_nonautosome and print STDERR "WARNING: Total of $count_nonautosome non-autosome CNV calls are skipped for conversion\n";
	print STDERR "WARNING: the numsnp field is set as 0 for all calls, since Birdsmerge output does not contain such information\n";
}

sub convertCanaryToPennCNV {
	my ($fh, $canarydeffile) = @_;
	
	$verbose and print STDERR "NOTICE: Converting Canary CNV calls to PennCNV format\n";
	
	my (%loc);
	
	open (DEF, $canarydeffile) or confess "Error: cannot read from canarydeffile: $!\n";
	$_ = <DEF>;
	s/[\r\n]+$//;
	m/^cnp_id\tchr\tstart\tend/ or confess "Error: invalid header line found in the canarydeffile ('cnp_id chr start end' expected): <$_>\n";
	while (<DEF>) {
		s/[\r\n]+$//;
		my @field = split (/\s+/, $_);
		$loc{$field[0]} = [$field[1], $field[2], $field[3]];
	}

	$_ = <$fh>;
	s/[\r\n]+$//;
	$_ =~ m/^cnp_id\t(.+)/ or confess "Error: invalid first line in the Canary call file: <$_>\n";
	my @sampleid = split (/\t/, $1);
	print STDERR "NOTICE: A total of ", scalar (@sampleid), " samples are included in the Canary file\n";
	
	my $count_nonautosome = 0;
	while (<$fh>) {
		s/[\r\n]+$//;
		my @field = split (/\t/, $_);
		@field == @sampleid+1 or confess "Error: invalid record found in Canary call file (expect to see ${\(scalar @sampleid)} CNV calls per line): <$_>\n";
		my ($cnpid, @cn) = @field;
		my ($startsnp, $endsnp) = ("NA", "NA");
		my ($conf) = ("NA");
		my ($chr, $start, $end) = @{$loc{$cnpid}};
		
		unless ($chr =~ m/^\d+$/ and $chr>=1 and $chr <= 22) {
			$count_nonautosome++;
			next;
		}
		

		for my $i (1 .. @field-1) {
			my $cn = $field[$i];
			if ($chr =~ m/^\d+$/) {
				if ($chr >= 1 and $chr <= 22) {
					$cn eq '2' and next;
				}
			}
			printPennCNVFormat ($chr, $start, $end, $cn, $sampleid[$i-1], $startsnp, $endsnp, $conf, 0);
		}
	}
	$count_nonautosome and print STDERR "WARNING: Total of $count_nonautosome non-autosome CNV calls are skipped for conversion\n";
	print STDERR "WARNING: the numsnp field is set as 0 for all calls, since Canary output does not contain such information\n";
}

sub convertTabToPennCNV {
	my ($infh) = @_;
	$verbose and print STDERR "NOTICE: Converting tab-delimited CNV calls to PennCNV format\n";
	
	while (<$infh>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		@record == 9 or confess "Error: invalid record found in tab-delimited CNV file (9 fields expected): <$_>\n";
		my ($chr, $start, $end, $cn, $sample, $startsnp, $endsnp, $conf, $numsnp) = @record;

		printPennCNVFormat ($chr, $start, $end, $cn, $sample, $startsnp, $endsnp, $conf, $numsnp);
	}
}

sub convertPennCNVToTab {
	my ($infh) = @_;
	$verbose and print STDERR "NOTICE: converting the CNV call file $cnvfile to tab-delimited format (Chr, Start, End, CN, Sample, StartSNP, EndSNP, Conf, NumSNP)\n";

	my $skipped_line;
	while (<$infh>) {
		if (m/^chr((\w+):(\d+)-(\d+))\s+numsnp=(\d+)\s+length=([\d\,]+)\s+state\d+,cn=(\d+)\s+(.+?)\s+startsnp=(\S+)\s+endsnp=(\S+)/) {
			print $2, "\t", $3, "\t", $4, "\t", $7, "\t", $8, "\t", $9, "\t", $10;
			my $numsnp = $5;
			if (m/conf=(\S+)/) {
				print "\t$1";
			} else {
				print "\tNA";
			}
			print "\t$numsnp\n";
		} else {
			$skipped_line++;
		}
	}
	$skipped_line and print STDERR "WARNING: $skipped_line lines in $cnvfile are skipped due to unrecognizable format\n";
}



sub convertXMLToPennCNV {
	my ($infh, $snplocfile) = @_;
	$verbose and print STDERR "NOTICE: Converting XML-formatted CNV calls (exported from BeadStudio/GenomeStudio) to PennCNV format\n";
	
	my $pos_name;
	
	$_ = <$infh>;
	defined $_ or confess "Error: cannot read anything from input\n";
	s/[\r\n]+$//;
	m/^<Project_Bookmarks>/ or confess "Error: invalid first line in XML file (<Project_Bookmarks> expected): <$_>\n";
	
	if ($snplocfile) {
		$pos_name = readSNPLocFile ($snplocfile);
	}
	
	my $countline = 0;
	while (<$infh>) {
		$countline++;
		if (m/^\s*<bookmark>/) {
			my ($sample_id, $chr_num, $cn, $base_start_pos, $base_end_pos, $conf);
			my ($startsnp, $endsnp) = qw/NA NA/;
			while (<$infh>) {
				if (m#^\s*<sample_id>(\S+)\s*\[\d+\]</sample_id>#) {
					$sample_id = $1;
				} elsif (m#\s*<bookmark_type>CNV Bin: Min (\S+) To Max (\S+)</bookmark_type>#) {
					$cn = int(($1+$2)/2);
				} elsif (m#^\s*<chr_num>(\w+)</chr_num>#) {
					$chr_num = $1;
				} elsif (m#^\s*<base_start_pos>(\d+)</base_start_pos>#) {
					$base_start_pos = $1;
				} elsif (m#^\s*<base_end_pos>(\d+)</base_end_pos>#) {
					$base_end_pos = $1;
				} elsif (m/CNV Confidence: (\S+)/) {
					$conf = $1;
				} elsif (m/^\s*<bookmark>/) {
					confess "Error: Reading a new Bookmark before finished reading the previous one at line $countline\n";
				}
				if (defined $sample_id and defined $chr_num and defined $base_start_pos and defined $base_end_pos and defined $conf) {
					if ($snplocfile) {
						$startsnp = $pos_name->{$chr_num, $base_start_pos} || 'NA';
						$endsnp = $pos_name->{$chr_num, $base_end_pos} || 'NA';
					}
					printPennCNVFormat ($chr_num, $base_start_pos, $base_end_pos, $cn, $sample_id, $startsnp, $endsnp, $conf, 0);
					last;
				}
			}
		}
	}
}

sub readSNPLocFile {
	my ($snplocfile) = @_;
	my (%pos_name);
	my ($header, $header_seg, $lrr_index, $pos_index, $chr_index, $name_index);
	
	open (SIG, $snplocfile) or confess "Error: cannot read from SNPLOC file $snplocfile: $!\n";

	$header = <SIG>;
	$header =~ m/(.+)Pos/ or confess "error: the header line of $snplocfile does not contain Pos annotation";
	$header_seg = $1;
	$pos_index = ($header_seg =~ tr/\t/\t/+0);
	$header =~ m/(.+)Chr/ or confess "error: the header file of $snplocfile does not contain Chr annotation";
	$header_seg = $1;
	$chr_index = ($header_seg =~ tr/\t/\t/+0);
	$header =~ m/(.*)Name/ or confess "error: the header file of $snplocfile does not contain Name annotation";
	$header_seg = $1;
	$name_index = ($header_seg =~ tr/\t/\t/+0);


	while (<SIG>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		my ($curname, $curchr, $curpos) = @record[$name_index, $chr_index, $pos_index];
		$pos_name{$curchr, $curpos} = $curname;
	}
	close (SIG);
	return (\%pos_name);
}


=head1 SYNOPSIS

 convert_cnv.pl [arguments] <cnvfile>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --intype <string>		input file format (default: penncnv)
 	    --outtype <strong>		output file format (default: tab)
 	    --output <file>		specify output file name (default: STDOUT)
 	    --snplocfile <file>		file annotating SNP locations (for Birdseye/XML format)
 	    --canarydeffile <file>	file containing CNP locations (for Canary calls)

 Function: convert formats between CNV calls generated from different 
 programs/algorithms

 Example: convert_cnv.pl -intype penncnv -outtype tab ex1.rawcnv > ex1.tabcnv
          convert_cnv.pl -intype birdseye -outtype penncnv -output ex1.penncnv -snploc gw6.hg18.pfb ex1.birdcnv 
          convert_cnv.pl -intype xml -outtype penncnv ex1.xml > ex1.penncnv
          convert_cnv.pl -intype canary -outtype penncnv -canarydef GenomeWideSNP_6.hg18.cnv_defs case.canary_calls 
          convert_cnv.pl -intype oncosnp -outtype penncnv input.cnv > output.penncnv

 Version: $LastChangedDate: 2011-06-16 22:13:39 -0700 (Thu, 16 Jun 2011) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--intype>

the input type can be penncnv, tab, xml (exported from BeadStudio/GenomeStudio), birdseye, canary

=item B<--outtype>

the outtype can be penncnv and tab.

=item B<--output>

specify the output file name (default: STDOUT)

=item B<--snplocfile>

specify a file that annotate chromosome coordinates for SNPs. This is useful for 
converting birdseye calls since the output does not give marker names.

=item B<--snplocfile>

specify a file that annotate chromosome coordinates for CNPs for Canary calls. 
The four tab-delimited fields are CNP ID, Chr, Start and End.

=back

=head1 DESCRIPTION

This program is used to convert copy number variation (CNV) calls between 
different formats. Currently, the input format includes penncnv, tab (penncnv-
specific tab output), xml (generated from BeadStudio/GenomeStudio) or birdseye.

Various software tools have been developed in recent years to generate CNV calls 
from various sources of data, while the basic information in CNV calls from each 
software is generally similar: usually the chromosome position is given, the 
input file name is given, the start and end marker names are given, and the 
confidence score for each CNV call is given. This program operates on these 
specific fields and convert the file formats to facilitate comparative analysis. 
For example, it can convert the PennCNV calls to tab-delimited formats to be 
loaded into Excel or other software that operates on tab-delimited data; it can 
also convert CNV calls generated by other software to PennCNV format, so that 
the compare_cnv.pl program can be used to perform comparative analysis on CNV 
call concordance.

=head2 PennCNV format

The PennCNV output format is a simple text format with one CNV call per line. An 
example is given below:

	chr3:37957465-37961253        numsnp=3      length=3,789       state2,cn=1 father.txt startsnp=rs9837352 endsnp=rs9844203 conf=15.133
	chr3:75511365-75650909        numsnp=7      length=139,545     state2,cn=1 father.txt startsnp=rs4677005 endsnp=rs2004089 conf=26.862
	chr11:81181640-81194909       numsnp=9      length=13,270      state2,cn=1 father.txt startsnp=rs7947005 endsnp=rs12293984 conf=35.083
	chr20:10440279-10511908       numsnp=10     length=71,630      state2,cn=1 father.txt startsnp=rs8114269 endsnp=rs682562 conf=39.410
	chr11:55127597-55204003       numsnp=11     length=76,407      state2,cn=1 mother.txt startsnp=rs2456022 endsnp=rs7934845 conf=23.106

The fields are chromosome coordinates, the number of markers (SNPs or CN 
markers), the length of the CNV call, the copy number (CN) estimate, the signal 
file name, the start marker, the end marker and the confidence score for CNV call.

=head2 Tab format

The Tab output format is a penncnv-specific format that cooresponds to the 
same fields in PennCNV output, albeit being tab-delimited. An example is given 
below:

	3       37957465        37961253        1       father.txt      rs9837352       rs9844203       15.133  3
	3       75511365        75650909        1       father.txt      rs4677005       rs2004089       26.862  7
	11      81181640        81194909        1       father.txt      rs7947005       rs12293984      35.083  9
	20      10440279        10511908        1       father.txt      rs8114269       rs682562        39.410  10
	11      55127597        55204003        1       mother.txt      rs2456022       rs7934845       23.106  11

=head2 XML format

The XML format refers to a specific format exported from the Illumina 
BeadStudio/GenomeStudio software. An example of the first a few lines and part 
of the file are given below:

	<Project_Bookmarks>
	  <Version>2.0.0</Version>
	  <Name>PennCNV 2009Aug27</Name>
	  <Author></Author>
	  <Comment><![CDATA[  ]]></Comment>
	  <CreateDate>8/26/2009 5:14:38 PM</CreateDate>
	  <Algorithm></Algorithm>
	  ...........
	  ...........
	  ...........
	    <bookmark>
	      <sample_id>99HI0697A [1]</sample_id>
	      <bookmark_type>CNV Bin: Min 0.5 To Max 1.5</bookmark_type>
	      <entry_date>5:14:38 PM</entry_date>
	      <chr_num>2</chr_num>
	      <base_start_pos>242565979</base_start_pos>
	      <base_end_pos>242656041</base_end_pos>
	      <author></author>
	      <value>1</value>
	      <comment>
	        <![CDATA[ CNV Confidence: 52.304  ]]>
	      </comment>
	    </bookmark>

Note that the start and end marker information is not provided in the XML output 
file. Therefore, one must use the --snplocfile argument to supply this 
information, otherwise the marker will be annotated as NA.

=head2 Birdseye format

This format is a birdseye-specific tab-delimited format generated by the 
BirdSuite software. The columns are sample, sample_index, copy_number, chr, 
start, end, per_probe_score, size, num_probes, lod_score, respectively. Similar 
to the XML format above, since the marker name is not annotaed in the file, one 
must use --snplocfile to supply this information. Additionally, since the 
birdseye file contain both CNV information as well as normal copy information, 
only the autosomal CNVs will be converted. (For sex chromosome CNVs, there is no 
way to tell whether it is a real CNV or not, based on CN estimate only, in the 
Birdseye output)

=head2 Birdmerge format

This format is used to represent that merged calls from Birdseye and 
Canary output files. The columns are sample, sample_index, copy_number, chr, 
start, end, confidence, respectively. Since the marker name is not annotaed in 
the file, one must use --snplocfile to supply this information. Additionally, 
since the birdseye file contain both CNV information as well as normal copy 
information, only the autosomal CNVs will be converted. (For sex chromosome 
CNVs, there is no way to tell whether it is a real CNV or not, based on CN 
estimate only, in the Birdmerge output)

sample	sample_index	copy_number	chr	start	end	confidence
NA06991.affy6.cel	1	2	1	51598	17063449	8576.0
NA06991.affy6.cel	1	3	1	17067742	17134834	10.0
NA06991.affy6.cel	1	2	1	17148828	25455928	5135.0
NA06991.affy6.cel	1	1	1	25465715	25534592	10.0
NA06991.affy6.cel	1	2	1	25534593	48052671	12098.0
NA06991.affy6.cel	1	1	1	48054319	48055107	0.62
NA06991.affy6.cel	1	2	1	48055171	61885690	8806.0
NA06991.affy6.cel	1	0	1	61886594	61890775	10.0
NA06991.affy6.cel	1	2	1	61892397	72522941	7032.0
NA06991.affy6.cel	1	1	1	72528701	72583736	60.81
NA06991.affy6.cel	1	2	1	72584492	105814899	20505.0

=head2 Canary format

This format is produced by Canary to infer CN numbers on common CNP regions. One 
example is given below:

cnp_id  NA06991.affy6.cel       NA06994.affy6.cel       NA07000.affy6.cel       NA07019.affy6.cel       NA07034.affy6.cel
CNP10000        2       3       2       2       2
CNP10010        2       2       2       3       2
CNP10026        2       4       2       2       2
CNP10028        2       1       1       2       2


For questions, comments or bug reports, please contact me at 
kai@openbioinformatics.org.

=back

