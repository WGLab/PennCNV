#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2010-06-29 14:13:11 -0700 (Tue, 29 Jun 2010) $';

our ($verbose, $help, $man);
our ($inputfile);
our ($output, $numsnp, $length, $maxnumsnp, $maxlength, $type, $confidence, $maxconfidence, $coversnp, $signalfile, $qclogfile, $qclrrsd, $qcbafdrift, $qcwf, $qcnumcnv, $qcpassout, $qcsumout, 
	$maxtotalcnvlength, $mintotalcnvlength, $chrx, $chroms, $nochroms);
our (%chr_include, %chr_exclude);

GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'output=s'=>\$output, 'numsnp=i'=>\$numsnp, 'length=s'=>\$length, 'maxnumsnp=i'=>\$maxnumsnp,
	'maxlength=s'=>\$maxlength, 'confidence=f'=>\$confidence, 'maxconfidence=f'=>\$maxconfidence, 'coversnp=i'=>\$coversnp, 'signalfile=s'=>\$signalfile, 'type=s'=>\$type,
	'qclogfile=s'=>\$qclogfile, 'qclrrsd=f'=>\$qclrrsd, 'qcbafdrift=f'=>\$qcbafdrift, 'qcwf=f'=>\$qcwf, 'qcnumcnv=i'=>\$qcnumcnv, 
	'qcpassout=s'=>\$qcpassout, 'qcsumout=s'=>\$qcsumout, 'maxtotalcnvlength=s'=>\$maxtotalcnvlength, 'mintotalcnvlength=s'=>\$mintotalcnvlength, 'chrx'=>\$chrx, 
	'chroms=s'=>\$chroms, 'nochroms=s'=>\$nochroms) or pod2usage ();


$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error:");

($inputfile) = @ARGV;

my ($marker_pos, $file_qc, $file_numcnv, $file_sum, $file_cnvlen);

if (defined $length) {
	$length =~ s/m$/000000/;
	$length =~ s/k$/000/;
	$length =~ m/^\d+$/ or pod2usage ("Error in argument: --length must be a positive inteer (suffix of k and m are okay)");
}

if (defined $maxlength) {
	$maxlength =~ s/m$/000000/;
	$maxlength =~ s/k$/000/;
	$maxlength =~ m/^\d+$/ or pod2usage ("Error in argument: --maxlength must be a positive inteer (suffix of k and m are okay)");
}

if (defined $maxtotalcnvlength) {
	$maxtotalcnvlength =~ s/m$/000000/;
	$maxtotalcnvlength =~ s/k$/000/;
	$maxtotalcnvlength =~ m/^\d+$/ or pod2usage ("Error in argument: --maxtotalcnvlength must be a positive inteer (suffix of k and m are okay)");
}
if (defined $mintotalcnvlength) {
	$mintotalcnvlength =~ s/m$/000000/;
	$mintotalcnvlength =~ s/k$/000/;
	$mintotalcnvlength =~ m/^\d+$/ or pod2usage ("Error in argument: --mintotalcnvlength must be a positive inteer (suffix of k and m are okay)");
}

if (defined $type) {
	$type eq 'del' or $type eq 'dup' or pod2usage ("Error in argument: --type must be either del or dup");
}

if ($coversnp) {
	defined $signalfile or pod2usage ("Error in argument: please specify the --signalfile argument when using --coversnp argument");
	$marker_pos = readMarkerPos ($signalfile);
}

if (defined $qclrrsd or defined $qcbafdrift or defined $qcwf) {
	$qclogfile or pod2usage ("Error in argument: please specify the --qclogfile argument for QC based CNV filtering");
}

if (defined $chroms) {
	defined $nochroms and pod2usage ("Error in argument: please do not specify --chroms and --nochroms together");
	my @element = split (/,/, $chroms);
	for my $element (@element) {
		if ($element =~ m/^(\d+)$/) {
			$chr_include{$1} = 1;
		} elsif ($element =~ m/^(\d+)\-(\d+)$/) {
			$1 < $2 or pod2usage ("Error in argument: the element $element in the --chroms argument '$chroms' is invalid (must start from smaller number to larger number)");
			for ($1 .. $2) {
				$chr_include{$_} = 1;
			}
		} elsif ($element =~ m/^([a-zA-Z]+)$/) {
			$chr_include{$1} = 1;
		} else {
			pod2usage ("Error in argument: the element $element in the --chroms argument '$chroms' is invalid");
		}
	}
}

if (defined $nochroms) {
	defined $chroms and pod2usage ("Error in argument: please do not specify --chroms and --nochroms together");
	my @element = split (/,/, $nochroms);
	for my $element (@element) {
		if ($element =~ m/^(\d+)$/) {
			$chr_exclude{$1} = 1;
		} elsif ($element =~ m/^(\d+)\-(\d+)$/) {
			$1 < $2 or pod2usage ("Error in argument: the element $element in the --chroms argument '$chroms' is invalid (must start from smaller number to larger number)");
			for ($1 .. $2) {
				$chr_exclude{$_} = 1;
			}
		} elsif ($element =~ m/^([a-zA-Z]+)$/) {
			$chr_exclude{$1} = 1;
		} else {
			pod2usage ("Error in argument: the element $element in the --chroms argument '$chroms' is invalid");
		}
	}
}

if ($qclogfile) {
	defined $qclrrsd or print STDERR "NOTICE: the --qclrrsd argument is set as 0.3 by default\n" and $qclrrsd = 0.3;
	defined $qcbafdrift or print STDERR "NOTICE: the --qcbafdrift argument is set as 0.01 by default\n" and $qcbafdrift = 0.01;
	defined $qcwf or print STDERR "NOTICE: the --qcwf argument is set as 0.05 by default\n" and $qcwf = 0.05;
	open (LOG, $qclogfile) or confess "Error: cannot read from qclogfile $qclogfile: $!\n";
	while (<LOG>) {
		if (not $chrx and m/([^"`\s]+|"[^"]+"|`[^`]+`): LRR_mean=(\S+) LRR_median=(\S+) LRR_SD=(\S+) BAF_mean=(\S+) BAF_median=(\S+) BAF_SD=(\S+) BAF_DRIFT=(\S+) WF=(\S+)/) {	#"
			my $abs_wf = $9;
			$abs_wf eq 'NA' and $abs_wf = 0;
			$abs_wf = abs ($abs_wf);
			if ($4 <= $qclrrsd and $8 <= $qcbafdrift and $abs_wf <= $qcwf) {
				$file_qc->{$1} = "pass";
			} else {
				$file_qc->{$1} = "fail";
				$verbose and print STDERR "WARNING: File $1 failes QC due to LRR_SD=$4, BAF_DRIFT=$8 and WF=$9\n";
			}
			$file_sum->{$1} = [$2, $3, $4, $5, $6, $7, $8, $9];
		} elsif ($chrx and m/([^"`\s]+|"[^"]+"|`[^`]+`): LRR_Xmean=(\S+) LRR_Xmedian=(\S+) LRR_XSD=(\S+) BAF_Xhet=(\S+)/) {	#"
			my $abs_wf = 0;
			if ($4 <= $qclrrsd) {
				$file_qc->{$1} = "pass";
			} else {
				$file_qc->{$1} = "fail";
				$verbose and print STDERR "WARNING: File $1 failes QC due to LRR_XSD=$4\n";
			}
			$file_sum->{$1} = [$2, $3, $4, $5];
		}
	}
}

if ($qcpassout or $qcsumout) {
	$qclogfile or pod2usage ("Error in argument: please specify --qclogfile for the --qcpassout or the --qcsumout argument");
}

if ($qcnumcnv or $qcsumout or $mintotalcnvlength or $maxtotalcnvlength) {
	if ($inputfile eq 'stdin') {
		open (CNV, ">STDIN.filtercnv") or confess "Error: the --qcnumcnv requires writting a temporary file STDIN.filtercnv but fails: $!\n";
		while (<STDIN>) {
			print CNV $_;
		}
		close (CNV);
		open (CNV, "STDIN.filtercnv") or confess "Error: cannot read from inputfile STDIN.filtercnv: $!\n";
	} else {
		open (CNV, $inputfile) or confess "Error: cannot read from inputfile $inputfile: $!\n";
	}
	while (<CNV>) {
		if (m/^(?:chr)?(\w+):(\d+)-(\d+)\s+numsnp=(\d+)\s+length=(\S+)\s+state(\d),cn=(\d)\s+(.+?)\s+startsnp=(\S+) endsnp=(\S+)/) {
			my ($len, $file) = ($5, $8);
			$file_numcnv->{$file}++;
			$len =~ s/,//g;
			$file_cnvlen->{$file} += $len;
			
		}
	}
	close (CNV);
	
	if ($qclogfile and defined $mintotalcnvlength || defined $maxtotalcnvlength) {		#check the QC based on max/min totalcnvlength
		for my $key (keys %$file_cnvlen) {
			if (defined $mintotalcnvlength) {
				$file_cnvlen->{$key} < $mintotalcnvlength and $file_qc->{$key} = 'fail';
			}
			if (defined $maxtotalcnvlength) {
				$file_cnvlen->{$key} > $maxtotalcnvlength and $file_qc->{$key} = 'fail';
			}
		}
	}
}



if (defined $output) {	#when --nfperfile is set, we need to write to a file (rather than STDOUT)
	open (STDOUT, ">$output") or confess "Error: cannot write to output file $output: $!\n";
}
filterCNV ($inputfile, $marker_pos, $file_qc, $file_sum, $file_numcnv);

sub filterCNV {
	my ($cnvfile, $marker_pos, $file_qc, $file_sum, $file_numcnv) = @_;
	my ($curchr, $curstart, $curend, $curnumsnp, $curlength, $curstate, $curcn, $curfile, $cursnpstart, $cursnpend);
	my ($skipped_line);
	my (%file_wo_qc);		#files without QC measures in the qclogfile
	
	if ($cnvfile eq 'stdin') {
		if ($qcnumcnv or $qcsumout or $mintotalcnvlength or $maxtotalcnvlength) {
			open (CNV, "STDIN.filtercnv") or confess "Error: cannot read from inputfile STDIN.filtercnv: $!\n";
		} else {
			*CNV = *STDIN;
		}
	} else {
		open (CNV, $cnvfile) or confess "Error: cannot read from inputfile $cnvfile: $!\n";
	}
	while (<CNV>) {
		if (m/^(?:chr)?(\w+):(\d+)-(\d+)\s+numsnp=(\d+)\s+length=(\S+)\s+state(\d),cn=(\d)\s+(.+?)\s+startsnp=(\S+) endsnp=(\S+)/) {
			($curchr, $curstart, $curend, $curnumsnp, $curlength, $curstate, $curcn, $curfile, $cursnpstart, $cursnpend) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10);
			$curlength =~ s/,//g;
			
			if (defined $qclogfile) {
				if (not exists $file_qc->{$curfile}) {
					$file_wo_qc{$curfile}++;
					next;
				} elsif ($file_qc->{$curfile} eq "fail") {
					next;
				}
			}
			
			if (%chr_include or %chr_exclude) {		#if --chroms or --nochroms are specified, check the chromosome of the CNV call to satisfy the criteria
				%chr_include and $chr_include{$curchr} || next;
				%chr_exclude and $chr_exclude{$curchr} && next;
			}
			
			if ($qcnumcnv) {
				if (exists $file_numcnv->{$curfile} and $file_numcnv->{$curfile} > $qcnumcnv) {
					next;
				}
			}
			
			if (defined $maxtotalcnvlength) {
				if ($file_cnvlen->{$curfile} > $maxtotalcnvlength) {	#the file has total CNV length larger than a pre-specified threshold
					next;
				}
			}
			if (defined $mintotalcnvlength) {
				if ($file_cnvlen->{$curfile} < $mintotalcnvlength) {	#the file has total CNV length larger than a pre-specified threshold
					next;
				}
			}
			
			defined $length and $curlength < $length and next;
			defined $numsnp and $curnumsnp < $numsnp and next;
			defined $maxlength and $curlength > $maxlength and next;
			defined $maxnumsnp and $curnumsnp > $maxnumsnp and next;
			
			if (defined $type) {
				if ($type eq 'del') {
					$curcn < 2 or next;
				} elsif ($type eq 'dup') {
					$curcn >= 2 or next;		#male chrY duplication, male chrX duplication, autosome duplication
				}
			}
			
			if (defined $confidence) {
				s/conf=nan/conf=-1/;
				s/conf=-Inf/conf=-1000000/;
				m/conf=([\d\.\-]+)/ or pod2usage ("Error in argument: the --confidence argument can only be applied to CNV calls with confidence annotations (unrecognized call: <$_>)");
				$1 < $confidence and next;
			}
			
			if (defined $maxconfidence) {
				s/conf=nan/conf=-1/;
				s/conf=-Inf/conf=-1000000/;
				m/conf=([\d\.\-]+)/ or pod2usage ("Error in argument: the --confidence argument can only be applied to CNV calls with confidence annotations (unrecognized call: <$_>)");
				$1 > $maxconfidence and next;
			}
			
			if ($marker_pos) {		#this is used to calculate the number of covered SNPs per individual
				my $curdata = $marker_pos->{$curchr};
				my $cover_count = 0;
				$curdata or next;
				for my $i (0 .. @$curdata-1) {
					if ($curdata->[$i][0] < $curstart) {
						next;
					} elsif ($curdata->[$i][0] <= $curend) {
						$cover_count++;
					} else {
						last;
					}
				}
				$cover_count < $coversnp and next;
			}
				
			
		} else {
			$skipped_line++;
		}
		print;
	}
	$qcnumcnv and $inputfile eq 'stdin' and unlink ("STDIN.filtercnv");
	$skipped_line and print STDERR "WARNING: $skipped_line lines in CNV file $cnvfile are skipped due to unrecognizable format\n";
	$qclogfile and %file_wo_qc and print STDERR "WARNING: CNV calls from ${\(scalar keys %file_wo_qc)} files are skipped due to lack of QC measure in qclogfile $qclogfile\n";

	if ($qcpassout) {
		writeQCPass ($qcpassout, $file_qc, $file_sum, $file_numcnv);
	}
	
	if ($qcsumout) {
		writeQCSum ($qcsumout, $file_qc, $file_sum, $file_numcnv);
	}
}

sub writeQCPass {
	my ($qcpassout, $file_qc, $file_sum, $file_numcnv) = @_;
	my $count_qcpass = 0;
	open (QCPASS, ">$qcpassout") or confess "Error: cannot write qcpass file $qcpassout: $!\n";
	for my $file (keys %$file_qc) {
		if ($file_qc->{$file} eq "pass") {
			if ($qcnumcnv) {
				if (exists $file_numcnv->{$file} and $file_numcnv->{$file} > $qcnumcnv) {
					next;
				}
			}
			print QCPASS $file, "\n";
			$count_qcpass++;
		}
	}
	close (QCPASS);
	print STDERR "NOTICE: Writting $count_qcpass file names that pass QC to qcpass file $qcpassout\n";
}

sub writeQCSum {
	my ($qcsumout, $file_qc, $file_sum, $file_numcnv) = @_;
	my $count_qcsum = 0;
	open (QCSUM, ">$qcsumout") or confess "Error: cannot write to qcsum file $qcsumout: $!\n";
	if ($chrx) {
		print QCSUM "File\tLRR_Xmean\tLRR_Xmedian\tLRR_XSD\tBAF_Xhet\tNumCNV\n";
	} else {
		print QCSUM "File\tLRR_mean\tLRR_median\tLRR_SD\tBAF_mean\tBAF_median\tBAF_SD\tBAF_drift\tWF\tNumCNV\n";
	}
	for my $file (keys %$file_qc) {
		print QCSUM $file, "\t", join ("\t", @{$file_sum->{$file}}), "\t", $file_numcnv->{$file} || 0, "\n";
		$count_qcsum++;
	}
	close (QCSUM);
	print STDERR "NOTICE: Writting $count_qcsum records of QC summary to qcsum file $qcsumout\n";
}


sub readMarkerPos {
	my ($signalfile) = @_;
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
	my $count;
	my (%snp_chr, %snp_pos);
	while (<SIG>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		my ($curname, $curchr, $curpos) = @record[$name_index, $chr_index, $pos_index];
		push @{$alldata{$curchr}}, [$curpos, $curname];
		$snp_chr{$curname} = $curchr;
		$snp_pos{$curname} = $curpos;
		$count++;
	}
	close (SIG);
	for my $key (keys %alldata) {
		@{$alldata{$key}} = sort {$a->[0] <=> $b->[0]} @{$alldata{$key}};
	}
	print STDERR "NOTICE: Total of $count records are read from $signalfile\n";
	return \%alldata;
}

=head1 SYNOPSIS

 filter_cnv.pl [arguments] <cnv-call-file>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	
 	    --output <file>		write output to this file (default=STDOUT)
 	    --numsnp <int>		minimum number of SNPs in CNV calls
 	    --length <int>		minimum length of CNV calls (suffix of k or m is okay)
 	    --maxnumsnp <int>		maximum number of SNPs in CNV calls
 	    --maxlength <int>		maximum length of CNV calls (suffix of k or m is okay)
 	    --type <dup|del>		type of CNVs (dup or del)
 	    --confidence <float>	minimum confidence score of CNV calls
 	    --maxconfidence <float>	maximum confidence score of CNV calls
 	    --coversnp <int>		minimum covered SNPs (in another signal file)
 	    --signalfile <file>		signal intensity file for checking coverage
 	    
 	    --qclogfile <file>		the log file with QC information for all samples
 	    --qclrrsd <float>		LRR_SD threshold for inclusion of the sample
 	    --qcbafdrift <float>	BAF_DRIFT threshold for inclusion of the sample
 	    --qcwf <float>		WF threshold for inclusion of the sample
 	    --qcnumcnv <int>		number of CNV call for inclusion of the sample
 	    --qcpassout <file>		write samples passing QC threshold to this file
 	    --qcsumout <file>		write QC summary for ALL samples to this file
 	    --chrx			process QC for chrX CNV calls (default: QC on autosome only)
 	    --mintotalcnvlength <int>	minimum total CNV length for an individual to be included in output
 	    --maxtotalcnvlength <int>	maximum total CNV length for an individual to be included in output
 	    --chroms <string>		a list of chromosomes to include in output (eg: 1-5,10)
 	    --nochroms <string>		a list of chromosomes to exclude in output (eg: 3-7,8,10,12,14)

 Function: filter CNV calls from input by various criteria

 Example: filter_cnv.pl cnvfile -numsnp 10 > newcnvfile
          filter_cnv.pl cnv_hh550 -coversnp 20 -signal affygw6.hg18.pfb -out new_cnvfile
          filter_cnv.pl cnvfile -qclogfile list.log -qclrrsd 0.28 -qcnumcnv 100 -qcpassout list.qcpass -qcsumout list.qcsum > newcnvfile

 Version: $LastChangedDate: 2010-06-29 14:13:11 -0700 (Tue, 29 Jun 2010) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--output>

specify the output file name (default is to print to STDOUT)

=item B<--numsnp>

specify the minimum number of SNPs in CNV calls

=item B<--length>

specify the minimum length of CNV calls

=item B<--maxnumsnp>

specify the maximum number of SNPs in CNV calls

=item B<--maxlength>

specify the maximum length of CNV calls

=item B<--confidence>

specify the minimum confidence score of CNV calls

=item B<--coversnp>

specify the mimimum number SNPs that can be covered by another signal file or PFB file

=item B<--signalfile>

a signal intensity file from another array for checking coverage of the specified CNVs

=item B<--qclogfile>

specify a log file generated by PennCNV for the automated QC procedure

=item B<--qclrrsd>

specify the maximum allowed LRR_SD measure for a sample to be included in output 
during the automated QC procedure

=item B<--qcbafdrift>

specify the maximum allowed BAF_drift measure for a sample to be included in output 
during the automated QC procedure

=item B<--qcwf>

specify the maximum allowed absolute value of WF measure for a sample to be 
included in output during the automated QC procedure. The WF measure could be 
positive or negative, so this argument uses the absolute value as a threshold.

=item B<--qcnumcnv>

specify the maximum allowed number of CNV calls for a sample, to include CNV 
calls for this sample in output. This argument can be used to exclude samples 
with exceedingly large number of CNV calls in subsequent analysis.

=item B<--qcpassout>

write all samples that passed QC procedure to this output file, with one file 
name per line.

=item B<--qcsumout>

write the QC summary measures for all samples to this output file, with 
information for one file per line. The summary information includes several QC 
measures such as LRR_SD, LRR_median, LRR_mean, etc. for each sample.

=item B<--chrx>

handle log files for chromosome X CNV calls (the log files for chromosome X 
differ from those for the autosomes)

=item B<--mintotalcnvlength>

the minimum total CNV length of an individual, for CNV calls on this individual 
to be included in output. A suffix of 'k' or 'm' is allowed.

=item B<--maxtotalcnvlength>

the maximum total CNV length of an individual, for CNV calls on this individual 
to be included in output. A suffix of 'k' or 'm' is allowed. This argument can 
be used to exclude samples with extremely large CNV calls (for example, trisomy 
21 samples) in the output.

=item B<--chroms>

specify a list of chromosomes to be included in output. It cannot be used with 
the --nochroms argument together. By default, CNV calls from all chromosomes 
will be included in output. The chromosome list needs to be comma separated, but 
a dash can be used to specify a range of chromosomes; for example, B<--chroms 1-3,4,6> 
is a valid argument. This argument has effects only on the output CNV 
calls, and it will not change the calculation of numCNV for each individual 
during the QC procedure.

=item B<--nochroms>

specify a list of chromosomes to be excluded in output. It cannot be used with 
the --chroms argument together. By default, CNV calls from all chromosomes will 
be included in output. The chromosome list needs to be comma separated, but a 
dash can be used to specify a range of chromosomes; for example, B<--nochroms 1-3,4,6> 
is a valid argument. This argument has effects only on the output CNV 
calls, and it will not change the calculation of numCNV for each individual 
during the QC procedure.

=back

=head1 DESCRIPTION

This program is designed for filter CNV calls generated by the PennCNV program 
by various criteria. Some major functionalities of this program are described in 
detail below:

=item * B<Filter CNV calls by specific criteria>

=over 4

The default parameter in PennCNV generate liberal CNV calls, but users may wish 
to impose strict criteria on the CNV calls for their own analysis. For example, 
it is generally useful to focus any analysis to CNV calls containing more than 
10 SNPs, although the default CNV calling parameter in PennCNV is 3 SNPs. In 
such cases, a "--numsnp 10" argument can be used in this program to identify a 
subset of CNV calls containing at least 10 SNPs.

Several other arguments, such as --length, --maxnumsnp, --maxlength, --type, --
confidence, --maxconfidence, are also provided and their argument names are 
self-evident. The --chroms and --nochroms argument can be used to selectively 
output CNVs in specific chromosomes, or not in specific chromosomes.

Additionally, if --qcnumcnv, --mintotalcnvlength or --maxtotalcnvlength is 
specified, the program will filter out all CNV calls for specific subjects, 
whose total number of CNV calls are higher than the threshold, or whose total 
CNV span (that is, length of all CNV calls added together) exceed the 
thresholds.

=back

=item * B<Identify a subset of CNV calls detectable by another technical platform>

=over 4

Another utility of this program is to identify a subset of CNV calls, that may 
potentially be detected by another technical platform. For example, a research may 
utilize Affymetrix arrays to detect CNVs for a specific subject, but this same 
person was also genotyped by Illumina arrays, so he/she wants to identify a 
subset of CNVs that are covered by at least 10 SNPs in the Illumin array, and 
then see how many of them can actually be detected in the Illumina array. 
Similarly, next-generation sequencing may detect thousands of CNVs in any given 
subject, but very often one wants to know whether a set of CNVs that cover at 
least 10 SNPs in Illumina arrays can be detected in the Illumina array to assess 
false negative rates of the Illumina arrays.

The following command can be used for this purpose

	filter_cnv.pl cnv_hh550 -coversnp 20 -signal affygw6.hg18.pfb -out new_cnvfile

the cnv_hh550 is the input file with CNV calls. Thhe --signal argument specify a 
random signal file for the other technical platform, which contains Name, Chr, 
Position in the columns. In some cases, users may provide a standard PFB file as 
the signal file, but  need to process the PFB file such that the file contains 
exactly the same markers as the desired "target" array.

=item * B<automatic quality control for CNV calls>

After running CNV calling, it is important to perform quality control measures 
to remove samples with low quality of signal intensity measures. Assuming that 
the log file is list.log and the CNV output file is cnvfile, the following 
command can be used.

	filter_cnv.pl cnvfile -qclogfile list.log -qclrrsd 0.28 -qcnumcnv 100 -qcpassout list.qcpass -qcsumout list.qcsum > newcnvfile

It will check all samples annotated in the cnvfile, then identify the subset of 
samples with LRR_SD less than 0.28, with number of CNV calls less than 100, then write 
this subset of subjects to the list.qcpass file. Additionally, it will write the 
QC summary measures for all samples to the list.qcsum file. Finally, it will 
print out CNV calls that pass the QC measures, and the output is written to the 
newcnvfile.

=back

For questions, comments or bug reports, please contact me at 
kai@openbioinformatics.org.

=cut