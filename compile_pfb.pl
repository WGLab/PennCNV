#!/usr/bin/perl -w
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2011-03-07 21:43:22 -0800 (Mon, 07 Mar 2011) $';

our ($verbose, $help, $man, $output, $snpposfile, $listfile);
our (@inputfile);

GetOptions ('verbose'=>\$verbose, 'help'=>\$help, 'man'=>\$man, 'output=s'=>\$output, 'snpposfile=s'=>\$snpposfile, 'listfile=s'=>\$listfile) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);

if ($listfile) {
	@ARGV and pod2usage ("Error in argument: please do not supply signal file names if the --listfile argument is specified");
	$listfile or pod2usage ("Error in argument: please either specify the --listfile argument or provide signal file names in command line\n");
	open (LIST, $listfile) or confess "\nERROR: cannot read from listfile $listfile: $!";
	while (<LIST>) {
		s/^\s+|\s*[\r\n]+$//g;		#get rid of leading and trailing spaces
		$_ or print STDERR "WARNING: blank lins in $listfile detected and skipped\n" and next;
		m/^([^\t]+)$/ or die "Error: the --listfile should contain one file name per line without tab character (invalid line found: <$_>)\n";
		push @inputfile, $1;
	}
	close (LIST);
	print STDERR "NOTICE: A total of ${\(scalar @inputfile)} input signal files is specified in $listfile\n";
} else {
	@ARGV > 0 or pod2usage ("Error in argument: please either specify the --listfile argument or provide signal file names in command line");
	(@inputfile) = @ARGV;
	print STDERR "NOTICE: A total of ${\(scalar @inputfile)} input signal files is specified in command line\n";
}

my (%filecount, @dupfile);
for my $nextfile (@inputfile) {
	$filecount{$nextfile}++;
	$filecount{$nextfile} == 2 and push @dupfile, $nextfile;
}
if (@dupfile) {
	my @exdupfile = splice (@dupfile, 0, 3);
	print STDERR "WARNING: Some signal files are provided multiple times (example: @exdupfile)\n";
}


my ($snp_chr, $snp_pos);
defined $snpposfile and ($snp_chr, $snp_pos) = readSNPPosFile ($snpposfile);

if ($output) {
	open (STDOUT, ">$output") or confess "Error: cannot write to output file $output: $!\n";
}

my (@fh, $header_seg, $baf_index, $chr_index, $pos_index);
for my $i (0 .. @inputfile-1) {
	my $fh;
	eval {
		open ($fh, $inputfile[$i]) or die $!;
	};
	if ($@) {
		if ($@ =~ m/Too many open files/) {
			print STDERR "NOTICE: File handle cannot be created by the operating system after reading ${\(scalar @fh)} files\n";
			last;
		} elsif ($@ =~ m/No such file or directory/) {
			print STDERR "WARNING: Skipping the file $inputfile[$i] that cannot be read by the current program\n";
			next;
		} else {
			print STDERR "WARNING: Skipping the file $inputfile[$i] because '$@'\n";
			next;
		}
	}
	
	push @fh, $fh;
	$_ = readline ($fh);
	s/[\r\n]+$//;
	m/(.*)B Allele Freq/ or confess "Error: invalid header line in signalfile $inputfile[$i] (B Allele Freq expected): <$_>\n";
	$header_seg = $1;
	if (not defined $baf_index) {
		$baf_index = ($header_seg =~ tr/\t/\t/+0);
	} else {
		$baf_index == ($header_seg =~ tr/\t/\t/+0) or die "Error: the column index for B Allele Freq is $baf_index in file ", $inputfile[$i-1], ", but becomes ", ($header_seg =~ tr/\t/\t/+0), " in file $inputfile[$i]\n";
	}
	if ($i == 0) {
		if (m/(.*)Chr/) {
			$header_seg = $1;
			$chr_index = ($header_seg =~ tr/\t/\t/+0);
			if (m/(.*)Pos/) {
				$header_seg = $1;
				$pos_index = ($header_seg =~ tr/\t/\t/+0);
			} else {
				$snpposfile or die "Error: Position is not annotated in header line in $inputfile[$i]. Use --snpposfile to supply this information.\n";
			}
		} else {
			$snpposfile or die "Error: Chr is not annotated in header line in $inputfile[$i]. Use --snpposfile to supply this information.\n";
		}			
	}
}
print STDERR "NOTICE: The B Allele Freq information is annotated as column ${\($baf_index+1)} in input files\n";
print STDERR "NOTICE: A total of ${\(scalar @fh)} input files will be used for compiling PFB values\n";



print "Name\tChr\tPosition\tPFB\n";

my ($count_noanno, $count_total) = (0, 0);
MAIN: while (1) {
	my (@baf, $snp, $chr, $pos);
	for my $i (0 .. @fh-1) {
		my $fh = $fh[$i];
		$_ = readline ($fh);
		defined $_ or last MAIN;
		s/[\r\n]+$//;
		my @field = split (/\t/, $_);
		
		$snp ||= $field[0];
		if (not defined $snpposfile) {
			$chr ||= $field[$chr_index];
			$pos ||= $field[$pos_index];
		}

		$field[$baf_index] eq 'NaN' and next;
		$field[$baf_index] eq 'NA' and next;
		$field[$baf_index] =~ m/^([\d\.\-eE]+)$/ or warn "WARNING: the BAF value ($field[$baf_index]) in file $inputfile[$i] is not a valid number in input line: <$_>\n" and next;	#sometimes decimal points become comma in the BAF files exported from BeadStudio
		push @baf, $field[$baf_index];
	}
	my $mean;
	if ($snp =~ m/^cnv/) {
		$mean = 2;
	} elsif (not @baf) {
		$mean = 2;
	} else {
		$mean = mean (\@baf);
		$mean = int($mean*1000+0.5)/1000;
	}
	if (defined $chr) {
		print "$snp\t$chr\t$pos\t$mean\n";
		$count_total++;
	} else {
		if (defined $snp_chr->{$snp}) {
			print "$snp\t$snp_chr->{$snp}\t$snp_pos->{$snp}\t$mean\n";
			$count_total++;
		} else {
			$count_noanno++;
		}
	}
}

print STDERR "NOTICE: PFB values for $count_total markers were written to output\n";
$count_noanno and print STDERR "WARNING: $count_noanno markers were not written in output due to lack of annotation in --snpposfile\n";


sub readSNPPosFile {
	my ($snpposfile) = @_;
	my (%name_pos, %name_chr);
	my ($header, $header_seg, $lrr_index, $pos_index, $chr_index, $name_index);
	
	open (SIG, $snpposfile) or confess "Error: cannot read from snpposfile $snpposfile: $!\n";
	print STDERR "NOTICE: Start reading snpposfile $snpposfile ...";
	$header = <SIG>;
	$header =~ m/(.+)Pos/ or confess "error: the header line of $snpposfile does not contain Pos annotation";
	$header_seg = $1;
	$pos_index = ($header_seg =~ tr/\t/\t/+0);
	$header =~ m/(.+)Chr/ or confess "error: the header file of $snpposfile does not contain Chr annotation";
	$header_seg = $1;
	$chr_index = ($header_seg =~ tr/\t/\t/+0);
	$header =~ m/(.*)Name/ or confess "error: the header file of $snpposfile does not contain Name annotation";
	$header_seg = $1;
	$name_index = ($header_seg =~ tr/\t/\t/+0);


	while (<SIG>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		my ($curname, $curchr, $curpos) = @record[$name_index, $chr_index, $pos_index];
		$name_chr{$curname} = $curchr;
		$name_pos{$curname} = $curpos;
	}
	close (SIG);
	print STDERR " Done with location information for ${\(scalar keys %name_chr)} markers\n";
	return (\%name_chr, \%name_pos);
}


sub mean {
	my ($score) = @_;
	@$score or confess "\nERROR: NO VALUES for calculating mean\n";
	my $sum;
	for (@$score) {
		$sum += $_;
	}
	return $sum/@$score;
}

=head1 SYNOPSIS

 compile_pfb.pl [arguments] <inputfiles | --listfile file>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --output <file>		specify output file (default: STDOUT)
 	    --snpposfile <file>		a file that contains Chr and Position for SNPs
 	    --listfile <file>		a listfile that contains signal file names

Function: compile PFB file from multiple signal intensity files containing BAF values

 Example: compile_pfb.pl input1 input2 input3 -output out.pfb
          compile_pfb.pl input1 input2 input3 -snpposfile affygw6.hg18.pfb -output out.pfb
          compile_pfb.pl -listfile signal_file_list -output out.pfb

 Version: $LastChangedDate: 2011-03-07 21:43:22 -0800 (Mon, 07 Mar 2011) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--output>

specify the output file name. By default, the output will be written to STDOUT.

=item B<--snpposfile>

specify a file that contains Chr and Position information for SNPs. Normally 
this information is provided in the input signal intensity file, but if they are 
not available, the users can provide the snpposfile. Markers that are not 
annotated in the snpposfile will NOT be printed in the output. Users could use a 
PFB file in PennCNV package for this argument.

=item B<--listfile>

specify a file that contains a list of signal file names, one per line. The 
program will process these files and calculate PFB values based on these files. 
When --listfile is specified, the signal file names cannot be additionally 
supplied in command line.

=back

=head1 DESCRIPTION

This program is used to compile a PFB file (used in PennCNV package for CNV 
calling), given a list of signal intensity files. The signal intensity file 
contains information for one marker per line, and the B Allele Freq column in 
the file will be used in the PFB calculation. It is recommended to ensure that 
these signal files represent the same ethnicity group before compiling a custom 
PFB file.

If the signal files contains the Chr and Position column, these columns will 
also show up in the output PFB file. Otherwise, the users are required to 
provided a --snpposfile specifying the position of SNPs. This file contains 
marker name, chr, position at each line for each marker. The PFB file included 
in the PennCNV package can be used as a --snpposfile (it may sound circular, but 
it is a convenient way for users to compile their own PFB files).

For questions, comments or bug reports, please contact me at 
kai@openbioinformatics.org.

=cut