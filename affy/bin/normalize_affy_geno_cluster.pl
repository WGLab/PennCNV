#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: 271 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2009-06-25 16:44:19 -0400 (Thu, 25 Jun 2009) $';

our ($verbose, $help, $man);
our ($locfile, $output, $power2);
our ($clusterfile, $sigfile);

GetOptions('verbose'=>\$verbose, 'help'=>\$help, 'man|m'=>\$man, 'locfile=s'=>\$locfile, 'output=s'=>\$output, 'power2!'=>\$power2) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error: <inputfile> missing");


($clusterfile, $sigfile) = @ARGV;
$output ||= "$sigfile.lrr_baf.txt";
defined $power2 or $power2 = 1;

my ($cluster, $snploc);

$cluster = readClusterFile ($clusterfile);
$locfile and $snploc = readLocFile ($locfile);

processSignalFile ($cluster, $sigfile, $output, $snploc);

sub processSignalFile {
	my ($cluster, $inputfile, $outputfile, $snploc) = @_;

	open (OUTPUT, ">$outputfile") or confess "Error: cannot write to outputfile $outputfile: $!";
	
	open (SIG, $inputfile) or confess "Error: cannot read from genotype file $inputfile: $!";
	print STDERR "NOTICE: Processing $inputfile...\n";
	while (<SIG>) {
		m/^#/ or last;
	}
	s/[\r\n]+$//;
	m/^probeset_id/ or confess "Error: unable to find the header line (which should start with 'probeset_id') for $inputfile";
	my @header = split (/\t/, $_);
	shift @header;
	my @output;
	
	#print the first line of the OUTPUT file
	print OUTPUT "Name";
	$snploc and print OUTPUT "\tChr\tPosition";
	for my $i (0 .. @header-1) {
		print OUTPUT "\t$header[$i].Log R Ratio\t$header[$i].B Allele Freq";
	}
	print OUTPUT "\n";
	
	my ($count_success, $count_nocluster, $count_noloc, $count_badchr) = qw/0 0 0 0/;
	while (<SIG>) {
		m/^#/ and next;
		s/[\r\n]+$//;

		my (@siga, @sigb, $sig_psida, $sig_psidb, $snp, $cn_flag);
		@siga = split (/\t/, $_);
		$sig_psida = shift @siga;
		if ($sig_psida =~ m/(\S+)\-A$/) {
			$snp = $1;
			$_ = <SIG>;
			s/[\r\n]+$//;
			@sigb = split (/\t/, $_);
			$sig_psidb = shift @sigb;
			$sig_psidb =~ m/$snp\-B$/ or confess "Error: expect to read B allele after reading A allele for $sig_psida, but instead read $sig_psidb";
		} else {
			$snp = $sig_psida;			#CN probe
			$cn_flag++;				#mark this probe as cn probe
		}
		
		if (not $cluster->{$snp}) {			#SNP or CN marker without cluster information
			$count_nocluster++;
			next;
		}
		
		if ($snploc) {
			$snploc->{$snp} or ++$count_noloc and next;
			$snploc->{$snp} =~ m/^(\-|0|MT)/ and ++$count_badchr and next;	#skip markers in these weird chromosomes
			print OUTPUT $snp, "\t", $snploc->{$snp};
		} else {
			print OUTPUT $snp;
		}
	
		if (@sigb) {		#SNP probe
			#transform signalA and signalB into exponential
			if ($power2) {
				@siga = map {2**$_} @siga;
				@sigb = map {2**$_} @sigb;
			}

			for my $i (0 .. @header-1) {
				my ($signala, $signalb) = ($siga[$i], $sigb[$i]);
				#print STDERR "NOTICE: Reading snp=$snp signala=$signala signalb=$signalb\n";
				my ($r, $theta) = ($signala+$signalb, atan2 ($signalb, $signala) / (3.1415926/2));
				my @snpinfo = split (/\t/, $cluster->{$snp});
				my @rmean = @snpinfo[0..2];
				my @thetamean = @snpinfo[3..5];
				my ($lrr, $baf);
				
				if ($theta < $thetamean[0]) {
					$baf = 0;
				} elsif ($theta < $thetamean[1]) {
					$baf = 0.5 * ($theta-$thetamean[0]) / ($thetamean[1]-$thetamean[0]);
				} elsif ($theta < $thetamean[2]) {
					$baf = 0.5 + 0.5 * ($theta-$thetamean[1]) / ($thetamean[2]-$thetamean[1]);
				} else {
					$baf = 1;
				}
				
				if ($theta < $thetamean[1]) {
					my $rexpected;
					if ($thetamean[1]-$thetamean[0]) {		#sometimes the two reference theta are identical (do not make sense but it may happen)
						$rexpected = $rmean[0] + ($theta-$thetamean[0]) * ($rmean[1]-$rmean[0]) / ($thetamean[1]-$thetamean[0]);
					} else {
						$rexpected = $rmean[0];
					}
					if ($rexpected < 0) {
						print STDERR "WARNING: R_expected < 0 because rmean=@rmean thetamean=@thetamean a=$signala b=$signalb r=$r theta=$theta\n";
						$rexpected = $r;		#so that LRR=0
					}
					$lrr = log ($r/$rexpected) / log(2);
				} else {
					my $rexpected;
					if ($thetamean[2]-$thetamean[1]) {
						$rexpected = $rmean[1] + ($theta-$thetamean[1]) * ($rmean[2]-$rmean[1]) / ($thetamean[2]-$thetamean[1]);
					} else {
						$rexpected = $rmean[1];
					}
					if ($rexpected < 0) {
						print STDERR "WARNING: R_expected < 0 due to rmean=@rmean thetamean=@thetamean a=$signala b=$signalb r=$r theta=$theta\n";
						$rexpected = $r;		#so that LRR=0
					}
					$lrr = log ($r/$rexpected) / log(2);
				}
				
				$lrr = sprintf ("%.4f", $lrr);
				$baf = sprintf ("%.4f", $baf);
				print OUTPUT "\t$lrr\t$baf";
			}
		} else {		#Non-polymorphic probe
			if (not $power2) {					#convert the signal to log2 based before subtracting the median values
				@siga = map {log($_)/log(2)} @siga;
			}
			my ($lrr, $baf);
			for my $i (0 .. @header-1) {
				my @npinfo = split (/\t/, $cluster->{$sig_psida});
				$lrr = sprintf ("%.4f", $siga[$i] - $npinfo[0]);
				$baf = 2;						#BAF=2 indicates a non-polymorphic marker!!!!!!
				print OUTPUT "\t$lrr\t$baf";
			}
		}
		print OUTPUT "\n";
		$count_success++;
	}
	close (OUTPUT);
	print STDERR "NOTICE: Finished writting to $outputfile with $count_success LRR and BAF values\n";
	$count_nocluster and print STDERR "WARNING: A total of $count_nocluster markers do not have cluster information and were skipped\n";
	$count_noloc and print STDERR "WARNING: A total of $count_noloc markers do not have chromosome location annotation and were skipped\n";
	$count_badchr and print STDERR "WARNING: A total of $count_badchr markers in unrecognizable chromosome (such as 0, ---, MT) were skipped\n";

}

sub readLocFile {
	my ($inputfile) = @_;
	my (%snploc);
	print STDERR "NOTICE: Reading snp-location-file $inputfile ...\n";
	open (LOC, $inputfile) or confess "Error: cannot read from inputfile $inputfile";
	while (<LOC>) {
		s/[\r\n]+$//;
		m/^([^\t]+)\t([^\t]+\t[^\t]+)/ or print STDERR "WARNING: unrecognizable record skipped in locfile $inputfile: <$_> (at least 3 tab-delimited fields expected)\n" and next;
		$snploc{$1} = $2;
	}
	close (LOC);
	print STDERR "NOTICE: Reading snp-locatoin-file $inputfile done with ${\(scalar keys %snploc)} markers!\n";
	return (\%snploc);
}

sub readClusterFile {
	my ($inputfile) = @_;
	my (%cluster, $count);
	print STDERR "NOTICE: Beginning reading canonical cluster information from $inputfile ...";
	open (CLUSTER, $inputfile) or confess "Error: cannot read from clusterfile $inputfile: $!";
	while (<CLUSTER>) {
		s/[\r\n]+$//;
		m/^(\S+)\t(.+)/ or confess "Error: invalid record found in clusterfile: <$_>";
		$cluster{$1} = $2;
		$count++;
	}
	close (CLUSTER);
	print STDERR "Done with $count markers!\n";
	return (\%cluster);
}


=head1 SYNOPSIS

 normalize_affy_geno_cluster.pl [arguments] <genotype-cluster-file> <signal-file>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --locfile <file>		a file containing genomic location of markers
 	    --output <file>		output file name (default: <signal-file>.lrr_baf.txt)
 	    --(no)power2		signal file contains log2 signal intensity (default:ON)

 Function: given a canonical genotype cluster file, convert raw signal intensity 
 for Affymetrix SNP arrays to Log R Ratio and B Allele Frequency values 

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--locfile>

specify a file that contains the genomic location annotation for markers. The 
file contains information for one marker per line, and each line contains three 
tab-delimited records, representing marker name, chromosome and position 
information.

=item B<--output>

specify the output file name.

=item B<--(no)power2>

specifies whether the signal file contains log2 signal intensity or not. By 
default it is on, and the values in the file are usually around 10. When the 
signal file has values around 1000-10000, then it means the the signal file 
contains raw signal values, and one should turn off this argument.


=back

=head1 DESCRIPTION

This program is used to convert raw signal intensity values from Affymetrix SNP 
arrays to Log R Ratio (LRR) and B Allele Frequency (BAF) measures. The signal 
intensity should be calculated by the Affymetrix Power Tools, and an appropriate 
canonical clustering file is necessary for the conversion. This clustering file 
can be generated from the generate_affy_geno_cluster.pl program.

For questions or comments on this program, please contact 
kai@mail.med.upenn.edu.