#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: 270 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2009-06-12 13:29:23 -0400 (Fri, 12 Jun 2009) $';

our ($verbose, $help, $man);
our ($genofile, $conffile, $sigfile);
our ($confidence_threshold, $clean, $sexfile, $locfile, $output, $median_cluster, $min_subject, $power2, $ignore_name_discord);

GetOptions('verbose'=>\$verbose, 'help'=>\$help, 'man|m'=>\$man, 'confidence=f'=>\$confidence_threshold, 'clean'=>\$clean, 'sexfile=s'=>\$sexfile,
	'locfile=s'=>\$locfile, 'output=s'=>\$output, 'median!'=>\$median_cluster, 'min_subject=i'=>\$min_subject, 'power2!'=>\$power2,
	'ignore_name_discord'=>\$ignore_name_discord) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 3 or pod2usage ("Syntax error: three <inputfiles> are required");

($genofile, $conffile, $sigfile) = @ARGV;
if (not defined $confidence_threshold) {
	print STDERR "NOTICE: The --confidence_threshold argument is automatically set as 0.01 (default for Affymetrix Power Tools calling)\n";
	$confidence_threshold = 0.01;
}
defined $median_cluster or $median_cluster = 1;
$min_subject ||= 10;
defined $power2 or $power2 = 1;			#by default turn on the power2 option

$locfile or pod2usage ("Error in argument: please specify the --locfile argument which provides the chromosome location of markers");
if ($output) {
	open (STDOUT, ">$output") or confess "Error: cannot write to output file $output: $!";
}

my $file_sex = {};
$sexfile and $file_sex = readSexFile ($sexfile);
my $snpchr = readLocFile ($locfile);
generateAffyGenoCluster ($genofile, $conffile, $sigfile, $snpchr, $file_sex);
generateAffyCNCluster ($sigfile, $snpchr, $file_sex);

sub generateAffyGenoCluster {
	my ($genofile, $conffile, $sigfile, $snpchr, $file_sex) = @_;
	open (GENO, $genofile) or confess "Error: cannot read from genotype file $genofile: $!\n";
	open (CONF, $conffile) or confess "Error: cannot read from confidence file $conffile: $!\n";
	open (SIG, $sigfile) or confess "Error: cannot read from signal intensity file $sigfile: $!\n";
	
	while (<GENO>) {
		m/^#/ or last;
	}
	m/^probeset_id\t/ or confess "Error: invalid record found in genotype file $genofile (probeset_id expected): <$_>";
	s/\S+[\\\/]//g;						#get rid of directory name in file path and only keep file name
	s/\s*[\r\n]+$//;
	my $geno_header_line = $_;
	$verbose and print STDERR "NOTICE: the genoytpe file $genofile is open and ready for processing\n";
	
	while (<CONF>) {
		m/^#/ or last;
	}
	m/^probeset_id\t/ or confess "Error: invalid record found in genotype file $conffile: <$_>";
	s/\S+[\\\/]//g;						#get rid of directory name in file path and only keep file name
	s/\s*[\r\n]+$//;
	$_ eq $geno_header_line or confess "Error: mismatch between the header line in genotype file (\n$geno_header_line) and in confidence file (\n$_)";
	$verbose and print STDERR "NOTICE: the conf file $conffile is open and ready for processing\n";
	
	while (<SIG>) {
		m/^#/ or last;
	}
	m/^probeset_id\t/ or confess "Error: invalid record found in genotype file $conffile (probeset_id expected at the beginning of the line): <$_>";
	s/\S+[\\\/]//g;						#get rid of directory name in file path and only keep file name
	s/\s*[\r\n]+$//;
	if ($_ ne $geno_header_line) {
		my @temp1 = split (/\t/, $_);
		my @temp2 = split (/\t/, $geno_header_line);
		for my $i (0 .. @temp1-1) {
			if ($temp1[$i] ne $temp2[$i]) {
				print STDERR "WARNING: Discordance in file name records in signalfile and genofile: index=$i signalfile=$temp1[$i] genofile=$temp2[$i]\n";
				if ($ignore_name_discord) {
					print STDERR "NOTICE: Since you have specified --ignore_name_discord, these discordances will be ignored\n" and last;
				} else {
					print STDERR "ERROR: Program stopped due to the mismatch between the header line in genotype file $genofile and in signal file $sigfile (use --ignore_name_discord to override)\n";
					exit (1);
				}
			}
		}
	}
		
	my ($processxy, %index_male, %index_female) = (0);
	my @header = split (/\t/, $_);
	shift @header;						#first element is "probeset_id"
	for my $i (0 .. @header-1) {
		exists $file_sex->{$header[$i]} or next;
		if ($file_sex->{$header[$i]} == 1) {
			$index_male{$i}++;
		} elsif ($file_sex->{$header[$i]} == 2) {
			$index_female{$i}++;
		}
	}
	if (scalar (keys %index_male) < $min_subject or scalar (keys %index_female) < $min_subject) {
		print STDERR "WARNING: unable to find enough males (count=${\(scalar keys %index_male)}) and females (count=${\(scalar keys %index_female)}) for chrX/chrY cluster generation (use --sexfile to supply this information)\n";
	} else {
		$processxy++;
		print STDERR "NOTICE: A total of  ${\(scalar keys %index_male)} males and ${\(scalar keys %index_female)} females are found in the signal file $sigfile\n";
	}

	print "probeset_id\tr_aa\tr_ab\tr_bb\ttheta_aa\ttheta_ab\ttheta_bb\tcount_aa\tcount_ab\tcount_bb";
	$clean and print "\ts_aa\ts_ab\ts_ba\ts_bb";			#the --clean is an OBSELETE argument and should NOT be used any more
	print "\n";

	my ($geno_psid, $conf_psid, $sig_psid, @geno, @conf, @sig, @siga, @sigb);
	my ($saa, $sab, $sba, $sbb);
	my ($count_genoline, $count_sigline, $count_success, $count_noanno, $count_genofail, $count_nopoly, $count_abnormal) = qw/0 0 0 0 0 0 0/;
	while (<GENO>) {
		s/\s*[\r\n]+$//;
		@geno = split (/\t/, $_);
		$geno_psid = shift @geno;
		$verbose and defined $snpchr->{$geno_psid} and print STDERR "NOTICE: Processing marker $geno_psid in chr$snpchr->{$geno_psid}\n";
		$count_genoline++;
		
		$_ = <CONF>;
		s/\s*[\r\n]+$//;
		@conf = split (/\t/, $_);
		$conf_psid = shift @conf;
		$conf_psid eq $geno_psid or confess "Error: probeset identifier mismatch between genotype file $genofile and confidence file $conffile: $geno_psid vs $conf_psid\n";
		@geno == @conf or confess "Error: probeset $conf_psid record number mismatch between genotype file $genofile and confidence file $conffile: ${\(scalar @geno)} versus ${\(scalar @conf)}";
		
		while (<SIG>) {
			$count_sigline++;
			s/\s*[\r\n]+$//;
			@siga = split (/\t/, $_);
			$sig_psid = shift @siga;
			$sig_psid eq $geno_psid.'-A' and last;
		}
		$_ = <SIG>;
		$count_sigline++;
		s/\s*[\r\n]+$//;
		@sigb = split (/\t/, $_);
		$sig_psid = shift @sigb;
		$sig_psid eq $geno_psid.'-B' or confess "Error: Found probeset $sig_psid A allele but not B allele in signal file $sigfile line $count_sigline (genofile line $count_genoline) for marker $geno_psid: <$_>";
		
		@geno == @siga or confess "Error: probeset $sig_psid record number mismatch between genotype file $genofile and signal file $sigfile: ${\(scalar @geno)} versus ${\(scalar @siga)}";
		@geno == @sigb or confess "Error: probeset $sig_psid record number mismatch between genotype file $genofile and signal file $sigfile: ${\(scalar @geno)} versus ${\(scalar @sigb)}";

		defined $snpchr->{$geno_psid} or ++$count_noanno and next;
		if ($snpchr->{$geno_psid} eq 'X' or $snpchr->{$geno_psid} eq 'Y') {				#skip X and Y markers
			$processxy or next;
		}
		
		if ($power2) {
			#transform signalA and signalB to real measurements !!!!!!!!!!!!!
			@siga = map {2**$_} @siga;
			@sigb = map {2**$_} @sigb;
		}
		
		if ($clean) {			#OBSELETE procedure
			($saa, $sab, $sba, $sbb) = cleanSignal (\@geno, \@conf, \@siga, \@sigb);		#if --clean argument is set, use a model that substract the cross-hybridization signal (eliminate background noise) to re-caluculate @siga and @sigb
		}
		
		my (@raa, @rab, @rbb, @thetaaa, @thetaab, @thetabb);
		my @valid_count = qw/0 0 0/;
		for my $i (0 .. @geno-1) {
			$conf[$i] > $confidence_threshold and next;
			if (defined $snpchr->{$geno_psid} and $snpchr->{$geno_psid} eq 'X') {			#chrX only
				$index_female{$i} or next;
			}
			if (defined $snpchr->{$geno_psid} and $snpchr->{$geno_psid} eq 'Y') {			#chrY only
				$index_male{$i} or next;
			}
			
			if ($geno[$i] eq '0') {				#AA genotype
				push @raa, $siga[$i]+$sigb[$i];
				push @thetaaa, atan2 ($sigb[$i], $siga[$i]) / (3.1415926/2);
				$valid_count[0]++;
			} elsif ($geno[$i] eq '1') {			#AB genotype
				push @rab, $siga[$i]+$sigb[$i];
				push @thetaab, atan2 ($sigb[$i], $siga[$i]) / (3.1415926/2);
				$valid_count[1]++;
			} elsif ($geno[$i] eq '2') {			#BB genotype
				push @rbb, $siga[$i]+$sigb[$i];
				push @thetabb, atan2 ($sigb[$i], $siga[$i]) / (3.1415926/2);
				$valid_count[2]++;
			} else {
				print STDERR "WARNING: No Call found for $geno_psid: <$geno[$i]> (this usually should NOT happen!)\n";
				next;
			}
		}
		
		my (@rmean, @thetamean);
		if ($median_cluster) {
			@rmean = (median (\@raa), median (\@rab), median (\@rbb));
			@thetamean = (median (\@thetaaa), median (\@thetaab), median (\@thetabb));
		} else {
			@rmean = (mean (\@raa), mean (\@rab), mean (\@rbb));
			@thetamean = (mean (\@thetaaa), mean (\@thetaab), mean (\@thetabb));
		}
		
		if (not defined $rmean[0] and not defined $rmean[1] and not defined $rmean[2]) {
			$count_genofail++;
			next;
		}
		
		if (not $rmean[0] and not $rmean[1] or not $rmean[0] and not $rmean[2] or not $rmean[1] and not $rmean[2]) {
			$count_nopoly++;
			next;
		}
			
		$rmean[0] ||= ($rmean[1] || $rmean[2]);
		$rmean[1] ||= ($rmean[0] || $rmean[2]);
		$rmean[2] ||= ($rmean[1] || $rmean[0]);
		if (not $thetamean[1]) {
			if ($thetamean[0] and $thetamean[2]) {
				$thetamean[1] = ($thetamean[0]+$thetamean[2])/2;	#make sure that AB is between AA and BB whenever possible
			} else {
				$thetamean[1] = 0.5;
			}
		}
		$thetamean[0] ||= ($thetamean[1]-0.3);			#empirically, I have found that the theta values typically are centered around 0.2, 0.5 and 0.8, respectively
		$thetamean[2] ||= ($thetamean[1]+0.3);
		if ($thetamean[0] > $thetamean[1] or $thetamean[2] < $thetamean[1]) {
			$count_abnormal++;
			next;
		}
		
		@rmean = map {sprintf ("%.4f", $_)} @rmean;
		@thetamean = map {sprintf ("%.4f", $_)} @thetamean;
		print join ("\t", $geno_psid, @rmean, @thetamean, @valid_count);
		if ($clean) {
			($saa, $sab, $sba, $sbb) = map {sprintf ("%.4f", $_)} ($saa, $sab, $sba, $sbb);
			print "\t", join ("\t", $saa, $sab, $sba, $sbb);
		}
		print "\n";
		$count_success++;
		$count_genoline =~ m/000$/ and print STDERR "NOTICE: Finished processing $count_genoline markers (cluster generation for $count_success markers have been successfully)\n";
		
	}
	close (GENO); close (CONF); close (SIG);

	$count_noanno and print STDERR "WARNING: A total of $count_noanno markers do not have chromosome annotation in locfile $locfile and were skipped\n";
	$count_genofail and print STDERR "WARNING: A total of $count_genofail markers have complete genotyping failure (no genotype calls pass the confidence score threshold) and were skipped\n";
	$count_nopoly and print STDERR "WARNING: A total of $count_nopoly markers do not have at least two types of genotypes (out of AA, AB, BB) and were skipped\n";
	$count_abnormal and print STDERR "WARNING: A total of $count_abnormal markers have abnormal thetamean patterns (mean theta value for AA, AB and BB cluster) and were skipped\n";
	print STDERR "NOTICE: A total of $count_success SNP markers have been analyzed to construct canonical clustering positions\n";
}

sub cleanSignal {
	my ($geno, $conf, $siga, $sigb) = @_;
	my ($e1, $e2, $e3, $e4, $e5, $e6) = qw/0 0 0 0 0 0/;
	my ($naa, $nab, $nbb) = qw/0 0 0/;
	my ($saa, $sba, $sab, $sbb);
	my (@cleansiga, @cleansigb);
	for my $i (0 .. @$geno-1) {
		$conf->[$i] > $confidence_threshold and next;
		if ($geno->[$i] == 0) {
			$e1 += $siga->[$i];
			$e2 += $sigb->[$i];
			$naa++;
		} elsif ($geno->[$i] == 1) {
			$e3 += $siga->[$i];
			$e4 += $sigb->[$i];
			$nab++;
		} elsif ($geno->[$i] == 2) {
			$e5 += $siga->[$i];
			$e6 += $sigb->[$i];
			$nbb++;
		}
	}
	$saa = (2*$nbb*($e1+$e3) + $nab*($e1-$e5)) / ((2*$naa+$nab)*(2*$nbb+$nab) - $nab*$nab);
	$sab = (2*$nbb*($e2+$e4) + $nab*($e2-$e6)) / ((2*$naa+$nab)*(2*$nbb+$nab) - $nab*$nab);
	$sba = (2*$naa*($e3+$e5) + $nab*($e5-$e1)) / ((2*$naa+$nab)*(2*$nbb+$nab) - $nab*$nab);
	$sbb = (2*$naa*($e4+$e6) + $nab*($e6-$e2)) / ((2*$naa+$nab)*(2*$nbb+$nab) - $nab*$nab);
	
	for my $i (0 .. @$geno-1) {
		my $ca = ($siga->[$i]*$sbb - $sigb->[$i]*$sba) / ($saa*$sbb - $sab*$sba);
		my $cb = ($siga->[$i]*$sab - $sigb->[$i]*$saa) / ($sba*$sab - $saa*$sbb);
		$ca < 0.001 and print "i=$i $ca less than 0.001\n" and $ca = 0.001 ;
		$cb < 0.001 and print "i=$i $cb less than 0.001\n" and $cb = 0.001 ;
		$cleansiga[$i] = $ca * $saa;
		$cleansigb[$i] = $cb * $sbb;
		#print STDERR "NOTICE: i=$i e1-6=$e1 $e2 $e3 $e4 $e5 $e6\n";
		#print STDERR "NOTICE: i=$i a=$siga->[$i] b=$sigb->[$i] saa=$saa sab=$sab sba=$sba sbb=$sbb cleana=$cleansiga[$i] cleanb=$cleansigb[$i]\n";
	}
	@$siga = @cleansiga;
	@$sigb = @cleansigb;
	return ($saa, $sab, $sba, $sbb);
}

sub readSexFile {
	my ($sexfile) = @_;
	my (%sex);
	open (SEX, $sexfile) or confess "Error: cannot read from sex file $sexfile: $!";
	while (<SEX>) {
		s/[\r\n]+$//;
		m/^([^\t]+)\t(\S+)\s*$/ or confess "Error: invalid record found in sexfile (two tab-delimited records expected): <$_>";
		my ($ind, $gen) = ($1, $2);
		$gen eq 'male' and $gen = 1;
		$gen eq 'female' and $gen = 2;
		$gen =~ m/^[012]$/ or confess "Error: Erraneous sex annotation ($gen): only 0, 1, 2, male or female are correct annotations in sexfile $sexfile\n";
		$sex{$ind} = $gen;
	}
	close (SEX);
	return (\%sex);
}

sub readLocFile {
	my ($locfile) = @_;
	my (%snpchr);
	print STDERR "NOTICE: Reading marker-location-file $locfile ...";
	open (LOC, $locfile) or confess "Error: cannot read from inputfile $locfile";
	while (<LOC>) {
		s/[\r\n]+$//;
		m/^([^\t]+)\t([^\t]+)/ or print STDERR "WARNING: unrecognizable record skipped in locfile $locfile: <$_> (at least 2 tab-delimited fields expected)\n" and next;
		$snpchr{$1} = $2;
	}
	close (LOC);
	print STDERR " Done with ${\(scalar keys %snpchr)} markers!\n";
	return (\%snpchr);
}

sub generateAffyCNCluster {
	my ($sigfile, $snpchr, $file_sex) = @_;

	my (@header, @index_male, @index_female);
	my $processxy = 0;			#by default do not process chrX and chrY
	open (SIG, $sigfile) or confess "Error: cannot read from signal file $sigfile: $!";
	while (<SIG>) {
		m/^probeset_id/ and last;
	}
	s/\s*[\r\n]+$//;
	
	@header = split (/\t/, $_);
	shift @header;
	for my $i (0 .. @header-1) {
		exists $file_sex->{$header[$i]} or next;
		if ($file_sex->{$header[$i]} == 1) {
			push @index_male, $i;
		} elsif ($file_sex->{$header[$i]} == 2) {
			push @index_female, $i;
		}
	}
	if (@index_male >= $min_subject and @index_female >= $min_subject) {
		$processxy++;
		#print STDERR "NOTICE: A total of ${\(scalar @index_male)} males and ${\(scalar @index_female)} females will be used in chrX and chrY cluster construction\n";
	}
	
	my ($count_success, $count_noanno) = qw/0 0/;
	my ($count_goodx, $count_badx, $count_goody, $count_bady) = qw/0 0 0 0/;
	my (@meandifx, @meandify);
	while (<SIG>) {
		m/^CN/ or next;
		s/\s*[\r\n]+$//;
		my @sig = split (/\t/, $_);
		my $psid = shift @sig;
		my $curchr = $snpchr->{$psid};
		defined $curchr or ++$count_noanno and next;	#this marker is not annotated in locfile so skipped
		
		if (not $power2) {				#by default, log2(signal) expected, unless --nopower2 is specified.
			@sig = map {log($_)/log(2)} @sig;
		}
		
		if ($curchr eq 'X' or $curchr eq 'Y') {
			$processxy or next;			#skip chrX and chrY
			my ($mean1, $mean2, $sd1, $sd2) = (mean ([@sig[@index_male]]), mean ([@sig[@index_female]]), sd ([@sig[@index_male]]), sd ([@sig[@index_female]]));
			if ($curchr eq 'X') {
				if ($mean2-$mean1 > $sd1+$sd2) {	#this is ad hoc formula to evaluate how "close" the signal intensity means are for male and female
					$count_goodx++;
				} else {
					$count_badx++;
				}
				push @meandifx, $mean2-$mean1;
				print $psid, "\t", sprintf ("%.4f", median ([@sig[@index_female]])), "\t", sprintf ("%.4f", $mean2), "\t", sprintf ("%.4f", $sd2), "\n";	#female only
			} elsif ($curchr eq 'Y') {
				if ($mean1-$mean2 > $sd1+$sd2) {
					$count_goody++;
				} else {
					$count_bady++;
				}
				push @meandify, $mean1-$mean2;
				print $psid, "\t", sprintf ("%.4f", median ([@sig[@index_male]])), "\t", sprintf ("%.4f", $mean1), "\t", sprintf ("%.4f", $sd1), "\n";		#male only
			}
		} else {
			print $psid, "\t", sprintf ("%.4f", median (\@sig)), "\t", sprintf ("%.4f", mean (\@sig)), "\t", sprintf ("%.4f", sd (\@sig)), "\n";
		}
		$count_success++;		
	}
	close (SIG);
	$count_noanno and print STDERR "WARNING: A total of $count_noanno CN markers do not have chromosome annotation and were skipped\n";
	print STDERR "NOTICE: A total of $count_success CN markers have been analyzed to construct canonical clustering positions\n";
	
	if ($processxy) {
		@meandifx and print STDERR "INFORMATION: For CN probes, meandif_X=", mean (\@meandifx), " and sd=", sd (\@meandifx), " count=", scalar (@meandifx), " (including $count_goodx high-quality markers)\n";
		@meandify and print STDERR "INFORMATION: For CN probes, meandif_Y=", mean (\@meandify), " and sd=", sd (\@meandify), " count=", scalar (@meandify), " (including $count_goody high-quality markers)\n";
	}
}


sub sum {
	my ($score) = @_;
	@$score or return undef;
	my $sum;
	for (@$score) {
		$sum += $_;
	}
	return $sum;
}

sub mean {
	my ($score) = @_;
	@$score or return undef;
	my $sum;
	for (@$score) {
		$sum += $_;
	}
	return $sum/@$score;
}

sub sd {
	my ($score) = @_;
	@$score or return undef;
	my $mean = mean ($score);
	my $sum;
	for my $i (0 .. @$score-1) {
		$sum += ($score->[$i]-$mean)*($score->[$i]-$mean);
	}
	$sum /= (@$score-1);
	return sqrt ($sum);
}

sub median {
	my ($score) = @_;
	@$score or return undef;
	my @newscore = sort {$a<=>$b} @$score;
	if (@newscore % 2 == 0) {
		return ($newscore[@newscore/2-1]+$newscore[@newscore/2])/2;
	} else {
		return $newscore[@newscore/2];
	}
}

=head1 SYNOPSIS

 generate_affy_geno_cluster.pl [arguments] <genotype-file> <confidence-file> <signal-file>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --sexfile <file>		a file containing sex of samples
 	    --locfile <file>		a file containing genomic location of markers
 	    --output <file>		specify output file name (default: STDOUT)
 	    --confidence <float>	threshold for genotyping calls confidence score (default:0.01)
 	    --min_subject <int>		minimum number of subject for each gender (default: 10)
 	    --(no)median		calculate median cluster position (default:ON)
 	    --clean			clean cross-hybridization signals (OBSELETE option)
 	    --(no)power2		use power(2,x) transformation for input signal (default:ON)
 	    --ignore_name_discord	ignore cases when file names are discordant in the three input files

 Function: based on Affymetrix Power Tools genotype calls and signal intensity 
 values, compute the R and theta values for three canonical genotype clusters 
 for each SNP or CN marker.

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--sexfile>

specify a file that contains the sex annotation of samples. The file contains 
information for one sample per line, and each line contains two tab-delimited 
records, representing file name and gender (1=male, 2=female, 0=unknown), 
respectively. (Alternatively, one can use male, female, unknown as the second 
field in the file explicitely, without using 1, 2 and 0 symbol.)

=item B<--locfile>

specify a file that contains the genomic location annotation for markers. The 
file contains information for one marker per line, and each line contains 
tab-delimited records, the first two of which represent marker name and 
chromosome information.

=item B<--output>

specify the output file names. By default the output will be printed to STDOUT.

=item B<--confidence>

confidence score threshold for genotype calling. The default score used in this 
program is 0.01, which is a highly stringent criteria and will usually results 
in high non-call rate. High stringency ensures high quality of genotype cluster 
positions.

=item B<--(no)median>

use median value rather than mean value when specifying the cluster positions.

=item B<--min_subject>

specify the minimum number of subjects for each gender for chrX and chrY cluster 
construction. If the number of subjects for each gender is less than this 
number, then chrX and chrY markers will not be processed for cluster generation.

=item B<--clean>

OBSELETE argument! use model-based method to adjust background noise and 
generate cleaned version of cross-hybridization signals.

=back

=head1 DESCRIPTION

This program is used for generating reference intensity clusters for Affymetrix 
SNP arrays. More specifically, given a large data set with genotype, confidence 
(of genotype calling) and signal intensity values for many individuals 
calculated by the Affymetrix Power Tools, this program computes the R and theta 
values for three canonical genotype clusters (AA, AB and BB genotype) for each 
SNP, and it also computes the reference median and mean values of each 
non-polymorphic probe.

The generated genotype cluster file will be useful for normalizing and 
converting signal intensity of new samples that were genotyped in similar 
fashion. This process can be performed by the B<normalize_affy_geno_cluster.pl> 
program. In the end the Log R Ratio (LRR) value and the B Allele Frequency (BAF) 
value of each SNP or NP marker will be computed for the new study samples. (For 
NP markers, the BAF value is set arbitrarily as 2). After the LRR and BAF values 
are generated, one can use various CNV detection program (such as PennCNV) for 
CNV inference and analysis.

This program takes three main input files, a genotype call file, a confidence 
file and a signal intensity file. They are briefly described below:

=over 8

=item B<Genotype call file>

the genotype call file can be generated by the apt-probeset-genotype program in 
the APT (Affymetrix Power Tools) package, which is composed of a few command 
line driven programs that automate several signal processing procedures for 
Affymetrix CEL files. See the documentation for APT for more detailed usage 
information. Note that all markers in all individuals will have genotype calls, 
but with different confidence scores; as we are only interested in using high-
quality genotypes for constructing clustering positions, we will need to use the 
confidence score file.

An example of command line usage is shown below:

 [kai@node126 apt]$ apt-probeset-genotype -c GenomeWideSNP_6.cdf -a birdseed --read-models-birdseed GenomeWideSNP_6.birdseed.models --special-snps GenomeWideSNP_6.specialSNPs --out-dir . *.CEL

This command request to run the birdseed genotyping calling algorithm on all CEL 
files in the current directory. Two output files will be generated, one 
containing genotyping calls and the other containing confidence scores.

=item B<Genotype confidence score file>

this file is generated at the same time as the genotype call file by the apt-
probeset-genotype program. It contains the confidence score for all genotype 
calls. By default only highly confident (score < 0.01) genotyping calls are used 
in the construction of clustering positions.

=item B<signal intensity file>

the signal intensity file can be generated by the apt-probeset-summarize program 
in the APT package. An example of command line usage is shown below:

 [kai@cytosine apt]$ apt-probeset-summarize --cdf-file GenomeWideSNP_6.cdf --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --out-dir . *.CEL

This command request processing all CEL files in the current directory with 
quantile normalization, followed by median polish, and then export the signal 
intensity values for A and B alleles. The pm-only option specifies that only PM 
probes be used in the data processing, but this does not matter for genome-wide 
6.0 arrays, which only contains PM probes.

=back

This program also needs a location file that contains the genomic location 
annotation for markers. The file contains information for one marker per line, 
and each line contains at least two tab-delimited records, representing marker 
name and chromosome information. Indeed, the only purpose of this file is to 
specify chrX and chrY specific markers for accurate processing of signal 
intensities and clustering of genotyping calls.

This program also needs a file that contains the sex annotation of samples. The 
file contains information for one sample per line, and each line contains two 
tab-delimited records, representing file name and gender (1=male, 2=female, 
0=unknown), respectively. It is important to have the correct sex annotation: 
for chrX, only females are used in the construction of genotyping clusters; for 
chrY, only males are used in the construction of genotyping clusters. Without 
specifying --sexfile, the chrX and chrY markers will not be processed!

The final output file is explained here: The first line of the output file is 
called header line, which contains column annotations for subsequent lines. The 
following lines contains the clustering positions for SNP markers and the signal 
intensity values for non-polymorphic markers.

An example output file is shown below:

	probeset_id     r_aa    r_ab    r_bb    theta_aa        theta_ab        theta_bb        count_aa        count_ab        count_bb
	SNP_A-2131660   2503.3793       2697.9633       2551.7423       0.1652  0.4829  0.8110  40      396     879
	SNP_A-1967418   819.9688        770.7794        734.2616        0.2188  0.4878  0.7735  30      342     934
	SNP_A-1969580   3297.4381       3297.4381       3524.6461       0.2644  0.5644  0.8437  0       14      1365
	SNP_A-4263484   2053.8534       1940.7644       1605.4141       0.1020  0.4008  0.8280  128     585     672
	SNP_A-1978185   1361.8169       2138.6936       3027.7066       0.1573  0.6466  0.8333  1342    32      3

The above shows the R values, the theta values and the counts for three 
genotypes for each SNP marker.

For NP probes, the format is a little bit different:

	CN_441519       10.51856        10.5260246530872        0.274807919008648
	CN_441512       10.82210        10.8057111839592        0.24389875689443
	CN_441516       10.26397        10.2513341947804        0.400280575164214

There are three columns after the marker name, representing median value, mean 
value and standard deviation of NP probes for presumably 2-copy marker (1-copy 
for chrY) in base 2 logarithm scale.

Note that this program first process all SNP markers and then process all NP 
markers (in fact the same file is read twice), so all the NP markers will be at 
the later half of the output file.

For questions or comments on this program, please contact 
kai@openbioinformatics.org.