#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2010-07-08 04:34:07 -0700 (Thu, 08 Jul 2010) $';

our ($verbose, $help, $man);
our ($adjust, $calwf, $refchr, $gcmodelfile, $suffix, $distance, $gap, $skip_marker, $listfile, $prefix, $output);
our (@inputfile);

GetOptions ('verbose'=>\$verbose, 'help'=>\$help, 'man|m'=>\$man, 'adjust'=>\$adjust, 'calwf'=>\$calwf, 'refchr=s'=>\$refchr, 'gcmodelfile=s'=>\$gcmodelfile,
	'suffix=s'=>\$suffix, 'distance=s'=>\$distance, 'gap=i'=>\$gap, 'skip_marker=i'=>\$skip_marker, 'listfile=s'=>\$listfile, 'prefix=s'=>\$prefix,
	'output=s'=>\$output) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);

(@inputfile) = @ARGV;


$refchr ||= '11';
$refchr eq '11' or pod2usage ("Error in argument: the --ref_chromosome must be 11 in the current version of the program for processing human genome");
$suffix ||= 'adjusted';
$distance ||= 1_000_000;			#regress on every 1Mb by default
$distance =~ s/m$/000000/;
$distance =~ s/k$/000/;
$distance =~ m/^\d+$/ or pod2usage ("Error in argument: --distance should be an integer (suffix of k or m is okay)");
$skip_marker ||= 0;

if (@ARGV) {
	$listfile and pod2usage ("Error in argument: please do not specify the --listfile argument when providing signal file names (@ARGV) in command line\n");
	@inputfile = @ARGV;
} elsif ($listfile) {
	open (LIST, $listfile) or confess "Error: cannot read from listfile $listfile: $!";
	while (<LIST>) {
		s/[\r\n]+$//;
		m/^\s*(.+?)\s*$/ or confess "Error: the --listfile should contain one file name per line: <$_>\n";
		push @inputfile, $1;
	}
	close (LIST);
} else {
	pod2usage ("Error in argument: please either specify the --listfile argument or provide signal file names in command line\n");
}

if ($adjust) {
	defined $gcmodelfile or pod2usage ("Error in argument: please specifiy the --snpgcfile argument for --adjust operation");
	my $snpgc = readSNPGC ($gcmodelfile);
	
	for my $i (0 .. @inputfile-1) {
		my ($b1, $b2) = trainGCModel ($inputfile[$i], $snpgc, $gap, $distance, $skip_marker);
		my $outputfile;
		if ($prefix) {
			-d $prefix or die "Error: the output directory $prefix (specified by --prefix) does not exist\n";
			my ($vol, $dir, $file) = File::Spec->splitpath ($inputfile[$i]);
			$outputfile = File::Spec->catfile ($prefix, "$file.$suffix");
		} else {
			$outputfile = "$inputfile[$i].$suffix";
		}
		print STDERR "NOTICE: Output file with adjusted signal intensity is written to $outputfile\n";
		adjustSignal ($inputfile[$i], $snpgc, $outputfile, $b1, $b2);
	}
} elsif ($calwf) {
	if (defined $output) {
		open (OUTPUT, ">$output");
		print STDERR "NOTICE: the waviness factor (WF) values will be written to $output\n";
	}
	for my $nextfile (@inputfile) {
		my ($wf, $cc) = calculateMADWF ($nextfile, $refchr);
		my $gcwf = $wf * abs ($cc);
		defined $output and print OUTPUT "$nextfile\t$gcwf\t$wf\t$cc\n";
		print STDERR "NOTICE: For $nextfile, GC-wave factor (GCWF) is $gcwf, WF is $wf, GC correlation is $cc\n";
	}
	defined $output and close (OUTPUT);
} else {
	pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
}
		
sub readSNPGC {
	my ($snpgcfile) = @_;
	my (%snpgc);
	open (SNPGC, $snpgcfile) or confess "Error: cannot read from SNP-GC file $snpgcfile: $!";
	while (<SNPGC>) {
		s/[\r\n]+$//;
		my @record = split (/\s+/, $_);
		$snpgc{$record[0]} = [@record[1..3]];
	}
	close (SNPGC);
	return \%snpgc;
}


sub adjustSignal {
	my ($signalfile, $snpgc, $outputfile, $b1, $b2) = @_;
	
	my ($lrr_index, $name_index, $header, $header_seg);
	my $count_nowave = 0;
	
	open (OUT, ">$outputfile") or confess "Error: cannot write to output file $outputfile: $!";
	open (SIG, $signalfile) or confess "Error: cannot read from file $signalfile: $!";
	$header = <SIG>;

	$header =~ m/(.+)Log R Ratio/ or confess "Error: the header file of signalfile $signalfile does not contain 'Log R Ratio' annotation";
	$header_seg = $1;
	$lrr_index = ($header_seg =~ tr/\t/\t/+0);

	$header =~ m/(.*)Name/ or confess "error: the header file of signalfile $signalfile does not contain 'Name' annotation";
	$header_seg = $1;
	$name_index = ($header_seg =~ tr/\t/\t/+0);
	
	print OUT $header;
	
	while (<SIG>) {
		s/[\r\n]+$//;
		my @record = split (/\s+/, $_);
		my ($curname, $curlrr) = @record[$name_index, $lrr_index];
		if ($snpgc->{$curname}) {
			if ($curlrr =~ m/^[\d\-\.Ee]+$/) {
				$record[$lrr_index]  = $curlrr - $b1 - $b2 * $snpgc->{$curname}[2];
			} else {
				$record[$lrr_index]  = $curlrr;
			}
		} else {
			$count_nowave++;
		}
		print OUT join ("\t", @record), "\n";
	}
	close (SIG);
	close (OUT);
	
	$count_nowave and print STDERR "NOTICE: Total of $count_nowave markers cannot be adjusted due to lack of GC model information\n";
}
	

sub trainGCModel {
	my ($signalfile, $snpgc, $gap, $distance, $skip_marker) = @_;
	my ($lrr_index, $name_index, %alldata, $count);
	my ($header, $header_seg);
	

	open (SIG, $signalfile) or confess "Error: cannot read from file $signalfile: $!";
	$header = <SIG>;

	$header =~ m/(.*)Log R Ratio/ or confess "Error: the header file of signalfile $signalfile does not contain log r ratio annotation";
	$header_seg = $1;
	$lrr_index = ($header_seg =~ tr/\t/\t/+0);

	$header =~ m/(.*)Name/ or confess "error: the header file of signalfile $signalfile does not contain Chr annotation";
	$header_seg = $1;
	$name_index = ($header_seg =~ tr/\t/\t/+0);
	
	while (<SIG>) {
		s/[\r\n]+$//;
		my @record = split (/\s+/, $_);
		my ($curname, $curlrr) = @record[$name_index, $lrr_index];
		$snpgc->{$curname} or next;			#no annotation for this SNP in the SNP-GC file
		
		my ($curchr, $curpos, $curgc) = @{$snpgc->{$curname}};
		push @{$alldata{$curchr}}, [$curname, $curlrr, $curpos, $curgc];
		$count++;
	}
	
	my (@x, @y);
	for my $key (keys %alldata) {
		
		$key =~ m/^\d+$/ or next;			#the regression is only performed on autosomes !!!!!
		
		my @data = @{$alldata{$key}};
		my $previous = 0;
		@data = sort {$a->[2] <=> $b->[2]} @data;
		for my $i ($skip_marker .. @data-1-$skip_marker) {		#exclude the first a few or the last a few markers in each chromosome from training the model
			if ($distance) {
				if ($data[$i]->[2] - $previous > $distance) {
					if ($data[$i]->[3] > 15 and $data[$i]->[3] < 80 and $data[$i]->[1] > -1 and $data[$i]->[1] < 1) {
						push @x, $data[$i]->[3];
						push @y, $data[$i]->[1];
					}
					$previous = $data[$i]->[2];
					#print STDERR "NOTICE: chr-$key previous=$previous current=$data[$i]->[2]\n";
				}
				
			} elsif ($gap and $i % $gap == 0) {
				if ($data[$i]->[3] > 15 and $data[$i]->[3] < 80 and $data[$i]->[1] > -2) {
					push @x, $data[$i]->[3];
					push @y, $data[$i]->[1];
				}
			} else {
				pod2usage ("Error in argument: please specify a non-zero value of --distance or --gap");
			}
		}
	}
	print STDERR "NOTICE: Collecting ${\(scalar @x)} autosome SNPs for GC-based regression\n";
	
	my ($b1, $b2) = reg_linear (\@x, \@y);
	$verbose and print STDERR "NOTICE: Regression coefficient a=$b1 b=$b2\n";
	return ($b1, $b2);
}

sub calculateMADWF {
	my ($sigfile, $refchr) = @_;
	my ($wf, $cc);
	
	#ref_median can be GC percentage in 135 1MB non-overlapping sliding windows in chr11 (2004 or 2006 human genome assembly)
	my @ref_median = qw/54.8207535282258 56.8381472081218 53.1218950320513 46.9484174679487 39.9367227359694 38.3365384615385 41.9867788461538 40.4431401466837 44.5320512820513 42.1979166666667 41.6984215561224 43.1598557692308 43.4388020833333 40.8104967948718 39.8444475446429 41.5357572115385 38.7496995192308 45.0213249362245 42.3251201923077 43.5287459935897 40.7440808354592 37.0492788461538 36.5006009615385 35.8518016581633 35.2767427884615 35.1972155448718 36.5286192602041 39.4890825320513 36.5779246794872 36.7275641025641 38.3256935586735 37.791266025641 41.1777844551282 41.950534119898 42.3639823717949 41.9208733974359 41.2061543367347 35.4974959935897 35.2123397435897 36.5101841517857 36.7135416666667 36.8268229166667 37.6945153061224 40.7453926282051 47.7049278846154 47.3233173076923 44.7361288265306 46.6585536858974 39.1593549679487 36.5684789540816 38.2718466806667 37.184425 37.184425 37.184425 37.184425 35.9227764423077 41.1157852564103 41.6662348533163 39.7402844551282 40.0149238782051 46.6417211415816 49.9136618589744 45.2016225961538 51.3019172512755 52.0818309294872 51.1320112179487 49.9807185102302 49.9807185102302 49.5874473187766 50.547349024718 50.7186498397436 45.6435347576531 46.3352363782051 42.4091546474359 46.6399274553571 43.7746394230769 45.0160256410256 41.8526642628205 43.8899075255102 38.5112179487179 36.1038661858974 36.1689851721939 39.8506610576923 37.0439703525641 36.8012595663265 40.2521033653846 39.661858974359 37.5013769093564 35.5448717948718 36.9039979272959 35.2046274038462 38.2195512820513 40.074537627551 40.7097355769231 40.5470753205128 38.4104380072343 36.131109775641 35.3915264423077 34.9693080357143 36.2953725961538 37.9602363782051 39.1942362882653 37.4464142628205 36.8879206730769 35.7242588141026 36.7556202168367 37.0639022435897 40.6929086538462 38.385084502551 39.4121594551282 40.2410857371795 42.0772879464286 43.2935697115385 43.2345753205128 40.9113919005102 44.9575320512821 46.2513020833333 46.4753069196429 48.3886217948718 47.8520633012821 43.8001802884615 39.808274872449 44.5042067307692 38.3835136217949 44.9097177933673 45.5366586538462 41.7346754807692 39.2198461415816 41.9489182692308 44.3351362179487 42.7910754145408 42.3190104166667 42.0425681089744 47.0514787946429 45.3482603740699/;

	open (SIG, $sigfile) or confess "Error: cannot read from signal file $sigfile: $!";
	my ($header, $header_temp, $lrr_index, $chr_index, $pos_index);
	$header = <SIG>;
	$header =~ m/(.*)Log R Ratio/ or confess "Error: invalid record in header line of $sigfile: <Log R Ratio> not found\n";
	$header_temp = $1;
	$lrr_index = ($header_temp =~ tr/\t/\t/ + 0);
	$header =~ m/(.*)Chr/ or confess "Error: invalid record in header file: <Chr> not found\n";
	$header_temp = $1;
	$chr_index = ($header_temp =~ tr/\t/\t/ + 0);
	$header =~ m/(.*)Position/ or confess "Error: invalid record in header file: <Position> not found\n";
	$header_temp = $1;
	$pos_index = ($header_temp =~ tr/\t/\t/ + 0);

	my (%chrsig, @allsig, @alllrr, @allblock);
	while (<SIG>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		if (not defined $record[$chr_index] or not defined $record[$pos_index] or not defined $record[$lrr_index]) {
			confess "Error: invalid record found in $sigfile (Chr, Pos, Log R Ratio not found): <$_>\n";
		}
		$record[$chr_index] =~ m/^\d+$/ or next;			#only consider autosome
		lc $record[$lrr_index] eq 'nan' and next;			#ignore missing data (annotated as NaN or possibly nan in signal intensity file)
		lc $record[$lrr_index] eq '-infinity' and next;			#ignore missing data (annotated as NaN or possibly nan in signal intensity file)
		$record[$pos_index] =~ m/^\d+$/ or confess "Error: invalid position in $sigfile (positive integer expected): <$_>\n";
		if (not $record[$lrr_index] =~ m/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
			print STDERR "WARNING: skipping invalid Log R Ratio in $sigfile (floating point number expected): <$_>\n";
			next;
		}
		push @{$chrsig{$record[$chr_index]}}, [$record[$pos_index], $record[$lrr_index]];
	}
	close (SIG);
	
	my (@cursignal, %chryi, @yi);
	for my $nextchr (keys %chrsig) {
		my @allsig = sort {$a->[0] <=> $b->[0]} @{$chrsig{$nextchr}};	#sort by position
		my $current_bin = 0;
		for my $i (0 .. @allsig-1) {
			if ($allsig[$i]->[0] >= $current_bin * 1_000_000 and $allsig[$i]->[0] < ($current_bin+1) * 1_000_000) {
				push @cursignal, $allsig[$i]->[1];
			} else {
				$current_bin++;
				if (@cursignal and @cursignal > 10) {
					push @{$chryi{$nextchr}}, median (\@cursignal);
				} else {
					push @{$chryi{$nextchr}}, 'NA';
				}
				@cursignal = ();
			}
		}
		if (@cursignal and @cursignal > 10) {
			push @{$chryi{$nextchr}}, median (\@cursignal);
		} else {
			push @{$chryi{$nextchr}}, 'NA';
		}
	}
	
	for my $nextchr (keys %chryi) {
		push @yi, @{$chryi{$nextchr}};
	}
	
	print STDERR "NOTICE: Detected ${\(scalar @yi)} 1MB sliding windows in $sigfile";
	@yi = grep {!m/NA/} @yi;
	my $yi_median = median (\@yi);
	print STDERR " (${\(scalar @yi)} were used in WF calculation with median=$yi_median)\n";
	@yi = map {abs ($_ - $yi_median)} @yi;
	$wf = median (\@yi);

	if (not exists $chryi{$refchr}) {
		confess "ERROR: Unable to find reference chromosome chr$refchr from inputfile so waviness factor cannot be calculated\n";
	}
	$cc = calCC (\@ref_median, [@{$chryi{$refchr}}]);
	$cc > 0 and $wf = -$wf;
	return ($wf, $cc);
}

sub calCC {
	my ($array1, $array2) = @_;
	my (@newarray1, @newarray2);
	@$array1 == @$array2 or print STDERR "WARNING: Unequal dimensions of arrays (${\(scalar @$array1)} vs ${\(scalar @$array2)})\n";
	my $dimension = scalar @$array1;
	$dimension > @$array2 and $dimension = scalar @$array2;
	for my $i (0 .. $dimension-1) {
		$array1->[$i] eq 'NA' and next;
		$array2->[$i] eq 'NA' and next;
		push @newarray1, $array1->[$i];
		push @newarray2, $array2->[$i];
	}
	return cc (\@newarray1, \@newarray2);
}

#the following subroutine calculates the correlation coefficient
sub cc {
	my ($score1, $score2) = @_;
	return ssr ($score1, $score2) / sqrt ((ssr ($score1, $score1) * ssr ($score2, $score2)));
}

#the following subroutine calculates the ssr score, which is used in cc (correlation coefficient) calculation
sub ssr {
	my @score1 = @{$_[0]};
	my @score2 = @{$_[1]};
	my $mean1 = mean ($_[0]);
	my $mean2 = mean ($_[1]);
	my $product = 0;
	for my $i (0 .. @score1-1) {
		$product += $score1[$i] * $score2[$i];
	}
	return ($product - @score1 * $mean1 * $mean2);
}

sub mean {
	my ($score) = @_;
	@$score or confess "Error: NO VALUES for calculating mean";
	my $sum;
	for (@$score) {
		$sum += $_;
	}
	return $sum/@$score;
}

sub sd {
	my ($score) = @_;
	@$score > 1 or confess "Error: Less than 2 values are given for SD calculation";
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
	@$score or confess "Error: NO VALUES for calculating median";
	my @newscore = sort {$a<=>$b} @$score;
	if (@newscore % 2 == 0) {
		return ($newscore[@newscore/2-1]+$newscore[@newscore/2])/2;
	} else {
		return $newscore[@newscore/2];
	}
}

sub reg_linear {
	my ($x, $y) = @_;
	my ($a, $b);					#regression coefficient alpha, beta
	
	my ($sx, $sy, $st2, $mss, $rss, $tss);		#mean x, mean y, sum, model sum of squares, residual sum of squares, total sum of squares
	
	@$x and @$y or die "Error: input to reg_linear() should be two reference to arrays\n";
	@$x and @$y or die "Error: input arrays to reg_linear() should have equal dimensions, but the current dimensions are ${\(scalar @$x)} and ${\(scalar @$y)}, respectively.\n";
	
	for my $i (0 .. @$x-1) {
		$sx += $x->[$i];
		$sy += $y->[$i];
	}
	
	for my $i (0 .. @$x-1) {
		my $t = $x->[$i] - $sx/@$x;
		$st2 += $t*$t;
		$b += $t*$y->[$i];
	}
	$b /= $st2;
	$a = ($sy-$sx*$b)/@$x;
	return ($a, $b);
}

=head1 SYNOPSIS

 genomic_wave.pl [arguments] <inputfiles | --listfile file | --wflistfile file>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --calwf			calculate Waviness Factor (WF) for input signal files
 	    --adjust			adjust signals based on a wave model file
 	    
 	    --listfile <file>		listfile contains file names to be processed (one file per line)
 	    --gcmodelfile <file>	specify GC model file to be used
 	    --suffix <string>		suffix of output file name (default=adjusted) (for -adjust operation)
 	    --prefix <string>		prefix of output file name (usually path-to-directory) (for -adjust operation)
 	    --output <file>		output file name (for -calwf and -train operation)
 	    
 	    --distance <int>		minimum marker-marker distance for training model (default=1Mb)
 	    --gap <int>			minimum marker-marker gap for training model

 Function: calculate or adjust for genomic waves in signal intensities for whole-genome SNP 
 genotyping arrays.

 Example: genomic_wave.pl -calwf -list list -output list.wf
 	  genomic_wave.pl -calwf inputfile
          genomic_wave.pl -adjust -gcmodel hh550.gcmodel inputfile

 Version: $LastChangedDate: 2010-07-08 04:34:07 -0700 (Thu, 08 Jul 2010) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--calwf>

calculate waviness factor for input signal files

=item B<--adjust>

adjust signal values for given input signal files, based on a given GC model 
file

=item B<--listfile>

a file containing the list of file names to be processed by this program, with 
one file name per line.

=item B<--gcmodelfile>

a GC model file containing the GC content of the genomic regions around each 
marker (columns are marker, chr, position and GC percentage)

=item B<--suffix>

specify the file name suffix for signal files after signal adjustment (by 
default the 'adjusted' will be appended to the input file name as the output file 
name)

=item B<--prefix>

specify the prefix of output file name for the --adjust operation. By default 
the output file name is the input file name with the 'adjusted' suffix; however, 
in many cases, one may want to write the output file to a new directory other 
than the original directory. In such cases the --prefix argument can be used to 
override the directory path for the input files, and generate output files in 
the desired output directories.

=item B<--output>

output file name for the --calwf operation. It has no effect on the -adjust 
operation, as the output file name will be determined by input file name as well 
as the --suffix argument

=item B<--ref_chromosome>

The reference chromosome for computing the sign of the WF. By default chromosome 
11 is used. User should not change this parameter at this moment. (in the future 
version this argument might be a useful)

=item B<--distance>

The minimum distance between any two markers in the regression model. This 
parameter is used to reduce the potential dependency between markers, and the 
default value is 1Mb, resulting in 2500-3000 points in a regression model for 
human genome.

=item B<--gap>

The minimum gap between any two markers in the regression model. The "gap" 
actually refers to the number of other markers between the two adjacent marker 
used in regression model. In general, one should consider using --distance 
argument rather than --gap, since many arrays do not contain evenly spaced 
markers.

=back

=head1 DESCRIPTION

This program is used to calculate the magnitude or adjust the effect of Genomic 
Waves. Genomic waves refers to the obvious wavy patterns that can be observed in 
many genomics applications that measure signal intensities along chromosomes; 
for example, hybridization intensity ratios for array-CGH experiments, log2Ratio 
signals from Affymetrix SNP genotyping assay and the log R Ratio values from the 
Illumina SNP genotyping assay.

The exact technical cause of genomic wave is not fully understood, however, 
empirically this is typically caused by too much or too little DNAs used in any 
hybridization experiments. The more deviated the quantity of DNA from normal 
amount (recommended amount, usually used in constructing canonical genotype 
clusters), the more obvious the wavy pattern is. Given a sample affected by 
genomic waves, we can observe that the wave patterns highly correlate with GC 
content, which forms the basis of the signal adjustment model used in this 
program.

=over 8

=item * B<calculation of waviness factor (WF) and GCWF>

The waviness factor (WF) is a value (typically between -0.15 and 0.1) that measure 
the strength of waviness in a given sample. The GCWF is a value that measures 
the proportion of waves that can be explained by genome-wide GC distribution 
patterns, which can be a good measure of the performance of the regression 
model. It is calculated simply as the median of the absolute values of median of 
signals for all 1Mb windows along the autosome. The sign of the WF is determined 
by comparing its signal patterns in a reference chromosome (chr11 by default) 
with a sample defined as positive sample.

The --calwf argument is used to specify that WF be calculated for several input 
files (either specify all their file names in the command line, or put their 
file names into a text file and use the --listfile argument). For example:

	genomic_wave.pl -calwf -list list -output list.wf

tells the program to calculate WF values for all signal file whose names are 
specified in the list1 file, and write the output to the list1.wf file.Each 
signal file should contain information for one marker per line, and each line 
contains tab-delimited records. For example, a typical signal file is below:

	Name    Chr     Position        1784.GType      1784.B Allele Freq      1784.Log R Ratio
	rs3094315       1       792429  AB      0.4678481       -0.03128142
	rs12562034      1       808311  BB      1       -0.0216024
	rs3934834       1       1045729 AA      0       -0.1426667
	rs9442372       1       1058627 AA      0.01218348      -0.1583388
	rs3737728       1       1061338 BB      1       0.07487191

This program will first read the first line of the signal file to determine 
which tab-delimited column corresponds to the 'Log R Ratio' values, and then use 
this column in subsequent lines for the calculation. The signal file need not be 
sorted by chromosome and positions, although the first line (header line) must 
contains a column ending with 'Log R Ratio' so that the program knows where to 
find the LRR values for calculatoin.

The output contains several tab-delimited records per line and each line 
corresponds to one file name. The three columns in each line indicate the file 
name, the GCWF value, the WF value, and the correlation coefficient of 1MB 
window signal in chr11 with a canonical positive samples (the sign of the 
correlation coefficient determines the sign of the WF value). An example file is 
shown below:

	bp3.5930.adj200kb       0.00531254037319981     0.021957369778586       -0.241947939428559
	bp3.5931.adj200kb       0.00488808609688766     0.0129159837157837      -0.378452482168609
	bp3.5932.adj200kb       0.000475109412672177    0.0176191800877771      -0.026965466628142
	bp3.5933.adj200kb       0.00194776631711844     0.0126938462557015      -0.153441776265692
	bp3.5934.adj200kb       0.0031053760715499      0.0155745026279633      -0.199388458542126

=item * B<adjust signal files based on a wave model file>

The -adjust argument is used for reducing the genomic waves, and it requires 
user to specify a GC model file using the --gcmodel argument. For example:

	genomic_wave.pl -adjust -gcmodel hhall.gcmodel -list list

The above command specifies that we want to apply the regression model in 
hhall.gcmodel file and adjust signal values for all files in the list file 
(one file name per line). The new file name will be the old file name plus a 
"adjusted" suffix, and the "Log R Ratio" column in the old file will be replaced 
by the new LRR values.

For example, the first a few lines in an old un-adjusted signal file is:

	Name    Chr     Position        1784.GType      1784.B Allele Freq      1784.Log R Ratio
	rs3094315       1       792429  AB      0.4678481       -0.03128142
	rs12562034      1       808311  BB      1       -0.0216024
	rs3934834       1       1045729 AA      0       -0.1426667
	rs9442372       1       1058627 AA      0.01218348      -0.1583388
	rs3737728       1       1061338 BB      1       0.07487191

The first a few lines in the adjusted signal file is:

	Name    Chr     Position        1784.GType      1784.B Allele Freq      1784.Log R Ratio
	rs3094315       1       792429  AB      0.4678481       0.0174676876557617
	rs12562034      1       808311  BB      1       -0.0573741457755519
	rs3934834       1       1045729 AA      0       0.0998165381118793
	rs9442372       1       1058627 AA      0.01218348      0.119878665107063
	rs3737728       1       1061338 BB      1       0.315737089180156

You can see that the only different between the old and the new file is the Log 
R Ratio column (last tab-delimited column).

After running the signal adjustment, one may re-run the --calwf operation on 
the new files, and then compare the GCWF values with those of the un-adjusted 
files to see whether there is really reduction in waviness patterns.

=back

For questions, bugs or comments, please email kai@openbioinformatics.org.

=cut                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                