#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;

our $VERSION = 			'$Revision: 0d69c478fd8b9e0a33c63ebda6bdc4ccf65bcdb3 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2013-04-22 21:03:06 -0700 (Mon, 22 Apr 2013) $';

our ($verbose, $help, $man);
our ($operator, $inputfile);
our ($output, $signalfile, $fraction, $bp);

GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'output=s'=>\$output, 'signalfile=s'=>\$signalfile, 'fraction=f'=>\$fraction,
	'bp'=>\$bp) or pod2usage ();


$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error:");

($operator, $inputfile) = @ARGV;

my ($marker_pos);


if (defined $output) {
	open (STDOUT, ">$output") or confess "Error: cannot write to output file $output: $!\n";
}

if (defined $fraction) {
	$fraction > 0 and $fraction < 1 or pod2usage ("Error in argument: the --fraction must be between 0 and 1");
} else {
	print STDERR "NOTICE: the --fraction argument is set as 0.2 by default\n";
	$fraction = 0.2;
}

if ($operator eq 'combineseg') {
	defined $signalfile or pod2usage ("Error in argument: please specify the --signalfile argument");
	my $posindex = readMarkerPos ($signalfile);
	combineSegment ($inputfile, $posindex);
}



sub combineSegment {
	my ($cnvfile, $posindex) = @_;
	my (%cnvcall);
	my ($skipped_line);
	
	if ($cnvfile eq 'stdin') {
		*CNV = *STDIN;
	} else {
		open (CNV, $cnvfile) or confess "Error: cannot read from inputfile $cnvfile: $!\n";
	}
	
	while (<CNV>) {
		if (m/^(?:chr)?(\w+):(\d+)-(\d+)\s+numsnp=(\d+)\s+length=(\S+)\s+state(\d),cn=(\d)\s+(.+?)\s+startsnp=(\S+)\s+endsnp=(\S+)/) {
			my ($curchr, $curstart, $curend, $curnumsnp, $curlength, $curstate, $curcn, $curfile, $cursnpstart, $cursnpend) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10);
			my $curconf;
			m/conf=([\-\.\deE]+)/ and $curconf=$1;
			$curlength =~ s/,//g;
			
			defined $posindex->{$cursnpstart} and defined $posindex->{$cursnpend} or die "Error: the index for SNPs ($cursnpstart and $cursnpend) are not found from signalfile\n";
			push @{$cnvcall{$curfile, $curchr}}, [$curstart, $curend, $posindex->{$cursnpstart}, $posindex->{$cursnpend}, $cursnpstart, $cursnpend, $curcn, $curconf];
			
		} else {
			$skipped_line++;
		}
	}
	$skipped_line and print STDERR "WARNING: $skipped_line lines were skipped due to unrecognizable formats\n";
	
	for my $key (keys %cnvcall) {
		my @call = sort {$a->[0] <=> $b->[0]} @{$cnvcall{$key}};
		my @newcall;
		
		my ($curfile, $curchr) = split (/$;/, $key);
		my ($prestart, $preend, $prestartindex, $preendindex, $presnpstart, $presnpend, $precn, $preconf) = @{$call[0]};
		
		
		
		push @newcall, [@{$call[0]}];
		for my $i (1 .. @call-1) {
			my ($curstart, $curend, $curstartindex, $curendindex, $cursnpstart, $cursnpend, $curcn, $curconf) = @{$call[$i]};
			if ($curcn eq $precn) {
				if ($bp and ($curstart-$preend-2)/($curend-$prestart+1) <= $fraction or not $bp and ($curstartindex-$preendindex-1) / ($curendindex-$prestartindex+1) <= $fraction) {
					if ($bp) {
						print STDERR "NOTICE: Merging chr$curchr:$prestart-$preend and chr$curchr:$curstart-$curend with base pair fraction=", sprintf("%.3f", ($curstart-$preend-1) / ($curend-$prestart)), " ($curstart-$preend-1) / ($curend-$prestart)\n";
					} else {
						print STDERR "NOTICE: Merging chr$curchr:$prestart-$preend and chr$curchr:$curstart-$curend with marker fraction=", sprintf("%.3f", ($curstartindex-$preendindex-1) / ($curendindex-$prestartindex)), " ($curstartindex-$preendindex-1) / ($curendindex-$prestartindex)  ($curstart-$preend-1) / ($curend-$prestart)\n";
					}
					pop @newcall;
					my $newconf;
					defined $preconf and defined $curconf and $newconf = $preconf+$curconf;
					push @newcall, [$prestart, $curend, $prestartindex, $curendindex, $presnpstart, $cursnpend, $curcn, $newconf];
					($prestart, $preend, $prestartindex, $preendindex, $presnpstart, $presnpend, $precn, $newconf) = @{$newcall[$#newcall]};
					next;
				}
			}
			
			push @newcall, [@{$call[$i]}];
			($prestart, $preend, $prestartindex, $preendindex, $presnpstart, $presnpend, $precn, $preconf) = @{$call[$i]};
		}
		for my $i (0 .. @newcall-1) {
			my ($curstart, $curend, $curstartindex, $curendindex, $cursnpstart, $cursnpend, $curcn, $curconf) = @{$newcall[$i]};
			my $cnvregion = "chr$curchr:$curstart-$curend";
			my $curnumsnp = $curendindex-$curstartindex+1;
			my $curstate = $curcn+1; $curcn >= 2 and $curstate++;
			if (length ($cnvregion) < 29) {
				$cnvregion = substr ("$cnvregion                              ", 0, 29);
			}
		
			if (length($curnumsnp) < 6) {
				$curnumsnp = substr ("$curnumsnp      ", 0, 6);
			}
			
			my $cnvlength = join ('', reverse split (//, $curend-$curstart+1)); $cnvlength =~ s/(...)/$1,/g; $cnvlength =~ s/,$//; $cnvlength = join ('', reverse split (//, $cnvlength));
			if (length($cnvlength) < 11) {
				$cnvlength = substr ("$cnvlength            ", 0, 11);
			}
			
			print "$cnvregion numsnp=$curnumsnp length=$cnvlength state$curstate,cn=$curcn $curfile startsnp=$cursnpstart endsnp=$cursnpend";
			print STDOUT (defined $curconf)?" conf=$curconf":"", "\n";
		}
	}
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
	$header =~ m/(.*)Name/ or confess "error: the header file of signalfile $signalfile does not contain Name annotation";
	$header_seg = $1;
	$name_index = ($header_seg =~ tr/\t/\t/+0);
	
	my %alldata;		#key is chromosome, value is a [pos, name] array
	my $count;
	my (%index);
	while (<SIG>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		my ($curname, $curchr, $curpos) = @record[$name_index, $chr_index, $pos_index];
		push @{$alldata{$curchr}}, [$curpos, $curname];
		$count++;
	}
	close (SIG);
	for my $key (keys %alldata) {
		@{$alldata{$key}} = sort {$a->[0] <=> $b->[0]} @{$alldata{$key}};
		for my $i (0 .. @{$alldata{$key}}-1) {
			$index{$alldata{$key}->[$i][1]} = $i; 
			#$alldata{$key}->[$i][1] eq 'rs375087' and print STDERR "NOTICE: found i=$i pos=$alldata{$key}->[$i][0] rs375087\n";
		}
	}
	print STDERR "NOTICE: Total of $count records are read from $signalfile\n";
	return \%index;
}

=head1 SYNOPSIS

 clean_cnv.pl [arguments] <operation> <cnv-call-file>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	
 	    --output <file>		write output to this file (default=STDOUT)
 	    --signalfile <file>		a file that contains Chr and Position for SNPs/markers
 	    --fraction <float>		maximum fraction of gap divided by combined CNV (default: 0.2)
 	    --bp			use base pairs to calculate fraction (default: number of markers)

 Function: post-process CNV calls. Current <operation> can be only 'combineseg'

 Example: clean_cnv.pl combineseg input.rawcnv > output.rawcnv

 Version: $LastChangedDate: 2013-04-22 21:03:06 -0700 (Mon, 22 Apr 2013) $

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

=item B<--signalfile>

a file that annotates the Chr and Position for each SNP markers. Normally, it is 
easiest to just use the signal intensity file that contains these measures (one 
marker per line).

=item B<--fraction>

maximum fraction of gap for merging two neighboring CNV calls into one CNV call. 
By default, the CNV length is measured by number of SNPs/probes. For example, 
there are three segments A, B and C. A and C are called as deletion. If the 
length of B divided by the length of A+B+C is less than this threshold, then 
A+B+C will be called as a single deletion.

=item B<--bp>

use base pairs (rather than number of SNPs/probes) to calculate the fraction 
measure.

=back

=head1 DESCRIPTION

This program is designed for clean CNV calls generated by the PennCNV program 
by various criteria.

the combineseg operation combined neighboring CNV calls if the gap between them 
is less than 20% (use --fraction to change the default) of the total length 
(either counted by number of mrakers or by number of base pairs) of combined 
CNVs.

For questions, comments or bug reports, please contact me at 
kai@openbioinformatics.org.

=cut