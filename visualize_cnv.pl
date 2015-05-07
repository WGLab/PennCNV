#!/usr/bin/env perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;
BEGIN {
	($_=$0)=~s{[^\\\/]+$}{};
	$_||="./";
	$ENV{PATH} .= ":$_";
}

our $VERSION = 			'$Revision: d973451fd0646e545c8537ef91cd10f3f5cc9eff $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2012-10-23 23:32:05 -0700 (Tue, 23 Oct 2012) $';

our ($verbose, $help, $man);
our ($format, $track_name, $bedstrand, $intype, $familyfile, $minoverlap, $commonfile, $highlightfile, $snpposfile, $output, $nfperfile, $idmapfile, $signalfile, $flankinglength, $pdfout, $keeptemp);
our ($inputfile);

GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'format=s'=>\$format, 'track_name=s'=>\$track_name, 'bedstrand=s'=>\$bedstrand,
	'intype=s'=>\$intype, 'familyfile=s'=>\$familyfile, 'minoverlap=f'=>\$minoverlap, 'commonfile=s'=>\$commonfile, 'highlightfile=s'=>\$highlightfile,
	'snpposfile=s'=>\$snpposfile, 'output=s'=>\$output, 'nfperfile=i'=>\$nfperfile, 'idmapfile=s'=>\$idmapfile, 'signalfile=s'=>\$signalfile, 'flankinglength=s'=>\$flankinglength,
	'pdfout'=>\$pdfout, 'keeptemp'=>\$keeptemp) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error: <inputfile> missing");


($inputfile) = @ARGV;
$format ||= 'bed';			#by default, convert CNV calls to Genome Browser track
$intype ||= 'cnv';
$bedstrand ||= '.';
$bedstrand eq '+' or $bedstrand eq '-' or $bedstrand eq '.' or pod2usage ("Error: the --bedstrand can be only + or - or .");
$minoverlap ||= 0;
$flankinglength ||= 50000;
$flankinglength =~ s/m$/000000/;
$flankinglength =~ s/k$/000/;
$flankinglength =~ m/^\d+$/ or pod2usage ("Error in argument: the --flankinglength argument must be a positive integer (suffix of 'k' or 'm' are fine)");


if ($output and not $nfperfile) {	#when --nfperfile is set, we need to write to a file (rather than STDOUT)
	open (STDOUT, ">$output") or die "Error: cannot write to output file $output: $!\n";
}

if ($intype eq 'cnv') {
	if ($format eq 'bed') {
		convertCNVToBED ($inputfile, $track_name, $bedstrand);
	} elsif ($format eq 'html') {
		$familyfile or pod2usage ("Error in argument: please specify the --familyfile argument for the HTML output");
		if (not defined $output and $nfperfile) {	#when --nfperfile is set, the output must be written to a file, rather than STDOUT
			$output = "$inputfile.output";
			print STDERR "NOTICE: The --output argument is automatically set as $output\n";
		}
		convertCNVToHTML ($inputfile, $familyfile, $output);
	} elsif ($format eq 'beadstudio') {
		$idmapfile or pod2usage ("Error in argument: please specify the --idmapfile argument for BeadStudio output");
		convertCNVToBeadStudio ($inputfile, $idmapfile);
	} elsif ($format eq 'tab') {
		convertCNVToTab ($inputfile);
	} elsif ($format eq 'plot') {
		$signalfile or pod2usage ("Error in argument: please specify the --signalfile argument (which contains LRR/BAF for all markers) for signal intensity plot output");
		generateCNVPlot ($inputfile, $signalfile, $signalfile);
	} elsif ($format eq 'idetab') {
		$signalfile or pod2usage ("Error in argument: please specify the --signalfile argument (which contains LRR/BAF/Chr/Location for all markers) for Idetab output");
		convertCNVToIdetab ($inputfile, $signalfile);
	} else {
		pod2usage ("Error in argument: for --intype of cnv, the valid --format arguments are bed, html, idetab and beadstudio only");
	}
} elsif ($intype eq 'assoc') {
	if ($format eq 'wig') {
		convertAssocToWIG ($inputfile, $snpposfile);
	} elsif ($format eq 'html') {
		pod2usage ("Error in argument: the functionality of converting association results to HTML file is to be implemented");
		convertAssocToHTML ($inputfile);
	}
} elsif ($intype eq 'tab') {
	if ($format eq 'cnv') {
		convertTabToCNV ($inputfile);
	}
} else {
	pod2usage ("Error in argument: the --intype can be only 'cnv', 'tab' or 'assoc'");
}

sub generateCNVPlot {
	my ($cnvfile, $signalfile, $sampleid) = @_;
	$verbose and print STDERR "NOTICE: plotting the CNV calls from file $cnvfile and generating JPG outputs (via R scripts)\n";

	if ($cnvfile eq 'stdin') {
		*REGION = *STDIN;
	} else {
		open (REGION, $cnvfile) or die "Error: cannot read from region file $cnvfile: $!\n";
	}
	
	my (%chrregion, $skipped_line);
	while (<REGION>) {
		if (m/^chr((\w+):(\d+)-(\d+))\s+numsnp=(\d+)\s+length=([\d\,]+)\s+state\d+,cn=(\d+)\s+(.+?)\s+startsnp=(\S+)\s+endsnp=(\S+)/) {
			$sampleid eq $8 or next;
			push @{$chrregion{$2}}, [$3, $4, $7];
		} else {
			$skipped_line++;
		}
	}
	$skipped_line and print STDERR "WARNING: $skipped_line lines in $cnvfile are skipped due to unrecognizable format\n";
	for my $chr (keys %chrregion) {
		@{$chrregion{$chr}} = sort {$a->[0] <=> $b->[0]} @{$chrregion{$chr}};
	}
		
	open (FH, $signalfile) or die "Error: cannot read from signalfile $signalfile: $!\n";
	$_ = <FH>;
	defined ($_) or die "\nERROR: NOTHING is found in signalfile $signalfile. Please check the file before proceeding\n";
	s/[\r\n]+$//;
	my @header = split (/\t/, $_);
	my ($name_index, $chr_index, $pos_index, $lrr_index, $baf_index);
	for my $i (0 .. @header-1) {
		if ($header[$i] eq 'Name') {
			$name_index = $i;
		} elsif ($header[$i] eq 'Chr' or $header[$i] eq 'Chromosome') {
			$chr_index = $i;
		} elsif ($header[$i] eq 'Position') {
			$pos_index = $i;
		} elsif ($header[$i] =~ m/Log R Ratio$/ or $header[$i] eq 'LRR') {
			$lrr_index = $i;
		} elsif ($header[$i] =~ m/B Allele Freq$/ or $header[$i] eq 'BAF') {
			$baf_index = $i;
		}
	}

	my (%signal);
	defined $name_index or die "Error: the signal file $signalfile does not contain the Name column in the header line\n";
	if (not defined $snpposfile) {
		defined $chr_index and defined $pos_index or die "Error: the signal file $signalfile does not contain Chr and Position column, use --snpposfile to supply this information\n";
	}
	defined $name_index and defined $lrr_index and defined $baf_index or die "Error: the signalfile $signalfile does not contain the Name, LRR and BAF columns in the first line\n";
	
	my ($name_chr, $name_pos);
	$snpposfile and ($name_chr, $name_pos) = readSNPPosFile ($snpposfile);
	
	while (<FH>) {
		s/[\r\n]+$//;
		my @field = split (/\t/, $_);
		my ($name, $chr, $pos, $lrr, $baf);
		if (defined $chr_index and defined $pos_index) {
			($name, $chr, $pos, $lrr, $baf) = @field[$name_index, $chr_index, $pos_index, $lrr_index, $baf_index];
		} else {
			($name, $lrr, $baf) = @field[$name_index, $lrr_index, $baf_index];
			($chr, $pos) = ($name_chr->{$name}, $name_pos->{$name});
			#$name eq 'rs12106502' and print STDERR "name=$name, $lrr, $baf, $chr, $pos\n";
		}
		defined $name and defined $chr and defined $pos and defined $lrr and defined $baf or next;		#this record does not contain required information

		if (defined $chrregion{$chr}) {
			for my $i (0 .. @{$chrregion{$chr}}-1) {
				my ($start, $end, $cn) = @{$chrregion{$chr}->[$i]};
				my $length = $end-$start+1;
				$length < $flankinglength and $length = $flankinglength;			#set the minimum length to 10kb (markers s10kb before or after CNV call will be printed)
				if ($pos >= $start-$length and $pos <= $end+$length) {
					push @{$signal{$chr, $start}}, [$pos, $lrr, $baf];
				}
			}
		}
	}
	close (FH);
	print STDERR "NOTICE: Signal values for ", scalar (keys %signal) > 1 ? (scalar (keys %signal) . " CNV regions are") : (scalar (keys %signal) . " CNV region is"), " found in $signalfile\n";
	
	for my $chr (sort keys %chrregion) {
		for my $i (0 .. @{$chrregion{$chr}}-1) {
			my ($start, $end, $cn) = @{$chrregion{$chr}->[$i]};
			my ($start1, $end1) = ($start, $end);		#used in plotting the two vertical grey lines
			defined $signal{$chr, $start} or print STDERR "NOTICE: Skipping chr$chr:$start region due to lack of signal information\n" and next;
			my @signal = sort {$a->[0] <=> $b->[0]} @{$signal{$chr, $start}};
			my $megabp = "bp";
			my $cex = 0.3;					#by default, 0.3 is used in plotting the dots
			my (@lrr);
			
			if (@signal > 4000) {				#decrease cex to avoid cluttering the screen, when whole-chromosome is being printed
				$cex = 0.01;
			} elsif (@signal>2000) {			#about >10M in 550K array
				$cex = 0.05;
			} elsif (@signal>400) {				#about >2M in 550K array
				$cex = 0.1;
			} elsif (@signal>100) {				#about >500kb in 550K array
				$cex = 0.2;
			}
			print STDERR "NOTICE: Processing sample $sampleid CNV chr$chr:$start-$end with copy number of $cn ... ";
			open (SIGNAL, ">$sampleid.chr$chr.$start.signal") or die "Error: cannot write to output file $sampleid.chr$chr.$start.signal: $!\n";
			if ($signal[$#signal]->[0] > 1_000_000) {
				$megabp = "Mb";
				$start1 /= 1_000_000;
				$end1 /= 1_000_000;
			}
			for my $j (0 .. @signal-1) {
				if ($megabp eq "Mb") {
					$signal[$j]->[0] /= 1_000_000;
				}
				$signal[$j]->[0] >= $start1 and $signal[$j]->[0] <= $end1 and $signal[$j]->[1] =~ m/^[\-\d\.e]+$/ and push @lrr, $signal[$j]->[1];		#store all LRR values
				print SIGNAL join ("\t", @{$signal[$j]}), "\n";
			}
			close (SIGNAL);
			
			@lrr = sort {$a<=>$b} @lrr;			#sort LRR from small to large values
			my $lrr_bound = @lrr[@lrr*0.01];		#if <100 markers, take the lowest values; otherwise, take the 1% percentile of LRR values as the lower bound
			$lrr_bound > -1 and $lrr_bound = -1;		#make it lower bound at -1.
			
			open (R, ">$sampleid.chr$chr.$start.R") or die "Error: cannot write to output file $sampleid.chr$chr.$start.R: $!\n";
			if ($pdfout) {
				print R <<END;
one=read.table("$sampleid.chr$chr.$start.signal",sep="\\t")
pdf(file="$sampleid.chr$chr.$start.PDF", width=3,height=3, pointsize=6)
par(mfrow=c(2,1),las=1,mar=c(4,4.2,1.5,1.3))
END
			} else {
			print R <<END;
one=read.table("$sampleid.chr$chr.$start.signal",sep="\\t")
bitmap(file="$sampleid.chr$chr.$start.JPG",type="jpeg",width=3,height=3,res=300,pointsize=6)
par(mfrow=c(2,1),las=1,mar=c(4,4.2,1.5,1.3))
END
			}

			print R <<END;
a=rep("blue",length(one\$V1))
for(i in 1:length(one\$V1)){
	if(one\$V1[i]>=$start1 & one\$V1[i]<=$end1){
		a[i]="red"
	}
}
plot(one\$V1,one\$V2,pch=20,cex=$cex,ylim=c($lrr_bound,1),ylab="Log R Ratio",xlab="Chr$chr Position ($megabp)",col=a)
title(paste("$sampleid","(CN=$cn)"))
abline(v=c($start1,$end1),col="grey",lty=1,lwd=0.3)
abline(h=0,col="grey",lty=2,lwd=0.3)
plot(one\$V1,one\$V3,pch=20,cex=$cex,ylim=c(0,1),ylab="B Allele Freq",xlab="Chr$chr Position ($megabp)",col=a)
abline(v=c($start1,$end1),col="grey",lty=1,lwd=0.3)
END

			if ($cn == 2) {
				print R q/abline(h=0.5,col="grey",lty=2,lwd=0.3)/, "\n";
			} elsif ($cn == 3) {
				print R q/abline(h=c(0.33,0.67),col="grey",lty=2,lwd=0.3)/, "\n";
			} elsif ($cn == 4) {
				print R q/abline(h=c(0.25,0.5,0.75),col="grey",lty=2,lwd=0.3)/, "\n";
			}
			print R "dev.off()\n";
			close (R);
			system ("R CMD BATCH $sampleid.chr$chr.$start.R") and die "Error: cannot execute system command R CMD BATCH $sampleid.chr$chr.$start.R\n";
			if (not $keeptemp) {
				unlink ("$sampleid.chr$chr.$start.R","$sampleid.chr$chr.$start.Rout", "$sampleid.chr$chr.$start.signal");
				my $routfile = "$sampleid.chr$chr.$start.Rout";		#Rout file is typically written to current directory, which may differ from the signal file directory
				if ($routfile =~ m/.*[\\\/]([\\\/]+)/) {
					unlink ($1);
				}
			}
			print STDERR " written to $sampleid.chr$chr.$start.", $pdfout?"PDF\n":"JPG\n";
		}
	}
}

sub readSNPPosFile {
	my ($snpposfile) = @_;
	my (%name_pos, %name_chr);
	my ($header, $header_seg, $lrr_index, $pos_index, $chr_index, $name_index);
	
	open (SIG, $snpposfile) or die "Error: cannot read from SNPLOC file $snpposfile: $!\n";
	print STDERR "NOTICE: Start reading snpposfile $snpposfile ...";
	$header = <SIG>;
	$header =~ m/(.+)Pos/ or die "error: the header line of $snpposfile does not contain Pos annotation";
	$header_seg = $1;
	$pos_index = ($header_seg =~ tr/\t/\t/+0);
	$header =~ m/(.+)Chr/ or die "error: the header file of $snpposfile does not contain Chr annotation";
	$header_seg = $1;
	$chr_index = ($header_seg =~ tr/\t/\t/+0);
	$header =~ m/(.*)Name/ or die "error: the header file of $snpposfile does not contain Name annotation";
	$header_seg = $1;
	$name_index = ($header_seg =~ tr/\t/\t/+0);

	my @header = split (/\t/, $header);
	@header>=3 or die "Error: the -snpposfile $snpposfile must be a tab-delimited file\n";
	$pos_index == $chr_index and die "Error: the -snpposfile $snpposfile must be tab-delimited with Pos and Chr field\n";
	$pos_index == $name_index and die "Error: the -snpposfile $snpposfile must be tab-delimited with Pos and Name field\n";
	$chr_index == $name_index and die "Error: the -snpposfile $snpposfile must be tab-delimited with Chr and Name field\n";

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
	
sub convertCNVToIdetab {
	my ($cnvfile, $signalfile) = @_;
	$verbose and print STDERR "NOTICE: converting the CNV call file $cnvfile to tab-delimited format for the IdeiogramBrowser software (Chr, Start, End, Name, Value)\n";

	if ($cnvfile eq 'stdin') {
		*REGION = *STDIN;
	} else {
		open (REGION, $cnvfile) or die "Error: cannot read from region file $cnvfile: $!";
	}
	
	my ($sampleid, $idwarning, %chrregion, $skipped_line);
	while (<REGION>) {
		if (m/^chr((\w+):(\d+)-(\d+))\s+numsnp=(\d+)\s+length=([\d\,]+)\s+state\d+,cn=(\d+)\s+(.+?)\s+startsnp=(\S+)\s+endsnp=(\S+)/) {
			defined $sampleid or $sampleid = $8;
			if ($sampleid ne $8) {
				$idwarning or print STDERR "WARNING: Multiple samples are included in inputfile, but only CNVs on the first sample ($sampleid) will be in output\n";
				$idwarning++;
			}
			push @{$chrregion{$2}}, [$3, $4, $7];
		} else {
			$skipped_line++;
		}
	}
	$skipped_line and print STDERR "WARNING: $skipped_line lines in $cnvfile are skipped due to unrecognizable format\n";
	for my $chr (keys %chrregion) {
		@{$chrregion{$chr}} = sort {$a->[0] <=> $b->[0]} @{$chrregion{$chr}};
	}
		
		
	open (FH, $signalfile) or die "Error: cannot read from signalfile $signalfile: $!\n";
	$_ = <FH>;
	defined ($_) or die "\nERROR: NOTHING is found in signalfile $signalfile. Please check the file before proceeding\n";
	s/[\r\n]+$//;
	my @header = split (/\t/, $_);
	my ($name_index, $chr_index, $pos_index);
	for my $i (0 .. @header-1) {
		if ($header[$i] eq 'Name') {
			$name_index = $i;
		} elsif ($header[$i] eq 'Chr' or $header[$i] eq 'Chromosome') {
			$chr_index = $i;
		} elsif ($header[$i] eq 'Position') {
			$pos_index = $i;
		}
	}

	defined $name_index and defined $chr_index and defined $pos_index or die "Error: the signalfile $signalfile does not contain the Name, Chr and Position columns\n";
	while (<FH>) {
		s/[\r\n]+$//;
		my @field = split (/\t/, $_);
		my ($name, $chr, $pos) = @field[$name_index, $chr_index, $pos_index];
		print $chr, "\t", $pos, "\t", $pos, "\t";
		my $foundcn;
		if (defined $chrregion{$chr}) {
			for my $i (0 .. @{$chrregion{$chr}}-1) {
				$foundcn and last;
				my ($start, $end, $cn) = @{$chrregion{$chr}->[$i]};
				if ($pos >= $start and $pos <= $end) {
					$foundcn++;
					if ($cn >= 2) {
						print "gain\t$cn\n";
					} else {
						print "loss\t$cn\n";
					}
				}
			}
		}
		if (not $foundcn) {
			print "no_aberration\t2\n";
		}
	}
}

sub convertCNVToTab {
	my ($cnvfile) = @_;
	print STDERR "NOTICE: converting the CNV call file $cnvfile to tab-delimited format (Chr, Start, End, CN, Sample, StartSNP, EndSNP, Conf, NumSNP)\n";
	my $skipped_line;
	if ($cnvfile eq 'stdin') {
		*REGION = *STDIN;
	} else {
		open (REGION, $cnvfile) or die "Error: cannot read from region file $cnvfile: $!";
	}
	while (<REGION>) {
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
	close (REGION);
	$skipped_line and print STDERR "WARNING: $skipped_line lines in $cnvfile are skipped due to unrecognizable format\n";
}

sub convertTabToCNV {
	my ($tabfile) = @_;
	if ($tabfile eq 'stdin') {
		*REGION = *STDIN;
	} else {
		open (REGION, $tabfile) or die "Error: cannot read from tab-delimited input file $tabfile: $!";
	}
	while (<REGION>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		@record == 9 or die "Error: invalid record found in tab-delimited CNV file (9 fields expected): <$_>\n";
		my ($chr, $start, $end, $cn, $sample, $startsnp, $endsnp, $conf, $numsnp) = @record;

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
		
		my $state=$cn+1; $cn>=2 and $state++;
		
		print "$cnvregion numsnp=$numsnp length=$cnvlength state$state,cn=$cn $sample startsnp=$startsnp endsnp=$endsnp";
		$conf ne 'NA' and print " conf=$conf";
		print "\n";
	}
}

		

sub convertCNVToBED {
	my ($cnvfile, $track_name, $bedstrand) = @_;
	print STDERR "NOTICE: converting the CNV call file $cnvfile to BED format to be displayed in UCSC genome browser\n";
	if ($cnvfile eq 'stdin') {
		*REGION = *STDIN;
	} else {
		open (REGION, $cnvfile) or die "Error: cannot read from region file $cnvfile: $!";
	}

	$track_name ||= "CNVs in $cnvfile";
	print qq/track name="Track: $track_name" description="$track_name" visibility=2 itemRgb="On"\n/;
	
	my %statecolor = (1=>'128,0,0', 2=>'255,0,0', 3=>'0,255,0', 4=>'0,255,0', 5=>'0,255,0', 6=>'0,128,0');
	my ($skipped_line);
	while (<REGION>) {
		if (m/^(?:chr)?(\w+):(\d+)-(\d+)\s+numsnp=(\d+)\s+length=\S+\s+state(\d)\S*\s+(\S+)/) {
			defined $1 and defined $2 and defined $4 and defined $5 and defined $6 and defined $statecolor{$5} or die "ERROR <$_>\n";
			print "chr$1\t", $2-1, "\t$3\t$6\t", $4*10, "\t$bedstrand\t0\t0\t$statecolor{$5}\n";
		} elsif (m/^(?:chr)?(\w+):(\d+)-(\d+)/) {			#the file contains only the chr_region, such as common CNV file
			print "chr$1\t$2\t$3\tCNV\n";
		} else {
			$skipped_line++;
		}
	}
	$skipped_line and print STDERR "WARNING: $skipped_line lines in $cnvfile are skipped due to unrecognizable format\n";
}

sub convertCNVToBeadStudio {
	my ($cnvfile, $idmapfile) = @_;
	print STDERR "NOTICE: converting the CNV call file $cnvfile into BeadStudio bookmark format to be imported into BeadStudio software\n";
	if ($cnvfile eq 'stdin') {
		*REGION = *STDIN;
	} else {
		open (REGION, $cnvfile) or die "Error: cannot read from region file $cnvfile: $!";
	}
	
	my (%sample_id);
	open (SAMPLEID, $idmapfile) or die "Error: cannot read from sampleidfile $idmapfile: $!\n";
	while (<SAMPLEID>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		@record == 2 or die "Error: invalid record in $idmapfile (2 tab-delimited records expected): <$_>\n";
		$sample_id{$record[0]} = $record[1];
	}

	
	my $time = scalar (localtime);
	my ($skipped_line, %skipped_file);
	my %bookmark_type = (0=>'CNV Bin: Min 0 To Max 0.5', 1=>'CNV Bin: Min 0.5 To Max 1.5', 2=>'CNV Bin: Min 1.5 To Max 2.5', 3=>'CNV Bin: Min 2.5 To Max 3.5', 4=>'CNV Bin: Min 3.5 To Max 4.5');

	print <<FILE_HEADER;
<Project_Bookmarks>
  <Version>2.0.0</Version>
  <Name>PennCNV calls in $cnvfile</Name>
  <Author></Author>
  <Comment><![CDATA[  ]]></Comment>
  <CreateDate>$time</CreateDate>
  <Algorithm>PennCNV</Algorithm>
  <AlgorithmVersion></AlgorithmVersion>
  <Module></Module>
  <GenomeSpecies></GenomeSpecies>
  <GenomeBuild ></GenomeBuild >
  <TableSource>GT Samples</TableSource>
  <Bookmark_Templates>
   <bookmark_template>
     <type>CNV Bin: Min 0 To Max 0.5</type>
     <fill_color>DarkRed</fill_color>
     <fill_style>Solid</fill_style>
     <fill_opacity>50</fill_opacity>
   </bookmark_template>
   <bookmark_template>
     <type>CNV Bin: Min 0.5 To Max 1.5</type>
     <fill_color>DarkOrange</fill_color>
     <fill_style>Solid</fill_style>
     <fill_opacity>50</fill_opacity>
   </bookmark_template>
   <bookmark_template>
     <type>CNV Bin: Min 1.5 To Max 2.5</type>
     <fill_color>DarkGreen</fill_color>
     <fill_style>Solid</fill_style>
     <fill_opacity>50</fill_opacity>
   </bookmark_template>
   <bookmark_template>
     <type>CNV Bin: Min 2.5 To Max 3.5</type>
     <fill_color>DarkBlue</fill_color>
     <fill_style>Solid</fill_style>
     <fill_opacity>50</fill_opacity>
   </bookmark_template>
   <bookmark_template>
     <type>CNV Bin: Min 3.5 To Max 4.5</type>
     <fill_color>BlueViolet</fill_color>
     <fill_style>Solid</fill_style>
     <fill_opacity>50</fill_opacity>
   </bookmark_template>
  </Bookmark_Templates>
  <Bookmarks>
FILE_HEADER

	while (<REGION>) {
		if (m/^(?:chr)?(\w+):(\d+)-(\d+)\s+numsnp=(\d+)\s+length=\S+\s+state\d,cn=(\d)\s+(\S+).*(conf=(\S+))?/) {
			my ($curchr, $curstart, $curend, $cursnp, $curcn, $curfile, $curconf) = ($1, $2, $3, $4, $5, $6, $8);
			$curconf ||= 0;
			
			if (not $sample_id{$curfile}) {
				$skipped_file{$curfile}++;
				next;
			}
			
			print <<BOOKMARK;
    <bookmark>
      <sample_id>$sample_id{$curfile}</sample_id>
      <bookmark_type>$bookmark_type{$curcn}</bookmark_type>
      <entry_date>$time</entry_date>
      <chr_num>$curchr</chr_num>
      <base_start_pos>$curstart</base_start_pos>
      <base_end_pos>$curend</base_end_pos>
      <author></author>
      <comment>
        <![CDATA[ CNV Value: $curcn; CNV Confidence: $curconf; CNV numsnp: $cursnp ]]>
      </comment>
    </bookmark>
BOOKMARK


		} else {
			$skipped_line++;
		}
	}
	print "  </Bookmarks>\n</Project_Bookmarks>\n";

	$skipped_line and print STDERR "WARNING: $skipped_line lines in $cnvfile are skipped due to unrecognizable format\n";
	%skipped_file and print STDERR "WARNING: ${\(scalar keys %skipped_file)} file names in $cnvfile are skipped due to lack of correspondence in idmapfile $idmapfile\n";

}




sub generateRandFile {
	my ($num) = @_;
	$num < 100 or die "Error: the maximum number of files are limited to 100";
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

#this subroutine takes a family file and a CNV call file as input to find match between individual identifiers and file names
#in many cases, the individual identifiers are part of the signal file name; for example, individual NA1234 corresponds to a signal file NA1234.lrr_baf.txt
#the subroutine is based on an assumption that the individual identifier is contained within the file name, separated by word boundaries or underlines
#the return vlues %ind_file is a hash, with key as individual identifier, and value as file name
#the return value @nf is an array, with each element being a reference to array of father, mother and offsprings
#the return value @ind is an array containing all individuals (sorted in the same order as in familyfile)
#the return value $non_family_cnv contains all CNV calls that does not belong to any individuals in the given familyfile
sub findIndFileMatch {
	my ($familyfile, $cnvfile, $include_non_family) = @_;
	my (@nf, @ind, %ind_file);
	my $non_family_cnv = '';
	open (FAM, $familyfile) or die "Error: cannot read from family file $familyfile: $!";
	while (<FAM>) {
		s/[\r\n]+$//;
		s/^\s+|\s+$//g;
		my @record = split (/\s+/, $_);
		push @nf, [@record];
		push @ind, @record;
	}
	close (FAM);
	my %ind = map {$_, 1} @ind;
	my $pattern = join ("|", keys %ind);
	
	open (CNV, $cnvfile) or die "Error: cannot read from CNV file $cnvfile: $!";
	while (<CNV>) {
		if (m/state\d(?:,cn=\d)?\s+(\S+)/) {
			my $filename = $1;
			if ($filename =~ m/(\b|_)($pattern)(\b|_)/) {
				if (exists $ind_file{$2} and $ind_file{$2} ne $filename) {
					die "Error: multiple file names ($1 and $ind_file{$3}) in CNV file $cnvfile contains the individual identifier $3";
				}
				$ind_file{$2} ||= $filename;
			} else {
				$non_family_cnv .= $_;
			}
		}
	}
	close (CNV);
	
	#make sure that all individual in family file has CNV calls in the cnv file
	for (@ind) {
		exists $ind_file{$_} or print STDERR "WARNING: cannot find CNV calls for individual $_ in cnvfile $cnvfile\n";
		$ind_file{$_} ||= 'NO_FILE_ASSOCIATED_WITH_INDIVIDUAL';
	}
	return (\%ind_file, \@nf, \@ind, $non_family_cnv);
}

sub convertCNVToHTML {
	my ($cnvfile, $familyfile, $output) = @_;
	my ($temp, $temp1, $temp2, $temp3, $temp4, $temp5, $temp6) = generateRandFile (7);

	my ($ind_file, $nf, $ind, $non_family_cnv) = findIndFileMatch ($familyfile, $cnvfile);
	my %ind_file = %$ind_file;
	my @nf = @$nf;
	my @ind = @$ind;
	my $num_processed_nf = 0;
	my $fh_output;
	
	if (not $nfperfile) {
		open ($fh_output, ">$output") or die "Error: cannot write to HTML file $output: $!\n";
	}
	
	my (@keyword, $keyword);
	if ($highlightfile) {
		@keyword = @{&readHighlightFile ($highlightfile)};
		$keyword = join ('|', @keyword);
	}
	
	for my $nextnf (@nf) {
		my @nextind = @$nextnf;
		my ($result);
		my @table;
		my $childname = 'child';
		
		if ($nfperfile and $num_processed_nf % $nfperfile == 0) {
			open ($fh_output, ">$output.part" . int ($num_processed_nf/$nfperfile+1) . ".html") or die "Error: cannot write to new HTML file: $!\n";
			print STDERR "NOTICE: Writting HTML file $output.part" . int ($num_processed_nf/50+1) . ".html\n";
		}

		@nextind >= 3 or die "Error: invalid record found in family file $familyfile: 3 fields (father, mother, child) expected in <$_>";
		print STDERR "NOTICE: Processing family @nextind\n";
		$result = qx/fgrep -w $ind_file{$nextind[0]} $cnvfile | tee $temp1/;		#father
		push @table, "father($nextind[0])", split (/\n/, $result);
		$result = qx/fgrep -w $ind_file{$nextind[1]} $cnvfile | tee $temp2/;		#mother
		push @table, "mother($nextind[1])", split (/\n/, $result);
		system ("cat $temp1 $temp2 > $temp3");						#father+mother
		
		for my $i (2 .. @nextind-1) {
			@nextind > 3 and $childname = 'child' . ($i-1);
			$childname .= "($nextind[$i])";
			system ("fgrep -w $ind_file{$nextind[$i]} $cnvfile > $temp4");				#child
			
			#print STDERR "next child $ind_file{$nextind[$i]} has", `wc -l $temp4`, "CNVs\n";
			$result = qx/scan_region.pl $temp4 $temp1 -minoverlap $minoverlap --quiet | uniq | tee $temp5/;		#overlap with father
			$result = qx/scan_region.pl $temp4 $temp2 -minoverlap $minoverlap --quiet | uniq | tee $temp6/;		#overlap with mother
			
			
			$result = qx/cat $temp5 $temp6 | sort | uniq -d | tee $temp/;			#overlap with both father and mother
			push @table, "${childname}_overlap_both_parents", split (/\n/, $result);
			
			$result = qx/fgrep -v -f $temp $temp5/;
			push @table, "${childname}_overlap_father", split (/\n/, $result);
			$result = qx/fgrep -v -f $temp $temp6/;
			push @table, "${childname}_overlap_mother", split (/\n/, $result);
			
			$result = qx/fgrep -v -f $temp5 $temp4 | fgrep -v -f $temp6/;
			push @table, "${childname}_denovo", split (/\n/, $result);
		}

		#when --commonfile argument is set, this will calculate rare inherited CNVs
		if ($commonfile) {
			$result = qx/scan_region.pl $temp4 $commonfile -minoverlap $minoverlap/;
			push @table, "${childname}_common", $result;
		}
		
		if ($highlightfile) {
			for my $nextline (@table) {
				$nextline =~ s/\b($keyword)\b/*$1*/g;
			}
		}
		
		$num_processed_nf++;
		writeTableAsHTML ($fh_output, \@table);
	}
	$keeptemp or unlink ($temp, $temp1, $temp2, $temp3, $temp4, $temp5, $temp6);
}

sub readHighlightFile {
	my ($inputfile) = @_;
	my (%keyword, @keyword);
	open (INPUT, $inputfile) or die "Error: cannot read from highlight file $inputfile";
	while (<INPUT>) {
		m/^(\S+)/ and $keyword{$1}++;
	}
	@keyword = keys %keyword;
	return (\@keyword);
}

sub writeTableAsHTML {
	my ($fh_output, $table) = @_;
	print $fh_output "<hr>\n<table>\n";
	my $bgcolor = '#FFFFFF';
	for my $i (0 .. @$table-1) {
		if ($table->[$i] =~ m/\s/) {
			my @td = split (/\s+/, $table->[$i]);
			defined $td[$#td] or pop @td;
			print $fh_output "<tr bgcolor=$bgcolor><td>", join ("</td><td>", @td), "</td></tr>\n";
		} else {
			if ($table->[$i] =~ m/parent/) {
				$bgcolor = '#CCCCFF';
			} elsif ($table->[$i] =~ m/father/) {
				$bgcolor = '#CCFFFF';
			} elsif ($table->[$i] =~ m/mother/) {
				$bgcolor = '#FFCCFF';
			} else {
				$bgcolor = '#FFFFFF';
			}
			print $fh_output "<tr bgcolor=$bgcolor><td><h3>$table->[$i]</h3></td></tr>\n";
		}
	}
	print $fh_output "</table>\n";
}

sub convertAssocToWIG {
	my ($assocfile, $snpposfile) = @_;
	my (%track);
	
	my (%snppos, %snpchr);
	if ($snpposfile) {
		print STDERR "NOTICE: Reading marker position from $snpposfile\n";
		open (PFB, $snpposfile) or die "Error: cannot read from pfb file $snpposfile: $!";
		while (<PFB>) {
			m/^(\S+)\t(\S+)\t(\S+)/ or die "Error: invald record found in position file";
			$snppos{$1} = $3;
			$snpchr{$1} = $2;
		}
	}
	
	if ($assocfile eq 'stdin') {
		*ASSOC = *STDIN;
	} else {
		open (ASSOC, $assocfile) or die "Error: cannot read association file $assocfile: $!";
	}
	print STDERR "NOTICE: Begin reading asso file $assocfile\n";
	my $header = <ASSOC>;
	$header =~ s/^\s*|\s*[\r\n]+$//g;
	my @record = split (/\s+/, $header);
	my ($markerindex, $pindex, $chrindex, $posindex);
	for my $i (0 .. @record-1) {
		if ($record[$i] eq "SNP" or $record[$i] eq "Name") {
			$markerindex = $i;
		} elsif ($record[$i] eq 'P_TDT' or $record[$i] eq 'CHI2_P' or $record[$i] eq 'P') {
			$pindex = $i;
		} elsif ($record[$i] eq 'Chr' or $record[$i] eq 'CHR') {
			$chrindex = $i;
		} elsif ($record[$i] eq 'Position' or $record[$i] eq 'Pos') {
			$posindex = $i;
		}
	}
	defined $markerindex and defined $pindex or die "Error: the first line of assocfile $assocfile does not contain information to determine the position of markers and P-values: <$header>\n";
	%snppos or defined $posindex && defined $chrindex or die "Error: the first line of association file $assocfile does not contains position information (you can use --snpposfile argument to supply a snp position file\n";
	
	my $skipped_marker = 0;
	while (<ASSOC>) {
		s/^\s*|\s*[\r\n]+$//g;				#this is important, as PLINK results may have heading or trailing spaces
		my @record = split (/\s+/, $_);
		my ($marker, $p) = @record[$markerindex, $pindex];
		defined $marker and defined $p or die "Error: unable to find marker and P-value from input line: <$_>\n";
		my $chr = (defined $chrindex)?$record[$chrindex]:$snpchr{$marker};
		defined $chr or next;
		$p =~ m/^[\d\.\-\+e]+$/ or ++$skipped_marker and next;
		$chr =~ m/^\d+$/ or $chr eq 'X' or $chr eq 'Y' or ++$skipped_marker and next;
		$p == 0 and $p = 1e-10;				#so that we can take a logarithm alter;
		$p < 0 and die "Error: P-value cannot be less than ZERO: p=$p in <$_>\n";
		if (%snppos) {					#SNP position information is read from a snppos file (as opposed to the associaion results, which may be in a different genome assembly)
			$snppos{$marker} or ++$skipped_marker and next;
			push @{$track{$snpchr{$marker}}}, [$snppos{$marker}, $p];
		} else {
			push @{$track{$chr}}, [$record[$posindex], $p];
		}
	}
	print STDERR "NOTICE: Finished reading assoc file $assocfile and skipped $skipped_marker markers\n";
	
	$track_name ||= "Genome-wide association -log10(P) in $assocfile";
	print qq{track type=wiggle_0 name="Track: $track_name" description="$track_name" visibility=full color=0,0,0 yLineMark=3 yLineOnOff=on autoScale=on\n};
	for my $chr (keys %track) {
		my @snpstat = sort {$a->[0] <=> $b->[0]} @{$track{$chr}};
		my $prepos = 0;
		$chr eq '23' and $chr = 'X';					#PLINK and some other software use 23 to represent chrX
		$chr eq '24' and $chr = 'Y';
		$chr =~ m/^\d+$/ and $chr >= 25 || $chr == 0 and next;
		$chr =~ m/^(\d+|X|Y|XY)$/ or next;				#skip this marker (probaly in unknown chromosome, or in MT)
		print qq{variableStep chrom=chr$chr span=1\n};
		for my $i (0 .. @snpstat-1) {
			$snpstat[$i]->[0] eq $prepos and next;			#rs4886982 and rs16969329 has identifical position in Broad Affy AGRE data set
			print $snpstat[$i]->[0], "\t", sprintf ("%.3f", -log($snpstat[$i]->[1])/log(10)), "\n";
			$prepos = $snpstat[$i]->[0];
		}
	}
}

sub convertAssocToHTML {
	1;
}


=head1 SYNOPSIS

 visualize_cnv.pl [arguments] <cnv-call-file>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --intype <string>		input file type: cnv (default) or assoc
 	    --output <file>		output root file name (default=inputfile name)
 	    --format <string>		out file format (bed, html, beadstudio, wig)
 	    --track_name <string>	name for UCSC genome browser track (default=inputfile name)
 	    --bedstrand <string>	strand orientation for UCSC genome browser track (+ or - or .)
 	    --familyfile <file>		a file containg identifiers of family members
 	    --minoverlap <float>	minimum overlap to determine concordance
 	    --commonfile <file>		a file containg common CNV regions
 	    --highlightfile <file>	a file containing keywords to be highlighted in HTML output
 	    --snpposfile <file>		a file whose first tab-delimited columns are marker, chr and pos
 	    --nfperfile <int>		number of nuclear families per file (default=50)
 	    --idmapfile <file>		a file continaing file name (in PennCNV call) and sample id (in BeadStudio) mapping
 	    --flankinglength <int>	minimum flanking distance in plot around the CNV (default: 50000)
 	    --pdfout			produce vector-based PDF file in plot (default: JPG)

 Function: reformat CNV calls by PennCNV for visualization in UCSC Genome 
 Browser, or in web browser, or in BeadStudio software

 Example: visualize_cnv.pl -format bed -out test.bed sampleall.cnv
          visualize_cnv.pl -format html -fam familyfile -out test.html sampleall.cnv
          visualize_cnv.pl -format beadstudio -idmap id_map.txt -out test.xml sampleall.cnv
          visualize_cnv.pl -intype assoc -format wig -out result.wig plink.tdt
          visualize_cnv.pl -format tab -out test.txt sampleall.cnv
          visualize_cnv.pl -format plot -signal father.txt sampleall.cnv

 Version: $LastChangedDate: 2012-10-23 23:32:05 -0700 (Tue, 23 Oct 2012) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--intype>

input file type can be cnv (default) or assoc, which is a file containing 
genome-wide association test results, with one marker per line and each line 
contains at least the marker name and P-value.

=item B<--format>

output format can be bed (the format for UCSC genome browser custom tracks), 
html (the format that can be displayed by a web browser), beadstudio (the format 
that can be loaded by beadstudio software to visualize cnv calls along the 
genome), wig (the format for UCSC genome browser custom tracks that displays 
quantitative values for individual bases/SNPs)

=item B<-track_name>

specify the name for UCSC genome browser custom tracks

=item B<--bedstrand>

the strand annotation for UCSC genome browser custom tracks. By default no 
directionality is attatched to CNV tracks.

=item B<--familyfile>

a file containng identifiers of samples within nuclear families, with one family 
per line and each line contains identifiers for father, mother, and one or more 
offsprings. This file is useful in generating PennCNV report for families in 
HTML format.

=item B<--minoverlap>

minimum overlap of CNV calls to determine the concordance of two CNVs 
(default=any overlap). When specifying a high number (such as 0.5), this 
requires that half of one CNV must overlap with another CNV to establish 
concordance.

=item B<--commonfile>

a file containing common CNV regions, which might be useful in generating HTML 
files for rare CNVs only. each line has the format like chr2:1000- 5000.

=item B<--highlightfile>

a file containing gene identifiers that should be highlighted in the HTML file. 
The gene list could be for example candidate genes for a disease. This helps 
spotting important CNVs in the CNV calls.

=item B<--snpposfile>

a file containing SNP genome coordinates, which is used to generate wig 
formatted file in UCSC genome browser.

=item B<--nfperfile>

the output HTML file are divided such that each file contains only a limited set 
of nuclear families to reduce the burden of web browsers. Default=50 nuclear 
families per HTML file.

=item B<--idmapfile>

a file continaing file name (in PennCNV call) and sample id (in BeadStudio) 
mapping. This file is useful when converting PennCNV calls to Illumina 
BeadStudio/GenomeStudio bookmarks.

=item B<--flankinglength>

when -format plot is used, this argument specifies the the minimum flanking 
distance in plot around the CNV (default: 50000) when the - -format plot is 
used. If users want to see whole- chromosome markers, use a very large value 
(for example, 200m) in this argument. The CNV region will be marked as red dots, 
so it can be still differentiated from the rest of the markers in blue dots.

=item B<--pdfout>

when -format plot is used, generate PDF output files for publication purposes. 
By default, low quality JPG files will be produced.


=back

=head1 DESCRIPTION

This program is designed for converting CNV calls generated by the PennCNV 
software to various other formats for better visualization in UCSC genome 
browser, or a web browser, or the BeadStudio software. An extra feature is to 
convert SNP-based genome-wide association results to UCSC genome browser wiggle 
format, so that one can visualize the association results on a genomic region 
easily.

The most common use of this program would be to convert CNV calls generated by 
PennCNV to BED files to be displayed in UCSC genome browser. The first a few 
lines of the PennCNV output file may look like this:

	chr1:59077355-59078584        numsnp=3      length=1,230       state5,cn=3 sample1.txt startsnp=rs942123 endsnp=rs3015321 father triostate=533
	chr1:147305744-147427061      numsnp=7      length=121,318     state5,cn=3 sample1.txt startsnp=rs11579261 endsnp=rs3853524 father triostate=535
	chr1:147305744-147427061      numsnp=7      length=121,318     state5,cn=3 sample3.txt startsnp=rs11579261 endsnp=rs3853524 offspring triostate=535

After running the command with all default options:

	visualize_cnv.pl sampleall.allcnv

The first a few lines of output look like this:

	track name="Track: CNVs in sampleall.hg18.allcnv" description="CNVs in sampleall.hg18.allcnv" visibility=2 itemRgb="On"
	chr1    59077354        59078584        sample1.txt     30      .       0       0       0,255,0
	chr1    147305743       147427061       sample1.txt     70      .       0       0       0,255,0
	chr1    147305743       147427061       sample3.txt     70      .       0       0       0,255,0

This file is in BED format, and it can be loaded into UCSC genome browser for 
visualization within the context of genome annotations.

This program can also convert CNV calls to the BookMark files to be loaded into 
BeadStudio for visualization as well. for example, after running the command:

	visualize_cnv.pl sampleall.hg18.allcnv -format beadstudio -id sampleidfile > cnvcall.xml

An XML file will be generated that can be loaded into the BeadStudio Genome 
Viewer (click View menu, select Bookmark Viewer, select Load Bookmark Analysis 
File).

Another very useful utility of this program is to generate signal intensity plots for all CNV calls for a given sample. This can be achieved by commands such as this one: 

	 visualize_cnv.pl -format plot -signal father.txt sampleall.cnv

This command essentially check the sampleall.cnv file, find all CNV calls for 
the father.txt file, then examine the father.txt file for signal intensity 
(LRR/BAF) values for all markers with each CNV, then plot them as JPG files. The 
plotting function requires R to be installed in the system. Several output files 
with suffix JPG will be generated, each corresponding to one CNV call. Users can 
open these files using any graphics viewer software and confirm whether the CNV 
looks real, or may represent some technical artifacts.

For questions, comments or bug reports, please contact me at 
kai@openbioinformatics.org.

=cut
