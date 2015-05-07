#!/usr/bin/perl -w
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;
BEGIN {($_=$0)=~s{[^\\\/]+$}{};$_||="."}
use lib $_, $_."kext";
use khmm;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2010-05-01 16:33:32 -0700 (Sat, 01 May 2010) $';

our ($verbose, $help, $man);
our ($startsnp, $endsnp, $denovocn, $allcn, $pfbfile, $directory, $outfile, $logfile, $hmmfile, $zdiff);
our (@nffile, @nfsig, $hmm);


GetOptions ('verbose'=>\$verbose, 'help'=>\$help, 'man|m'=>\$man, 'startsnp=s'=>\$startsnp, 'endsnp=s'=>\$endsnp, 'denovocn=i'=>\$denovocn, 'allcn=i'=>\$allcn, 'pfbfile=s'=>\$pfbfile,
	'directory=s'=>\$directory, 'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile, 'hmmfile=s'=>\$hmmfile, 'zdiff=f'=>\$zdiff) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV >= 1 or pod2usage ("Syntax error");

defined $pfbfile or pod2usage ("Error in argument: please specify the --pfbfile argument");
defined $denovocn or defined $allcn or pod2usage ("Error in argument: please specify either the --denovocn argument (for de novo CNV) or the --allcn argument");
defined $denovocn and defined $allcn and pod2usage ("Error in argument: pelase specify either the --denovocn argument (for de novo CNV) or the --allcn argument but not both");
$allcn and length ($allcn) == @ARGV || pod2usage ("Error in argument: please specify the --allcn with ${\(scalar @ARGV)} digits, which corresponds to ${\(scalar @ARGV)} signal files");
$allcn and $allcn =~ m/^\d+$/ || pod2usage ("Error in argument: the --allcn argument must be composed of digits, corresponding to copy number in ${\(scalar @ARGV)} signal files");
defined $denovocn and $denovocn =~ m/^[13]$/ || pod2usage ("Error in argument: the --denovocn argument must be one single digit (current version only supports 1 or 3)");
defined $denovocn and @ARGV==3 || pod2usage ("Error in argument: when --denovocn is set, only three signal files (father, mother, offspring) can be supplied in command line");
defined $startsnp and defined $endsnp or pod2usage ("Error in argument: please specify the --startsnp and the --endsnp argument");
defined $hmmfile or pod2usage ("Error in argument: please supply the --hmmfile argument");
$hmm = readHMMFile ($hmmfile);
$zdiff ||= 0.5;

(@nffile) = @ARGV;
my ($snp_chr, $snp_pos) = readPFB ($pfbfile);
my ($tchr, $tstart, $tend);		#target chromoosme, target start bp, target end bp

$tchr = $snp_chr->{$startsnp};
defined $tchr or confess "Error: the --startsnp $startsnp is not found in PFB file $pfbfile\n";
defined $snp_chr->{$endsnp} or confess "Error: the --endsnp $endsnp is not found in PFB file $pfbfile\n";
$tchr eq $snp_chr->{$endsnp} or confess "Error: the --startsnp is located in chr $tchr but --endsnp is located in chr $snp_chr->{$endsnp}\n";

$tstart = $snp_pos->{$startsnp};
$tend = $snp_pos->{$endsnp};

$outfile and open (STDOUT, ">$outfile") || confess "Error: cannot write to output file $outfile: $!\n";
$logfile and open (STDERR, ">$logfile") || confess "Error: cannot write to log file $logfile: $!\n";

for my $i (0 .. @nffile-1) {
	push @nfsig, readSignalFile ($nffile[$i], $tchr, $tstart, $tend, $snp_chr, $snp_pos, $directory);
}

if ($denovocn) {
	analyzeDeNovoCNV (\@nfsig, $denovocn);
} elsif ($allcn) {
	analyzeCNV (\@nfsig, $allcn);
}

sub analyzeCNV {
	my ($nfsig, $allcn) = @_;
	my @allcn = split(//, $allcn);
	
	for my $i (0 .. @{$nfsig->[0]}-1) {
		my (@lrr, @baf, @cngt);
		
		for my $j (0 .. @$nfsig-1) {
			push @lrr, $nfsig->[$j]->[$i][2];
			push @baf, $nfsig[$j]->[$i][3];
			push @cngt, bafToGeno ($nfsig->[$j][$i][3], $allcn[$j], $hmm);
		}
		
		if ($i == 0) {
			print "Name";
			for (1 .. @$nfsig) {
				print "\tLRR($_)";
			}
			for (1 .. @$nfsig) {
				print "\tBAF($_)";
			}
			for (1 .. @$nfsig) {
				print "\tGENO($_)";
			}
			print "\n";
		}
		
		print $nfsig->[0][$i][1];
		print "\t", join ("\t", @lrr);				#LRR values
		print "\t", join ("\t", @baf);				#BAF values
		print "\t", join ("\t", @cngt);				#GT values
		print "\n";
	}
}
	

sub analyzeDeNovoCNV {
	my ($triosig, $denovocn) = @_;
	my %origin = (	AA_AB_B=>'F',
			AB_AA_B=>'M',
			AB_BB_A=>'M',
			BB_AB_A=>'F',
			AA_BB_B=>'F',
			AA_BB_A=>'M',
			BB_AA_A=>'F',
			BB_AA_B=>'M',
			AA_AB_AAB=>'F',
			AA_AB_ABB=>'M',
			AB_AA_AAB=>'M',
			AB_AA_ABB=>'F',
			BB_AB_ABB=>'F',
			BB_AB_AAB=>'M',
			AB_BB_ABB=>'M',
			AB_BB_AAB=>'F'
			);
	
	my ($sigf, $sigm, $sigo) = @$triosig;
	my ($cf, $cm) = qw/0 0/;
	
	if (@$sigf != @$sigm) {
		confess "Error: marker counts discordance found: father=${\(scalar @$sigf)} mother=${\(scalar @$sigm)}\n";
	}
	if (@$sigf != @$sigo) {
		confess "Error: marker counts discordance found: father=${\(scalar @$sigf)} offspring=${\(scalar @$sigo)}\n";
	}
	
	print STDERR "NOTICE: Analyzing trio @nffile\n";
	print "Name\tLRR_F\tLRR_M\tLRR_O\tBAF_F\tBAF_M\tBAF_O\tGENO_F\tGENO_M\tGENO_O\tOrigin\n";
	for my $i (0 .. @$sigf-1) {
		print $sigf->[$i][1];
		print "\t", join ("\t", $sigf->[$i][2], $sigm->[$i][2], $sigo->[$i][2]);		#LRR values
		print "\t", join ("\t", $sigf->[$i][3], $sigm->[$i][3], $sigo->[$i][3]);		#BAF values
		my ($gtf, $gtm, $gto) = (bafToGeno ($sigf->[$i][3], 2, $hmm), bafToGeno ($sigm->[$i][3], 2, $hmm), bafToGeno ($sigo->[$i][3], $denovocn, $hmm));
		defined $gto or confess;
		print "\t", join ("\t", $gtf, $gtm, $gto);
		my ($origin) = $origin{join("_", $gtf, $gtm, $gto)};
		if ($origin) {
			if ($origin eq 'F') {
				$cf++;
			} else {
				$cm++;
			}
			print "\t", $origin, "\n";
		} else {
			print "\t?\n";
		}
	}
	
	if ($cf+$cm) {
		my $p = khmm::bitest ($cf+$cm, $cf, 0.5) * 2;
		print STDERR "NOTICE: Evidence for parental origin for the putative de novo CNVs (de novo CN=$denovocn in trio @nffile ): Marker= ${\(scalar @$sigf)} Paternal_origin(F)= $cf Maternal_origin(M)= $cm P-value= $p\n";
	} else {
		print STDERR "NOTICE: No evidence of parental origin for the putative de novo CNVs (de novo CN=$denovocn in trio @nffile )\n";
	}
}


sub bafToGeno {
	my ($baf, $cn, $hmm) = @_;
	my $gt;
	if (not $baf =~ m/^(([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?)$/) {		#BAF is not a number
		return 'NC';
	}
	if ($baf > 1 or $baf < 0) {
		return 'NC';				#probably a non-polymorphic marker, or an error in data entry
	}
	
	my @thres = @{$hmm->{B2_sd}};			#SD for 0, 0.25, 0.33, 0.5
	my ($z1, $z2, $z3, $z4, $z5);
	
	if ($cn == 4) {
		$z1 = abs ($baf-1)/$thres[0];
		$z2 = abs ($baf-0.75)/$thres[1];
		$z3 = abs ($baf-0.5)/$thres[3];
		$z4 = abs ($baf-0.25)/$thres[1];
		$z5 = abs ($baf-0)/$thres[0];
		
		if ($baf <= 1 and $baf >0.75) {
			if ($z1-$z2 > $zdiff) {
				$gt = "ABBB";
			} elsif ($z2-$z1 > $zdiff) {
				$gt = "BBBB";
			} else {
				$gt = "NC";
			}
		} elsif ($baf <= 0.75 and $baf > 0.5) {
			if ($z2-$z3 > $zdiff) {
				$gt = "AABB";
			} elsif ($z3-$z2 > $zdiff) {
				$gt = "ABBB";
			} else {
				$gt = "NC";
			}
		} elsif ($baf <= 0.5 and $baf > 0.25) {
			if ($z3-$z4 > $zdiff) {
				$gt = "AAAB";
			} elsif ($z4-$z3 > $zdiff) {
				$gt = "AABB";
			} else {
				$gt = "NC";
			}
		} elsif ($baf <= 0.25 and $baf >= 0) {
			if ($z4-$z5 > $zdiff) {
				$gt = "AAAA";
			} elsif ($z5-$z4 > $zdiff) {
				$gt = "AAAB";
			} else {
				$gt = "NC";
			}
		}
	} elsif ($cn == 3) {
		$z1 = abs ($baf-1)/$thres[0];
		$z2 = abs ($baf-0.667)/$thres[2];
		$z3 = abs ($baf-0.333)/$thres[2];
		$z4 = abs ($baf-0)/$thres[0];
		
		if ($baf <= 1 and $baf > 0.667) {
			if ($z1-$z2 > $zdiff) {
				$gt = "ABB";
			} elsif ($z2-$z1 > $zdiff) {
				$gt = "BBB";
			} else {
				$gt = "NC";
			}
		} elsif ($baf <= 0.667 and $baf > 0.333) {
			if ($z2-$z3 > $zdiff) {
				$gt = "AAB";
			} elsif ($z3-$z2 > $zdiff) {
				$gt = "ABB";
			} else {
				$gt = "NC";
			}
		} elsif ($baf <= 0.333 and $baf >= 0) {
			if ($z3-$z4 > $zdiff) {
				$gt = "AAA";
			} elsif ($z4-$z3 > $zdiff) {
				$gt = "AAB";
			} else {
				$gt = "NC";
			}
		}
	} elsif ($cn == 2) {
		$z1 = abs ($baf-1)/$thres[0];
		$z2 = abs ($baf-0.5)/$thres[3];
		$z3 = abs ($baf-0)/$thres[0];
		
		if ($baf <= 1 and $baf > 0.5) {
			if ($z1-$z2 > $zdiff) {
				$gt = "AB";
			} elsif ($z2-$z1 > $zdiff) {
				$gt = "BB";
			} else {
				$gt = "NC";
			}
		} elsif ($baf <= 0.5 and $baf >= 0) {
			if ($z2-$z3 > $zdiff) {
				$gt = "AA";
			} elsif ($z3-$z2 > $zdiff) {
				$gt = "AB";
			} else {
				$gt = "NC";
			}
		}
	} elsif ($cn == 1) {			#this is similar to the cn=2 above, 
		$z1 = abs ($baf-1)/$thres[0];
		$z2 = abs ($baf-0.5)/$thres[3];
		$z3 = abs ($baf-0)/$thres[0];
		
		if ($baf <= 1 and $baf > 0.5) {
			if ($z2-$z1 > $zdiff) {
				$gt = "B";
			} else {
				$gt = "NC";
			}
		} elsif ($baf <= 0.5 and $baf >= 0) {
			if ($z2-$z3 > $zdiff) {
				$gt = "A";
			} else {
				$gt = "NC";
			}
		}
	} elsif ($cn == 0) {
		$gt = 'NC';
	}
	return $gt;
}

sub readPFB {
	my ($pfbfile) = @_;
	my (%snp_chr, %snp_pos, %snp_pfb);
	my (%skip_snp, $skip_snp_count, $count_np);
	open (PFB, $pfbfile) or confess "Error: cannot read from pfb file $pfbfile: $!\n";
	print STDERR "NOTICE: Reading marker coordinates and population frequency of B allele (PFB) from $pfbfile ...";	
	while (<PFB>) {
		s/[\r\n]+$//;							#delete line feed and return characters
		m/^(\S+)\t(\S+)\t(\S+)\t(\S+)$/ or confess "Error: invalid record found in PFB file $pfbfile (4 tab-delimited records expected): <$_>\n";
		my ($name, $chr, $pos, $pfb) = ($1, $2, $3, $4);
		$chr eq 'Chr' and next;						#this is the header line in a regular PFB file
		if (not $chr =~ m/^(\d+|X|Y)$/ or $chr eq '0') {
			$skip_snp{$chr}++;
			$skip_snp_count++;
			next;
		}
		$snp_chr{$name} = $chr;
		$snp_pos{$name} = $pos;
		$snp_pfb{$name} = $pfb;
		if ($snp_pfb{$name} > 1) {
			$count_np++;
		}
		$snp_pfb{$name} >= 0 and $snp_pfb{$name} < 0.01 and $snp_pfb{$name} = 0.01;	#prevent underflow/overflow errors in subsequent calculations.
		$snp_pfb{$name} <= 1 and $snp_pfb{$name} > 0.99 and $snp_pfb{$name} = 0.99;	#prevent underflow/overflow errors in subsequent calculations.
	}
	print STDERR " Done with ${\(scalar keys %snp_chr)} records";
	if ($skip_snp_count) {
		print STDERR " ($skip_snp_count records in chr ", join (",", keys %skip_snp), " were discarded)\n";
	} else {
		print STDERR "\n";
	}
	$verbose and print STDERR "NOTICE: FPB file $pfbfile contains $count_np non-polymorphic markers\n";
	close (PFB);
	return (\%snp_chr, \%snp_pos, \%snp_pfb);
}

sub readSignalFile {
	my ($signalfile, $tchr, $tstart, $tend, $snp_chr, $snp_pos, $directory) = @_;
	my ($header, $header_seg, $name_index, $lrr_index, $baf_index, $pos_index, $chr_index);
	
	$directory and $signalfile = "$directory/$signalfile";
	
	open (SIG, $signalfile) or confess "Error: cannot read from file $signalfile: $!";
	$header = <SIG>;
	$header =~ m/(.*)Name/ or confess "Error: the header file of signalfile $signalfile does not contain log r ratio annotation";
	$header_seg = $1;
	$name_index = ($header_seg =~ tr/\t/\t/+0);
	$header =~ m/(.*)Log R Ratio/ or confess "Error: the header file of signalfile $signalfile does not contain log r ratio annotation";
	$header_seg = $1;
	$lrr_index = ($header_seg =~ tr/\t/\t/+0);
	$header =~ m/(.*)B Allele Freq/ or confess "Error: the header file of signalfile $signalfile does not contain B Allele Freq annotation";
	$header_seg = $1;
	$baf_index = ($header_seg =~ tr/\t/\t/+0);
	
	
	my (@data, $count);
	my $skipped = 0;
	while (<SIG>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		my ($curlrr, $curbaf, $curname) = @record[$lrr_index, $baf_index, $name_index];
		my ($curchr, $curpos) = ($snp_chr->{$curname}, $snp_pos->{$curname});
		
		defined $curchr or ++$skipped and next;
		
		$curchr eq $tchr or next;
		$curpos >= $tstart or next;
		$curpos <= $tend or next;
		
		push @data, [$curpos, $curname, $curlrr, $curbaf];
		$count++;
	}
	@data = sort {$a->[0]<=>$b->[0]} @data;
	close (SIG);
	print STDERR "NOTICE: For the region chr$tchr:$tstart-$tend, $count markers were identified from $signalfile\n";
	return (\@data);
}

sub readHMMFile {
	my ($inputfile) = @_;
	my (%hmm, @cell);
	open (HMM, $inputfile) or confess "\nERROR: cannot read from HMM file $hmmfile: $!\n";
	my @line = <HMM>;
	map {s/[\r\n]+$//} @line;
	$line[0] eq 'M=6' or confess "\nERROR: invalid record found in HMM file: <$_> ('M=6' expected)\n";
	$line[1] eq 'N=6' or confess "\nERROR: invalid record found in HMM file: <$_> ('N=6' expected)\n";
	$line[2] eq 'A:' or confess "\nERROR: invalid record found in HMM file: <$_> ('A:' expected)\n";
	$line[9] eq 'B:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B:' expected)\n";
	$line[16] eq 'pi:' or confess "\nERROR: invalid record found in HMM file: <$_> ('pi:' expected)\n";
	$line[18] eq 'B1_mean:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B1_mean:' expected)\n";
	$line[20] eq 'B1_sd:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B1_sd:' expected)\n";
	$line[22] eq 'B1_uf:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B1_uf:' expected)\n";
	$line[24] eq 'B2_mean:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B2_mean:' expected)\n";
	$line[26] eq 'B2_sd:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B2_sd:' expected)\n";
	$line[28] eq 'B2_uf:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B2_uf:' expected)\n";
	if (@line > 30) {
		#print STDERR "NOTICE: HMM model file $inputfile contains parameters for non-polymorphic (NP) probes\n";
		@line == 36 or @line == 38 or confess "\nERROR: invalid number of records found in HMM file: 30 or 36 or 38 lines expected\n";
		$line[30] eq 'B3_mean:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B3_mean:' expected)\n";
		$line[32] eq 'B3_sd:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B3_sd:' expected)\n";
		$line[34] eq 'B3_uf:' or confess "\nERROR: invalid record found in HMM file: <$_> ('B3_uf:' expected)\n";
		if (@line == 38) {
			$line[36] eq 'DIST:' or confess "\nERROR: invalid record found in HMM file: <$_> ('DIST:' expected)\n";
		}
	}

	for my $i (3 .. 8) {
		@cell = split (/\s+/, $line[$i]);
		abs (sum (\@cell) - 1) < 1e-5 or confess "\nERROR: invalid line ${\($i+1)} in HMM file: <$_> (sum of line should be 1)\n";
		push @{$hmm{'A'}}, [@cell];
	}
	
	@cell = split (/\s+/, $line[17]);
	abs (sum (\@cell) - 1) < 1e-5 or confess "\nERROR: invalid line in HMM file: <$line[17]> (sum of line should be 1)\n";
	push @{$hmm{'pi'}}, @cell;
	
	@cell = split (/\s+/, $line[19]);
	@cell == 6 or confess "\nERROR: invalid line in HMM file: <@cell> (6 fields expected)\n";
	push @{$hmm{'B1_mean'}}, @cell;
	
	@cell = split (/\s+/, $line[21]);
	@cell == 6 or confess "\nERROR: invalid line in HMM file: <@cell> (6 fields expected)\n";
	push @{$hmm{'B1_sd'}}, @cell;
	grep {$_>0} @cell == 6 or confess "\nERROR: invalid line in HMM file: <@cell> (all values should be between greater than zero)\n";
	
	@cell = split (/\s+/, $line[23]);
	@cell == 1 or confess "\nERROR: invalid line in HMM file: <@cell> (1 fields expected)\n";
	push @{$hmm{'B1_uf'}}, @cell;
	
	@cell = split (/\s+/, $line[25]);
	@cell == 5 or confess "\nERROR: invalid line in HMM file: <@cell> (5 fields expected)\n";
	push @{$hmm{'B2_mean'}}, @cell;
	
	@cell = split (/\s+/, $line[27]);
	@cell == 5 or confess "\nERROR: invalid line in HMM file: <@cell> (5 fields expected)\n";
	push @{$hmm{'B2_sd'}}, @cell;
	grep {$_>0} @cell == 5 or confess "\nERROR: invalid line in HMM file: <@cell> (all values should be between greater than zero)\n";

	@cell = split (/\s+/, $line[29]);
	@cell == 1 or confess "\nERROR: invalid line in HMM file: <@cell> (1 fields expected)\n";
	push @{$hmm{'B2_uf'}}, @cell;
	
	if (@line > 30) {
		@cell = split (/\s+/, $line[31]);
		@cell == 6 or confess "\nERROR: invalid line in HMM file: <@cell> (6 fields exepcted)\n";
		push @{$hmm{'B3_mean'}}, @cell;

		@cell = split (/\s+/, $line[33]);
		@cell == 6 or confess "\nERROR: invalid line in HMM file: <@cell> (6 fields expected)\n";
		push @{$hmm{'B3_sd'}}, @cell;
		grep {$_>0} @cell == 6 or confess "\nERROR: invalid line in HMM file: <@cell> (all values should be between greater than zero)\n";

		@cell = split (/\s+/, $line[35]);
		@cell == 1 or confess "\nERROR: invalid line in HMM file: <@cell> (1 fields expected)\n";
		push @{$hmm{'B3_uf'}}, @cell;
	}
	
	if (@line == 38) {
		$line[37] =~ m/^\d+$/ or confess "Error: invalid line in HMM file: <$line[37]> (an integer expected for DIST)\n";
		$hmm{dist} = $line[37];
	}
	
	close (HMM);
	return (\%hmm);
}		

sub sum {
	my ($score) = @_;
	@$score or confess "\nERROR: NO VALUES for calculating sum\n";
	my $sum;
	for (@$score) {
		$sum += $_;
	}
	return $sum;
}

=head1 SYNOPSIS

 infer_snp_allele.pl [arguments] <fatherfile> <motherfile> <offspringfile | ...>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --startsnp <string>		start SNP of the CNV
 	    --endsnp <string>		end SNP of the CNV
 	    --denovocn <int>		copy number of offspring (for de novo CNV)
 	    --allcn <int>		copy number of nuclear family (for inherited CNV)
 	    --pfbfile <file>		population frequency of B allele file
 	    --directory <dir>		the directory where the signal files are stored
 	    --outfile <file>		the output file (default: STDOUT)
 	    --logfile <file>		the log file (default: STDERR)
 	    --hmmfile <file>		the HMM file used in PennCNV calling
 	    --zdiff <float>		a parameter controlling genotype call rate (default:0.5)
 	    
 Function: infer SNP genotypes in CNV regions, or assign P-values to putative de novo CNV 
 calls

 Example: #assign P-value to de novo CNV calls:
          infer_snp_allele.pl -pfbfile hh550.pfb -hmmfile hh550.hmm -denovocn 3 -start rs100 -end rs2000 father.txt mother.txt offspring.txt
          
          #infer SNP allele composition within CNV regions for one or a few subjects
          infer_snp_allele.pl -pfbfile hh550.pfb -hmmfile hh550.hmm -allcn 32333 -start rs300 -end rs500 sample1.txt sample2.txt sample3.txt sample4.txt sample5.txt

 Version: $LastChangedDate: 2010-05-01 16:33:32 -0700 (Sat, 01 May 2010) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--startsnp>

specify the start SNP in the CNV region

=item B<--endsnp>

specify the end SNP in the CNV region

=item B<--denovocn>

specify the copy number of the de novo CNV call (the parents are assumed to have copy number of 2)

=item B<--allcn>

specify the copy number of a nuclear family (father, mother and offspring). For 
example, a --allcn of 323 means that father has 3 copies, mother has two copies 
and child has 3 copies

=item B<--pfbfile>

specify the PFB file, which is used to infer all the SNPs between the start SNP 
and end SNP.

=item B<--hmmfile>

specify the HMM file, which is used to calculate the Z-scores for the genotype 
calls (Z-score measures the distance between an observed BAF value versus 
expected cluster centers, so it controls the genotype confidence scores).

=item B<--directory>

specify the directory that stores the signal files.

=item B<--outfile>

specify the output file that contains all the inferred SNP genotypes. By default 
this information is printed to STDOUT.

=item B<--logfile>

specify the log file that contains all the counts of markers and their 
associated P-values. By default this information is printed to STDERR.

=item B<--zdiff>

a parameter controlling the genotype call rate. The parameter refers to the 
difference of the two Z-scores for two competing genotype calls. Smaller Z-score 
will usually get a genotype call, but if two Z-scores are close together 
(absolute value less than the --zdiff), then a no call (NC) will be assigned as 
the genotype call.

=back

=head1 DESCRIPTION

This program is part of the PennCNV package, and it is used to validate the de 
novo CNV calls based on genotype information, and assigns a P-value to the 
validated de novo CNVs. Alternatively, the program can be used to generate CNV-
based genotype calls, given the signal file and given the known copy number.

=head2 validate de novo CNVs and assign P-value

This program should be used after calling CNVs and after generating a putative 
de novo CNVs (with known startsnp and endsnp). It provides an additional level 
of validation to de novo CNV calls, and filter out spurious de novo CNV calls.

For example, suppose the following CNV call is detected from offspring, but not in father or mother:

	chr6:79029920-79088461        numsnp=25     length=58,542      state2,cn=1 offspring.txt startsnp=rs100 endsnp=rs200 offspring triostate=332

To further validate it by genotype information, we can use this program:

	infer_snp_allele.pl -pfbfile hh500.pfb -denovocn 2 -hmm hh550.hmm -start rs100 -end rs2000 father.txt mother.txt offspring.txt

Some results will be printed: All SNPs genotypes and signal intensities within 
the CNV regions will be printed. For father and mother, this is the regular 
genotype with two allels, but for offspring, this will be single-allele genoytpe 
due to the copy number of 1.

In addition, the program will print how many informative markers are in the 
region, how many support paternal origin and how many support maternal origin, 
then give a P-value to this analysis.

For deletion de novo CNVs, this genotype inference is usually quite accurate: 
for example, I have cases where hundreds of markers support paternal origin and 
ZERO markers support maternal origin for a de novo deletion. For duplication de 
novo CNVs, the data is somewhat more heterogeneous, due to the inherent 
difficulty in inferring tri-allelic genotypes. For example, if 20 SNPs support 
paternal origin, maybe several will support maternal as well. But in any case, 
the P-value can be a good indicator on whether the de novo CNV call is genuine 
or not.

=head2 infer SNP allele composition within CNV region

If the --allcn argument is supplied to the program, then the program will try to 
generate CNV-based genotype calls, that is, genoype calls that are not composed 
of two alleles at CNV regions. For a deletion region, example of calls are A, B. 
For a duplication region, example of calls are AAB, ABB and AABB. The number of 
digits supplied by the --allcn argument must correspond to the number of signal 
files supplied in the command line. For example, suppose we want to infer the 
CNV-based genotype calls for four people:

	infer_snp_allele.pl -pfbfile hh500.pfb -allcn 2121 -hmm hh550.hmm -start rs100 -end rs2000 father.txt mother.txt offspring1.txt offspring2.txt

This means that among the four people, the mother and offspring2 has CN=1 in the 
specified region, whereas all other people has normal copy in the specified 
region.

For questions, comments and bug reports, contact me at kai@openbioinformatics.org.

=cut