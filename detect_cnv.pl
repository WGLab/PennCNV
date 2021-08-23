#!/usr/bin/env perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;
BEGIN {($_=$0)=~s{[^\\\/]+$}{};$_||="./"}
use lib $_, $_."kext";
eval {
	require "khmm.pm";
};

if ($@ and $@ =~ m/^Can't locate loadable object/ || $@ =~ m/^Floating point exception/ || $@ =~ m/^Core dumped/) {	#'
	require Config;
	my $arch = $Config::Config{archname} || 'UNKNOWN';
	print STDERR "PennCNV compilation error: $@\n";
	print STDERR "PennCNV compilation error: Your system architecture is '$arch', which is not compatible with pre-compiled executables in PennCNV package.\n";
	print STDERR "PennCNV compilation error: Please download source code from http://www.openbioinformatics.org/penncnv and compile executable program.\n";
	exit (100);
}
if ($@ and $@ =~ m/^Can't locate khmm.pm/) {	#'
	print STDERR "PennCNV module error: Can't locate khmm.pm in $_/kext\n";
	print STDERR "PennCNV module error: Please make sure that $_/kext/khmm.pm file exists after decompressing the PennCNV package\n";
	exit (200);
}
if ($@) {
	print STDERR "PennCNV error: $@\n";
	print STDERR "PennCNV compilation error: Please download source code from http://www.openbioinformatics.org/penncnv/ and compile executable program.\n";
	exit (400);
}

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2013-02-08 11:10:50 -0800 (Fri, 08 Feb 2013) $';

our ($verbose, $help, $man);
our ($inputfile, @inputfile);
our ($train, $ctrain, $test, $hmmfile, $region, $loh, $minsnp, $minlength, $minconf, $pfbfile, $trio, $quartet, $chrx, $chry, $medianadjust, $cnvfile, $heterosomic_threshold, $denovo_rate, $sdadjust, $bafadjust, $sexfile, $fmprior, $listfile, $exclude_heterosomic, $output, $joint, $record_count, $logfile, $summary, $confidence, $tabout, $coordinate_from_input, $gcmodelfile, $cctest, $phenofile, $onesided, $control_label, $type_filter, $directory, $flush, $validate, $delfreq, $dupfreq, $backfreq, $startsnp, $endsnp, $candlist, $bafxhet_threshold, $valilog, $cnprior, $tseq, $gamma, $paramfile, $lastchr, $cache, $stroma, $refchr, $refgcfile, @ref_median);
our (@fmprior, @cnprior, $gamma_k, $gamma_theta);


processArgument ();

if ($train) {
	trainCHMM (\@inputfile, $hmmfile, $pfbfile, $gcmodelfile, $directory);
} elsif ($ctrain) {
	ctrainCHMM ($inputfile[0], $hmmfile);
} elsif ($test) {
	my $hmm = readHMMFile ($hmmfile);
	$loh and $hmm->{B1_mean}[3] > 1 and pod2usage ("Error: the supplied -hmmfile $hmmfile is not suitable for LOH inference, but for CNV detection only");
	if ($stroma) {
		newtumorCHMM (\@inputfile, $hmmfile, $pfbfile, $sexfile, $gcmodelfile, $directory);
	} else {
		newtestCHMM (\@inputfile, $hmmfile, $pfbfile, $sexfile, $gcmodelfile, $directory);
	}
} elsif ($validate) {
	my $candregion;
	if ($candlist) {
		$candregion = retrieveCandidateRegion ($candlist);
	} else {
		if ($cnprior) {
			$candregion = [[$startsnp, $endsnp, @cnprior]];
			$verbose and print STDERR "NOTICE: For validating CNV calls of $startsnp-$endsnp, the prior distribution is given as @cnprior\n";
		} else {
			$candregion = [[$startsnp, $endsnp, convertDelDupFreqToPrior ($delfreq, $dupfreq, $backfreq)]];
			$verbose and print STDERR "NOTICE: For validating CNV calls of $startsnp-$endsnp, the prior distribution is set at ", join (",", convertDelDupFreqToPrior ($delfreq, $dupfreq, $backfreq)), "\n";
		}
	}
	newvalidateCNVCall (\@inputfile, $hmmfile, $pfbfile, $sexfile, $gcmodelfile, $directory, $candregion);
} elsif ($trio) {
	scalar (@inputfile) % 3 == 0 or pod2usage ("Error in argument: the inputfiles do not appear to form one or more complete trios (${\(scalar @inputfile)} inputifles provided)");
	for my $i (0 .. @inputfile/3-1) {
		my @trioinputfile = @inputfile[$i*3..($i*3+2)];
		print STDERR "NOTICE: Recalling CNVs for trio: @trioinputfile\n";
		newtestTrioCNVFile (\@trioinputfile, $hmmfile, $pfbfile, $cnvfile, $sexfile, $gcmodelfile, $directory);
	}
} elsif ($quartet) {
	scalar (@inputfile) % 4 == 0 or pod2usage ("Error in argument: the inputfiles do not appear to form one or more complete quartets (${\(scalar @inputfile)} inputifles provided)");
	for my $i (0 .. @inputfile/4-1) {
		my @quartetinputfile = @inputfile[$i*4..($i*4+3)];
		print STDERR "NOTICE: Recalling CNVs for quartet: @quartetinputfile\n";
		newtestQuartetCNVFile (\@quartetinputfile, $hmmfile, $pfbfile, $cnvfile, $sexfile, $gcmodelfile, $directory);
	}
} elsif ($exclude_heterosomic) {
	excludeHeterosomic ($cnvfile);
} elsif ($summary) {
	calculateSampleSummary (\@inputfile, $pfbfile, $gcmodelfile, $directory);
} elsif ($joint) {
	scalar (@inputfile) % 3 == 0 or pod2usage ("Error in argument: the inputfiles do not appear to form one or more complete trios (${\(scalar @inputfile)} inputifles provided)");
	for my $i (0 .. @inputfile/3-1) {
		my @familyfile = @inputfile[$i*3..($i*3+2)];
		jointCNVCall (\@familyfile, $hmmfile, $pfbfile, $sexfile, $chrx||0, $chry||0, $gcmodelfile, $directory);
	}
} elsif ($cctest) {
	cctestCNV ($cnvfile, $phenofile, $pfbfile);
} elsif ($tseq) {
	testSEQ (\@inputfile, $hmmfile);
}

#check the validity of arguments supplied to the program, and assign default values to various arguments

sub read_refgcfile {
	print $refgcfile."\n";
	open (REFGC, $refgcfile) or confess "\nERROR: cannot read from refgcfile $refgcfile: $!\n";
	chomp (@ref_median = <REFGC>);
	close REFGC;
}

sub processArgument {
	my @command_line = @ARGV;		#command line argument

		GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 
				'train'=>\$train, 'ctrain'=>\$ctrain, 'test|wgs'=>\$test, 'trio'=>\$trio, 'quartet'=>\$quartet, 'joint'=>\$joint, 'summary'=>\$summary, 'cctest'=>\$cctest, 'validate'=>\$validate, 
				'hmmfile=s'=>\$hmmfile, 'pfbfile=s'=>\$pfbfile, 'cnvfile=s'=>\$cnvfile, 'output=s'=>\$output, 'sexfile=s'=>\$sexfile, 'logfile=s'=>\$logfile, 
				'minsnp=i'=>\$minsnp, 'minlength=s'=>\$minlength, 'minconf=f'=>\$minconf,
				'chrx'=>\$chrx, 'chry'=>\$chry, 'medianadjust!'=>\$medianadjust, 'bafadjust!'=>\$bafadjust, 'sdadjust!'=>\$sdadjust, 
				'heterosomic_threshold=i'=>\$heterosomic_threshold, 'denovo_rate=f'=>\$denovo_rate,
				'fmprior=s'=>\$fmprior, 'listfile=s'=>\$listfile, 'exclude_heterosomic'=>\$exclude_heterosomic, 'record_count=i'=>\$record_count,
				'confidence'=>\$confidence, 'tabout'=>\$tabout, 'loh'=>\$loh, 
				'coordinate_from_input'=>\$coordinate_from_input, 'phenofile=s'=>\$phenofile,
				'onesided'=>\$onesided, 'control_label=s'=>\$control_label, 'type_filter=s'=>\$type_filter, 'gcmodelfile|gcwavefile=s'=>\$gcmodelfile,
				'directory=s'=>\$directory, 'flush!'=>\$flush, 'startsnp=s'=>\$startsnp, 'endsnp=s'=>\$endsnp, 'delfreq=f'=>\$delfreq, 'dupfreq=f'=>\$dupfreq,
				'backfreq=f'=>\$backfreq, 'candlist=s'=>\$candlist, 'region=s'=>\$region, 'bafxhet=f'=>\$bafxhet_threshold, 'valilog=s'=>\$valilog, 'cnprior=s'=>\$cnprior,
				'tseq'=>\$tseq, 'gamma=s'=>\$gamma, 'paramfile=s'=>\$paramfile, 'lastchr=i'=>\$lastchr, 'cache'=>\$cache, 'stroma=f'=>\$stroma,
				'refchr=s' => \$refchr, 'refgcfile=s' => \$refgcfile,
				) or pod2usage ();

	$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
	$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);

#CHECK THE NUMBER OF OPERATIONS SPECIFIED IN COMMAND LINE; ONE AND ONLY ONE OPERATION IS ALLOWED
	my $num_operation = 0;
	$train and $num_operation++;
	$ctrain and $num_operation++;
	$test and $num_operation++;
	$validate and $num_operation++;
	$trio and $num_operation++;
	$quartet and $num_operation++;
	$exclude_heterosomic and $num_operation++;
	$summary and $num_operation++;
	$joint and $num_operation++;
	$cctest and $num_operation++;
	$tseq and $num_operation++;
	$num_operation == 1 or pod2usage ("Error in argument: please specify one and only one operation such as --test, --trio, --quartet, --joint, --summary or --cctest");


#CHECK THE REQUIRED ARGUMENTS FOR EACH OPERATION; WARN THE USER ABOUT MISSING REQUIRED ARGUMENTS
	if ($train) {
		defined $hmmfile and defined $pfbfile and defined $output or pod2usage ("Error in argument: please specify the -hmmfile, --pfbfile and --output arguments for the --train operation");
	} elsif ($ctrain) {
		defined $hmmfile and defined $output or pod2usage ("Error in argument: please specify the -hmmfile and --output arguments for the --train operation");
	} elsif ($test) {
		defined $hmmfile and defined $pfbfile or pod2usage ("Error in argument: please specify the -hmmfile and --pfbfile arguments for the --test operation");
	} elsif ($validate) {
		defined $hmmfile and defined $pfbfile or pod2usage ("Error in argument: please specify the -hmmfile and --pfbfile arguments for the --validate operation");
		if ($candlist) {
			if ($cnprior) {
				pod2usage ("Error in argument: please do NOT specify --cnprior when --candlist is provided");
			}
			if ($startsnp or $endsnp or defined $delfreq or defined $dupfreq) {
				pod2usage ("Error in argument: please do NOT specify --startsnp or --endsnp or --delfreq or --dupfreq for --validate operation when --candlist is provided");
			}
		} elsif ($cnprior) {
			defined $startsnp and defined $endsnp or pod2usage ("Error in argument: please specify --startsnp and --endsnp for the --validate operation when --candlist is not provided");
			if (defined $delfreq or defined $dupfreq) {
				pod2usage ("Error in argument: please do NOT specify --delfreq or --dupfreq for --validate operation when --cnprior is provided");
			}
		} else {
			defined $delfreq || defined $dupfreq or pod2usage ("Error in argument: please specify --delfreq or --dupfreq or --cnprior for the --validate operation when --candlist is not provided");
			defined $startsnp and defined $endsnp or pod2usage ("Error in argument: please specify --startsnp and --endsnp for the --validate operation when --candlist is not provided");

			defined $delfreq and $delfreq <= 1 && $delfreq >=0 || pod2usage ("Error in argument: the --delfreq ($delfreq) must be between 0 and 1");
			defined $dupfreq and $dupfreq <= 1 && $dupfreq >=0 || pod2usage ("Error in argument: the --delfreq ($dupfreq) must be between 0 and 1");
			defined $delfreq and defined $dupfreq and $delfreq+$dupfreq<1 || pod2usage ("Error in argument: the sum of --delfreq ($delfreq) and --dupfreq ($dupfreq) must be less than 1");
		}
	} elsif ($trio) {
		defined $hmmfile and defined $pfbfile and defined $cnvfile or pod2usage ("Error in argument: please specify the --hmmfile and --pfbfile and --cnvfile arguments for the --trio operation");
	} elsif ($quartet) {
		defined $hmmfile and defined $pfbfile and defined $cnvfile or pod2usage ("Error in argument: please specify the --hmmfile and --pfbfile and --cnvfile arguments for the --quartet operation");
	} elsif ($exclude_heterosomic) {
		$cnvfile or pod2usage ("Error in argument: please specify --cnvfile for the --exclude_heterosomic operation");
	} elsif ($summary) {
		defined $pfbfile or pod2usage ("Error in argument: please specify the --pfbfile arguments for the --summary operation");
	} elsif ($joint) {
		defined $hmmfile and defined $pfbfile or pod2usage ("Error in argument: please specify the --hmmfile and --pfbfile arguments for the --joint operation");
	} elsif ($cctest) {
		defined $pfbfile and defined $cnvfile and $phenofile or pod2usage ("Error in argument: please specify the --pfbfile and --cnvfile and --phenofile arguments for the --cctest opearation");
	} elsif ($tseq) {
		defined $hmmfile or pod2usage ("Error in argument: please specify the -hmmfile arguments for the --tseq operation");
	}


#CHECK OPERATION-SPECIFIC ARGUMENTS: PREVENT USER FROM SPECIFYING USELESS ARGUMENTS TO THE PROGRAM
	if ($startsnp or $endsnp or defined $delfreq or defined $dupfreq or defined $backfreq or defined $candlist) {
		$validate or pod2usage ("Error in argument: the --startsnp, --endsnp, --delfreq, --dupfreq, --backfreq or --candlist argument is only supported for the --validate operation");
	}
	if ($cnvfile) {
		$trio or $quartet or $exclude_heterosomic or $cctest or pod2usage ("Error in argument: the --cnvfile argument is only supported for the --trio or --quartet or --exclude_heterosomic or --cctest operation");
	}
	if (defined $minconf or $confidence) {
		$test or $validate or pod2usage ("Error in argument: the --confidence or --minconf argument is supported only for the --test and --validate operation");
	}
	if ($loh) {
		print STDERR "WARNING in argument: the --loh argument is obselete and may be removed in future releases!!!\n";
		$test or pod2usage ("Error in argument: the --loh argument is supported only for the --test operation");
	}
	if (defined $type_filter or defined $control_label or defined $phenofile or $onesided) {
		$cctest or pod2usage ("Error in argument: the --type_filter, --control_label, --phenofile or --onesided argument is supported only for the --cctest operation");
	}
	if (defined $fmprior or defined $denovo_rate) {
		$trio or $quartet or pod2usage ("Error in argument: the --fmprior or --denovo_rate argument is supported only for the --trio and --quartet operation");
	}
	if (defined $cnprior) {
		$validate or pod2usage ("Error in argument: the --cnprior argument is supported only for the --validate operation");
	}
	if ($hmmfile) {
		$train or $ctrain or $test or $trio or $quartet or $joint or $tseq or $validate or pod2usage ("Error in argument: the --hmmfile argument is supported only for the --train, --ctrain, --test, --trio, --quartet, --joint, --validation or --tseq operations");
	}
	if (defined $directory) {
		$test or $trio or $quartet or $joint or $validate or $summary or pod2usage ("Error in argument: the --directory argument is supported only for the --test, --trio, --quartet, --joint, --validation or --summary operations");
	}
	if (defined $heterosomic_threshold) {
		$exclude_heterosomic or pod2usage ("Error in argument: the --heterosomic_threshold argument is only supported for the --exclude_heterosomic operation");
	}
	if (defined $record_count) {
		$ctrain or pod2usage ("Error in argument: the --record_count argument is only supported for the --ctrain operation");
	}
#CNV-calling related arguments should not be specified when performing non-calling operations
	if ($chrx or $chry) {
		$test or $trio or $quartet or $joint or $validate or pod2usage ("Error in argument: the --chrx or --chry argument is supported only for the --test, --trio, --quartet, --joint or --validation operations");
	}
	if ($medianadjust or $bafadjust or $sdadjust) {
		$test or $trio or $quartet or $joint or $validate or pod2usage ("Error in argument: the --medianadjust, --bafadjust or --sdadjust argument is supported only for the --test, --trio, --quartet, --joint or --validation operations");
	}
	if ($gcmodelfile) {
		$test or $trio or $quartet or $joint or $validate or pod2usage ("Error in argument: the --gcmodelfile argument is supported only for the --test, --trio, --quartet, --joint or --validation operations");
	}
	if (defined $refchr){
		defined $refgcfile or pod2usage("Error in argument: the --refgcfile argument must be specified if --refchr is specified.")
	}
	if (defined $refgcfile){
		defined $refchr or pod2usage("Error in argument: the --refchr argument must be specified if --refgcfile is specified.")
	}
	if (defined $refchr and defined $refgcfile){
		$test or $trio or $quartet or $joint or $validate or pod2usage ("Error in argument: the --refchr and --refgcfile argument is supported only for the --test, --trio, --quartet, --joint or --validation operations");
	}

	if (defined $refgcfile) {
		read_refgcfile (); 
	}else{
		$refchr = '11';
		@ref_median = qw/54.8207535282258 56.8381472081218 53.1218950320513 46.9484174679487 39.9367227359694 38.3365384615385 41.9867788461538 40.4431401466837 44.5320512820513 42.1979166666667 41.6984215561224 43.1598557692308 43.4388020833333 40.8104967948718 39.8444475446429 41.5357572115385 38.7496995192308 45.0213249362245 42.3251201923077 43.5287459935897 40.7440808354592 37.0492788461538 36.5006009615385 35.8518016581633 35.2767427884615 35.1972155448718 36.5286192602041 39.4890825320513 36.5779246794872 36.7275641025641 38.3256935586735 37.791266025641 41.1777844551282 41.950534119898 42.3639823717949 41.9208733974359 41.2061543367347 35.4974959935897 35.2123397435897 36.5101841517857 36.7135416666667 36.8268229166667 37.6945153061224 40.7453926282051 47.7049278846154 47.3233173076923 44.7361288265306 46.6585536858974 39.1593549679487 36.5684789540816 38.2718466806667 37.184425 37.184425 37.184425 37.184425 35.9227764423077 41.1157852564103 41.6662348533163 39.7402844551282 40.0149238782051 46.6417211415816 49.9136618589744 45.2016225961538 51.3019172512755 52.0818309294872 51.1320112179487 49.9807185102302 49.9807185102302 49.5874473187766 50.547349024718 50.7186498397436 45.6435347576531 46.3352363782051 42.4091546474359 46.6399274553571 43.7746394230769 45.0160256410256 41.8526642628205 43.8899075255102 38.5112179487179 36.1038661858974 36.1689851721939 39.8506610576923 37.0439703525641 36.8012595663265 40.2521033653846 39.661858974359 37.5013769093564 35.5448717948718 36.9039979272959 35.2046274038462 38.2195512820513 40.074537627551 40.7097355769231 40.5470753205128 38.4104380072343 36.131109775641 35.3915264423077 34.9693080357143 36.2953725961538 37.9602363782051 39.1942362882653 37.4464142628205 36.8879206730769 35.7242588141026 36.7556202168367 37.0639022435897 40.6929086538462 38.385084502551 39.4121594551282 40.2410857371795 42.0772879464286 43.2935697115385 43.2345753205128 40.9113919005102 44.9575320512821 46.2513020833333 46.4753069196429 48.3886217948718 47.8520633012821 43.8001802884615 39.808274872449 44.5042067307692 38.3835136217949 44.9097177933673 45.5366586538462 41.7346754807692 39.2198461415816 41.9489182692308 44.3351362179487 42.7910754145408 42.3190104166667 42.0425681089744 47.0514787946429 45.3482603740699/;
	}

	if ($tabout or $coordinate_from_input) {
		$test or $trio or $quartet or $joint or $validate or pod2usage ("Error in argument: the --tabout or --coordinate_from_input argument is supported only for the --test, --trio, --quartet, --joint or --validation operations");
	}
	if (defined $minsnp or defined $minlength) {
		$test or $trio or $quartet or $joint or $validate or pod2usage ("Error in argument: the --minsnp or --minlength argument is supported only for the --test, --trio, --quartet, --joint or --validation operations");
	}
	if (defined $region) {
		$test or $trio or $quartet or $joint or $validate or pod2usage ("Error in argument: the --region argument is supported only for the --test, --trio, --quartet, --joint or --validation operations");
	}
	if (defined $gamma) {
		$tseq or pod2usage ("Error in argument: the --gamma argument is supported only for the --tseq operation");
		$gamma =~ m/^[\d\.]+,[\d\.]+$/ or pod2usage ("Error in argument: the --gamma argument should be specified as two floating point numbers separated by comma");
		($gamma_k, $gamma_theta) = split (/,/, $gamma);
	}
	if (defined $paramfile) {
		defined $gamma and pod2usage ("Error in argument: please do not specify --gamma and --paramfile simultaneously");
		open (PARAM, $paramfile) or confess "Error: cannot read from --paramfile $paramfile: $!\n";
		while (<PARAM>) {
			if (m/k=([\d\.]+)\s+theta\d?=([\d\.]+)/) {
				($gamma_k, $gamma_theta) = ($1, $2);
			}
		}
		close (PARAM);
		$gamma_k or confess "Error: unable to read gamma values from the --paramfile $paramfile\n";
	}

#GET THE LIST OF INPUT FILE NAMES USED IN THE SUBSEQUENT ANALYSIS
	if (not $cctest and not $exclude_heterosomic) {		#these two operations do NOT require signal files as input
		if (@ARGV) {
			$listfile and pod2usage ("Error in argument: please do not specify the --listfile argument ($listfile) and provide signal file names (@ARGV) in command line simultaneously");
			@inputfile = @ARGV;
		} else {
			$listfile or pod2usage ("Error in argument: please either specify the --listfile argument or provide signal file names in command line\n");
			open (LIST, $listfile) or confess "\nERROR: cannot read from listfile $listfile: $!";
			while (<LIST>) {
				s/^\s+|\s*[\r\n]+$//g;		#get rid of leading and trailing spaces
					$_ or print STDERR "WARNING: blank lins in $listfile detected and skipped\n" and next;
				if ($trio or $joint) {
					m/^([^\t]+)\t([^\t]+)\t([^\t]+)$/ or confess "\nERROR: the --listfile should contain 3 tab-delimited file names per line for the --trio or --joint operation (invalid line found: <$_>)\n";
					push @inputfile, $1, $2, $3;
				} elsif ($quartet) {
					m/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/ or confess "\nERROR: the --listfile should contain 4 tab-delimited file names per line for the --quartet operation (invalid line found: <$_>)\n";
					push @inputfile, $1, $2, $3, $4;
				} else {
					m/^([^\t]+)$/ or confess "\nERROR: the --listfile should contain one file name per line without tab character (invalid line found: <$_>)\n";
					push @inputfile, $1;
				}
			}
			close (LIST);
		}
	}


	defined $flush or $flush = 1;
	$flush and $| = 1;					#flush input/output buffer to minitor program progress in real time (when using remote connections)

		if (defined $minlength) {
			$minlength =~ m/^\d+(k|m)?$/i or pod2usage ("Error in argument: the --minlength argument should be a positive integer (suffix of k or m is okay)");
			$minlength =~ s/k$/000/i;
			$minlength =~ s/m$/000000/i;
		}

	defined $minconf and $confidence = 1;							#automatically set the --confidence argument to 1 for confidence score calculation

		defined $hmmfile and readHMMFile ($hmmfile);						#check the validity of HMM file

		$minsnp ||= 3;										#by default only call CNVs containing more than or equal to 3 SNPs for Illumina HumanHap550 array
		$heterosomic_threshold ||= 5;								#by default if there are more than 5 CNVs in the same chromosome with same CN state, it is likely to be caused by cell-line artifacts
		$denovo_rate ||= 0.0001;								#by default de novo rate is 0.0001 (the rationale is that 1% individuals have denovo CNVs, but each individual may have 20-100 detectable CNVs, so the prior is somewhere around 0.0001 for any given CNV call)
		defined $medianadjust or $medianadjust = 1;						#by default --medianadjust is turned on, since it usually improve calling
		defined $sdadjust or $sdadjust = 1;							#as of Jan 2008, this is on by default (to reduce false positve rate for low-quality samples empirically)
		defined $bafadjust or $bafadjust = 1;							#as of June 2008, this is on by default (to reduce false positve rate for low-quality samples empirically)
		defined $bafxhet_threshold or $bafxhet_threshold = 0.1;					#as of Jan 2009, this becomes an argumen that user can change (by default is 10% markers have BAF around 0.5 the same is predicted as female; it works for Illumina SNP arrays well but not for Affy arrays).
		$bafxhet_threshold > 0 and $bafxhet_threshold < 1 or pod2usage ("Error in argument: the --bafxhet argument ($bafxhet_threshold) should be between 0 and 1");

	$chrx and $chry and pod2usage ("Error in argument: please do not specify --chrx and --chry simultaneously");
	if ($chrx or $chry) {
		$chrx and $chry and pod2usage ("Error in argument: please do NOT specify both --chrx and --chry argument");
		$region and pod2usage ("Error in argument: the --region argument should NOT be specified when --chrx or --chry is provided");
	}
	if (defined $sexfile) {
		$chrx or $chry or pod2usage ("Error in argument: the --sexfile argument should NOT be specified when neither --chrx nor --chry is provided");
	}
	if (defined $lastchr) {
		$lastchr > 0 or pod2usage ("Error in argument: the --lastchr argument must be a positive integer");
	} else {
		$lastchr = 22;
	}
	if ($chrx) {
		$region = 'X';
	} elsif ($chry) {
		$region = 'Y';
	} else {
		$region ||= "1-$lastchr";
	}

	if (defined $type_filter) {
		$type_filter eq 'dup' or $type_filter eq 'del' or pod2usage ("Error in argument: the --type_filter can be only 'del' or 'dup' but not '$type_filter'");
	}
	$control_label ||= 'control';								#default text label for control is "control" (or the number 1)

		if ($fmprior) {
			@fmprior = split (/,/, $fmprior);
			map {s/^$/0/} @fmprior;								#handle situations where no number is given for a state (blank between consecutive commas)
				@fmprior == 6 or pod2usage ("Error in argument: the --fmprior argument should be 6 comma-separated numbers that sum up to 1");
			abs (sum (\@fmprior) - 1) < 0.01 or pod2usage ("Error in argument: the --fmprior argument should be 6 comma-separated numbers that sum up to 1");
			unshift @fmprior, 0;								#fill in the first element of array so hidden state starts from subscript 1 to 6
		} else {
			@fmprior = (0, 0.01*0.5, 0.49*0.5, 0.5, 0, 0.49*0.5, 0.01*0.5);
		}
	if ($cnprior) {
		@cnprior = split (/,/, $cnprior);
		map {s/^$/0/} @cnprior;
		map {$_||=1e-9} @cnprior;
		@cnprior == 5 or pod2usage ("Error in argument: the --cnprior argument should be 5 comma-separated numbers, corresponding to probability for CN=0 to CN=4");
	}

	if (defined $backfreq) {
		$backfreq > 0 and $backfreq < 0.5 or pod2usage ("Error in argument: the --backfreq argument ($backfreq) should be less than 0.5");
	} else {
		$backfreq = 0.0001;								#by default, the background freq means that 0.01% of whole genome markers are in CNV
	}
	if ($validate) {
		defined $delfreq or $delfreq = $backfreq;
		defined $dupfreq or $dupfreq = $backfreq;
	}

	if (defined $output) {									#default output is STDOUT
		open (STDOUT, ">$output") or confess "\nERROR: cannot write to output file $output";
	}
	if ($logfile) {
		open (LOG, ">$logfile") or die "Error: cannot write to log file $logfile: $!\n";
		print LOG "PennCNV Version:\n\t", q/$LastChangedDate: 2013-02-08 11:10:50 -0800 (Fri, 08 Feb 2013) $/, "\n";
		print LOG "PennCNV Information:\n\tFor questions, comments, documentation, bug reports and program update, please visit http://www.openbioinformatics.org/penncnv/\n";
		print LOG "PennCNV Command:\n\t$0 @command_line\n";
		print LOG "PennCNV Started:\n\t", scalar (localtime), "\n";
		close (LOG);

		my $notee = 0;
		my $msg = qx/perl -v | tee 2>&1/;						#check whether tee is installed in the system
			if ($msg =~ m/tee/ or not $msg) {						#system error message should contain "tee" word
				$notee++;
			} else {
				eval {
					open(STDERR, " | tee $logfile 1>&2");				#redirect STDERR to a logfile, but also print out the message to STDERR
				};
				$@ and $notee++;
			}
		if ($notee) {
			$verbose and print STDERR "WARNING: cannot set up dual-output of notification/warning message to STDERR and to log file $logfile (unable to execute 'tee' command in your system)\n";
			print STDERR "NOTICE: All program notification/warning messages will be written to log file $logfile\n";
			open (STDERR, ">>$logfile") or confess "\nERROR: unable to write to log file $logfile: $!\n";
		} else {
			print STDERR "NOTICE: All program notification/warning messages that appear in STDERR will be also written to log file $logfile\n";
		}
	}
}

#this subroutine is used to generate an offspring_state array that contains all possible combinations of offspring states. Right now it only support up to two offsprings
sub generateOstateArray {
	my ($num_element, $symbol) = @_;
	my @return;

	if ($num_element == 2) {
		for my $i (@$symbol) {
			for my $j (@$symbol) {
				push @return, [$i, $j];
			}
		}
	} else {
		confess "\nERROR: the generateOstateArray subroutine can only handle two offsprings at this time (but you gave $num_element)\n";
	}
	return @return;
}

#this subroutine converts an index into three states (state range from 1 to 6)
sub convertTrioIndex2State {
	my ($index) = @_;
	my $fstate = int ($index / 49);
	my $mstate = int (($index - $fstate*49) / 7);
	my $ostate = ($index - $fstate*49) % 7;
	return ($fstate, $mstate, $ostate);
}

#this subroutine converts an index into four states (state range from 1 to 6)
sub convertQuartetIndex2State {
	my ($index) = @_;
	my $o2state = int ($index/343);
	my $o1state = int (int ($index - $o2state * 343) / 49);

	my $mstate = int (($index - $o2state*343 - $o1state*49) / 7);
	my $fstate = $index - $o2state*343 - $o1state*49 - $mstate*7;
	return ($fstate, $mstate, $o1state, $o2state);
}

#read transition matrix from the HMM file
sub readTransitionMatrix {
	my ($hmmfile) = @_;
	my ($found_transition, $transition);
	push @$transition, 0;					#convert 0-based index to 1-based index
		open (HMM, $hmmfile) or confess "\nERROR: cannot read from HMM file $hmmfile: $!\n";
	while (<HMM>) {
		m/^A:/ and ++$found_transition and next;
		m/^B:/ and last;
		if ($found_transition) {
			s/^\s+|\s*[\r\n]+$//g;
			my @trans = split (/\s+/, $_);
			unshift @trans, 0;			#convert 0-based index to 1-based index
				push @$transition, [@trans];
		}
	}
	close (HMM);
	return ($transition);
}

#from siginfo, get the LRR/BAF values for a specific genomic region
sub getRegionInfo {
	my ($siginfo, $pfbinfo, $nextregion, $startname, $endname) = @_;
	my ($region_name, $region_chr, $region_pos, $region_lrr, $region_baf, $region_pfb, $region_snpdist);
	my ($c, $s, $e) = @$nextregion;
	my ($found, $last);

	my $info = $siginfo->{$c};
	my ($name, $pos, $lrr, $baf) = ($info->{name}, $info->{pos}, $info->{lrr}, $info->{baf});
	for my $i (0 .. @$name-1) {
		if ($name->[$i] eq $s) {
			$found = 1;
		} elsif ($name->[$i] eq $e) {
			$last = 1;
		}
		if ($found) {
			push @$region_name, $name->[$i];
			push @$region_chr, $c;
			push @$region_pos, $pos->[$i];
			push @$region_lrr, $lrr->[$i];
			push @$region_baf, $baf->[$i];
			my @temp = split (/\t/, $pfbinfo->{$name->[$i]}); 
			push @$region_pfb, $temp[2];
			if ($i < @$name-1) {
				push @$region_snpdist, ($pos->[$i+1]-$pos->[$i]) || 1;
			} else {
				push @$region_snpdist, 100_000_000;
			}
		}
		$last and last;
	}

#generally speaking, we try to make the start of the block open (shift the startsnp), but the end of the block solid (keep the endsnp)
	if ($endname->{$s}) {
		shift @$region_name; shift @$region_chr; shift @$region_pos; shift @$region_lrr; shift @$region_baf; shift @$region_pfb; shift @$region_snpdist;
	}
	if ($startname->{$e}) {
		if (not $endname->{$e}) {

#when there is only one element (when startsnp is the end of previous block, and endsnp is begin of next block), arbitrarily create a single-SNP block
			if (@$region_name == 1) {
				$endname->{$e}++;		#mark this SNP as constituting an end SNP, even though it is not a endSNP in the CNV file
			} else {
				pop @$region_name; pop @$region_chr; pop @$region_pos; pop @$region_lrr; pop @$region_baf; pop @$region_pfb; pop @$region_snpdist;
			}
		}
	}
	return ($region_name, $region_chr, $region_pos, $region_lrr, $region_baf, $region_pfb, $region_snpdist);
}

#read the HMM file
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

sub trainCHMM {
	my ($ref_inputfile, $hmmfile, $pfbfile, $gcmodelfile, $directory) = @_;

	my ($pfbinfo) = newreadPFB ($pfbfile);
	my $hmm = readHMMFile ($hmmfile);
	my $gcmodel;
	defined $gcmodelfile and $gcmodel = newreadGCModel ($gcmodelfile, $pfbinfo);	#read the gc model for SNPs defined in snp_chr

		open (LOGR_B_PFB, ">$output.lrr_baf_pfb") or confess "\nERROR: cannot write to temporary file $output.lrr_baf_pfb: $!\n";

	my ($counter, $counter_file, $siginfo, $sigdesc);
	for my $inputfile (@$ref_inputfile) {
		$region ||= "1-$lastchr";

		($siginfo, $sigdesc) = readLRRBAF ($inputfile, $region, $pfbinfo, $gcmodel, $directory);

#if the file processing fails and no signal data is read from file, move on to the next file
		%$siginfo or next;

		print STDERR "NOTICE: quality summary for $inputfile: LRR_mean=", sprintf ("%.4f", $sigdesc->{lrr_mean}), " LRR_median=", sprintf ("%.4f", $sigdesc->{lrr_median}), " LRR_SD=", sprintf ("%.4f", $sigdesc->{lrr_sd}), 
			  " BAF_mean=", sprintf ("%.4f", $sigdesc->{baf_mean}), " BAF_median=", sprintf ("%.4f", $sigdesc->{baf_median}), " BAF_SD=", sprintf ("%.4f", $sigdesc->{baf_sd}), " BAF_DRIFT=", sprintf ("%.6f", $sigdesc->{baf_drift}), " WF=", sprintf("%.4f", $sigdesc->{wf}), "\n";
		$sigdesc->{lrr_sd} > $hmm->{B1_sd}[2]+0.05 || $sigdesc->{lrr_sd} < $hmm->{B1_sd}[2]-0.05 and print STDERR "NOTICE: Sample does not pass quality control for training due to its LRR_SD ($sigdesc->{lrr_sd})!\n" and next;
		$sigdesc->{baf_median} >= 0.55 || $sigdesc->{baf_median} <= 0.45 and print STDERR "NOTICE: Sample does not pass quality control for training due to its BAF_MEDIAN ($sigdesc->{baf_median})!\n" and next;
		$sigdesc->{baf_drift} >= $hmm->{B1_uf}/10 and print STDERR "NOTICE: Sample does not pass quality control for training due to its large BAF_DRIFT ($sigdesc->{baf_drift})!\n" and next;

		if ($sigdesc->{cn_count}) {
			if (not $hmm->{B3_mean}) {
				print STDERR "WARNING: The signal file $inputfile contains Non-Polymorphic markers but the HMM file $hmmfile is for SNP markers only. It is recommended to use a HMM file containing parameters specifically for non-polymorphic markers.\n";
			}
		}

#perform signal pre-processing to adjust the median LRR and BAF values
		if ($medianadjust) {				#THIS IS A MANDATORY ARGUMENT NOW, SINCE IT ALWAYS IMPROVE PERFORMANCE!!! (use "--nomedianadjust" to disable this feature)
			print STDERR "NOTICE: Median-adjusting LRR values for all autosome markers from $inputfile by ", sprintf ("%.4f", $sigdesc->{lrr_median}), "\n";
			adjustLRR ($siginfo, $sigdesc->{lrr_median});
			$sigdesc->{lrr_mean} -= $sigdesc->{lrr_median};
			$sigdesc->{lrr_median} = 0;
		}
		if ($bafadjust) {
			print STDERR "NOTICE: Median-adjusting BAF values for all autosome markers from $inputfile by ", sprintf ("%.4f", $sigdesc->{baf_median}-0.5), "\n";
			adjustBAF ($siginfo, $sigdesc->{baf_median}-0.5);
			$sigdesc->{baf_mean} -= $sigdesc->{baf_median}-0.5;
			$sigdesc->{baf_median} = 0.5;
		}

		for my $curchr (keys %$siginfo) {
			$curchr =~ m/^(\d+)$/ or next;			#only consider autosomes
				my $pos = $siginfo->{$curchr}{pos};
			my $lrr = $siginfo->{$curchr}{lrr};
			my $baf = $siginfo->{$curchr}{baf};
			my $pfb = $siginfo->{$curchr}{pfb};
			my @snpdist = @$pos;
			my $logprob = 0;
			for my $i (1 .. @snpdist-2) {
				$snpdist[$i] = $snpdist[$i+1]-$snpdist[$i];
			}
			$snpdist[@snpdist-1] = 100_000_000;
			for my $i (1 .. @$lrr-1) {				#first element in the array is ignored
				print LOGR_B_PFB $lrr->[$i], "\t", $baf->[$i], "\t", $pfb->[$i], "\t", $snpdist[$i], "\n";
				$counter++;
			}
		}
		$counter_file++;
	}
	$counter or confess "\nERROR: none of the files (@$ref_inputfile) pass quality control threshold for training HMM model\n";
	print STDERR "NOTICE: Total of $counter markers in $counter_file files will be used in training HMM models\n";

#the following paragraph sets up pseudocounts in the array to handle cases where the state is extremely rare, and it is unlikely that the training set contains such state
	my (@alllogr, @allbaf, @allpfb, @allsnpdist);
	push @alllogr, ($hmm->{'B1_mean'}[0]) x 100; push @allbaf, ('0.5') x 100; push @allsnpdist, ('5000')x99, 10_000_000;				#state1
		push @alllogr, ('0') x 1000; for (1..250) {push @allbaf, 0.001, 0.5, 0.5, 0.999}
	push @allpfb, ('0.5') x 1100; push @allsnpdist, ('5000')x999, 10_000_000;

	push @alllogr, ($hmm->{'B1_mean'}[1]) x 100; for (1..50) {push @allbaf, 0.001, 0.999} push @allsnpdist, ('5000')x99, 10_000_000;		#state2
		push @alllogr, ('0') x 1000; for (1..250) {push @allbaf, 0.001, 0.5, 0.5, 0.999}
	push @allpfb, ('0.5') x 1100; push @allsnpdist, ('5000')x999, 10_000_000;

	push @alllogr, ('0') x 100; for (1..50) {push @allbaf, 0, 1} push @allsnpdist, ('5000')x99, 10_000_000;						#state4
		push @alllogr, ('0') x 1000; for (1..250) {push @allbaf, 0.001, 0.5, 0.5, 0.999}
	push @allpfb, ('0.5') x 1100; push @allsnpdist, ('5000')x999, 10_000_000;

	push @alllogr, ($hmm->{'B1_mean'}[4]) x 100; for (1..25) {push @allbaf, 0.001, 0.33, 0.66, 0.999} push @allsnpdist, ('5000')x99, 10_000_000;	#state5
		push @alllogr, ('0') x 1000; for (1..250) {push @allbaf, 0.001, 0.5, 0.5, 0.999}
	push @allpfb, ('0.5') x 1100; push @allsnpdist, ('5000')x999, 10_000_000;

	push @alllogr, ($hmm->{'B1_mean'}[5]) x 100; for (1..25) {push @allbaf, 0.001, 0.75, 0.25, 0.5} push @allsnpdist, ('5000')x99, 10_000_000;	#state6
		push @alllogr, ('0') x 1000; for (1..250) {push @allbaf, 0.001, 0.5, 0.5, 0.999}
	push @allpfb, ('0.5') x 1100; push @allsnpdist, ('5000')x999, 10_000_000;

	@alllogr==@allsnpdist and @alllogr==@allbaf and @alllogr==@allbaf or confess "FATAL ERROR: discordance in array elements: ${\(scalar @alllogr)} vs ${\(scalar @allsnpdist)} vs ${\(scalar @allbaf)}\n";
	for my $i (0 .. @alllogr-1) {
		print LOGR_B_PFB $alllogr[$i], "\t", $allbaf[$i], "\t", $allpfb[$i], "\t", $allsnpdist[$i], "\n";
		$counter++;
	}
	close (LOGR_B_PFB);

	my $system_command = "$0 --ctrain --record_count $counter --hmmfile $hmmfile $output.lrr_baf_pfb --output $output";
	print STDERR "NOTICE: Handing over the command '$system_command' to your operating system to continue Baum-Welch training\n";
	exec ($system_command) or confess "FATAL ERROR: Unable to execute system command <$system_command> via your operating system\n";
}

sub testSEQ {
	my ($ref_inputfile, $hmmfile) = @_;

	for my $inputfile (@$ref_inputfile) {
		my @mlstate;
		my $fh_temp = khmm::fopen ($inputfile, "r");
		my $record_count = 0;
		my ($curchr, $name);
		push @$name, 'padding_element';

		open (SIG, $inputfile) or confess "Error: cannot read from signal file $inputfile: $!\n";
		$_ = <SIG>;
		while (<SIG>) {
			$record_count++;
			m/^(\S+)/ or confess "Error: invalid record found: <$_>\n";
			push @$name, $1;
			push @mlstate, 0;
		}
		close (SIG);

		$name->[1] =~ m/^(\w+)/ or confess "Error: invalid marker name";
		$curchr = $1;

		$record_count or confess "\nERROR: the signal file $inputfile does not contain any signal values for CNV calling\n";
		print STDERR "NOTICE: Calling CNVs on $record_count records\n";

		my $hmm_model = khmm::ReadCHMM ($hmmfile);

		khmm::callCNVFromFile_SEQ ($hmm_model, $record_count, $fh_temp, \@mlstate, $gamma_k||0, $gamma_theta||0);
		print STDERR "mlstate(1000-1500)=@mlstate[1000..1500]\n";

		my ($normal_state, $found_signal, $cnv, $startpos, $endpos, $stretch_start_i) = (3);
		for my $i (1 .. @mlstate-1) {
			if ($mlstate[$i] != $normal_state) {
				if ($found_signal and $found_signal ne $mlstate[$i]) {
					$name->[$stretch_start_i] =~ m/^\w+:(\d+)/ or confess "Error: invalid name";
					$startpos = $1;
					$name->[$i-1] =~ m/^\w+:\d+\-(\d+)/ or confess "Error: invalid name";		#"
						$endpos = $1;
					push @{$cnv->{$found_signal}}, [$curchr, $startpos, $endpos, $i-$stretch_start_i, $name->[$stretch_start_i], $name->[$i-1]];
					$stretch_start_i = $i;
					$found_signal = $mlstate[$i];
				} elsif ($found_signal) {
					1;
				} else {
					$found_signal = $mlstate[$i];
					$stretch_start_i = $i;
				}
			} else {
				if ($found_signal) {
					$name->[$stretch_start_i] =~ m/^\w+:(\d+)/ or confess "Error: invalid name";
					$startpos = $1;
					$name->[$i-1] =~ m/^\w+:\d+\-(\d+)/ or confess "Error: invalid name";		#"
						$endpos = $1;
					push @{$cnv->{$found_signal}}, [$curchr, $startpos, $endpos, $i-$stretch_start_i, $name->[$stretch_start_i], $name->[$i-1]];
					$found_signal = 0;
				}
			}
		}
		if ($found_signal) {
			$name->[$stretch_start_i] =~ m/^\w+:(\d+)/ or confess "Error: invalid name";
			$startpos = $1;
			$name->[@$name-1] =~ m/^\w+:\d+\-(\d+)/ or confess "Error: invalid name";		#"
				$endpos = $1;
			push @{$cnv->{$found_signal}}, [$curchr, $startpos, $endpos, @$name-1-$stretch_start_i, $name->[$stretch_start_i], $name->[@$name-1]];
		}

#print all CNV calls, sorted by copy numbers first, then by chromosomes
		if ($tabout) {
			printTabbedCNV ($cnv, $inputfile);				#use tab-delimited output
		} else {
			printFormattedCNV ($cnv, $inputfile);
		}
	}
}

sub ctrainCHMM {
	my ($sigfile, $hmmfile, $record_count) = @_;

	my $fh_stdout = khmm::fh_stdout ();
	my $hmm_model = khmm::ReadCHMM ($hmmfile);

	print "<--------------------ORIGINAL HMM-------------------------\n";
	khmm::PrintCHMM ($fh_stdout, $hmm_model);
	print "--------------------------------------------------------->\n";

	if (not $record_count) {		#when record_count is not given, open the file and count the number of lines
		$record_count = 0;
		open (SIG, $sigfile) or confess "\nERROR: cannot read from signal file $sigfile: $!\n";
		while (<SIG>) {
			$record_count++;
		}
		close (SIG);
		$record_count or confess "\nERROR: the signal file $sigfile does not contain any signal values for Baum-Welch training\n";
	}

	my ($logprobinit, $logprobfinal, $num_iteration) = (0, 0, 0, 0);
	my $fh1 = khmm::fopen ($sigfile, "r");
	khmm::estHMMFromFile_CHMM ($hmm_model, $record_count, $fh1, \$num_iteration, \$logprobinit, \$logprobfinal);

	print STDERR "NOTICE: Baum-Welch estimation done with initp=$logprobinit endp=$logprobfinal delta=${\($logprobfinal-$logprobinit)} num_iteration=$num_iteration\n";
	print STDERR "NOTICE: The final HMM model is below (the model is also written to $output.hmm file):\n";
	khmm::PrintCHMM ($fh_stdout, $hmm_model);

	my $fh2 = khmm::fopen ("$output.hmm", 'w');
	khmm::PrintCHMM ($fh2, $hmm_model);

	khmm::FreeCHMM ($hmm_model);
	khmm::fclose ($fh1);
	khmm::fclose ($fh2);
}

sub calculateSampleSummary {
	my ($ref_inputfile, $pfbfile, $gcmodelfile, $directory) = @_;
	my ($pfbinfo) = newreadPFB ($pfbfile);
	my $gcmodel;

	defined $gcmodelfile and $gcmodel = newreadGCModel ($gcmodelfile, $pfbinfo);	#read the gc model for SNPs defined in snp_chr

		for my $inputfile (@$ref_inputfile) {
			my ($siginfo, $sigdesc) = readLRRBAF ($inputfile, "1-$lastchr", $pfbinfo, $gcmodel, $directory);

#if the file processing fails and no signal data is read from file, move on to the next file
			%$siginfo or next;

			print "$inputfile";
			print "\tLRR_MEAN=", sprintf("%.4f", $sigdesc->{lrr_mean});
			print "\tLRR_MEDIAN=", sprintf("%.4f", $sigdesc->{lrr_median});
			print "\tLRR_SD=", sprintf ("%.4f", $sigdesc->{lrr_sd});
			print "\tLRR_XMEAN=", sprintf("%.4f", $sigdesc->{lrr_xmean});
			print "\tLRR_XMEDIAN=", sprintf("%.4f", $sigdesc->{lrr_xmedian});
			print "\tLRR_XSD=", sprintf ("%.4f", $sigdesc->{lrr_xsd});
			print "\tBAF_MEAN=", sprintf ("%.4f", $sigdesc->{baf_mean});
			print "\tBAF_MEDIAN=", sprintf ("%.4f", $sigdesc->{baf_median});
			print "\tBAF_SD=", sprintf ("%.4f", $sigdesc->{baf_sd});
			print "\tBAF_DRIFT=", sprintf ("%.6f", $sigdesc->{baf_drift});		#this number is usually small
				print "\tBAF_XHET=", sprintf ("%.4f", $sigdesc->{baf_xhet});
			print "\tWF=", sprintf ("%.4f", $sigdesc->{wf});
			print "\n";
		}
}

#the new testCHMM subroutine treat each chromosome separately, to reduce the burden of system memory, when analyzing a large number of markers
sub newtestCHMM {
	my ($ref_inputfile, $hmmfile, $pfbfile, $sexfile, $gcmodelfile, $directory) = @_;
	my $pfbinfo = newreadPFB ($pfbfile);
	my $hmm = readHMMFile ($hmmfile);
	my ($file_sex, $gcmodel) = ({});

	defined $sexfile and $file_sex = readSexFile ($sexfile);
	defined $gcmodelfile and $gcmodel = newreadGCModel ($gcmodelfile, $pfbinfo);	#read the GC model for SNPs defined in snp_chr

		for my $inputfile (@$ref_inputfile) {
			my ($siginfo, $sigdesc, $sample_sex);

			($siginfo, $sigdesc) = readLRRBAF ($inputfile, $region, $pfbinfo, $gcmodel, $directory);

#if the file processing fails and no signal data is read from file, move on to the next file
			%$siginfo or print STDERR "WARNING: Skipping $inputfile since no signal values can be retrieved from the file\n" and next;

			$sample_sex = QCSignal ($inputfile, $siginfo, $sigdesc, $hmm, $file_sex->{$inputfile});

#read HMM model file
			my $hmm_model = khmm::ReadCHMM ($hmmfile);

#Now a mandatory step to handle low-quality samples and reduce false-positve calls, through matching the SD measure for sample and model
			if ($sdadjust) {
				if ($chrx) {
					khmm::adjustBSD ($hmm_model, $sigdesc->{lrr_xsd});
				} elsif ($chry) {
					khmm::adjustBSD ($hmm_model, $sigdesc->{lrr_ysd});
				} else {
					khmm::adjustBSD ($hmm_model, $sigdesc->{lrr_sd});
				}
			}

#do the CNV calling for each chromosome separately, save the CNV calls in $cnvcall hash
			my $cnvcall = {};
			for my $curchr (keys %$siginfo) {
				my $name = $siginfo->{$curchr}{name};
				my $pos = $siginfo->{$curchr}{pos}; 
				my $lrr = $siginfo->{$curchr}{lrr}; 
				my $baf = $siginfo->{$curchr}{baf};

				@$pos >= 10 or print STDERR "WARNING: Skipping chromosome $curchr due to insufficient data points (<10 markers)!\n" and next;

				my $pfb = [@{$siginfo->{$curchr}{pfb}}];	#we will use PFB again later (for assigning confidence scores)
					my $snpdist = [@$pos];
				my $logprob = 0;
				my $curcnvcall = {};
				for my $i (1 .. @$snpdist-2) {
					$snpdist->[$i] = ($snpdist->[$i+1]-$snpdist->[$i]) || 1;	#sometimes two markers have the same chromosome location (in Affymetrix array annotation)
				}

#generate CNV calls
				my $probe_count = scalar (@$lrr)-1;
				khmm::testVit_CHMM ($hmm_model, $probe_count, $lrr, $baf, $pfb, $snpdist, \$logprob);
				analyzeStateSequence ($curcnvcall, $curchr, $pfb, $name, $pos, $sample_sex);

#assign confidence scores (if the --conf argument is given)
				if ($confidence) {
					$pfb = $siginfo->{$curchr}{pfb};			#reassign value to PFB as it has been changed
						assignConfidence ($curcnvcall, $hmm_model, $name, $lrr, $baf, $pfb, $snpdist);
				}

#add the new CNV calls for one chromosome to the total CNV calls
				for my $key (keys %$curcnvcall) {
					push @{$cnvcall->{$key}}, @{$curcnvcall->{$key}};
				}
			}

#print all CNV calls, sorted by copy numbers first, then by chromosomes
			if ($tabout) {
				printTabbedCNV ($cnvcall, $inputfile);				#use tab-delimited output
			} else {
				printFormattedCNV ($cnvcall, $inputfile);
			}

			khmm::FreeCHMM  ($hmm_model);
		}
}

#the new tumorCHMM is copied from newtestCHMM, with tiny changes.
sub newtumorCHMM {
	my ($ref_inputfile, $hmmfile, $pfbfile, $sexfile, $gcmodelfile, $directory) = @_;
	my $pfbinfo = newreadPFB ($pfbfile);
	my $hmm = readHMMFile ($hmmfile);
	my ($file_sex, $gcmodel) = ({});

	defined $sexfile and $file_sex = readSexFile ($sexfile);
	defined $gcmodelfile and $gcmodel = newreadGCModel ($gcmodelfile, $pfbinfo);	#read the GC model for SNPs defined in snp_chr

		for my $inputfile (@$ref_inputfile) {
			my ($siginfo, $sigdesc, $sample_sex);

			($siginfo, $sigdesc) = readLRRBAF ($inputfile, $region, $pfbinfo, $gcmodel, $directory);

#if the file processing fails and no signal data is read from file, move on to the next file
			%$siginfo or print STDERR "WARNING: Skipping $inputfile since no signal values can be retrieved from the file\n" and next;

			$sample_sex = QCSignal ($inputfile, $siginfo, $sigdesc, $hmm, $file_sex->{$inputfile});

#read HMM model file
			my $hmm_model = khmm::ReadCHMM ($hmmfile);

#Now a mandatory step to handle low-quality samples and reduce false-positve calls, through matching the SD measure for sample and model
			if ($sdadjust) {
				if ($chrx) {
					khmm::adjustBSD ($hmm_model, $sigdesc->{lrr_xsd});
				} elsif ($chry) {
					khmm::adjustBSD ($hmm_model, $sigdesc->{lrr_ysd});
				} else {
					khmm::adjustBSD ($hmm_model, $sigdesc->{lrr_sd});
				}
			}

#do the CNV calling for each chromosome separately, save the CNV calls in $cnvcall hash
			my $cnvcall = {};
			for my $curchr (keys %$siginfo) {
				my $name = $siginfo->{$curchr}{name};
				my $pos = $siginfo->{$curchr}{pos}; 
				my $lrr = $siginfo->{$curchr}{lrr}; 
				my $baf = $siginfo->{$curchr}{baf};

				@$pos >= 10 or print STDERR "WARNING: Skipping chromosome $curchr due to insufficient data points (<10 markers)!\n" and next;

				my $pfb = [@{$siginfo->{$curchr}{pfb}}];	#we will use PFB again later (for assigning confidence scores)
					my $snpdist = [@$pos];
				my $logprob = 0;
				my $curcnvcall = {};
				for my $i (1 .. @$snpdist-2) {
					$snpdist->[$i] = ($snpdist->[$i+1]-$snpdist->[$i]) || 1;	#sometimes two markers have the same chromosome location (in Affymetrix array annotation)
				}

#generate CNV calls
				my $probe_count = scalar (@$lrr)-1;


#------------------------------------------------------------------------------------------
#khmm::testVit_CHMM ($hmm_model, $probe_count, $lrr, $baf, $pfb, $snpdist, \$logprob);
				khmm::tumorVit_CHMM ($hmm_model, $probe_count, $lrr, $baf, $pfb, $snpdist, \$logprob, $stroma);
#------------------------------------------------------------------------------------------


				analyzeStateSequence ($curcnvcall, $curchr, $pfb, $name, $pos, $sample_sex);

#assign confidence scores (if the --conf argument is given)
				if ($confidence) {
					$pfb = $siginfo->{$curchr}{pfb};			#reassign value to PFB as it has been changed
						assignConfidence ($curcnvcall, $hmm_model, $name, $lrr, $baf, $pfb, $snpdist);
				}

#add the new CNV calls for one chromosome to the total CNV calls
				for my $key (keys %$curcnvcall) {
					push @{$cnvcall->{$key}}, @{$curcnvcall->{$key}};
				}
			}

#print all CNV calls, sorted by copy numbers first, then by chromosomes
			if ($tabout) {
				printTabbedCNV ($cnvcall, $inputfile);				#use tab-delimited output
			} else {
				printFormattedCNV ($cnvcall, $inputfile);
			}

			khmm::FreeCHMM  ($hmm_model);
		}
}

#this subroutine is used to assign a confidence score to each predicted CNV
#this is an experimental feature. so far it seems to work but in the future the derivation of confidence score may change significantly
sub assignConfidence {
	my ($curcnvcall, $hmm_model, $name, $lrr, $baf, $pfb, $snpdist) = @_;
	for my $key (keys %$curcnvcall) {
		$key == 4 and next;				#do not consider LOH
			for my $stretch (@{$curcnvcall->{$key}}) {
				my ($curchr, $curstart, $curend, $curnumsnp, $curnamestart, $curnameend, $curid, $curpath, $curconf) = @$stretch;
				my (@region_lrr, @region_baf, @region_pfb, @region_snpdist);
				push @region_lrr, 0;
				push @region_baf, 0;
				push @region_pfb, 0;
				push @region_snpdist, 0;
				my ($curlogprob, $logprob, @logprob) = (0, 0);	#curlogprob is the value at desired state, logprob are temp values
					my ($found, $last) = qw/0 0/;
				for my $i (0 .. @$name-1) {
					if ($name->[$i] eq $curnamestart) {
						$found = 1;
					}
					if ($name->[$i] eq $curnameend) {	#use "if" not "elsif", because the CNV may contain only one SNP
						$last = 1;
					}
					if ($found) {
						push @region_lrr, $lrr->[$i];
						push @region_baf, $baf->[$i];
						push @region_pfb, $pfb->[$i];
						push @region_snpdist, $snpdist->[$i];
					}
					$last and last;
				}
				@region_lrr or confess "UNKNOWN ERROR: unable to find CNV region startsnp=$curnamestart endsnp=$curnameend\n";

				for my $nextstate (1, 2, 3, 5, 6) {		#do not consider LOH
					khmm::GetStateProb_CHMM ($hmm_model, @region_lrr-1, \@region_lrr, \@region_baf, \@region_pfb, \@region_snpdist, \$logprob, $nextstate);
					push @logprob, $logprob;
					$key == $nextstate and $curlogprob = $logprob;
				}
				@logprob = sort {$b<=>$a} @logprob;
				@$stretch = ($curchr, $curstart, $curend, $curnumsnp, $curnamestart, $curnameend, $curid, $curpath, sprintf ("%.3f", $curlogprob-$logprob[1]));
			}
	}
}

#print out formatted CNV calls (one CNV per line, with chromosome coordinates, length, copy number, file name, start and end SNP name and other supplementary information)
sub printFormattedCNV {
	my ($cnvcall, $inputfile) = @_;
	for my $key (sort {$a<=>$b} keys %$cnvcall) {
		$key == 4 and $loh || next;						#skip LOH regions unless --loh argument is set
			my $curcn = $key-1;
		$curcn >= 3 and $curcn--;						#calculate current copy number
			for my $stretch (sort {sortChr ($a->[0], $b->[0])} @{$cnvcall->{$key}}) {
				my ($curchr, $curstart, $curend, $curnumsnp, $curnamestart, $curnameend, $trioid, $triopath, $confscore) = @$stretch;

				if (defined $minsnp) {						#minimum number of SNPs within a CNV to be printed out
					$curnumsnp >= $minsnp or next;
				}
				if (defined $minlength) {					#minimum length of CNVs to be printed out
					$curend-$curstart+1 >= $minlength or next;
				}
				if (defined $minconf) {						#minimum confidence of CNVs to be printed out
					$confscore >= $minconf or next;
				}

				my $cnvregion = "chr$curchr:$curstart-$curend";
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

				print "$cnvregion numsnp=$curnumsnp length=$cnvlength state$key,cn=$curcn $inputfile startsnp=$curnamestart endsnp=$curnameend";
				print $trioid?" $trioid":"", $triopath?" statepath=$triopath":"", (defined $confscore)?" conf=$confscore":"", "\n";
			}
	}
}

#print out formatted CNV calls (one CNV per line, with chromosome coordinates, length, copy number, file name, start and end SNP name and other supplementary information)
sub printTabbedCNV {
	my ($cnvcall, $inputfile) = @_;
	for my $key (sort {$a<=>$b} keys %$cnvcall) {
		$key == 4 and $loh || next;						#skip LOH regions unless --loh argument is set
			my $curcn = $key-1;
		$curcn >= 3 and $curcn--;						#calculate current copy number
			for my $stretch (sort {sortChr ($a->[0], $b->[0])} @{$cnvcall->{$key}}) {
				my ($curchr, $curstart, $curend, $curnumsnp, $curnamestart, $curnameend, $trioid, $triopath, $confscore) = @$stretch;
				my $cnvlength = $curend-$curstart+1;

				if (defined $minsnp) {						#minimum number of SNPs within a predicted CNV to be printed out
					$curnumsnp >= $minsnp or next;
				}
				if (defined $minlength) {					#minimum length of CNVs to be printed out
					$cnvlength >= $minlength or next;
				}
				if (defined $minconf) {
					$confscore >= $minconf or next;
				}
				print "$curchr\t$curstart\t$curend\t$curcn\t$inputfile\t$curnamestart\t$curnameend\t$confscore\n";
			}
	}
}


sub analyzeStateSequence {
	my ($cnv, $newchr, $symbol, $name, $pos, $sample_sex) = @_;
	my ($normal_state, $found_signal, $stretch_start_i);

	$sample_sex and $sample_sex eq 'male' || $sample_sex eq 'female' || confess "\nERROR: sample_sex argument to the subroutine can be only 'male' or 'female'";

#set up the normal state in case of sex chromosome analysis
	if ($newchr eq 'X' and $sample_sex eq 'male') {
		$normal_state = 2;
	} elsif ($newchr eq 'Y') {
		if ($sample_sex eq 'male') {
			$normal_state = 2;
		} else {
			$normal_state = 1;
		}
	} else {
		$normal_state = 3;
	}

	$cnv ||= {};									#cnv is a reference to a hash (key=state value=array of info)
		for my $i (1 .. @$symbol-1) {							#the first element is zero (arrays start from 1 in khmm module)
#found a new CNV or continue within a previously identified stretch of CNV
			if ($symbol->[$i] != $normal_state) {
				if ($found_signal and $found_signal ne $symbol->[$i]) {		#transition to one CNV to another CNV with different copy number
					push @{$cnv->{$found_signal}}, [$newchr, $pos->[$stretch_start_i], $pos->[$i-1], $i-$stretch_start_i, $name->[$stretch_start_i], $name->[$i-1]];
					$stretch_start_i = $i;
					$found_signal = $symbol->[$i];
				} elsif ($found_signal) {					#do nothing, continue within a previously identified stretch of CNV
					1;
				} else {							#found the start of a new CNV
					$found_signal = $symbol->[$i];
					$stretch_start_i = $i;
				}
			} else {
				if ($found_signal) {						#end a previously identified CNV
					push @{$cnv->{$found_signal}}, [$newchr, $pos->[$stretch_start_i], $pos->[$i-1], $i-$stretch_start_i, $name->[$stretch_start_i], $name->[$i-1]];
					$found_signal = 0;
				}
			}
		}
	if ($found_signal) {								#finish the last CNV
		push @{$cnv->{$found_signal}}, [$newchr, $pos->[$stretch_start_i], $pos->[@$symbol-1], @$symbol-$stretch_start_i, $name->[$stretch_start_i], $name->[@$symbol-1]];
	}
}

sub cctestCNV {
	my ($cnvfile, $phenofile, $pfbfile) = @_;

	my (%chr_marker, %snp_index, @snpid, $index);
	open (PFB, $pfbfile) or confess "\nERROR: cannot read from PFB file $pfbfile: $!\n";
	while (<PFB>) {
		s/\s*[\r\n]+$//;
		my @record = split (/\s+/, $_);			#snp name, chr, position, pfb
			$record[1] =~ m/^\d+$/ or $record[1] =~ m/^[XY]+$/ or next;
		push @{$chr_marker{$record[1]}}, [$record[2], $record[0]];
	}
	close (PFB);
	for my $nextchr (sort {sortChr ($a, $b)} keys %chr_marker) {
		@{$chr_marker{$nextchr}} = sort {$a->[0] <=> $b->[0]} @{$chr_marker{$nextchr}};
		push @snpid, map {$_->[1]} @{$chr_marker{$nextchr}};
	}
	for my $i (0 .. @snpid-1) {
		$snp_index{$snpid[$i]} = $i;
	}

	my ($caseid, $controlid) = readPhenoFile ($phenofile, $control_label);
	my ($casecount, $controlcount) = (scalar keys %$caseid, scalar keys %$controlid);
	my (%snpcase, %snpcontrol, %skipped_ind);

	$casecount and $controlcount or confess "\nERROR: unable to find sufficient case/control from the phenotype file $phenofile (case=$casecount; control=$controlcount)\n";

	open (CNV, $cnvfile) or confess "\nERROR: cannot read from CNV file $cnvfile: $!\n";
	while (<CNV>) {
		m/(?:chr)?(\w+).+?state\d,cn=(\d)\s+(.+?)\s+startsnp=([\w\-\.\:]+)\s+endsnp=([\w\-\.\:]+)/ or confess "\nERROR: invalid record found in cnvfile $cnvfile: <$_>\n";
		my ($cn, $indfile, $startsnp, $endsnp) = ($2, $3, $4, $5);

		if (not $caseid->{$indfile} and not $controlid->{$indfile}) {
			$skipped_ind{$indfile}++;
			next;
		}

		my ($starti, $endi) = ($snp_index{$startsnp}, $snp_index{$endsnp});
		defined $starti and defined $endi or confess "\nERROR: the marker $startsnp or $endsnp cannot be found in PFB file $pfbfile\n";
		if ($type_filter) {
			if ($type_filter eq 'dup') {
				$cn >= 3 or next;
			} elsif ($type_filter eq 'del') {
				$cn < 3 or next;
			}
		}
		for my $index ($starti .. $endi) {
			if ($caseid->{$indfile}) {
				$snpcase{$snpid[$index]}++;
				$snpcontrol{$snpid[$index]}||=0;
			} elsif ($controlid->{$indfile}) {
				$snpcontrol{$snpid[$index]}++;
				$snpcase{$snpid[$index]}||=0;
			}
		}
	}
	print STDERR "NOTICE: Finished reading CNV file $cnvfile";
	%skipped_ind and print STDERR " (${\(scalar keys %skipped_ind)} individuals are skipped due to lack of phenotype annotation)";
	print STDERR "\n";

#perform case-control association tests
	print STDERR "NOTICE: Performing case-control comparison by Fisher's Exact Test (8 output columns are: SNP, Chr, Position, case-cnv, case-normal, control-cnv, control-normal, P-value)\n";
	for my $nextchr (sort {sortChr ($a, $b)} keys %chr_marker) {
		for my $i (0 .. @{$chr_marker{$nextchr}}-1) {
			my ($nextpos, $nextsnp) = @{$chr_marker{$nextchr}->[$i]};
			exists $snpcase{$nextsnp} or next;		#it could be zero in cases!!!
				my $fisherp;
			if ($onesided) {
				if ($snpcase{$nextsnp} / $casecount > $snpcontrol{$nextsnp} / $controlcount) {
					$fisherp = khmm::fisher_exact_1sided ($snpcase{$nextsnp}, $casecount-$snpcase{$nextsnp}, $snpcontrol{$nextsnp}, $controlcount-$snpcontrol{$nextsnp});
				} else {
					$fisherp = 1;
				}
			} else {
				$fisherp = khmm::fisher_exact_2sided ($snpcase{$nextsnp}, $casecount-$snpcase{$nextsnp}, $snpcontrol{$nextsnp}, $controlcount-$snpcontrol{$nextsnp});
			}
			print STDOUT join ("\t", $nextsnp, $nextchr, $nextpos, $snpcase{$nextsnp}, $casecount-$snpcase{$nextsnp}, $snpcontrol{$nextsnp}, $controlcount-$snpcontrol{$nextsnp}, $fisherp), "\n";
		}
	}
}

sub readPhenoFile {
	my ($phenofile, $control_label) = @_;
	my (%case, %control);
	open (PHENO, $phenofile) or confess "\nERROR: cannot read from phenotype file $phenofile: $!\n";
	while (<PHENO>) {
		s/[\r\n]+$//;
		s/^\s+//;
		my @record = split (/\s+/, $_);
		@record == 2 or confess "\nERROR: invalid record found in $phenofile (2 records expected): <$_>\n";
		if ($record[1] eq '1' or $record[1] eq $control_label) {
			$case{$record[0]} and confess "\nERROR: individual $record[0] occur more than once in phenotype file $phenofile with discordant annotation\n";			
			$control{$record[0]}++;
		} else {
			$control{$record[0]} and confess "\nERROR: individual $record[0] occur more than once in phenotype file $phenofile with discordant annotation\n";
			$case{$record[0]}++;
		}
	}
	close (PHENO);
	return (\%case, \%control);
}		
sub newreadPFB {
	my ($pfbfile) = @_;
	my (%pfbinfo);
	my (%skip_snp, $skip_snp_count, $count_np);
	open (PFB, $pfbfile) or confess "\nERROR: cannot read from pfb file $pfbfile: $!\n";
	print STDERR "NOTICE: Reading marker coordinates and population frequency of B allele (PFB) from $pfbfile ...";	
	while (<PFB>) {
		s/[\r\n]+$//;							#delete line feed and return characters
			m/^(\S+)\t(\S+)\t(\S+)\t(\S+)$/ or confess "\nERROR: invalid record found in PFB file $pfbfile (4 tab-delimited records expected): <$_>\n";
		my ($name, $chr, $pos, $pfb) = ($1, $2, $3, $4);
		$chr =~ m/^chr/i and next;					#this is the header line in a regular PFB file
			if (not $chr =~ m/^(\d+|X|Y)$/ or $chr eq '0') {
				$skip_snp{$chr}++;
				$skip_snp_count++;
				next;
			}
		if ($pfb > 1) {
			$count_np++;
		} elsif ($pfb < 0) {						#when PFB<0, this marker is marked as being excluded from analysis
			next;
		}
		$pfb < 0.01 and $pfb = 0.01;
		$pfb > 0.99 and $pfb = 0.99;
		length ($pfb) > 6 and $pfb = sprintf("%.4f", $pfb);
		$pfbinfo{$name} = "$chr\t$pos\t$pfb";
	}
	print STDERR " Done with ${\(scalar keys %pfbinfo)} records";
	if ($skip_snp_count) {
		print STDERR " ($skip_snp_count records in chr ", join (",", keys %skip_snp), " were discarded)\n";
	} else {
		print STDERR "\n";
	}
	$verbose and $count_np and print STDERR "NOTICE: FPB file $pfbfile contains $count_np non-polymorphic markers\n";
	close (PFB);
	return (\%pfbinfo);
}

sub readLRRBAF  {
	my ($inputfile, $region, $pfbinfo, $gcmodel, $directory) = @_;
	print STDERR "NOTICE: Reading LRR and BAF values for from $inputfile ...";

	my (%chr_include);
	my @element = split (/,/, $region);
	for my $element (@element) {
		if ($element =~ m/^(\d+)$/) {
			$1 >= 1 and $1 <= $lastchr or pod2usage ("Error in argument: invalid --region argument that contains chromosome $1 (valid chromosomes are 1-$lastchr,X,Y,XY,M");
			$chr_include{$1} = 1;
		} elsif ($element =~ m/^(\d+)\-(\d+)$/) {
			$1 >= 1 and $1 <= $lastchr or pod2usage ("Error in argument: invalid --region argument that contains chromosome $1 (valid chromosomes are 1-$lastchr,X,Y,XY,M");
			$2 >= 1 and $2 <= $lastchr or pod2usage ("Error in argument: invalid --region argument that contains chromosome $2 (valid chromosomes are 1-$lastchr,X,Y,XY,M");
			$1 < $2 or pod2usage ("Error in argument: chromosome designation of --region argument '$element' (in $region) is invalid");
			for ($1 .. $2) {
				$chr_include{$_} = 1;
			}
		} elsif ($element =~ m/^(X|Y|XY|M)$/) {
			$chr_include{$1} = 1;
		} else {
			pod2usage ("Error in argument: valid --region argument should contain chromosome names (valid chromosomes are 1-$lastchr,X,Y,XY,M) but you specified $region");
		}
	}

	my (@header, $name_index, $chr_index, $pos_index, $logr_index, $b_index);
	if ($directory) {
		require File::Spec;
		$directory =~ s#\\\/$##;
		$inputfile = File::Spec->catfile ($directory, $inputfile);
	}
	if ($inputfile =~ m/^`(.+)`$/) {
		open (SIG, "$1 |") or confess "\nERROR: cannot read from system command: $inputfile: $!\n";
	} else {
		open (SIG, $inputfile) or confess "\nERROR: cannot read from inputfile $inputfile: $!\n";
	}

	$_ = <SIG>;
	defined ($_) or die "\nERROR: NOTHING is found in inputfile $inputfile. Please check the file before proceeding\n";
	s/[\r\n]+$//;
	@header = split (/\t/, $_);
	my ($sample_count, $sample_name, $found_sample);
	for my $i (0 .. @header-1) {
		if ($header[$i] eq 'Name' or $header[$i] eq 'SNP' or $header[$i] eq 'SNP ID' or $header[$i] eq 'ProbeID') {		#some files have "SNP ID" or "SNP"
			$name_index = $i;
		} elsif ($header[$i] eq 'Chr' or $header[$i] eq 'Chromosome') {
			$chr_index = $i;
		} elsif ($header[$i] eq 'Position') {
			$pos_index = $i;
		} elsif ($header[$i] =~ m/(.*)\.?Log R Ratio$/ or $header[$i] =~ m/(.*)\.?LRR$/) {
			if ($found_sample) {
				$sample_name eq $1 or confess "\nERROR: previously read B Allele Freq for sample $sample_name but now encountered sample $1\n";
				$logr_index = $i;
				last;
			} else {
				$sample_count++;
				if ($sample_count == 1) {		#previous version use --sample_index here
					$logr_index = $i;
					$sample_name = $1;
					$found_sample = 1;
				}
			}
		} elsif ($header[$i] =~ m/(.*)\.?B Allele Freq(uency)?$/ or $header[$i] =~ m/(.*)\.?BAF$/) {
			if ($found_sample) {
				$sample_name eq $1 or confess "\nERROR: previously read Log R Ratio for sample $sample_name but now encountered sample $1\n";
				$b_index = $i;
				last;
			} else {
				$sample_count++;
				if ($sample_count == 1) {		#previous version use --sample_index here
					$b_index = $i;
					$sample_name = $1;
					$found_sample = 1;
				}
			}
		}
	}
	defined $name_index and defined $logr_index and defined $b_index or die "\nERROR: inputfile $inputfile do not contain Name, *.B Allele Freq and *.Log R Ratio column header (@header)";

	if ($coordinate_from_input) {
		defined $chr_index and defined $pos_index or confess "\nERROR: signal file $inputfile does not contain Chr and Position information in header line (but you specified --coordinate_from_input argument)\n";
	}

#read the logr values for specified individual
	my (%siginfo, %sigdesc, @alllrr, @xlrr, @ylrr, @chr_baf_drift);
	my ($skipped_line, $record_count, $cn_count, $line_count, $nan_count) = qw/0 0 0 0 0/;
	while (<SIG>) {
		s/[\r\n]+$//;
		$line_count++;
		my @record = split (/\t/, $_);
		my ($name, $logr, $baf) = @record[$name_index, $logr_index, $b_index];
		my ($chr, $pos, $pfb);

#in non-English version of BeadStudio, sometimes weird characters may appear in the signal intensity file: they can be ignored
		if (not defined $name) {
			print STDERR "WARNING: Unable to read Name information from line ${\($line_count+1)} in signal file $inputfile: <$_>\n" and next;
		}
		if (not defined $logr or not defined $baf) {
			print STDERR "WARNING: Unable to read LRR and BAF information from line ${\($line_count+1)} in signal file $inputfile: <$_>\n";
			$logr = 0;
			$baf = 0;
			$nan_count++;
		}

		if (not exists $pfbinfo->{$name}) {
			$skipped_line++;				#no PFB information is provided for this probe
				next;
		}

		if ($coordinate_from_input) {		#use chromosome coordinates from the input file (rather than PFB file)
			($chr, $pos) = ($record[$chr_index], $record[$pos_index]);
			defined $chr and defined $pos or confess "\nERROR: Unable to read Chr and Position information from line ${\($line_count+1)} in signal file $inputfile: <$_>\n";
			my @temp = split (/\t/, $pfbinfo->{$name});
			@temp == 3 or die "ERROR: Fatal mistake encountered when handling PFB values when --coord_from_input is specified. Contact PennCNV authors for solutions\n";
			$pfb = $temp[2];
		} else {
			($chr, $pos, $pfb) = split (/\t/, $pfbinfo->{$name});
		}

		if ($logr =~ m/^NaN/ or $baf =~ m/^NaN/) {		#assign the logr and B values to zero, when there are annotated as 'Not a Number' (in non-English version of BeadStudio, some additional strings may be appended after NaN
				$logr = 0;
				$baf = 0;
				$nan_count++;
				}
				$logr eq '-Infinity' and $logr = -5;			#this usually means homozygous deletion (zero copy number)
				$logr eq '-inf' and $logr = -5;				#same as above: this usually means homozygous deletion (zero copy number)
				if ($logr =~ m/[^\d\.eE-]/ or $baf =~ m/[^\d\.eE-]/) {			#current (2008Sep04) to handle "weird character" problems caused by non-English version of BeadStudio
				print STDERR "WARNING: Unrecognizable LRR or BAF values are treated as zero: LRR=$logr BAF=$baf\n";
				$logr = 0;
				$baf = 0;
				$nan_count++;
				}

				if ($name =~ m/^CN/ or $name =~ m/^cnv/) {		#this step force the BAF for CN marker (non-polymorphic marker) to be 2 (normally it should be in 0-1 range for a SNP marker)
				$baf = 2;
				$cn_count++;
				} elsif ($pfb > 1) {					#if the PFB value is annotated as >1, it is also treated as a NP marker (this is used for custom array design)
				$baf = 2;
				$cn_count++;
				}
				$record_count++;

				if (not $siginfo{$chr}) {
					push @{$siginfo{$chr}->{name}}, 0;
					push @{$siginfo{$chr}->{pos}}, [0, 0];
					push @{$siginfo{$chr}->{lrr}}, 0;
					push @{$siginfo{$chr}->{baf}}, 0;
					push @{$siginfo{$chr}->{pfb}}, 0;
				}
				push @{$siginfo{$chr}->{name}}, $name;
				push @{$siginfo{$chr}->{pos}}, [$pos, @{$siginfo{$chr}->{name}}-1];
				push @{$siginfo{$chr}->{lrr}}, $logr;
				push @{$siginfo{$chr}->{baf}}, $baf;
				push @{$siginfo{$chr}->{pfb}}, $pfb;
	}
	close (SIG);

#make sure that the SNPs in the inputfile are sorted based on positions from the PFB file
	for my $nextchr (keys %siginfo) {
		@{$siginfo{$nextchr}->{pos}} = sort {$a->[0]<=>$b->[0]} @{$siginfo{$nextchr}->{pos}};
		my @index = map {$_->[1]} @{$siginfo{$nextchr}->{pos}};
		@{$siginfo{$nextchr}->{pos}} = map {$_->[0]} @{$siginfo{$nextchr}->{pos}};

		@{$siginfo{$nextchr}->{name}} = @{$siginfo{$nextchr}->{name}}[@index];
		@{$siginfo{$nextchr}->{lrr}} = @{$siginfo{$nextchr}->{lrr}}[@index];
		@{$siginfo{$nextchr}->{baf}} = @{$siginfo{$nextchr}->{baf}}[@index];
		@{$siginfo{$nextchr}->{pfb}} = @{$siginfo{$nextchr}->{pfb}}[@index];

		if ($nextchr =~ m/^\d+$/ and $nextchr>=1 and $nextchr<=$lastchr) {	#baf_drift is a measure of BAF in "unlikely regions" (except when 4-copy duplication is present): high baf_drift indicate sample quality issue
			my @baf_drift = grep {$_>0.2 and $_<0.25 or $_>0.75 and $_<0.8} @{$siginfo{$nextchr}->{baf}};
			push @chr_baf_drift, @baf_drift/@{$siginfo{$nextchr}->{baf}};
		}
	}

	print STDERR " Done with $record_count records in ${\(scalar keys %siginfo)} chromosomes";
	if ($skipped_line) {
		print STDERR " ($skipped_line records are discarded due to lack of PFB information for the markers)\n";
	} else {
		print STDERR "\n";
	}
	$verbose and $cn_count and print STDERR "NOTICE: Total of $cn_count non-polymorphic probes detected in the $record_count records from $inputfile\n";

#handle situations where nothing was in the file
	$record_count or print STDERR "ERROR: NO SIGNAL DATA FOUND IN INPUTFILE $inputfile\n" and return ({}, {});

	$sigdesc{num_record} = $record_count;					#total number of records (not including the leading zero in each array for each chromosome)
		$sigdesc{cn_count} = $cn_count;						#total number of CN probes
		$sigdesc{nocall_rate} = $nan_count/$record_count;

	($sigdesc{wf}, $sigdesc{cc}) = calWF (\%siginfo);
	$sigdesc{gcwf} = $sigdesc{cc} eq 'NA' ? 'NA' : ($sigdesc{wf} * abs ($sigdesc{cc}));

#adjust genomic wave using GC model file
	if ($gcmodel) {
		adjustLRRByGCModel (\%siginfo, \%sigdesc, $gcmodel);
	}

	for my $nextchr (keys %siginfo) {
		if ($nextchr =~ m/^\d+$/) {
			push @alllrr, grep {$_>-2 and $_<2} @{$siginfo{$nextchr}->{lrr}};	#must use PUSH here, for each chromsome
		} elsif ($nextchr eq 'X') {
			@xlrr = grep {$_>-2 and $_<2} @{$siginfo{$nextchr}->{lrr}};
		} elsif ($nextchr eq 'Y') {
			@ylrr = grep {$_>-2 and $_<2} @{$siginfo{$nextchr}->{lrr}};
		}
	}
	if (@alllrr <= 1) {							#for example, when the signal file does not contain any autosomes
		@sigdesc{'lrr_mean', 'lrr_median', 'lrr_sd'} = qw/NA NA NA/;
	} else {
		$sigdesc{lrr_mean} = mean (\@alllrr);				#lrr mean, sd and median are for LRR values in (-2, 2) region only!
			$sigdesc{lrr_median} = median (\@alllrr);
		$sigdesc{lrr_sd} = sd (\@alllrr);
	}
	if (@xlrr <= 1) {
		@sigdesc{'baf_xhet', 'lrr_xmean', 'lrr_xmedian', 'lrr_xsd'} = qw/NA NA NA NA/;
	} else {
		@sigdesc{'lrr_xmean', 'lrr_xmedian', 'lrr_xsd'} = (mean (\@xlrr), median (\@xlrr), sd (\@xlrr));
	}

	if (@ylrr <= 1) {
		@sigdesc{'baf_yhet', 'lrr_ymean', 'lrr_ymedian', 'lrr_ysd'} = qw/NA NA NA NA/;
	} else {
		@sigdesc{'lrr_ymean', 'lrr_ymedian', 'lrr_ysd'} = (mean (\@ylrr), median (\@ylrr), sd (\@ylrr));
	}

	@alllrr = ();
	@xlrr = ();
	@ylrr = ();
	for my $nextchr (keys %siginfo) {
		if ($nextchr =~ m/^\d+$/) {
			push @alllrr, grep {$_>0.25 and $_<0.75} @{$siginfo{$nextchr}->{baf}};	#must use PUSH here, for each chromosome
		} elsif ($nextchr eq 'X') {
			@xlrr = @{$siginfo{$nextchr}->{baf}};		#unfilterd BAF values here
		} elsif ($nextchr eq 'Y') {
			@ylrr = @{$siginfo{$nextchr}->{baf}};		#unfilterd BAF values here
		}
	}

	if (@alllrr <= 1) {
		@sigdesc{'baf_mean', 'baf_median', 'baf_sd'} = qw/NA NA NA/;
	} else {
		$sigdesc{baf_mean} = mean (\@alllrr);				#baf mean, sd and median are for BAF values in (0.25,0.75) region only!
			$sigdesc{baf_median} = median (\@alllrr);
		$sigdesc{baf_sd} = sd (\@alllrr);
		if (@chr_baf_drift) {
			$sigdesc{baf_drift} = median (\@chr_baf_drift);		#baf_drift measure the random drift of BAF values to unlikely regions (we use median here, to eliminate effects of large baf_drift caused by one chromosome having large 4-copy duplication)
		} else {
			$sigdesc{baf_drift} = 0;				#for chrX, there is no baf_drift value
		}
	}

	if (@xlrr > 1) {								#occasionally the file do not contain SNPs for X chromosome
		if (scalar (grep {$_>=0 and $_<=1} @xlrr)) {
			$sigdesc{baf_xhet} = scalar (grep {$_>0.25 and $_<0.75} @xlrr)/scalar (grep {$_>=0 and $_<=1} @xlrr);	#the heterozygosity fraction (this is especially useful for chrX, as this value can be used to infer sex)
		} else {
			$sigdesc{baf_xhet} = 'NA';
		}
	}

	if (@ylrr > 1) {								
		if (scalar (grep {$_>=0 and $_<=1} @ylrr)) {
			$sigdesc{baf_yhet} = scalar (grep {$_>0.25 and $_<0.75} @ylrr)/scalar (grep {$_>=0 and $_<=1} @ylrr); # fraction of heterozygous SNPs in chrY. this can not be used to predict the gender
		} else {
			$sigdesc{baf_yhet} = 'NA';  # no SNPs on chrY
		}
	}


	my @dropchr;
	for my $nextchr (keys %siginfo) {
		if (not $chr_include{$nextchr}) {
			delete $siginfo{$nextchr};
			push @dropchr, $nextchr;
		}
	}
	@dropchr and print STDERR "NOTICE: Data from chromosome ", join (",", sort {sortChr($a, $b)} @dropchr), " will not be used in analysis\n";
	return (\%siginfo, \%sigdesc);
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

sub adjustLRRByGCModel {
	my ($siginfo, $sigdesc, $model) = @_;
	my (@x, @y);
	my ($b1, $b2, $f, $p) = qw/0 0 0 0/;
	my ($prepos, $distance) = (0, 1_000_000);

	my ($wf, $cc, $gcwf) = ($sigdesc->{wf}, $sigdesc->{cc}, $sigdesc->{gcwf});

	if ($cc eq 'NA') {
		print STDERR "WARNING: Unable to adjust LRR values by GC model due to lack of GCWF measure\n";
		return;
	}

	for my $nextchr (keys %$siginfo) {
		$nextchr =~ m/^\d+$/ or next;			#model training on autosome only
			$prepos = 0;					#reset prepos for the new chromosome
			my ($name, $pos, $lrr) = ($siginfo->{$nextchr}{name}, $siginfo->{$nextchr}{pos}, $siginfo->{$nextchr}{lrr});
		for my $i (0 .. @$name-1) {
			$model->{$name->[$i]} or next;		#this SNP is not in GC model
				if ($pos->[$i] - $prepos > $distance) {
					my ($x, $y) = ($model->{$name->[$i]}, $lrr->[$i]);
					if ($x > 15 and $x < 80 and $y > -1 and $y < 1) {
						push @x, $x;
						push @y, $y;
					}
					$prepos = $pos->[$i];
				}
		}
	}
	($b1, $b2) = reg_linear (\@x, \@y);
	$verbose and print STDERR "NOTICE: Collecting ${\(scalar @x)} autosome SNPs for GC-based regression ... GC model regression parameter: a=$b1 b=$b2\n";

	my ($count_skipped_marker, @alllrr, @xlrr, @ylrr);
	for my $nextchr (keys %$siginfo) {
		my ($name, $lrr) = ($siginfo->{$nextchr}{name}, $siginfo->{$nextchr}{lrr});
		for my $i (0 .. @$name-1) {
			if ($model->{$name->[$i]}) {
				$lrr->[$i] -= ($b1 + $b2 * $model->{$name->[$i]});
				if ($nextchr =~ m/^\d+$/) {
					push @alllrr, $lrr->[$i];
				} elsif ($nextchr eq 'X') {
					push @xlrr, $lrr->[$i];
				} elsif ($nextchr eq 'Y') {
					push @ylrr, $lrr->[$i];
				}
			} else {
				$count_skipped_marker++;
			}
		}
	}
	($sigdesc->{wf}, $sigdesc->{cc}) = calWF ($siginfo);
	$sigdesc->{gcwf} = $sigdesc->{wf} * abs ($sigdesc->{cc});

	@alllrr = grep {$_>-2 and $_<2} @alllrr;
	@xlrr = grep {$_>-2 and $_<2} @xlrr;
	if (@alllrr <= 1) {
		@{$sigdesc->{'lrr_mean', 'lrr_median', 'lrr_sd'}} = qw/NA NA NA/;
	} else {
		@{$sigdesc->{'lrr_mean', 'lrr_median', 'lrr_sd'}} = (sprintf("%.4f", mean (\@alllrr)), sprintf("%.4f", median (\@alllrr)), sprintf("%.4f", sd (\@alllrr)));
	}

	@xlrr = grep {$_>-2 and $_<2} @xlrr;
	if (@xlrr <= 1) {
		@{$sigdesc->{'lrr_xmean', 'lrr_xmedian', 'lrr_xsd'}} = qw/NA NA NA/;
	} else {
		@{$sigdesc->{'lrr_xmean', 'lrr_xmedian', 'lrr_xsd'}} = (sprintf("%.4f", mean (\@xlrr)), sprintf("%.4f", median (\@xlrr)), sprintf("%.4f", sd (\@xlrr)));
	}
	$verbose and $count_skipped_marker and print STDERR "NOTICE: A total of $count_skipped_marker markers are not adjusted by GC model due to lack of model parameters\n";
	print STDERR "NOTICE: Adjusting LRR by GC model: WF changes from ", sprintf ("%.4f", $wf), " to ", sprintf ("%.4f", $sigdesc->{wf}), ", GCWF changes from ", sprintf ("%.4f", $gcwf), " to ", sprintf ("%.4f", $sigdesc->{gcwf}), "\n";
}

sub newreadGCModel {
	my ($gcmodelfile, $pfbinfo) = @_;
	my (%snp_gc);
	open (GC, $gcmodelfile) or confess "\nERROR: cannot read from gcmodelfile $gcmodelfile: $!\n";
	while (<GC>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		@record >= 4 or confess "\nERROR: invalid record found in gcmodelfile (four tab-delimited fields expected): <$_>\n";
		if ($pfbinfo->{$record[0]}) {
			$snp_gc{$record[0]} = $record[3];
		}
	}
	close (GC);
	return (\%snp_gc);
}

sub identifyCNVRegion {
	my ($cnvfile, $ref_inputfile, $pfbinfo) = @_;

	my $snp_pos;


	my %inputsignalfile = map {$_, 1} @$ref_inputfile;
	my (@region, %region);				#@region contains the region (3 elements in array, representing chr, start, end) %region stores files that have this region
		open (TRIOCNV, $cnvfile) or confess "\nERROR: cannot read from cnvfile $cnvfile: $!";
	while (<TRIOCNV>) {
		m/^chr(\w+):(\d+)-(\d+)\s+numsnp=(\d+)\s+length=([\d\,]+)\s+state(\d+),cn=(\d+)\s+(.+?)\s+startsnp=(\S+)\s+endsnp=(\S+)/ or die "Error: invalid record found in cnvfile $cnvfile: <$_>\n";
		my ($curchr, $curstate, $curfile, $curstartsnp, $curendsnp) = ($1, $6, $8, $9, $10);

		$inputsignalfile{$curfile} or next;		#only consider regions that are related to the inputfiles
			$curstartsnp eq $curendsnp and next;		#ignore CNV that contain only one SNP
			if ($chrx) {
				$curchr eq 'X' or next;
			} else {
				$curchr eq 'X' and next;
			}
		if ($chry) {
			$curchr eq 'Y' or next;
		} else {
			$curchr eq 'Y' and next;
		}
		$region{"$curchr:$curstartsnp-$curendsnp"} and $region{"$curchr:$curstartsnp-$curendsnp"} .= "$curfile(state$curstate) " and next;		#this region is already processed (for example, when it occur in both father and son)
			$region{"$curchr:$curstartsnp-$curendsnp"} .= "$curfile(state$curstate) ";		#key=CNV value=filename(state)
			push @region, [$curchr, $curstartsnp, $curendsnp];					#array containing CNV regions

			my @temp = split (/\t/, $pfbinfo->{$curstartsnp});
		$snp_pos->{$curstartsnp} = $temp[1];
		@temp = split (/\t/, $pfbinfo->{$curendsnp});
		$snp_pos->{$curendsnp} = $temp[1];
	}
	close (TRIOCNV);
	print STDERR "NOTICE: Total of ${\(scalar @region)} CNV regions will be validated by family information (@$ref_inputfile)\n";
	@region or print STDERR "WARNING: NO CNV CAN BE CALLED BY TRIO-BASED ALGORITHM FOR THE TRIO: @$ref_inputfile\n";
	@region or return;								#no region will be used for family validation (none of the parent or offspring has CNV calls)

#the goal of the following paragraph is to combine overlapping regions into one single region
#for example, the father may have a region of 100-200 as CNV, but the child may have a region of 105-205 as CNV region, so which region is more reasonable?
#to solve this problem, the best approach is to directly test the four regions via Viterbi algorithm: 100-104, 205-200, 201-205 and then see which one has the highest prob
#after the following paragraph:
#@region contains all regions, or arrays of [chr, start1, end1, end2, end3], etc.
#@subregion contains all subregions, or arrays of  [chr, start, end]
		@region = sort {sortChr ($a->[0], $b->[0]) || $snp_pos->{$a->[1]}<=>$snp_pos->{$b->[1]} || $snp_pos->{$a->[2]} <=> $snp_pos->{$b->[2]}} @region;
	my (@newregion, @subregion);
	my ($prechr, $prestart, $preend) = (0);
	my (%startname, %endname);			#key is snp name, value is the count that it is start or it is end

		for my $nextregion (@region) {
			my ($curchr, $curstart, $curend) = @$nextregion;
			if ($curchr ne $prechr or $snp_pos->{$curstart} > $snp_pos->{$preend}) {
				push @newregion, $nextregion;
			} else {
				push @{$newregion[$#newregion]}, $curstart, $curend;
			}
			$startname{$curstart}++;
			$endname{$curend}++;
			($prechr, $prestart, $preend) = ($curchr, $curstart, ((!$preend || $snp_pos->{$curend} > $snp_pos->{$preend} || $curchr ne $prechr) ? $curend : $preend));
		}
	@region = ();
	for my $nextregion (@newregion) {
		my ($curchr, @curpos) = @$nextregion;
		my %seen;
		@curpos = grep {!$seen{$_}++} @curpos;
		@curpos = sort {$snp_pos->{$a} <=> $snp_pos->{$b}} @curpos;
		push @region, [$curchr, @curpos];
		for my $i (0 .. @curpos-2) {
			push @subregion, [$curchr, $curpos[$i], $curpos[$i+1]];
		}
		$verbose and print STDERR "NOTICE: adding region for analysis: chr$curchr, @curpos\n";
	}
	@region = sort {sortChr ($a->[0], $b->[0]) || $snp_pos->{$a->[1]} <=> $snp_pos->{$b->[1]}} @region;
	@subregion = sort {sortChr ($a->[0], $b->[0]) || $snp_pos->{$a->[1]} <=> $snp_pos->{$b->[1]}} @subregion;
	return (\@subregion, \@region, \%region, \%startname, \%endname);
}

#read the signal files for the family and retrieve the signal information for the candidate CNV regions
sub retrieveRegionSignal {
	my ($ref_inputfile, $file_sex, $hmmfile, $subregion, $pfbinfo, $gcmodel, $directory, $hash_startname, $hash_endname) = @_;

#the following paragraph read the three signal intensity files and calculate the probability for each subregion
#the hash region_* all have key as subregion (chr,start,end)
	my (%region_alllogprob, %region_numsnp, %region_cnvlength, %region_chr, %region_name_start, %region_name_end, %region_pos_start, %region_pos_end, %region_snpdist);
	my (@sample_sex, $sample_sex);

	my $hmm = readHMMFile ($hmmfile);
	my ($siginfo, $sigdesc);			#the siginfo for the last individual will be returned by the program (for calculating CNV length when printing CNV calls)

		for my $i (0 .. @$ref_inputfile-1) {
			my $inputfile = $ref_inputfile->[$i];
			my $sample_sex = $file_sex->{$inputfile} || undef;
			if (defined $sample_sex) {
				if ($i == 0) {
					$sample_sex eq 'male' or print STDERR "WARNING: sample_sex is annotated as $sample_sex in $sexfile but it is assumed as father (male) in analysis\n" and $sample_sex = 'male';
				} elsif ($i == 1) {
					$sample_sex eq 'female' or print STDERR "WARNING: sample_sex is annotated as $sample_sex in $sexfile but it is assumed as mother (female) in analysis\n" and $sample_sex = 'female';
				}
			}

			($siginfo, $sigdesc) = readLRRBAF ($inputfile, $region, $pfbinfo, $gcmodel, $directory);

#if the file processing fails and no signal data is read from file, move on to the next file
			%$siginfo or confess "\nERROR: Unable to read signal intensity information for chr$region from $inputfile for family-based CNV calling. Program halted.\n";

			$sample_sex = QCSignal ($inputfile, $siginfo, $sigdesc, $hmm, $sample_sex);
			push @sample_sex, $sample_sex;

#read HMM files
			my $hmm_model = khmm::ReadCHMM ($hmmfile);

#Now a mandatory step to handle low-quality samples and reduce false-positve calls, through matching the SD measure for sample and model
			if ($sdadjust) {
				if ($chrx) {
					khmm::adjustBSD ($hmm_model, $sigdesc->{lrr_xsd});
				} elsif ($chry) {
					khmm::adjustBSD ($hmm_model, $sigdesc->{lrr_ysd});
				} else {
					khmm::adjustBSD ($hmm_model, $sigdesc->{lrr_sd});
				}
			}

			my ($pre_chr, $pre_snpdist);
			for my $nextsubregion (@$subregion) {
				$verbose and print STDERR "NOTICE: Processing subregion @$nextsubregion\n";
				my $subregionstring = join (",", @$nextsubregion);
				my ($region_name, $region_chr, $region_pos, $region_logr, $region_baf, $region_pfb, $region_snpdist) = getRegionInfo ($siginfo, $pfbinfo, $nextsubregion, $hash_startname, $hash_endname);

				my ($logprob, $scale) = (0);
				my @logprob;
				$verbose and print STDERR "NOTICE: processing @$nextsubregion containing ", scalar (@$region_name), " SNPs, where start=", $hash_startname->{$nextsubregion->[1]}||0, " and end=", $hash_endname->{$nextsubregion->[2]}||0, "\n";
				unshift @$region_logr, 0; unshift @$region_baf, 0; unshift @$region_pfb, 0; unshift @$region_snpdist, 0;
				for my $nextstate (1..6) {
					khmm::GetStateProb_CHMM ($hmm_model, @$region_logr-1, $region_logr, $region_baf, $region_pfb, $region_snpdist, \$logprob, $nextstate);
					push @logprob, $logprob;
					$scale or $scale = $logprob;
					$scale < $logprob and $scale = $logprob;
				}
				@logprob = map {$_-$scale} @logprob;		#make the highest probability as zero arbitrarily
					unshift @logprob, 0;				#fill in the state0 (dummy state)

					$verbose and print STDERR "logprob for $inputfile in region ($region_name) in 6 states are @logprob\n";
				push @{$region_alllogprob{$subregionstring}}, \@logprob;
				$region_numsnp{$subregionstring} = @$region_logr-1;
				$region_cnvlength{$subregionstring} = $region_pos->[@$region_pos-1]-$region_pos->[0]+1;
				$region_chr{$subregionstring} = $region_chr->[0];
				$region_name_start{$subregionstring} = $region_name->[0];
				$region_name_end{$subregionstring} = $region_name->[@$region_name-1];
				$region_pos_start{$subregionstring} = $region_pos->[0];
				$region_pos_end{$subregionstring} = $region_pos->[@$region_pos-1];
				if ($pre_chr and $pre_chr eq $region_chr->[0]) {			#encountered new chromosome
					$region_snpdist{$subregionstring} = $region_pos->[0]-$pre_snpdist;
					$verbose and print STDERR "NOTICE: region $subregionstring snpdist=$region_snpdist{$subregionstring} ($region_pos->[0]-$pre_snpdist) curend=$region_pos_end{$subregionstring}\n";
				}
				$pre_snpdist = $region_pos_end{$subregionstring};
				$pre_chr = $region_chr->[0];
			}
		}

	return ($siginfo, \@sample_sex, \%region_alllogprob, \%region_numsnp, \%region_cnvlength, \%region_chr, \%region_name_start, \%region_name_end, \%region_pos_start, \%region_pos_end, \%region_snpdist);
}

#validate CNV calls using father-mother-offspring trio informatoin and reconcile boundary discordance using family information.
sub newtestTrioCNVFile {
	my ($ref_inputfile, $hmmfile, $pfbfile, $cnvfile, $sexfile, $gcmodelfile, $directory) = @_;
#my ($name, $chr, $pos, $logr, $baf, $pfb, $snpdist);
	my ($pfbinfo) = newreadPFB ($pfbfile);
	my $file_sex;
	$sexfile and $file_sex = readSexFile ($sexfile);
	my $hmm = readHMMFile ($hmmfile);
	my $gcmodel;

	defined $gcmodelfile and $gcmodel = newreadGCModel ($gcmodelfile, $pfbinfo);	#read the gc model for SNPs defined in snp_chr

		my ($array_subregion, $array_region, $hash_region, $hash_startname, $hash_endname) = identifyCNVRegion ($cnvfile, $ref_inputfile, $pfbinfo);
	$array_subregion or return;		#no CNV region can be identified for posterior calling algorithm

		my @subregion = @$array_subregion;
	my @region = @$array_region;
	my %region = %$hash_region;
	my %startname = %$hash_startname;
	my %endname = %$hash_endname;


	my ($siginfo, $sample_sex, $region_alllogprob, $region_numsnp, $region_cnvlength, $region_chr, $region_name_start, $region_name_end, $region_pos_start, $region_pos_end, $region_snpdist) 
		= retrieveRegionSignal ($ref_inputfile, $file_sex, $hmmfile, $array_subregion, $pfbinfo, $gcmodel, $directory, $hash_startname, $hash_endname);


#now apply pedigree information to validate CNV calls and fine-map CNV boundaries
#first define the CNV inheritance matrix, and read the state transition matrix

	my $transition = readTransitionMatrix ($hmmfile);
	my ($e, %oprior) = ($denovo_rate);
	if ($chrx and $sample_sex->[2] eq 'male') {
		%oprior = (111=>1-$e, 112=>0.25*$e, 113=>0.25*$e, 115=>0.25*$e, 116=>0.25*$e, 121=>0.5*(1-$e), 122=>0.5*(1-$e), 123=>0.3333*$e, 125=>0.3333*$e, 126=>0.3333*$e, 131=>0.25*$e, 132=>1-$e, 133=>0.25*$e, 135=>0.25*$e, 136=>0.25*$e, 151=>0.3333*$e, 152=>0.5*(1-$e), 153=>0.5*(1-$e), 155=>0.3333*$e, 156=>0.3333*$e, 161=>0.5*$e, 162=>0.25*(1-$e), 163=>0.5*(1-$e), 165=>0.25*(1-$e), 166=>0.5*$e, 211=>0.25*$e, 212=>1-$e, 213=>0.25*$e, 215=>0.25*$e, 216=>0.25*$e, 221=>0.3333*$e, 222=>0.5*(1-$e), 223=>0.5*(1-$e), 225=>0.3333*$e, 226=>0.3333*$e, 231=>0.25*$e, 232=>0.25*$e, 233=>1-$e, 235=>0.25*$e, 236=>0.25*$e, 251=>0.3333*$e, 252=>0.3333*$e, 253=>0.5*(1-$e), 255=>0.5*(1-$e), 256=>0.3333*$e, 261=>0.5*$e, 262=>0.5*$e, 263=>0.25*(1-$e), 265=>0.5*(1-$e), 266=>0.25*(1-$e), 311=>0.25*$e, 312=>0.25*$e, 313=>1-$e, 315=>0.25*$e, 316=>0.25*$e, 321=>0.3333*$e, 322=>0.3333*$e, 323=>0.5*(1-$e), 325=>0.5*(1-$e), 326=>0.3333*$e, 331=>0.25*$e, 332=>0.25*$e, 333=>0.25*$e, 335=>1-$e, 336=>0.25*$e, 351=>0.3333*$e, 352=>0.3333*$e, 353=>0.3333*$e, 355=>0.5*(1-$e), 356=>0.5*(1-$e), 361=>0.3333*$e, 362=>0.3333*$e, 363=>0.3333*$e, 365=>0.25*(1-$e), 366=>0.75*(1-$e), 511=>0.25*$e, 512=>0.25*$e, 513=>0.25*$e, 515=>1-$e, 516=>0.25*$e, 521=>0.3333*$e, 522=>0.3333*$e, 523=>0.3333*$e, 525=>0.5*(1-$e), 526=>0.5*(1-$e), 531=>0.25*$e, 532=>0.25*$e, 533=>0.25*$e, 535=>0.25*$e, 536=>1-$e, 551=>0.25*$e, 552=>0.25*$e, 553=>0.25*$e, 555=>0.25*$e, 556=>1-$e, 561=>0.25*$e, 562=>0.25*$e, 563=>0.25*$e, 565=>0.25*$e, 566=>1-$e, 611=>0.25*$e, 612=>0.25*$e, 613=>0.25*$e, 615=>0.25*$e, 616=>1-$e, 621=>0.25*$e, 622=>0.25*$e, 623=>0.25*$e, 625=>0.25*$e, 626=>1-$e, 631=>0.25*$e, 632=>0.25*$e, 633=>0.25*$e, 635=>0.25*$e, 636=>1-$e, 651=>0.25*$e, 652=>0.25*$e, 653=>0.25*$e, 655=>0.25*$e, 656=>1-$e, 661=>0.25*$e, 662=>0.25*$e, 663=>0.25*$e, 665=>0.25*$e, 666=>1-$e);
	} elsif ($chrx and $sample_sex->[2] eq 'female') {
		%oprior = (111=>1-$e, 112=>0.25*$e, 113=>0.25*$e, 115=>0.25*$e, 116=>0.25*$e, 121=>0.5*(1-$e), 122=>0.5*(1-$e), 123=>0.3333*$e, 125=>0.3333*$e, 126=>0.3333*$e, 131=>0.25*$e, 132=>1-$e, 133=>0.25*$e, 135=>0.25*$e, 136=>0.25*$e, 151=>0.3333*$e, 152=>0.5*(1-$e), 153=>0.5*(1-$e), 155=>0.3333*$e, 156=>0.3333*$e, 161=>0.5*$e, 162=>0.25*(1-$e), 163=>0.5*(1-$e), 165=>0.25*(1-$e), 166=>0.5*$e, 211=>1-$e, 212=>0.25*$e, 213=>0.25*$e, 215=>0.25*$e, 216=>0.25*$e, 221=>0.5*(1-$e), 222=>0.5*(1-$e), 223=>0.3333*$e, 225=>0.3333*$e, 226=>0.3333*$e, 231=>0.25*$e, 232=>1-$e, 233=>0.25*$e, 235=>0.25*$e, 236=>0.25*$e, 251=>0.3333*$e, 252=>0.5*(1-$e), 253=>0.5*(1-$e), 255=>0.3333*$e, 256=>0.3333*$e, 261=>0.5*$e, 262=>0.25*(1-$e), 263=>0.5*(1-$e), 265=>0.25*(1-$e), 266=>0.5*$e, 311=>1-$e, 312=>0.25*$e, 313=>0.25*$e, 315=>0.25*$e, 316=>0.25*$e, 321=>0.5*(1-$e), 322=>0.5*(1-$e), 323=>0.3333*$e, 325=>0.3333*$e, 326=>0.3333*$e, 331=>0.25*$e, 332=>1-$e, 333=>0.25*$e, 335=>0.25*$e, 336=>0.25*$e, 351=>0.3333*$e, 352=>0.5*(1-$e), 353=>0.5*(1-$e), 355=>0.3333*$e, 356=>0.3333*$e, 361=>0.5*$e, 362=>0.25*(1-$e), 363=>0.5*(1-$e), 365=>0.25*(1-$e), 366=>0.5*$e, 511=>1-$e, 512=>0.25*$e, 513=>0.25*$e, 515=>0.25*$e, 516=>0.25*$e, 521=>0.5*(1-$e), 522=>0.5*(1-$e), 523=>0.3333*$e, 525=>0.3333*$e, 526=>0.3333*$e, 531=>0.25*$e, 532=>1-$e, 533=>0.25*$e, 535=>0.25*$e, 536=>0.25*$e, 551=>0.3333*$e, 552=>0.5*(1-$e), 553=>0.5*(1-$e), 555=>0.3333*$e, 556=>0.3333*$e, 561=>0.5*$e, 562=>0.25*(1-$e), 563=>0.5*(1-$e), 565=>0.25*(1-$e), 566=>0.5*$e, 611=>1-$e, 612=>0.25*$e, 613=>0.25*$e, 615=>0.25*$e, 616=>0.25*$e, 621=>0.5*(1-$e), 622=>0.5*(1-$e), 623=>0.3333*$e, 625=>0.3333*$e, 626=>0.3333*$e, 631=>0.25*$e, 632=>1-$e, 633=>0.25*$e, 635=>0.25*$e, 636=>0.25*$e, 651=>0.3333*$e, 652=>0.5*(1-$e), 653=>0.5*(1-$e), 655=>0.3333*$e, 656=>0.3333*$e, 661=>0.5*$e, 662=>0.25*(1-$e), 663=>0.5*(1-$e), 665=>0.25*(1-$e), 666=>0.5*$e);
	}elsif ($chry and $sample_sex->[2] eq 'male'){
		%oprior = (111=>1.000-$e, 112=>1.000*$e, 113=>0.500*$e, 115=>0.333*$e, 116=>0.250*$e, 121=>1.000-$e, 122=>1.000*$e, 123=>0.500*$e, 125=>0.333*$e, 126=>0.250*$e, 131=>1.000-$e, 132=>1.000*$e, 133=>0.500*$e, 135=>0.333*$e, 136=>0.250*$e, 151=>1.000-$e, 152=>1.000*$e, 153=>0.500*$e, 155=>0.333*$e, 156=>0.250*$e, 161=>1.000-$e, 162=>1.000*$e, 163=>0.500*$e, 165=>0.333*$e, 166=>0.250*$e, 211=>1.000*$e, 212=>1.000-$e, 213=>1.000*$e, 215=>0.500*$e, 216=>0.333*$e, 221=>1.000*$e, 222=>1.000-$e, 223=>1.000*$e, 225=>0.500*$e, 226=>0.333*$e, 231=>1.000*$e, 232=>1.000-$e, 233=>1.000*$e, 235=>0.500*$e, 236=>0.333*$e, 251=>1.000*$e, 252=>1.000-$e, 253=>1.000*$e, 255=>0.500*$e, 256=>0.333*$e, 261=>1.000*$e, 262=>1.000-$e, 263=>1.000*$e, 265=>0.500*$e, 266=>0.333*$e, 311=>0.500*$e, 312=>1.000*$e, 313=>1.000-$e, 315=>1.000*$e, 316=>0.500*$e, 321=>0.500*$e, 322=>1.000*$e, 323=>1.000-$e, 325=>1.000*$e, 326=>0.500*$e, 331=>0.500*$e, 332=>1.000*$e, 333=>1.000-$e, 335=>1.000*$e, 336=>0.500*$e, 351=>0.500*$e, 352=>1.000*$e, 353=>1.000-$e, 355=>1.000*$e, 356=>0.500*$e, 361=>0.500*$e, 362=>1.000*$e, 363=>1.000-$e, 365=>1.000*$e, 366=>0.500*$e, 511=>0.333*$e, 512=>0.500*$e, 513=>1.000*$e, 515=>1.000-$e, 516=>1.000*$e, 521=>0.333*$e, 522=>0.500*$e, 523=>1.000*$e, 525=>1.000-$e, 526=>1.000*$e, 531=>0.333*$e, 532=>0.500*$e, 533=>1.000*$e, 535=>1.000-$e, 536=>1.000*$e, 551=>0.333*$e, 552=>0.500*$e, 553=>1.000*$e, 555=>1.000-$e, 556=>1.000*$e, 561=>0.333*$e, 562=>0.500*$e, 563=>1.000*$e, 565=>1.000-$e, 566=>1.000*$e, 611=>0.250*$e, 612=>0.333*$e, 613=>0.500*$e, 615=>1.000*$e, 616=>1.000-$e, 621=>0.250*$e, 622=>0.333*$e, 623=>0.500*$e, 625=>1.000*$e, 626=>1.000-$e, 631=>0.250*$e, 632=>0.333*$e, 633=>0.500*$e, 635=>1.000*$e, 636=>1.000-$e, 651=>0.250*$e, 652=>0.333*$e, 653=>0.500*$e, 655=>1.000*$e, 656=>1.000-$e, 661=>0.250*$e, 662=>0.333*$e, 663=>0.500*$e, 665=>1.000*$e, 666=>1.000-$e);
	}elsif ($chry and $sample_sex->[2] eq 'female'){
		%oprior = (111=>$e*$e, 112=>$e*$e, 113=>$e*$e, 115=>$e*$e, 116=>$e*$e, 121=>$e*$e, 122=>$e*$e, 123=>$e*$e, 125=>$e*$e, 126=>$e*$e, 131=>$e*$e, 132=>$e*$e, 133=>$e*$e, 135=>$e*$e, 136=>$e*$e, 151=>$e*$e, 152=>$e*$e, 153=>$e*$e, 155=>$e*$e, 156=>$e*$e, 161=>$e*$e, 162=>$e*$e, 163=>$e*$e, 165=>$e*$e, 166=>$e*$e, 211=>$e*$e, 212=>$e*$e, 213=>$e*$e, 215=>$e*$e, 216=>$e*$e, 221=>$e*$e, 222=>$e*$e, 223=>$e*$e, 225=>$e*$e, 226=>$e*$e, 231=>$e*$e, 232=>$e*$e, 233=>$e*$e, 235=>$e*$e, 236=>$e*$e, 251=>$e*$e, 252=>$e*$e, 253=>$e*$e, 255=>$e*$e, 256=>$e*$e, 261=>$e*$e, 262=>$e*$e, 263=>$e*$e, 265=>$e*$e, 266=>$e*$e, 311=>$e*$e, 312=>$e*$e, 313=>$e*$e, 315=>$e*$e, 316=>$e*$e, 321=>$e*$e, 322=>$e*$e, 323=>$e*$e, 325=>$e*$e, 326=>$e*$e, 331=>$e*$e, 332=>$e*$e, 333=>$e*$e, 335=>$e*$e, 336=>$e*$e, 351=>$e*$e, 352=>$e*$e, 353=>$e*$e, 355=>$e*$e, 356=>$e*$e, 361=>$e*$e, 362=>$e*$e, 363=>$e*$e, 365=>$e*$e, 366=>$e*$e, 511=>$e*$e, 512=>$e*$e, 513=>$e*$e, 515=>$e*$e, 516=>$e*$e, 521=>$e*$e, 522=>$e*$e, 523=>$e*$e, 525=>$e*$e, 526=>$e*$e, 531=>$e*$e, 532=>$e*$e, 533=>$e*$e, 535=>$e*$e, 536=>$e*$e, 551=>$e*$e, 552=>$e*$e, 553=>$e*$e, 555=>$e*$e, 556=>$e*$e, 561=>$e*$e, 562=>$e*$e, 563=>$e*$e, 565=>$e*$e, 566=>$e*$e, 611=>$e*$e, 612=>$e*$e, 613=>$e*$e, 615=>$e*$e, 616=>$e*$e, 621=>$e*$e, 622=>$e*$e, 623=>$e*$e, 625=>$e*$e, 626=>$e*$e, 631=>$e*$e, 632=>$e*$e, 633=>$e*$e, 635=>$e*$e, 636=>$e*$e, 651=>$e*$e, 652=>$e*$e, 653=>$e*$e, 655=>$e*$e, 656=>$e*$e, 661=>$e*$e, 662=>$e*$e, 663=>$e*$e, 665=>$e*$e, 666=>$e*$e);
	} else {
		%oprior = (111=>1-$e, 112=>0.25*$e, 113=>0.25*$e, 115=>0.25*$e, 116=>0.25*$e, 121=>0.5*(1-$e), 122=>0.5*(1-$e), 123=>0.3333*$e, 125=>0.3333*$e, 126=>0.3333*$e, 131=>0.25*$e, 132=>1-$e, 133=>0.25*$e, 135=>0.25*$e, 136=>0.25*$e, 151=>0.3333*$e, 152=>0.5*(1-$e), 153=>0.5*(1-$e), 155=>0.3333*$e, 156=>0.3333*$e, 161=>0.5*$e, 162=>0.25*(1-$e), 163=>0.5*(1-$e), 165=>0.25*(1-$e), 166=>0.5*$e, 211=>0.5*(1-$e), 212=>0.5*(1-$e), 213=>0.3333*$e, 215=>0.3333*$e, 216=>0.3333*$e, 221=>0.25*(1-$e), 222=>0.5*(1-$e), 223=>0.25*(1-$e), 225=>0.5*$e, 226=>0.5*$e, 231=>0.3333*$e, 232=>0.5*(1-$e), 233=>0.5*(1-$e), 235=>0.3333*$e, 236=>0.3333*$e, 251=>0.5*$e, 252=>0.25*(1-$e), 253=>0.5*(1-$e), 255=>0.25*(1-$e), 256=>0.5*$e, 261=>$e, 262=>0.125*(1-$e), 263=>0.375*(1-$e), 265=>0.375*(1-$e), 266=>0.125*(1-$e), 311=>0.25*$e, 312=>1-$e, 313=>0.25*$e, 315=>0.25*$e, 316=>0.25*$e, 321=>0.3333*$e, 322=>0.5*(1-$e), 323=>0.5*(1-$e), 325=>0.3333*$e, 326=>0.3333*$e, 331=>0.25*$e, 332=>0.25*$e, 333=>1-$e, 335=>0.25*$e, 336=>0.25*$e, 351=>0.3333*$e, 352=>0.3333*$e, 353=>0.5*(1-$e), 355=>0.5*(1-$e), 356=>0.3333*$e, 361=>0.5*$e, 362=>0.5*$e, 363=>0.25*(1-$e), 365=>0.5*(1-$e), 366=>0.25*(1-$e), 511=>0.3333*$e, 512=>0.5*(1-$e), 513=>0.5*(1-$e), 515=>0.3333*$e, 516=>0.3333*$e, 521=>0.5*$e, 522=>0.25*(1-$e), 523=>0.5*(1-$e), 525=>0.25*(1-$e), 526=>0.5*$e, 531=>0.3333*$e, 532=>0.3333*$e, 533=>0.5*(1-$e), 535=>0.5*(1-$e), 536=>0.3333*$e, 551=>0.5*$e, 552=>0.5*$e, 553=>0.25*(1-$e), 555=>0.5*(1-$e), 556=>0.25*(1-$e), 561=>0.5*$e, 562=>0.5*$e, 563=>0.125*(1-$e), 565=>0.375*(1-$e), 566=>0.5*(1-$e), 611=>0.5*$e, 612=>0.25*(1-$e), 613=>0.5*(1-$e), 615=>0.25*(1-$e), 616=>0.5*$e, 621=>$e, 622=>0.125*(1-$e), 623=>0.375*(1-$e), 625=>0.375*(1-$e), 626=>0.125*(1-$e), 631=>0.5*$e, 632=>0.5*$e, 633=>0.25*(1-$e), 635=>0.5*(1-$e), 636=>0.25*(1-$e), 651=>0.5*$e, 652=>0.5*$e, 653=>0.125*(1-$e), 655=>0.375*(1-$e), 656=>0.5*(1-$e), 661=>0.5*$e, 662=>0.5*$e, 663=>0.0625*(1-$e), 665=>0.25*(1-$e), 666=>0.6875*(1-$e));
	}

#dichotomize the CNV validation procedure:
#if there is no boundary discordance, do a simple Bayesian computation
#otherwise, use viterbi algorithm to infer the most likely path
	for my $nextregion (@region) {
		my ($nextchr, @nextpos) = @$nextregion;				#nextchr is the chromosome, nextname is the SNPs that constitute the CNV boundary
			my ($delta, $psi);						#delta: maximum prob at each time, psi: maximum previous ind that generate the delta at this time
			my $nextsubregion = join (",", $nextchr, $nextpos[0], $nextpos[1]);
		my ($maxval, $maxvalind);
		$verbose and print "processing regions at chr$nextchr pos=@nextpos\n";

#the following paragraph applies when there is no boundary discordance (a simple Bayesian method to calculate most likely a posterior state combintations)
		if (@nextpos == 2) {
			validateTrioCNVCall ($ref_inputfile, \@fmprior, \%oprior, $nextsubregion, $sample_sex, $region_alllogprob, $region_numsnp, $region_cnvlength, $region_chr, $region_name_start, $region_name_end, $region_pos_start, $region_pos_end);
			next;
		}

#proceed to the following only if there are multiple overlapped subregions in the region. We will use a Viterbi algorithm to decode the best path
#delta is a matrix that record for each time point, for each fmostate, the maxlogprob of the path
#psi is a matrix that record for each time point, for each fmostate, the previous fmostate (that reach the current one with highest prob)

#step1: initialization
		for my $fstate (1, 2, 3, 5, 6) {
			for my $mstate (1, 2, 3, 5, 6) {
				for my $ostate (1, 2, 3, 5, 6) {
					my $fmoindex = $fstate*7*7+$mstate*7+$ostate;

					$psi->[0][$fmoindex] = 0;
					$delta->[0][$fmoindex] += log ($fmprior[$fstate]) + log ($fmprior[$mstate]) + log ($oprior{$fstate.$mstate.$ostate});
					$delta->[0][$fmoindex] += $region_alllogprob->{$nextsubregion}->[0][$fstate] + $region_alllogprob->{$nextsubregion}->[1][$mstate] + $region_alllogprob->{$nextsubregion}->[2][$ostate];
				}
			}
		}
		$verbose and print STDERR "NOTICE: initial scanning at $nextchr: @nextpos[0..1] found $maxvalind $maxval\n";

#step 2: recursion
		for my $i (1 .. @nextpos-2) {				#this is equivalent to "from time 1 to time t" (time 0 is the initialization stage)
			my $nextsubregion = join (",", $nextchr, $nextpos[$i], $nextpos[$i+1]);
			for my $fstate (1, 2, 3, 5, 6) {
				for my $mstate (1, 2, 3, 5, 6) {
					for my $ostate (1, 2, 3, 5, 6) {
						my $fmoindex = $fstate*7*7+$mstate*7+$ostate;
						my $old_fmoindex;
						undef $maxval;

						for my $old_fstate (1, 2, 3, 5, 6) {
							for my $old_mstate (1, 2, 3, 5, 6) {
								for my $old_ostate (1, 2, 3, 5, 6) {
									$old_fmoindex = $old_fstate*7*7+$old_mstate*7+$old_ostate;
									$region_snpdist->{$nextsubregion} or confess "\nERROR: For nextregion $nextsubregion, pos=@nextpos @$ref_inputfile: no snpdist defined!!!";
									my $val = $delta->[$i-1][$old_fmoindex] + log ($transition->[$old_fstate][$fstate] * (1-exp(-$region_snpdist->{$nextsubregion}/100_000)) / (1-exp(-5000/100_000))) + log ($transition->[$old_mstate][$mstate]* (1-exp(-$region_snpdist->{$nextsubregion}/100_000)) / (1-exp(-5000/100_000)));
									if ($fstate == $old_fstate and $mstate == $old_mstate and $ostate != $old_ostate and $oprior{$fstate.$mstate.$ostate} > $e) {
										my @tempoprior = sort {$a<=>$b} @oprior{$fstate.$mstate.'1', $fstate.$mstate.'2', $fstate.$mstate.'3', $fstate.$mstate.'5', $fstate.$mstate.'6'};
										$val += log ($tempoprior[0]);			#recombination happens without de novo event!
									} else {
										$val += log ($oprior{$fstate.$mstate.$ostate});
									}

#for the last block, arbitrarily set a "end" state that returns to state 3 (since all fstate and mstate return to 3, it can be directly added here as if there are two transitions)
									if ($i == @nextpos-2) {
										$val += log ($transition->[$fstate][3]) + log ($transition->[$mstate][3]);
									}

									if (not defined $maxval or $val > $maxval) {
										$maxval = $val;
										$maxvalind = $old_fmoindex;
									}
								}
							}
						}
						my $subpost = $region_alllogprob->{$nextsubregion}->[0][$fstate] + $region_alllogprob->{$nextsubregion}->[1][$mstate] + $region_alllogprob->{$nextsubregion}->[2][$ostate];

						$delta->[$i][$fmoindex] = $maxval + $subpost;
						$psi->[$i][$fmoindex] = $maxvalind;
					}
				}
			}
		}

		my ($pprob, $q, $qstate);			#pprob is the maximum prob, $q is the fmoindex that reach maximum prob, qstate is the [fstate,mstate,ostate]
#step 3: termination
			for my $fstate (1, 2, 3, 5, 6) {
				for my $mstate (1, 2, 3, 5, 6) {
					for my $ostate (1, 2, 3, 5, 6) {
						my $fmoindex = $fstate*7*7+$mstate*7+$ostate;
						if (not defined $pprob or $delta->[@nextpos-2][$fmoindex] > $pprob) {
							$pprob = $delta->[@nextpos-2][$fmoindex];
							$q->[@nextpos-2] = $fmoindex;
							$qstate->[@nextpos-2] = [convertTrioIndex2State ($fmoindex)];
						}
					}
				}
			}
		$verbose and print STDERR "q[", @nextpos-2, "]=", $q->[@$q-1], " (", convertTrioIndex2State ($q->[@$q-1]), ")\n";

#step 4: path backtracing
		for (my $i = @nextpos-3; $i >= 0; $i--) {
			$q->[$i] = $psi->[$i+1][$q->[$i+1]];
			$qstate->[$i] = [convertTrioIndex2State ($q->[$i])];
			$verbose and print STDERR "q[$i] = $q->[$i] (", convertTrioIndex2State ($q->[$i]), ")\n";
		}

#finally, compile all state sequences together to get a concensus sequence
		my @cnv;					#@cnv has 3 element, corresponding to father, mother, offspring CNV calls
			for my $i (0 .. 2) {
				my ($stretch_start_j, $found_signal);
				my $normal_state = 3;
				$chrx and $sample_sex->[$i] eq 'male' and $normal_state = 2;		#for male chrx
					$chry and $sample_sex->[$i] eq 'male' and $normal_state = 2;		#for male chry
					$chry and $sample_sex->[$i] eq 'female' and $normal_state = 1;		#for female chry
					for my $j (0 .. @$qstate-1) {
						my $nextstate = $qstate->[$j][$i];
						if ($nextstate ne $normal_state) {
							if ($found_signal and $found_signal ne $nextstate) {	#transition between different CNV states
								push @{$cnv[$i]}, [$stretch_start_j, $j, $found_signal];
								$stretch_start_j = $j;
								$found_signal = $nextstate;
							} elsif ($found_signal) {
								1;						#do nothing, still in the same state
							} else {
								$found_signal = $nextstate;
								$stretch_start_j = $j;
							}
						} else {
							if ($found_signal) {
								push @{$cnv[$i]}, [$stretch_start_j, $j, $found_signal];
								$found_signal = 0;
							}
						}
					}
#finish the last stretch
				if ($found_signal) {
					push @{$cnv[$i]}, [$stretch_start_j, scalar (@$qstate), $found_signal];
				}
			}

		printFamilyCNVCall (\@cnv, $nextchr, \@nextpos, $region_name_start, $region_name_end, $qstate, $siginfo, $ref_inputfile);

	}
}

sub printFamilyCNVCall {
	my ($cnv, $nextchr, $nextpos, $region_name_start, $region_name_end, $qstate, $siginfo, $ref_inputfile) = @_;

	for my $i (0 .. @$cnv-1) {
		for my $stretch (@{$cnv->[$i]}) {
			my ($starti, $endi, $nextstate) = @$stretch;
			my $nextsubregion1 = join (",", $nextchr, $nextpos->[$starti], $nextpos->[$starti+1]);
			my $nextsubregion2 = join (",", $nextchr, $nextpos->[$endi-1], $nextpos->[$endi]);
			my $namestart = $region_name_start->{$nextsubregion1};			#the actual namestart might be different from the start_snp
				my $nameend = $region_name_end->{$nextsubregion2};			#the actual nameend might be different from the end_snp
				my $identity = ("father", "mother", "offspring", "offspring")[$i];
			my $nextcn = $nextstate-1;
			$nextcn >= 3 and $nextcn--;						#calculate current copy number
				my $triostateseq;

			for my $tempi ($starti .. $endi-1) {			#endi is usually the last element of nextpos, but qstate has only endi-1 elements
				$triostateseq .= join ("", @{$qstate->[$tempi]})."-";
			}
			$triostateseq =~ s/-$//;
			my ($snp_starti, $snp_endi, $numsnp, $cnvlength);
			for my $j (0 .. @{$siginfo->{$nextchr}{name}}-1) {
				$siginfo->{$nextchr}{name}[$j] eq $namestart and $snp_starti = $j;
				$siginfo->{$nextchr}{name}[$j] eq $nameend and $snp_endi = $j;		#endi+1 must be used because we are asking for end SNP in a stretch
					defined $snp_starti and defined $snp_endi and last;
			}

			$numsnp = $snp_endi-$snp_starti+1;
			if (defined $minsnp) {						#minimum number of SNPs within a predicted CNV to be printed out
				$numsnp >= $minsnp or next;
			}
			if (length($numsnp) < 6) {
				$numsnp = substr ("$numsnp      ", 0, 6);
			}

			$cnvlength = $siginfo->{$nextchr}{pos}[$snp_endi]-$siginfo->{$nextchr}{pos}[$snp_starti]+1;
			if (defined $minlength) {					#minimum length of CNVs to be printed out
				$cnvlength >= $minlength or next;
			}
			$cnvlength = join ('', reverse split (//, $cnvlength)); $cnvlength =~ s/(...)/$1,/g; $cnvlength =~ s/,$//; $cnvlength = join ('', reverse split (//, $cnvlength));
			if (length($cnvlength) < 11) {
				$cnvlength = substr ("$cnvlength            ", 0, 11);
			}

			my $cnvregion = "chr" . $nextchr . ":" . $siginfo->{$nextchr}{pos}[$snp_starti] . "-" . $siginfo->{$nextchr}{pos}[$snp_endi];
			if (length ($cnvregion) < 29) {
				$cnvregion = substr ("$cnvregion                              ", 0, 29);
			}

			print "$cnvregion numsnp=$numsnp length=$cnvlength state$nextstate,cn=$nextcn $ref_inputfile->[$i] ";
			print "startsnp=$namestart endsnp=$nameend $identity boundary_reconciled=$triostateseq\n";			
		}
	}
}

sub retrieveCandidateRegion {
	my ($candlist) = @_;
	my @candregion;
	print STDERR "NOTICE: Reading candidate regions to be validated from candlist file $candlist ...";
	open (CAND, $candlist) or confess "Error: cannot read from candlist file $candlist: $!\n";

	while (<CAND>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		@record == 5 or confess "Error: the --candlist file $candlist is not in valid format (five tab-delimited records expected in each line): <$_>\n";
		$record[0] eq 'region' and next;		#header line in the candlist file

			my ($chr_region, $startsnp, $endsnp, $delfreq, $dupfreq) = @record;
		$delfreq eq 'null' and $dupfreq eq 'null' and next;		#skip this region
			$delfreq eq 'null' and $delfreq = $backfreq;
		$dupfreq eq 'null' and $dupfreq = $backfreq;
		$delfreq =~ m/^[\d\.]+$/ or confess "Error: invalid record found in candlist file $candlist (delfreq is specified as $delfreq): <$_>\n";
		$dupfreq =~ m/^[\d\.]+$/ or confess "Error: invalid record found in candlist file $candlist (dupfreq is specified as $dupfreq): <$_>\n";

		$delfreq <= 1 and $delfreq >=0 or confess "Error: the --candlist file $candlist contains invalid delfreq (must be between 0 and 1) in line: <$_>\n";
		$dupfreq <= 1 and $dupfreq >=0 or confess "Error: the --candlist file $candlist contains invalid dupfreq (must be between 0 and 1) in line: <$_>\n";
		$delfreq+$dupfreq<1 or confess "Error: the --candlist file $candlist contains invalid delfreq and dupfreq (the sum must be less than 1) in line: <$_>\n";
		$delfreq ||= 1e-30;		#convert zero to a small number to prevent underflow error
			$dupfreq ||= 1e-30;		#convert zero to a small number to prevent underflow error
			push @candregion, [$startsnp, $endsnp, convertDelDupFreqToPrior ($delfreq, $dupfreq, $backfreq)];
		$verbose and print STDERR "NOTICE: For validating CNV calls of $startsnp-$endsnp, the prior distribution is set at ", join (",", convertDelDupFreqToPrior ($delfreq, $dupfreq, $backfreq)), "\n";
	}
	close (CAND);
	print STDERR " Done with ${\(scalar @candregion)} regions\n";
	return (\@candregion);
}

sub convertDelDupFreqToPrior {
	my ($delfreq, $dupfreq, $backfreq) = @_;
	my @prior;

	if (not defined $delfreq or not defined $dupfreq) {
		defined $backfreq or confess "ERROR: convertDelDupFreqToPrior: no backfreq defined\n";				#by default, the background freq means that 0.01% of whole genome markers are in CNV
	}
	defined $delfreq or $delfreq = $backfreq;
	defined $dupfreq or $dupfreq = $backfreq;
	my $normfreq = 1-$delfreq-$dupfreq;					#allele frequqency of the normal-copy allele
		$delfreq > 0 and $delfreq < 1 or confess "ERROR: the deletion frequency ($delfreq) must be specified between 0 and 1\n";
	$dupfreq > 0 and $dupfreq < 1 or confess "ERROR: the duplication frequency ($dupfreq) must be specified between 0 and 1\n";
	$delfreq+$dupfreq<1 or confess "ERROR: the deletion frequency ($delfreq) and duplication frequency ($dupfreq) must sum up less than 1\n";

	$prior[0] = $delfreq*$delfreq;
	$prior[1] = 2*$delfreq*$normfreq;
	$prior[3] = 2*$dupfreq*$normfreq;
	$prior[4] = $dupfreq*$dupfreq;

	$prior[2] = 1-$prior[0]-$prior[1]-$prior[3]-$prior[4];
	return @prior;
}

sub QCSignal {
	my ($inputfile, $siginfo, $sigdesc, $hmm, $known_sex) = @_;
	my ($sample_sex, $use_sex, $qc_flag);			#predicted gender, the gender that should be used, a QC-flag that specifies problems with QC

		if ($chrx) {
			if ($sigdesc->{baf_xhet} eq 'NA') {
				print STDERR "\nERROR: skipping inputfile $inputfile since it does not contain signal information for X chromosome markers.\n";
				confess "FATAL ERROR: chrX markers not found in inputfile $inputfile for chrX operation (but the program reaches the QC step)\n";
			}
			print STDERR "NOTICE: quality summary for $inputfile: LRR_Xmean=", sprintf ("%.4f", $sigdesc->{lrr_xmean}), " LRR_Xmedian=", sprintf ("%.4f", $sigdesc->{lrr_xmedian}), " LRR_XSD=", sprintf ("%.4f", $sigdesc->{lrr_xsd}), " BAF_Xhet=", sprintf ("%.4f", $sigdesc->{baf_xhet}), "\n";

			$sample_sex = $sigdesc->{baf_xhet} > $bafxhet_threshold ? 'female' : 'male';
			if (not $known_sex) {
				print STDERR "NOTICE: Sample sex for $inputfile is predicted as '$sample_sex' based on BAF heterozygosity rate (this is different from genotype heterozygosity rate!) for chrX ($sigdesc->{baf_xhet})\n";
			}

			$use_sex = $known_sex||$sample_sex;

			if ($medianadjust) {					#for XXX, XXY and so on, the user may not want to use --medianadjust
				if ($use_sex eq 'male') {			#for males, the mean LRR should be equal to expected LRR for state 2 (one-copy)
					$verbose and print STDERR "NOTICE: Adjusting LRR values in chrX from $inputfile by ", sprintf ("%.4f", $sigdesc->{lrr_xmedian}-$hmm->{B1_mean}[1]), "\n";
					adjustLRR ($siginfo, $sigdesc->{lrr_xmedian}-$hmm->{B1_mean}[1]);
				} else {
					$verbose and print STDERR "NOTICE: Adjusting LRR values in chrX from $inputfile by ", sprintf ("%.4f", $sigdesc->{lrr_xmedian}), "\n";
					adjustLRR ($siginfo, $sigdesc->{lrr_xmedian});
				}
			}

#quality control: examine the variation of Log R Ratio values (>0.2 is generally treated as bad genotyping quality)
			if ($sigdesc->{lrr_xsd} > 0.2) {
				print STDERR "WARNING: Sample from $inputfile does not pass quality control criteria due to its large SD for LRR ($sigdesc->{lrr_xsd})!\n";
				$qc_flag++;
			}

#quality control: examine the waviness factor values (>0.04 or <-0.04 is generally treated as bad genotyping quality)
			if ($sigdesc->{wf} ne 'NA' and $sigdesc->{wf} > 0.04 || $sigdesc->{wf} < -0.04) {
				print STDERR "WARNING: Sample from $inputfile does not pass quality control criteria due to its waviness factor values ($sigdesc->{wf})!\n";
				$qc_flag++;
			}
		} elsif ($chry) {
			if ($sigdesc->{baf_yhet} eq 'NA') {
				print STDERR "\nERROR: skipping inputfile $inputfile since it does not contain signal information for Y chromosome markers.\n";
				confess "FATAL ERROR: chrY markers not found in inputfile $inputfile for chrY operation (but the program reaches the QC step)\n";
			}
			print STDERR "NOTICE: quality summary for $inputfile: LRR_Ymean=", sprintf ("%.4f", $sigdesc->{lrr_ymean}), " LRR_Ymedian=", sprintf ("%.4f", $sigdesc->{lrr_ymedian}), " LRR_YSD=", sprintf ("%.4f", $sigdesc->{lrr_ysd}), " BAF_Yhet=", sprintf ("%.4f", $sigdesc->{baf_yhet}), "\n";

## sex is also infered from chrX heterozygosity rate
			$sample_sex = $sigdesc->{baf_xhet} > $bafxhet_threshold ? 'female' : 'male';
			if (not $known_sex) {
				print STDERR "NOTICE: Sample sex for $inputfile is predicted as '$sample_sex' based on BAF heterozygosity rate (this is different from genotype heterozygosity rate!) for chrX ($sigdesc->{baf_xhet})\n";
			}
			$use_sex = $known_sex||$sample_sex;

			if ($medianadjust) {		
#for XXX, XXY and so on, the user may not want to use --medianadjust			
#for males, the mean LRR should be equal to expected LRR for state 2 (one-copy)
#for females, the mean LRR should be equal to expected LRR for state 1 (zero-copy)
				if ($use_sex eq 'male') {			
					$verbose and print STDERR "NOTICE: Adjusting LRR values in chrY from $inputfile by ", sprintf ("%.4f", $sigdesc->{lrr_ymedian}-$hmm->{B1_mean}[1]), "\n";
					adjustLRR ($siginfo, $sigdesc->{lrr_ymedian}-$hmm->{B1_mean}[1]);
				} else {
					$verbose and print STDERR "NOTICE: Adjusting LRR values in chrY from $inputfile by ", sprintf ("%.4f", $sigdesc->{lrr_ymedian}-$hmm->{B1_mean}[0]), "\n";
					adjustLRR ($siginfo, $sigdesc->{lrr_ymedian}-$hmm->{B1_mean}[0]);
				}
			}

#quality control: examine the variation of Log R Ratio values (>0.2 is generally treated as bad genotyping quality)
			if ($sigdesc->{lrr_ysd} > 0.2 and $use_sex eq 'male') {
				print STDERR "WARNING: Sample from $inputfile does not pass quality control criteria due to its large SD for LRR ($sigdesc->{lrr_ysd})!\n";
				$qc_flag++;
			}

#quality control: examine the waviness factor values (>0.04 or <-0.04 is generally treated as bad genotyping quality)
			if ($sigdesc->{wf} ne 'NA' and $sigdesc->{wf} > 0.04 || $sigdesc->{wf} < -0.04) {
				print STDERR "WARNING: Sample from $inputfile does not pass quality control criteria due to its waviness factor values ($sigdesc->{wf})!\n";
				$qc_flag++;
			}
		} else {
#perform signal pre-processing to adjust the median LRR and BAF values
			if ($medianadjust) {				#THIS IS A MANDATORY ARGUMENT NOW, SINCE IT ALWAYS IMPROVE PERFORMANCE!!! (use "--nomedianadjust" to disable this feature)
				print STDERR "NOTICE: Median-adjusting LRR values for all autosome markers from $inputfile by ", sprintf ("%.4f", $sigdesc->{lrr_median}), "\n";
				adjustLRR ($siginfo, $sigdesc->{lrr_median});
				$sigdesc->{lrr_mean} -= $sigdesc->{lrr_median};
				$sigdesc->{lrr_median} = 0;
			}
			if ($bafadjust) {
				if ($sigdesc->{baf_median} ne "NA") {	#for oligonucleotide arrays, there is no BAF information
					print STDERR "NOTICE: Median-adjusting BAF values for all autosome markers from $inputfile by ", sprintf ("%.4f", $sigdesc->{baf_median}-0.5), "\n";
					adjustBAF ($siginfo, $sigdesc->{baf_median}-0.5);
					$sigdesc->{baf_mean} -= $sigdesc->{baf_median}-0.5;
					$sigdesc->{baf_median} = 0.5;
				} 
			}

			if ($sigdesc->{baf_median} eq "NA") {
				$sigdesc->{baf_mean} = 0.5;
				$sigdesc->{baf_median} = 0.5;
				$sigdesc->{baf_sd} = 0;
				$sigdesc->{baf_drift} = 0;
			}

			print STDERR "NOTICE: quality summary for $inputfile: LRR_mean=", sprintf ("%.4f", $sigdesc->{lrr_mean}), " LRR_median=", sprintf ("%.4f", $sigdesc->{lrr_median}), " LRR_SD=", sprintf ("%.4f", $sigdesc->{lrr_sd}), 
				  " BAF_mean=", sprintf ("%.4f", $sigdesc->{baf_mean}), " BAF_median=", sprintf ("%.4f", $sigdesc->{baf_median}), " BAF_SD=", sprintf ("%.4f", $sigdesc->{baf_sd}), " BAF_DRIFT=", sprintf ("%.6f", $sigdesc->{baf_drift}), " WF=", $sigdesc->{wf} eq "NA"?"NA":sprintf("%.4f", $sigdesc->{wf}), " GCWF=", $sigdesc->{gcwf} eq "NA"?"NA":sprintf("%.4f", $sigdesc->{gcwf}), "\n";

#quality control: examine the variation of Log R Ratio values (>0.2 is generally treated as bad genotyping quality)
			if ($sigdesc->{lrr_sd} > 0.2) {
				print STDERR "WARNING: Sample from $inputfile does not pass default quality control criteria due to its large SD for LRR ($sigdesc->{lrr_sd})!\n";
				$qc_flag++;
			}

#quality control: examine the median of BAF values (>0.6 or <0.4 is treated as bad clustering quality)
			if ($sigdesc->{baf_median} < 0.4 or $sigdesc->{baf_median} > 0.6) {
				print STDERR "WARNING: Sample from $inputfile does not pass default quality control criteria due to its shifted BAF values (median=$sigdesc->{baf_median})!\n";
				$qc_flag++;
			}

#quality control: examine the drifting of BAF values (>0.002 is generally treated as bad genotyping quality)
			if ($sigdesc->{baf_drift} > 0.002) {
				print STDERR "WARNING: Sample from $inputfile does not pass default quality control criteria due to its drifting BAF values (drift=$sigdesc->{baf_drift})!\n";
				$qc_flag++;
			}
		}

#quality control: examine the waviness factor values (>0.04 or <-0.04 is generally treated as bad genotyping quality)
	if ($sigdesc->{wf} ne 'NA' and $sigdesc->{wf} > 0.04 || $sigdesc->{wf} < -0.04) {
		print STDERR "WARNING: Sample from $inputfile does not pass default quality control criteria due to its waviness factor values (wf=$sigdesc->{wf})!\n";
		$qc_flag++;
	}
	$sigdesc->{nocall_rate} > 0.1 and print STDERR "WARNING: Sample from $inputfile does not pass default quality control criteria due to its NoCall rate ($sigdesc->{nocall_rate})!\n" and $qc_flag++;
	$qc_flag and print STDERR "WARNING: Small-sized CNV calls may not be reliable and should be interpreted with caution!\n";

#check whether the signal file contains non-polymorphic markers and whether the HMM file contains parameters to handle them
	if ($sigdesc->{cn_count}) {
		if (not $hmm->{B3_mean}) {
			print STDERR "WARNING: The signal file $inputfile contains $sigdesc->{cn_count} non-polymorphic markers but the HMM file $hmmfile is for SNP markers only!\n";
			print STDERR "WARNING: To get the most accurate CNV calls, please use a HMM file containing parameters specifically for non-polymorphic markers.\n";
		}
	}
	if ($known_sex) {
		if ($known_sex ne $sample_sex) {
			print STDERR "WARNING: sample sex for $inputfile is predicted as $sample_sex based on chrX, but is specified as $known_sex by user. It is treated as $known_sex in analysis\n";
		}
		return $known_sex;
	} else {
		return $sample_sex;
	}
}

sub newvalidateCNVCall {
	my ($ref_inputfile, $hmmfile, $pfbfile, $sexfile, $gcmodelfile, $directory, $candregion) = @_;
	my ($pfbinfo) = newreadPFB ($pfbfile);
	my $hmm = readHMMFile ($hmmfile);
	my ($file_sex, $gcmodel) = ({});

	defined $sexfile and $file_sex = readSexFile ($sexfile);
	defined $gcmodelfile and $gcmodel = newreadGCModel ($gcmodelfile, $pfbinfo);	#read the gc model for SNPs defined in snp_chr

		my $valioutfh;
	defined $valilog and open ($valioutfh, ">$valilog") || confess "Error: cannot write to valilog file $valilog: $!\n";

	for my $inputfile (@$ref_inputfile) {
		my ($siginfo, $sigdesc, $sample_sex);

		($siginfo, $sigdesc) = readLRRBAF ($inputfile, $region, $pfbinfo, $gcmodel, $directory);

#if the file processing fails and no signal data is read from file, move on to the next file
		%$siginfo or print STDERR "WARNING: Skipping $inputfile since no signal values can be retrieved from the file\n" and next;

		$sample_sex = QCSignal ($inputfile, $siginfo, $sigdesc, $hmm, $file_sex->{$inputfile});

#read HMM model file
		my $hmm_model = khmm::ReadCHMM ($hmmfile);

#Now a mandatory step to handle low-quality samples and reduce false-positve calls, through matching the SD measure for sample and model
		if ($sdadjust) {
			if ($chrx) {
				khmm::adjustBSD ($hmm_model, $sigdesc->{lrr_xsd});
			} elsif ($chry) {
				khmm::adjustBSD ($hmm_model, $sigdesc->{lry_xsd});
			} else {
				khmm::adjustBSD ($hmm_model, $sigdesc->{lrr_sd});
			}
		}

		for my $i (0 .. @$candregion-1) {
			my ($startsnp, $endsnp, @prior) = @{$candregion->[$i]};
			my $snp_chr;
			my @temp = split (/\t/, $pfbinfo->{$startsnp});
			$snp_chr->{$startsnp} = $temp[0];
			@temp = split (/\t/, $pfbinfo->{$endsnp});
			$snp_chr->{$endsnp} = $temp[0];

			my $curchr = $snp_chr->{$startsnp} or confess "Error: the --startsnp $startsnp is not annotated in PFB file $pfbfile\n";
			$snp_chr->{$endsnp}  or confess "Error: the --endsnp $endsnp is not annotated in PFB file $pfbfile\n";
			$curchr eq $snp_chr->{$endsnp} or confess "Error: the --startsnp $startsnp is located in chr$curchr but --endsnp $endsnp is located in chr$snp_chr->{$endsnp}\n";

			my $name = $siginfo->{$curchr}{name};
			my $pos = $siginfo->{$curchr}{pos};
			my $lrr = $siginfo->{$curchr}{lrr};
			my $baf = $siginfo->{$curchr}{baf};

			my $cnvcall = validateRegion ($name, $pos, $lrr, $baf, $startsnp, $endsnp, $pfbinfo, $curchr, \@prior, $hmm_model, $sample_sex, $valioutfh);

#print all CNV calls, sorted by copy numbers first, then by chromosomes
			if ($tabout) {
				printTabbedCNV ($cnvcall, $inputfile);				#use tab-delimited output
			} else {
				printFormattedCNV ($cnvcall, $inputfile);
			}
		}
	}
}

#assign the score for each region (similar to assignConfidence subroutine)
sub validateRegion {
	my ($name, $pos, $lrr, $baf, $startsnp, $endsnp, $pfbinfo, $curchr, $prior_prob, $hmm_model, $sample_sex, $valioutfh) = @_;
	my (@region_lrr, @region_baf, @region_pfb, @region_snpdist);
	push @region_lrr, 0;
	push @region_baf, 0;
	push @region_pfb, 0;
	push @region_snpdist, 0;
	my (@logprob);
	my ($found, $last) = qw/0 0/;
	my $cnvcall = {};

	my ($snp_pos);

	for my $i (0 .. @$name-1) {
		if ($name->[$i] eq $startsnp) {
			$found = 1;
		}
		if ($name->[$i] eq $endsnp) {
			$last = 1;
		}
		if ($found) {
			push @region_lrr, $lrr->[$i];
			push @region_baf, $baf->[$i];
			my @temp = split (/\t/, $pfbinfo->{$name->[$i]});
			push @region_pfb, $temp[2];
			$snp_pos->{$name->[$i]} = $temp[1];
			push @region_snpdist, ($pos->[$i+1]||($pos->[$i]+1))-$pos->[$i];
		}
		$last and last;
	}

	@region_lrr or confess "UNKNOWN ERROR: unable to find CNV region startsnp=$startsnp endsnp=$endsnp\n";

	my @rawll;
	for my $nextstate (1, 2, 3, 5, 6) {		#do not consider LOH
		my $logprob = 0;
		my $stateindex = $nextstate-1;
		$nextstate > 4 and $stateindex--;
		khmm::GetStateProb_CHMM ($hmm_model, @region_lrr-1, \@region_lrr, \@region_baf, \@region_pfb, \@region_snpdist, \$logprob, $nextstate);
		push @logprob, [$logprob + log ($prior_prob->[$stateindex]), $nextstate];
#push @rawll, $logprob;
		push @rawll, $logprob + log ($prior_prob->[$stateindex]);
	}

	$valioutfh and print $valioutfh "$startsnp\t$endsnp\t", join ("\t", @rawll), "\n";

	@logprob = sort {$b->[0]<=>$a->[0]} @logprob;
	my $beststate = $logprob[0]->[1];	#best state should be the first state! 20160621
		$verbose and print STDERR "NOTICE: best state is $beststate\n";

	if ($chrx) {
		$sample_sex or confess "Error: the sample_sex is not defined for chrX CNV calling\n";
		if ($sample_sex eq 'male') {
			$beststate == 2 and return $cnvcall;
		}
	} elsif ($chry) {
		$sample_sex or confess "Error: the sample_sex is not defined for chrY CNV calling\n";
		if ($sample_sex eq 'male') {
			$beststate == 2 and return $cnvcall;
		} else {
			$beststate == 1 and return $cnvcall;
		}
	} else {
		$beststate == 3 and return $cnvcall;
	}

	if ($confidence) {
		if ($minconf) {
			$logprob[0]->[0]-$logprob[1]->[0] >= $minconf or return $cnvcall;
		}
		push @{$cnvcall->{$beststate}}, [$curchr, $snp_pos->{$startsnp}, $snp_pos->{$endsnp}, scalar (@region_lrr)-1, $startsnp, $endsnp, undef, undef, sprintf ("%.3f", $logprob[0]->[0]-$logprob[1]->[0])];
	} else {
		push @{$cnvcall->{$beststate}}, [$curchr, $snp_pos->{$startsnp}, $snp_pos->{$endsnp}, scalar (@region_lrr)-1, $startsnp, $endsnp];
	}
	return $cnvcall;
}

sub validateTrioCNVCall {
	my ($ref_inputfile, $fmprior, $oprior, $nextsubregion, $sample_sex, $region_alllogprob, $region_numsnp, $region_cnvlength, $region_chr, $region_name_start, $region_name_end, $region_pos_start, $region_pos_end) = @_;

	my ($maxfstate, $maxmstate, $maxostate, $maxfmostate, $maxpost);
	for my $fstate (1, 2, 3, 5, 6) {
		for my $mstate (1, 2, 3, 5, 6) {
			for my $ostate (1, 2, 3, 5, 6) {
				my $cnvgeno = $fstate . $mstate . $ostate;
				my $bayesposterior = log ($fmprior->[$fstate]) + log ($fmprior->[$mstate]) + log ($oprior->{$cnvgeno}) + $region_alllogprob->{$nextsubregion}->[0][$fstate] + $region_alllogprob->{$nextsubregion}->[1][$mstate] + $region_alllogprob->{$nextsubregion}->[2][$ostate];
				if (not defined $maxpost or $bayesposterior > $maxpost) {
					($maxpost, $maxfstate, $maxmstate, $maxostate) = ($bayesposterior, $fstate, $mstate, $ostate);
				}
				$verbose and print STDERR "NOTICE: cnvgeno=$cnvgeno F=$fstate M=$mstate C=$ostate P=$bayesposterior", ",", log ($fmprior->[$fstate]), ",", log ($fmprior->[$mstate]), ",", log ($oprior->{$cnvgeno}), "\n";
			}
		}
	}
	$maxfmostate = $maxfstate . $maxmstate . $maxostate;

#$verbose and print STDERR "NOTICE: States with highest logprob in region ", $region{"$nextregion->[0]:$nextregion->[1]-$nextregion->[2]"}, "is $maxfstate,$maxmstate,$maxostate (rawlogprob=$maxpost)\n";
	for my $i (0 .. 2) {
		my $nextstate = ($maxfstate, $maxmstate, $maxostate)[$i];
		if ($chrx) {
			if ($sample_sex->[$i] eq 'male') {
				$nextstate == 2 and next;
			} else {
				$nextstate == 3 and next;
			}
		} elsif ($chry) {
			if ($sample_sex->[$i] eq 'male') {
				$nextstate == 2 and next;
			} else {
				$nextstate == 1 and next;
			}
		} else {
			$nextstate == 3 and next;
		}
		my $nextcn = $nextstate-1;
		$nextcn >= 3 and $nextcn--;

		my $identity = ("father", "mother", "offspring")[$i];
		my $numsnp = substr ($region_numsnp->{$nextsubregion} . "       ", 0, 6);
		my $cnvlength = $region_cnvlength->{$nextsubregion};
		$cnvlength = join ('', reverse split (//, $cnvlength)); $cnvlength =~ s/(...)/$1,/g; $cnvlength =~ s/,$//; $cnvlength = join ('', reverse split (//, $cnvlength));
		$cnvlength = substr ("$cnvlength            ", 0, 12);
		print substr ("chr" . $region_chr->{$nextsubregion} . ":" . $region_pos_start->{$nextsubregion} . "-" . $region_pos_end->{$nextsubregion} . "                              ", 0, 30), "numsnp=$numsnp length=$cnvlength";
		print "state$nextstate,cn=$nextcn $ref_inputfile->[$i] ";
		print "startsnp=", $region_name_start->{$nextsubregion}, " endsnp=", $region_name_end->{$nextsubregion}, " $identity triostate=$maxfmostate\n";
	}
}

sub validateQuartetCNVCall {
	my ($ref_inputfile, $fmprior, $oprior, $nextsubregion, $sample_sex, $region_alllogprob, $region_numsnp, $region_cnvlength, $region_chr, $region_name_start, $region_name_end, $region_pos_start, $region_pos_end) = @_;

	my ($maxfstate, $maxmstate, $max_o1state, $max_o2state, $maxfmostate, $maxpost);
	for my $fstate (1, 2, 3, 5, 6) {
		for my $mstate (1, 2, 3, 5, 6) {
			my @ostate_array = generateOstateArray (@$ref_inputfile-2, [1,2,3,5,6]);

			for my $o12states (@ostate_array) {
				my ($o1state, $o2state) = @$o12states;
				my $cnvgeno1 = $fstate . $mstate . $o1state;
				my $cnvgeno2 = $fstate . $mstate . $o2state;
#$verbose and print "f=$fmprior[$fstate] m=$fmprior[$mstate] geno1=$cnvgeno1 o1=$oprior{$cnvgeno1}} $region_alllogprob{@$nextregion}->[0][$fstate] and $region_alllogprob{@$nextregion}->[1][$mstate] and $region_alllogprob{@$nextregion}->[2][$o1state]\n";
				my $bayesposterior = log ($fmprior->[$fstate]) + log ($fmprior->[$mstate]) + $region_alllogprob->{$nextsubregion}->[0][$fstate] + $region_alllogprob->{$nextsubregion}->[1][$mstate];
				$bayesposterior += log ($oprior->{$cnvgeno1})  + $region_alllogprob->{$nextsubregion}->[2][$o1state];
				$bayesposterior += log ($oprior->{$cnvgeno2})  + $region_alllogprob->{$nextsubregion}->[3][$o2state];

				$verbose and print STDERR "NOTICE: cnvgeno=$cnvgeno1 F=$fstate M=$mstate C=$o1state P=$bayesposterior", ",", log ($fmprior->[$fstate]), ",", log ($fmprior->[$mstate]), ",", log ($oprior->{$cnvgeno1}), "\n";
				if (not defined $maxpost or $bayesposterior > $maxpost) {
					($maxpost, $maxfstate, $maxmstate, $max_o1state, $max_o2state) = ($bayesposterior, $fstate, $mstate, $o1state, $o2state);
				}
				$verbose and print STDERR "NOTICE: max=$maxpost, $maxfstate, $maxmstate, $max_o1state, $max_o2state, current=$bayesposterior, $fstate, $mstate, $o1state, $o2state\n";

			}					
		}
	}
	$maxfmostate = join ('', $maxfstate, $maxmstate, $max_o1state, $max_o2state);

#$verbose and print STDERR "NOTICE: States with highest logprob in region ", $region{"$nextregion->[0]:$nextregion->[1]-$nextregion->[2]"}, "is $maxfstate,$maxmstate,$max_o1state,max_o2state (rawlogprob=$maxpost)\n";
	for my $i (0 .. @$ref_inputfile-1) {
		my $nextstate = ($maxfstate, $maxmstate, $max_o1state, $max_o2state)[$i];
		if ($chrx) {
			if ($sample_sex->[$i] eq 'male') {
				$nextstate == 2 and next;
			} else {
				$nextstate == 3 and next;
			}
		} elsif ($chry) {
			if ($sample_sex->[$i] eq 'male') {
				$nextstate == 2 and next;
			} else {
				$nextstate == 1 and next;
			}
		} else {
			$nextstate == 3 and next;
		}

		my $identity;
		if ($i >= 2) {
			$identity = 'offspring';
		} else {
			$identity = ("father", "mother")[$i];
		}

		my $nextcn = $nextstate-1;
		$nextcn >= 3 and $nextcn--;

		my $numsnp = substr ($region_numsnp->{$nextsubregion} . "       ", 0, 6);
		my $cnvlength = $region_cnvlength->{$nextsubregion};
		$cnvlength = join ('', reverse split (//, $cnvlength)); $cnvlength =~ s/(...)/$1,/g; $cnvlength =~ s/,$//; $cnvlength = join ('', reverse split (//, $cnvlength));
		$cnvlength = substr ("$cnvlength            ", 0, 12);
		print substr ("chr" . $region_chr->{$nextsubregion} . ":" . $region_pos_start->{$nextsubregion} . "-" . $region_pos_end->{$nextsubregion} . "                              ", 0, 30), "numsnp=$numsnp length=$cnvlength";
		print "state$nextstate,cn=$nextcn $ref_inputfile->[$i] ";
		print "startsnp=", $region_name_start->{$nextsubregion}, " endsnp=", $region_name_end->{$nextsubregion}, " $identity quartetstate=$maxfmostate\n";
	}
}

#validate CNV calls using father-mother-offspring1-offspring2 quartet information and reconcile boundary discordance using family information.
sub newtestQuartetCNVFile {
	my ($ref_inputfile, $hmmfile, $pfbfile, $cnvfile, $sexfile, $gcmodelfile, $directory) = @_;
#my ($name, $chr, $pos, $logr, $baf, $pfb, $snpdist);
	my ($pfbinfo) = newreadPFB ($pfbfile);
	my $file_sex;
	$sexfile and $file_sex = readSexFile ($sexfile);
	my $hmm = readHMMFile ($hmmfile);
	my $gcmodel;

	defined $gcmodelfile and $gcmodel = newreadGCModel ($gcmodelfile, $pfbinfo);	#read the gc model for SNPs defined in snp_chr

		my ($array_subregion, $array_region, $hash_region, $hash_startname, $hash_endname) = identifyCNVRegion ($cnvfile, $ref_inputfile, $pfbinfo);
	$array_subregion or return;		#no CNV region can be identified for posterior calling algorithm

		my @subregion = @$array_subregion;
	my @region = @$array_region;
	my %region = %$hash_region;
	my %startname = %$hash_startname;
	my %endname = %$hash_endname;


	my ($siginfo, $sample_sex, $region_alllogprob, $region_numsnp, $region_cnvlength, $region_chr, $region_name_start, $region_name_end, $region_pos_start, $region_pos_end, $region_snpdist) 
		= retrieveRegionSignal ($ref_inputfile, $file_sex, $hmmfile, $array_subregion, $pfbinfo, $gcmodel, $directory, $hash_startname, $hash_endname);

	$chrx and confess "ERROR: THIS VERSION OF PENNCNV DOES NOT SUPPORT CHRX CALLING FOR A QUARTET YET\n";

#now apply pedigree information to validate CNV calls and fine-map CNV boundaries
#first define the CNV inheritance matrix, and read the state transition matrix
	my $transition = readTransitionMatrix ($hmmfile);
	my ($e, %oprior1, %oprior2) = ($denovo_rate);
	if ($chrx and $sample_sex->[2] eq 'male') {		#first child in the family
		%oprior1 = (111=>1-$e, 112=>0.25*$e, 113=>0.25*$e, 115=>0.25*$e, 116=>0.25*$e, 121=>0.5*(1-$e), 122=>0.5*(1-$e), 123=>0.3333*$e, 125=>0.3333*$e, 126=>0.3333*$e, 131=>0.25*$e, 132=>1-$e, 133=>0.25*$e, 135=>0.25*$e, 136=>0.25*$e, 151=>0.3333*$e, 152=>0.5*(1-$e), 153=>0.5*(1-$e), 155=>0.3333*$e, 156=>0.3333*$e, 161=>0.5*$e, 162=>0.25*(1-$e), 163=>0.5*(1-$e), 165=>0.25*(1-$e), 166=>0.5*$e, 211=>0.25*$e, 212=>1-$e, 213=>0.25*$e, 215=>0.25*$e, 216=>0.25*$e, 221=>0.3333*$e, 222=>0.5*(1-$e), 223=>0.5*(1-$e), 225=>0.3333*$e, 226=>0.3333*$e, 231=>0.25*$e, 232=>0.25*$e, 233=>1-$e, 235=>0.25*$e, 236=>0.25*$e, 251=>0.3333*$e, 252=>0.3333*$e, 253=>0.5*(1-$e), 255=>0.5*(1-$e), 256=>0.3333*$e, 261=>0.5*$e, 262=>0.5*$e, 263=>0.25*(1-$e), 265=>0.5*(1-$e), 266=>0.25*(1-$e), 311=>0.25*$e, 312=>0.25*$e, 313=>1-$e, 315=>0.25*$e, 316=>0.25*$e, 321=>0.3333*$e, 322=>0.3333*$e, 323=>0.5*(1-$e), 325=>0.5*(1-$e), 326=>0.3333*$e, 331=>0.25*$e, 332=>0.25*$e, 333=>0.25*$e, 335=>1-$e, 336=>0.25*$e, 351=>0.3333*$e, 352=>0.3333*$e, 353=>0.3333*$e, 355=>0.5*(1-$e), 356=>0.5*(1-$e), 361=>0.3333*$e, 362=>0.3333*$e, 363=>0.3333*$e, 365=>0.25*(1-$e), 366=>0.75*(1-$e), 511=>0.25*$e, 512=>0.25*$e, 513=>0.25*$e, 515=>1-$e, 516=>0.25*$e, 521=>0.3333*$e, 522=>0.3333*$e, 523=>0.3333*$e, 525=>0.5*(1-$e), 526=>0.5*(1-$e), 531=>0.25*$e, 532=>0.25*$e, 533=>0.25*$e, 535=>0.25*$e, 536=>1-$e, 551=>0.25*$e, 552=>0.25*$e, 553=>0.25*$e, 555=>0.25*$e, 556=>1-$e, 561=>0.25*$e, 562=>0.25*$e, 563=>0.25*$e, 565=>0.25*$e, 566=>1-$e, 611=>0.25*$e, 612=>0.25*$e, 613=>0.25*$e, 615=>0.25*$e, 616=>1-$e, 621=>0.25*$e, 622=>0.25*$e, 623=>0.25*$e, 625=>0.25*$e, 626=>1-$e, 631=>0.25*$e, 632=>0.25*$e, 633=>0.25*$e, 635=>0.25*$e, 636=>1-$e, 651=>0.25*$e, 652=>0.25*$e, 653=>0.25*$e, 655=>0.25*$e, 656=>1-$e, 661=>0.25*$e, 662=>0.25*$e, 663=>0.25*$e, 665=>0.25*$e, 666=>1-$e);
	} elsif ($chrx and $sample_sex->[2] eq 'female') {
		%oprior1 = (111=>1-$e, 112=>0.25*$e, 113=>0.25*$e, 115=>0.25*$e, 116=>0.25*$e, 121=>0.5*(1-$e), 122=>0.5*(1-$e), 123=>0.3333*$e, 125=>0.3333*$e, 126=>0.3333*$e, 131=>0.25*$e, 132=>1-$e, 133=>0.25*$e, 135=>0.25*$e, 136=>0.25*$e, 151=>0.3333*$e, 152=>0.5*(1-$e), 153=>0.5*(1-$e), 155=>0.3333*$e, 156=>0.3333*$e, 161=>0.5*$e, 162=>0.25*(1-$e), 163=>0.5*(1-$e), 165=>0.25*(1-$e), 166=>0.5*$e, 211=>1-$e, 212=>0.25*$e, 213=>0.25*$e, 215=>0.25*$e, 216=>0.25*$e, 221=>0.5*(1-$e), 222=>0.5*(1-$e), 223=>0.3333*$e, 225=>0.3333*$e, 226=>0.3333*$e, 231=>0.25*$e, 232=>1-$e, 233=>0.25*$e, 235=>0.25*$e, 236=>0.25*$e, 251=>0.3333*$e, 252=>0.5*(1-$e), 253=>0.5*(1-$e), 255=>0.3333*$e, 256=>0.3333*$e, 261=>0.5*$e, 262=>0.25*(1-$e), 263=>0.5*(1-$e), 265=>0.25*(1-$e), 266=>0.5*$e, 311=>1-$e, 312=>0.25*$e, 313=>0.25*$e, 315=>0.25*$e, 316=>0.25*$e, 321=>0.5*(1-$e), 322=>0.5*(1-$e), 323=>0.3333*$e, 325=>0.3333*$e, 326=>0.3333*$e, 331=>0.25*$e, 332=>1-$e, 333=>0.25*$e, 335=>0.25*$e, 336=>0.25*$e, 351=>0.3333*$e, 352=>0.5*(1-$e), 353=>0.5*(1-$e), 355=>0.3333*$e, 356=>0.3333*$e, 361=>0.5*$e, 362=>0.25*(1-$e), 363=>0.5*(1-$e), 365=>0.25*(1-$e), 366=>0.5*$e, 511=>1-$e, 512=>0.25*$e, 513=>0.25*$e, 515=>0.25*$e, 516=>0.25*$e, 521=>0.5*(1-$e), 522=>0.5*(1-$e), 523=>0.3333*$e, 525=>0.3333*$e, 526=>0.3333*$e, 531=>0.25*$e, 532=>1-$e, 533=>0.25*$e, 535=>0.25*$e, 536=>0.25*$e, 551=>0.3333*$e, 552=>0.5*(1-$e), 553=>0.5*(1-$e), 555=>0.3333*$e, 556=>0.3333*$e, 561=>0.5*$e, 562=>0.25*(1-$e), 563=>0.5*(1-$e), 565=>0.25*(1-$e), 566=>0.5*$e, 611=>1-$e, 612=>0.25*$e, 613=>0.25*$e, 615=>0.25*$e, 616=>0.25*$e, 621=>0.5*(1-$e), 622=>0.5*(1-$e), 623=>0.3333*$e, 625=>0.3333*$e, 626=>0.3333*$e, 631=>0.25*$e, 632=>1-$e, 633=>0.25*$e, 635=>0.25*$e, 636=>0.25*$e, 651=>0.3333*$e, 652=>0.5*(1-$e), 653=>0.5*(1-$e), 655=>0.3333*$e, 656=>0.3333*$e, 661=>0.5*$e, 662=>0.25*(1-$e), 663=>0.5*(1-$e), 665=>0.25*(1-$e), 666=>0.5*$e);
	} elsif ($chry and $sample_sex->[2] eq 'male') {
		%oprior1 = (111=>1.000-$e, 112=>1.000*$e, 113=>0.500*$e, 115=>0.333*$e, 116=>0.250*$e, 121=>1.000-$e, 122=>1.000*$e, 123=>0.500*$e, 125=>0.333*$e, 126=>0.250*$e, 131=>1.000-$e, 132=>1.000*$e, 133=>0.500*$e, 135=>0.333*$e, 136=>0.250*$e, 151=>1.000-$e, 152=>1.000*$e, 153=>0.500*$e, 155=>0.333*$e, 156=>0.250*$e, 161=>1.000-$e, 162=>1.000*$e, 163=>0.500*$e, 165=>0.333*$e, 166=>0.250*$e, 211=>1.000*$e, 212=>1.000-$e, 213=>1.000*$e, 215=>0.500*$e, 216=>0.333*$e, 221=>1.000*$e, 222=>1.000-$e, 223=>1.000*$e, 225=>0.500*$e, 226=>0.333*$e, 231=>1.000*$e, 232=>1.000-$e, 233=>1.000*$e, 235=>0.500*$e, 236=>0.333*$e, 251=>1.000*$e, 252=>1.000-$e, 253=>1.000*$e, 255=>0.500*$e, 256=>0.333*$e, 261=>1.000*$e, 262=>1.000-$e, 263=>1.000*$e, 265=>0.500*$e, 266=>0.333*$e, 311=>0.500*$e, 312=>1.000*$e, 313=>1.000-$e, 315=>1.000*$e, 316=>0.500*$e, 321=>0.500*$e, 322=>1.000*$e, 323=>1.000-$e, 325=>1.000*$e, 326=>0.500*$e, 331=>0.500*$e, 332=>1.000*$e, 333=>1.000-$e, 335=>1.000*$e, 336=>0.500*$e, 351=>0.500*$e, 352=>1.000*$e, 353=>1.000-$e, 355=>1.000*$e, 356=>0.500*$e, 361=>0.500*$e, 362=>1.000*$e, 363=>1.000-$e, 365=>1.000*$e, 366=>0.500*$e, 511=>0.333*$e, 512=>0.500*$e, 513=>1.000*$e, 515=>1.000-$e, 516=>1.000*$e, 521=>0.333*$e, 522=>0.500*$e, 523=>1.000*$e, 525=>1.000-$e, 526=>1.000*$e, 531=>0.333*$e, 532=>0.500*$e, 533=>1.000*$e, 535=>1.000-$e, 536=>1.000*$e, 551=>0.333*$e, 552=>0.500*$e, 553=>1.000*$e, 555=>1.000-$e, 556=>1.000*$e, 561=>0.333*$e, 562=>0.500*$e, 563=>1.000*$e, 565=>1.000-$e, 566=>1.000*$e, 611=>0.250*$e, 612=>0.333*$e, 613=>0.500*$e, 615=>1.000*$e, 616=>1.000-$e, 621=>0.250*$e, 622=>0.333*$e, 623=>0.500*$e, 625=>1.000*$e, 626=>1.000-$e, 631=>0.250*$e, 632=>0.333*$e, 633=>0.500*$e, 635=>1.000*$e, 636=>1.000-$e, 651=>0.250*$e, 652=>0.333*$e, 653=>0.500*$e, 655=>1.000*$e, 656=>1.000-$e, 661=>0.250*$e, 662=>0.333*$e, 663=>0.500*$e, 665=>1.000*$e, 666=>1.000-$e);
	} elsif ($chry and $sample_sex->[2] eq 'female'){
		%oprior1 = (111=>$e*$e, 112=>$e*$e, 113=>$e*$e, 115=>$e*$e, 116=>$e*$e, 121=>$e*$e, 122=>$e*$e, 123=>$e*$e, 125=>$e*$e, 126=>$e*$e, 131=>$e*$e, 132=>$e*$e, 133=>$e*$e, 135=>$e*$e, 136=>$e*$e, 151=>$e*$e, 152=>$e*$e, 153=>$e*$e, 155=>$e*$e, 156=>$e*$e, 161=>$e*$e, 162=>$e*$e, 163=>$e*$e, 165=>$e*$e, 166=>$e*$e, 211=>$e*$e, 212=>$e*$e, 213=>$e*$e, 215=>$e*$e, 216=>$e*$e, 221=>$e*$e, 222=>$e*$e, 223=>$e*$e, 225=>$e*$e, 226=>$e*$e, 231=>$e*$e, 232=>$e*$e, 233=>$e*$e, 235=>$e*$e, 236=>$e*$e, 251=>$e*$e, 252=>$e*$e, 253=>$e*$e, 255=>$e*$e, 256=>$e*$e, 261=>$e*$e, 262=>$e*$e, 263=>$e*$e, 265=>$e*$e, 266=>$e*$e, 311=>$e*$e, 312=>$e*$e, 313=>$e*$e, 315=>$e*$e, 316=>$e*$e, 321=>$e*$e, 322=>$e*$e, 323=>$e*$e, 325=>$e*$e, 326=>$e*$e, 331=>$e*$e, 332=>$e*$e, 333=>$e*$e, 335=>$e*$e, 336=>$e*$e, 351=>$e*$e, 352=>$e*$e, 353=>$e*$e, 355=>$e*$e, 356=>$e*$e, 361=>$e*$e, 362=>$e*$e, 363=>$e*$e, 365=>$e*$e, 366=>$e*$e, 511=>$e*$e, 512=>$e*$e, 513=>$e*$e, 515=>$e*$e, 516=>$e*$e, 521=>$e*$e, 522=>$e*$e, 523=>$e*$e, 525=>$e*$e, 526=>$e*$e, 531=>$e*$e, 532=>$e*$e, 533=>$e*$e, 535=>$e*$e, 536=>$e*$e, 551=>$e*$e, 552=>$e*$e, 553=>$e*$e, 555=>$e*$e, 556=>$e*$e, 561=>$e*$e, 562=>$e*$e, 563=>$e*$e, 565=>$e*$e, 566=>$e*$e, 611=>$e*$e, 612=>$e*$e, 613=>$e*$e, 615=>$e*$e, 616=>$e*$e, 621=>$e*$e, 622=>$e*$e, 623=>$e*$e, 625=>$e*$e, 626=>$e*$e, 631=>$e*$e, 632=>$e*$e, 633=>$e*$e, 635=>$e*$e, 636=>$e*$e, 651=>$e*$e, 652=>$e*$e, 653=>$e*$e, 655=>$e*$e, 656=>$e*$e, 661=>$e*$e, 662=>$e*$e, 663=>$e*$e, 665=>$e*$e, 666=>$e*$e);
	} else {
		%oprior1 = (111=>1-$e, 112=>0.25*$e, 113=>0.25*$e, 115=>0.25*$e, 116=>0.25*$e, 121=>0.5*(1-$e), 122=>0.5*(1-$e), 123=>0.3333*$e, 125=>0.3333*$e, 126=>0.3333*$e, 131=>0.25*$e, 132=>1-$e, 133=>0.25*$e, 135=>0.25*$e, 136=>0.25*$e, 151=>0.3333*$e, 152=>0.5*(1-$e), 153=>0.5*(1-$e), 155=>0.3333*$e, 156=>0.3333*$e, 161=>0.5*$e, 162=>0.25*(1-$e), 163=>0.5*(1-$e), 165=>0.25*(1-$e), 166=>0.5*$e, 211=>0.5*(1-$e), 212=>0.5*(1-$e), 213=>0.3333*$e, 215=>0.3333*$e, 216=>0.3333*$e, 221=>0.25*(1-$e), 222=>0.5*(1-$e), 223=>0.25*(1-$e), 225=>0.5*$e, 226=>0.5*$e, 231=>0.3333*$e, 232=>0.5*(1-$e), 233=>0.5*(1-$e), 235=>0.3333*$e, 236=>0.3333*$e, 251=>0.5*$e, 252=>0.25*(1-$e), 253=>0.5*(1-$e), 255=>0.25*(1-$e), 256=>0.5*$e, 261=>$e, 262=>0.125*(1-$e), 263=>0.375*(1-$e), 265=>0.375*(1-$e), 266=>0.125*(1-$e), 311=>0.25*$e, 312=>1-$e, 313=>0.25*$e, 315=>0.25*$e, 316=>0.25*$e, 321=>0.3333*$e, 322=>0.5*(1-$e), 323=>0.5*(1-$e), 325=>0.3333*$e, 326=>0.3333*$e, 331=>0.25*$e, 332=>0.25*$e, 333=>1-$e, 335=>0.25*$e, 336=>0.25*$e, 351=>0.3333*$e, 352=>0.3333*$e, 353=>0.5*(1-$e), 355=>0.5*(1-$e), 356=>0.3333*$e, 361=>0.5*$e, 362=>0.5*$e, 363=>0.25*(1-$e), 365=>0.5*(1-$e), 366=>0.25*(1-$e), 511=>0.3333*$e, 512=>0.5*(1-$e), 513=>0.5*(1-$e), 515=>0.3333*$e, 516=>0.3333*$e, 521=>0.5*$e, 522=>0.25*(1-$e), 523=>0.5*(1-$e), 525=>0.25*(1-$e), 526=>0.5*$e, 531=>0.3333*$e, 532=>0.3333*$e, 533=>0.5*(1-$e), 535=>0.5*(1-$e), 536=>0.3333*$e, 551=>0.5*$e, 552=>0.5*$e, 553=>0.25*(1-$e), 555=>0.5*(1-$e), 556=>0.25*(1-$e), 561=>0.5*$e, 562=>0.5*$e, 563=>0.125*(1-$e), 565=>0.375*(1-$e), 566=>0.5*(1-$e), 611=>0.5*$e, 612=>0.25*(1-$e), 613=>0.5*(1-$e), 615=>0.25*(1-$e), 616=>0.5*$e, 621=>$e, 622=>0.125*(1-$e), 623=>0.375*(1-$e), 625=>0.375*(1-$e), 626=>0.125*(1-$e), 631=>0.5*$e, 632=>0.5*$e, 633=>0.25*(1-$e), 635=>0.5*(1-$e), 636=>0.25*(1-$e), 651=>0.5*$e, 652=>0.5*$e, 653=>0.125*(1-$e), 655=>0.375*(1-$e), 656=>0.5*(1-$e), 661=>0.5*$e, 662=>0.5*$e, 663=>0.0625*(1-$e), 665=>0.25*(1-$e), 666=>0.6875*(1-$e));
	}
	if ($chrx and $sample_sex->[3] eq 'male') {		#second child in the family
		%oprior2 = (111=>1-$e, 112=>0.25*$e, 113=>0.25*$e, 115=>0.25*$e, 116=>0.25*$e, 121=>0.5*(1-$e), 122=>0.5*(1-$e), 123=>0.3333*$e, 125=>0.3333*$e, 126=>0.3333*$e, 131=>0.25*$e, 132=>1-$e, 133=>0.25*$e, 135=>0.25*$e, 136=>0.25*$e, 151=>0.3333*$e, 152=>0.5*(1-$e), 153=>0.5*(1-$e), 155=>0.3333*$e, 156=>0.3333*$e, 161=>0.5*$e, 162=>0.25*(1-$e), 163=>0.5*(1-$e), 165=>0.25*(1-$e), 166=>0.5*$e, 211=>0.25*$e, 212=>1-$e, 213=>0.25*$e, 215=>0.25*$e, 216=>0.25*$e, 221=>0.3333*$e, 222=>0.5*(1-$e), 223=>0.5*(1-$e), 225=>0.3333*$e, 226=>0.3333*$e, 231=>0.25*$e, 232=>0.25*$e, 233=>1-$e, 235=>0.25*$e, 236=>0.25*$e, 251=>0.3333*$e, 252=>0.3333*$e, 253=>0.5*(1-$e), 255=>0.5*(1-$e), 256=>0.3333*$e, 261=>0.5*$e, 262=>0.5*$e, 263=>0.25*(1-$e), 265=>0.5*(1-$e), 266=>0.25*(1-$e), 311=>0.25*$e, 312=>0.25*$e, 313=>1-$e, 315=>0.25*$e, 316=>0.25*$e, 321=>0.3333*$e, 322=>0.3333*$e, 323=>0.5*(1-$e), 325=>0.5*(1-$e), 326=>0.3333*$e, 331=>0.25*$e, 332=>0.25*$e, 333=>0.25*$e, 335=>1-$e, 336=>0.25*$e, 351=>0.3333*$e, 352=>0.3333*$e, 353=>0.3333*$e, 355=>0.5*(1-$e), 356=>0.5*(1-$e), 361=>0.3333*$e, 362=>0.3333*$e, 363=>0.3333*$e, 365=>0.25*(1-$e), 366=>0.75*(1-$e), 511=>0.25*$e, 512=>0.25*$e, 513=>0.25*$e, 515=>1-$e, 516=>0.25*$e, 521=>0.3333*$e, 522=>0.3333*$e, 523=>0.3333*$e, 525=>0.5*(1-$e), 526=>0.5*(1-$e), 531=>0.25*$e, 532=>0.25*$e, 533=>0.25*$e, 535=>0.25*$e, 536=>1-$e, 551=>0.25*$e, 552=>0.25*$e, 553=>0.25*$e, 555=>0.25*$e, 556=>1-$e, 561=>0.25*$e, 562=>0.25*$e, 563=>0.25*$e, 565=>0.25*$e, 566=>1-$e, 611=>0.25*$e, 612=>0.25*$e, 613=>0.25*$e, 615=>0.25*$e, 616=>1-$e, 621=>0.25*$e, 622=>0.25*$e, 623=>0.25*$e, 625=>0.25*$e, 626=>1-$e, 631=>0.25*$e, 632=>0.25*$e, 633=>0.25*$e, 635=>0.25*$e, 636=>1-$e, 651=>0.25*$e, 652=>0.25*$e, 653=>0.25*$e, 655=>0.25*$e, 656=>1-$e, 661=>0.25*$e, 662=>0.25*$e, 663=>0.25*$e, 665=>0.25*$e, 666=>1-$e);
	} elsif ($chrx and $sample_sex->[2] eq 'female') {
		%oprior2 = (111=>1-$e, 112=>0.25*$e, 113=>0.25*$e, 115=>0.25*$e, 116=>0.25*$e, 121=>0.5*(1-$e), 122=>0.5*(1-$e), 123=>0.3333*$e, 125=>0.3333*$e, 126=>0.3333*$e, 131=>0.25*$e, 132=>1-$e, 133=>0.25*$e, 135=>0.25*$e, 136=>0.25*$e, 151=>0.3333*$e, 152=>0.5*(1-$e), 153=>0.5*(1-$e), 155=>0.3333*$e, 156=>0.3333*$e, 161=>0.5*$e, 162=>0.25*(1-$e), 163=>0.5*(1-$e), 165=>0.25*(1-$e), 166=>0.5*$e, 211=>1-$e, 212=>0.25*$e, 213=>0.25*$e, 215=>0.25*$e, 216=>0.25*$e, 221=>0.5*(1-$e), 222=>0.5*(1-$e), 223=>0.3333*$e, 225=>0.3333*$e, 226=>0.3333*$e, 231=>0.25*$e, 232=>1-$e, 233=>0.25*$e, 235=>0.25*$e, 236=>0.25*$e, 251=>0.3333*$e, 252=>0.5*(1-$e), 253=>0.5*(1-$e), 255=>0.3333*$e, 256=>0.3333*$e, 261=>0.5*$e, 262=>0.25*(1-$e), 263=>0.5*(1-$e), 265=>0.25*(1-$e), 266=>0.5*$e, 311=>1-$e, 312=>0.25*$e, 313=>0.25*$e, 315=>0.25*$e, 316=>0.25*$e, 321=>0.5*(1-$e), 322=>0.5*(1-$e), 323=>0.3333*$e, 325=>0.3333*$e, 326=>0.3333*$e, 331=>0.25*$e, 332=>1-$e, 333=>0.25*$e, 335=>0.25*$e, 336=>0.25*$e, 351=>0.3333*$e, 352=>0.5*(1-$e), 353=>0.5*(1-$e), 355=>0.3333*$e, 356=>0.3333*$e, 361=>0.5*$e, 362=>0.25*(1-$e), 363=>0.5*(1-$e), 365=>0.25*(1-$e), 366=>0.5*$e, 511=>1-$e, 512=>0.25*$e, 513=>0.25*$e, 515=>0.25*$e, 516=>0.25*$e, 521=>0.5*(1-$e), 522=>0.5*(1-$e), 523=>0.3333*$e, 525=>0.3333*$e, 526=>0.3333*$e, 531=>0.25*$e, 532=>1-$e, 533=>0.25*$e, 535=>0.25*$e, 536=>0.25*$e, 551=>0.3333*$e, 552=>0.5*(1-$e), 553=>0.5*(1-$e), 555=>0.3333*$e, 556=>0.3333*$e, 561=>0.5*$e, 562=>0.25*(1-$e), 563=>0.5*(1-$e), 565=>0.25*(1-$e), 566=>0.5*$e, 611=>1-$e, 612=>0.25*$e, 613=>0.25*$e, 615=>0.25*$e, 616=>0.25*$e, 621=>0.5*(1-$e), 622=>0.5*(1-$e), 623=>0.3333*$e, 625=>0.3333*$e, 626=>0.3333*$e, 631=>0.25*$e, 632=>1-$e, 633=>0.25*$e, 635=>0.25*$e, 636=>0.25*$e, 651=>0.3333*$e, 652=>0.5*(1-$e), 653=>0.5*(1-$e), 655=>0.3333*$e, 656=>0.3333*$e, 661=>0.5*$e, 662=>0.25*(1-$e), 663=>0.5*(1-$e), 665=>0.25*(1-$e), 666=>0.5*$e);
	} elsif ($chry and $sample_sex->[2] eq 'male') {
		%oprior2 = (111=>1.000-$e, 112=>1.000*$e, 113=>0.500*$e, 115=>0.333*$e, 116=>0.250*$e, 121=>1.000-$e, 122=>1.000*$e, 123=>0.500*$e, 125=>0.333*$e, 126=>0.250*$e, 131=>1.000-$e, 132=>1.000*$e, 133=>0.500*$e, 135=>0.333*$e, 136=>0.250*$e, 151=>1.000-$e, 152=>1.000*$e, 153=>0.500*$e, 155=>0.333*$e, 156=>0.250*$e, 161=>1.000-$e, 162=>1.000*$e, 163=>0.500*$e, 165=>0.333*$e, 166=>0.250*$e, 211=>1.000*$e, 212=>1.000-$e, 213=>1.000*$e, 215=>0.500*$e, 216=>0.333*$e, 221=>1.000*$e, 222=>1.000-$e, 223=>1.000*$e, 225=>0.500*$e, 226=>0.333*$e, 231=>1.000*$e, 232=>1.000-$e, 233=>1.000*$e, 235=>0.500*$e, 236=>0.333*$e, 251=>1.000*$e, 252=>1.000-$e, 253=>1.000*$e, 255=>0.500*$e, 256=>0.333*$e, 261=>1.000*$e, 262=>1.000-$e, 263=>1.000*$e, 265=>0.500*$e, 266=>0.333*$e, 311=>0.500*$e, 312=>1.000*$e, 313=>1.000-$e, 315=>1.000*$e, 316=>0.500*$e, 321=>0.500*$e, 322=>1.000*$e, 323=>1.000-$e, 325=>1.000*$e, 326=>0.500*$e, 331=>0.500*$e, 332=>1.000*$e, 333=>1.000-$e, 335=>1.000*$e, 336=>0.500*$e, 351=>0.500*$e, 352=>1.000*$e, 353=>1.000-$e, 355=>1.000*$e, 356=>0.500*$e, 361=>0.500*$e, 362=>1.000*$e, 363=>1.000-$e, 365=>1.000*$e, 366=>0.500*$e, 511=>0.333*$e, 512=>0.500*$e, 513=>1.000*$e, 515=>1.000-$e, 516=>1.000*$e, 521=>0.333*$e, 522=>0.500*$e, 523=>1.000*$e, 525=>1.000-$e, 526=>1.000*$e, 531=>0.333*$e, 532=>0.500*$e, 533=>1.000*$e, 535=>1.000-$e, 536=>1.000*$e, 551=>0.333*$e, 552=>0.500*$e, 553=>1.000*$e, 555=>1.000-$e, 556=>1.000*$e, 561=>0.333*$e, 562=>0.500*$e, 563=>1.000*$e, 565=>1.000-$e, 566=>1.000*$e, 611=>0.250*$e, 612=>0.333*$e, 613=>0.500*$e, 615=>1.000*$e, 616=>1.000-$e, 621=>0.250*$e, 622=>0.333*$e, 623=>0.500*$e, 625=>1.000*$e, 626=>1.000-$e, 631=>0.250*$e, 632=>0.333*$e, 633=>0.500*$e, 635=>1.000*$e, 636=>1.000-$e, 651=>0.250*$e, 652=>0.333*$e, 653=>0.500*$e, 655=>1.000*$e, 656=>1.000-$e, 661=>0.250*$e, 662=>0.333*$e, 663=>0.500*$e, 665=>1.000*$e, 666=>1.000-$e);
	} elsif ($chry and $sample_sex->[2] eq 'female'){
		%oprior2 = (111=>$e*$e, 112=>$e*$e, 113=>$e*$e, 115=>$e*$e, 116=>$e*$e, 121=>$e*$e, 122=>$e*$e, 123=>$e*$e, 125=>$e*$e, 126=>$e*$e, 131=>$e*$e, 132=>$e*$e, 133=>$e*$e, 135=>$e*$e, 136=>$e*$e, 151=>$e*$e, 152=>$e*$e, 153=>$e*$e, 155=>$e*$e, 156=>$e*$e, 161=>$e*$e, 162=>$e*$e, 163=>$e*$e, 165=>$e*$e, 166=>$e*$e, 211=>$e*$e, 212=>$e*$e, 213=>$e*$e, 215=>$e*$e, 216=>$e*$e, 221=>$e*$e, 222=>$e*$e, 223=>$e*$e, 225=>$e*$e, 226=>$e*$e, 231=>$e*$e, 232=>$e*$e, 233=>$e*$e, 235=>$e*$e, 236=>$e*$e, 251=>$e*$e, 252=>$e*$e, 253=>$e*$e, 255=>$e*$e, 256=>$e*$e, 261=>$e*$e, 262=>$e*$e, 263=>$e*$e, 265=>$e*$e, 266=>$e*$e, 311=>$e*$e, 312=>$e*$e, 313=>$e*$e, 315=>$e*$e, 316=>$e*$e, 321=>$e*$e, 322=>$e*$e, 323=>$e*$e, 325=>$e*$e, 326=>$e*$e, 331=>$e*$e, 332=>$e*$e, 333=>$e*$e, 335=>$e*$e, 336=>$e*$e, 351=>$e*$e, 352=>$e*$e, 353=>$e*$e, 355=>$e*$e, 356=>$e*$e, 361=>$e*$e, 362=>$e*$e, 363=>$e*$e, 365=>$e*$e, 366=>$e*$e, 511=>$e*$e, 512=>$e*$e, 513=>$e*$e, 515=>$e*$e, 516=>$e*$e, 521=>$e*$e, 522=>$e*$e, 523=>$e*$e, 525=>$e*$e, 526=>$e*$e, 531=>$e*$e, 532=>$e*$e, 533=>$e*$e, 535=>$e*$e, 536=>$e*$e, 551=>$e*$e, 552=>$e*$e, 553=>$e*$e, 555=>$e*$e, 556=>$e*$e, 561=>$e*$e, 562=>$e*$e, 563=>$e*$e, 565=>$e*$e, 566=>$e*$e, 611=>$e*$e, 612=>$e*$e, 613=>$e*$e, 615=>$e*$e, 616=>$e*$e, 621=>$e*$e, 622=>$e*$e, 623=>$e*$e, 625=>$e*$e, 626=>$e*$e, 631=>$e*$e, 632=>$e*$e, 633=>$e*$e, 635=>$e*$e, 636=>$e*$e, 651=>$e*$e, 652=>$e*$e, 653=>$e*$e, 655=>$e*$e, 656=>$e*$e, 661=>$e*$e, 662=>$e*$e, 663=>$e*$e, 665=>$e*$e, 666=>$e*$e);
	} else {
		%oprior2 = (111=>1-$e, 112=>0.25*$e, 113=>0.25*$e, 115=>0.25*$e, 116=>0.25*$e, 121=>0.5*(1-$e), 122=>0.5*(1-$e), 123=>0.3333*$e, 125=>0.3333*$e, 126=>0.3333*$e, 131=>0.25*$e, 132=>1-$e, 133=>0.25*$e, 135=>0.25*$e, 136=>0.25*$e, 151=>0.3333*$e, 152=>0.5*(1-$e), 153=>0.5*(1-$e), 155=>0.3333*$e, 156=>0.3333*$e, 161=>0.5*$e, 162=>0.25*(1-$e), 163=>0.5*(1-$e), 165=>0.25*(1-$e), 166=>0.5*$e, 211=>0.5*(1-$e), 212=>0.5*(1-$e), 213=>0.3333*$e, 215=>0.3333*$e, 216=>0.3333*$e, 221=>0.25*(1-$e), 222=>0.5*(1-$e), 223=>0.25*(1-$e), 225=>0.5*$e, 226=>0.5*$e, 231=>0.3333*$e, 232=>0.5*(1-$e), 233=>0.5*(1-$e), 235=>0.3333*$e, 236=>0.3333*$e, 251=>0.5*$e, 252=>0.25*(1-$e), 253=>0.5*(1-$e), 255=>0.25*(1-$e), 256=>0.5*$e, 261=>$e, 262=>0.125*(1-$e), 263=>0.375*(1-$e), 265=>0.375*(1-$e), 266=>0.125*(1-$e), 311=>0.25*$e, 312=>1-$e, 313=>0.25*$e, 315=>0.25*$e, 316=>0.25*$e, 321=>0.3333*$e, 322=>0.5*(1-$e), 323=>0.5*(1-$e), 325=>0.3333*$e, 326=>0.3333*$e, 331=>0.25*$e, 332=>0.25*$e, 333=>1-$e, 335=>0.25*$e, 336=>0.25*$e, 351=>0.3333*$e, 352=>0.3333*$e, 353=>0.5*(1-$e), 355=>0.5*(1-$e), 356=>0.3333*$e, 361=>0.5*$e, 362=>0.5*$e, 363=>0.25*(1-$e), 365=>0.5*(1-$e), 366=>0.25*(1-$e), 511=>0.3333*$e, 512=>0.5*(1-$e), 513=>0.5*(1-$e), 515=>0.3333*$e, 516=>0.3333*$e, 521=>0.5*$e, 522=>0.25*(1-$e), 523=>0.5*(1-$e), 525=>0.25*(1-$e), 526=>0.5*$e, 531=>0.3333*$e, 532=>0.3333*$e, 533=>0.5*(1-$e), 535=>0.5*(1-$e), 536=>0.3333*$e, 551=>0.5*$e, 552=>0.5*$e, 553=>0.25*(1-$e), 555=>0.5*(1-$e), 556=>0.25*(1-$e), 561=>0.5*$e, 562=>0.5*$e, 563=>0.125*(1-$e), 565=>0.375*(1-$e), 566=>0.5*(1-$e), 611=>0.5*$e, 612=>0.25*(1-$e), 613=>0.5*(1-$e), 615=>0.25*(1-$e), 616=>0.5*$e, 621=>$e, 622=>0.125*(1-$e), 623=>0.375*(1-$e), 625=>0.375*(1-$e), 626=>0.125*(1-$e), 631=>0.5*$e, 632=>0.5*$e, 633=>0.25*(1-$e), 635=>0.5*(1-$e), 636=>0.25*(1-$e), 651=>0.5*$e, 652=>0.5*$e, 653=>0.125*(1-$e), 655=>0.375*(1-$e), 656=>0.5*(1-$e), 661=>0.5*$e, 662=>0.5*$e, 663=>0.0625*(1-$e), 665=>0.25*(1-$e), 666=>0.6875*(1-$e));
	}

#dichotomize the CNV validation procedure:
#if there is no boundary discordance, do a simple Bayesian computation
#otherwise, use viterbi algorithm to infer the most likely path
	for my $nextregion (@region) {
		my ($nextchr, @nextpos) = @$nextregion;				#nextchr is the chromosome, nextname is the SNPs that constitute the CNV boundary
			my ($delta, $psi);						#delta: maximum prob at each time, psi: maximum previous ind that generate the delta at this time
			my $nextsubregion = join (",", $nextchr, $nextpos[0], $nextpos[1]);
		my ($maxval, $maxvalind);
		$verbose and print "processing regions at chr$nextchr pos=@nextpos\n";

#the following paragraph applies when there is no boundary discordance (a simple Bayesian method to calculate most likely a posterior state combintations)
		if (@nextpos == 2) {
			validateQuartetCNVCall ($ref_inputfile, \@fmprior, \%oprior1, $nextsubregion, $sample_sex, $region_alllogprob, $region_numsnp, $region_cnvlength, $region_chr, $region_name_start, $region_name_end, $region_pos_start, $region_pos_end);
			next;
		}

#proceed to the following only if there are multiple overlapped subregions in the region. We will a Viterbi algorithm to decode the best path
#delta is a matrix that record for each time point, for each fmostate, the maxlogprob of the path
#psi is a matrix that record for each time point, for each fmostate, the previous fmostate (that reach the current one with highest prob)

#step1: initialization
		for my $fstate (1, 2, 3, 5, 6) {
			for my $mstate (1, 2, 3, 5, 6) {
				my @ostate_array = generateOstateArray (@$ref_inputfile-2, [1,2,3,5,6]);

				for my $o12states (@ostate_array) {
					my ($o1state, $o2state) = @$o12states;
					my $fmoindex = $fstate+$mstate*7+$o1state*7*7+$o2state*7*7*7;

					$psi->[0][$fmoindex] = 0;
					$delta->[0][$fmoindex] += log ($fmprior[$fstate]) + log ($fmprior[$mstate]) + log ($oprior1{$fstate.$mstate.$o1state})+log ($oprior2{$fstate.$mstate.$o2state});
					$delta->[0][$fmoindex] += $region_alllogprob->{$nextsubregion}->[0][$fstate] + $region_alllogprob->{$nextsubregion}->[1][$mstate] + $region_alllogprob->{$nextsubregion}->[2][$o1state] + $region_alllogprob->{$nextsubregion}->[3][$o2state];
				}
			}
		}
		$verbose and print STDERR "NOTICE: initial scanning at $nextchr: @nextpos[0..1] found $maxvalind $maxval\n";

#step 2: recursion
		for my $i (1 .. @nextpos-2) {	#this is equivalent to "from time 1 to time t" (time 0 is the initialization stage)
			my $nextsubregion = join (",", $nextchr, $nextpos[$i], $nextpos[$i+1]);

			for my $fstate (1, 2, 3, 5, 6) {
				for my $mstate (1, 2, 3, 5, 6) {
					my @ostate_array = generateOstateArray (@$ref_inputfile-2, [1,2,3,5,6]);
					for my $o12states (@ostate_array) {
						my ($o1state, $o2state) = @$o12states;
						my $fmoindex = $fstate+7*$mstate+7*7*$o1state+7*7*7*$o2state;
						undef $maxval;

						for my $old_fstate (1, 2, 3, 5, 6) {
							for my $old_mstate (1, 2, 3, 5, 6) {
								my @old_ostate_array = generateOstateArray (@$ref_inputfile-2, [1,2,3,5,6]);
								for my $index (0 .. @old_ostate_array-1) {
									my $old_o12states = $old_ostate_array[$index];
									my ($old_o1state, $old_o2state) = @$old_o12states;
									my $old_fmoindex = $old_fstate+7*$old_mstate+7*7*$old_o1state+7*7*7*$old_o2state;
									$region_snpdist->{$nextsubregion} or confess "\nERROR: for nextregion $nextsubregion, pos=@nextpos @$ref_inputfile: no snpdist defined!!!";

#val is calculated as previous prob + transition prob of parent + transition prob of offspring
									my $val = $delta->[$i-1][$old_fmoindex];
									$val += log ($transition->[$old_fstate][$fstate] * (1-exp(-$region_snpdist->{$nextsubregion}/100_000)) / (1-exp(-5000/100_000))) + log ($transition->[$old_mstate][$mstate]* (1-exp(-$region_snpdist->{$nextsubregion}/100_000)) / (1-exp(-5000/100_000)));

									if ($fstate == $old_fstate and $mstate == $old_mstate and $o1state != $old_o1state and $oprior1{$fstate.$mstate.$o1state} > $e) {
										my @tempoprior = sort {$a<=>$b} @oprior1{$fstate.$mstate.'1', $fstate.$mstate.'2', $fstate.$mstate.'3', $fstate.$mstate.'5', $fstate.$mstate.'6'};
										$val += log ($tempoprior[0]);			#recombination happens without de novo event!
									} else {
										$val += log ($oprior1{$fstate.$mstate.$o1state});
									}
									if ($fstate == $old_fstate and $mstate == $old_mstate and $o2state != $old_o2state and $oprior2{$fstate.$mstate.$o2state} > $e) {
										my @tempoprior = sort {$a<=>$b} @oprior2{$fstate.$mstate.'1', $fstate.$mstate.'2', $fstate.$mstate.'3', $fstate.$mstate.'5', $fstate.$mstate.'6'};
										$val += log ($tempoprior[0]);			#recombination happens without de novo event!
									} else {
										$val += log ($oprior2{$fstate.$mstate.$o2state});
									}

#for the last block, arbitrarily set a "end" state that returns to state 3 (since all fstate and mstate return to 3, it can be directly added here as if there are two transitions)
									if ($i == @nextpos-2) {
										$val += log ($transition->[$fstate][3]) + log ($transition->[$mstate][3]);
									}

									if (not defined $maxval or $val > $maxval) {
										$maxval = $val;
										$maxvalind = $old_fmoindex;
									}
								}
							}
						}
						my $subpost = $region_alllogprob->{$nextsubregion}->[0][$fstate] + $region_alllogprob->{$nextsubregion}->[1][$mstate] + $region_alllogprob->{$nextsubregion}->[2][$o1state] + $region_alllogprob->{$nextsubregion}->[3][$o2state];
						$delta->[$i][$fmoindex] = $maxval + $subpost;
						$psi->[$i][$fmoindex] = $maxvalind;
					}
				}
			}
		}

		my ($pprob, $q, $qstate);	#pprob is the maximum prob, $q is the fmoindex that reach maximum prob, qstate is the [fstate,mstate,ostate]
#step 3: termination
			for my $fstate (1, 2, 3, 5, 6) {
				for my $mstate (1, 2, 3, 5, 6) {
					my @ostate_array = generateOstateArray (@$ref_inputfile-2, [1,2,3,5,6]);

					for my $o12states (@ostate_array) {
						my ($o1state, $o2state) = @$o12states;
						my $fmoindex = $fstate+7*$mstate+7*7*$o1state+7*7*7*$o2state;
						if (not defined $pprob or $delta->[@nextpos-2][$fmoindex] > $pprob) {
							$pprob = $delta->[@nextpos-2][$fmoindex];
							$q->[@nextpos-2] = $fmoindex;
							$qstate->[@nextpos-2] = [convertQuartetIndex2State ($fmoindex)];
						}
					}
				}
			}
		$verbose and print STDERR "q[", @nextpos-2, "]=", $q->[@$q-1], " (", convertTrioIndex2State ($q->[@$q-1]), ")\n";

#step 4: path backtracing
		for (my $i = @nextpos-3; $i >= 0; $i--) {
			$q->[$i] = $psi->[$i+1][$q->[$i+1]];
			$qstate->[$i] = [convertQuartetIndex2State ($q->[$i])];
			$verbose and print STDERR "q[$i] = $q->[$i] (", convertTrioIndex2State ($q->[$i]), ")\n";
		}

#finally, compile all state sequences together to get a concensus sequence
		my @cnv;			#@cnv has 4 element, corresponding to father, mother, offspring CNV calls
			for my $i (0 .. 3) {
				my ($stretch_start_j, $found_signal);
				my $normal_state = 3;
				$chrx and $sample_sex->[$i] eq 'male' and $normal_state = 2;		#for male chrx
					$chry and $sample_sex->[$i] eq 'male' and $normal_state = 2;		#for male chry
					$chry and $sample_sex->[$i] eq 'female' and $normal_state = 1;		#for female chry
					for my $j (0 .. @$qstate-1) {
						my $nextstate = $qstate->[$j][$i];
						if ($nextstate ne $normal_state) {
							if ($found_signal and $found_signal ne $nextstate) {	#transition between different CNV states
								push @{$cnv[$i]}, [$stretch_start_j, $j, $found_signal];
								$stretch_start_j = $j;
								$found_signal = $nextstate;
							} elsif ($found_signal) {
								1;						#do nothing, still in the same state
							} else {
								$found_signal = $nextstate;
								$stretch_start_j = $j;
							}
						} else {
							if ($found_signal) {
								push @{$cnv[$i]}, [$stretch_start_j, $j, $found_signal];
								$found_signal = 0;
							}
						}
					}
#finish the last stretch
				if ($found_signal) {
					push @{$cnv[$i]}, [$stretch_start_j, scalar (@$qstate), $found_signal];
				}
			}

		printFamilyCNVCall (\@cnv, $nextchr, \@nextpos, $region_name_start, $region_name_end, $qstate, $siginfo, $ref_inputfile);

	}
}

sub jointCNVCall {
	my ($ref_inputfile, $hmmfile, $pfbfile, $sexfile, $chrx, $chry, $gcmodelfile, $directory) = @_;
	my ($pfbinfo) = newreadPFB ($pfbfile);
	my $file_sex = {};
	$sexfile and $file_sex = readSexFile ($sexfile);		# if --sexfile is set, read gender information from sexfile for each sample
		my $hmm = main::readHMMFile ($hmmfile);				#check the validity of the HMM file
		my $region;

	my (%fmolrr, %fmobaf);						#key=chr value=ref_array (fat, mot, off)
		my (%name, %pos, %pfb, %snpdist, %probe_count);			#key=chr value=name, pos, pfb, snpdist, probe_count
		my (@trio_sample_sex, @trio_lrr_sd);
	my $logprob = 0;
	my @cnvcall = ();

	my $gcmodel;
	defined $gcmodelfile and $gcmodel = newreadGCModel ($gcmodelfile, $pfbinfo);

#read HMM model file
	my $hmm_model = khmm::ReadCHMM ($hmmfile);

	my $medianadjust = $main::medianadjust;
	my $sdadjust = $main::sdadjust;

#read intensity signals for all members
	for my $i (0 .. @$ref_inputfile-1) {
		my $inputfile = $ref_inputfile->[$i];
		my $sample_sex = $file_sex->{$inputfile} || 'unknown';	# if the sample is not in sexfile, or if it is 'unknown' in sex file
			my ($siginfo, $sigdesc);				# the hash contains all the signal intensity values for each chromosome (key)


			if ($chrx) {
				$region = 'X';
			} elsif ($chry) {
				$region = 'Y';
			} else {
				$region ||= "1-$lastchr";					# by default --region argument should be all autosomes
			}
		($siginfo, $sigdesc) = readLRRBAF ($inputfile, $region, $pfbinfo, $gcmodel, $directory);
#print STDERR "NOTICE: Descriptive summary of $sigdesc->{num_record} records: mean=", sprintf ("%.4f", $sigdesc->{lrr_mean}), " sd=", sprintf ("%.4f", $sigdesc->{lrr_sd}), "\n";
########################################################################################################################
		if ($chrx) {
			if ($sigdesc->{baf_xhet} eq 'NA') {
				print STDERR "ERROR: skipping inputfile $inputfile since it does not contain signal information for X chromosome markers.\n";
				return;
			}
			print STDERR "NOTICE: Descriptive summary of $sigdesc->{num_record} records in $inputfile: LRR_Xmean=", sprintf ("%.4f", $sigdesc->{lrr_xmean}), " LRR_Xmedian=", sprintf ("%.4f", $sigdesc->{lrr_xmedian}), " LRR_XSD=", sprintf ("%.4f", $sigdesc->{lrr_xsd}), " BAF_Xhet=", sprintf ("%.4f", $sigdesc->{baf_xhet}), "\n";

#examine whether sample_sex is correct or predict sample_sex when it is not specified in --sexfile
			if ($sigdesc->{baf_xhet} > 0.1) {		#if >10% of BAF are heterozygotes in chrX, it must be a female sample
				if (defined $sample_sex) {
					$sample_sex eq 'male' and print STDERR "WARNING: Sample sex for $inputfile is specified as 'male' in $sexfile but BAF heterozygous rate for chrX ($sigdesc->{baf_xhet}) indicates 'female'\n";
				} else {
					$sample_sex = 'female';
					print STDERR "NOTICE: Sample sex for $inputfile is predicted as 'female' based on BAF heterozygous rate for chrX ($sigdesc->{baf_xhet})\n";
				}
			} else {
				if (defined $sample_sex) {
					$sample_sex eq 'female' and print STDERR "WARNING: Sample sex for $inputfile is specified as 'female' in $sexfile but BAF heterozygous rate for chrX ($sigdesc->{baf_xhet}) indicates 'male'\n";
				} else {
					$sample_sex = 'male';
					print STDERR "NOTICE: Sample sex for $inputfile is predicted as 'male' based on BAF heterozygous rate for chrX ($sigdesc->{baf_xhet})\n";
				}
			}

			if ($sample_sex eq 'male') {			#for males, the mean LRR should be equal to expected LRR for state 2 (one-copy)
				print STDERR "NOTICE: Adjusting LRR values in chrX from $inputfile by ", sprintf ("%.4f", $sigdesc->{lrr_xmedian}-$hmm->{B1_mean}[1]), "\n";
				main::adjustLRR ($siginfo, $sigdesc->{lrr_xmedian}-$hmm->{B1_mean}[1]);
			} else {
				print STDERR "NOTICE: Adjusting LRR values in chrX from $inputfile by ", sprintf ("%.4f", $sigdesc->{lrr_xmedian}), "\n";
				main::adjustLRR ($siginfo, $sigdesc->{lrr_xmedian});
			}

#quality control: examine the variation of Log R Ratio values (>0.2 is generally treated as bad genotyping quality)
			if ($sigdesc->{lrr_xsd} > 0.2) {
				print STDERR "WARNING: Sample from $inputfile does not pass quality control criteria due to its large SD for LRR ($sigdesc->{lrr_xsd})!\n";
				print STDERR "WARNING: Small-sized CNV calls may not be reliable and should be interpreted with caution!\n";
			}
########################################################################################################################
########################################################################################################################
		} elsif ($chry) {
			if ($sigdesc->{baf_yhet} eq 'NA') {
				print STDERR "ERROR: skipping inputfile $inputfile since it does not contain signal information for Y chromosome markers.\n";
				return;
			}
			print STDERR "NOTICE: Descriptive summary of $sigdesc->{num_record} records in $inputfile: LRR_Ymean=", sprintf ("%.4f", $sigdesc->{lrr_ymean}), " LRR_Ymedian=", sprintf ("%.4f", $sigdesc->{lrr_ymedian}), " LRR_YSD=", sprintf ("%.4f", $sigdesc->{lrr_ysd}), "\n";

#examine whether sample_sex is correct or predict sample_sex when it is not specified in --sexfile
			if ($sigdesc->{baf_xhet} > 0.1) {		#if >10% of BAF are heterozygotes in chrX, it must be a female sample
				if (defined $sample_sex) {
					$sample_sex eq 'male' and print STDERR "WARNING: Sample sex for $inputfile is specified as 'male' in $sexfile but BAF heterozygous rate for chrX ($sigdesc->{baf_xhet}) indicates 'female'\n";
				} else {
					$sample_sex = 'female';
					print STDERR "NOTICE: Sample sex for $inputfile is predicted as 'female' based on BAF heterozygous rate for chrX ($sigdesc->{baf_xhet})\n";
				}
			} else {
				if (defined $sample_sex) {
					$sample_sex eq 'female' and print STDERR "WARNING: Sample sex for $inputfile is specified as 'female' in $sexfile but BAF heterozygous rate for chrX ($sigdesc->{baf_xhet}) indicates 'male'\n";
				} else {
					$sample_sex = 'male';
					print STDERR "NOTICE: Sample sex for $inputfile is predicted as 'male' based on BAF heterozygous rate for chrX ($sigdesc->{baf_xhet})\n";
				}
			}

			if ($sample_sex eq 'male') {			#for males, the mean LRR should be equal to expected LRR for state 2 (one-copy)
				print STDERR "NOTICE: Adjusting LRR values in chrY from $inputfile by ", sprintf ("%.4f", $sigdesc->{lrr_ymedian}-$hmm->{B1_mean}[1]), "\n";
				main::adjustLRR ($siginfo, $sigdesc->{lrr_ymedian}-$hmm->{B1_mean}[1]);
			} else {
#for females, the mean LRR should be equal to expected LRR for state 1 (zero-copy)
				print STDERR "NOTICE: Adjusting LRR values in chrY from $inputfile by ", sprintf ("%.4f", $sigdesc->{lrr_ymedian}-$hmm->{B1_mean}[0]), "\n";
				main::adjustLRR ($siginfo, $sigdesc->{lrr_ymedian}-$hmm->{B1_mean}[0]);
			}

#quality control: examine the variation of Log R Ratio values (>0.2 is generally treated as bad genotyping quality)
			if ($sigdesc->{lrr_ysd} > 0.2 and $sample_sex eq 'male') {
				print STDERR "WARNING: Sample from $inputfile does not pass quality control criteria due to its large SD for LRR ($sigdesc->{lrr_ysd})!\n";
				print STDERR "WARNING: Small-sized CNV calls may not be reliable and should be interpreted with caution!\n";
			}
########################################################################################################################
########################################################################################################################
		} else {
#quality control: examine the variation of Log R Ratio values (>0.2 is generally treated as bad genotyping quality)
			if ($sigdesc->{lrr_sd} > 0.2) {
				print STDERR "WARNING: Sample from $inputfile does not pass quality control criteria due to its large SD for LRR ($sigdesc->{lrr_sd})!\n";
				print STDERR "WARNING: Small-sized CNV calls may not be reliable and should be interpreted with caution!\n";
			}

#quality control: examine the median of BAF values (>0.6 or <0.4 is treated as bad clustering quality)
			if ($sigdesc->{baf_median} < 0.4 or $sigdesc->{baf_median} > 0.6) {
				print STDERR "WARNING: Sample from $inputfile does not pass quality control criteria due to its shifted BAF values (median=$sigdesc->{baf_median})!\n";
				print STDERR "WARNING: Small-sized CNV calls may not be reliable and should be interpreted with caution!\n";
			}

#quality control: examine the drifting of BAF values (>0.002 is generally treated as bad genotyping quality)
			if ($sigdesc->{baf_drift} > 0.002) {
				print STDERR "WARNING: Sample from $inputfile does not pass quality control criteria due to its drifting BAF values (drift=$sigdesc->{baf_drift})!\n";
				print STDERR "WARNING: Small-sized CNV calls may not be reliable and should be interpreted with caution!\n";
			}

			if ($medianadjust) {				#THIS IS A MANDATORY ARGUMENT NOW, SINCE IT ALWAYS IMPROVE PERFORMANCE!!! (use "--medianadjust 0" to disable this feature)
				if ($sigdesc->{num_record} >= 5_000) {			#just in case one is only interested in one small region (<5000 SNPs) that may already have large chunk of CNVs
					main::adjustLRR ($siginfo, $sigdesc->{lrr_median});
					print STDERR "NOTICE: Median-adjusting LRR values for all markers by ", sprintf ("%.4f", $sigdesc->{lrr_median}), "\n";
				}
			}
		}
		if ($i == 0) {
			$sample_sex eq 'female' and print STDERR "ERROR: sample_sex for $inputfile should male since it corresponds to father signals\n";
		} elsif ($i == 1) {
			$sample_sex eq 'male' and print STDERR "ERROR: sample_sex for $inputfile should be female since it corresponds to mother signals\n";
		}
		push @trio_sample_sex, $sample_sex;			# sample_sex is important information for CNV calling in chrX and chrY, but not for autosomes

			if ($sdadjust) {
				if ($chrx) {
					push @trio_lrr_sd, $sigdesc->{lrr_xsd};
				} elsif ($chry) {
					push @trio_lrr_sd, $sigdesc->{lrr_ysd};
				} else {
					push @trio_lrr_sd, $sigdesc->{lrr_sd};
				}
			} else {
				push @trio_lrr_sd, $hmm->{B1_sd}[2];
				print STDERR "NOTICE: Setting trio_lrr_sd as $hmm->{B1_sd}[2], effectively eliminating the sdadjust step in trio CNV detection\n";
			}

		for my $curchr (keys %$siginfo) {
			my $pos = $siginfo->{$curchr}{pos};
			my $lrr = $siginfo->{$curchr}{lrr};
			my $baf = $siginfo->{$curchr}{baf};
			@$pos >= 10 or print STDERR "WARNING: Skipping chromosome $curchr due to insufficient data points (<10 markers)!\n" and next;

			my @snpdist = @$pos;
			for my $i (1 .. @snpdist-2) {
				$snpdist[$i] = ($snpdist[$i+1]-$snpdist[$i]) || 1;	#sometimes two markers have the same chromosome location (in Affymetrix array annotation)
			}

			if ($snpdist{$curchr}) {
				if (@{$snpdist{$curchr}} != @snpdist) {
					confess "Error: unequal amounts of markers are read from $ref_inputfile->[0] (${\(scalar @{$snpdist{$curchr}})} markers) and $ref_inputfile->[$i] (${\(scalar @snpdist)} markers)";
				}
			}
			$snpdist{$curchr} ||= \@snpdist;
			$name{$curchr} ||= $siginfo->{$curchr}{name};
			$pfb{$curchr} ||= $siginfo->{$curchr}{pfb};
			$pos{$curchr} ||= $siginfo->{$curchr}{pos};
			push @{$fmolrr{$curchr}}, $lrr;
			push @{$fmobaf{$curchr}}, $baf;
			$probe_count{$curchr} ||= @snpdist-1;
		}
	}

#generate CNV calls
	for my $curchr (sort keys %name) {
		print STDERR "NOTICE: Calling CNVs in chromosome $curchr with $probe_count{$curchr} markers\n";
		khmm::testVitTrio_CHMM ($hmm_model, $probe_count{$curchr}, $fmolrr{$curchr}->[0], $fmobaf{$curchr}->[0], $fmolrr{$curchr}->[1], $fmobaf{$curchr}->[1], $fmolrr{$curchr}->[2], $fmobaf{$curchr}->[2], $pfb{$curchr}, $snpdist{$curchr}, \$logprob, \@trio_lrr_sd, $trio_sample_sex[2] eq 'male'?1:2, $chrx);
		my @trioid= qw/father mother offspring/;
		for my $i (0 .. 2) {
			analyzeTrioStateSequence (\@cnvcall, $i, $curchr, $snpdist{$curchr}, $name{$curchr}, $pos{$curchr}, $trio_sample_sex[$i], $trioid[$i]);
		}
	}
	for my $i (0 .. 2) {
		if ($tabout) {
			printTabbedCNV ($cnvcall[$i], $ref_inputfile->[$i]);
		} else {
			printFormattedCNV ($cnvcall[$i], $ref_inputfile->[$i]);
		}
	}
}

sub analyzeTrioStateSequence {
	my ($triocnvcall, $curindex, $curchr, $symbol, $name, $pos, $sample_sex, $trioid) = @_;
	my ($normal_state, $found_signal, $stretch_start_i);

#set up the normal state in case of sex chromosome analysis
	if ($curchr eq 'X' and $sample_sex eq 'male') {
		$normal_state = 2;
	} elsif ($curchr eq 'Y') {
		if ($sample_sex eq 'male') {
			$normal_state = 2;
		} else {
			$normal_state = 1;
		}
	} else {
		$normal_state = 3;
	}

	$triocnvcall->[$curindex] ||= {};						#cnv is a reference to a hash (key=state value=array of info)
		my (@curstate, @triopath, $currentdn);
	for my $i (1 .. @$symbol-1) {							#the first element is zero (arrays start from 1 in khmm module)
		$currentdn = int ($symbol->[$i] / 125);
		$curstate[0] = int (($symbol->[$i] - $currentdn*125) / 25);
		$curstate[1] = int (($symbol->[$i] - $currentdn*125 - $curstate[0]*25) / 5);
		$curstate[2] = $symbol->[$i] - $currentdn*125 - $curstate[0]*25 - $curstate[1]*5;

		if ($i <= 10) {
#print "i=$i dn=$currentdn state=@curstate\n";
		}

		for my $j (0 .. 2) {
			$curstate[$j]++;
			$curstate[$j] >= 4 and $curstate[$j]++;		#convert CN to state
		}

		my $curstate = $curstate[$curindex];

#found a new CNV or continue within a previously identified stretch of CNV
		if ($curstate != $normal_state) {
			if ($found_signal and $found_signal ne $curstate) {		#transition to one CNV to another CNV with different copy number
				push @{$triocnvcall->[$curindex]{$found_signal}}, [$curchr, $pos->[$stretch_start_i], $pos->[$i-1], $i-$stretch_start_i, $name->[$stretch_start_i], $name->[$i-1], $trioid, join ('-', @triopath)];
				$stretch_start_i = $i;
				$found_signal = $curstate;
				@triopath = (join ('', @curstate));
			} elsif ($found_signal) {					#do nothing, continue within a previously identified stretch of CNV
				$triopath[$#triopath] eq join ('', @curstate) or push @triopath, join ('', @curstate);
			} else {							#found the start of a new CNV
				$found_signal = $curstate;
				$stretch_start_i = $i;
				@triopath = (join ('', @curstate));			#start a new trio path (each node in the path is the trio state combination at a block of SNPs)
			}
		} else {
			if ($found_signal) {						#end a previously identified CNV
				push @{$triocnvcall->[$curindex]{$found_signal}}, [$curchr, $pos->[$stretch_start_i], $pos->[$i-1], $i-$stretch_start_i, $name->[$stretch_start_i], $name->[$i-1], $trioid, join ('-', @triopath)];
				$found_signal = 0;
			}
		}
	}
	if ($found_signal) {								#finish the last CNV
		push @{$triocnvcall->[$curindex]{$found_signal}}, [$curchr, $pos->[$stretch_start_i], $pos->[@$symbol-1], @$symbol-$stretch_start_i, $name->[$stretch_start_i], $name->[@$symbol-1], $trioid, join ('-', @triopath)];
	}
}

#this is a completely empirical method to exclude possible heterosomic regions (chromosomes) from the CNV calling results. It usually results in over-prediction and may exclude real rearrangements!
#Therefore, it is always important to check the signal intensities in visualization software to make sure that there are indeed in vitro chromosome deletion or duplication, since usually you can observe clear non-canonical clustering patterns of B Allelel Frequency and slightly increased or decreased log R Ratio values.
sub excludeHeterosomic {
	my ($cnvfile) = @_;
	my (%cnvchrcount, %cnvchrlength);
	open (CNV, $cnvfile) or confess "\nERROR: cannot read from cnvfile $cnvfile: $!";
	while (<CNV>) {
		m/^chr(\w+):(\d+)-(\d+)\s+numsnp=(\d+).+?state(\d)(?:,cn=\d)?\s+(\S+)\s+startsnp=(\S+)\s+endsnp=(\S+)/ or confess "\nERROR: invalid record found in cnvfile $cnvfile: <$_>";
		$5 == 4 and next;
		$4 > 5 or next;
		$cnvchrcount{$6, $1, $5}++;
		$cnvchrlength{$6, $1, $5} += ($3-$2+1);
	}
	open (CNV, $cnvfile) or confess "\nERROR: cannot read from cnvfile $cnvfile: $!";
	while (<CNV>) {
		m/^chr(\w+).+?state(\d)(?:,cn=\d)?\s+(\S+)\s+startsnp=(\S+)\s+endsnp=(\S+)/ or confess "\nERROR: invalid record found in cnvfile $cnvfile: <$_>";
		$cnvchrcount{$3, $1, $2} and $cnvchrcount{$3, $1, $2} == -1 and next;		#this has already been processed
			if ($cnvchrcount{$3, $1, $2} and ($cnvchrcount{$3, $1, $2} >= $heterosomic_threshold and $cnvchrlength{$3, $1, $2} > 1_000_000 or $cnvchrcount{$3, $1, $2} >= 40)) {
				print STDERR "NOTICE: chr$1 from $3 is excluded due to possible problems with heterosomic abberations ($cnvchrcount{$3, $1, $2} CNVs $cnvchrlength{$3, $1, $2} bases found in one chromosome)\n";
				$cnvchrcount{$3, $1, $2} = -1;
			} else {
				print;
			}
	}
	close (CNV);
}

sub calWF {
	my ($siginfo) = @_;
	my ($wf, $cc);
	my (%chrsig, @allsig, @alllrr, @allblock);
	my (@cursignal, %chryi, @yi);

	for my $nextchr (keys %$siginfo) {
		for my $i (0 .. @{$siginfo->{$nextchr}{pos}}-1) {
			push @{$chrsig{$nextchr}}, [$siginfo->{$nextchr}{pos}[$i], $siginfo->{$nextchr}{lrr}[$i]];
		}
	}

	for my $nextchr (keys %chrsig) {
		$nextchr =~ m/^\d+$/ or next;					#only consider autosome in WF calculation
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

	$verbose and print STDERR "NOTICE: Detected ${\(scalar @yi)} 1MB sliding windows in autosomes\n";
	@yi = grep {!m/NA/} @yi;
	if (not @yi) {
		print STDERR "(WARNING: cannot calculate waviness factor) ...";
		return ('NA', 'NA');
	}
	my $yi_median = median (\@yi);
#print STDERR " (${\(scalar @yi)} sliding windows were used in WF calculation with median=$yi_median)\n";
	@yi = map {abs ($_ - $yi_median)} @yi;
	$wf = median (\@yi);

	if (not exists $chryi{$refchr}) {
		print STDERR "(WARNING: cannot calculate waviness factor) ...";
		return ('NA', 'NA');
	}
	$cc = calCC ([@ref_median[0 .. @{$chryi{$refchr}}-1]], [@{$chryi{$refchr}}]);
	$cc eq 'NA' and return ($wf, 'NA');
	$cc > 0 and $wf = -$wf;
	return ($wf, $cc);
}

sub readSexFile {
	my ($inputfile) = @_;
	my (%file_sex);
	my %sexcode = (1=>'male', 2=>'female', 'm'=>'male', 'f'=>'female', 'male'=>'male', 'female'=>'female');
	open (SEX, $inputfile) or confess "\nERROR: cannot read from sexfile $sexfile: $!\n";
	while (<SEX>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		@record == 2 or confess "\nERROR: invalid record found in sexfile: <$_> (2 tab-delimited fields expected)";
		$record[1] = lc $record[1];
		$sexcode{$record[1]} eq 'male' or $sexcode{$record[1]} eq 'female' or confess "\nERROR: unrecognizable sex code in sexfile $sexfile: $record[1]\n";
		$file_sex{$record[0]} = $sexcode{$record[1]};
	}
	return (\%file_sex);
}

sub sortChr {
	my ($a1, $b1) = @_;
	if ($a1 =~ m/^(\d+)$/) {
		if ($b1 =~ m/^(\d+)$/) {
			return $a1<=>$b1;
		} else {
			return -1;
		}
	} elsif ($b1 =~ m/^(\d+)$/) {
		return 1;
	} else {
		return $a1 cmp $b1;
	}
}

sub adjustLRR {
	my ($siginfo, $adjustment) = @_;
	for my $curchr (keys %$siginfo) {
		map {$_-=$adjustment} @{$siginfo->{$curchr}{lrr}};
	}
}

sub adjustBAF {
	my ($siginfo, $adjustment) = @_;
	for my $curchr (keys %$siginfo) {
		map {$_>0.25 and $_<0.75 and $_-=$adjustment} @{$siginfo->{$curchr}{baf}};
	}
}

sub calCC {
	my ($array1, $array2) = @_;
	my (@newarray1, @newarray2);
	@$array1 == @$array2 or print STDERR "WARNING: Unequal dimensions of arrays: ${\(scalar @$array1)} vs ${\(scalar @$array2)}\n";
	for my $i (0 .. @$array1-1) {
		defined $array1->[$i] or next;
		defined $array2->[$i] or next;
		$array1->[$i] eq 'NA' and next;
		$array2->[$i] eq 'NA' and next;
		push @newarray1, $array1->[$i];
		push @newarray2, $array2->[$i];
	}
	$verbose and print STDERR "NOTICE: CC calculated on ", scalar (@newarray1), " elements\n";
	@newarray1 <= 2 and return 'NA';
	return cc (\@newarray1, \@newarray2);
}

#the following subroutine calculates the correlation coefficient
sub cc {
	my ($score1, $score2) = @_;
	my ($ssr12, $ssr11, $ssr22) = (ssr ($score1, $score2), ssr ($score1, $score1), ssr ($score2, $score2));
	$ssr11*$ssr22 or return "NA";
	return $ssr12 / sqrt ($ssr11 * $ssr22);
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

sub sum {
	my ($score) = @_;
	@$score or confess "\nERROR: NO VALUES for calculating sum\n";
	my $sum;
	for (@$score) {
		$sum += $_;
	}
	return $sum;
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

sub sd {
	my ($score) = @_;
	@$score > 1 or confess "\nERROR: NO sufficient VALUES for calculating SD\n";
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
	@$score or confess "\nERROR: NO VALUES for calculating median\n";
	my @newscore = sort {$a<=>$b} @$score;
	if (@newscore % 2 == 0) {
		return ($newscore[@newscore/2-1]+$newscore[@newscore/2])/2;
	} else {
		return $newscore[@newscore/2];
	}
}

=head1 SYNOPSIS

 detect_cnv.pl [arguments] <inputfile [inputfile2] ... | -listfile file>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation

        Analysis Type:
 	    --train			train optimized HMM model (not recommended to use)
 	    --test			test HMM model to identify CNV
 	    --wgs			test HMM model to identify CNV from transformed wgs data
 	    --trio			posterior CNV calls for father-mother-offspring trio
 	    --quartet			posterior CNV calls for quartet
 	    --joint			joint CNV calls for trio
 	    --cctest			case-control comparison of per-marker CNV frequency
 	    --validate			validate copy number at a pre-specified region

        Input/output Files:
 	    --listfile <file>		a list file containing path to files to be processed
 	    --output <file>		specify output root filename
 	    --logfile <file>		write notification/warning messages to this file
 	    --hmmfile <file>		HMM model file
 	    --pfbfile <file>		population frequency for B allelel file
 	    --sexfile <file>		a 2-column file containing filename and sex (male/female) for chrx CNV calling
 	    --cnvfile <file>		specify CNV call file for use in family-based CNV calling by -trio or -quartet
 	    --directory <string>	specify the directory where signal files are located
 	    --refchr <string>		specify a chromosome for wave adjustment (default: 11 for human)
 	    --refgcfile <file>		a file containing GC percentage of each 1M region of the chromosome specified by --refchr (default: 11 for human)
 	    --gcmodelfile <file>	a file containing GC model for wave adjustment
 	    --phenofile <file>		a file containing phenotype information for each input file for -cctest operation

        CNV output control:
 	    --minsnp <int>		minimum number of SNPs within CNV (default=3)
 	    --minlength <int>		minimum length of bp within CNV
 	    --minconf <float>		minimum confidence score of CNV
 	    --confidence		calculate confidence for each CNV
 	    --chrx			use chrX-specific treatment
 	    --chry			use chrY-specific treatment (available soon)

        Validation-calling arguments:
 	    --startsnp <string>		start SNP of a pre-specified region for --validate operation
 	    --endsnp <string>		end SNP of a pre-specified region for --validate operation
 	    --delfreq <float>		prior deletion frequency of a pre-specified region for --validate operation
 	    --dupfreq <float>		prior duplication frequency of a pre-specified region for --validate operation
 	    --backfreq <float>		background CNV probability for any loci (default: 0.0001)
 	    --candlist <file>		a file containing all candidate CNV regions to be validated
 	    
        Misc options
 	    --loh			detect copy-neutral LOH (obselete argument; for SNP arrays only!)
 	    --exclude_heterosomic	empirically exclude CNVs in heterosomic chromosomes
 	    --fmprior <numbers>		prior belief on CN state for regions with CNV calls
 	    --denovo_rate <float>	prior belief on genome-wide de novo event rate (default=0.0001)
 	    --tabout			use tab-delimited output
 	    --coordinate_from_input	get marker coordindate information from input signal file
 	    --control_label <string>	the phenotype label for control subjects in the phenotype file (default=control)
 	    --onesided			performed one-sided test for --cctest operation
 	    --type_filter <dup|del>	used together with --cctest to specify types of CNVs to be tested
 	    --(no)medianadjust		adjust genome-wide LRR such that median=0 (default=ON)
 	    --(no)bafadjust		adjust genome-wide BAF such that median=0.5 (default=ON)
 	    --(no)sdadjust		adjust SD of hidden Markov model based on input signal (default=ON)
 	    --(no)flush			flush input/output buffer (default=ON)
 	    --bafxhet <float>		minimum BAF het rate to predict female gender when -sexfile is not supplied (default=0.1)
 	    --lastchr <integer>		the number of the last autosomal chromosome, change for non human data (default=22).

 Function: generate CNV calls from high-density SNP genotyping data that 
 contains Log R Ratio and B Allele Frequency for each SNP or CN marker. Use -m 
 argument to read the complete manual.

 Example: detect_cnv.pl -test -hmm hhall.hmm -pfb hhall.pfb file1 file2 file3 -log logfile -out outfile
          detect_cnv.pl -trio -hmm hhall.hmm -pfb hhall.pfb -cnv outfile file1 file2 file3
          detect_cnv.pl -validate -hmm hhall.hmm -pfb hhall.pfb -startsnp rs100 -endsnp rs200 -delfreq 0.2 file1 file2 file3

 Version: $LastChangedDate: 2013-02-08 11:10:50 -0800 (Fri, 08 Feb 2013) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--train>

train an optimized HMM model, given a HMM model file, a PFB file, and several 
signal intensity files.

=item B<--ctrain>

this argument is NOT intended to be directly usable by the user. It trains a HMM 
model by calling a C subroutine directly, without reading data in Perl, thus 
providing a less memory-intensive way for training HMM model (compared to 
previous versions of the program).

=item B<--test>

generate CNV calls given a list of signal intensity files, a HMM model file and 
a PFB file.

=item B<--wgs>

generate CNV calls given a list of signal intensity files, a HMM model file and
a PFB file, which are transformed from whole-genome sequence data using PennCNV-Seq protocols.

=item B<--trio>

generate CNV calls for a father-mother-offspring trio, given a CNV file 
containing calls generated on each individual separately, a HMM model file, a 
PFB file, and the three signal intensity files.

=item B<--quartet>

jointly generate CNV calls for a father-mother-offspring1-offpspring2 quartet, 
given a CNV file containing calls generated on each individual separately, a HMM 
model file, a PFB file, and the four signal intensity files.

=item B<--joint>

New in July 2008: generate CNV calls for a father-mother-offspring trio via a 
one-step procedure. It is considerably slower than the --trio argument, but 
generates more accurate CNV calls with reduced false negative rates in 
simulation studies. One major issue for joint-calling is that it is far more 
memory intensive than trio-calling algorithm; nevertheless, it works well for a 
550K array in a 64-bit computer with 8GB memory. For arrays with much higher 
density, users are advised to divide PFB files into multiple files for the CNV 
calling (for example, each PFB file may contain 1 million consecutive markers).

=item B<--summary>

generate summary statistics on signal quality for each input file. Usually the 
summary is provided when calling CNVs and can be written to a log file via the 
--log argument; however, sometimes users forget to use --log, such that the 
signal quality information is lost. The --summary argument can calculate the 
signal quality again quickly without calling CNVs.

=item B<--cctest>

perform a case-control test on the frequency of having CNVs for each marker 
within CNVs. A separate phenotype file must be specified via the --phenofile 
argument for this to work. The actual test is a two-sided Fisher exact test. The 
--onesided argument can be specified for performning one-sided test, and the 
--type_filter argument can be specified so that only "dup" or "del" is compared 
between cases and controls.

=item B<--validate>

generate validation-based CNV calls. Unlike HMM algorithm that attempts to 
identify consecutive markers with copy number changes, in validation-based CNV 
calling, the user must specify a specific region (start-end markers), and the 
program attempts to identify the most likely copy number for all markers in the 
region, assuming that all markers have identical copy number.

=item B<--listfile>

a file containg a list of file names (path to other files) to be processed, one 
per line (for --test operation) or three per line (for --trio operation) or four 
per line (for --quartet operation). Normally you can also list all files to be 
processed in the command line; however, due to the system/shell limitations, you 
can only list maybe ~1000 file names in command line in most Linux systems (use 
ulimit command to check this). Using --listfile bypass this limitation so that 
many more files can be processed to give CNV calls.

As of July 2008, the file name seperator for list file changes from space to tab 
characters. Therefore, one can use spaces in the file name (such as "C:\program 
files\data\signal file 1") in the list file, and they can be processed correctly.

As of July 2008, the file names can be a piped output from a system command, if 
the file name starts with "`" and ends with "`". See B<Advanced program options: 
input filtering> below for more details.

=item B<--output>

specify an output file for CNV calls. By default, the CNV calls will be written 
to STDOUT.

=item B<--logfile>

specify a file name so the notification/warning messages from the program are 
written to this file (in addition to being printed in STDERR). The file can be 
further parsed to figure out the information for each file (such as LRR_SD 
values). Use the filter_cnv.pl program in PennCNV package for automatic analysis 
of logfile.

=item B<--hmmfile>

specify a HMM model file containing elements necessary for specifying the hidden 
Markov model for CNV calling.

=item B<--pfbfile>

a population frequency of B allele file containing chromosome coordinates of 
each SNP, as well as the frequency of B allele in a large reference population 
for this SNP.

In CNV calling, the Chr and Position in PFB will be used in annotating the 
position of SNPs. This makes it easy for users to generate CNV calls using 
different reference genome assemblies, without changing the input signal 
intensity files. However, when --coord_from_input is specified, the Chr and 
Position information for all markers will be read from the input signal 
intensity file instead.

=item B<--sexfile>

a 2-column file containing filename and sex (male/female) for sex chromosome 
calling with -chrx argument. The first tab-delimited column should be the input 
signal file name, while the second tab-delimited column should be male or 
female. Alternatively, abbreviations including m (male), f (female), 1 (male) or 
2 (female) are also fine.

=item B<--cnvfile>

a file containing CNV calls, that could be generated by the --test operation of 
this program, and be used subsequently in --trio operation or --quartet 
operation of this program.

=item B<--directory>

specify a directory to read signal files from. The user can either give the full 
path to signal files (in the command line or a list file), or only give file 
names but use --directory argument to specify where to find the files. For 
example, if the file name is inputdir/subdir/inputfile1, and the directory is 
dir1/dir2, the program will read signal intensity from the file 
dir1/dir2/inputdir/subdir/inputfile1, which is created by concatenating the dir 
and the file names.

=item B<--gcmodelfile>

a file that contains the GC percentage in the genomic region surrounding each 
marker for the GC-model based signal adjustment.

=item B<--phenofile>

a file containing phenotype information for each individual, so that --cctest 
can be used to compare the frequency between cases and controls. Each line has 
two tab-delimited fields: file name and the phenotype. By default, "control" 
means control subjects, and other words means cases; however, the user can use 
--control_label argument to change the phenotype label for controls.

=item B<--minsnp>

the minimum number of SNPs that a CNV call must contain to be in output

=item B<--minlength>

the minimum length of base pairs that a CNV call must contain to be in output

=item B<--minconf>

minimum confidence score for a CNV call to be in output. This is an experimental 
feature, and the actual definition of "confidence score" may change in the 
future.

=item B<--confidence>

specify that confidence scores be calculated for each CNV call in the -test 
operation. The confidence score is calculated as the log likelihood of the 
called copy number state minus the second most likely copy number state. This is 
an experimental feature, and the actual definition of "confidence score" may 
change in the future.

=item B<--chrx>

process chromosome X specifically. By default only autosomes will be processed 
by this program. If the species is not human, the PFB file must be changed that 
use X and Y to specify the two sex chromosomes. Otherwise, the program will not 
understand which chromosome is the sex chromosome.

=item B<--chry>

process chromosome Y specifically. By default only autosomes will be processed 
by this program. If the species is not human, the PFB file must be changed that 
use X and Y to specify the two sex chromosomes. Otherwise, the program will not 
understand which chromosome is the sex chromosome.

=item B<--startsnp>

specify the start SNP of a pre-specified region used in --validate operation

=item B<--endsnp>

specify the end SNP of a pre-specified region used in --validate operation

=item B<--delfreq>

specify the prior deletion allele frequency of a pre-specified region used in 
--validate operation (this frequency can be estimated from CNV calls by --test 
operation)

=item B<--dupfreq>

specify the prior duplication allele frequency of a pre-specified region used in 
--validate operation (this frequency can be estimated from CNV calls by --test 
operation)

=item B<--backfreq>

background CNV probability for any loci, with default value as 0.0001. This 
argument is useful in validation calling. When -delfreq/-dupfreq is not 
specified, the background frequency is used to calculate the prior probability 
of different copy number states.

=item B<--candlist>

specify a candidate region list file. In the file, multiple candidate regions, 
and their a prior deletion frequency, duplication frequency can be specified, 
and the -validate operation will attemp to generate CNV calls on these specific 
genomic regions, based on prior information on deletion/duplication frequencies.

=item B<--loh>

output copy-neutral LOH CNV calls. This is an obselete feature; one can 
alternatively use SNP genotype data for the LOH inference. This argument does 
NOT work for arrays with non-polymorphic markers.

=item B<--exclude_heterosomic>

exclude CNV calls from chromosomes showing evidence of heterosomic abberations 
from a given file containing CNV calls. An purely empirical method is applied in 
this procedure, although I recommended always manually examine the patterns of 
BAF to determine whether heterosomic abberation is present in a particular 
sample, if the sample size is relatively small (<100).

=item B<--fmprior>

the prior probability of 6 hidden states a given CNV call in father or mother. 
This is used for joint calling of trios or quartets. It is specified as six 
numbers separated by a comma. The empirically derived default values actually 
work well.

=item B<--denovo_rate>

specify the probability that a given CNV is a de novo event for family-based CNV 
calling. The default is 0.0001. (In previous version this rate is set as 0.01)

=item B<--tabout>

use tab-delimited output format. This is useful for parsing the output by third-
party programs.

=item B<--coordinate_from_input>

specifies that the genome coordinate information for the markers are retrieved 
from the input signal file, rather than the PFB file (default).

=item B<--control_label>

specify the text label for control subjects in the phenotype file specified by 
the --phenofile argument. Normally the "control" is used to specify controls, 
and all other individuals are treated as cases. However, some times users may 
use 1 to denote controls and 2 to denote cases; in such situations the 
"--control_label 1" should be used for the --cctest operation.

=item B<--onesided>

perform one-sided Fisher exact test in the --cctest operation.

=item B<--type_filter>

specify the particular types of CNVs to be used in the --cctest operation. By 
default both duplications and deletions are treated as a single group of CNVs 
and be used to compare cases and controls.

=item B<--(no)medianadjust>

this option is turned on by default. It adjust the log R Ratio 
values of the entire genome by a constant so that the median is zero.

=item B<--(no)bafadjust>

this option is turned ON by default (new July 2008): it adjust the BAF values 
genome-wide such that the median value is 0.5.

=item B<--(no)sdadjust>

this option is turned ON by default: it adjust the SD values in HMM model such 
that the model fits the signal quality of the testing sample to reduce false 
positive calls

=item B<--(no)flush>

this argument is turned ON by default. It requires the input/output buffer to 
flush immediately (that is, no input/output is buffered). When PennCNV is 
running remotely (for example, through a SSH connection) or when the output is 
redirected, this argument cause the program to report progress in real-time. 
When running PennCNV in parallel with many processes accessing disks 
simultaneously, this option should be turned off to decrease system overhead.

=item B<--bafxhet>

this argument specifies the BAF heterozygosity rate in chrX to predict the sex 
for a sample. Note that this rate is based on BAF values so it is not genotype 
heterozygosity rate and indeed quite different/smaller than that genotype 
heterozygosity rate. By default if >10% chrX markers have BAF values around 0.5, 
the sample is predicted as female. For chrX CNV calling, rather than relying on 
PennCNV prediction of gender, it is always best to explicitely specify the 
sample sex using the -sexfile argument.

=back

=head1 DESCRIPTION

This program is designed for detecting Copy Number Variation (CNV) based on Log 
R Ratio and B Allele Frequency data from high-density SNP genotyping platform, 
such as Illumina SNP arrays, Affymetrix SNP arrays, Perlegen SNP arrays, as well 
as custom-made SNP arrays. Users have also reported success in analyzing Agilent 
tiling arrays without SNP markers. Please cite: Wang K, Li M, Hadley D, Liu R, 
Glessner J, Grant S, Hakonarson H, Bucan M. PennCNV: an integrated hidden Markov 
model designed for high-resolution copy number variation detection in whole-
genome SNP genotyping data. Genone Research, 17:1665-1674, 2007.

=over 8

=item * B<Background>

Copy number variations (CNVs) refers to genomic segments of at least 1 kb in 
size, for which copy number differences have been observed in comparison to 
reference genome assemblies. Several high-throughput technologies, including 
array-CGH with large-insert clones, SNP genotyping arrays, whole-genome 
oligonucleotide arrays, paired-end sequencing, next-generation resequencing and 
so on, have been developed to detect CNVs in mammalian genome. The detect_cnv.pl 
program in PennCNV package is originally designed for inferring CNVs from marker 
signal intensities in Illumina SNP genotyping arrays, but has been extended to 
other types of arrays.

To use PennCNV, the user needs to supply Log R Ratio (LRR) and B Allele 
Frequency (BAF) values for all SNPs. (If Agilent/Nimblegen array is used, the BAF values 
still need to be specified, but it can be a random number like zero.) These two 
measures were originally developed by Illumina but can be applied to other 
platform as well. The program outputs the chromosome region of the CNV calls, 
the copy number of calls, the start and end SNP in the CNV as well as some other 
auxillary information such as the confidence score and the input file name. For 
SNP genotyping platforms, PennCNV can give calls on the exact copy number 
ranging from zero to four, while higher copy number cannot be inferred from the 
signal intensity confidently. The CNV calls can be further processed by other 
programs in the PennCNV package; for example, by the visualize_cnv.pl program to 
display the calls in genome browser, by filter_cnv.pl program to filter out CNV 
calls meeting particular criteria, or by the scan_region.pl program to infer the 
neighboring genes for CNVs.

=item * B<File formats>

Several input files are required for CNV calling: signal intensity files (which 
need to be generated by the user), a population frequency file (which is 
provided together this the program package, such as as hh550.hg17.pfb and 
hh550.hg18.pfb files, corresponding to genome coordinates for the 2004 and 2006 
human genome assembly for Illumina HumanHap550 arrays, respectively), a HMM 
model file (which is included in the program package).

I used a similar (but longer and more complicated) file format for HMM model 
file as the one used in the UMDHMM package. Actually it is a little weird to do 
that, since many entries in the file have either no meaning at all (such as 
emission probabilities), or has very different meanings (such as transition 
probabilities). Nevertheless, I still decide to use a similar format, solely to 
indicate that this program adopted the UMDHMM framework. A sample HMM model 
file is listed below:

	M=6
	N=6
	A:
	0.932302 0.004004 0.013696 0.050027 0.000010 0.000010 
	0.000010 0.917243 0.067506 0.013988 0.001292 0.000010 
	0.000011 0.000031 0.986449 0.013531 0.000015 0.000013 
	0.000018 0.000021 0.039392 0.960596 0.000012 0.000011 
	0.000010 0.000011 0.026750 0.000071 0.973197 0.000011 
	0.000010 0.000011 0.059087 0.006112 0.000012 0.934818 
	B:
	0.950000 0.000001 0.050000 0.000001 0.000001 0.000001 
	0.000001 0.950000 0.050000 0.000001 0.000001 0.000001 
	0.000001 0.000001 1.000000 0.000001 0.000001 0.000001 
	0.000001 0.000001 0.050000 0.950000 0.000001 0.000001 
	0.000001 0.000001 0.050000 0.000001 0.950000 0.000001 
	0.000001 0.000001 0.050000 0.000001 0.000001 0.950000 
	pi:
	0.001000 0.001000 0.995000 0.001000 0.001000 0.001000 
	B1_mean:
	-3.93000 -0.675382 0.000000 0.000000 0.361171 0.680540 
	B1_sd:
	1.000000 0.328682 0.222929 0.222929 0.244476 0.274889 
	B1_uf:
	0.010000
	B2_mean:
	0.000000 0.250000 0.333333 0.500000 
	B2_sd:
	0.016371 0.046211 0.050623 0.044834 
	B2_uf:
	0.010000

The signal intensity file is a tab-delimited file containing signal values 
(LRR,BAF) for each SNP (one SNP per line). The first a few lines for a signal 
intensity file is shown below:

	Name    Chr     Position        sample1.Allele Calls  sample1.B Allele Freq sample1.Log R Ratio
	rs3094315	1	792429	AG	0.375	0.340
	rs12562034	1	808311	GG	1.000	-0.008
	rs3934834	1	1045729	GG	1.000	-0.016
	rs9442372	1	1058627	AG	0.599	-0.155
	rs3737728	1	1061338	AG	0.555	-0.525
	rs11260588	1	1061581	GG	1.000	0.176
	rs6687776	1	1070488	GG	1.000	-0.055
	rs9651273	1	1071463	AG	0.494	0.012
	rs4970405	1	1088878	AA	0.000	-0.089


Each line should contain at least the information on SNP name, the LRR and the 
BAF. The relative positions of columns are determined by the annotation in the 
first line of the file (so-called header line). So the following file will be 
also recognized just fine:

	Name	sample1.Log R Ratio	sample1.B Allele Freq
	rs3094315	0.340	0.375
	rs12562034	-0.008	1.000
	rs3934834	-0.016	1.000
	rs9442372	-0.155	0.599
	rs3737728	-0.525	0.555
	rs11260588	0.176	1.000
	rs6687776	-0.055	1.000
	rs9651273	0.012	0.494
	rs4970405	-0.089	0.000

The population frequency file (PFB) is also a tab-delimited file containing 
frequency values, chromosome, and positions for each SNP (one SNP per line). The 
first a few lines for a sample PFB file is shown below:

	Name    Chr     Position        PFB
	rs3094315	1	792429	0.168
	rs12562034	1	808311	0.866
	rs3934834	1	1045729	0.808
	rs9442372	1	1058627	0.559
	rs3737728	1	1061338	0.748
	rs11260588	1	1061581	0.994
	rs6687776	1	1070488	0.811
	rs9651273	1	1071463	0.740
	rs4970405	1	1088878	0.132

Each line should contain SNP name, chr, position and PFB values, separated by 
tab. Besides providing the marker allele frequency information, the PFB file 
also specifies the chromosome coordinates of markers in a particular array. 
Therefore, to generate CNV calls for different genome assemblies, one can simply 
switch to a different PFB file without regenerating signal intensity files. If 
the user has his/her own genome assembly (different from NCBI human genome 
assembly), the user must supply the coordinate information in the input signal 
intensity file, then use --coord_from_input argument when running PennCNV.

=item * B<Program arguments: an overview>

For a simple help message, the -h or --help argument can be used; for the 
complete manual on program usage, the -m or --man argument can be used.

Both single dash and double dash in the argument name in the command line can be 
accepted, so -test and --test have the same effect. In addition, you can even 
omit some letters from the argument, as long as there is no ambiguity. For 
example, --test and --tes and --te have the same effect.

=item * B<Program usage: individual-based CNV calling>

The most typical usage of the program is to give CNV calls for a new sample, 
given a HMM model file and a PFB file. With the provided example, you can do:

	detect_cnv.pl -test -hmm hh550.hmm -pfb hh550.hg18.pfb -log outfile.log -out outfile file1.txt file2.txt file3.txt

The above command process three files, including file1.txt, file2.txt and 
file3.txt and generates CNV calls for them. The hh550.hmm and hh550.hg18.pfb 
file are constructed using the Illumina HumanHap550 array for the 2006 human 
genome assembly. The --test argument tells the program to call CNVs from the 
files. Note that if other types of Illumina arrays (such as HumanHap300, 
HumanHap1M and CNV370 arrays) are used, the hhall.hmm and hhall.hg18.pfb file 
should be used. When Affymetrix arrays are used, you can ues array-specific HMM 
and PFB files in the PennCNV-Affy package (for example, the affygw6.hmm and 
affygw6.hg18.pfb files for genome-wide 6.0 arrays).

Users can put multiple file names in the command line for CNV calling, or can 
also generate a listfile, with one file name per line, and use the --list 
argument to process each file and generate CNV calls.

To generate CNV calls specifically for chrX, the user need to add the --chrx 
argument. Without this argument, only CNVs in autosomes will be detected!

=item * B<Program usage: family-based posterior CNV calling>

After generating individual-based CNV calls, one may use family relationship to 
re-call CNVs and reconcile the boundary discordance in CNV calls in family 
members. This algorithm is described in the original PennCNV paper in 2007. For 
example, suppose that the father, mother and offspring CNV call information are 
all stored in the rawcnvfile file. The following command:

	detect_cnv.pl -trio -hmm hh550.hmm -pfb hh550.hg18.pfb -log newout.log -out newout -cnv rawcnvfile father.txt mother.txt child.txt

The signal intensity files for father, mother and child need to be specified in 
the command line consecutively for the trio-based CNV calling. The results will 
be written to the newout file.

=item * B<Program usage: trio-based joint CNV calling>

As of July 2008, a new trio-based joint-calling algorithm is implemented in 
the PennCNV package as beta version. Unlike the posterior-calling algorithm that 
use two steps, the joint calling algorithm is a one-step algorithm that directly 
generate CNV calls for father, mother and offspring simultaneously, but the 
speed of the program is considerably slower than the posterior calling 
algorithm. The advantage is that the joint-calling algorithm has lower false 
negative rates, and better power in resolving the correct CNV boundaries for 
trios.

To use the joint-calling algorithm, issue the following command: 

	detect_cnv.pl -joint -hmm hh550.hmm -pfb hh550.hg18.pfb -log newout.log -out newout father.txt mother.txt child.txt

In the above command, the --joint argument specifies that the joint-calling 
algorithm should be used for CNV calling, and the three input files are 
specified in the command line. (Alternatively, a list file can be specified by 
the --listfile argument, and the file contains three file names for one trio in 
each line.)

=item * B<Program usage: case-control comparisons>

The generated CNV calls for two phenotype groups can be used to test whether one 
group (such as cases) have more CNVs at a particular marker than the other group 
(such as controls). The --cctest operation provides a convenient way of doing 
this comparison, however, the functionality is rudimentary, and does not 
consider the various signal qualities of different sample groups. This operation 
therefore is purely used for discovery purposes and should NOT be used to report 
association. The P-values are generated by two-sided Fisher exact test, unless 
--onesided argument is specified. The users can use --type_filter to test only 
on duplications or deletions, since by default all CNVs are used in the 
comparison between cases and controls.

To use the case-control comparison, the user should first run generate CNV calls 
for all individuals and save the calls into a file such as outfile. The user 
need to create another file that specify the phenotype for each individual: in 
the file, each line contains a file name and a phenotype label separated by tab 
character.

The following command can then be issued:

	detect_cnv.pl --cctest -cnv outfile -phenofile phenotype -pfb pfbfile

For each marker that was located in a CNV region in the outfile, a result will 
be printed out. It tells the number of cases with CNV in this marker, number of 
cases without CNV in this marker, number of controls with CNV in this marker and 
number of controls without CNV in this marker, and also the P-values from Fisher 
exact test.

Again these P-values and the comparisons should be used for discovery purposes 
only. Different signal qualities and different batches of genotyping could 
result in artificially high level of CNVs at particular regions in cases but not 
in controls.

=item * B<Advanced program options: GC-model based signal adjustment>

As of July 2008, a new GC-model adjustment procedure is implemented in the 
detect_cnv.pl program as well as the genomic_wave.pl program. This method is 
used to handle genotyped samples that show wavy patterns in signal intensities. 
Our analysis demonstrate that these wavy patterns usually correlate with 
genome-wide GC distributions, such that a simple regression model based on GC 
content can be used for adjustment of signal intensities to increase the 
signal-to-noise ratio. For eac array, a separate GC-model file is provided.

For example, to perform an individual-based CNV calling for the file1.txt file:

	detect_cnv.pl -test -hmm hh550.hmm -pfb hh550.hg18.pfb -log outfile.log -out outfile file1.txt -gcmodel hh550.hg18.gcmodel

In the above command, the --gcmodel argument is used to specify the particular 
GC-model file to be used in the signal adjustment. For other types of Illumina 
arrays, one can use the hhall.* files.

=item * B<Advanced program options: empirical signal adjustment>

As of July 2008, several default arguments are set on in the signal 
pre-processing for CNV calling. The --bafadjust, --medianadjust and the 
--sdadjust are all set on as default. When they are in effect, the input file 
signal intensities will be adjusted such that the median of autosome BAF values 
will be 0.5, the median LRR values will be 0.5, and the SD of LRR/BAF values for 
normal copy number in the HMM model will match that in the input file. These 
options usually improve the program performance and reduce false positive CNV 
calls, so they are turned on by default. However, if the user is dealing with a 
regional array (such as a custom array from the GoldenGate platform, these 
options should be probabably turned off).

For custome-made arrays on particular genomic regions, these options should be 
turned OFF! This is because genome-wide data is NOT available to assess how to 
adjust the signal.

=item * B<Advanced program options: input filtering>

Although most times users create their own input signal files and then ask 
PennCNV to process these files sequentially, there may be situations where the 
signals are not stored in the desired format (for example, the signal 
information may be stored in a MySQL server). The PennCNV --list argument 
supports input filtering, such that piped input from a system command can be 
read as signal intensities and be processed. This functionality eliminate the 
need to create temporary files, so this method avoids using excessive disk 
storage for large-scale projects, and this method avoids all the overhead 
associated with disk writting and reading to achieve the fastest CNV calling 
speed, especially when running PennCNV on dozens of cluster machines 
simultaneously.

For example, suppose the user writes a program called retrieve-signal that 
retrieves the desired signal intensity for a given identifier from a database, 
then print to STDOUT. Now for CNV calling, instead of writting the signal to a 
file and then analyze the file, the user can simply put these lines into a list 
file:

	`retrieve-signal sample1 databaseA userB`
	`retrieve-signal sample2 databaseC userD`

then using the --list argument, PennCNV automatically recognize that these 
commands are not real file names, but piped output from some system commands, 
and can read from the pipe directly.

It is important to note that these commands above should reproduce the same 
output as a real signal file, that is, with Name, Log R Ratio and B Allele 
Frequency in the header line, and with the actual values in subsequent lines.

It is important to note that the PennCNV output will contain the system command, 
rather than a signal file name. This can be seen by the "`" character in the CNV 
calls. For example:

	chr20:51031721-51052942       numsnp=4      length=21,222      state2,cn=1 `extract_lrr_baf /mnt/1b/call/20070302/1773195392_A.int` startsnp=rs6022196 endsnp=rs6013594

Some text processing can be used to remove the "`" character and keep only a 
file name or sample identifier in the CNV output.

=item * B<Advanced program options: output filtering>

The --minsnp, --minlength and --minconf arguments can be used to filter out some 
CNV calls such that the output contains a subset of more confident calls. By 
default the --minsnp threshold is set as 3, while no --minlength or --minconf is 
in effect.

It is recommended to just leave these argument as is so that many CNV calls are 
generated. One can then post-process the CNV calls by custom scripts, or by 
using the filter_cnv.pl program provided in the PennCNV package to get a more 
confident list of calls for subsequent analysis.

=item * B<Advanced program options: handling genome coordinates>

The PennCNV program uses a PFB (population frequency of B allele) file for 
providing the genome coordinates for markers intended to be used in CNV calling. 
The advantage is that by switching to different PFB files, one can generate CNV 
calls for different genome corodinates (the most popular ones are 2004 and 2006 
human genome assembly). In the PennCNV package, I adapted the nomenclature of 
UCSC Genome Browser for file naming: I use hg17 and hg18 to denote 2004 and 
2006 human genome assembly, respectively, and may use similar nomenclature for 
other mammalian genomes in the future.

The user should be careful when using the PFB files, and make sure that the file 
corresponds to the particular array that the user is interested in. A very 
common mistake for many users is to use hh550.hg18.pfb file (designed for 
Illumina HumanHap550 array) on Illumina 1M chips. Essentially only half the 
markers in the 1M chip can be used for CNV inference in this case. In the latest 
version of PennCNV, I have not provided the hhall.* files, that should work on 
all Illumina pre-made arrays, including HumanHap300, Human1, Human610, 
HumanHap550, HumanHap650, HumanHap1M.

There are cases when one wants to use his or her own genome assembly, without 
resorting to PFB files. In such cases, the --coordinate_from_input argument can 
be used, which specifies that the genome coordinates for each marker is read 
from the input file. The input file must contain the Chr and the Position column 
in this case, otherwise an error message will be issued. (When using the 
BeadStudio plug-in version of PennCNV, the --coordinate_from_input argument is 
turned on by default.)

=item * B<Version history>

2010may01 version: 

=over 4

reduce memory usage by rewritting PFB-related subroutines

rewrite the genomic wave adjustment subroutine to solve compatibility issues 
during compilation at some systems

=back

2009Aug27 version:

=over 4

added --lastchr argument for handling non-human arrays

slight changes to C code to make it more efficient

=back

2008Nov19 version: 

=over 4

eliminated the --wavemodelfile and --sample_index arguments

change the architecture of multiple subroutines (should not affect CNV call output)

better argument error handling (if user supply useless argument, a warning message will be printed)

=back

2008jun26 version:

=over 4

beta version of GC-model based signal pre-processing to handle low-quality 
samples showing genomic waves

beta version of joint-calling algorithm for family data (slower than 
posterior-calling but more accurate)

changed default arguments: -bafadjust is turned ON by default; -denovo is 
changed to 0.0001 by default

provided the hhall.hmm, hhall.hg18.pfb and hhall.hg18.gcmodel file for all 
commonly used Illumina arrays

the previous wave-adjustment procedure (--wavemodel) is now obselete!!!

added genomic_wave.pl program for generating new signal files with signal 
adjustment

added filter_cnv.pl program for filtering CNV calls meeting user-specified 
criteria

=back

2008mar11 version: 

=over 4

improved version of BeadStudio plug-in for generating CNV calls from BeadStudio

beta-version of CNV visualization functions

=back

2008mar03 version: beta-version for using PennCNV with BeadStudio plug-in, 
beta-version for signal adjustment for wavy samples (via --wavemodel argument)

2008feb15 version: minor changes with more QC summary

2007dec14 version: support 64bit system and support cygwin now, print sample 
quality measure, produce confidence score by --conf argument (experimental 
feature)

2007nov13 version: cumulative minor bug fix and function enhancement

2007oct31 version: fix bugs for chrX processing; fix bugs in kcolumn.pl for 
splitting huge files

2007sep28 version: reimplement calling algorithm to process each chromosome 
separately; change output format

=item * B<Acknowledgements>

The software framework is adapted from the UMDHMM software package 
(http://www.kanungo.us/software/software.html), which has excellent 
implementation of the framework for Forward, Backward, Viterbi and Baum-Welch 
algorithm for discrete state HMMs. All the codes are rewritten, since the 
PennCNV uses continuous likelihood calculation, that is, there is no discrete 
emission probabilities for each HMM state in PennCNV. However, I tried to keep 
the names of the subroutines/variables similar to the original package, to help 
readers understand and rewrite the source code, if they need to.

The PennCNV package is mostly written in Perl, except that the core HMM 
computation subroutines are written in C solely for speed purposes. I have 
indeed implemented the HMM in Perl for the --trio and the --quartet arguments 
that uses family information for CNV calling. However, for the --joint argument 
that generates family-based CNV calls, the code is again implemented in C, again 
purely due to speed reasons.

The development of PennCNV is a collaborative project that involves Maja Bucan 
in Penn, Mingyao Li in Penn, Hakon Hakonarson in CHOP. Many users of PennCNV 
have provided many valuable feedbacks on the program, and I appreciate these 
comments for making the program better.

All code in the PennCNV package (except those under the include/ directory) is 
copiable, distributable, modifiable, and usable without any restrictions and 
without further permission from me.

Please contact Kai (kai@openbioinformatics.org) for questions or concerns on the 
PennCNV program.

=back

=cut                                                                                                                                                                                                                                                                                                                                                                 
