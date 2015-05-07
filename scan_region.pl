#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2011-01-17 01:16:59 -0800 (Mon, 17 Jan 2011) $';

our ($verbose, $help, $man);
our ($queryfile, $dbfile);
our ($snp_flag, $score_threshold, $normscore_threshold, $overlap, $dbregion, $append, $condense_query, $mce_flag, $phastcons_flag, $evofold_flag, $tfbs_flag, 
	$wgrna_flag, $refgene_flag, $segdup_flag, $knowngene_flag, $test_flag, $anno_flag, $refcds_flag, $refexon_flag,
	$autoexpand, $expandleft, $expandright, $expandmax, $bothside, $strand, $minoverlap, $kgxref, $reflink, $name2_flag, $quiet, $maxoverlap, $expanddb,
	$maxquerydbratio, $minquerydbratio, $mindbfrac, $maxdbfrac, $minqueryfrac, $maxqueryfrac, $queryinfo);
	
GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'snp_flag'=>\$snp_flag, 'score_threshold=f'=>\$score_threshold, 
	'normscore_threshold=i'=>\$normscore_threshold, 'overlap'=>\$overlap, 'dbregion'=>\$dbregion,
	'append'=>\$append, 'condense_query'=>\$condense_query, 
	'mce_flag'=>\$mce_flag, 'phastcons_flag'=>\$phastcons_flag, 'evofold_flag'=>\$evofold_flag, 'tfbs_flag'=>\$tfbs_flag, 'wgrna_flag'=>\$wgrna_flag, 
	'refgene_flag'=>\$refgene_flag, 'segdup_flag'=>\$segdup_flag, 'knowngene_flag'=>\$knowngene_flag, 'test_flag'=>\$test_flag, 'anno_flag'=>\$anno_flag,
	'refcds_flag'=>\$refcds_flag, 'refexon_flag'=>\$refexon_flag,
	'expandleft=s'=>\$expandleft, 'expandright=s'=>\$expandright, 'expandmax=s'=>\$expandmax, 'bothside'=>\$bothside, 'expanddb=s'=>\$expanddb,
	'strand=s'=>\$strand, 'kgxref=s'=>\$kgxref, 'reflink=s'=>\$reflink, 'name2_flag'=>\$name2_flag, 'quiet'=>\$quiet,
	'minoverlap=f'=>\$minoverlap, 'maxoverlap=f'=>\$maxoverlap, 'maxquerydbratio=f'=>\$maxquerydbratio, 'minquerydbratio=f'=>\$minquerydbratio, 
	'mindbfrac=f'=>\$mindbfrac, 'maxdbfrac=f'=>\$maxdbfrac, 'minqueryfrac=f'=>\$minqueryfrac, 'maxqueryfrac=f'=>\$maxqueryfrac,
	'queryinfo'=>\$queryinfo) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");

($queryfile, $dbfile) = @ARGV;

#preparing default parameters
$snp_flag and pod2usage ("Error in argument: the --snp_flag argument is in development");
$strand and $strand =~ m/^(plus|minus|forward|reverse)$/ || pod2usage ("Error in argument: the --strand argument can be set as only plus (forward) or minus (reverse)");

$minoverlap ||= 0;
$expandleft ||= 0; $expandright ||= 0; $expandmax ||= 0;
if ($expandmax) {
	$expandmax =~ s/k$/000/i; $expandmax =~ s/m$/000000/i;
	$expandmax =~ m/^\d+$/ or pod2usage ("Error in argument: the --expandmax argument should be a positive integer (suffix of k and m are acceptable)");
	$expandleft || $expandright and pod2usage ("Error in argument: please do not specify --expandleft or --expandright when --expandmax is specified");
}
if ($expandleft) {
	$expandleft =~ s/k$/000/i; $expandleft =~ s/m$/000000/i;
	$expandleft =~ m/^\d+$/ or pod2usage ("Error in argument: the --expandleft argument should be a positive integer (suffix of k and m are acceptable)");
}
if ($expandright) {
	$expandright =~ s/k$/000/i; $expandright =~ s/m$/000000/i;
	$expandright =~ m/^\d+$/ or pod2usage ("Error in argument: the --expandright argument should be a positive integer (suffix of k and m are acceptable)");
}
if ($expanddb) {
	$expanddb =~ s/k$/000/i; $expanddb =~ s/m$/000000/i;
	$expanddb =~ m/^\d+$/ or pod2usage ("Error in argument: the --expanddb argument should be a positive integer (suffix of k and m are acceptable)");
}
$name2_flag and $refgene_flag || pod2usage ("Error in argument: please do not specify --name2_flag unless --refgene_flag is specified");
$kgxref and $knowngene_flag || pod2usage ("Error in argument: please do not specify --kgxref unless --knowngene_flag is specified");
$reflink and $refgene_flag || $refexon_flag || $refcds_flag || pod2usage ("Error in argument: please do not specify --reflink unless --refgene_flag/refexon_flag/refcds_flag is specified");
$queryinfo and $dbregion || pod2usage "Error in argument: the --qinfo argument can be used only if --dbregion is specified";

if ($minoverlap or $maxoverlap or $minquerydbratio or $maxquerydbratio or $mindbfrac or $maxdbfrac or $minqueryfrac or $maxqueryfrac) {
	if ($mce_flag or $evofold_flag or $tfbs_flag or $wgrna_flag or $segdup_flag or $anno_flag or $knowngene_flag or $refgene_flag or $refcds_flag or $refexon_flag or $phastcons_flag) {
		pod2usage ("Error in argument: the --minoverlap, --maxoverlap, --minquerydbratio, --maxquerydbratio, --mindbfrac, --maxdbfrac, --minqueryfrac, --maxqueryfrac arguments are applicable only to plain DB-file format");
	}
	my @arg = sort {$a<=>$b} ($minoverlap||0, $maxoverlap||0, $minquerydbratio||0, $mindbfrac||0, $maxdbfrac||0, $minqueryfrac||0, $maxqueryfrac||0);
	$arg[0] >=0 and $arg[$#arg]<=1 or pod2usage ("Error in argument: the --minoverlap, --maxoverlap, --minquerydbratio, --maxquerydbratio, --mindbfrac, --maxdbfrac, --minqueryfrac, --maxqueryfrac arguments must be between 0 and 1");
}

if ($score_threshold or $normscore_threshold) {
	if ($mce_flag or $evofold_flag or $tfbs_flag or $wgrna_flag or $segdup_flag or $anno_flag) {
		pod2usage ("Error in argument: the --score_threshold or --normscore_threshold arguments are applicable only to selected UCSC annotation formats");
	}
}

if ($expandmax or $expandleft or $expandright or $expanddb) {
	unless ($mce_flag or $evofold_flag or $tfbs_flag or $wgrna_flag or $segdup_flag or $anno_flag or $knowngene_flag or $refgene_flag or $refcds_flag or $refexon_flag or $phastcons_flag) {
		pod2usage ("Error in argument: the --expandmax, --expandleft, --expandright or --expanddb arguments are NOT applicable to plain DB-file format");
	}
}


if ($mce_flag or $evofold_flag or $tfbs_flag or $wgrna_flag or $segdup_flag or $anno_flag) {
	my ($chr_loc, $count_query_bp) = readChrLoc ($queryfile);
	$condense_query and ($chr_loc, $count_query_bp) = condenseOverlap ($chr_loc);
	$count_query_bp or print STDERR "WARNING: there is NOTHING in query to scan\n" and exit(0);	#there is no query regions in the query-location-file
	
	if ($mce_flag) {
		scanUCSCRegion ($chr_loc, $dbfile, $count_query_bp, 'mce');
	} elsif ($evofold_flag) {
		scanUCSCRegion ($chr_loc, $dbfile, $count_query_bp, 'evofold');
	} elsif ($tfbs_flag) {
		scanUCSCRegion ($chr_loc, $dbfile, $count_query_bp, 'tfbs');
	} elsif ($wgrna_flag) {
		scanUCSCRegion ($chr_loc, $dbfile, $count_query_bp, 'wgrna');
	} elsif ($segdup_flag) {
		scanUCSCRegion ($chr_loc, $dbfile, $count_query_bp, 'segdup');
	} elsif ($anno_flag) {
		scanUCSCRegion ($chr_loc, $dbfile, $count_query_bp, 'anno');
	}
} elsif ($knowngene_flag or $refgene_flag or $refcds_flag or $refexon_flag) {
	if ($knowngene_flag) {
		scanUCSCGene ($queryfile, $dbfile, 0, 'knowngene', $kgxref, $name2_flag);
	} elsif ($refgene_flag) {
		scanUCSCGene ($queryfile, $dbfile, 0, 'refgene', $kgxref, $name2_flag);
	} elsif ($refcds_flag) {
		scanUCSCGene ($queryfile, $dbfile, 0, 'refcds', $kgxref, $name2_flag);
	} elsif ($refexon_flag) {
		scanUCSCGene ($queryfile, $dbfile, 0, 'refexon', $kgxref, $name2_flag);
	}
} elsif ($phastcons_flag) {							#special treatment for phsatcons file
	my ($chr_loc, $count_query_bp) = readChrLoc ($queryfile);
	$count_query_bp or print STDERR "WARNING: there is NOTHING in query to scan\n" and exit(0);	#there is no query regions in the query-location-file
	print STDERR "NOTICE: the --phastcons flag requires condensing query regions into a set of non-overlapping regions\n";
	($chr_loc, $count_query_bp) = condenseOverlap ($chr_loc);
	scanPhastCons ($chr_loc, $dbfile, $count_query_bp);
} elsif ($test_flag) {
	scanRegion ($queryfile, $dbfile);
} else {
	scanRegion ($queryfile, $dbfile);
}





sub readChrLoc {
	my ($queryfile) = @_;
	my ($count_query_bp, $count_query_region, %chr_loc);	
	if ($queryfile eq 'stdin') {
		*QUERY = *STDIN;
	} else {
		open (QUERY, $queryfile) or confess "Error: cannot read from query file $queryfile: $!";		#!
	}
	while (<QUERY>) {
		/\S/ or next;												#skip empty lines
		s/[\r\n]+$//;												#get rid of newline characters
		m/^(chr)?(\w+):(\d+)(\-(\d+))?(.*)/ or confess "Error: invalid record in chromosome file: <$_>";	#might contain only one locatoin (SNP) or a start-end (for chromosome regions)
		my ($chr, $loc_start, $loc_end, $other_info) = ($2, $3, $5, $6);
		my $switch = 0;												#switch start and end position
		
		defined $loc_end or $loc_end = $loc_start;
		defined $other_info or $other_info = '';
		if ($loc_start > $loc_end) {										#in-del mutation such as rs34899723 in human
			$switch = 1;
			($loc_start, $loc_end) = ($loc_end, $loc_start);
		}
		push @{$chr_loc{$chr}}, [$loc_start, $loc_end, $other_info, $switch];
	
		$count_query_bp += ($loc_end-$loc_start+1);
		$count_query_region++;
	}
	for my $chr (keys %chr_loc) {
		@{$chr_loc{$chr}} = sort {$a->[0] <=> $b->[0]} @{$chr_loc{$chr}};				#sort the chr_location by the start position for each chr
	}
	$verbose and print STDERR "NOTICE: A total of $count_query_bp (possibly overlapped) base pairs found in $queryfile in ${\(scalar keys %chr_loc)} chromosomes: " , join (" ", sort keys %chr_loc), "\n";
	return (\%chr_loc, $count_query_bp, $count_query_region);
}

sub scanPhastCons {
	my ($chr_loc, $dbfile, $count_query_bp) = @_;
	my ($chr, $start, $current_loc, $step, $qloc, $qpointer);
	my ($count_db_bp, $count_db_bp_found) = (0, 0);
	my (%chr_finish, %chr_pointer);
	
	open (DB, $dbfile) or confess "Error: cannot read from database file $dbfile: $!";
	while (<DB>) {
		if (m/^fixedStep chrom=chr(\w+) start=(\d+) step=(\d+)$/) {
			#$verbose and print STDERR "NOTICE: next header $_";
			($chr, $start, $current_loc, $step) = ($1, $2-1, $2-1, $3);	#the next line is the real start (so first minus 1 to get the real start)
			$qloc = $chr_loc->{$chr} || undef;
			$qpointer = $chr_pointer{$chr} || 0;
			$step == 1 or confess "Error: multiple steps found in phastcons file";
		} else {
			length ($_) == 6 or die "Error: invalid record in db file: <$_>";
			$current_loc++;
			$count_db_bp++;

			$verbose and $count_db_bp =~ m/000000$/ and print STDERR "NOTICE: Processing chr $chr in $count_db_bp base of dbfile $dbfile with $count_db_bp_found/$count_query_bp interested base pairs found in overlapping regions\n";
			$qloc or next;					#the $qloc is the query location array for this chromosome
			$chr_finish{$chr} and next;

			for my $i ($qpointer .. @$qloc-1) {
				my ($qstart, $qend, $qinfo, $qswitch) = @{$qloc->[$i]};
				if ($qend < $current_loc) {
					$chr_pointer{$chr} = $i;
					$i == @$qloc-1 and $chr_finish{$chr} = 1;	#this chromosome is DONE
					next;
				} elsif ($qstart > $current_loc) {
					last;
				} else {
					$count_db_bp_found++;
					print $overlap?"$chr:$current_loc-$current_loc$qinfo":"$chr:$qstart-$qend$qinfo";
					print $append?"\t$_":"\n";
				}
			}
		}
	}
}

sub scanRegion {
	my ($queryfile, $dbfile) = @_;
	my ($db_chr_loc, $count_db_bp, $count_db_region) = readChrLoc ($dbfile);
	($db_chr_loc, $count_db_bp, $count_db_region) = condenseOverlap ($db_chr_loc);			#condense overlapping regions in DB file into consensus regions

	my ($chr_loc, $count_query_bp, $count_query_region) = readChrLoc ($queryfile);
	$count_query_bp or print STDERR "WARNING: there is NOTHING in queryfile $queryfile to scan\n" and exit(0);	#there is no query regions in the query-location-file
	$condense_query and ($chr_loc, $count_query_bp, $count_query_region) = condenseOverlap ($chr_loc);		#--condense_query is rarely necessary to be used

	my ($count_db_bp_found, %chr_finish, %chr_pointer) = (0);
	for my $chr (sort sortChr keys %$db_chr_loc) {
		for my $nextregion (@{$db_chr_loc->{$chr}}) {
			my ($start, $end) = @$nextregion;
			
			my $qloc = $chr_loc->{$chr} or next;				#skip this chromosome because query does not have anything on it
			$chr_finish{$chr} and next;					#skip this chromosome because it has been finished processing by query
	
			my $pointer = $chr_pointer{$chr} || 0;				#current pointer for this chromosome
			for my $i ($pointer .. @$qloc-1) {
				my ($qstart, $qend, $qinfo, $qswitch) = @{$qloc->[$i]};
	
				if (defined $maxquerydbratio and ($qend-$qstart+1)/($end-$start+1) > $maxquerydbratio) {
					next;			#maximum allowable query/db size ratio difference (query is too big to claim overlap/similarity with db)
				}
				if (defined $minquerydbratio and ($qend-$qstart+1)/($end-$start+1) < $minquerydbratio) {
					next;			#minimum allowable query/db size ratio difference (query is too small to claim overlap/similarity with db)
				}
	
				if ($qend < $start) {					#query end is before the start position of this db region
					#db:            <------------------------->
					#query: <--->
					$chr_pointer{$chr} = $i;			#move the pointer only when qend<start (I thought this for a long time)
					$i == @$qloc-1 and $chr_finish{$chr} = 1;	#this chromosome is DONE
					next;
				} elsif ($qend <= $end) {				#move pointer to this region? (probably NOT, especially when one query overlap with another)
					if ($qstart >= $start) {			#query contained completely within db region
						#db:      <-------------------------->
						#query:       <------------------>
						$count_db_bp_found += ($qend-$qstart+1);
						if (defined $minoverlap) {
							1;				#query is contained within hit
						}
						if (defined $maxoverlap) {
							($qend-$qstart+1)/($end-$start+1) > $maxoverlap and next;
						}
						
						if (defined $mindbfrac) {
							($qend-$qstart+1)/($end-$start+1) < $mindbfrac and next;
						}
						if (defined $maxdbfrac) {
							($qend-$qstart+1)/($end-$start+1) > $maxdbfrac and next;
						}
						if (defined $minqueryfrac) {
							1;
						}
						if (defined $maxqueryfrac) {
							1;
						}

						if ($dbregion) {
							print "chr$chr:$start-$end", $queryinfo?$qinfo:"";	#there is no qinfo for db region
						} else {
							print "chr$chr:$qstart-$qend$qinfo";
						}
						print "\n";
					} else {					#query overlap but upstream of db region
						#db:       <------------------------->
						#query: <---------------------->
						$count_db_bp_found += ($qend-$start+1);
						if (defined $minoverlap) {
							if (($qend-$start+1)/($qend-$qstart+1) < $minoverlap and ($qend-$start+1)/($end-$start+1) < $minoverlap) {
								next;
							}
						}
						if (defined $maxoverlap) {
							if (($qend-$start+1)/($qend-$qstart+1) > $maxoverlap and ($qend-$start+1)/($end-$start+1) > $maxoverlap) {
								next;
							}
						}
						
						if (defined $mindbfrac) {
							($qend-$start+1)/($end-$start+1) < $mindbfrac and next;
						}
						if (defined $maxdbfrac) {
							($qend-$start+1)/($end-$start+1) > $maxdbfrac and next;
						}
						if (defined $minqueryfrac) {
							($qend-$start+1)/($qend-$qstart+1) < $minqueryfrac and next; 
							#print "$queryinfo ($qend-$start+1)/($qend-$qstart+1) < $minqueryfrac\n";
						}
						if (defined $maxqueryfrac) {
							($qend-$start+1)/($qend-$qstart+1) > $maxqueryfrac and next;
						}
						
						if ($dbregion) {
							print $overlap?"chr$chr:$start-$qend$qinfo":"chr$chr:$start-$end", $queryinfo?$qinfo:"";
						} else {
							print $overlap?"chr$chr:$start-$qend$qinfo":"chr$chr:$qstart-$qend$qinfo";
						}
						print "\n";
					}
				} elsif ($qstart <= $end) {
					if ($qstart >= $start) {			#query overlap but downstream of db region
						#db:      <------------------------>
						#query:        <----------------------->
						$count_db_bp_found += ($end-$qstart+1);
						if (defined $minoverlap) {
							if (($end-$qstart+1)/($qend-$qstart+1) < $minoverlap and ($end-$qstart+1)/($end-$start+1) < $minoverlap) {
								next;
							}
						}
						if (defined $maxoverlap) {
							if (($end-$qstart+1)/($qend-$qstart+1) > $maxoverlap and ($end-$qstart+1)/($end-$start+1) > $maxoverlap) {
								next;
							}
						}
						
						if (defined $mindbfrac) {
							($end-$qstart+1)/($end-$start+1) < $mindbfrac and next;
						}
						if (defined $maxdbfrac) {
							($end-$qstart+1)/($end-$start+1) > $maxdbfrac and next;
						}
						if (defined $minqueryfrac) {
							($end-$qstart+1)/($qend-$qstart+1) < $minqueryfrac and next;
						}
						if (defined $maxqueryfrac) {
							($end-$qstart+1)/($qend-$qstart+1) > $maxqueryfrac and next;
						}
						#print "$queryinfo ($end-$start+1)/($qend-$qstart+1)\n";
						if ($dbregion ) {
							print $overlap?"chr$chr:$qstart-$end$qinfo":"chr$chr:$start-$end", $queryinfo?$qinfo:"";
						} else {
							print $overlap?"chr$chr:$qstart-$end$qinfo":"chr$chr:$qstart-$qend$qinfo";
						}
						print "\n";
					} else {					#db region completely contained within query
						#db:      <------------------------->
						#query: <------------------------------>
						$count_db_bp_found += ($end-$start+1);
						if (defined $minoverlap) {
							1;
						}
						if (defined $maxoverlap) {
							($end-$start+1)/($qend-$qstart+1) > $maxoverlap and next;
						}
						
						if (defined $mindbfrac) {
							1;
						}
						if (defined $maxdbfrac) {
							1;
						}
						if (defined $minqueryfrac) {
							($end-$start+1)/($qend-$qstart+1) < $minqueryfrac and next; 
							#print "$queryinfo ($end-$start+1)/($qend-$qstart+1)\n";
						}
						if (defined $maxqueryfrac) {
							($end-$start+1)/($qend-$qstart+1) > $maxqueryfrac and next;
							#print "$queryinfo ($end-$start+1)/($qend-$qstart+1)\n";
						}
						
						if ($dbregion) {
							print $overlap?"chr$chr:$start-$end$qinfo":"chr$chr:$start-$end", $queryinfo?$qinfo:"";
						} else {
							print $overlap?"chr$chr:$start-$end$qinfo":"chr$chr:$qstart-$qend$qinfo";
						}
						print "\n";
					}
				} else {
					#db:      <---------->
					#query:       	         <----------------------->
					last;						#should examine next db region
				}
			}
		}
		$verbose and $count_db_region =~ m/00000$/ and print STDERR "NOTICE: Processing chr $chr in $count_db_region line of template file $dbfile with $count_db_bp_found/$count_query_bp interested base pairs found in overlapping regions\n";
	}
	$verbose and print STDERR "NOTICE: Total of $count_db_region template chr_regions (total of $count_db_bp) examined and found $count_db_bp_found / $count_query_bp (${\($count_db_bp_found / $count_query_bp)}) interested base pairs that belong to these template regions\n";

}

#this subroutine is used to scan various files downloaded from UCSC table browser (usually tab-delimited files). Many of them have the same or similar file formats
sub scanUCSCRegion {
	my ($chr_loc, $dbfile, $count_query_bp, $dbtype) = @_;
	my (%chr_warning, %chr_finish, %chr_pointer);			#warning for no chr information in query file, chr in query has been finished processing, pointer to current query chr info
	my ($count_db_bp, $count_db_region, $count_db_bp_found, $count_db_region_found) = (0, 0, 0, 0);
	my @record;
	my ($chr, $start, $end, $score, $normscore);			#score is defined as LOD for MCE, Z for TFBS and conf for Evofold; normscore is between 0-1000
	
	
	open (DB, $dbfile) or confess "Error: cannot read from database file $dbfile: $!";
	while (<DB>) {
		s/[\r\n]+$//;							#deleting the newline characters
		@record = split (/\t/, $_);
		$record[0] eq '#bin' and next;					#comment line
		if ($dbtype eq 'mce') {
			@record == 6 or confess "Error: invalid record in dbfile $dbfile (expecting 10 fields in evofold file): <$_>";
			($chr, $start, $end, $score, $normscore) = ($record[1], $record[2], $record[3], $record[4], $record[5]);
			$score =~ s/^lod=// or confess "Error: invalid lod score designation (no 'lod=' found) in dbfile $dbfile: <$_>";
		} elsif ($dbtype eq 'evofold') {
			@record == 10 or confess "Error: invalid record in dbfile $dbfile (expecting 10 fields in evofold file): <$_>";
			($chr, $start, $end, $score, $normscore) = ($record[1], $record[2], $record[3], $record[5], $record[5]);
		} elsif ($dbtype eq 'tfbs') {
			@record == 8 or confess "Error: invalid record in dbfile $dbfile (expecting 10 fields in evofold file): <$_>";
			($chr, $start, $end, $score, $normscore) = ($record[1], $record[2], $record[3], $record[7], $record[5]);
		} elsif ($dbtype eq 'wgrna') {
			@record == 10 or confess "Error: invalid record in dbfile $dbfile (expecting 10 fields in evofold file): <$_>";
			($chr, $start, $end, $score, $normscore) = ($record[1], $record[2], $record[3], $record[5], $record[5]);
		} elsif ($dbtype eq 'segdup') {
			@record == 30 or confess "Error; invalid record in dbfile $dbfile (expecting 30 fields in segdup file): <$_>";
			($chr, $start, $end, $score, $normscore) = @record[1, 2, 3, 5, 5];
		} elsif ($dbtype eq 'anno') {
			($chr, $start, $end, $score, $normscore) = @record[1, 2, 3, 5, 5];
		} else {
			confess "Error: the dbtype $dbtype cannot be handled by the current program";
		}

		if ($dbtype =~ m/^(mce|evofold|tfbs|wgrna|segdup|anno)$/) {
			$score_threshold and $score < $score_threshold and next;			#if --score_threshold is set, the low scoring segment will be skipped
			$normscore_threshold and $normscore < $normscore_threshold and next;		#if --normscore_threshold is set, the low scoring segment will be skipped
			$start++;									#due to the zero-opening coordinate system in UCSC
		}

		$chr =~ s/^chr// or confess "Error: invalid chromosome designation (no 'chr' found) in dbfile $dbfile: <$_>";
		$count_db_region++;
		$count_db_bp += ($end-$start+1);
		
		my $qloc = $chr_loc->{$chr} or next;				#skip this chromosome because query does not have anything on it
		$chr_finish{$chr} and next;					#skip this chromosome because it has been finished processing by query

		my $pointer = $chr_pointer{$chr} || 0;				#current pointer for this chromosome
		for my $i ($pointer .. @$qloc-1) {
			my ($qstart, $qend, $qinfo, $qswitch) = @{$qloc->[$i]};

			if ($qend < $start) {					#query end is before the start position of this db region
				#db:            <------------------------->
				#query: <--->
				$chr_pointer{$chr} = $i;			#move the pointer only when qend<start (I thought this for a long time)
				$i == @$qloc-1 and $chr_finish{$chr} = 1;	#this chromosome is DONE
				next;
			} elsif ($qend <= $end) {				#move pointer to this region? (probably NOT, especially when one query overlap with another)
				if ($qstart >= $start) {			#query contained completely within db region
					#db:      <-------------------------->
					#query:       <------------------>
					$count_db_bp_found += ($qend-$qstart+1);
					if ($dbregion) {
						print "chr$chr:$start-$end";	#there is no qinfo for db region
					} else {
						print "chr$chr:$qstart-$qend$qinfo";
					}
					$append and print "\t$score\t$normscore";
					print "\n";
				} else {					#query overlap but upstream of db region
					#db:       <------------------------->
					#query: <---------------------->
					$count_db_bp_found += ($qend-$start+1);
					if ($minoverlap) {
						if (($qend-$start+1)/($qend-$qstart+1) < $minoverlap and ($qend-$start+1)/($end-$start+1) < $minoverlap) {
							next;
						}
					}
					if ($dbregion) {
						print $overlap?"chr$chr:$start-$qend$qinfo":"chr$chr:$start-$end";
					} else {
						print $overlap?"chr$chr:$start-$qend$qinfo":"chr$chr:$qstart-$qend$qinfo";
					}
					$append and print "\t$score\t$normscore";
					print "\n";
				}
			} elsif ($qstart <= $end) {
				if ($qstart >= $start) {			#query overlap but downstream of db region
					#db:      <------------------------>
					#query:        <----------------------->
					$count_db_bp_found += ($end-$qstart+1);
					if ($minoverlap) {
						if (($end-$qstart+1)/($qend-$qstart+1) < $minoverlap and ($end-$qstart+1)/($end-$start+1) < $minoverlap) {
							next;
						}
					}
					if ($dbregion ) {
						print $overlap?"chr$chr:$qstart-$end$qinfo":"chr$chr:$start-$end";
					} else {
						print $overlap?"chr$chr:$qstart-$end$qinfo":"chr$chr:$qstart-$qend$qinfo";
					}
					$append and print "\t$score\t$normscore";
					print "\n";
				} else {					#db region completely contained within query
					#db:      <------------------------->
					#query: <------------------------------>
					$count_db_bp_found += ($end-$start+1);
					if ($dbregion) {
						print $overlap?"chr$chr:$start-$end$qinfo":"chr$chr:$start-$end";
					} else {
						print $overlap?"chr$chr:$start-$end$qinfo":"chr$chr:$qstart-$qend$qinfo";
					}
					$append and print "\t$score\t$normscore";
					print "\n";
				}
			} else {
				last;						#should examine next db region
			}
		}
		$verbose and $count_db_region =~ m/00000$/ and print STDERR "NOTICE: Processing chr $chr in $count_db_region line of template file $dbfile with $count_db_bp_found/$count_query_bp interested base pairs found in overlapping regions\n";
	}
	$quiet or print STDERR "NOTICE: Total of $count_db_region template chr_regions (total of $count_db_bp) examined and found $count_db_bp_found / $count_query_bp (${\($count_db_bp_found / $count_query_bp)}) interested base pairs that belong to these template regions\n";

}

#this subroutine scan a queryfile against a knownGene file or refGene file from UCSC to identify the nearby genes for each query_region in queryfile
#instead of scanning each line in the template file, I decided to just read everything in memory, then process each query sequentially. this is conceptually different and opposite the scanDBRegion subroutine
sub scanUCSCGene {
	my ($queryfile, $dbfile, $count_query_bp, $dbtype, $kgxref, $name2_flag) = @_;
	my (%chr_warning, %chr_finish, %chr_pointer);				#warning for no chr information in query file, chr in query has been finished processing, pointer to current query chr info
	my ($count_db_bp, $count_db_region) = (0, 0);
	my (%db_chr_loc);
	my ($name, $chr, $dbstrand, $start, $end, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend);
	my ($gene_xref);

	open (DB, $dbfile) or confess "Error: cannot read from database file $dbfile: $!";
	while (<DB>) {
		m/^#/ and next;							#comment line
		s/[\r\n]+$//;							#deleting the newline characters
		my @record = split (/\t/, $_);

		if ($dbtype eq 'knowngene') {
			@record == 12 or confess "Error: invalid record in dbfile $dbfile (expecting 12 fields in knownGene file): <$_>";
			($name, $chr, $dbstrand, $start, $end, $cdsstart, $cdsend) = @record[0..6];				#human hg17
		} elsif ($dbtype eq 'refgene' or $dbtype eq 'refcds' or $dbtype eq 'refexon') {
			if (@record == 16) {
				($name, $chr, $dbstrand, $start, $end, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend) = @record[1..10];		#human hg18, mouse
				$name2_flag and $name = $record[12];
			} elsif (@record == 10) {
				($name, $chr, $dbstrand, $start, $end, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend) = @record;			#rat, human hg17
				$name2_flag and pod2usage ("Error in argument: the --name2_flag argument can be used only for refGene table with name2 annotations (such as those in hg18 and mm8 database)");
			} else {
				confess "Error: invalid record in template-location-file $dbfile (expecting 16 or 10 tab-delimited fields in refGene file): <$_>";
			}
		} else {
			confess "Error: unknown dbtype flat: $dbtype\n";
		}
		
		if (defined $strand) {						#only consider a particular strand (by default both strand will be considered)
			$strand eq 'plus' || $strand eq 'forward' and $dbstrand eq '+' || next;
			$strand eq 'minus' || $strand eq 'reverse' and $dbstrand eq '-' || next;
		}

		if ($expanddb) {			#expand the definition of gene
			$start -= $expanddb;
			$end += $expanddb;
		}

		if ($dbtype eq 'knowngene' or $dbtype eq 'refgene') {
			$start++;							#due to the zero-opening coordinate system in UCSC
			$chr =~ s/^chr// or confess "Error: invalid chromosome designation (no 'chr' found) in dbfile $dbfile: <$_>";
			$count_db_region++;
			$count_db_bp += ($end-$start+1);
			push @{$db_chr_loc{$chr}}, [$start, $end, $name];
		} elsif ($dbtype eq 'refcds') {
			$cdsstart++;							#due to the zero-opening coordinate system in UCSC
			$chr =~ s/^chr// or confess "Error: invalid chromosome designation (no 'chr' found) in dbfile $dbfile: <$_>\n";
			$count_db_region++;
			$count_db_bp += ($cdsend-$cdsstart+1);
			push @{$db_chr_loc{$chr}}, [$cdsstart, $cdsend, $name];
		} elsif ($dbtype eq 'refexon') {
			$chr =~ s/^chr// or confess "Error: invalid chromosome designation (no 'chr' found) in dbfile $dbfile: <$_>\n";
			my @exonstart = split (/,/, $exonstart);
			my @exonend = split (/,/, $exonend);
			@exonstart == @exonend or confess "Error: invalid exon annotation in dbfile $dbfile (exonstart=@exonstart exonend=@exonend): <$_>\n";
			map {$_++} @exonstart;
			for my $i (0 .. @exonstart-1) {
				$count_db_region++;
				$count_db_bp += ($exonend[$i]-$exonstart[$i]+1);
				push @{$db_chr_loc{$chr}}, [$exonstart[$i], $exonend[$i], $name];
			}
		}
	}
	for my $chr (sort keys %db_chr_loc) {					#sort dbregion to make sure that smaller regions occur before bigger ones
		@{$db_chr_loc{$chr}} = sort {$a->[0] <=> $b->[0] or $a->[1] <=> $b->[1]} @{$db_chr_loc{$chr}};
	}
	$verbose and print STDERR "NOTICE: Total of $count_db_region dbregions are found in $dbfile with $count_db_bp base pairs\n";
	
	#now scan each query_loc against all sorted db entry in the template file
	#by default, any gene/transcript/cds/exon overlapping with query will be printed out
	#if no gene/transcript/cds/exon overlap, then the closest gene (as specified by --expandmax) will be printed out
	#db: ---------          -------------------                  -----        ---------------
	#        --------    -------------------------    ---          -------------                ----------
	#query:    <------->                     <-------------------------->
	#                       <--------->                                             <-->
	#                                                       <-->

	if ($kgxref) {		#in the output, rather than using gene identifier (name), use the corresponding gene symbol in the known gene cross-reference file
		$gene_xref = readKgXref ($kgxref);
	} elsif ($reflink) {
		$gene_xref = readRefLink ($reflink);
	}
		
	if ($queryfile eq 'stdin') {
		*QUERY = *STDIN;
	} else {
		open (QUERY, $queryfile) or confess "Error: cannot read from queryfile $queryfile: $!";
	}
	while (<QUERY>) {
		s/[\r\n]+$//;												#get rid of newline characters
		/\S/ or next;												#skip empty lines
		m/^(chr)?(\w+):(\d+)(\-(\d+))?(.*)/ or confess "Error: invalid record in chromosome file: <$_>";	#might contain only one locatoin (SNP) or a start-end (for chromosome regions)
		my ($chr, $qstart, $qend, $qinfo, $qswitch) = ($2, $3, $5||$3, $6||'');
		if ($qstart > $qend) {
			$qswitch = 1;
			($qstart, $qend) = ($qend, $qstart);
		}
		
		my (@closest_gene, $closest_dist, @overlap);
		my ($qstart_expand, $qend_expand) = ($qstart - $expandmax, $qend + $expandmax);
		$expandleft and $qstart_expand = $qstart - $expandleft;
		$expandright and $qend_expand = $qend + $expandright;
		
		my $dbloc = $db_chr_loc{$chr};
		if (not $dbloc) {								#for example, when chr is 22_random, there is NO genes there
			print "$chr:$qstart-$qend$qinfo\tNOT_FOUND\tNOT_FOUND\n";
			next;
		}
		for my $j (0 .. @$dbloc-1) {
			my ($dbstart, $dbend, $dbname) =@{$dbloc->[$j]};
			$qstart_expand > $dbend and next;					#the dbend has not reached expanded query yet
			if ($qend >= $dbstart and $qstart <= $dbend) {				#overlap found between query and db!
				push @overlap, $dbname;
			}
			if (not @overlap and $expandleft || $expandright || $expandmax) {
				if ($qend_expand >= $dbstart and $qstart_expand <= $dbend) {	#overlap found between expanded query and db
					my $dist;
					if ($qend <= $dbstart) {
						$dist = $dbstart-$qend;
					} else {
						$dist = $qstart-$dbend;
					}

					if ($closest_dist and $dist < $closest_dist) {
						$closest_dist = $dist;
						@closest_gene = $dbname;
					} elsif (not defined $closest_dist or $dist == $closest_dist) {
						$closest_dist = $dist;
						push @closest_gene, $dbname;

					}
				}
			}
			$qend_expand < $dbstart and last;	#the dbstart has passed expanded query
		}
		
		if ($kgxref or $reflink) {
			@overlap = map {$gene_xref->{$_} || $_} @overlap;
			@closest_gene = map {$gene_xref->{$_} || $_} @closest_gene;
		}
			
		
		my %seen;
		@overlap = sort grep {!$seen{$_}++} @overlap;
		@closest_gene = sort grep {!$seen{$_}++} @closest_gene;
		
		if (@overlap) {
			print "chr$chr:$qstart-$qend$qinfo\t", join(",", @overlap), "\t0\n";
		} elsif (@closest_gene) {
			print "chr$chr:$qstart-$qend$qinfo\t", join(",", @closest_gene), "\t$closest_dist\n";
		} else {
			print "chr$chr:$qstart-$qend$qinfo\tNOT_FOUND\tNOT_FOUND\n";
		}
	}
}

#sometimes two exons enclose with each other, sometimes two regions locate in two strands so have overlaps, so it is important to eliminate overlap positions
#this subroutine reduce overlaps, when feeded with a series of chromosome regions in two arrays, one for start position and one for end position
sub condenseOverlap {
	my ($chr_loc) = @_;
	my ($length, $newlength, $newcount) = (0, 0, 0);
	for my $chr (keys %$chr_loc) {
		my @loc = @{$chr_loc->{$chr}};
		my @newloc;
		my $pointer = 0;
		my ($prestart, $preend) = ($loc[0]->[0], $loc[0]->[1]);
		$length += $loc[0]->[1]-$loc[0]->[0]+1;
		
		while (1) {
			$pointer++;					#process next segment
			$pointer == @loc and last;			#finish processing the loc array
			$length += ($loc[$pointer]->[1]-$loc[$pointer]->[0]+1);
			if ($loc[$pointer]->[0] <= $preend) {		#start of next element less than end of previous element
				if ($loc[$pointer]->[1] >= $preend) {	#make sure the next element is not contained within the previous element
					$preend = $loc[$pointer]->[1];
				}
			} else {
				push @newloc, [$prestart, $preend, '', ''];
				$newlength += ($preend-$prestart+1);
				($prestart, $preend) = ($loc[$pointer]->[0], $loc[$pointer]->[1]);
			}
		}
		push @newloc, [$prestart, $preend, '', ''];		#process the last element
		$newlength += ($preend-$prestart+1);
		$newcount += @newloc;

		$chr_loc->{$chr} = \@newloc;
	}
	$verbose and print STDERR "NOTICE: After cleaning query, the length of query changes from $length to $newlength\n";
	return ($chr_loc, $newlength, $newcount);
}

sub sortChr {
	my ($a, $b) = @_;		#this line is required for functioning properly and I do not know why
	if ($a =~ m/^(\d+)$/) {
		if ($b =~ m/^(\d+)$/) {
			return $a<=>$b;
		} else {
			return -1;
		}
	} elsif ($b =~ m/^(\d+)$/) {
		return 1;
	} else {
		return $a cmp $b;
	}
}

sub readKgXref {
	my ($inputfile) = @_;
	my (%gene_xref);
	open (XREF, $inputfile) or confess "Error: cannot read from kgxref file $inputfile: $!";
	while (<XREF>) {
		m/^#/ and next;
		my @record = split (/\t/, $_);
		if ($gene_xref{$record[0]}) {			#BC003168 occur twice in kgxref file (OSBPL10, BC003168)
			if ($gene_xref{$record[0]} =~ m/^(BC|AK)\d+$/) {
				$gene_xref{$record[0]} = $record[4];
			}
		} else {
			$gene_xref{$record[0]} = $record[4];
		}
	}
	close (XREF);
	return (\%gene_xref);
}

sub readRefLink {
	my ($inputfile) = @_;
	my (%gene_xref);
	open (REFLINK, $inputfile) or confess "Error: cannot read from refLink file $inputfile: $!";
	while (<REFLINK>) {
		m/^#/ and next;
		my @record = split (/\t/, $_);
		$gene_xref{$record[2]} = $record[0];
	}
	close (REFLINK);
	return (\%gene_xref);
}


=head1 SYNOPSIS

 scan_region.pl [arguments] <query-file> <DB-file>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
        
        Flags specifying type of db-file
            --mce_flag			dbfile is UCSC MCE annotation file
            --evofold_flag		dbfile is UCSC EvoFold annotation file
            --tfbs_flag			dbfile is UCSC TFBS annotation file
            --wgrna_flag		dbfile is UCSC wgRna annotation file
            --knowngene_flag		dbfile is UCSC knownGene annotation file
            --refgene_flag		dbfile is UCSC refGene annotation file
            --phastcons_flag		dbfile is UCSC phastcons conservation score file
            --anno_flag			dbfile is in generic UCSC annotation database format
            --refexon_flag		dbfile is UCSC refGene with exon annotation
            --refcds_flag		dbfile is UCSC refGene with coding sequence annotation
        
        Criteria for defining query-db match
            --condense_query		condense and eliminate overlapping regions in query
            --score_threshold <float>	score threshold for db in UCSC annotation file
            --normscore_threshold <float>	normalized score threshold for db in UCSC annotation file
            --minoverlap <float>	minimum portion of overlap (either query or db) to decide query-db match
            --maxoverlap <float>	maximum portion of overlap (either query or db) to decide query-db match
            --minquerydbratio <float>	minimum query to db length ratio to decide match
            --maxquerydbratio <float>	maximum query to db length ratio to decide match
            --mindbfrac <float>		minimum fraction of db in overlap
            --maxdbfrac <float>		maximum fraction of db in overlap
            --minqueryfrac <float>	minimum fraction of query in overlap
            --maxqueryfrac <float>	maximum fraction of query in overlap
            
         Expansion of query to find match:
            --expandleft <int>		expand left side of query regions (overwrite --expandmax)
            --expandright <int>		expand right side of query regions (overwrite --expandmax)
            --expandmax <int>		size of maximum expansion for query region to find overlap
         
         Expansion of db to find match:
            --expanddb <int>		expand definition of gene/cds/exon at both sides

         dbfile-specific arguments:
            --kgxref <file>		use gene symbol in known gene xref file in output
            --reflink <file>		use gene symbol in refGene link file in output
            --name2_flag		use name2 annotation in refGene file in output
         
         Input/output options:
            --quiet			suppress printing progress messages
            --append			append extra information from annotation file to output
            --dbregion			print db region (default is to print query region)
            --queryinfo			force to print query info when --dbregion is used
            --overlap			print overlapped portion of region only

 Function: scan genomic regions in a query-file against a DB-file which contains 
 chromosome locations for various genomics features
 
 Example: scan_region.pl query.txt refGene.txt -refgene -reflink refLink.txt
          scan_region.pl query.txt db.txt
          scan_region.pl cnvcall.txt centromere.txt -minqueryfrac 0.5

 Version: $LastChangedDate: 2011-01-17 01:16:59 -0800 (Mon, 17 Jan 2011) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--overlap>

instead of printing the query region, only print the overlapped portion of the 
query region and template region.

=item B<--dbregion>

print the region in DB file, rather than query file, when an overlapped 
hit is found

=item B<--append>

append the score and normscore for the overlapped template region to the output 
for DB files downloaded as UCSC tables.

=item B<--condense_query>

condense overlapped regions in the query file into non-overlapped regions. When 
this argument is set, the annotation for each query (the strings after the 
chromosome location in each line of the query file) will not in the output.

=item B<--minoverlap>

the minimum percentage of overlap between query and DB chromosome region to 
infer a matching event. By default, even a single base pair overlap will be 
considered as overlap, but a lot of times people prefer to use something like 
0.5 (50% overlap) to make sure that the query and DB regions have high 
concordance. This argument operates by OR operation: if overlap for either query 
or DB exceed the specified criteria, it will be treated as passing the criteria.

=item B<--maxoverlap>

the maximum percentage of overlap between query and DB chromosome region to 
infer a matching event. 

=item B<--minquerydbratio>

minimum query to db length ratio to infer a matching event

=item B<--maxquerydbratio>

maximum query to db length ratio to infer a matching event

=item B<--mindbfrac>

minimum fraction of db in overlap to infer a matching event

=item B<--maxdbfrac>

maximum fraction of db in overlap to infer a matching event

=item B<--minqueryfrac>

minimum fraction of query in overlap to infer a matching event

=item B<--maxdbfrac>

maximum fraction of query in overlap to infer a matching event

=item B<--score_threshold>

specify the score threshold in the DB file to include in the search for 
overlaps. This argument is file format dependent.

=item B<--normscore_threshold>

specify the normalized score threshold in the DB file to include in the 
search for overlaps. This argument is file format dependent.

=item B<--mce_flag>

specify that the DB file is in MCE format from UCSC genome browser

=item B<--evofold_flag>

specify that the DB file is in EvoFold format from UCSC genome browser

=item B<--tfbs_flag>

specify that the DB file is in TFBS format from UCSC genome browser

=item B<--wgrna_flag>

specify that the DB file is in WGRNA format from UCSC genome browser

=item B<--anno_flag>

generic UCSC genome browser annotation database format (first a few columns are 
bin, chrom, chromStart, chromEnd, name)

=item B<--knowngene_flag>

specify that the DB file is in knownGene format from UCSC genome browser

=item B<--refgene_flag>

specify that the DB file is in refGene format from UCSC genome browser

=item B<--refcds_flag>

specify that the DB file is in refGene format from UCSC genome browser, but user 
is only interested in the overlap of coding region (first exon to last exon).

=item B<--refexon_flag>

specify that the DB file is in refGene format from UCSC genome browser, but user 
is only interested in the overlap of query with exons.

=item B<--phastcons_flag>

specify that the DB file is in phastcons format from UCSC genome browser

=item B<--expandleft>

expand the query region on the left side (5' in forward strand, 3' in reverse 
strand) to find overlap (used in conjunction with --refgene or --knowngene 
argument)

=item B<--expandright>

expand the query region on the right side (3' in forward strand, 5' in reverse 
strand) to find overlap (used in conjunction with --refgene or --knowngene 
argument)

=item B<--expandmax>

maximum expansion size of the query region on both side to find at least one 
overlap (used in junction with --refgene or --knowngene argument). After query 
expansion, only the closet gene will be printed; other genes, even if 
overlapping with the query after expansion, will not be printed.

=item B<--expanddb>

expand the chromosome region specified in the DB-file to find overlap with the 
query regions.

=item B<--reflink>

specify a cross-reference file for the RefGene track in UCSC genome browser, 
so that in the output, the gene identifier (gene name or refseq id) are replaced 
by the gene symbol specified in the link file. (If not found in the reflink 
file, the gene identifiers are still used)

=item B<--kgxref>

specify a cross-reference file for the knownGene track in UCSC genome browser, 
so that in the output, the gene identifier (gene name or refseq id) are replaced 
by the gene symbol specified in the kgxref file. (If not found in the kgxref 
file, the gene identifiers are still used)

=item B<--name2_flag>

this argument is used in conjunction with the --refgene argument, to specify 
that the alternative gene symbol in the "name2" field in the refGene file be 
printed in the output.

=back

=head1 DESCRIPTION

This program is used to scan chromosome regions from a query-location-file 
against a db-location-file that has been sorted by chromosome and start 
location. Both the query-location-file and the db-location-file contains one 
chromosome location per line. GENERALLY SPEAKING, THE DB FILE MUST BE SORTED BY 
CHROMOSOME AND START POSITIONS! Some exceptions exist, for example, when db-file 
is in UCSC RefGene format.

The various file formats for query file and db file are described below:

=over 8

=item * B<query file>

A sample query-location-file is shown below (The I<chr> prefix is 
optional but highly desirable in the beginning of each line):

	chr3:2000-349990
	chr19:32333-52333

If the query-location-file is for SNPs, then it is okay to use only the start 
location in the line (for example, chr3:2000).

=item * B<DB-file: plain format>

DB-file can be in multiple different formats. The plain format is the simplest 
format, in the same form of "chr3:2000-349990" (same as query-file).

=item * B<DB-file: General overview of UCSC genome browser tables>

When processing DB files as UCSC tables, this program works by first reading all 
information from the query file and store them in memory, then scan each line of 
the DB file to find overlaps. THEREFORE, THE DB FILE CANNOT CONTAIN TWO 
LOCATIONS THAT SHARE OVERLAPS. IN ADDITION, THE TEMPLATE FILE MUST BE SORTED 
FIRST BY CHROMOSOME THEN BY START LOCATION. The chromosome can be sorted 
alphabetically (since there will usually be X, Y and MT chromosomes), but start 
location should be always sorted numerically.

These DB files can be downloaded from the UCSC Genome Browser Annotation 
Database. Generally speaking, I recommended saving the DB file with the same 
file name as the original database file, but prefixing with the version and 
organism of the corresponding genome. For example, for human genome, the 
DB-file name should start with hg17_ or hg18_, to eliminate confusions 
when dealing with different genome builds. For example, after downloading the 
file evofold.txt.gz for human May 2004 assembly from the UCSC genome browser 
database, you can use the sort command to rename the file as hg17_evofold.sorted 
(see details below).

=item * B<MCE format>

The --mce_flag argument specifies that template file is in MCE (most conserved element) format from the 
UCSC genome browser. There are two ways to download the file: You can download a 
MCE template file using a simple command:

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/phastConsElements17way.txt.gz

After downloading the file, next decompress it:

	gunzip phastConsElements17way.txt.gz

This will generate the phastConsElements17way.txt file. Next perform sorting:

	sort -k 2,2 -k 3,3n phastConsElements17way.txt > hg18_phastConsElements17way.sorted

A sample sorted MCE file is shown below:

	585     chr1    1865    1943    lod=31  307
	585     chr1    2039    2096    lod=104 444
	585     chr1    2473    2566    lod=179 505
	585     chr1    2873    2917    lod=107 447
	585     chr1    3081    3135    lod=58  378

The --score_threshold argument operates on the fifth field in the line (such as 
"lod=31"), and the --normscore_threshold argumet operates on the sixth field 
(such as 307).

To scan the queryfile against a MCE file:

	scan_region.pl queryfile hg18_phastConsElements17way.sorted -mce

Note that in the above command, we specified the -mce argument, indicating that 
the DB-file is in MCE format. This argument can be given as single dash (-mce) 
or double dash (--mce), as partial name (-mce) or full name (-mce_flag) as long 
as there is no ambiguity with another argument.

Using the above query file with two lines, the output will be multiple 
identical lines, each representing one hit in the MCE file. The first 
a few lines are:

	19:32333-52333
	19:32333-52333
	19:32333-52333
	19:32333-52333
	19:32333-52333

This means that the query has multiple hits in the DB-file, and for each match, 
the query is printed out.

If only the region in MCE DB-file that overlap with query is desired, you can use:

	scan_region.pl queryfile hg18_phastConsElements17way.sorted -mce -dbregion

The first a few lines in output file are:

	19:35715-35751
	19:37371-37466
	19:38216-38406
	19:38437-38504
	19:38551-38678

If only the overlapped region between query-file and DB-file is desired, you can use:

	scan_region.pl queryfile hg18_phastConsElements17way.sorted -mce -overlap

The output will be the overlapped chromosome regions between query and MCE file. 
The first a few output lines are:

	19:35715-35751
	19:37371-37466
	19:38216-38406
	19:38437-38504
	19:38551-38678

In this case, the query region is very big, so that the overlapped regions are 
usually identical to the dbregion (region in the MCE file).

=item * B<EvoFold format>

The EvoFold file from UCSC genome browser contains genomic segments that are 
predicted to retain stable RNA structural folds in different species. You can 
download a EvoFold DB-file for human NCBI36 genome assembly using

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/evofold.txt.gz
	gunzip evofold.txt.gz
	sort -k 2,2 -k 3,3n evofold.txt > hg18_evofold.sorted

The first a few lines of a sorted EvoFold table is shown below:

	591     chr1    886773  886798  608_0_-_96      96      -       25      (((((.....(((...))).)))))       0.97,0.98,0.99,0.99,0.9,0.78,0.87,0.99,0.98,0.98,0.21,0.22,0.11,0.99,0.95,0.91,0.11,0.22,0.21,0.98,0.9,0.99,0.99,0.98,0.97
	591     chr1    888417  888435  617_0_+_156     156     +       18      ((((((......))))))      0.97,0.99,0.99,1.0,0.99,0.92,1.0,1.0,0.99,1.0,1.0,1.0,0.92,0.99,1.0,0.99,0.99,0.97

The --score_threshold and --normscore_threshold arguments both operate on the 
sixth field in the line (such as "96").

=item * B<TFBS format>

The TFBS file from UCSC genome browser contains predicted transcription factor 
binding sites (TFBS), based on positional weight matrices (PWM) from the TRANFAC 
database. You can download a TFBS DB-file for human NCBI36 genome assembly 
using:

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/tfbsConsSites.txt.gz
	gunzip tfbsConsSites.txt.gz
	sort -k 2,2 -k 3,3n tfbsConsSites.txt > hg18_tfbsConsSites.sorted

The first a few lines of a sorted TFBS table is shown below:

	585     chr1    1832    1848    V$ARP1_01       818     -       2.01
	585     chr1    3211    3232    V$NRSF_01       749     -       2.19
	585     chr1    3603    3623    V$YY1_02        790     -       1.93
	585     chr1    3949    3956    V$NKX25_01      1000    -       2.48
	585     chr1    3996    4014    V$CART1_01      798     -       1.75

The --score_threshold argument operates on the eighth field in the line (such as 
"2.01"), and the --normscore_threshold argumet operates on the sixth field 
(such as 818).

To decode the strange values like V$ARP1_01 and V$NRSF_01, you can additionally 
download the TF annotation file at 
http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/tfbsConsFactors.txt.gz 
and examine the file for detailed information on the transcription factors 
binding to the sites.

=item * B<wgRna format>

The wgRna table from UCSC genome browser contains microRNA and small nucleolar 
RNA information. You can download a wgRna DB-file for human NCBI36 genome 
assembly using:

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/wgRna.txt.gz
	gunzip wgRna.txt.gz
	sort -k 2,2 -k 3,3n wgRna.txt > hg18_wgRna.sorted

The first a few lines of a sorted wgRna table is shown below:

	593     chr1    1092346 1092441 hsa-mir-200b    960     +       1092402 1092425 miRna
	593     chr1    1093105 1093195 hsa-mir-200a    960     +       1093120 1093142 miRna
	593     chr1    1093105 1093195 hsa-mir-200a    960     +       1093158 1093180 miRna
	593     chr1    1094247 1094330 hsa-mir-429     960     +       1094297 1094319 miRna
	611     chr1    3467118 3467214 hsa-mir-551a    480     -       3467133 3467154 miRna

The --score_threshold and --normscore_threshold arguments both operate on the 
sixth field in the line (such as "960").

=item * B<SegDup format>

The segDup table from UCSC genome browser contains chromosome regions with 
segmental duplications, as well as the "target" regions that match these 
duplications. You can download a segdup DB-file for human NCBI36 genome assembly 
using:

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/genomicSuperDups.txt.gz
	gunzip genomicSuperDups.txt.gz
	sort -k 2,2 -k 3,3n genomicSuperDups.txt.gz > hg18_genomicSuperDups.sorted

The first a few lines of a sorted segdup table is shown below:

	585     chr1    465     30596   No.1139,chr2:114046528  570211  -       chr2    114046528       114076216       243018229       1139    1000    N/A     Filtered        N/A     N/A     build35/align_both/0012//both064929     30176   43      531     29645   29282   363     128     235     0.987755        0.986324        0.012346        0.0123574
	585     chr1    486     30596   No.2251,chr9:844        582086  +       chr9    844     30515   138429268       2251    1000    N/A     Filtered        N/A     N/A     build35/align_both/0013//both065185     30137   14      491     29646   29474   172     64      108     0.994198        0.993729        0.00582435      0.00582657

The --score_threshold and --normscore_threshold arguments both operate on the 
sixth field in the line (such as "570211").

=item * B<knownGene format>

The knownGene table from UCSC genome browser contains gene identifier (name), 
location, exon count and location, Swiss-Prot identifiers. It does not contain 
the gene symbol information, so an additional step is necessary to interrogate 
the kgXref table to find the gene symbol for each gene identifier. You can 
download a knownGene DB-file for NCBE36 human genome assembly using:

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/knownGene.txt.gz
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/kgXref.txt.gz
	gunzip knownGene.txt.gz
	gunzip kgXref.txt.gz
	mv knownGene.txt hg18_knownGene.txt
	mv kgXref.txt hg18_kgXref.txt

Unlike previous sections, a sorting of chrosome and location is not necessary 
for knownGene or refGene tables.

To scan a query file against the knownGene table, use:

	scan_region.pl --knowngene queryfile knownGene.txt

The output is:

	chr3:2000-349990        uc003bot.1,uc003bou.1,uc003bov.1,uc003bow.1     0
	chr19:32333-52333       NOT_FOUND       NOT_FOUND

Unless the --condense_query argument is set, the output will be in the same 
order as the query file. Those query that does not match any knownGene will 
still be printed, although "NOT_FOUND" will show as the last 2 columns in the 
output. If there is a match with the DB-file, the output line contains gene 
identifiers and their distances (should be ZERO here) that overlap with query 
location.

In the above output, four genes all overlap with the query region, so they are 
printed out together but separated by comma. In comparison, no genes are 
present that overlap with the second region.

To use gene symbols in the generated output, use:

	scan_region.pl --knowngene --kgxref hg18_kgXref.txt queryfile hg18_knownGene.txt

The output is:

	chr3:2000-349990        CALL,CHL1       0
	chr19:32333-52333       NOT_FOUND       NOT_FOUND

To identify surrounding genes to the query location, use:

	scan_region.pl --knowngene --kgxref hg18_kgXref.txt -expandmax 1m queryfile hg18_knownGene.txt

The output is:

	chr3:2000-349990        CALL,CHL1       0
	chr19:32333-52333       F379    3647

This command will expand the query location by 1m (1 million base pairs) in both 
sides (use --expandleft and --expandright if you only want to find overlapping 
genes to the left or right of the query location), and print out closest 
overlapping genes (in this case, F379) as well as their distance (in this case, 
3647bp) to the query. It is very possible that after expanding the query for 1m 
base pairs, multiple genes will be found, but only the CLOSEST gene will be 
printed in the output.

If you want to print out all genes within 1m of the query, then the --expanddb 
argument should be used, which will expand the chromosome region in the DB-file 
to find query-DB match. The command is:

	scan_region.pl --knowngene --kgxref hg18_kgXref.txt queryfile hg18_knownGene.txt -expanddb 1m

The output is:

	chr3:2000-349990        AK126307,BC065754,CALL,CHL1,CNTN6,hNB-3 0
	chr19:32333-52333       ABCA7,ABCA7/ABCA-SSN,AK311622,ARID3A,AZU1,BC048429,BSG,BX647595,C19orf20,C19orf21,C19orf22,C19orf6,CDC34,CFD,CNN2,DQ574670,DQ578983,DQ599872,DQ600587,ELA2,F379,FAM148C,FGF22,FLJ00038,FLJ00277,FSTL3,GRIN3B,GZMM,HCN2,HMHA1,KISS1R,LPPR3,MADCAM1,MED16,MIER2,ODF3L2,OR4F17,PALM,PHP2,POLR2E,POLRMT,PPAP2C,PRG2,PRSSL1,PRTN3,PTBP1,RNF126,SHC2,THEG,UNQ674,WDR18,hEMMPRIN       0

In the above result, you can alternatively think that all genes within the 1mb 
region of the query is printed: this is the major differnece between --expandmax 
and --expanddb.

=item * B<refGene format>

The refGene table from UCSC Genome Browser contains genome location of RefSeq 
mRNA (usually with NM_ suffix in its identifier). You can download a refGene 
template file using:

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/refGene.txt.gz
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/refLink.txt.gz
	gunzip refGene.txt.gz; mv refGene.txt hg18_refGene.txt
	gunzip refLink.txt.gz; mv refLink.txt hg18_refLink.txt

The syntax for scanning refGene table is similar to that scanning the knownGene 
table. However, some refgene table (such as hg18) contains more information than 
other refgene table from earlier genome assembly (such as hg17). In that case, 
the --name2 argument can be used to replace the gene name (the "name" field in 
the refGene table) in output by the alternative gene name ("name2" field in the 
refGene table). Having said that, to reproduce the output gene names in the UCSC 
genome browser, it is best to still use the refLink.txt file.

To scan a query file against the refGene table:

	scan_region.pl --refgene --name2 queryfile hg18_refgene -expandmax 100k

The output is:

	chr3:2000-349990        CHL1    0
	chr19:32333-52333       OR4F17  9346

To use the gene name in the refLink.txt file to annotate the output:

	scan_region.pl --refgene --reflink hg18_refLink.txt queryfile hg18_refGene.txt -expandmax 100k

The output is the same in this case:

	chr3:2000-349990        CHL1    0
	chr19:32333-52333       OR4F17  9346

Generally speaking, refGene table contains less genes than the knownGene table, 
but the genes in refGene is better annotated and usually have gene symbol 
associated with them.

=item * B<refCds format>

In some cases, one might be particularly interested in finding overlap between a 
query and the coding sequence of a gene (from the first exon to the last exon), rather than the entire transcript. In 
this case, a refGene file can be provided, but the --refcds argument should be 
specified.

=item * B<refExon format>

In many cases, especially when analyzing copy number variations, one might be 
particularly interested in finding overlap between a query and the exon of a 
gene, excluding all introns. In this case, a refGene file can be provided, but 
the --refexon argument should be specified.

=item * B<phastCons format>

Despite its name, do not confuse it with the MCE (most conserved element) track 
of the genome browser. The phastcons file merely contains conservation scores 
for each genome position that can be aligned to other species (in comparison, 
the MCE file contains the conservation scores for the top 5% most conserved 
genomic regions in the genome). For human genome, much more than 5% of the 
genome (probably ~50%) can be actually aligned, although only 5% of the genome 
can be considered as conserved during evolution.

In other word, the MCE annotates conserved sequence segments, but the phastcons 
annotates individual base pairs that can be aligned.

Since the file contains a score for each genomic position, the file size is 
extremely large (typically >10GB), and the scanning takes a lot of time (up to a 
whole day). In most circumstances, I recommended using the MCE file for sequence 
conservation analysis, since the phastcons score is less meaningful than the MCE 
score, which use a hidden Markov model to identify and summarize a genomic 
region that is conserved during evolution.

=item * B<generic format>

The program can also take a generic UCSC genome browser format, where the first 
four columns are index, chromosome, start and end position, with zero-based 
(half-open) position counting schemes. This is indeed how the MCE, EvoFold and 
other tracks are organized in the database for the newest genome build (but not 
necessarily the case for older genome build, and not necessarily true for 
non-human species).

=back

Please contact Kai (kai@openbioinformatics.org) for questions or concerns on the 
PennCNV program.

=cut                                       