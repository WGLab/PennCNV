#!/usr/bin/perl -w
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: 35 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2008-10-03 08:57:47 -0700 (Fri, 03 Oct 2008) $';

our ($verbose, $help, $man);
our ($inputfile, $operator, @number, %number, $numfile, $reverse, $tolerate, $tab, $heading, $output, $start_split, $end_split, $name_by_header, 
	$topleft, $wordchar, $stopchar, $beforestring, $dryrun, $filenameout);

GetOptions ('verbose'=>\$verbose, 'help'=>\$help, 'man'=>\$man, 'numfile=s'=>\$numfile, 'reverse'=>\$reverse, 'tolerate'=>\$tolerate, 'tab'=>\$tab, 
	'heading=i'=>\$heading, 'output=s'=>\$output, 'start_split=i'=>\$start_split, 'end_split=i'=>\$end_split, 'name_by_header'=>\$name_by_header,
	'topleft=s'=>\$topleft, 'wordchar=s'=>\$wordchar, 'stopchar=s'=>\$stopchar, 'dryrun'=>\$dryrun, 'filenameout=s'=>\$filenameout, 
	'beforestring=s'=>\$beforestring) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV >= 2 or pod2usage ("Syntax error");

($inputfile, $operator, @number) = @ARGV;

if (defined $numfile) {
	@number and pod2usage ("Error in argument: number already defined in command line, so no need to use --numfile option");
	open (NUMBER, $numfile) or confess "Error opening numfile $numfile: $!";
	$_ = <NUMBER>;
	s/[\r\n]+$//;
	s/^\s+|\s+$//g;
	@number = split (/\s+/, $_);
	close (NUMBER);
} else {
	@number or pod2usage ("Error in argument: you must specify at least one column");
}

if (defined $wordchar) {
	length ($wordchar) == 1 or pod2usage ("Error in argument: --wordchar <$wordchar> can only be a single character/letter");
	$stopchar || $beforestring and pod2usage ("Error in argument: please do not specify --stopchr or --beforestring when --wordchar is specified");
}

if (defined $stopchar) {
	$beforestring and pod2usage ("Error in argument: please do not specify --beforestring when --wordchar is specified");
	length ($stopchar) == 1 or pod2usage ("Error in argument: --stopchar <$stopchar> can only be a single character/letter");
	defined $wordchar and pod2usage ("Error in argument: --stopchar override the --wordchar option");
	if ($stopchar eq '.' or $stopchar eq '(' or $stopchar eq ')' or $stopchar eq '[' or $stopchar eq ']' or $stopchar eq '{' or $stopchar eq '}' or $stopchar eq '$' or $stopchar eq '^' or $stopchar eq '|' or $stopchar eq '*' or $stopchar eq '+' or $stopchar eq '?' or $stopchar eq "\\") {
		$stopchar = "\\" . $stopchar;
	}
}

if (defined $beforestring) {
	$beforestring =~ s/\./\\./g;
	$beforestring =~ s/\-/\\-/g;
	$beforestring =~ s/\?/\\?/g;
	$beforestring =~ s/\+/\\+/g;
	$beforestring =~ s/\*/\\*/g;
	$beforestring =~ s/\$/\\\$/g;		#special handling here required
	$beforestring =~ s/\^/\\^/g;
	$beforestring =~ s/\|/\\|/g;
	$beforestring =~ s/\(/\\(/g;
	$beforestring =~ s/\)/\\)/g;
	$beforestring =~ s/\[/\\[/g;
	$beforestring =~ s/\]/\\]/g;
	$beforestring =~ s/\{/\\{/g;
	$beforestring =~ s/\}/\\}/g;
	$verbose and print STDERR "NOTICE: Changing --beforestring argument to $beforestring\n";
}

$operator eq 'f' and $operator = 'fetch';
$operator eq 'a' and $operator = 'append';
$operator eq 's' and $operator = 'split';
$operator eq 'as' and $operator = 'appendstring';
$heading ||= 0;

my @old_number;
for (@number) {
	if (/^(\d+)-(\d+)$/) {
		push @old_number, ($1 .. $2);
	} else {
		push @old_number, $_;
	}
}
@number = @old_number;

if ($operator eq 'split') {
	@number == 1 or pod2usage ("Error in argument: when operator is split, only one number is allowed that specifies the number of columns in each output");
	$output or pod2usage ("Error in argument: please specify the --output argument for the split operation");
}

#get the maximum column so that later we can check whether the input line has sufficient number of columns
@old_number = sort {$b <=> $a} @number;

if ($operator ne 'appendstring') {
	for (@number) {
		m/^\d+$/ or pod2usage ("<number> should be a positive integer or zero, but your specified <$_>");
		$number{$_} = 1;		#a hash storing all the column numbers
		$reverse and $_ eq '0' and pod2usage ("Error in argument: you cannot use column 0 together with the --reverse option");
	}
}

if ($inputfile eq 'stdin') {
	$operator eq 'append' and pod2usage ("Error in argument: you cannot append to stdin using input from stdin");
	*FH = *STDIN;
} else {
	open (FH, $inputfile) || confess "Error: cannot read from inputfile $inputfile: $!";
}

my (@fhout, @fileout, @heading);								#file handle for a bunch of output files (when operator is 'split')
while (my $nextline = <FH>) {
	$nextline =~ s/[\r\n]+$//;

	my @content;
	if ($tab) {
		@content = split (/\t/, $nextline);			#split the line by the tab character
	} else {
		$nextline =~ s/^\s+|\s+$//g;
		@content = split (/\s+/, $nextline);			#split the line by any blank character, and delete all heading and trailing spaces
	}
	my $newline = '';
	if ($operator eq 'append') {
		$newline = <STDIN>;
		defined $newline or confess "Error: STDIN has less lines than those in $inputfile";
		$newline =~ s/\s*[\r\n]+$//;
	}
	if ($operator eq 'append' or $operator eq 'fetch') {
		if ($reverse) {		#only print those unspecified column (reverse column selection)
			for my $i (0 .. @content-1) {
				$number{$i+1} or $newline .= "\t$content[$i]";
			}
		} else {
			for my $i (0 .. @number-1) {
				if ($number[$i]) {
					if ($number[$i] > @content) {
						$tolerate or confess "Error: too few columns (${\(scalar @content)}) in input line <$nextline> for column $number[$i]";
						$newline .= "\t";
					} else {
						$newline .= "\t" . $content[$number[$i]-1];
					}
				} else {
					$newline .= "\t" . join ("\t", @content);
				}
			}
		}
		$newline =~ s/^\t//;
		$newline .= "\n";
		print $newline;
	} elsif ($operator eq 'appendstring') {
		print $nextline, "\t", join ("\t", @number), "\n";
	} elsif ($operator eq 'split') {
		if ($topleft and not @fhout) {				#replace the top left cell of the file by this string
			$content[0] = $topleft;
		}
		
		if ($heading) {
			@heading = splice (@content, 0, $heading);
		}
		if ($end_split) {					#delete anything after --end_split
			if (@content >= $number[0]*$end_split) {
				splice (@content, $number[0]*$end_split);
			}
		}

		if ($start_split) {					#delete anything before --start_split
			if (@content <= $number[0] * ($start_split-1)) {
				confess "Error: --start_split argument is too small, as there are only ${\(scalar @content)} columns in the row";
			}
			splice (@content, 0, $number[0]*($start_split-1));
		}
		if (not @fhout) {					#reading the first line of the input stream
			@content % $number[0] == 0 or print STDERR "WARNING: the input contains ", scalar (@content)+$heading, " columns, by excluding $heading heading columns, the rest cannot be split to $number[0] whole parts. The last file will contain fewer columns!!!\n";
			my $num_fh = int (@content / $number[0]);
			@content % $number[0] == 0 or $num_fh++;	#one additional file will be printed but it will have only partial columns
			my (%opened_filename);				#to prevent files from having same name (sometimes header line contains same samples ID)
			for my $i (1 .. $num_fh) {
				#if ($start_split) {
				#	$i += ($start_split-1);
				#}
				my ($temp_fh, $temp_output);
				if ($name_by_header) {			#output file name is based on the corresponding column in header line
					my $firstword;
					if ($wordchar) {
						$content[($i-1)*$number[0]] =~ m/^([\w$wordchar]+)/ or confess "Error: --name_by_header argument cannot be applied, since the header (" . $content[($i-1)*$number[0]] . ") must start with words or --wordchar <$wordchar>\n";
						$firstword = $1;
					} elsif ($stopchar) {
						$content[($i-1)*$number[0]] =~ m/^([^$stopchar]+)/ or confess "Error: --name_by_header argument cannot be applied, since the header (" . $content[($i-1)*$number[0]] . ") must start with letters other than --stopchar <$stopchar>\n";
						$firstword = $1;
						if ($firstword =~ m/\s/) {
							print STDERR "WARNNING: Sample name $firstword contains space: change to underline instead\n";
							$firstword =~s/\s/_/g;
						}
					} elsif ($beforestring) {
						$content[($i-1)*$number[0]] =~ m/(.+?)$beforestring/ or confess "Error: --beforestring argument cannot be applied, since the header (" . $content[($i-1)*$number[0]] . ") must contain the string ($beforestring) that you specified\n";
						$firstword = $1;
						if ($firstword =~ m/\s/) {
							print STDERR "WARNNING: Sample name $firstword contains space: change to underline instead\n";
							$firstword =~s/\s/_/g;
						}
					} else {
						$content[($i-1)*$number[0]] =~ m/^(\w+)/ or confess "Error: --name_by_header argument cannot be applied, since the header (" . $content[($i-1)*$number[0]] . ") must start with words";
						$firstword = $1;
					}
					if ($opened_filename{$firstword}) {
						while (1) {
							-f "$output.$firstword" . "_R$opened_filename{$firstword}" or last;	#make sure such file does not exist in current directory (otherwise they will be overwritten)
							$opened_filename{$firstword}++;
						}
						$temp_output = "$output.$firstword" . "_R$opened_filename{$firstword}";
						print STDERR "WARNING: Multiple samples with same identifier $firstword found: file renamed as $temp_output\n";
					} else {
						$temp_output = "$output.$firstword";
						if (-f $temp_output) {
							$opened_filename{$firstword} = 1;
							while (1) {
								-f "$output.$firstword" . "_R$opened_filename{$firstword}" or last;	#make sure such file does not exist in current directory (otherwise they will be overwritten)
								$opened_filename{$firstword}++;
							}
							print STDERR "WARNING: A file named $temp_output already exists: rename the output file to $output.${firstword}_R$opened_filename{$firstword}\n";
							$temp_output = "$output.$firstword" . "_R$opened_filename{$firstword}";
						}
					}
					
					$opened_filename{$firstword}++;			#add one to the opened file names
				} else {
					my $file_index = $i;
					$start_split and $file_index += ($start_split-1);
					$temp_output = "$output.split$file_index";
				}
				
				eval {
					open ($temp_fh, ">$temp_output") or die "$!\n";
				};
				if ($@) {
					print STDERR "ERROR: cannot open $temp_output for writting after opening ${\(scalar @fhout)} files: $@";
					print STDERR "POSSIBLE SOLUTION: try to use --start_split and --end_split argument to limit the number of files to be processed at once\n";
					for my $i (0 .. @fhout-1) {
						close ($fhout[$i]);				#close all opened file handles
						unlink ($fileout[$i]);				#delete all opened files
					}
					exit (100);
				}
				push @fhout, $temp_fh;
				push @fileout, $temp_output;
			}
			if ($name_by_header) {
				print STDERR "NOTICE: splitted content will be written to ${\(scalar @fhout)} files with name starting as $output and ending as the sample name (as specified in first line of input file)\n";
			} elsif ($start_split) {
				print STDERR "NOTICE: splitted content will be written to files $output.split$start_split to $output.split", $start_split+$num_fh-1, "\n";
			} else {
				print STDERR "NOTICE: splitted content will be written to files $output.split1 to $output.split$num_fh\n";
			}
			
			if ($dryrun) {
				for my $i (0 .. @fhout-1) {
					close ($fhout[$i]);				#close all opened file handles
					unlink ($fileout[$i]);				#delete all opened files
				}
			}
			
			if ($filenameout) {
				open (FILENAME, ">$filenameout") or confess "Error: cannot write to file name out $filenameout\n";
				print FILENAME join ("\n", @fileout), "\n";
				close (FILENAME);
			}
			
			if ($dryrun) {
				print STDERR "NOTICE: dry run of the program finished with ${\(scalar @fileout)} output files\n";
				exit (200);
			}
		}
		
		for my $i (0 .. @fhout-1) {
			my $temp_fh = $fhout[$i];
			$heading and print $temp_fh join ("\t", @heading), "\t";
			if ($i*$number[0]+$number[0] <= @content) {
				print $temp_fh join ("\t", @content[($i*$number[0]) .. ($i*$number[0]+$number[0]-1)]), "\n";
			} else {
				print $temp_fh join ("\t", @content[($i*$number[0]) .. (@content-1)]), "\n";
			}
		}
	} else {
		pod2usage ("Error in argument: operator can only be append or fetch");
	}
}

if ($operator eq 'append') {
	$_ = <STDIN>;
	defined $_ and print STDERR "WARNING: STDIN has more lines than those in $inputfile";
}

=head1 SYNOPSIS

 kcolumn.pl <inputfile> <operator> <operated-numbers>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	-n, --numfile <string>		file containing operated number
 	-r, --reverse			reverse the column selection
 	    --tolerate			tolerate lines with too few columns
 	    --tab			use tab as field separator (default=space+tab)
 	    --heading <int>		number of heading columns
 	    --output <string>		output filename for split operation
 	    --start_split <int>		start processing file splitting from this file index
 	    --end_split <int>		end processing file splitting after this file index
 	    --name_by_header		name outputfiles by the first word in header columns
 	    --topleft <string>		replace the top left cell by specified string
 	    --wordchar <characters>	treat the characters as word in file naming
 	    --stopchar <characters>	all characters before this char is used in file naming
 	    --beforestring <string>	all characters before this string is used in file naming
 	    --dryrun			test program for number of simultaneous files without output
 	    --filenameout <file>	write all splitted file names to this file sequentially

 Function: perform column extraction operations (fetch or append or split) on 
 inputfile. All column index start from 1, not zero. Zero means "all columns" in this 
 program.
 
 Example: kcolumn.pl inputfile3 fetch 1 3
          cat inputfile2 | kcolumn.pl inputfile3 append 0
          kcolumn.pl inputfile split 3 -head 3 -tab

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--numfile>

specify a file containing numbers that is operated by the operator

=item B<--reverse>

reverse the column selection, so only those columns that are not 
selected are printed

=item B<--tolerate>

Tolerate input lines that have insufficient columns (normally this 
program will throw an error when such line is encountered).

=item B<--heading>

number of heading columns. "heading column" refers to the first a few columns in 
each line that will be retained in all splitted output files.

=item B<--start_slit>

specify the starting index of splitting files. When there are too much columns 
in the inputfile so that split file is more than 1020, many operating systems 
will fail to open new file handles. Specifying the start_split and end_split can 
deal with this problem. For example, an inputfile contains 300 columns, with 3 
heading columns and 10 columns per splitted file, specifying --start_split as 
100 and --end_split as 200 will generate splitted files for the 1004-20003 
columns.

=item B<--end_split>

specify the ending index of splitting files.

=item B<--name_by_header>

name the splitted output files by the corresponding column in the header line, 
as opposed to name the outputfiles by numbers.

=item B<--topleft>

replace the top left cell of the file by this specified string, since in many 
cases an input file does not have a top-left cell (the first word in a file)

=item B<--wordchar>

treat the specified characters/letters as word characters, for use together with 
the -- name_by_header argument in output file naming. By default, only 
alphabetical characters and underline characters are treated as word characters, 
and only the first word in each column is used as the output file name.

=item B<--stopchar>

treat the specified characters as stop characters, for use together with the --
name_by_header argument in output file naming. By default, only alphabetical 
characters and underline characters are treated as word characters and are used 
in file naming, but once --stopchar is specified, anything other than this 
character will be treated as a word, and such word will be used as part of the 
output file name.

=back

=head1 DESCRIPTION

This program is used to perform various column operations on a text file. Column 
is defined as text strings separated by empty character (either space or tab or 
return) or tab character in each line (when -tab argument is set). Currently three 
operations are supported: fetch, append and split. In the future, more column 
operations will be implemented and added to this program. 

The "fetch" operation retrieves specific columns from a file and displays the 
content. The "append" operation fetched columns from an inputfile and then 
append the content to STDIN stream. The "split" operation split the inputfile 
into several output files, each containing certain columns in the input file.

The operated-numbers (column specifiers) can be expressed in two formats: (1) a 
list of numbers such as "1 3 7" (2) a range of numbers such as "1 3-5 8-10" 
(equivalent to "1 3 4 5 8 9 10"). Alternatively, you can put all numbers in a 
single line in a file separated by spaces, and use the --numfile option to 
specify the file to the program. ALL NUMBERS MUST BE POSITIVE INTEGER OR ZERO, 
AND ZERO REPRESENTS ALL COLUMNS.

=head2 B<fetch operation>

The fetch operation is identical to the "cut" program in many Linux 
distributions. However, by default the blank character (such as space and tab) 
is used in defining columns, instead of the tab character which is used in the 
cut program. In addition, the cut program can only display columns in sorted 
order, but kcolumn.pl program can be very flexible to display desired columns in 
specified order. For example, the content of a file called "inputfile" is below:

	11       10     red
	12       22     yellow
	13       37     green
	14       50     black
	15       70     blue

When using the command

	kcolumn.pl inputfile f 3 1

the output will be:

	red     11
	yellow  12
	green   13
	black   14
	blue    15

=head2 B<append operation>

The append operation is used to fetch and append to a data stream received from 
STDIN. For example, using the command

	perl -e 'print "good\n" for (1..5)' | kcolumn.pl inputfile a 3

will generate the output:

	good	red
	good	yellow
	good	green
	good	black
	good	blue

The first column in each line (good) is printed from the first command, while 
the second column in each line is retrieved from the third column of the 
inputfile and appended to the output line.

=head2 B<split operation>

The split operation is used to split an input file to multiple files, each file 
contains certain columns from the input file. For example, suppose the input 
file contains 10 columns, and the first column represent the "heading" column 
that should be kept in each individual file. Now if we want to split the 
inputfile into three separate files, we can use

	kcolumn.pl inputfile split 3 -heading 1 -out newfile

The resulting three files (file names are newfile.split1, newfile.split2 and 
newfile.split3, respectively) will each contain four columns from the inputfile, 
including one heading column.

When the --name_by_header argument is set, the output file names will start with 
newfile, but with a prefix that is determined by the header line (first line) of 
the inputfile. For the example above, the output file names will be newfile.11, 
newfile.10 and newfile.red.

When the inputfile contains many columns to be splitted into >1000 individual 
files, usually the operating system will complain that there are too many open 
file handles. In this case, one can use the --start_split and --end_split 
argument to overcome this limitation and process selected columns. For example, 
for an input file that will generate 2500 individual files after splitting, one 
can use:

	kcolumn.pl inputfile split 3 -start 1 -end 1000 -out newfile
	kcolumn.pl inputfile split 3 -start 1001 -end 2000 -out newfile
	kcolumn.pl inputfile split 3 -start 2001 -end 3000 -out newfile

The above three commands will generate the 2500 individual files (note that it 
is okay to use higher number of --end_split than the actual number of columns, 
as the program will figure the correct column automatically).

=head2 B<other column operations>

In the future, more column operations will be added to this program.

=cut
                                                                                                                                                                     