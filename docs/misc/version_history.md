Older versions of PennCNV are given below together with one-sentence description of changes. It is highly recommended to use only the latest version. More detailed change log is recorded in the program code per se.

- Latest version (2011Jun16): cumulative bug fixes and function enhancement. Added clean_cnv.pl and cal_gc_snp.pl scripts to clean CNV calls and for calculate GC content for user-supplied genome builds. Additionally, pre-compiled executables are included for 32-bit Perl 5.8.9, 5.10.1, 5.12.3 and 5.14.0 in Windows XP, Vista and 7.
- 2010May01 version: reduce memory usage for PennCNV such that Affy6 array requires <2GB memory for CNV calling. Added compile_pfb.pl program to generate users' own PFB files given a list of signal intensity files. Added functionality to plot the signal intensity values for each CNV call for visual validation of reliability of CNV calls. Re-write genomic wave adjustment procedure to solve compatibility issues in certain system architectures during compilation.
- 2009Aug27 verion: minor bug fix, added --lastchr argument to detect_cnv.pl to handle non-human arrays. Added infer_snp_allele.pl program to infer CNV-based SNP genotypes, or to validate de novo CNVs and assign P-values. Reorganize kext/ directory structure to accormodate different Perl versions. Fix the missing "-" before "minsnp" problem in the BeadStudio/GenomeStudio plug-in. Enhanced functionality of convert_cnv.pl, which now handles XML files exported from BeadStudio/GenomeStudio. Added -reciprocal argumen to the compare_cnv.pl program to fine-tune -minoverlap argument. Updated scan_region.pl program for better functionality and more accurate control of overlapping criteria. See "PennCNV main package" section above for link to download files.
- 2008Nov19 version: adding functionality in the filter_cnv.pl and compare_cnv.pl program.
- 2008Jun26 version: GC-model signal pre-processing to handle low-quality samples; family-based CNV calls by joint-calling algorithm; other enhancements
- 2008Mar11 version: better compatibility with BeadStudio, beta-version of CNV visualization
- 2008Mar03 version: beta-version for using PennCNV with BeadStudio plug-in, beta-version for signal adjustment for wavy samples
- 2008Feb15 version: minor changes with more QC summary
- 2007Dec14 version: support 64bit system and support cygwin now, print sample quality measure, produce confidence score by --conf argument (experimental feature)
- 2007Nov13 version: cumulative minor bug fix and function enhancement
- 2007Oct31 version: fix bugs for chrX processing; fix bugs in kcolumn.pl for splitting huge files
- 2007Sep28 version: re-implement calling algorithm to process each chromosome separately; change output format

