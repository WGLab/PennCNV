1. **How to ask a question?**

    Drop me an email as which should contain (1) command line argument (2) error message in screen or the LOG file (3) sometimes example inputfile (4) in case you use Mac or Windows, let me know. The reason is that I can fix something or diagnose something only if I can understand the question and reproduce the results.

1. **How to process arrays other than HumanHap550?**

    The PennCNV was originally developed for Illumina HumanHap550 array (in the Genome Research paper) but has been extended to other arrays. In PennCNV, the hhall.hg18.pfb, hhall.hmm and hhall.hg18.gcmodel files have been provided in the lib/ folder. They can be used to process arrays such as HumanHap300, HumanHap550, HumanHap650, HumanHap1M, HumanCNV370, Human610, Human1. For Affymetrix arrays, the PennCNV-Affy package contains library files for commonly used marker sets. In a nut shell, PennCNV can process an array if markers for this array are annotated in the PFB file.

    For other Illumina&Affymetrix SNP arrays (for example, the Illumina Omni 2.5M array, or the Affy cytogenetics array), the user needs to compile your own PFB file that specifies the chromosome coordinates and PFB value for each marker in the array. Use the compile_pfb.pl program included in the PennCNV package for doing that.

    For other types of arrays, such as oligonucleotide arrays, and for situations where user does not have access to the raw data but can access processed signal values, please check Input Files section in the tutorial for potential solutions. Besides commonly used SNP arrays, users have also reported success on Agilent non-SNP array, Perlegen SNP array, Affymetrix 100K SNP arrays.

1. **How to make a PFB file?**

    The PFB file is described in the Input Files section in the tutorial. If you already have hundreds of signal files (with LRR/BAF values) generated on non-human species, or on non-European populations, then you can use the compile_pfb.pl program to generate PFB file from a collection of signal files.

1. **Does PennCNV handle non-polymorphic probes?**

    Yes. In principle, it can also handle oligonucleotide arrays without SNP markers. There are two ways to represent non-polymorphic probes: (1) if the marker name contains the word cnv (for Illumina arrays) or CN (for Affymetrix arrays), they will be treated as non-polymorphic markers by PennCNV. (2) if the allele frequency annotation in the PFB file for a marker is more than 1, then this marker will be treated as a non-polymorphic marker by PennCNV program. The second situation is suitable for custom designed arrays.

1. **Is male chromosome X similar to heterozygous deletions in autosomes?**

    No. One may think that by replacing the signal intensity of a segment of autosome by intensity of male chrX, he/she can create an artificial heterozygous deletion in autosomes, which is not correct. Autosome and chrX have different signal intensity distributions. When building reference clustering files, Illumina used both males and females: as a result, the expected ZERO value of LRR in chrX corresponds to ~1.5 copy of chrX, not 2 copies!

    To further understand this, one can open BeadStudio and examine by eye the mean intensity of chrX from multiple female individuals. The mean is not ZERO, but higher than ZERO, indicating that female chrX is not even similar to 2-copy autosomes. For more deails and discussions, refer to Sup Figure 1 of the 2007 PennCNV paper in Genome Research.

1. ** How to use genomic wave adjustment independent of CNV calling?**

    Some users just want to adjust signal intensity values, without generating CNV calls by PennCNV. The genomic_wave.pl program in PennCNV package can be used to adjust signal intensity values. The input file must have a field in the header line that says "*.Log R Ratio". The -adjust argument can be used to generate a new file with updated Log R Ratio measures. This procedure can be also used in Agilent arrays or Nimblegen arrays for adjustment. Email me for a script to generate GC model file for these custom arrays.

1. **A sample generates >1000 CNV calls, whats wrong?**

    Most probably it has low signal quality. Please check the QC&Annotation section in the tutorial to see how to take advantage of the filter_cnv.pl program for automatic QC analysis.

    It could be also due to heterosomatic CNVs where a fraction of cells lose/gain one chromosomal regions (such as one chromosome arm), so PennCNV gives lots of CNV calls within this region specifically. The 2007 Genome Research paper on PennCNV have some discussions on that, with some illustrations in supplementary figures. In this case, just delete this chromosome in this individual from analysis.

1. **How to identify large CNV calls (>10SNPs, >100kb)?**

    This can be done by the following command: `filter_cnv.pl sampleall.rawcnv -numsnp 10 -length 100k -out sampleall.largecnv`

1. **How to compare PennCNV calls with other algorithms?**

    If a user also generates CNV calls by cnvPartition or QuantiSNP in BeadStudio/GenomeStudio, you can export the CNV calls to a XML file (see illustration described in this tutorial), then use the convert_cnv.pl program (with "-intype xml -outtype penncnv" argument) to convert the XML file to PennCNV format. If a user also generates CNV calls by Birdseye program, the convert_cnv.pl program can also convert BirdSeed calls to standard PennCNV format. Next, you can use the compare_cnv.pl program (use -m argument to read the manual) to compare the CNV calls generated by different programs that are all in PennCNV format. A two-column list file should be prepared that specifies the file name correspondence in PennCNV call file and the call file generated by the other algorithm (but converted to PennCNV format).

1. **How to remove CNV calls in immunoglobulin regions?**

    The CNV calls in immunoglobulin regions are most likely cell line artifact, so they should be removed as part of the QC procedure. The `scan_region.pl` program can be used to do this conveniently: `scan_region.pl cnvcall imm_region -minqueryfrac 0.5 > cnvcall.imm; fgrep -v -f cnvcall.imm cnvcall > cnvcall.clean`

    This command first scan the cnvcall file against known immunoglobulin regions, and any CNV call that overlap with immunoglobulin regions are written to the cnvcall.imm file (the --minqueryfrac means that at least 50% of the length in the CNV call must overlap with the immunoglobulin region, to exclude cases where a very large CNV call happens to encompass the immunoglobulin regions). Then the fgrep program is used to remove these regions from the file and generate a cleaned cnvcall.clean file. The imm_region file contains immunoglobulin regions. For the 2006 human genome assembly, these four regions can be put into the file:

    
    ```
 chr22:20715572-21595082
 chr14:105065301-106352275
 chr2:88937989-89411302
 chr14:21159897-22090937
    ```
    

1. **How to remove CNV calls in centromeric and telomeric regions?**

    The same techniques described above can be used. For telomeric regions, one can treat the 100kb or 500kb region within start or end of chromosome as telomeric region. For example, for the 500kb threshold, you can put the following regions ino a file and then use scan_region.pl to remove CNV calls:

    ```
chr1:1-500000
chr1:246749719-247249719
......
chr22:1-500000
chr22:49191432-49691432
    ```

    For centromeric regions, the following definition can be used (NCBI36 2006 human genome assembly). In fact, you may want to add 100kb (or 500kb) to both the left and right of these regions, just to make sure that centromeric CNVs are identified comprehensively.

    ```
    chr1:121100001-128000000
    chr2:91000001-95700000
    chr3:89400001-93200000
    chr4:48700001-52400000
    chr5:45800001-50500000
    chr6:58400001-63400000
    chr7:57400001-61100000
    chr8:43200001-48100000
    chr9:46700001-60300000
    chr10:38800001-42100000
    chr11:51400001-56400000
    chr12:33200001-36500000
    chr13:13500001-18400000
    chr14:13600001-19100000
    chr15:14100001-18400000
    chr16:34400001-40700000
    chr17:22100001-23200000
    chr18:15400001-17300000
    chr19:26700001-30200000
    chr20:25700001-28400000
    chr21:10000001-13200000
    chr22:9600001-16300000
    chrX:56600001-65000000
    chrY:11200001-12500000
    ```

    For centromeric regions, the following definition can be used (NCBI37/hg19 human genome assembly). In fact, you may want to add 100kb (or 500kb) to both the left and right of these regions, just to make sure that centromeric CNVs are identified comprehensively.

    ```
    chr1:121500000-128900000
    chr2:90500000-96800000
    chr3:87900000-93900000
    chr4:48200000-52700000
    chr5:46100000-50700000
    chr6:58700000-63300000
    chr7:58000000-61700000
    chr8:43100000-48100000
    chr9:47300000-50700000
    chr10:38000000-42300000
    chr11:51600000-55700000
    chr12:33300000-38200000
    chr13:16300000-19500000
    chr14:16100000-19100000
    chr15:15800000-20700000
    chr16:34600000-38600000
    chr17:22200000-25800000
    chr18:15400000-19000000
    chr19:24400000-28600000
    chr20:25600000-29400000
    chr21:10900000-14300000
    chr22:12200000-17900000
    chrX:58100000-63000000
    chrY:11600000-13400000
    ```

1. **Does chromosome X requires special handling?**

    Yes. By default PennCNV only works on autosomes, without generating CNV calls for chrX and chrY. Several important issues for sex chromosome CNV calling are:

    The `-chrx` argument should be used for chrX CNV calling. Whenever possible, the `--sexfile` argument should be used to specify sample sex. This is a tab-delimited two-column file, containing signal file names and sex (male or female).
PennCNV tries to predict sample sex based on BAF values in chrX markers. The 0.1 threshold is used by default and it works well for Illumina but not Affy arrays. The `--bafxhet`argument can be used to change this default value. This value is calculated as fraction of chrX markers with BAF values between 0.25 and 0.75, so it is not the same (but smaller) as genotype heterozygosity rate from your genotype calling software.

    Right now chrY calling is not supported yet. If you want to call CNV for chrY, then you can remove all chrX markers from PFB file, then annotate chrY markers as located in chromosome X in PFB file, then use the -chrx argument instead.

1. **How to handle "weirld characters" in signal intensity files?**

    Sometimes, Illumina BeadStudio/Genome studio may export signal intensity values with weird characters. This could occur in non-English version of BeadStudio, in non-English version of Windows, in non-human SNP arrays, or any other reasons. For example, the LRR values for several markers in a file may display as "ABCDE" rather than a number, and PennCNV will ignore these values and ignore these markers in analysis. If they are "-inf" instead, PennCNV will treat them as -5. This should usually affect only a few markers for each sample.

    But some other times, all the signal intensity values are wrong so PennCNV will not work at all. For example, all the decimal points in LRR/BAF become "comma", so they are not valid numbers. In that case, users can do "perl -pe 's/,/./g' < inputfile > outputfile" to generate a new signal intensity file for CNV calling. One example is shown below:

    ```
[kaiwang@cc ~]$ head -n 3 sample.split1
Name    Chr     Position        4622780469_F.GType 4622780469_F.Log R Ratio        4622780469_F.B Allele Freq
rs109696 3       108367588       AB      -0,1223288      0,5079593
rs109701 13      39306521        BB      -0,1577153      1
rs109702 16      6508957 BB      -0,001872403    0,9723207

[kaiwang@cc ~]$ perl -pe 's/,/./g' < sample.split1 | head -n 3 
Name    Chr     Position        4622780469_F.GType 4622780469_F.Log R Ratio        4622780469_F.B Allele Freq
rs109696 3       108367588       AB      -0.1223288      0.5079593
rs109701 13      39306521        BB      -0.1577153      1
rs109702 16      6508957 BB      -0.001872403    0.9723207
    ```

1. **How to call CNVs if I have signal data for 50 SNPs in a candidate region?**

    Whenever whole-genome data is available, it is always best to use them for CNV calling even if the interest is on one particular gene or region. Sometimes, if you do not have access to whole-genome data, and you know the region or genes that you are interested in, it's still best to acquire data for the region plus flanking markers (for example, if the region per se has 50 markers, try to ask for 150 markers surrounding the region).

    If you want to find all CNVs (possibly with different sizes) in this region, use the --test argument for HMM-based CNV calling, as well as the --nomedianadjust argument. The latter argument is very important because by default a median adjustment procedure is used (to make the median LRR of whole-genome markers to be zero), and this procedure obviously should not apply here to candidate regions. If you know there exist a common CNV in this candidate region with known start and end marker and known deletion frequency and duplication frequency, you can also use the --validate argument for validation-based CNV calling to reduce false negative rates (but again the --nomedianadjust argument should be used for the operation).

1. **How to identify a subset of the most confident de novo CNV calls?**

    The 2009Aug27 vesion of PennCNV added a script for validating de novo CNVs and assigning P-values to de novo calls. If you want to know whether a particularly interesting de novo CNV is real or not, or if you want to select a set of most confident de novo CNVs for experimental validation, then this program should definitely be used. Check it out here.

1. **Can PennCNV give allele-specific CNV calls?**

    In practive, a user can multiply BAF with CN estimate to get allele-specific CNV calls. Believe it or not, it is as simple as that. People treat it as a big deal simply because most other software do not have the concept of BAF. Please read the web page on infer_snp_allele.pl for a more thorough description on this issue.

1. **How to interpret "state" in "state2,cn=1" in PennCNV output?**

    The word "state" is used internally in HMM algorithm, and its six-state definition is borrowed from the original quantisnp paper. You can simply ignore the "state", and focus on the CN=1 part, which means a deletion with copy number of 1.

1. **How to detect CNVs in Bovine genome?**

    First, compile your own PFB file (see format description here). Next, use the "--lastchr 29" argument in the latest version of PennCNV, which means that the last autosomal chromosome is 29 (rather than 22).

1. **Why I get "PennCNV compilation error: Can't locate loadable object" in Windows?**

    The PennCNV package includes pre-built executables for Windows, but only if you install 32-bit Perl (version 5.8.8 or 5.10.1). See the kext/ directory: it has two sub-directories 5.8.8 and 5.10.1, each containing appropriate DLL files for PennCNV. So if you install say 5.12.1, or 5.8.9, PennCNV won't work. You have to install the correct Perl version, or compile the executable yourself.

1. **Is there other similar software for CNV calling?**

    Some of the free software for CNV calling from SNP arrays include QuantiSNP, cnvPartition, BirdSuite, dChip, CNAT, CNAG, GenoCN and CORKEN. There are two recent free review articles on SNP arrays and CNVs here and here. There are quite a few companies that sell CNV calling software as well, including GoldenHelix, Partek, Nexus Copy Number software.

1. **How to call CNVs from tumor samples?**

    PennCNV considers normal copy to be 2 copies, and only gives integer estimate of copy number. For tumors, it is best to use an algorithm that specifically handles tumor samples, that gives continuous estimates of copy number and that accounts for aneuploidy levels as well as sample heterogeneity levels. This paper gives an overview of these issues and provides a software that only works in BeadStudio though. Other similar software include GAP, SOMATICs, GenoCN.

1. **What published study used PennCNV?**

    Some citations can be found [here](http://scholar.google.com/scholar?q=penncnv&hl=en&btnG=Search&as_sdt=800000000001).

1. **What databases are available for checking CNVs?**

    The most widely used CNV database is the Database for Genomic Variants (DGV). Other resources include the UCSC Genome Browser (with Structural Var track turned on), CNVVdb, CHOPPY, DECIPHER database and the new NCBI dbVar database.

1. **What tools can annotate CNV calls?**

    PennCNV package contains `scan_region.pl` program to annotate CNV calls (see tutorial here). Other free software or web server include WGAViewer, SCAN, PLINK CNV function .

