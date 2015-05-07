Starting from the November 2008 version of PennCNV, several example data sets and scripts have been included to test that the programs are installed correctly, and to demonstrate some quick examples to use the various programs in PennCNV. In the penncnv/example/ directory, we will see several files there. Among them, the father.txt, mother.txt and offspring.txt are three signal intensive files with signal values (to keep the file size small, only a few chromosomes are included in these files). In addition, there are several list files that contain file names to be processed by PennCNV.
                                
If we run the program ./runex.pl in this directory, we will see a description of available options:

[kaiwang@cc ~/project/penncnv/example]$ ./runex.pl 
Usage:
     runex.pl <1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16>

     Optional arguments:
            -v, --verbose                   use verbose output
            -h, --help                      print help message
            -m, --man                       print complete documentation
                --path_detect_cnv <string>          path to detect_cnv.pl
                --path_visualize_cnv <string>       path to visualize_cnv.pl
                --path_convert_cnv <string>         path to convert_cnv.pl
                --path_filter_cnv <string>          path to filter_cnv.pl
                --path_compare_cnv <string>         path to compare_cnv.pl

     Function: test-drive PennCNV and related scripts in your system

     Example: runex.pl 1 (run PennCNV to call CNV on three signal files)
              runex.pl 2 (run posterior CNV calling algorithm on a trio)
              runex.pl 3 (run PennCNV with GCmodel adjustment of signal intensity)
              runex.pl 4 (validation-based CNV calling on a candidate region)
              runex.pl 5 (validation-based CNV calling on all candidate regions in a file)
              runex.pl 6 (run joint CNV calling on a trio)
              runex.pl 7 (convert CNV call to BED format)
              runex.pl 8 (convert CNV call to tab-delimited format)
              runex.pl 9 (convert tab-delimited calls to PennCNV format)
              runex.pl 10 (filter CNV and print out a subset of calls)
              runex.pl 11 (compare CNV calls on duplicated samples)
              runex.pl 12 (compare CNV calls on same sample called by different algorithms)
              runex.pl 13 (generate CNV-based genotype calls)
              runex.pl 14 (validate de novo CNV calls and assign a P-value)
              runex.pl 15 (convert Canary CNV calls to PennCNV format)
              runex.pl 16 (plot signal intensity for CNV calls in JPG formats)

 

The user can try to run these examples one by one and get some idea on what PennCNV can do and how to use the command line options. For example, let’s first try the first example:

[kaiwang@cc ~/project/penncnv/example]$ runex.pl 1
Exercise 1: individual-based calling and write the output to ex1.rawcnv
        Running command <detect_cnv.pl -test -hmm ../lib/hh550.hmm -pfb ../lib/hh550.hg18.pfb -conf -log ex1.log -out ex1.rawcnv father.txt mother.txt offspring.txt>

******************************************************************************
NOTICE: All program notification/warning messages that appear in STDERR will be also written to log file ex1.log
NOTICE: Reading marker coordinates and population frequency of B allele (PFB) from ../lib/hh550.hg18.pfb ... Done with 566108 records (178 records in chr M,XY were discarded)
NOTICE: Reading LRR and BAF values for from father.txt ... Done with 93129 records in 4 chromosomes
NOTICE: Data from chromosome X will not be used in analysis
NOTICE: Median-adjusting LRR values for all autosome markers from father.txt by 0.0221
NOTICE: Median-adjusting BAF values for all autosome markers from father.txt by 0.0029
NOTICE: quality summary for father.txt: LRR_mean=0.0027 LRR_median=0.0000 LRR_SD=0.1335 BAF_mean=0.5063 BAF_median=0.5000 BAF_SD=0.0390 BAF_DRIFT=0.000037 WF=0.0184 GCWF=0.0136
NOTICE: Reading LRR and BAF values for from mother.txt ... Done with 93129 records in 4 chromosomes
NOTICE: Data from chromosome X will not be used in analysis
NOTICE: Median-adjusting LRR values for all autosome markers from mother.txt by -0.0199
NOTICE: Median-adjusting BAF values for all autosome markers from mother.txt by 0.0324
NOTICE: quality summary for mother.txt: LRR_mean=0.0039 LRR_median=0.0000 LRR_SD=0.1374 BAF_mean=0.5044 BAF_median=0.5000 BAF_SD=0.0418 BAF_DRIFT=0.000140 WF=0.0100 GCWF=0.0028
NOTICE: Reading LRR and BAF values for from offspring.txt ... Done with 93129 records in 4 chromosomes
NOTICE: Data from chromosome X will not be used in analysis
NOTICE: Median-adjusting LRR values for all autosome markers from offspring.txt by -0.0087
NOTICE: Median-adjusting BAF values for all autosome markers from offspring.txt by 0.0260
NOTICE: quality summary for offspring.txt: LRR_mean=0.0028 LRR_median=0.0000 LRR_SD=0.1263 BAF_mean=0.5045 BAF_median=0.5000 BAF_SD=0.0429 BAF_DRIFT=0.000293 WF=-0.0171 GCWF=-0.0100

******************************************************************************

Note: when running the examples, the program suppose that the PennCNV executables are already in your system path. If not, then you can use"runex.pl 1 -path_detect_cnv ../detect_cnv.pl" instead to specify the path to the executable (the executable is located at the upper directory). 

To add the penncnv directory to PATH environment variable: In Unix/Linux system with BASH shell, suppose the penncnv is installed in the /home/user1/penncnv directory, then you can either execute "export PATH=$PATH:/home/user1/penncnv" to add the directory to the system PATH environment variable, or directly add the above command to the ~/.bashrc file. In Windows system, go to "Control Panel", click "System", click "Advanced" tab, click "Environment variables" button, then select "Path" and click "Edit" button, manually type in the path and append it to the current path environment variables.

In the above command, the program gives the first example, which use individual-based CNV calling algorithm on three input files (father.txt, mother.txt and offspring.txt), and then write the output to the ex1.rawcnv file. The actual command line arguments are printed after “Running command”. Some LOG information was printed out in the screen between the two ******* lines, but they are also written to the ex1.log file as well. We can check the content of the ex1.rawcnv file:

[kaiwang@cc ~/project/penncnv/example]$ cat ex1.rawcnv
chr3:37957465-37961253        numsnp=3      length=3,789       state2,cn=1 father.txt startsnp=rs9837352 endsnp=rs9844203 conf=15.133
chr3:75511365-75650909        numsnp=7      length=139,545     state2,cn=1 father.txt startsnp=rs4677005 endsnp=rs2004089 conf=26.862
chr11:81181640-81194909       numsnp=9      length=13,270      state2,cn=1 father.txt startsnp=rs7947005 endsnp=rs12293984 conf=35.083
chr20:10440279-10511908       numsnp=10     length=71,630      state2,cn=1 father.txt startsnp=rs8114269 endsnp=rs682562 conf=39.410
chr11:55127597-55204003       numsnp=11     length=76,407      state2,cn=1 mother.txt startsnp=rs2456022 endsnp=rs7934845 conf=23.106
chr11:539119-548884           numsnp=4      length=9,766       state5,cn=3 mother.txt startsnp=rs4963136 endsnp=rs2061586 conf=7.342
chr11:55127597-55193702       numsnp=8      length=66,106      state1,cn=0 offspring.txt startsnp=rs2456022 endsnp=rs17498926 conf=48.159
chr3:3974670-4071644          numsnp=50     length=96,975      state2,cn=1 offspring.txt startsnp=rs11716390 endsnp=rs17039742 conf=219.722
chr11:81181640-81194909       numsnp=9      length=13,270      state2,cn=1 offspring.txt startsnp=rs7947005 endsnp=rs12293984 conf=41.264
chr20:10440279-10511908       numsnp=10     length=71,630      state2,cn=1 offspring.txt startsnp=rs8114269 endsnp=rs682562 conf=41.645

These fields are chromosome coordinates, number of markers (SNPs markers and sometimes CN markers as well) in the region, the CNV length, the copy number estimate, the signal file name, the start and end SNP and the confidence score. See more in-depth description on the default CNV calling algorithm here.

The example 2 illustrates trio-based calling algorithm, which requires the output file from example 1 (ex1.rawcnv) as one of the input files. The example 6 illustrates a joint-calling algorithm for trios that uses one step only, and does not require the ex1.rawcnv as input files. See more in-depth description on the trio-based algorithm here.

The example 3 illustrates the use of GC-model adjustement of signal intensity values for CNV calling. The algorithm was previously published and described here. The signal values for the 3 files are not really affected by genomic waves, so the adjustment has little effects on CNV calls.

The example 4 and 5 illustrate the validation-calling algorithm. This is not based on HMM, but instead based on a validation subroutine that takes prior probability, calculate likelihood of the region being various copy numbers, and then select the most likely copy number. See more in-depth description on the validation-calling algorithm here.

The example 6 illustrates the joint-calling algorithm for trios. This is a HMM-based algorithm that simultaneously model the signal measures for a trio (the detailed algorithm is described in this paper), so it is computationally expensive and may take quite some time to run.

The example 7, 8 and 9 illustrate the use of convert_cnv.pl program to convert CNV call formats.

The example 10 illustrates the use of filter_cnv.pl program to select a subset of CNV calls.

The example 11 and 12 illustrate the use of compare_cnv.pl program to compare CNV calls on same file given by different algorithms, or on duplicated samples given by the same algorithm.

The example 13 illustrates the use of infer_snp_allele.pl program to infer CNV-based genotype calls in CNV regions for three subjects. See more in-depth description on the program here.

The example 14 illustrates the use of infer_snp_allele.pl program to validate putative de novo CNV calls and assign P-values. See more in-depth description on the program here.

The example 15 illustrates the use of convert_cnv.pl program to convert CNV calls from other algorithms to penncnv output format. The resulting calls are useful for comparative analysis of calls between algorithms, or can be used in visualize_cnv.pl program to plot the actual signal intensity values.

The example 16 illustrates the use of visualize_cnv.pl program to plot signal intensity values (LRR and BAF) for all CNV calls for a given individual, so that users can visually examine and decide whether the calls are reliable or not (without relying on manual examination of GenomeStudio). The program requires calling R subroutines to work.

Note that the output_expected/ directory under the example/ directory contains the expected output files. If you see a difference between your output and the expected output, then maybe there is a problem with the PennCNV installation and a re-compilation is necessary.

