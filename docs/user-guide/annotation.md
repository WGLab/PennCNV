## Filtering CNV calls by user-specified criteria

The raw CNV calls often need to be filtered to keep a specific subset of calls for further analysis. In the PennCNV package, the filter_cnv.pl program can filter CNV calls based on various criteria, including both sample-level and on call-level criteria.

For example, `-numsnp 10 -length 50k` can be used to select CNV calls containing 10 SNPs and larger than 50kb.

```
[kaiwang@cc penncnv]$ filter_cnv.pl -numsnp 10 -length 50k sampleall.rawcnv 
chr2:242565979-242656041      numsnp=16     length=90,063      state2,cn=1 sample1.txt startsnp=rs12987376 endsnp=rs6740738
chr5:17671545-17780897        numsnp=15     length=109,353     state2,cn=1 sample1.txt startsnp=rs6894435 endsnp=rs6451161
chr11:55127597-55204003       numsnp=11     length=76,407      state2,cn=1 sample1.txt startsnp=rs2456022 endsnp=rs7934845
chr6:79034386-79088461        numsnp=22     length=54,076      state5,cn=3 sample1.txt startsnp=rs818258 endsnp=rs818280
chr3:3974670-4071644          numsnp=50     length=96,975      state2,cn=1 sample2.txt startsnp=rs11716390 endsnp=rs17039742
chr5:17671545-17780897        numsnp=15     length=109,353     state2,cn=1 sample2.txt startsnp=rs6894435 endsnp=rs6451161
chr20:10440279-10511908       numsnp=10     length=71,630      state2,cn=1 sample2.txt startsnp=rs8114269 endsnp=rs682562
```

As we can see, in the output, only CNV calls meeting specified criteria are printed out to the screen. Adding the `-output` argument to the above command will write the CNV calls to a file.

## Merging adjacent CNV calls

For recently developed SNP arrays with high-density markers, PennCNV may tends to split large CNVs (such as those >500kb) into smaller parts (such as two or three 150kb CNV calls). It is possible to examine the CNV calls and merge adjacent calls together into one single call, using user-specified threshold (for example, gap\<20% of total length). In other word, suppose there are three genomic segments, A, B and C, whereas A and C are called as deletion by PennCNV. If you divide the length of the gap B (measured by #probes) by the length of A+B+C, and if this fraction if <=20%, then it is recommended to merge A+B+C as a single deletion call.

In the latest version of PennCNV, the `clean_cnv.pl` program can do this automatically for you, with the -fraction argument to control the fraction threshold. The -bp argument can be used so that the fraction is calculated as base pair length, rather than number of SNPs/markers withint he call.

## Automatic Quality Control of CNV calls

The `filter_cnv.pl` program also has functionality to perform sample-level QC, in addition to call-level QC. Therefore, we can use the program to identify low-quality samples from a genotyping experiment, and eliminate them from future analysis. This analysis requires the LOG file used in CNV calling.

The three samples (`sample1.txt`, `sample2.txt` and `sample3.txt`) used in the tutorial are actually of very good data quality. However, to illustrate how to apply the QC procedure, next we specify an extraordinarily stringent LRR_SD criteria to show that one sample (`sample3.txt`) cannot meet this criteria and is dropped during the QC procedure.

```
[kaiwang@cc penncnv]$ filter_cnv.pl sampleall.rawcnv -qclogfile sampleall.log -qclrrsd 0.135 -qcpassout sampleall.qcpass -qcsumout sampleall.qcsum -out sampleall.goodcnv
NOTICE: the --qcbafdrift argument is set as 0.01 by default
NOTICE: the --qcwf argument is set as 0.05 by default
NOTICE: Writting 2 file names that pass QC to qcpass file sampleall.qcpass
NOTICE: Writting 3 records of QC summary to qcsum file sampleall.qcsum
```

The above command asked to analyze the LOG file (`sampleall.log`), find all samples with LRR_SD less than 0.135, then write these samples to the `sampleall.qcpass` file, write the CNV calls for these samples to the sampleall.goodcnv file, and write the QC summary for all samples to the sampleall.qcsum file.

Now we can examine the `sampleall.qcsum` file:

```
[kaiwang@cc penncnv]$ cat sampleall.qcsum 
File    LRR_mean        LRR_median      LRR_SD  BAF_mean        BAF_median      BAF_SD  BAF_drift       WF      NumCNV
sample3.txt     0.0022  0.0000  0.1260  0.5045  0.5000  0.0431  0.000186        -0.0177 17
sample2.txt     0.0034  0.0000  0.1338  0.5063  0.5000  0.0391  0.000070        0.0194  18
sample1.txt     0.0038  0.0000  0.1383  0.5046  0.5000  0.0425  0.000151        0.0120  16
```

This is a tab-delimited text file that can be easily loaded into Excel for plots and histograms. For example, it is often informative to plot the number of CNV calls versus the LRR_SD measure to diagnose what's a good LRR_SD threshold to use in QC for a particular data set.

Now examine the sampleall.qcpass file:

```
[kaiwang@cc ~/project/penncnv/tutorial]$ cat sampleall.qcpass 
sample3.txt
sample2.txt
```

This means that sample1.txt does not pass the QC threshold. Checking the qcsum file, we can see that the sample1.txt has LRR_SD of 0.138, which is higher than the 0.135 as specified in the command line.

Now check the sampleall.goodcnv file: we will see that only sample2.txt and sample3.txt are included in this file, that is, the CNV calls from “bad” sample is filtered from the output file.

Again, the above example is used only for illustrating how the procedure works. In reality the 0.135 threshold is too stringent. Using a 0.24, or even 0.3 or 0.35, would be advised. (it does not hurt to change the threshold and see how many samples pass QC). Finally, it is highly recommended to also set the -qcnumcnv argument in the command line in practical settings; for example, “-qcnumcnv 50” would treat any samples with >50 CNV calls as low quality samples and eliminate them from analysis. (taking the qcsum file and draw a histogram to see the distribution and decide on a good threshold to use; it is not unusual for PennCNV to generate 2,000 CNV calls for a sample with very low quality even if LRR_SD appear to be normal!).

## Removing spurious CNV calls in specific genomic regions

Several genomic regions are known to harbor spurious CNV calls that represent cell-line artifacts, so the calls should be eliminated before analysis. For example, the original PennCNV paper in Genome Research has done a comparison to show that immunoglobulin regions are especially likely to have deletions in cell line samples than whole-blood samples. Additionally, centromeric and telomeric regions are especially likely to harbor spurious CNV calls as well. The scan_region.pl program can easily help remove CNV calls in specific genomic regions.

Some examples were given in the FAQ section since many people asked about this issue. One thing to note is that people may use different thresholds (like 100kb, 500kb, 1Mb) to define telomeric/centromeric regions. Another thing to note is that CNVs in some telomeric regions may be indeed functionally important rather than simple artifacts. So always double check when deleting things from the CNV call files to make sure that they make sense.

## Finding overlapping/neighboring genes for CNV calls

One of the most common tasks for CNV annotation is to identify overlapping or neighboring genes. If the CNV calls are generated using hg18 (Mar 2006, NCBI build 36) human genome assembly, we can download UCSC known gene annotation (`knownGene.txt.gz` and `kgXref.txt.gz`) or refGene annotation (`refGene.txt.gz` and `refLink.txt.gz`) for CNV annotation. (After downloading these files, it is recommended to first unzip the files (for example, gunzip refGene.txt.gz), then rename them (for example, `mv refGene.txt hg18_refGene.txt`) to add the genome build information to the file name to remind yourself about the version of genome assembly).

For example, to identify UCSC known genes that overlap with CNV calls (which were generated using the hg18 genome coordinates), we can run:

```
[kai@adenine penncnv]$ scan_region.pl sampleall.cnv hg18_refGene.txt -refgene -reflink hg18_refLink.txt > sampleall.cnv.rg18
```

The output file contains two additional columns to each line of the sampleall.cnv file. The first column represents the gene symbols and the second column indicates the distance between CNV and gene. In our case, the distance is always zero since we will only find CNVs that overlap with a gene. If the CNV does not overlap with any gene, the “NOT_FOUND” notation will be shown for the corresponding CNVs.

For example, several lines of the `sampleall.cnv.rg18` file is shown below:

```
chr2:242565979-242656041      numsnp=16     length=90,063      state2,cn=1 sample1.txt startsnp=rs12987376 endsnp=rs6740738 father triostate=233        NOT_FOUND       NOT_FOUND
chr3:3974670-4071644          numsnp=50     length=96,975      state2,cn=1 sample3.txt startsnp=rs11716390 endsnp=rs17039742 offspring triostate=332    NOT_FOUND       NOT_FOUND
chr3:37957465-37961253        numsnp=3      length=3,789       state2,cn=1 sample2.txt startsnp=rs9837352 endsnp=rs9844203 mother triostate=323 CTDSPL  0
chr3:75511365-75650909        numsnp=7      length=139,545     state2,cn=1 sample2.txt startsnp=rs4677005 endsnp=rs2004089 mother triostate=323 NOT_FOUND       NOT_FOUND
chr4:25166145-25184635        numsnp=8      length=18,491      state5,cn=3 sample2.txt startsnp=rs7686265 endsnp=rs7689639 mother boundary_reconciled=355       NOT_FOUND       NOT_FOUND
```

For the above five CNV calls, only one overlaps with a gene (CTDSPL), while others are in intergenic regions. The second CNV call (50 SNPs, 97kb) is a de novo CNV.

It is usually more useful to find neighboring genes for an intergenic CNV, we can use the `--expandmax` argument:

```
[kai@adenine penncnv]$ scan_region.pl sampleall.cnv hg18_refGene.txt -refgene -reflink hg18_refLink.txt -expandmax 5m > sampleall.cnv.rg18
```

This will expand the CNV up to 5 megabases in both direction and then try to find neighboring genes. Only the closest gene to the CNV will be written to output, while this closest gene might be located to the left or right side of the CNV (note that we use “left” and “right” here since CNVs occur in both forward and reverse chains without a predefined direction). To find only left genes, we can use the `--expandleft 5m` argument. To find both left and right genes, we have to run the program twice with `--expandleft` and `--expandright` argument respectively. This can be done easily in one single step:

```
[kai@adenine penncnv]$ scan_region.pl sampleall.cnv hg18_refGene.txt -refgene -reflink hg18_refLink.txt -expandleft 5m | scan_region.pl sampleall.cnv hg18_refGene.txt -refgene -reflink hg18_refLink.txt -expandright 5m > sampleall.cnv.rg18
```

Note that when <inputfile> is “stdin”, the scan_region.pl program will read input from STDIN (standard input, which can be a piped output from a previous command). The output will contain four extra columns, representing closest left gene, left distance, closest right gene and right distance, respectively. For example, the same five CNVs mentioned above were annotated as:

```
chr2:208064035-208066083      numsnp=5      length=2,049       state2,cn=1 sample2.txt startsnp=rs918843 endsnp=rs959668 mother triostate=323   KLF7    325176  CREB1   36848
chr2:242565979-242656041      numsnp=16     length=90,063      state2,cn=1 sample1.txt startsnp=rs12987376 endsnp=rs6740738 father triostate=233        FLJ33590        101824  NOT_FOUND       NOT_FOUND
chr3:3974670-4071644          numsnp=50     length=96,975      state2,cn=1 sample3.txt startsnp=rs11716390 endsnp=rs17039742 offspring triostate=332    LRRN1   110283  SETMAR  248372
chr3:37957465-37961253        numsnp=3      length=3,789       state2,cn=1 sample2.txt startsnp=rs9837352 endsnp=rs9844203 mother triostate=323 CTDSPL  0       CTDSPL  0
chr3:75511365-75650909        numsnp=7      length=139,545     state2,cn=1 sample2.txt startsnp=rs4677005 endsnp=rs2004089 mother triostate=323 CNTN3   858332  ROBO2   1521075
chr4:25166145-25184635        numsnp=8      length=18,491      state5,cn=3 sample2.txt startsnp=rs7686265 endsnp=rs7689639 mother boundary_reconciled=355       ANAPC4  136928  SLC34A2 81898
```

Tips: the UCSC known gene annotation is more comprehensive than RefSeq annotation, but it contains too many uncharacterized transcripts (with names such as AK091772 or BC015880). In most cases, it is probably a better idea to use RefSeq annotation to annotate CNVs.

## Finding overlapping exons for CNV calls

We often want to know which CNVs severely affects genes, so examination of exon deletions may be one way to ensure that the CNV call indeed disrupt gene function. The `-refexon` argument, rather than `-refgene` argument, can be used in the above command line to find out exonic overlaps.

Those CNV calls without exonic overlap with have “NOT_FOUND” appended to the end of the line, so a `fgrep -v NOT_FOUND` can be added to get rid of these CNVs affecting non-exonic regions.

## Finding overlapping functional elements for CNV calls

The `scan_region.pl` program can perform these tasks using the UCSC genome annotation file. Please use the “-m” argument for the scan_region.pl program to read a complete manual on how to achieve these goals.

## CNV case-control comparison

The PennCNV program implements a very preliminary functionality of case-control association analysis, to identify a stretch of SNPs that tend to have copy number changes in cases versus controls using Fisher’s Exact Test. This function is very preliminary and very rough, but can be a first-pass effort to identify potentially interesting regions in a case-control setting.

A more formal way to do case-control association analysis can be accomplished by ParseCNV, which is availalbe [here](http://parsecnv.sourceforge.net/). ParseCNV takes CNV calls as input and creates probe based statistics for CNV occurrence in (cases and controls, families, or population with quantitative trait) then calls CNVRs based on neighboring SNPs of similar significance. CNV calls may be from aCGH, SNP array, Exome Sequencing, or Whole Genome Sequencing.

A phenotype file needs to be supplied for this analysis. The file contains two tab-delimited columns, representing file name and the disease label, respectively. The disease label can be 1 (indicating non-affected) or 2 (indicating affected), or can be “control” and “case”. In fact, The disease label can be anything: if using the “-control_label” argument, then any string that match the control label will be treated as controls, otherwise treated as cases. The following is an example phenofile:

```
[kaiwang@cc penncnv]$ cat phenotype 
sample1.txt     chronic_disease
sampel2.txt     control
sample3.txt     acute_disease
```

So the program will treat `sample2.txt` as control, and the other samples as cases. If a file name is not annotated in the phenofile, it will be treated as unknown and not used in the association analysis.

For illustration purposes, using the CNV call file generated in the tutorial, now we can perform a simple case-control comparison between the case group (two samples) and the control group (one sample), using the CNV call file `sampleall.rawcnv`:

```
[kaiwang@cc ~/project/penncnv/tutorial]$ detect_cnv.pl -cctest -pfb ../lib/hh550.hg18.pfb -phenofile phenotype -cnv sampleall.rawcnv
NOTICE: Finished reading CNV file sampleall.rawcnv (1 individuals are skipped due to lack of phenotype annotation)
NOTICE: Performing case-control comparison by Fisher's Exact Test (8 output columns are: SNP, Chr, Position, case-cnv, case-normal, control-cnv, control-normal, P-value)
rs942123        1       59077355        1       1       0       1       1
rs942124        1       59077498        1       1       0       1       1
rs3015321       1       59078584        1       1       0       1       1
rs11579261      1       147305744       2       0       0       1       0.333333333333333
rs17162082      1       147306690       2       0       0       1       0.333333333333333
```

The association results will be printed in the screen in tab-delimited text format (using -out argument to write to a file). For example, for the SNP rs11579261, two cases (out of two total cases) have CNV, but zero control (out of one control) has CNV. The P-value is 0.33 by Fisher’s Exact Test. Note that the results are printed in a per-marker basis: obviously if a stretch of neighboring markers all have good P-values it may indicates that this entire CNV region is associated with disease status.

Some additional useful arguments include `-type del` or `-type dup` to examine deletions or duplications only. (By default both deletions and duplications are counted as CNV in the above calculation). The `-onesided` argument can be used to performed one-sided association test that only cares about regions where cases have more CNVs than controls.
