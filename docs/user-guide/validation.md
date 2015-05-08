The validation CNV calling algorithm in the PennCNV package is designed for calling CNVs in known common CNV regions, or CCNVRs. It does not use the HMM algorithm. Instead, it directly compares the posterior likelihood of copy number 0, 1, 2, 3 and 4 for a given genomic region, and selects the most likely copy number for this region. A “prior” probability of deletion and duplication should be provided for this region, which specifies the population prevalence of having deletion or duplication allele at this region.

This algorithm is somewhat similar to the “validation” procedure in trio-calling algorithm, which operates on a parent-offspring trio by validating the called CNV regions and selecting the most likely posterior copy number combination for a trio. Therefore, this algorithm is termed as validation-calling algorithm in PennCNV. Conceptually, it is also similar to a component called Canary in the BirdSuite program for the Affymetrix arrays. The BirdSuite package also provide a CCNVR file containing ~1000 genomic regions likely to harbor CNVs, and PennCNV validation-calling algorithm can operates on these ~1000 regions (for Affymetrix arrays). 
                                                            
An example of running the validation algorithm is given below:

```
[kai@adenine penncnv]$ detect_cnv.pl -validate -hmm lib/hh550.hmm -pfb lib/hh550.hg18.pfb -startsnp rs6894435 -endsnp rs6451161 -delfreq 0.15 -dupfreq 1e-6 -out sampleall.valicnv -conf sample1.txt sample2.txt sample3.txt
```

In the above command line, the -validation argument specifies that the validation-calling algorithm, as opposed to other HMM-based algorithm, should be used for CNV calling. The -startsnp and -endsnp combination tells the genomic region that should be examined. The -delfreq and -dupfreq arguments specify the prior probability of having a deletion allele and a duplication allele in this region. The result is the CNV call in this specified region, and we can see that both father and child has a deletion in this region with high confidence:

```
[kaiwang@cc penncnv]$ cat sampleall.valicnv 
chr5:17671545-17780897        numsnp=15     length=109,353     state2,cn=1 sample1.txt startsnp=rs6894435 endsnp=rs6451161 conf=63.042
chr5:17671545-17780897        numsnp=15     length=109,353     state2,cn=1 sample3.txt startsnp=rs6894435 endsnp=rs6451161 conf=55.186
```

When calling CNVs, the actual prior probability for copy numbers at this region is calculated by Hardy Weinberg equilibrium, based on the deletion/duplication allele frequency specified in the command line. For example, given a -delfreq of 0.01 the copy number 0 would have a prior probability of 0.01x0.01=0.0001, while copy number 1 would have a prior probability of 2x0.01x0.99=0.0198.

In most cases, rather than validating a specific region, we might want to validate all annotated “common CNV” in a given file, so-called candidate list file. This file can be supplied to the program by the -candlist argument. For example,

```
[kai@adenine penncnv]$ detect_cnv.pl -validate -hmm lib/hh550.hmm -pfb lib/hh550.hg18.pfb -candlist ccnvr -out sampleall.valicnv sample1.txt sample2.txt sample3.txt
```

The candidate list file contains five columns, representing chromosome region, start SNP, end SNP, deletion allele frequency and duplication allele frequency, respectively. All the regions specified in the file will be examined for CNV calls in the three files given in command line. For example,

```
[kaiwang@cc penncnv]$ cat ccnvr 
region  startsnp        endsnp  delfreq dupfreq
chr3:3974670-4071644    rs11716390      rs17039742      0.01    0
chr20:10440279-10511908 rs8114269       rs682562        0.01    0
```

The first line in the file is optional, since it will be ignored if starting with “region” in the line. PennCNV will score all candidate regions listed in this file.

In some cases, the user may want to know the exact log likelihood of all CNV state for a given region, rather than just the most likely copy number calls. This can be achieved by the -valilog argument. We can use the files in the "example/" directory in PennCNV package as an example below.

```
[kaiwang@cc ~/project/penncnv/example]$ detect_cnv.pl -validate -hmm example.hmm -pfb example.pfb -log ex4.log -out ex4.rawcnv -startsnp rs8114269 -endsnp rs682562 -delfreq 0.005 -list inputlist -valilog ex4.valilog
[kaiwang@cc ~/project/penncnv/example]$ cat ex4.valilog
rs8114269 rs682562 -57.329156975414 15.445268359525 -23.9674595603709 -37.735089228846 -41.4653035085859
rs8114269 rs682562 -60.7850764353335 -29.2158486346641 14.4433152191971 -24.4924596405965 -42.1807421936771
rs8114269 rs682562 -61.0796365132782 12.3599714086757 -29.2887422018293 -39.7798906904628 -44.3235225985156
```

Each of the numbers represent log likilihood of the CN=0, 1, 2, 3, and 4. As we can see from the above output, the three samples should have most likely copy number of 1, 2, 1, respectively (in other word, father and child has the deletion). Haveing these numbers (as opposed to the most likely CN state) may help the development of CNV association methods that can take into account of CNV calling uncertainty.

 