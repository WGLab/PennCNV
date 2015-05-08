Suppose we already generated the individual-based CNV calls for `sample1.txt`, `sample2.txt` and `sample3.txt`, as described in the previous tutorial section in CNV Calling. Below is a description of the procedure for trio-based (and quartet-based) CNV calling for each of the three individuals using PennCNV.
                            
The family structure can be used for generating more accurate CNV calls, since we can borrow and correlate CNV information from related family members that are very likely to share the same CNV region. To achieve this, we can run the following command:

```
[kai@adenine penncnv]$ detect_cnv.pl -trio -hmm lib/hh550.hmm -pfb lib/hh550.hg18.pfb -cnv sampleall.rawcnv sample1.txt sample2.txt sample3.txt -out sampleall.triocnv
```

In the above command, the `--trio` argument specify that we want to use family-based CNV detection algorithm to jointly update CNV status for a father-mother-offspring trio. The `--cnvfile` argument specify the prior CNV calls generated in individual-based calling step. The three files in command line represent signal data for father, mother and offspring, respectively. The output will be redirected and written to sampleall.cnv. Note that we can also generate a listfile, which contains 3 file names per line, to process multiple trios simultaneously.

The PennCNV trio-based calling algorithm analyzes the fifth column (file name column) of each line in the sampleall.rawcnv file, and then checks all individual-based CNV calls generated on any member of a trio (`sample1.txt`, `sample2.txt` and `sample3.txt`), and then try to re-call these regions on the trio and fine map boundaries. Therefore, it is important that the file names listed in the sampleall.rawcnv file is identical to the names in the command line, otherwise the program won’t work, since it cannot figure out the correct individual-based calls to use.

The first a few lines of the output file is listed below:

```
[kai@adenine penncnv]$ cat sampleall.triocnv
chr1:59077355-59078584        numsnp=3      length=1,230       state5,cn=3 sample1.txt startsnp=rs942123 endsnp=rs3015321 father triostate=533
chr1:147305744-147427061      numsnp=7      length=121,318     state5,cn=3 sample1.txt startsnp=rs11579261 endsnp=rs3853524 father triostate=535
chr1:147305744-147427061      numsnp=7      length=121,318     state5,cn=3 sample3.txt startsnp=rs11579261 endsnp=rs3853524 offspring triostate=535
chr1:153461604-153467859      numsnp=3      length=6,256       state5,cn=3 sample2.txt startsnp=rs2049805 endsnp=rs1045253 mother triostate=355
chr1:153461604-153467859      numsnp=3      length=6,256       state5,cn=3 sample3.txt startsnp=rs2049805 endsnp=rs1045253 offspring triostate=355
chr1:156783977-156788016      numsnp=6      length=4,040       state2,cn=1 sample1.txt startsnp=rs16840314 endsnp=rs10489835 father triostate=233
chr1:232415025-232419522      numsnp=4      length=4,498       state5,cn=3 sample1.txt startsnp=rs556585 endsnp=rs4333882 father triostate=535
chr1:232415025-232419522      numsnp=4      length=4,498       state5,cn=3 sample3.txt startsnp=rs556585 endsnp=rs4333882 offspring triostate=535
chr2:4191253-4200019          numsnp=3      length=8,767       state2,cn=1 sample1.txt startsnp=rs1175867 endsnp=rs1175854 father triostate=233
chr2:40075710-40100220        numsnp=9      length=24,511      state2,cn=1 sample2.txt startsnp=rs10865162 endsnp=rs2192721 mother triostate=323
chr2:183794033-183797494      numsnp=3      length=3,462       state2,cn=1 sample1.txt startsnp=rs17758247 endsnp=rs1462530 father triostate=232
chr2:183794033-183797494      numsnp=3      length=3,462       state2,cn=1 sample3.txt startsnp=rs17758247 endsnp=rs1462530 offspring triostate=232
chr2:208064035-208066083      numsnp=5      length=2,049       state2,cn=1 sample2.txt startsnp=rs918843 endsnp=rs959668 mother triostate=323
chr2:242565979-242656041      numsnp=16     length=90,063      state2,cn=1 sample1.txt startsnp=rs12987376 endsnp=rs6740738 father triostate=233
```

After using family information, we now generate a total of 62 CNVs for three members in this family.

Among the 21 CNVs in offspring (`sample3.txt`), 20 are inherited CNVs and 1 is de novo CNV. All inherited CNVs have identical boundaries as the corresponding CNVs from the father or mother. The de novo CNV is this one:

```
chr3:3974670-4071644          numsnp=50     length=96,975      state2  sample3.txt startsnp=rs11716390 endsnp=rs17039742 offspring triostate=332
```

As we can see, the new CNV file contains two extra fields: the eighth field indicates that sample3.txt is offspring in the trio-based CNV calling, while the ninth field tells us that the HMM states for the trio are 3 (normal), 3 (normal) and 2 (one-copy deletion) at this genomic region, respectively.

If the family has two children, then the `-quartet` argument can be used for CNV calling. Accordingly, four file names should be supplied in the command line, or given in each line of the list file, representing father, mother, child 1 and child 2, respectively.

PennCNV cannot generate calls on a pair of parents and 3 or more children; instead, the user need to split the family into trios and quartets for CNV calling, and then combine the CNV calls together into consensus calls.

