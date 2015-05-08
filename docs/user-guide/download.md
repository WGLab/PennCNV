## PennCNV main package

The link to latest stable version of PennCNV main package is given below. The package includes both source codes and pre-compiled executables for Windows. (Due to the problem of several unzipping software under Windows in handling *.tar.gz files, Windows users should try to download the zip file instead. The contents of these two files are identical).

- [tar.gz file](https://github.com/WangGenomicsLab/PennCNV/archive/v1.0.0.tar.gz)
- [zip file](https://github.com/WangGenomicsLab/PennCNV/archive/v1.0.0.zip)

Alternatively, users can download the cutting edge version from GitHub by `git clone git@github.com:WangGenomicsLab/PennCNV.git`.

## User-supplied files and programs

Users can use the `compile_pfb.pl` program included in PennCNV to generate a PFB file for their specific arrays. A few additional PFB files are provided below:

- The [hc12v1.hg18.pfb.gz](http://www.openbioinformatics.org/penncnv/download/hc12v1.hg18.pfb.gz) and [hc12v1.hg18.gcmodel.gz](http://www.openbioinformatics.org/penncnv/download/hc12v1.hg18.gcmodel.gz) file are provided for HumanCytoSNP12 V1 array with 300K markers on that array specifically (based on Caucasian populations).

- The [ho1v1.hg18.pfb.gz](http://www.openbioinformatics.org/penncnv/download/ho1v1.hg18.pfb.gz) and [ho1v1.hg18.gcmodel.gz](http://www.openbioinformatics.org/penncnv/download/ho1v1.hg18.gcmodel.gz) files are provided for HumanOmni1 QuadV1 array specifically (based on Caucasian populations).

- **Updated 2011Apr25**: A user (Denise Kay) submitted an updated [ho1v1.hg18.pfb.gz](http://www.openbioinformatics.org/penncnv/download/ho1v1.hg18.MOD.zip) file in which 32,290 markers have modified PFB=2; these are originally designed as SNP markers by Illumina but their quality is too low so it is best to treat them as intensity-only markers by flagging PFB=2.

- **Updated 2011Mar18**: Dr. Luis M. Franco from Baylor College of Medicine provided [PFB file](http://www.openbioinformatics.org/penncnv/download/OmniExpress_hg18.pfb.gz) for the Illumina OmniExpress array compiled from 96 unrelated Caucasian individuals. The number of SNPs in the files current form is 703,965 ( pseudoautosomal and Y-chromosome SNPs are excluded).

- **Updated 2011Aug16**: Dr. Colm O'Dushlaine lifted over the PFB file for Illumina OmniExpress array to [hg19 coordinate](http://www.openbioinformatics.org/penncnv/download/OmniExpress_hg19.pfb.gz).

- **Updated 2013Nov06**: Dr. Lijian Yu from Harvard University provided Affymetrix PFB files in [hg19 coordinate](http://www.openbioinformatics.org/penncnv/download/affygw6.hg19.pfb.gz), using version na33 annotations from Affymetrix. Fie updated on 2014Jul23 to fix white space problem in the file.

- **Updated 2014Sep13**: Stephen Sanders from UCSF provided [PFB file](http://www.openbioinformatics.org/penncnv/download/YALE_Merged_PFB_hg19.pfb) for Illumina Omni 5 array, generated from 600 controls from the Simon Simplex Collection.

A few user contributed programs/scripts are provided below:

- **User contributed programs (Uploaded 2010Dec01)**: The [penncnv_to_plink.pl](http://www.openbioinformatics.org/penncnv/download/penncnv_to_plink.pl) program can be used to convert PennCNV output into PLINK input format. Run it with -h argument to read the documentation. The program is written and provided by Matthew Gillman at the Wellcome Trust Sanger Institute.

- **User contributed programs (Uploaded 2011Feb27)**: The [plot_raw_PennCNV.R](http://www.openbioinformatics.org/penncnv/download/plot_raw_PennCNV.R) program can plot from PennCNV rawcnv file on screen or to a high solution png file. Example screen shot [1](img/PennCNV_76neu_3unique_chr18.png), [2](img/PennCNV_76neu_3unique_chr22.png), [3](img/test_chr02.png) and [4](img/test_chr18.png). This program is written band provided by Dr. Bowang Chen. See updated script below.

- **User contributed programs (Uploaded 2012Jul16)**: The [plot_raw_PennCNV2PDF.R](http://www.openbioinformatics.org/penncnv/download/plot_raw_PennCNV2PDF.R) program can plot from PennCNV rawcnv file to a high solution PDF file. This new script exports all plots to a pdf file, by default 70 samples/page (can be changed by users). It is much fast now, it plots a 50000+ line raw file to a 400+ page pdf in ~2 minutes. Example output is here. Example command line is "R CMD BATCH --no-save --no-restore "--args filename=\"$filename\" PDFname=\"$PDFname\"" $path_R/plot_raw_PennCNV2PDF.R ". This program is written band provided by Dr. Bowang Chen@German Cancer Research Center.

## PennCNV example data sets

This example data set contains genotyping data for a father-mother-offspring trio, genotyped on the Illumina HumanHap550 array. The file is in BeadStudio project format, and can be opened by the Illumina BeadStudio software. The file is used in the tutorial in the web site.

- [tutorial_beadstudio.zip](http://www.openbioinformatics.org/penncnv/download/tutorial_beadstudio.zip)

This file is an exported text file from the Illumina BeadStudio software using the above project file. It is useful for users without access to BeadStudio to test PennCNV.

- [tutorial_text.zip](http://www.openbioinformatics.org/penncnv/download/tutorial_text.zip)

This serial dilution data set contains genotyping data on a sample genotyped five times, each time with ~2-fold dilutions. For more details please see the Diskin et al paper. The file is in BeadStudio project format, and can be opened by the Illumina BeadStudio software.

- [tutorial2.zip](http://www.openbioinformatics.org/penncnv/download/tutorial2.zip)
 

## PennCNV-Affy protocols

This package contains PennCNV-Affy protocols and helper scripts for handling Affymetrix Mapping 500K, genome-wide 5.0 and genome-wide 6.0 data as raw CEL files. For instructions on how to use the PennCNV-Affy package, please refer to the PennCNV-Affy link in the menu bar to the left of the page. The package below is updated on 2009Aug27, with very minor bug fixes only. It should produce identical results as previous version.

- [gw6.tar.gz](http://www.openbioinformatics.org/penncnv/download/gw6.tar.gz)
 

## Genome Browser tracks for SNP arrays

Click any of the link below will load the corresponding track in the UCSC Genome Browser on hg18 coordinate, or 2006 human genome assembly. It may take several minutes to load each of the track. The user can visualize the location of the marker coverage in each specific arrays and compare them with each other. If the user does not want to see the marker names in the browser, the dense drop-down menu option can be used in the browser page to display red vertical bars for each marker.

All the genome coordinates are also annotated in the PFB file supplied with the PennCNV package. The Illumina data were based on supplied annotation from Illumina's sample files. The Affymetrix data were based on na26 (July 2008) annotation from Affymetrix website.

- [Illumina HumanHap300](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/hh300.hg18.bed.gz)
- [Illumina HumanHap550v1](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/hh550v1.hg18.bed.gz)
- [Illumina HumanHap550v3](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/hh550v3.hg18.bed.gz)
- [Illumina Human660-Quad v1](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/hh660v1.hg18.bed.gz)
- [Illumina HumanHap650Y](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/hh650.hg18.bed.gz)
- [Illumina HumanHap1M](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/hh1m.hg18.bed.gz)
- [Illumina HumanCNV370](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/hc370.hg18.bed.gz)
- [Illumina Human-1](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/human1.hg18.bed.gz)
- [Illumina Human-610](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/human610.hg18.bed.gz)
- [Illumina IBC v1](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/ibcv1.hg18.bed.gz)
- [Illumina OmniExpress](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/omniexpress.bed.gz)
- [Affymetrix Mapping 500K](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/affy500k.hg18.bed.gz)
- [Affymetrix Genomewide 5.0](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/affygw5.hg18.bed.gz)
- [Affymetrix Genomewide 6.0](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hgt.customText=http://www.openbioinformatics.org/penncnv/gbtrack/affygw6.hg18.bed.gz)

Please report annotation mistakes/bugs or request for additional arrays.

