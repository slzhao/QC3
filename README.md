Table of Content
================
* [Introduction](#Introduction)
* [Download](#download)
* [Change log](#Change)
 * [Release Version 1.33](#R133)
 * [Release Version 1.32](#R132)
 * [Release Version 1.31](#R131)
 * [Release Version 1.30](#R130)
 * [Release Version 1.26](#R126)
 * [Release Version 1.25](#R125)
 * [Release Version 1.24](#R124)
 * [Release Version 1.23](#R123) 
 * [Release Version 1.22](#R122)
 * [Release Version 1.21](#R121)
 * [Release Version 1.20](#R120)
 * [Release Version 1.15](#R115)
 * [Release Version 1.10](#R110)
 * [Release Version 1.00](#R100)
 * [Release candidate (RC) version 1.3](#RC13)
 * [Release candidate (RC) version 1.2](#RC12)
 * [Release candidate (RC) version 1.1](#RC11)
 * [Release candidate (RC) version 1.0](#RC10)
* [Prerequisites](#Prerequisites)
 * [Install required perl packages](#irpp)
 * [Install required software](#irs)
 * [Download required database](#drd)
* [Usage](#Usage)
 * [input file](#if)
 * [fastq QC](#fqc)
 * [bam QC](#bqc)
 * [vcf QC](#vqc)
* [A simple example](#sExample)
 * [Download example data](#ded)
 * [fastq QC example](#fqcE)
 * [bam QC example](#bqcE)
 * [vcf QC example](#vqcE)
* [An example from TCGA data](#TExample)
 * [fastq QC example](#TfqcE)
 * [bam QC example](#TbqcE)
 * [vcf QC example](#TvqcE)
* [Results](#Results)
 * [Introduction for fastq QC result](#ifr)
 * [Introduction for bam QC result](#ibr)
 * [Introduction for vcf QC result](#ivr)

<a name="Introduction"/>
# Introduction #

High throughput sequencing is the most effective way to screen for non-specific germline variants, somatic mutations, and structural variants. Some of the most popular sequencing paradigms in DNA sequencing are whole genome sequencing, exome sequencing, and target panel sequencing. While vastly informative, sequencing data poses significant bioinformatics challenges in areas such as data storage, computation time, and variant detection accuracy.  One of the easily overlooked challenges associated with sequencing is quality control Quality control (QC) for DNA sequencing can be categorized based on three stages: raw data, alignment, and variant calling. QC on raw sequencing data has been given more attention than QC on alignment and the variant calling. There are many QC tools aimed at raw data such as FastQC, FastQ Screen FastX-Toolkit, NGS QC Toolkit, RRINSEQ and QC-Chain. However, few tools have been developed for conducting quality control on alignment and variant calling.

We present QC3, a quality control tool designed for DNA sequencing data for all  the aforementioned three stages. QC3 provides both graphic and tabulated reports for quality control results. It also offers several unique features such as separation of bad and good reads (based on Illumina’s filter), detection of batch effect, and cross contamination. The input of QC3 takes three types of data: FASTQ, Binary Alignment Map (BAM), and Variant Calling Format, respectively, corresponding to the three stages of raw data, alignment, and variant detection. QC3 is written with Perl and R and is freely available for public use. It can be downloaded from [QC3 website on github](https://github.com/slzhao/QC3).

<a name="Download"/>
# Download #

You can directly download QC3 from [github](https://github.com/slzhao/QC3) by the following commands (If git has already been installed in your computer).

	#The source codes of QC3 software will be downloaded to your current directory
	git clone https://github.com/slzhao/QC3.git

Or you can also download the zip file of QC3 from [github](https://github.com/slzhao/QC3/archive/master.zip).

	#The zip file of QC3 software will be downloaded to your current directory
	wget https://github.com/slzhao/QC3/archive/master.zip
	#A directory named QC3-master will be generated and the source codes will be extracted there
	unzip master

<a name="Change"/>
# Change log #

<a name="R133">
## Release version 1.33 on October 19, 2015
Release version 1.33
 * Bam QC supports bam files with more than one alignments for one query;

<a name="R132">
## Release version 1.32 on June 26, 2015
Release version 1.32
 * Fix bug in of Bam QC;

<a name="R131">
## Release version 1.31 on May 01, 2015
Release version 1.31
 * Changes in methods of Vcf QC;
 * More robust for different formats of vcf files;

<a name="R130">
## Release version 1.30 on April 26, 2015
Release version 1.30
 * Vcf QC supports more different formats of vcf files;
 * Improvement of Vcf QC report;

<a name="R126">
## Release version 1.26 on April 23, 2015
Release version 1.26
 * More robust for different formats of vcf files;
 * new parameter -up for Vcf QC, Only use PASS variants in vcf;

<a name="R125">
## Release version 1.25 on April 16, 2015
Release version 1.25
 * More robust for different formats of vcf files;

<a name="R124">
## Release version 1.24 on April 15, 2015
Release version 1.24
 * More robust for different formats of vcf files;

<a name="R123">
## Release version 1.23 on October 20, 2014
Release version 1.23
 * A table with counts for different SNPs was added in vcf QC;
 * The SNPs on X, Y chromosome and Mitochondrion can be selected to keep or not;
 * More robust for different formats of vcf files;

<a name="R122">
## Release version 1.22 on July 30, 2014
Release version 1.22
 * Fix a bug when only one fastq or bam file is provided to fastq QC or Bam QC;

<a name="R121">
## Release version 1.21 on February 19, 2014
Release version 1.21
 * Vcf QC now supports vcf files generated by latest version of GATK;

<a name="R120">
## Release version 1.20 on January 24, 2014
Release version 1.20
 * Bam QC provided a statistic for alignment score generated by aligner(AS) if it was available in bam files;
 * Vcf QC provided statistics for reads on chromosome Y and heterozygous, Non-reference homozygous SNPs on chromosome X, which can be used for sex check;
 * The codes were improved to fix some issues for paper revision;
 * Documents were improved, an example from TCGA data was provided;

<a name="R115">
## Release version 1.15 on December 25, 2013
Release version 1.15
 * Input file for fastq QC and bam QC now supports label for each file;
 * Vcf QC now supports vcf files generated by latest version of GATK;
 * Vcf QC now supports vcf files with only two samples;
 * Documents were improved;

<a name="R110">
## Release version 1.10 on September 24, 2013
Release version 1.10
 * A configure file was provided so that the parameters of QC3 could be easily modified;
 * Some codes were updated to improve the performance;
 * Documents were improved;

<a name="R100">
## Release version 1.00 on August 09, 2013
Release version 1.00
 * Documents were improved;
 * The codes in bam QC were improved;

<a name="RC13">
## Release candidate (RC) version 1.3 on August 06, 2013
Release candidate version 1.3 for test
 * The ANNOVAR annotation function in vcf QC was changed;
 * The codes in bam QC and fastq QC were improved;

<a name="RC12">
## Release candidate (RC) version 1.2 on July 31, 2013
Release candidate version 1.2 for test
 * The algorithm in bam QC was improved and the memory usage was highly decreased;

<a name="RC11">
## Release candidate (RC) version 1.1 on July 29, 2013
Release candidate version 1.1 for test
 * Documents were improved;
 * Some bugs were fixed;
 * Example files were provided;

<a name="RC10">
## Release candidate (RC) version 1.0 on July 26, 2013
Release candidate version 1.0 for test

<a name="Prerequisites"/>
# Prerequisites #
<a name="irpp"/>
## Install Perl and required Perl packages ##

Perl is a highly capable, widely used, feature-rich programming language. It could be downloaded [Perl website](http://www.perl.org/get.html).

If Perl has already been installed on your computer, no other Perl module was needed to run QC3 in most cases. And you can run the following commands to make sure all the required modules have been installed.

	#Go to the folder where your QC3 software is.
	#And test whether all the required modules have been installed.
	bash test.modules

The successful output would look like this

	ok   File::Basename
	ok   FindBin
	ok   Getopt::Long
	ok   HTML::Template
	ok   source::bamSummary
	ok   source::fastqSummary
	ok   source::makeReport
	ok   source::vcfSummary
	ok   threads
	ok   threads::shared

Otherwise, for example, if HTML::Template package was missing, it may look like this

	ok   File::Basename
	ok   FindBin
	ok   Getopt::Long
	fail HTML::Template
	ok   source::bamSummary
	ok   source::fastqSummary
	ok   source::makeReport
	ok   source::vcfSummary
	ok   threads
	ok   threads::shared

Then you need to install the missing packages from [CPAN](http://www.cpan.org/). A program was also provided to make the package installation more convenient.

	#if HTML::Template was missing
	bash install.modules HTML::Template

<a name="irs"/>
## Install required software ##

### R ###

R is a free software environment for statistical computing and graphics. It could be downloaded from [R website](http://www.r-project.org/).

After you install R and add R bin file to your Path, the software can find and use R automatically. Or you can modify the config.txt file in the software directory and tell the program where the R is on your computer. Here is the line you need to modify.

	#where the R bin file is
	RBin="R"

### samtools ###

SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format. It could be downloaded from [SAM Tools website](http://samtools.sourceforge.net/).

After you install SAMtools and add SAMtools bin file to your Path, the software can find and use SAMtools automatically. Or you can modify the config.txt file in the software directory and tell the program where the SAMtools is on your computer. Here is the line you need to modify.

	#where the SAMtools bin file is
	samtoolsBin="samtools"


### annovar ###

**Please note: ANNOVAR and ANNOVAR database are not essential for vcf QC. If not provided or not found, ANNOVAR annotation will not be performed but other functions in vcf QC work well**.

ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes. We used ANNOVAR for annotation in vcf QC. It could be downloaded from [ANNOVAR website](http://www.openbioinformatics.org/annovar/).

After you install ANNOVAR and add ANNOVAR bin file to your Path, the software can find and use ANNOVAR automatically. Or you can modify the config.txt file in the software directory  and tell the program where the ANNOVAR is on your computer. Here is the line you need to modify.

	#where the ANNOVAR bin file is
	annovarBin="table_annovar.pl"
	#where the ANNOVAR convert bin file is
	annovarConvert="convert2annovar.pl"
	#other options for ANNOVAR
	annovarOption="-buildver hg19 -protocol refGene,snp137 -operation g,f --remove"

The ANNOVAR database files were also needed for ANNOVAR annotation. Please refer to the [ANNOVAR document](http://www.openbioinformatics.org/annovar/annovar_db.html) for the  database preparation. Or you can simply use the following commands to prepare ANNOVAR database files.

	#assume ANNOVAR was installed and added to Path. Then you want to download ANNOVAR database in annovarDatabase directory.
	mkdir annovarDatabase
	cd annovarDatabase
	annotate_variation.pl -downdb -buildver hg19 -webfrom annovar snp137 humandb/
	annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene humandb/

<a name="drd"/>
## Download required database ##
A gtf and/or a bed file was taken as database files in bam QC, and the position annotation for each sequence was exported from them. The chromosomes and positions information in these files should be exactly same with the bam files (for example, chromosome 1 can't be represented as chr1 in one file but 1 in another file). At least one of gtf/bed should be provided in bam QC. These files could be downloaded at [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables?command=start). The format of these files were:

	#bed file:
	#Column1	    Column2         Column3
	Chromosome	StartPosition	EndPosition

	#gtf file:
	#Column1   	...	Column3	Column4	        Column5	...
	Chromosome	...	exon	StartPosition	EndPosition	...

<a name="Usage"/>
# Usage #
The usage of QC3 software could be:

	perl qc3.pl -m module -i inputFile -o outputDir [-t threads -rp] [some other parameters]

	-m [module]         Required. QC module used. It should be f (fastq QC), b (bam QC), or v (vcf QC).
	-i [inputFile]      Required. Input file. In fastq QC or bam QC, it should be a file listing all analyzed files (supports .fastq, .fastq.gz, and .bam files). To analyze the pair-end fastq files in fastq QC, the two files for the same sample should be listed together in this file. In vcf QC, it should be a vcf file.
	-o [outputDir]      Required. Output directory for QC result. If the directory doesn't exist, it would be created.
	-t [int]            Optional. Threads used in analysis. The default value is 4. This parameter only valid for fastq and bam QC. Only one thread would be used in vcf QC.
	-rp					Optional. Re-plot the figures. The program will not re-analyze the input files but the result files in output directory will be used to re-plot the figures and re-generate the report if this parameter is used. You can use this parameter when you want to slightly modify the report, such as changing sample names in the figures.
	-h	                Optional. Show help information for each module.

<a name="if"/>
## input file
In fastq QC or bam QC, the input file should be a file listing all analyzed files (supports .fastq, .fastq.gz for fastq QC, and .bam for bam QC). To analyze the pair-end fastq files in fastq QC, the two files for the same sample should be listed together in this file. In vcf QC, the input file should be a vcf file.

Input file for fastq QC or bam QC also supports label for each file listed (The labels are optional). The label should follow the file and separated by a Tab. The labels will be used in the report instead of the file names. So that the report will be much easier to understand. An example was listed here:

	#An example of input file for fastq QC. The labels are optional.
	#Column1	    Column2
	sample1Pair1File	labelForThisFile
	sample1Pair2File	labelForThisFile
	Sample2Pair1File	labelForThisFile

Here is more details for each module:

<a name="fqc"/>
## fastq QC

	perl qc3.pl -m f -i inputFileList -o outputDir -t threads

	-se				Optional. Indicating the fastq files were single-end data. QC3 takes fastq files as pair-end data by default.

<a name="bqc"/>
## bam QC

	perl qc3.pl -m b -i inputFileList -o outputDir -t threads

	-r  [database]	Optional. A targetregion file. At least one targetregion file or gtf file should be provided.
	-g  [database]	Optional. A gtf file. At least one targetregion file or gtf file should be provided.
	-cm [int]	    Optional. Calculation method for data summary, should be 1 or 2. Method 1 means mean and method 2 means median. The default value is 1.
	-d				Optional. The depth in on-/off-target regions will be calculated. QC3 will not calculate depth by default, because it may take a long time.
	-nod            Optional. Do not compute the depth in off-target regions for saving time (no off-target depth) If set, the user must provide a bed file for target regions.
	-no_batch       Optional. The batch effect part of figures will not be plotted (considering BAM file not informative).
	-use_SM         Optional. Extract SM field from bam file and use it as sample name.
	-d_cumul        Optional. Depth values for cumulative distribution computation. Put 3 values if format: -d_cumul 0,10,30 (default values).

**Please note: At least one targetregion file or gtf file should be provided for bam QC**. So the user can use (-r) or (-g) or (-r and -g) option for bam QC.

<a name="vqc"/>
## vcf QC

	perl qc3.pl -m v -i inputFile -o outputDir

	-s [int]			Optional. Method used in consistence calculation, should be 1 or 2. In method 1, only the two samples with completely same allele will be taken as consistent. In method 2, the two samples satisfy the criterion in method 1, or if the two samples were 0/0 vs 0/1 but 0/0 sample has some read counts in alternative allele, they will be taken as consistent. Besides, consistence calculation in either method 1 or 2 only uses the alleles with at least 10 read depths in any of the two samples. The default value is 1.
	-c [database]	Optional. A file indicating the filter arguments for vcf files. If not specified, the default file 'GATK.cfg' in QC3 directory with GATK best practices recommended arguments will be used.
	-a [database]	Optional. Directory of annovar database.

<a name="sExample"/>
# A short Example #

The example files can be downloaded at [sourceforge](http://sourceforge.net/projects/qc3/files/).

You need to download and extract it to a directory. Then the example code for running QC3 with given example data set could be:

<a name="ded"/>
## download example data ##
	#download and extract example data into exampleDir
	mkdir exampleDir
	cd exampleDir
	wget http://sourceforge.net/projects/qc3/files/example.tar.gz/download
	tar zxvf example.tar.gz
	ls

<a name="fqcE"/>
## fastq QC ##

	#assume qc3.pl in qc3Dir, examples in exampleDir
	cd exampleDir
	cd fastq
	perl qc3Dir/qc3.pl -m f -i example_fastq_list.txt -o example_fastq_result


<a name="bqcE"/>
## bam QC

	#assume qc3.pl in qc3Dir, examples in exampleDir
	cd exampleDir
	cd bam

	#a simple example, with only bed file
	perl qc3Dir/qc3.pl -m b -i example_bam_list.txt -r ../bamQcDatabase/hg19_protein_coding.bed.Part -o example_bam_result1

	#a simple example, with only gtf file
	perl qc3Dir/qc3.pl -m b -i example_bam_list.txt -g ../bamQcDatabase/Homo_sapiens.GRCh37.63_protein_coding_chr1-22-X-Y-M.gtf.Part -o example_bam_result2

	#a example with both bed and gtf file, and will calculate the depth
	perl qc3Dir/qc3.pl -m b -i example_bam_list.txt -r ../bamQcDatabase/hg19_protein_coding.bed.Part -g ../bamQcDatabase/Homo_sapiens.GRCh37.63_protein_coding_chr1-22-X-Y-M.gtf.Part -o example_bam_result3

<a name="vqcE"/>
## vcf QC

	#assume qc3.pl in qc3Dir, examples in exampleDir
	cd exampleDir
	cd vcf

	#a simple example, will not perform ANNOVAR annotation
	perl qc3Dir/qc3.pl -m v -i CV-7261.vcf.top40000.txt -o example_vcf_result1

	#assume ANNOVAR database in annovarDatabase, will perform ANNOVAR annotation
	perl qc3Dir/qc3.pl -m v -i CV-7261.vcf.top40000.txt -s 2 -a annovarDatabase -o example_vcf_result2

<a name="TExample"/>
# An example from TCGA data #

In this part, we used 11 bam files from TCGA BRCA data as examples. These files was downloaded by [GeneTorrent](https://cghub.ucsc.edu/software/downloads.html). The file names were:

	TCGA-A7-A0D9-01A-31W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam
	TCGA-BH-A0B3-01A-11W-A071-09_IlluminaGA-DNASeq_exome.bam
	TCGA-BH-A0B8-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam
	TCGA-BH-A0BJ-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam
	TCGA-BH-A0BM-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam
	TCGA-BH-A0C0-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam
	TCGA-BH-A0DK-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam
	TCGA-BH-A0DP-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam
	TCGA-BH-A0DP-10A-01W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam
	TCGA-BH-A0E0-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam
	TCGA-BH-A0H7-01A-13W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam

These 11 bam files including 10 primary solid tumor samples and 1 blood derived normal samples. The tumor sample TCGA-BH-A0DP-01A and normal sample TCGA-BH-A0DP-10A were from same patient.

<a name="TfqcE"/>
## fastq QC ##
The fastq files were first generated based on the 11 bam files. Then fastq QC results were generated with following codes:

	#assume qc3.pl in qc3Dir, fastq files in fastqDir and listed in 11fastq.list file
	cd fastqDir
	perl qc3Dir/qc3.pl -m f -t 11 -i 11fastq.list -o 11BRCA_fastq_result

The fastq QC results can be visited [Here](http://htmlpreview.github.io/?https://github.com/slzhao/QC3/blob/master/reportExample/fastq/fastqReport.html)

<a name="TbqcE"/>
## bam QC
The results were generated with following codes:

	#assume qc3.pl in qc3Dir, bam files in bamDir and listed in 11bam.list file
	cd bamDir
	perl qc3Dir/qc3.pl -m b -t 11 -i 11bam.list -r hg19BedFile -g hg19GtfFile -o 11BRCA_bam_result

The bam QC results can be visited [Here](http://htmlpreview.github.io/?https://github.com/slzhao/QC3/blob/master/reportExample/bam/bamReport.html)

<a name="TvqcE"/>
## vcf QC
The vcf file was first generated based on the 11 bam files by GATK. Then vcf QC results were generated with following codes:

	#assume qc3.pl in qc3Dir, vcf files in vcfDir
	cd vcfDir
	#assume ANNOVAR database in annovarDatabase, will perform ANNOVAR annotation
	perl qc3Dir/qc3.pl -m v -i vcfFile -a annovarDatabase -o 11BRCA_vcf_result

The vcf QC results can be visited [Here](http://htmlpreview.github.io/?https://github.com/slzhao/QC3/blob/master/reportExample/vcf/vcfReport.html)

<a name="Results"/>
# Results #
The first section in all reports was the date and command for report generation. So that the user can easily reproduce the report.

<a name="ifr"/>
## Introduction for fastq QC result ##

An example report could be accessed at [Here](http://htmlpreview.github.io/?https://github.com/slzhao/QC3/blob/master/reportExample/fastq/fastqReport.html)

fastq QC report was constituted by three sections.

The "Result Table" section displayed the instrument, run number, flowcell, lane, total reads, base quality score and GC content for each fastq file.

The "Sequence quality and content" section indicated the per base sequence quality and per base nucleotide content for each fastq file (For pair end results the two files indicating one sample were plotted in the same figure).

The "Batch effect" section indicated the statistics of reads, BQ and GC by different runs/machines/Flowcells/Lanes to demonstrate if there were batch effects in all experiments. Boxplot was used to visualize the distributions of statistics in different groups. Kruskal-Wallis rank sum test and Fligner-Killeen test was used to determine the significance of distribution in different groups.

<a name="ibr"/>
## Introduction for bam QC result ##

An example report could be accessed at [Here](http://htmlpreview.github.io/?https://github.com/slzhao/QC3/blob/master/reportExample/bam/bamReport.html)

bam QC report was constituted by three sections.

The "Result Table" section displayed the instrument, run number, flowcell, lane, total reads, on/off target reads, intron/intergenic/mito reads,  for each bam file.

The "Distribution" section visualized the distributions of statistics in all experiments.

The "Batch effect" section indicated the statistics of reads, on-target and off-target by different runs/machines/Flowcells/Lanes to demonstrate if there were batch effects in all experiments. Boxplot was used to visualize the distributions of statistics in different groups. Kruskal-Wallis rank sum test and Fligner-Killeen test was used to determine the significance of distribution in different groups.

<a name="ivr"/>
## Introduction for vcf QC result ##

An example report could be accessed at [Here](http://htmlpreview.github.io/?https://github.com/slzhao/QC3/blob/master/reportExample/vcf/vcfReport.html)

vcf QC report was constituted by seven sections.

The "Filter" section displayed the arguments for the filter of vcf files.

The "Statistics" section displayed the Transitions:Transversions ratio and Heterozygous:Non-reference homozygous ratio in all samples.

The "Consistency" section visualized the consistency in all sample to sample pairs.

The "SNP count" section displayed the SNP count of each sample in different chromesomes.

The "Statistics for sex check" section displayed the Heterozygous and Non-reference Homozygous SNP counts in chromosome X, and total reads in chromosome Y, which can be used for sex check.

The "Score" section visualized the position and score of SNPs before and after filter.

The "Annovar annotation" section displayed the counts of SNPs with different functions, such as synonymous/nonsynonymous SNV, and whether they were in snp137 database.
