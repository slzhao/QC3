Table of Content
================
* [Introduction](#Introduction)
* [Download](#download)
* [Change log](#Change)
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
 * [fastq QC](#fqc)
 * [bam QC](#bqc)
 * [vcf QC](#vqc)
* [Example](#Example)
 * [Download example data](#ded)
 * [fastq QC example](#fqcE)
 * [bam QC example](#bqcE)
 * [vcf QC example](#vqcE)
* [Results](#Results)
 * [Introduction for fastq QC result](#ifr)
 * [Introduction for bam QC result](#ibr)
 * [Introduction for vcf QC result](#ivr)
* [Contacts](#Contacts)

<a name="Introduction"/>
# Introduction #

High throughput sequencing is the most effective way to screen for non-specific germline variants, somatic mutations, and structural variants. Some of the most popular sequencing paradigms in DNA sequencing are whole genome sequencing, exome sequencing, and target panel sequencing. While vastly informative, sequencing data poses significant bioinformatics challenges in areas such as data storage, computation time, and variant detection accuracy.  One of the easily overlooked challenges associated with sequencing is quality control Quality control for DNA sequencing data can be categorized into three stages: raw data, alignment, and variant calling. QC on raw sequencing data has been given more attention than QC on alignment and the variant calling. There are many QC tools aimed at raw data such as FastQC, FastQ Screen FastX-Toolkit, NGS QC Toolkit, RRINSEQ and QC-Chain. However, few tools have been developed for conducting quality control on alignment and variant calling.

We present QC3, a quality control tool designed for DNA sequencing data for all aforementioned three stages. QC3 provides both graphic and tabulated reports for quality control results. It also offers several unique features such as separation of bad and good reads (based on Illumina’s filter), detection of batch effect, and cross contamination. The input of QC3 can be three types of data: FASTQ, Binary Alignment Map (BAM), and Variant Calling Format, respectively, corresponding to the three stages of raw data, alignment, and variant detection. QC3 is written with Perl and R and is freely available for public use. It can be downloaded from [github](https://github.com/slzhao/QC3).

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

Perl is a highly capable, widely used, feature-rich programming language. It could be downloaded [here](http://www.perl.org/get.html).

If Perl has already been installed on your computer, no other Perl module was needed to run QC3 in most cases. And you can run the following commands to make sure all the required modules have been installed.

	#go the the folder where your QC3 software is.
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

Otherwise, it may look like this

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

R is a free software environment for statistical computing and graphics. It could be downloaded [here](http://www.r-project.org/).

After you install R and add R bin file to your Path, the software can find and use R automatically. Or you can modify the config.txt file in the software directory and tell the program where the R is on your computer. Here is the line you need to modify.

	#where the R bin file is
	RBin="R"

### samtools ###

SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format. It could be downloaded [here](http://samtools.sourceforge.net/).

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

	bed file:
	Column1	    Column2         Column3
	Chromosome	StartPosition	EndPosition

	gtf file:
	Column1   	...	Column3	Column4	        Column5	...
	Chromosome	...	exon	StartPosition	EndPosition	...

<a name="Usage"/>
# Usage #
The usage of QC3 software could be:

	perl qc3.pl -m module -i inputFile -o outputDir [-t threads -rp] [some other parameters]
	
	-m [module]         Required. QC module used. It should be f (fastq QC), b (bam QC), or v (vcf QC).
	-i [inputFile]      Required. Input file. It should be a file list including all analyzed files in fastq QC and bam QC. To analyze the pair-end fastq files in fastq QC, the two files for the same sample should be listed together in this file. In vcf QC, it should be a vcf file.
	-o [outputDir]      Required. Output directory for QC result. If the directory doesn't exist, it would be created.
	-t [int]            Optional. Threads used in analysis. The default value is 4. This parameter only valid for fastq and bam QC. Only one thread would be used in vcf QC.
	-rp					Optional. Re-plot the figures. The program will not re-analyze the input files but the result files in output directory will be used to re-plot the figures and re-generate the report if this parameter is used. You can use this parameter when you want to slightly modify the report, such as changing sample names in the figures.
	-h	                Optional. Show help information for each module.

Here is more details for each module:

<a name="fqc"/>
## fastq QC 

	perl qc3.pl -m f -i inputFileList -o outputDir -t threads

	-se				Optional. whether the fastq files were pair-end data. They were taken as pair-end by default, -se = single-end.

<a name="bqc"/>
## bam QC

	perl qc3.pl -m b -i inputFileList -o outputDir -t threads

	-r  [database]	Optional. A targetregion file. At least one targetregion file or gtf file should be provided.
	-g  [database]	Optional. A gtf file. At least one targetregion file or gtf file should be provided.
	-cm [int]	    Optional. Calculation method for data summary, should be 1 or 2. Method 1 means mean and method 2 means median. The default value is 1.
	-d				Optional. whether the depth in on-/off-target regions will be calculated, It will not be calculated by default, -d = will be calculated.

<a name="vqc"/>
## vcf QC

	perl qc3.pl -m v -i inputFile -o outputDir

	-s [int]			Optional. Method used in consistence calculation, should be 1 or 2. First of all, if any of the two samples had less than 10 read depths in an allele, it will not be used in consistence calculation. In method 1, only the two samples with completely same allele will be taken as consistent. In method 2, the two samples satisfy the criterion in method 1, or if the two samples were 0/0 vs 0/1 but 0/0 sample has some read counts in alternative allele, they will be taken as consistent. The default value is 1.
	-c [database]	Optional. A file indicating the filter arguments for vcf files. If not specified, the default file 'GATK.cfg' in QC3 directory with GATK best practices recommended arguments will be used.
	-a [database]	Optional. Directory of annovar database.

<a name="Example"/>
# Example #

The example files can be downloaded at [sourceforge](http://sourceforge.net/projects/qc3/files/).

You need to download and extract it to a directory. Then the example code for running QC with given example data set could be:

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
	
	#a example with both bed and gtf file, and will calculate the depth
	perl qc3Dir/qc3.pl -m b -i example_bam_list.txt -r ../bamQcDatabase/hg19_protein_coding.bed.Part -g ../bamQcDatabase/Homo_sapiens.GRCh37.63_protein_coding_chr1-22-X-Y-M.gtf.Part -o example_bam_result2

<a name="vqcE"/>
## vcf QC
	
	#assume qc3.pl in qc3Dir, examples in exampleDir
	cd exampleDir
	cd vcf

	#a simple example, will not perform ANNOVAR annotation
	perl qc3Dir/qc3.pl -m v -i CV-7261.vcf.top40000.txt -o example_vcf_result1
	
	#assume ANNOVAR database in annovarDatabase, will perform ANNOVAR annotation
	perl qc3Dir/qc3.pl -m v -i CV-7261.vcf.top40000.txt -s 2 -a annovarDatabase -o example_vcf_result2

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

vcf QC report was constituted by six sections. 

The "Filter" section displayed the arguments for the filter of vcf files.

The "Statistics" section displayed the Transitions:Transversions ratio and Heterozygous:Non-reference homozygous ratio in all samples.

The "Consistency" section visualized the consistency in all sample to sample pairs.

The "SNP count" section displayed the SNP count of each sample in different chromesomes.

The "Score" section visualized the position and score of SNPs before and after filter.

The "Annovar annotation" section displayed the counts of SNPs with different functions, such as synonymous/nonsynonymous SNV, and whether they were in snp137 database.

<a name="Contacts"/>
# Contacts #
Shilin Zhao:	shilin.zhao@vanderbilt.edu

Yan Guo:	    yan.guo@vanderbilt.edu