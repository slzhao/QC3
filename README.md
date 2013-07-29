Table of Content
================
* [Overview](#Overview)
* [Change log](#Change)
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
 * [Introduction for vcf QC result](#iqr)
* [Others](#Others)

<a name="Overview"/>
# Overview #
Say something

<a name="Change"/>
# Change log #

<a name="RC11">
## Release candidate (RC) version 1.1 on July 29, 2013
Release version 1.1 for test
 * Documents were improved;
 * Some bugs were fixed;
 * Example files were provided;

<a name="RC10">
## Release candidate (RC) version 1.0 on July 26, 2013
Release version 1.0 for test

<a name="Prerequisites"/>
# Prerequisites #
<a name="irpp"/>
## Install required perl packages ##

There are several other packages needed to be installed on your computer.

	#go the the folder where your software is.
	#And test whether all the required modules have been installed.
	./test.modules

The successful output would look like this

	ok   File::Basename
	ok   FindBin
	ok   Getopt::Long
	ok   HTML::Template
	ok   source::bamSummary
	ok   source::fastqSummary
	ok   source::makeReport
	ok   source::vcfSummary

Otherwise, it may look like this

	ok   File::Basename
	ok   FindBin
	ok   Getopt::Long
	fail HTML::Template is not usable (it or a sub-module is missing)
	ok   source::bamSummary
	ok   source::fastqSummary
	ok   source::makeReport
	ok   source::vcfSummary

Then you need to install the missing packages from [CPAN](http://www.cpan.org/). 

<a name="irs"/>
## Install required software ##

### R ###

R is a free software environment for statistical computing and graphics. It could be downloaded [here](http://www.r-project.org/).

After you install R and add R bin file to your Path, the software can find and use R automatically. Or you can change the qc3.pl script and tell the program where the R is on your computer. Here is the line you need to modify.

	$config{'RBin'}        = "R";       #where the R bin file is

### samtools ###

SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format. It could be downloaded [here](http://samtools.sourceforge.net/).

After you install SAMtools and add SAMtools bin file to your Path, the software can find and use SAMtools automatically. Or you can change the qc3.pl script and tell the program where the SAMtools is on your computer. Here is the line you need to modify.

	$config{'samtoolsBin'} = "samtools";	#where the SAMtools bin file is


### annovar ###

ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes. We used ANNOVAR for annotation in vcf summary. It could be downloaded from [ANNOVAR website](http://www.openbioinformatics.org/annovar/).

The ANNOVAR database files were also needed for ANNOVAR annotation. Please refer to the [ANNOVAR document](http://www.openbioinformatics.org/annovar/annovar_db.html) for the  database preparation.

After you install ANNOVAR and add ANNOVAR bin file to your Path, the software can find and use ANNOVAR automatically. Or you can change the qc3.pl script and tell the program where the ANNOVAR is on your computer. Here is the line you need to modify.

	$config{'annovarBin'}  = "annotate_variation.pl";		#where the ANNOVAR bin file is

**Please note: ANNOVAR and ANNOVAR database are not essential for vcf QC. If not provided or not found, ANNOVAR annotation will not be performed in vcf QC**.

<a name="drd"/>
## Download required database ##
A gtf and/or a bed file was taken as database files in bam QC, and the position annotation for each sequence was exported from them. At least one of them should be provided in bam QC. The format of these files were: 

	hg19_protein_coding.bed
	Column1	    Column2         Column3
	Chromosome	StartPosition	EndPosition

	Homo_sapiens.GRCh37.63_protein_coding_chr1-22-X-Y-M.gtf
	Column1   	...	Column3	Column4	        Column5	...
	Chromosome	...	exon	StartPosition	EndPosition	...

<a name="Usage"/>
# Usage #
Example code for running QC with given example dataset could be:

	perl qc3.pl -m module -i inputFile -o outputDir (some other parameters)
	
	-m [module]         QC module used. It should be f (fastq QC), b (bam QC), or v (vcf QC).
	-i [inputFile]      Input file. It should be a file list including all analyzed files in fastq QC and bam QC. In vcf QC, it should be the vcf file.
	-o [outputDir]      Output directory for QC result. If the directory doesn't exist, it would be created. If the directory already existed, the files in it would be deleted.
	

Here is more details for each module:

<a name="fqc"/>
## fastq QC 

	perl qc3.pl -m f -i inputFileList -o outputDir

	-t [int]        Threads used in analysis. The default value is 4. This parameter only valid for fastq QC. Only one thread would be used in bam and vcf QC.
	-p				whether the fastq files were pair-end data. -p = pair-end.

<a name="bqc"/>
## bam QC

	perl qc3.pl -m b -i inputFileList -o outputDir

	-r  [database]	A targetregion file.
	-g  [database]	A gtf file.
	-cm [int]	    Calculation method for data summary, should be 1 or 2. Method 1 means mean and method 2 means median.
	-d				whether the depth in on-/off-target regions will be calculated, -d = will be calculated.

<a name="vqc"/>
## vcf QC

	perl qc3.pl -m v -i inputFile -o outputDir

	-s [int]		    Method used in consistence calculation, should be 1 or 2. In method 1, only the two samples will completely same allele will be taken as consist. In method 2, the two samples satisfy the criterion in method 1 or ~~~ will be taken as consist.
	-a  [database]	Directory of annovar database.

<a name="Example"/>
# Example #

The example files can be downloaded at [sourceforge](http://sourceforge.net/projects/qc3/files/).

You need to download and extract it to a directory, such as exampleDir. Then the example code for running QC with given example data set could be:

<a name="ded"/>
## download example data ##
	#download and extract example data
	cd exampleDir
	wget http://sourceforge.net/projects/qc3/files/example.tar.gz/download
	tar zxvf example.tar.gz
	ls
	
	#make file list file for fastq and bam QC
	bash makeFileList.sh
	head example_bam_fileList.txt
	head example_fastq_fileList.txt
	#example_bam_fileList.txt and example_fastq_fileList.txt would be generated in exampleDir

<a name="fqcE"/>
## fastq QC ##

	#assume qc3.pl in qc3Dir, examples in exampleDir
	cd qc3Dir
	perl qc3.pl -m f -t 8 -i exampleDir/example_fastq_fileList.txt -o exampleDir/example_fastq_result -p


<a name="bqcE"/>
## bam QC

	#assume qc3.pl in qc3Dir, examples in exampleDir
	cd qc3Dir
	perl qc3.pl -m b -i exampleDir/example_bam_fileList.txt -r exampleDir/bamQcDatabase/hg19_protein_coding.bed.Part -g exampleDir/bamQcDatabase/Homo_sapiens.GRCh37.63_protein_coding_chr1-22-X-Y-M.gtf.Part -o exampleDir/example_bam_result -d

<a name="vqcE"/>
## vcf QC
	
	#assume qc3.pl in qc3Dir, examples in exampleDir
	perl qc3.pl -m v -i exampleDir/vcf/CV-7261.vcf.top40000.txt -s 1 -o exampleDir/example_vcf_result

<a name="Results"/>
# Results #
The first section in all reports was the date and command for report generation. So that the user can easily reproduce the report.
 
<a name="ifr"/>
## Introduction for fastq QC result ##
fastq QC report was constituted by two sections. 

The "Base" section indicated the per base sequence quality and per base nucleotide content for each fastq files (For pair end results the two files indicating one sample were plotted in the same figure).

The "Batch effect" section indicated the statistics of reads, BQ and GC by different runs/machines/Flowcells/Lanes to demonstrate if there were batch effects in all experiments. Boxplot was used to visualize the distributions of statistics in different groups. Kruskal-Wallis rank sum test and Fligner-Killeen test was used to determine the significance of distribution in different groups.

<a name="ibr"/>
## Introduction for bam QC result ##

The "Distribution" section visualized the distributions of statistics in all experiments.

The "Batch effect" section indicated the statistics of reads, on-target and off-target by different runs/machines/Flowcells/Lanes to demonstrate if there were batch effects in all experiments. Boxplot was used to visualize the distributions of statistics in different groups. Kruskal-Wallis rank sum test and Fligner-Killeen test was used to determine the significance of distribution in different groups.

<a name="ivr"/>
## Introduction for vcf QC result ##
The "Statistics" section visualized the Transitions:Transversions ratio and 0/1:1/1 ratio in all samples.

The "Consistence" section visualized the consistence in all samples.

The "Score" section visualized the score before and after filter.

The "Annovar" section listed the annovar annotation result files.

<a name="Others"/>
# Others #
