Table of Content
================
* [Overview](#Overview)
* [Change log](#Change)
 * [Release version 1.0](#RV10)
* [Prerequisites](#Prerequisites)
 * [Install required perl packages](#irpp)
 * [Install required software](#irs)
 * [Download required database](#drd)
* [Usage](#Usage)
 * [fastq QC](#fqc)
 * [bam QC](#bqc)
 * [vcf QC](#vqc)
* [Example](#Example)
 * [fastq QC example](#fqcE)
 * [bam QC example](#bqcE)
 * [vcf QC example](#vqcE)
* [Results](#Results)
 * [Introduction for fastq QC result](#ifr)
 * [Introduction for bam QC result](#ibr)
 * [Introduction for vcf QC result](#iqr)
* [Others](#Others)


# Overview #
Say something

# Change log #
Say something

# Prerequisites #
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

## Install required software ##

### R ###

R is a free software environment for statistical computing and graphics. It could be downloaded [here](http://www.r-project.org/).

After you install R and add R bin file to your Path, the software can find and use R automatically. Or you can change the main.pl script and tell the program where the R is on your computer. Here is the line you need to modify.

	$config{'RBin'}        = "R";       #where the R bin file is

### samtools ###

SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format. It could be downloaded [here](http://samtools.sourceforge.net/).

After you install SAMtools and add SAMtools bin file to your Path, the software can find and use SAMtools automatically. Or you can change the main.pl script and tell the program where the SAMtools is on your computer. Here is the line you need to modify.

	$config{'samtoolsBin'} = "samtools";	#where the SAMtools bin file is


### annovar ###

ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes. We used ANNOVAR for annotation in vcf summary. It could be downloaded from[ANNOVAR website](http://www.openbioinformatics.org/annovar/).

The ANNOVAR database files were also needed for ANNOVAR annotation. Please refer to the [ANNOVAR document](http://www.openbioinformatics.org/annovar/annovar_db.html) for the  database preparation.

After you install ANNOVAR and add ANNOVAR bin file to your Path, the software can find and use ANNOVAR automatically. Or you can change the main.pl script and tell the program where the ANNOVAR is on your computer. Here is the line you need to modify.

	$config{'annovarBin'}  = "~/bin/annovar/annotate_variation.pl";		#where the ANNOVAR bin file is

Please note: ANNOVAR and ANNOVAR database are not essential for vcf QC. If not provided or not found, ANNOVAR annotation will not be performed in vcf QC.

## Download required database ##
These database files were needed in bam QC. At least one of them should be provided in bam QC. They could be downloaded from:

Homo\_sapiens.GRCh37.63\_protein\_coding\_chr1-22-X-Y-M.gtf

hg19\_protein\_coding.bed


# Usage #
Example code for running QC with given example dataset could be:

	perl main.pl -m module -i inputFile -o outputDir (some other parameters)
	
	-m [module]         QC module used. It should be f (fastq QC), b (bam QC), or v (vcf QC).
	-i [inputFile]      Input file. It should be a file list including all analyzed files in fastq QC and bam QC. In vcf QC, it should be the vcf file.
	-o [outputDir]      Output directory for QC result. If the directory doesn't exist, it would be created. If the directory already existed, the files in it would be deleted.
	

Here is more details for each module:

## fastq QC 

	perl main.pl -m f -i inputFileList -o outputDir

	-t [int]        Threads used in analysis. The default value is 4. This parameter only valid for fastq QC. Only one thread would be used in bam and vcf QC.


## bam QC

	perl main.pl -m b -i inputFileList -o outputDir

	-r [database]	A targetregion file
	-g [database]	A gtf file
	-d				whether the depth in on-/off-target regions will be calculated, -d = will be calculated

## vcf QC

	perl main.pl -m v -i inputFile -o outputDir

	-s [int]		Method used in consistence calculation, should be 1 or 2. In method 1, only the two samples will completely same allele will be taken as consist. In method 2, the two samples satisfy the criterion in method 1 or ~~~ will be taken as consist.

# Example #

The example files can be downloaded at 121221.

Example code for running QC with given example data set could be:

## fastq QC ##

	perl main.pl -m f -i example_fastq_fileList.txt -o example_fastq_result



## bam QC

	perl main.pl -m b -i example_bam_fileList.txt -t /data/cqs/guoy1/reference/annotation/hg19/hg19_protein_coding.bed -g /data/cqs/guoy1/reference/annotation/hg19/Homo_sapiens.GRCh37.63_protein_coding_chr1-22-X-Y-M.gtf -o example_bam_result -d 1

## vcf QC

	perl main.pl -m v -i /scratch/cqs/zhaos/QC/test/example/vcf/HNSC_TCGA-CV-7261.vcf.top4000.txt -s 1 -o example_vcf_result

# Results #
The first section in all reports was the date and command for report generation. So that the user can easily reproduce the report.
 
## Introduction for fastq QC result ##
fastq QC report was constituted by two sections. 

The "Base" section indicated the per base sequence quality and per base nucleotide content for each fastq files (For pair end results the two files indicating one sample were plotted in the same figure).

The "Batch effect" section indicated the statistics of reads, BQ and GC by different runs/machines/Flowcells/Lanes to demonstrate if there were batch effects in all experiments. Boxplot was used to visualize the distributions of statistics in different groups. Kruskal-Wallis rank sum test and Fligner-Killeen test was used to determine the significance of distribution in different groups.

## Introduction for bam QC result ##

The "Distribution" section visualized the distributions of statistics in all experiments.

The "Batch effect" section indicated the statistics of reads, on-target and off-target by different runs/machines/Flowcells/Lanes to demonstrate if there were batch effects in all experiments. Boxplot was used to visualize the distributions of statistics in different groups. Kruskal-Wallis rank sum test and Fligner-Killeen test was used to determine the significance of distribution in different groups.

## Introduction for vcf QC result ##
The "Statistics" section visualized the Transitions:Transversions ratio and 0/1:1/1 ratio in all samples.

The "Consistence" section visualized the consistence in all samples.

The "Score" section visualized the score before and after filter.

The "Annovar" section listed the annovar annotation result files.
