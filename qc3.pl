#!/usr/bin/perl
use strict;
use warnings;
use threads;
use threads::shared;
use File::Basename;
use Getopt::Long;

use FindBin;
use lib $FindBin::Bin;
use HTML::Template;
use source::makeReport;
use source::fastqSummary;
use source::bamSummary;
use source::vcfSummary;

our $version = "1.34";

my $qc3ConfigFile = dirname($0) . "/config.txt";
my %config;

my $usageModule = "
Program: qc3.pl (a quality control tool for DNA sequencing data in raw data, alignment, and variant calling stages)
Version: $version

Usage:   perl qc3.pl -m module -i listfile -o outputDirectory [-t threads] [other options]

Options:

	-m	module             Required. QC module used. It should be f (fastq QC), b (bam QC), or v (vcf QC).
	-i	input filelist     Required. Input file. In fastq QC or bam QC, it should be a file listing all analyzed files (supports .fastq, .fastq.gz, and .bam files). To analyze the pair-end fastq files in fastq QC, the two files for the same sample should be listed together in this file. In vcf QC, it should be a vcf file.
	-o	output directory   Required. Output directory for QC result. If the directory doesn't exist, it would be created.
	-t	threads            Optional. Threads used in analysis. The default value is 4. This parameter only valid for fastq and bam QC. Only one thread will be used in vcf QC.
	-h	help               Optional. Show this information.

To see help information for each module, please use -m parameter to specify module.
For more information, please refer to the readme file in QC3 directory. Or visit the QC3 website at https://github.com/slzhao/QC3

";

$| = 1;
my $commandline = "perl $0 "
  . join( " ", @ARGV )
  ;    #Output you input command in the report to reproduce you result
my (
	$module,           $filelist,    $resultDir,  $singleEnd,
	$targetregionfile, $gtffile,     $isdepth,    $nod,
        $caculateMethod,   $vcfCfgFile,  $method,     $maxThreads,
        $annovarDb,        $usePASS,     $xym,        $rePlot,
        $no_batch,         $use_SM,      $d_cumul,		$showHelp
) = ();
our @log : shared;

GetOptions(
	"t=i" => \$maxThreads,
	"m=s" => \$module,
	"i=s" => \$filelist,
	"o=s" => \$resultDir,
	"se"  => \$singleEnd,

	"r=s"  => \$targetregionfile,
	"g=s"  => \$gtffile,
	"d"    => \$isdepth,
	"nod"  => \$nod,
	"d_cumul=s" => \$d_cumul,
	"no_batch"  => \$no_batch,
	"use_SM"	=>	\$use_SM,
	"cm=i" => \$caculateMethod,

	"c:s" => \$vcfCfgFile,
	"s=i" => \$method,
	"a=s" => \$annovarDb,
	"up"    => \$usePASS,
	"xym"    => \$xym,

	"rp" => \$rePlot,
	"h"  => \$showHelp,
);
if ( defined $module and $showHelp ) {
	$config{'showHelp'} = $showHelp;
	if ( $module eq "f" ) {
		my $rResult = &fastqSummary( $filelist, \%config );
	}
	elsif ( $module eq "b" ) {
		my $rResult = &bamSummary( $filelist, \%config );
	}
	elsif ( $module eq "v" ) {
		my $rResult = &vcfSummary( $filelist, \%config );
	}
}
elsif ($showHelp) { die "$usageModule"; }
open QC3CFG, "<$qc3ConfigFile" or die "Can't read $qc3ConfigFile\n$!";
while (<QC3CFG>) {
	chomp;
	if (/^#/) {
		next;
	}
	my @lines = ( split /[="]/, $_ );
	$config{ $lines[0] } = $lines[2];
}
my $detectR = `which $config{"RBin"}`;
if ( $detectR eq '' ) {
	die
"Can't find R. Please install R or modify the RBin in config.txt.\n$usageModule";
}
if ( !defined $module
	or ( $module ne 'f' and $module ne 'b' and $module ne 'v' ) )
{
	die(
"Module (-m) is required and must be f (fastq), b (bam), and v (vcf)\n$usageModule"
	);
}
elsif ( !defined $filelist or !defined $resultDir ) {
	die(
"Input file (-i) and Output directory (-o) must be provided\nYou can use '$commandline -h' to see help information for the module\n$usageModule"
	);
}
elsif ( !( -s $filelist ) ) {
	die(
"Input file (-i) didn't exist or size equal 0\nYou can use '$commandline -h' to see help information for the module\n$usageModule"
	);
}
if ( !defined $targetregionfile and defined $nod ) {
	die(
"Target region file (-r) must be provided if option -nod used.\nYou can use '$commandline -h' to see help information for the module\n$usageModule"
	);
}
if ( !defined $isdepth ) { $isdepth = 0; }
if ( !defined $nod ) { $nod = 0; }
if ( !defined $d_cumul ) { $d_cumul = "0,10,30"; }
my @d_cumul_array = split(',', $d_cumul);

if ( !defined $no_batch ) { $no_batch = 0; }
if ( !defined $use_SM ) { $use_SM = 0; }
if ( !defined $usePASS ) { $usePASS = 0; }
if ( !defined $xym ) { $xym = 0; }
if ( !defined $method )  { $method  = 1; }
if ( !( defined $vcfCfgFile ) or $vcfCfgFile eq '' ) {
	$vcfCfgFile = dirname($0) . '/GATK.cfg';
}
if ( !defined $maxThreads ) { $config{'maxThreads'} = 4; }
else {
	$config{'maxThreads'} = $maxThreads;
}
if ( !defined $singleEnd ) {
	$config{'singleEnd'} = 0;
}
else {
	$config{'singleEnd'} = $singleEnd;
}
if ( !defined $caculateMethod ) {
	$config{'caculateMethod'} = 1;
}    #1 as mean 2 as median
else {
	$config{'caculateMethod'} = $caculateMethod;
}
$config{'resultDir'} = $resultDir;
if ( !defined $rePlot ) {
	$config{'rePlot'} = 0;
}
else {
	$config{'rePlot'} = $rePlot;
}
if ( !( -e $resultDir ) ) {
	if ( mkdir $resultDir ) {
	}
	else {
		die "can't make result dir. $!\n";
	}
}

my $reportHash;
${$reportHash}{'COMMAND'}   = $commandline;
${$reportHash}{'CREATTIME'} = localtime;
my $template;
my $reportFileName;

if ( defined $resultDir ) {
	open LOG, ">$resultDir/" . $config{'logFileName'}
	  or die "can't open log file. $!\n";
}

if ( $module eq "f" ) {
	pInfo( "Start fastq summary for $filelist", \@log );
	my $rResult = &fastqSummary( $filelist, \%config );
	if ( $rResult != 0 ) {
		pInfo(
"Something wrong in plotting figure. Please check the fastqSummary.rLog file!",
			\@log
		);
	}

	#report
	my $figureList1 =
	  &dir2list( $resultDir, "/fastqFigure/", "^batch_", "FIGUREBATCH" );
	my $table1 =
	  &file2table( "$resultDir/fastqResult/fastqSummary.txt", '', 1 );
	${$reportHash}{'FIGUREBASESCOREYN'} = "./fastqFigure/fastqScoreYN.png";
	${$reportHash}{'FIGUREBASESCOREYNLEGEND'} =
	  "./fastqFigure/fastqScoreYN.png.legend.png";
	${$reportHash}{'FIGUREBASENUCYN'} = "./fastqFigure/nucPercentyYN.png";
	if ( defined($figureList1) ) {
		${$reportHash}{'FIGURELOOP1'} = $figureList1;
	};
	${$reportHash}{'MAKETABLE1'}      = $table1;
	$template =
	  &build_template( dirname($0) . "/template/report_tmpl_fastq.tmpl",
		$reportHash );
	$reportFileName = "fastqReport.html";
}
elsif ( $module eq "b" ) {
	pInfo( "Start bam summary for $filelist", \@log );
	$config{'targetregionfile'} = $targetregionfile;
	$config{'gtffile'}          = $gtffile;
	$config{'isdepth'}          = $isdepth;
	$config{'nod'}              = $nod;
	$config{'d_cumul1'}         = $d_cumul_array[0];
	$config{'d_cumul2'}         = $d_cumul_array[1];
	$config{'d_cumul3'}         = $d_cumul_array[2];
	$config{'no_batch'}         = $no_batch;
	$config{'use_SM'}           = $use_SM;
	my $rResult = &bamSummary( $filelist, \%config );
	if ( $rResult != 0 ) {
		pInfo(
"Something wrong in plotting figure. Please check the bamSummary.rLog file!",
			\@log
		);
	}

	#report
	my $figureList1 =
	  &dir2list( $resultDir, "/bamFigure/", "^summary_lines", "FIGUREDIS" );
	${$reportHash}{'FIGURELOOP1'} = $figureList1;
	my $figureList2 =
	  &dir2list( $resultDir, "/bamFigure/", "^batch_", "FIGUREBATCH" );
	if ( defined($figureList2) ) {
		${$reportHash}{'FIGURELOOP2'} = $figureList2;
	};
	my $table1 = &file2table( "$resultDir/bamResult/bamSummary.txt", '', 1 );
	${$reportHash}{'MAKETABLE1'} = $table1;

	$template =
	  &build_template( dirname($0) . "/template/report_tmpl_bam.tmpl",
		$reportHash );
	$reportFileName = "bamReport.html";
}
elsif ( $module eq "v" ) {
	pInfo( "Start vcf summary for $filelist", \@log );
	$config{'vcfCfgFile'} = $vcfCfgFile;
	$config{'method'}     = $method;
	$config{'annovarDb'}  = $annovarDb;
	$config{'usePASS'}  =$usePASS;
	$config{'xym'}          = $xym;
	my $rResult = &vcfSummary( $filelist, \%config );
	if ( $rResult != 0 ) {
		pInfo(
"Something wrong in plotting figure. Please check the vcfSummary.rLog file!",
			\@log
		);
	}

	#report
	my $vcfFileName    = basename($filelist);
	my $cfgFileContent = &file2text($vcfCfgFile);
	my $table1 =
	  &file2table( "$resultDir/vcfResult/$vcfFileName.SampleNumber.txt", '',
		1 );
	my $table2 =
	  &file2table( "$resultDir/vcfResult/$vcfFileName.Method$method.txt",
		[ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 ], 1 );
	my $table3 =
	  &file2table( "$resultDir/vcfResult/$vcfFileName.snpCount.txt", '', 1 );
	my $table4;
	if (
		-e "$resultDir/vcfAnnovarResult/$vcfFileName.pass.avinput.annovar.countTable.txt"
	  )
	{
		$table4 = &file2table(
"$resultDir/vcfAnnovarResult/$vcfFileName.pass.avinput.annovar.countTable.txt",
			'', 1
		);
	}
	my $table5 =
	  &file2table( "$resultDir/vcfResult/$vcfFileName.sexCheck.txt", '', 1 );
	my $table6 =
	  &file2table( "$resultDir/vcfResult/$vcfFileName.snpNucCount.txt", '', 1 );
	my $figureList1 =
	  &dir2list( $resultDir, "/vcfFigure/", "scoreCompare", "FIGURE2" );
	my $figureList2 =
	  &dir2list( $resultDir, "/vcfAnnovarResult/", ".variant_function", "FS2",
		1 );

	${$reportHash}{'MAKETABLE1'}        = $table1;
	${$reportHash}{'MAKETABLE2'}        = $table2;
	${$reportHash}{'MAKETABLE3'}        = $table3;
	${$reportHash}{'MAKETABLE4'}        = $table4;
	${$reportHash}{'MAKETABLE5'}        = $table5;
	${$reportHash}{'MAKETABLE6'}        = $table6;
	${$reportHash}{'FILTERFILE'}        = $vcfCfgFile;
	${$reportHash}{'FILTERFILECONTENT'} = $cfgFileContent;
	${$reportHash}{'FIG'} = "./vcfFigure/$vcfFileName.Method$method.txt.png";
	${$reportHash}{'FIGNUMBER'} = "./vcfFigure/$vcfFileName.sampleNumber.png";
	${$reportHash}{'FIGURE1'}   = $figureList1;

	if ( defined($figureList2) ) {
		${$reportHash}{'FL1'} = $figureList2;
	}

	$template =
	  &build_template( dirname($0) . "/template/report_tmpl_vcf.tmpl",
		$reportHash );
	$reportFileName = "vcfReport.html";
}
else {
	die
"Module (-m) is required and must be f (fastq), b (bam), and v (vcf)\n$usageModule";
}

open REPORT, ">$resultDir/$reportFileName" or die "can't open $!\n";
print REPORT $template->output;

pInfo( "Success!", \@log );
foreach my $log (@log) {
	print LOG $log;
}

sub pInfo {
	my $s = shift;
	print "[", scalar(localtime), "] $s\n";
	push @{ $_[1] }, "[" . scalar(localtime) . "] $s\n";
}
