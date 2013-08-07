#!/usr/bin/perl
use strict;
use warnings;
use threads;
use File::Basename;
use HTML::Template;
use Getopt::Long;

use FindBin;
use lib $FindBin::Bin;
use source::makeReport;
use source::fastqSummary;
use source::bamSummary;
use source::vcfSummary;

my $version="RC 1.3";
my %config;
$config{'RBin'}        = "R";       #where the R bin file is
$config{'annovarBin'}  = "table_annovar.pl";		#where the ANNOVAR bin file is
$config{'annovarConvert'}  = "convert2annovar.pl";		#where the ANNOVAR convert bin file is
$config{'annovarOption'}  = "-buildver hg19 -protocol refGene,snp137 -operation g,f --remove"; #other options for ANNOVAR
$config{'samtoolsBin'} = "samtools";	#where the SAMtools bin file is
$config{'logFileName'} = "qc3Log.txt";	#name of the log file

my $commandline = "perl $0 "
  . join( " ", @ARGV )
  ;    #Output you input command in the report to reproduce you result
my (
	$module,           $filelist, $resultDir, $pairend,
	$targetregionfile, $gtffile,  $isdepth,   $caculateMethod,$cfgFile,
	$method,           $maxThreads,$annovarDb,$rePlot
) = ();
GetOptions(
	"t=i" => \$maxThreads,
	"m=s" => \$module,
	"i=s" => \$filelist,
	"o=s" => \$resultDir,
	"p"   => \$pairend,

	"r=s" => \$targetregionfile,
	"g=s" => \$gtffile,
	"d"   => \$isdepth,
	"cm=i"   => \$caculateMethod,

	"c:s" => \$cfgFile,
	"s=i" => \$method,
	"a=s" => \$annovarDb,
	
	"rp"   => \$rePlot,
);
if ( !defined $module) {die ("Module (-m) must be specified\n");}
if ( !defined $filelist or !defined $resultDir) {die ("Input file (-i) and Output directory (-o) must be provided\n");}

if ( !defined $pairend )    { $pairend              = 0; }
if ( !defined $isdepth )    { $isdepth              = 0; }
if ( !defined $method ) { $method = 1; }
if ( !( defined $cfgFile ) or $cfgFile eq '' ) {
		$cfgFile = dirname($0) . '/GATK.cfg';
}
if ( !defined $maxThreads ) { $config{'maxThreads'} = 4; }
else {
	$config{'maxThreads'} = $maxThreads;
}
if ( !defined $caculateMethod ) { $config{'caculateMethod'} = 1; } #1 as mean 2 as median
else {
	$config{'caculateMethod'} = $caculateMethod;
}
$config{'resultDir'} = $resultDir;
if ( !defined $rePlot )    { 
	$config{'rePlot'}          = 0; } else {
	$config{'rePlot'}          = $rePlot;
}
if (!( -e $resultDir ) ) { 
	if (mkdir $resultDir) {
	} else {	
		die "can't make result dir. $!\n";
	}
}

my $reportHash;
${$reportHash}{'COMMAND'}   = $commandline;
${$reportHash}{'CREATTIME'} = localtime;
my $template;
my $reportFileName;

if (defined $resultDir) {
	open $config{'log'}, ">$resultDir/".$config{'logFileName'} or die "can't open log file. $!\n";
}

if ( $module eq "f" ) {
	pInfo("Start fastq summary for $filelist",$config{'log'});
	$config{'pairEnd'} = $pairend;
	my $rResult = &fasqSummary( $filelist, \%config );
	if ( $rResult != 0 ) {
		print(
"Something wrong in plotting figure. Please check the fastqSummary.rLog file!\n"
		);
	}

	#report
	my $figureList1 =
	  &dir2list( $resultDir, "/fastqFigure/", "^batch_", "FIGUREBATCH" );
	my $table1 =
	  &file2table("$resultDir/fastqResult/fastqSummary.txt",'',1);
	${$reportHash}{'FIGUREBASESCOREYN'} = "./fastqFigure/fastqScoreYN.png";
	${$reportHash}{'FIGUREBASESCOREYNLEGEND'} = "./fastqFigure/fastqScoreYN.png.legend.png";
	${$reportHash}{'FIGUREBASENUCYN'}   = "./fastqFigure/nucPercentyYN.png";
	${$reportHash}{'FIGURELOOP1'}     = $figureList1;
	${$reportHash}{'MAKETABLE1'} = $table1;
	$template =
	  &build_template( dirname($0) . "/template/report_tmpl_fastq.tmpl",
		$reportHash );
	$reportFileName = "fastqReport.html";
}
elsif ( $module eq "b" ) {
	pInfo("Start bam summary for $filelist",$config{'log'});
	$config{'targetregionfile'} = $targetregionfile;
	$config{'gtffile'}          = $gtffile;
	$config{'isdepth'}          = $isdepth;
	my $rResult = &bamSummary( $filelist, \%config );
	if ( $rResult != 0 ) {
		print(
"Something wrong in plotting figure. Please check the bamSummary.rLog file!\n"
		);
	}

	#report
	my $figureList1 =
	  &dir2list( $resultDir, "/bamFigure/", "^summary_lines", "FIGUREDIS" );
	${$reportHash}{'FIGURELOOP1'} = $figureList1;
	my $figureList2 =
	  &dir2list( $resultDir, "/bamFigure/", "^batch_", "FIGUREBATCH" );
	${$reportHash}{'FIGURELOOP2'} = $figureList2;
	my $table1 =
	  &file2table("$resultDir/bamResult/bamSummary.txt",'',1);
	${$reportHash}{'MAKETABLE1'} = $table1;
	
	$template =
	  &build_template( dirname($0) . "/template/report_tmpl_bam.tmpl",
		$reportHash );
	$reportFileName = "bamReport.html";
}
elsif ( $module eq "v" ) {
	pInfo("Start vcf summary for $filelist",$config{'log'});
	$config{'cfgFile'} = $cfgFile;
	$config{'method'}  = $method;
	$config{'annovarDb'}   = $annovarDb;
	my $rResult = &vcfSummary( $filelist, \%config );
	if ( $rResult != 0 ) {
		print(
"Something wrong in plotting figure. Please check the vcfSummary.rLog file!\n"
		);
	}

	#report
	my $vcfFileName = basename($filelist);
	my $cfgFileContent=&file2text($cfgFile);
	my $table1 =
	  &file2table("$resultDir/vcfResult/$vcfFileName.SampleNumber.txt",'',1);
	my $table2 =
	  &file2table( "$resultDir/vcfResult/$vcfFileName.Method$method.txt",
		[ 1, 2, 3,4,5,6,7, 8,9,10, 11 ],1 );
	my $table3 =
	  &file2table("$resultDir/vcfResult/$vcfFileName.snpCount.txt",'',1);
	my $table4;
	 if (-e "$resultDir/vcfAnnovarResult/$vcfFileName.pass.avinput.annovar.countTable.txt") {
	 		$table4 =
	  &file2table("$resultDir/vcfAnnovarResult/$vcfFileName.pass.avinput.annovar.countTable.txt",'',1);
	 }
	my $figureList1 =
	  &dir2list( $resultDir, "/vcfFigure/", "scoreCompare", "FIGURE2" );
	my $figureList2 =
	  &dir2list( $resultDir, "/vcfAnnovarResult/", ".variant_function", "FS2",
		1 );
#	${$reportHash}{'TABLE1'} = $table1;
	${$reportHash}{'MAKETABLE1'} = $table1;
	${$reportHash}{'MAKETABLE2'} = $table2;
	${$reportHash}{'MAKETABLE3'} = $table3;
	${$reportHash}{'MAKETABLE4'} = $table4;
	${$reportHash}{'FILTERFILE'}   = $cfgFile;
	${$reportHash}{'FILTERFILECONTENT'}   = $cfgFileContent;
	${$reportHash}{'FIG'}    = "./vcfFigure/$vcfFileName.Method$method.txt.png";
	${$reportHash}{'FIGNUMBER'} = "./vcfFigure/$vcfFileName.sampleNumber.png";
	${$reportHash}{'FIGURE1'}   = $figureList1;
	if (defined ($figureList2)) {
		${$reportHash}{'FL1'}       = $figureList2;
	}
	
	$template =
	  &build_template( dirname($0) . "/template/report_tmpl_vcf.tmpl",
		$reportHash );
	$reportFileName = "vcfReport.html";
}
else {
	die "Module (-m) must be f (fastq), b (bam), and v (vcf)\n";
}

open REPORT, ">$resultDir/$reportFileName" or die "can't open $!\n";
print REPORT $template->output;

pInfo ("Success!",$config{'log'});

sub pInfo {
	my $s = shift;
	my $logFile = shift;
	print "[", scalar(localtime), "] $s\n";
	print $logFile "[", scalar(localtime), "] $s\n";
}
