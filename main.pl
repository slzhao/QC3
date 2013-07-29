use strict;
use warnings;
#use forks;
use threads;
#use threads::shared;
use File::Basename;
use HTML::Template;
use Getopt::Long;

use FindBin;
use lib $FindBin::Bin;
use source::makeReport;
use source::fastqSummary;
use source::bamSummary;
use source::vcfSummary;

my %config;
$config{'RBin'}        = "R";       #where the R bin file is
$config{'annovarBin'}  = "annotate_variation.pl";		#where the ANNOVAR bin file is
$config{'samtoolsBin'} = "samtools";	#where the SAMtools bin file is
$config{'logFileName'} = "qcLog.txt";	#name of the log file

my $commandline = "perl $0 "
  . join( " ", @ARGV )
  ;    #Output you input command in the report to reproduce you result
my (
	$module,           $filelist, $resultDir, $pairend,
	$targetregionfile, $gtffile,  $isdepth,   $caculateMethod,$cfgFile,
	$method,           $maxThreads,$annovarDb
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
);
if ( !defined $pairend )    { $pairend              = 0; }
if ( !defined $isdepth )    { $isdepth              = 0; }
if ( !defined $maxThreads ) { $config{'maxThreads'} = 4; }
else {
	$config{'maxThreads'} = $maxThreads;
}
if ( !defined $caculateMethod ) { $config{'caculateMethod'} = 1; }
else {
	$config{'caculateMethod'} = $caculateMethod;
}
$config{'resultDir'} = $resultDir;
if ( !( -e $resultDir ) ) { mkdir $resultDir; }

my $reportHash;
${$reportHash}{'COMMAND'}   = $commandline;
${$reportHash}{'CREATTIME'} = localtime;
my $template;
my $reportFileName;

open $config{'log'}, ">$resultDir/".$config{'logFileName'} or die "can't open log file. $!\n";

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
#	${$reportHash}{'FIGUREBASESCORE'} = "./fastqfigure/fastqScore.png";
#	${$reportHash}{'FIGUREBASENUC'}   = "./fastqfigure/nucPercenty.png";
#	${$reportHash}{'FIGUREBASESCOREN'} = "./fastqfigure/fastqScoreN.png";
	${$reportHash}{'FIGUREBASESCOREYN'} = "./fastqfigure/fastqScoreYN.png";
	${$reportHash}{'FIGUREBASESCOREYNLEGEND'} = "./fastqfigure/fastqScoreYN.png.legend.png";
#	${$reportHash}{'FIGUREBASENUCN'}   = "./fastqfigure/nucPercentyN.png";
	${$reportHash}{'FIGUREBASENUCYN'}   = "./fastqfigure/nucPercentyYN.png";
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
	$config{'annovarDb'}   = $annovarDb; #"/scratch/cqs/shengq1/references/annovar/humandb/"
	my $rResult = &vcfSummary( $filelist, \%config );
	if ( $rResult != 0 ) {
		print(
"Something wrong in plotting figure. Please check the vcfSummary.rLog file!\n"
		);
	}

	#report
	my $vcfFileName = basename($filelist);
	my $table1 =
	  &file2table("$resultDir/vcfResult/$vcfFileName.SampleNumber.txt",'',1);
	my $table2 =
	  &file2table( "$resultDir/vcfResult/$vcfFileName.Method$method.txt",
		[ 1, 2, 3,4,5,6,7, 8,9,10, 11 ],1 );
	my $table3 =
	  &file2table("$resultDir/vcfResult/$vcfFileName.snpCount.txt",'',1);
	my $figureList1 =
	  &dir2list( $resultDir, "/vcfFigure/", "scoreCompare", "FIGURE2" );
	my $figureList2 =
	  &dir2list( $resultDir, "/vcfAnnovarResult/", ".variant_function", "FS2",
		1 );
	${$reportHash}{'TABLE1'} = $table1;
#	${$reportHash}{'TABLE2'} = $table2;
	${$reportHash}{'MAKETABLE1'} = $table3;
	${$reportHash}{'MAKETABLE2'} = $table2;
	${$reportHash}{'MAKETABLE3'} = $table1;
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
