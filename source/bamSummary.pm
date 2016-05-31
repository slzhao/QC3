package source::bamSummary;

use strict;
use warnings;

#use forks;
use threads;
use threads::shared;
use Exporter;
use File::Basename;

#use List::Util qw(sum);

our @ISA    = qw(Exporter);
our @EXPORT = qw(bamSummary);

1;

my $rSourceLocation = '/source/bam_plot_byPerl.R';
my $current : shared;
my @resultOut : shared;

my $logRef = \@main::log;

sub bamSummary {

	# Get the metrics from a bam file
	my ( $filelist, $config ) = @_;
	my $targetregionfile = $config->{'targetregionfile'};
	my $gtffile          = $config->{'gtffile'};
	my $resultDir        = $config->{'resultDir'};
	my $isdepth          = $config->{'isdepth'};
	my $nod              = $config->{'nod'};
	my $d_cumul1         = $config->{'d_cumul1'};
	my $d_cumul2         = $config->{'d_cumul2'};
	my $d_cumul3         = $config->{'d_cumul3'};
	my $no_batch         = $config->{'no_batch'};
	my $use_SM           = $config->{'use_SM'};
	my $RBin             = $config->{'RBin'};
	my $samtoolsBin      = $config->{'samtoolsBin'};

	my $logFile        = $config->{'log'};
	my $caculateMethod = $config->{'caculateMethod'};    #1 as mean 2 as median

	my $maxThreads = $config->{'maxThreads'};
	my $rePlot     = $config->{'rePlot'};
	my $showHelp   = $config->{'showHelp'};

	my $usage =
"Program: qc3.pl (a quality control tool for DNA sequencing data in raw data, alignment, and variant calling stages)
Version: $main::version

Module:  bam QC
Usage:   perl qc3.pl -m b -i bamFileList -o outputDirectory [-t threads] [other options]

Options:

	-i	input filelist        Required. A list file for bam files to be analyzed.
	-o	output directory      Required. Output directory for QC result. If the directory doesn't exist, it would be created.
	-t	threads               Optional. Threads used in analysis. The default value is 4.
	-r	region file           Optional. A targetregion file. At least one targetregion file or gtf file should be provided.
	-g	region file           Optional. A gtf file. At least one targetregion file or gtf file should be provided.
	-cm	method	              Optional. Calculation method for data summary, should be 1 or 2. Method 1 means mean and method 2 means median. The default value is 1.
	-d	                      Optional. The depth in on-/off-target regions will be calculated. QC3 will not calculate depth by default, because it may take a long time.
	-nod	no off-target depth Optional. The depth in off-target regions will not be calculated to save time.
	-no_batch                 Optional. The batch effect part of figures will not be plotted (considering BAM file not informative).
	-use_SM                   Optional. Extract SM field from bam file and use it as sample name.
	-d_cumul                  Optional. Depth values for cumulative distribution computation. Put 3 values if format: -d_cumul 0,10,30 (default values).

	-h	                      Optional. Show this information.

For more information, please refer to the readme file in QC3 directory. Or visit the QC3 website at https://github.com/slzhao/QC3

";

	die "\n$usage" if ($showHelp);
	die
"Can't find targetregion file(-r) or gtf file(-g). You need to specify at least one of them\n$usage"
	  if ( ( defined($targetregionfile) and !-e $targetregionfile )
		or ( defined($gtffile) and !-e $gtffile )
		or ( !defined($targetregionfile) and !defined($gtffile) ) );
	if ( !( -e $resultDir . '/bamResult/' ) ) {
		mkdir $resultDir . '/bamResult/';
	}
	my $outputFile = $resultDir . '/bamResult/bamSummary.txt';
	my $methodText;
	if ( $caculateMethod == 1 ) { $methodText = "Mean"; }
	else {
		$methodText = "Median";
	}

	$| = 1;
	if ($rePlot) {    #re plot by R only, comment below codes
		pInfo(
"Parameter -rp was used. The program will use the result files in output directory to re-generate the report",
			$logRef
		);
	}
	else {
		if ( $isdepth == 1 ) {
			my $off_on;
			if ( $nod == 1 ) { $off_on = "on"} else { $off_on = "on-/off"; $nod = 0 }
			pInfo(
"Will calculate the depth in $off_on-target regions. It will take a long time. You can turn it off without -d in your command line",
				$logRef
			);
		}
		else {
			$isdepth = 0;
		}

		my %regionDatabase;
		open( BAMFILE1, $filelist ) or die $!;
		my $bamFile1 = <BAMFILE1>;
		&initializeDatabase( $bamFile1, \%regionDatabase, $samtoolsBin );
		close(BAMFILE1);
		my $inBedSign = 1;    #1=10 in bed; 2=01 in exon; 3=11 in intron

		# 1) Load targetregionfile
		if ( defined($targetregionfile) ) {
			pInfo( "Load target region file '$targetregionfile'", $logRef );
			loadbed( $targetregionfile, \%regionDatabase );
		}
		else {
			$inBedSign = 2;
			pInfo(
"No targetregion file, will use exon regions in gtf file instead",
				$logRef
			);
		}

		# 2) Load gtf file
		if ( defined($gtffile) ) {
			pInfo( "Load gtf file '$gtffile'", $logRef );
			loadgtf( $gtffile, \%regionDatabase );
		}
		else {
			pInfo( "No gtf file", $logRef );
		}

		# 3) count total exon length
		my $totalExonLength = 0;
		if ($isdepth) {

			#			pInfo( "Count total exome length", $logFile );
			#			foreach my $chr ( keys %chrLength ) {
			#				for ( my $pos = 0 ; $pos <= $chrLength{$chr} ; $pos++ ) {
			#					if ( vec( $regionDatabase{$chr}, $pos, 2 ) == $inBedSign ) {
			#						$totalExonLength++;
			#					}
			#				}
			#			}
			#			pInfo( "Method1 Length: $totalExonLength", $logFile );
			$totalExonLength = 0;
			my $example = "";
			vec( $example, 0, 2 ) = $inBedSign;
			vec( $example, 1, 2 ) = $inBedSign;
			my ($find11) = unpack( 'h', $example );
			my $find01 = "";

			foreach my $x ( 0 .. 3 ) {
				if ( $x eq $inBedSign ) { next; }
				my $example = "";
				vec( $example, 0, 2 ) = $x;
				vec( $example, 1, 2 ) = $inBedSign;
				my ($hex) = unpack( 'h', $example );
				$find01  = $find01 . $hex;
				$example = "";
				vec( $example, 0, 2 ) = $inBedSign;
				vec( $example, 1, 2 ) = $x;
				($hex) = unpack( 'h', $example );
				$find01 = $find01 . $hex;
			}
			foreach my $chr ( keys %regionDatabase ) {
				my ($hex) = unpack( 'h*', $regionDatabase{$chr} );
				my $num1 = eval "\$hex =~ tr/[$find11]//";
				my $num2 = eval "\$hex =~ tr/[$find01]//";
				$totalExonLength = $totalExonLength + 2 * $num1 + $num2;
			}

			#			pInfo( "Method2 Length: $totalExonLength", $logFile );
		}

		# 4) read bam file
		pInfo( "Read bam file and write out", $logRef );
		open( IN,  $filelist )      or die $!;
		open( OUT, ">$outputFile" ) or die $!;

		my @fileList;    #get file list
		while ( my $f = <IN> ) {
			$f =~ s/\r|\n//g;
			my @fileLable = ( split /\t/, $f );
			if ( !-e $fileLable[0] ) {
				pInfo( "$fileLable[0] doesn't exist", $logRef );
				next;
			}
			push @fileList, $f;
		}
		close(IN);

		my @threads;     #open threads and get results
		$current = 0;
		foreach my $x ( 1 .. $maxThreads ) {
			push @threads,
			  threads->new(
				\&bamProcess,     $x,              \@fileList,
				\%regionDatabase, $inBedSign,      $isdepth,
				$samtoolsBin,     $caculateMethod, $totalExonLength, $nod,      $targetregionfile,
				$no_batch,        $use_SM,         $d_cumul1,        $d_cumul2, $d_cumul3,          $logRef
			  );
		}

		foreach my $thread (@threads) {
			my $threadNum = $thread->join;
			pInfo( "Thread $threadNum finished", $logRef );
		}
		print OUT join "\t",
		  (
			"Sample",
			"SM",
			"Instrument",
			"Run",
			"Flowcell",
			"Lane",
			"Unmapped",
			"Total Mapped",
			"On-target",
			"Off-target",
			"Off-target-intron",
			"Off-target-intergenic",
			"Off-target-mito",
			"Total Mapped($methodText MQ)",
			"On-target($methodText MQ)",
			"Off-target($methodText MQ)",
			"Off-target-intron($methodText MQ)",
			"Off-target-intergenic($methodText MQ)",
			"Off-target-mito($methodText MQ)",
			"Total Mapped($methodText InsertSize)",
			"On-target($methodText InsertSize)",
			"Off-target($methodText InsertSize)",
			"Off-target-intron($methodText InsertSize)",
			"Off-target-intergenic($methodText InsertSize)",
			"Off-target-mito($methodText InsertSize)"
		  );

		my $colNumber = $resultOut[0] =~ y/\t//;
		if (   ( $isdepth and $colNumber > 32 )
			or ( !$isdepth and $colNumber > 23 ) )
		{    #result with AS
			print OUT "\t";
			print OUT join "\t",
			  (
				"Total Mapped($methodText AS)",
				"On-target($methodText AS)",
				"Off-target($methodText AS)",
				"Off-target-intron($methodText AS)",
				"Off-target-intergenic($methodText AS)",
				"Off-target-mito($methodText AS)"
			  );
		}
		if ($isdepth) {
			print OUT "\t";
			print OUT join "\t",
			  (
				"Total Reads($methodText Depth)",
				"On-target($methodText Depth)",
				"Off-target($methodText Depth)",
				"Off-target-intron($methodText Depth)",
				"Off-target-intergenic($methodText Depth)",
				"Off-target-mito($methodText Depth)",
				"On-target percent (Depth larger than $d_cumul1)",
				"On-target percent (Depth larger than $d_cumul2)",
				"On-target percent (Depth larger than $d_cumul3)\n"
			  );
		}
		else {
			print OUT "\n";
		}
		foreach my $result (@resultOut) {
			print OUT "$result\n";
		}
		close OUT;

		#end comment by test R
	}

	#plot by R
	my $Rsource =
	  dirname($0) . "/source/rFunctions.R " . dirname($0) . $rSourceLocation;
	my $rResult = system(
"cat $Rsource | $RBin --vanilla --slave --args $resultDir $no_batch $use_SM 1>$resultDir/bamResult/bamSummary.rLog 2>$resultDir/bamResult/bamSummary.rLog"
	);
	pInfo( "Finish bam summary!", $logRef );
	return ($rResult);
}

sub getbammetric {
	my ( $in, $regionDatabaseRef, $inBedSign, $isdepth, $samtoolsBin,
		$caculateMethod, $totalExonLength, $nod, $targetregionfile, $use_SM, $d_cumul1, $d_cumul2, $d_cumul3 )
	  = @_;
	my ( $total, $ontarget, $offtarget, $ummapped, $offtargetintron,
		$offtargetintergenic, $offtargetmito )
	  = ( 0, 0, 0, 0, 0, 0, 0 );
	my ( $totalNorm, $ontargetNorm, $offtargetNorm, $ummappedNorm, $offtargetintronNorm,
        $offtargetintergenicNorm, $offtargetmitoNorm )
      = ( 0, 0, 0, 0, 0, 0, 0 );
	my (
		$totalMQflag,               $totalInsertflag,
		$ontargetMQflag,            $ontargetInsertflag,
		$offtargetMQflag,           $offtargetInsertflag,
		$offtargetintronMQflag,     $offtargetintronInsertflag,
		$offtargetintergenicMQflag, $offtargetintergenicInsertflag,
		$offtargetmitoMQflag,       $offtargetmitoInsertflag,
		$totalASflag,               $ontargetASflag,
		$offtargetASflag,           $offtargetintronASflag,
		$offtargetintergenicASflag, $offtargetmitoASflag
	) = ( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
	my (
		$totalMQ,                   $ontargetMQ,
		$offtargetMQ,               $offtargetintronMQ,
		$offtargetintergenicMQ,     $offtargetmitoMQ,
		$totalAS,                   $ontargetAS,
		$offtargetAS,               $offtargetintronAS,
		$offtargetintergenicAS,     $offtargetmitoAS,
		$totalinsert,               $ontargetinsert,
		$offtargetinsert,           $offtargetintroninsert,
		$offtargetintergenicinsert, $offtargetmitoinsert
	);
	my @fileLable = ( split /\t/, $in );
	$in = $fileLable[0];
	my $label = $in;
	if ( defined $fileLable[1] ) {
		$label = $fileLable[1];
	}
	open( BAMLINE1, "$samtoolsBin view $in|" ) or die $!;
	my $firstLine = <BAMLINE1>;
	my @headers = split( ":", ( split( "\t", $firstLine ) )[0] );
	my ( $instrument, $runNumber, $flowcell, $lane ) =
	  ( $headers[0], $headers[1], $headers[2], $headers[3] );

	my $flag_read_unmapped           = 0x0004;
	my $flag_read_mapped_proper_pair = 0x0002;

#	my $mutipleAlignInd=0;
	my $mutipleAlignKey='';
	while(<BAMLINE1>) {
		my @line   = split "\t", $_;
		my $flag   = $line[1];
		if (!($flag & $flag_read_unmapped)) { #mapped
			#find $mutipleAlignKey here
			foreach my $i (11..(scalar(@line)-1)) {
#				print $line[$i]."\n";
				if ($line[$i]=~/(NH:i:)\d+/ or $line[$i]=~/(IH:i:)\d+/) {
#					$mutipleAlignInd=$i;
					$mutipleAlignKey=$1;
#					print $mutipleAlignKey."\n";
				}
			}
			last;
		}
	}

	close(BAMLINE1);
	open( BAM, "$samtoolsBin view $in|" ) or die $!;


	while (<BAM>) {
		my @line   = split "\t", $_;
		my $flag   = $line[1];
		my $MQ     = $line[4];
		my $MQflag = 0;
		if ( $MQ =~ /^\d+$/ ) {
			$MQflag = 1;
		}
		else { print "$_\n"; }
		my $insert     = $line[8];
		my $insertflag = 0;
		if ( $insert =~ /^\d+$/
			and ( $flag & $flag_read_mapped_proper_pair ) )
		{
			$insertflag = 1;
		}
		my $chr = $line[2];
		$chr =~ s/^chr//i;
		my $pos = $line[3];

#		my $AS     = $line[11];
		my $AS     = 0;
		my $ASflag = 0;
		if (/\tAS:i:([-]*\d+)\t/ ) {
			$AS     = $1;
			$ASflag = 1;
		}

		if (   !exists( $regionDatabaseRef->{$chr} )
			or !defined vec( $regionDatabaseRef->{$chr}, $pos, 2 ) )
		{
			next;
		}

		my $mutipleAlignCount=1;
		if ($mutipleAlignKey ne "") {
			#find the $mutipleAlignCount by $mutipleAlignInd
			if (/\t$mutipleAlignKey(\d+)\t/) {
				$mutipleAlignCount=1/$1;
				#print("Find:".$mutipleAlignInd.":".$mutipleAlignKey)
			}
		}

		if ( $flag & $flag_read_unmapped ) {
			#read unmapped
#			$ummapped++;
			$ummappedNorm+=$mutipleAlignCount;
		}
		else {
#			$total++;
			$totalNorm+=$mutipleAlignCount;
			if ($MQflag) {
				&storeData( $MQ, $totalMQ, $caculateMethod );
				$totalMQflag++;
			}
			if ($ASflag) {
				&storeData( $AS, $totalAS, $caculateMethod );
				$totalASflag++;
			}
			if ($insertflag) {
				&storeData( $insert, $totalinsert, $caculateMethod );
				$totalInsertflag++;
			}

			if ( vec( $regionDatabaseRef->{$chr}, $pos, 2 ) == $inBedSign )
			{    #in bed file
#			    $ontarget++;
				$ontargetNorm+=$mutipleAlignCount;
				if ($MQflag) {
					&storeData( $MQ, $ontargetMQ, $caculateMethod );
					$ontargetMQflag++;
				}
				if ($ASflag) {
					&storeData( $AS, $ontargetAS, $caculateMethod );
					$ontargetASflag++;
				}
				if ($insertflag) {
					&storeData( $insert, $ontargetinsert, $caculateMethod );
					$ontargetInsertflag++;
				}

			}
			else {
#				$offtarget++;
				$offtargetNorm+=$mutipleAlignCount;
				if ($MQflag) {
					&storeData( $MQ, $offtargetMQ, $caculateMethod );
					$offtargetMQflag++;
				}
				if ($ASflag) {
					&storeData( $AS, $offtargetAS, $caculateMethod );
					$offtargetASflag++;
				}
				if ($insertflag) {
					&storeData( $insert, $offtargetinsert, $caculateMethod );
					$offtargetInsertflag++;
				}

				if ( vec( $regionDatabaseRef->{$chr}, $pos, 2 ) == 3 )
				{    #in intron
#				    $offtargetintron++;
					$offtargetintronNorm+=$mutipleAlignCount;
					if ($MQflag) {
						&storeData( $MQ, $offtargetintronMQ, $caculateMethod );
						$offtargetintronMQflag++;
					}
					if ($ASflag) {
						&storeData( $AS, $offtargetintronAS, $caculateMethod );
						$offtargetintronASflag++;
					}
					if ($insertflag) {
						&storeData( $insert, $offtargetintroninsert,
							$caculateMethod );
						$offtargetintronInsertflag++;
					}
				}
				elsif ( !vec( $regionDatabaseRef->{$chr}, $pos, 2 ) )
				{    #==0, then not in exon and intron, must be intergenic
#				    $offtargetintergenic++;
					$offtargetintergenicNorm+=$mutipleAlignCount;
					if ($MQflag) {
						&storeData( $MQ, $offtargetintergenicMQ,
							$caculateMethod );
						$offtargetintergenicMQflag++;
					}
					if ($ASflag) {
						&storeData( $AS, $offtargetintergenicAS,
							$caculateMethod );
						$offtargetintergenicASflag++;
					}
					if ($insertflag) {
						&storeData( $insert, $offtargetintergenicinsert,
							$caculateMethod );
						$offtargetintergenicInsertflag++;
					}
				}
				if ( $chr =~ /M/i ) {
					$offtargetmito+=$mutipleAlignCount;
					if ($MQflag) {
						&storeData( $MQ, $offtargetmitoMQ, $caculateMethod );
						$offtargetmitoMQflag++;
					}
					if ($ASflag) {
						&storeData( $AS, $offtargetmitoAS, $caculateMethod );
						$offtargetmitoASflag++;
					}
					if ($insertflag) {
						&storeData( $insert, $offtargetmitoinsert,
							$caculateMethod );
						$offtargetmitoInsertflag++;
					}
				}
			}
		}
	}
	close BAM;

	$totalMQ = findMedianMean( $totalMQ, $totalMQflag, $caculateMethod );
	$totalAS = findMedianMean( $totalAS, $totalASflag, $caculateMethod );
	$totalinsert =
	  findMedianMean( $totalinsert, $totalInsertflag, $caculateMethod );

	$ontargetMQ =
	  findMedianMean( $ontargetMQ, $ontargetMQflag, $caculateMethod );
	$ontargetAS =
	  findMedianMean( $ontargetAS, $ontargetASflag, $caculateMethod );
	$ontargetinsert =
	  findMedianMean( $ontargetinsert, $ontargetInsertflag, $caculateMethod );

	$offtargetMQ =
	  findMedianMean( $offtargetMQ, $offtargetMQflag, $caculateMethod );
	$offtargetAS =
	  findMedianMean( $offtargetAS, $offtargetASflag, $caculateMethod );
	$offtargetinsert =
	  findMedianMean( $offtargetinsert, $offtargetInsertflag, $caculateMethod );

	$offtargetintronMQ =
	  findMedianMean( $offtargetintronMQ, $offtargetintronMQflag,
		$caculateMethod );
	$offtargetintronAS =
	  findMedianMean( $offtargetintronAS, $offtargetintronASflag,
		$caculateMethod );
	$offtargetintroninsert =
	  findMedianMean( $offtargetintroninsert, $offtargetintronInsertflag,
		$caculateMethod );

	$offtargetintergenicMQ =
	  findMedianMean( $offtargetintergenicMQ, $offtargetintergenicMQflag,
		$caculateMethod );
	$offtargetintergenicAS =
	  findMedianMean( $offtargetintergenicAS, $offtargetintergenicASflag,
		$caculateMethod );
	$offtargetintergenicinsert =
	  findMedianMean( $offtargetintergenicinsert,
		$offtargetintergenicInsertflag,
		$caculateMethod );

	$offtargetmitoMQ =
	  findMedianMean( $offtargetmitoMQ, $offtargetmitoMQflag, $caculateMethod );
	$offtargetmitoAS =
	  findMedianMean( $offtargetmitoAS, $offtargetmitoASflag, $caculateMethod );
	$offtargetmitoinsert =
	  findMedianMean( $offtargetmitoinsert, $offtargetmitoInsertflag,
		$caculateMethod );

	# if need to calculate the depth
	if ($isdepth) {

		# use samtools depth to get the depth file
		#		pInfo("Getting depth of '$in'",$logFile);
		# add a variable to the samtools depth line for adapting it in case of  -nod option
		my $nod_part;
		if ($nod) { $nod_part = "-b $targetregionfile" } else { $nod_part = "" }
		open( DEPTH, "$samtoolsBin depth -d 1000000 -a $nod_part $in|" ) or die $!;
		my ( $totaldepth, $ontargetdepth, $offtargetdepth,
			$offtargetintrondepth, $offtargetintergenicdepth,
			$offtargetmitodepth );
		my (
			$totalnum,               $ontargetnum,   $ontargetnum1,
			$ontargetnum2,           $ontargetnum3,
			$offtargetnum,           $offtargetintronnum,
			$offtargetintergenicnum, $offtargetmitonum
		) = ( 0, 0, 0, 0, 0, 0, 0, 0 );    # used as denominator
		while (<DEPTH>) {
			s/\r|\n//g;
			my ( $chr, $pos, $depth ) = split "\t";
			$chr =~ s/^chr//i;
			$totalnum++;
			&storeData( $depth, $totaldepth, $caculateMethod );

			if (   !exists( $regionDatabaseRef->{$chr} )
				or !defined vec( $regionDatabaseRef->{$chr}, $pos, 2 ) )
			{
				next;
			}

			if ( vec( $regionDatabaseRef->{$chr}, $pos, 2 ) == $inBedSign ) {
				$ontargetnum++;
				&storeData( $depth, $ontargetdepth, $caculateMethod );

				if ( $depth > $d_cumul1 ) { $ontargetnum1++; }
				if ( $depth > $d_cumul2 ) { $ontargetnum2++; }
				if ( $depth > $d_cumul3 ) { $ontargetnum3++; }
			}
			else {
				$offtargetnum++;
				&storeData( $depth, $offtargetdepth, $caculateMethod );
				if ( vec( $regionDatabaseRef->{$chr}, $pos, 2 ) == 3 ) {
					$offtargetintronnum++;
					&storeData( $depth, $offtargetintrondepth,
						$caculateMethod );
				}
				elsif ( !vec( $regionDatabaseRef->{$chr}, $pos, 2 ) ) {
					$offtargetintergenicnum++;
					&storeData( $depth, $offtargetintergenicdepth,
						$caculateMethod );
				}
				if ( $chr =~ /M/i ) {
					$offtargetmitonum++;
					&storeData( $depth, $offtargetmitodepth, $caculateMethod );
				}
			}
		}
		close DEPTH;

		$totaldepth = findMedianMean( $totaldepth, $totalnum, $caculateMethod );
		$ontargetdepth =
		  findMedianMean( $ontargetdepth, $ontargetnum, $caculateMethod );
		$offtargetdepth =
		  findMedianMean( $offtargetdepth, $offtargetnum, $caculateMethod );
		$offtargetintrondepth =
		  findMedianMean( $offtargetintrondepth, $offtargetintronnum,
			$caculateMethod );
		$offtargetintergenicdepth =
		  findMedianMean( $offtargetintergenicdepth, $offtargetintergenicnum,
			$caculateMethod );
		$offtargetmitodepth =
		  findMedianMean( $offtargetmitodepth, $offtargetmitonum,
			$caculateMethod );

		my $exonRatio1 = $ontargetnum1 / $totalExonLength;
		my $exonRatio2 = $ontargetnum2 / $totalExonLength;
		my $exonRatio3 = $ontargetnum3 / $totalExonLength;

		my $SM="";
		$SM=`samtools view -H $in | grep \@RG | head -1 | sed "s/.*SM:\\([^\\t]*\\).*/\\1/" | tr -d '[:space:]'`;

		my @returnValue = (
			$label,                     $SM, $instrument,
			$runNumber,                 $flowcell,
			$lane,                      $ummappedNorm,
			$totalNorm,                     $ontargetNorm,
			$offtargetNorm,                 $offtargetintronNorm,
			$offtargetintergenicNorm,       $offtargetmitoNorm,
			$totalMQ,                   $ontargetMQ,
			$offtargetMQ,               $offtargetintronMQ,
			$offtargetintergenicMQ,     $offtargetmitoMQ,
			$totalinsert,               $ontargetinsert,
			$offtargetinsert,           $offtargetintroninsert,
			$offtargetintergenicinsert, $offtargetmitoinsert
		);

		if ( $totalAS ne "NA" ) {
			@returnValue = (
				@returnValue, $totalAS, $ontargetAS, $offtargetAS,
				$offtargetintronAS, $offtargetintergenicAS, $offtargetmitoAS
			);
		}
		@returnValue = (
			@returnValue,          $totaldepth,
			$ontargetdepth,        $offtargetdepth,
			$offtargetintrondepth, $offtargetintergenicdepth,
			$offtargetmitodepth,   $exonRatio1,
			$exonRatio2,           $exonRatio3
		);
		return ( \@returnValue );
	}
	else {

		my $SM="";
		$SM=`samtools view -H $in | grep \@RG | head -1 | sed "s/.*SM:\\([^\\t]*\\).*/\\1/" | tr -d '[:space:]'`;

		my @returnValue = (
			$label,                     $SM,	$instrument,
			$runNumber,                 $flowcell,
			$lane,                      $ummappedNorm,
			$totalNorm,                     $ontargetNorm,
			$offtargetNorm,                 $offtargetintronNorm,
			$offtargetintergenicNorm,       $offtargetmitoNorm,
			$totalMQ,                   $ontargetMQ,
			$offtargetMQ,               $offtargetintronMQ,
			$offtargetintergenicMQ,     $offtargetmitoMQ,
			$totalinsert,               $ontargetinsert,
			$offtargetinsert,           $offtargetintroninsert,
			$offtargetintergenicinsert, $offtargetmitoinsert
		);
		if ( $totalAS ne "NA" ) {
			@returnValue = (
				@returnValue, $totalAS, $ontargetAS, $offtargetAS,
				$offtargetintronAS, $offtargetintergenicAS, $offtargetmitoAS
			);
		}
		return ( \@returnValue );
	}
}

sub bamProcess {
	my (
		$threadNum,      $fileListRef,     $regionDatabaseRef,
		$inBedSign,      $isdepth,         $samtoolsBin,
		$caculateMethod, $totalExonLength, $nod, $targetregionfile, $no_batch, $use_SM,
		$d_cumul1, $d_cumul2, $d_cumul3, $logRef
	) = @_;
	pInfo( "Thread $threadNum started", $logRef );
	while (1) {
		my $file;
		my $current_temp;
		{
			lock $current;
			$current_temp = $current;
			++$current;
		}
		return $threadNum unless $current_temp < @{$fileListRef};
		$file = ${$fileListRef}[$current_temp];

		pInfo( "Thread $threadNum processing $file", $logRef );
		my ($metric) =
		  &getbammetric( $file, $regionDatabaseRef, $inBedSign, $isdepth,
			$samtoolsBin, $caculateMethod, $totalExonLength, $nod, $targetregionfile, $use_SM,
			$d_cumul1, $d_cumul2, $d_cumul3 );

		#return result here
		{
			lock @resultOut;
			$resultOut[$current_temp] = join "\t", ( @{$metric} );
		}
	}
}

sub initializeDatabase {
	my ( $in, $databaseRef, $samtoolsBin ) = @_;
	open( BAMHEAD, "$samtoolsBin view -H $in|" ) or die $!;
	while (<BAMHEAD>) {
		chomp;
		if (/^\@SQ\tSN:(\w+)\tLN:(\d+)/) {
			my $chr = $1;
			my $chrLength = $2;
			$chr =~ s/^chr//i;
			$databaseRef->{$chr} = "";
			vec( $databaseRef->{$chr}, $chrLength, 2 ) = 0;
		}
	}
}

sub loadbed {
	my ( $in, $ref ) = @_;
	open( IIN, $in ) or die $!;
	while (<IIN>) {
		s/\r|\n//g;
		my ( $chr, $start, $end ) = split "\t";
		$chr =~ s/^chr//i;
		if ( exists $ref->{$chr} ) {
			foreach ( my $i = $start + 1 ; $i <= $end ; $i++ ) {
				vec( $ref->{$chr}, $i * 2, 1 ) = 1;    #10 in bed
			}
		}
	}
	close IIN;
}

sub loadgtf {
	my ( $in, $ref ) = @_;
	open( IIN, $in ) or die $!;
	my %transcripts;
	while (<IIN>) {
		if (/^#/) {
			next;
		}
		s/\r|\n//g;
		my ( $chr, $feature, $start, $end, $attributes ) =
		  ( split( /\t/, $_ ) )[ ( 0, 2, 3, 4, 8 ) ];
		$chr =~ s/^chr//i;

		if ( $feature eq "exon" ) {
			if ( exists $ref->{$chr} ) {
				foreach my $temp ( $start .. $end ) {
					if ( !vec( $ref->{$chr}, $temp, 2 ) ) {    #not in bed file
						vec( $ref->{$chr}, $temp, 2 ) = 2;     #01 in exon
					}
				}

				my $gene_id = ( split( / \"|\"; /, $attributes ) )[1];
				$transcripts{$gene_id}->{'chr'} = $chr;
				if ( !exists( $transcripts{$gene_id}->{'start'} )
					or $start < $transcripts{$gene_id}->{'start'} )
				{
					$transcripts{$gene_id}->{'start'} = $start;
				}
				if ( !exists( $transcripts{$gene_id}->{'end'} )
					or $end > $transcripts{$gene_id}->{'end'} )
				{
					$transcripts{$gene_id}->{'end'} = $end;
				}
			}
		}
	}
	close IIN;

	#	pInfo("finish read gtf, now assign intron");
	foreach my $gene_id ( keys %transcripts ) {
		my $start = $transcripts{$gene_id}->{'start'};
		my $end   = $transcripts{$gene_id}->{'end'};
		my $chr   = $transcripts{$gene_id}->{'chr'};
		foreach my $temp ( $start .. $end ) {
			if ( !vec( $ref->{$chr}, $temp, 2 ) )
			{    #not in bed file and not in exon
				vec( $ref->{$chr}, $temp, 2 ) = 3;    #11 in intron
			}
		}
	}
	undef(%transcripts);
}

#sub pInfo {
#	my $s       = shift;
#	my $logFile = shift;
#	print "[", scalar(localtime), "] $s\n";
#	print $logFile "[", scalar(localtime), "] $s\n";
#}
sub pInfo {
	my $s = $_[0];
	print "[", scalar(localtime), "] $s\n";

	lock $_[1];
	{
		lock @{ $_[1] };
		push @{ $_[1] }, "[" . scalar(localtime) . "] $s\n";
	}
}

sub storeData {
	my $method = $_[2];
	if ( $method == 1 ) {    #mean
		$_[1] += $_[0];
	}
	else {                   #median
		my $input = $_[0];

		#		if ($input>=1000) {$input=1000;}
		${ $_[1] }{$input}++;
	}
}

sub findMedianMean {
	my $arrayLength = $_[1];
	if ( !( defined $arrayLength ) or $arrayLength == 0 ) {
		return ('NA');
	}
	my $method = $_[2];
	if ( $method == 1 ) {    #mean
		return $_[0] / $arrayLength;
	}
	else {                   #median
		my $totalCount = 0;

#use numbers in array as length, not numbers of count (including none MQ or none insert)
#		my $arryCount=0;
#		foreach my $temp (@{ $_[0] }) {
#			if (defined $temp) {
#				$arryCount += $temp;
#			}
#		}
#		$arrayLength=$arryCount;

		#		$arrayLength=sum (@{ $_[0] });

		if ( $arrayLength % 2 == 0 ) {
			my $medianPre;
			foreach my $key ( sort { $a <=> $b } keys %{ $_[0] } ) {
				$totalCount += ${ $_[0] }{$key};
				if ( $totalCount == $arrayLength / 2 ) {
					$medianPre = $key;
				}
				elsif ( $totalCount > $arrayLength / 2 ) {
					if ( defined $medianPre ) {
						return ( ( $medianPre + $key ) / 2 );
					}
					else {
						return ($key);
					}
				}
			}

		  #			for ( my $count = 0 ; $count < scalar( @{ $_[0] } ) ; $count++ ) {
		  #				if ( !defined( ${ $_[0] }[$count] ) ) { next; }
		  #				$totalCount += ${ $_[0] }[$count];
		  #				if ( $totalCount == $arrayLength / 2 ) {
		  #					$medianPre = $count;
		  #				}
		  #				elsif ( $totalCount > $arrayLength / 2 ) {
		  #					if ( defined $medianPre ) {
		  #						return ( ( $medianPre + $count ) / 2 );
		  #					}
		  #					else {
		  #						return ($count);
		  #					}
		  #				}
		  #			}
		}
		else {
			foreach my $key ( sort { $a <=> $b } keys %{ $_[0] } ) {
				$totalCount += ${ $_[0] }{$key};
				if ( $totalCount > $arrayLength / 2 ) {
					return ($key);
				}
			}

		  #			for ( my $count = 0 ; $count < scalar( @{ $_[0] } ) ; $count++ ) {
		  #				if ( !defined( ${ $_[0] }[$count] ) ) { next; }
		  #				$totalCount += ${ $_[0] }[$count];
		  #				if ( $totalCount > $arrayLength / 2 ) {
		  #					return ($count);
		  #				}
		  #			}
		}
	}
}
