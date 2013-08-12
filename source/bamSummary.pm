package source::bamSummary;

use strict;
use warnings;

#use forks;
use threads;

#use threads::shared;
use Exporter;
use File::Basename;

#use List::Util qw(sum);

our @ISA    = qw(Exporter);
our @EXPORT = qw(bamSummary);

1;

my $rSourceLocation = '/source/bam_plot_byPerl.R';

sub bamSummary {

	# Get the metrics from a bam file
	my ( $filelist, $config ) = @_;
	my $targetregionfile = $config->{'targetregionfile'};
	my $gtffile          = $config->{'gtffile'};
	my $resultDir        = $config->{'resultDir'};
	my $isdepth          = $config->{'isdepth'};
	my $RBin             = $config->{'RBin'};
	my $samtoolsBin      = $config->{'samtoolsBin'};

	my $logFile        = $config->{'log'};
	my $caculateMethod = $config->{'caculateMethod'};    #1 as mean 2 as median

	my $maxThreads = $config->{'maxThreads'};
	my $rePlot     = $config->{'rePlot'};
	my $showHelp     = $config->{'showHelp'};

	my $usage =
"Program: qc3.pl (a quality control tool for DNA sequencing data in raw data, alignment, and variant calling stages)
Version: $main::version

Module:  bam QC
Usage:   perl qc3.pl -m b -i bamFileList -o outputDirectory [-t threads] [other options]

Options:

	-i	input filelist      Required. A list file for bam files to be analyzed.
	-o	output directory    Required. Output directory for QC result. If the directory doesn't exist, it would be created.
	-t	threads             Optional. Threads used in analysis. The default value is 4.
	-r	region file         Optional. A targetregion file. At least one targetregion file or gtf file should be provided.
	-g	region file         Optional. A gtf file. At least one targetregion file or gtf file should be provided.
	-cm	method	            Optional. Calculation method for data summary, should be 1 or 2. Method 1 means mean and method 2 means median. The default value is 1.
	-d	                    Optional. whether the depth in on-/off-target regions will be calculated, It will not be calculated by default, -d = will be calculated.
	
	-h	                    Optional. Show this information.
		
For more information, please refer to the readme file in QC3 directory. Or visit the QC3 website at https://github.com/slzhao/QC3

";
	
	die "$!\n$usage"
	  if ($showHelp or ( defined($targetregionfile)  and !-e $targetregionfile )
		or ( defined($gtffile)           and !-e $gtffile )
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

	}
	else {
		if ( $isdepth == 1 ) {
			pInfo(
"Will calculate the depth in on-/off-target regions. It will take a long time. You can turn it off without -d in your command line",
				$logFile
			);
		}
		else {
			$isdepth = 0;
		}

		my %regionDatabase;
		open( BAMFILE1, $filelist ) or die $!;
		my $bamFile1 = <BAMFILE1>;
		&initializeDatabase( $bamFile1, \%regionDatabase, $samtoolsBin);
		close(BAMFILE1);
		my $inBedSign = 1;    #1=10 in bed; 2=01 in exon; 3=11 in intron

		# 1) Load targetregionfile
		if ( defined($targetregionfile) ) {
			pInfo( "Load target region file '$targetregionfile'", $logFile );
			loadbed( $targetregionfile, \%regionDatabase );
		}
		else {
			$inBedSign = 2;
			pInfo(
"No targetregion file, will use exon regions in gtf file instead",
				$logFile
			);
		}

		# 2) Load gtf file
		if ( defined($gtffile) ) {
			pInfo( "Load gtf file '$gtffile'", $logFile );
			loadgtf( $gtffile, \%regionDatabase );
		}
		else {
			pInfo( "No gtf file", $logFile );
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
		pInfo( "Read bam file and write out", $logFile );
		open( IN,  $filelist )      or die $!;
		open( OUT, ">$outputFile" ) or die $!;

		print OUT join "\t",
		  (
			"Sample",
			"Instrument",
			"Run",
			"Flowcell",
			"Lane",
			"Total Reads",
			"On-target",
			"Off-target",
			"Unmapped",
			"Off-target-intron",
			"Off-target-intergenic",
			"Off-target-mito",
			"Total Reads($methodText MQ)",
			"On-target($methodText MQ)",
			"Off-target($methodText MQ)",
			"Off-target-intron($methodText MQ)",
			"Off-target-intergenic($methodText MQ)",
			"Off-target-mito($methodText MQ)",
			"Total Reads($methodText InsertSize)",
			"On-target($methodText InsertSize)",
			"Off-target($methodText InsertSize)",
			"Off-target-intron($methodText InsertSize)",
			"Off-target-intergenic($methodText InsertSize)",
			"Off-target-mito($methodText InsertSize)"
		  );
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
				"On-target percent (Depth larger than 0)",
				"On-target percent (Depth larger than 10)",
				"On-target percent (Depth larger than 30)\n"
			  );
		}
		else {
			print OUT "\n";
		}
		while ( my $f = <IN> ) {
			$f =~ s/\r|\n//g;
			if ( !-e $f ) {
				pInfo( "$f doesn't exist", $logFile );
				next;
			}

			#		#single threads
			#		pInfo( "Processing $f ", $logFile );
			#		my @metric =
			#		  &getbammetric( $f, \%regionDatabase, $inBedSign, $isdepth,
			#			$samtoolsBin, $caculateMethod );
			#		print OUT join "\t", (@metric);
			#		print OUT "\n";

			#multi threads
			if ( scalar( threads->list() ) < $maxThreads ) {
				pInfo( "Processing $f ", $logFile );
				my ($t) = threads->new(
					\&getbammetric,  $f,         \%regionDatabase,
					$inBedSign,      $isdepth,   $samtoolsBin,
					$caculateMethod, $totalExonLength
				);
			}
			else {
				foreach my $thread ( threads->list() ) {
					my @metric = $thread->join;
					print OUT join "\t", (@metric);
					print OUT "\n";
					pInfo( "Processing $f ", $logFile );
					my ($t) = threads->new(
						\&getbammetric,  $f,         \%regionDatabase,
						$inBedSign,      $isdepth,   $samtoolsBin,
						$caculateMethod, $totalExonLength
					);
					last;
				}
			}
		}

		#join all left threads
		foreach my $thread ( threads->list() ) {
			my @metric = $thread->join;
			print OUT join "\t", (@metric);
			print OUT "\n";
		}
		close OUT;

		#end comment by test R
	}

	#plot by R
	my $Rsource =
	  dirname($0) . "/source/rFunctions.R " . dirname($0) . $rSourceLocation;
	my $rResult = system(
"cat $Rsource | $RBin --vanilla --slave --args $resultDir 1>$resultDir/bamResult/bamSummary.rLog 2>$resultDir/bamResult/bamSummary.rLog"
	);
	pInfo( "Finish bam summary!", $logFile );
	return ($rResult);
}

sub getbammetric {
	my ( $in, $regionDatabaseRef, $inBedSign, $isdepth, $samtoolsBin,
		$caculateMethod, $totalExonLength )
	  = @_;
	my ( $total, $ontarget, $offtarget, $ummapped, $offtargetintron,
		$offtargetintergenic, $offtargetmito )
	  = ( 0, 0, 0, 0, 0, 0, 0 );
	my (
		$totalMQflag,               $totalInsertflag,
		$ontargetMQflag,            $ontargetInsertflag,
		$offtargetMQflag,           $offtargetInsertflag,
		$offtargetintronMQflag,     $offtargetintronInsertflag,
		$offtargetintergenicMQflag, $offtargetintergenicInsertflag,
		$offtargetmitoMQflag,       $offtargetmitoInsertflag
	) = ( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
	my (
		$totalMQ,                   $ontargetMQ,
		$offtargetMQ,               $offtargetintronMQ,
		$offtargetintergenicMQ,     $offtargetmitoMQ,
		$totalinsert,               $ontargetinsert,
		$offtargetinsert,           $offtargetintroninsert,
		$offtargetintergenicinsert, $offtargetmitoinsert
	);
	open( BAMLINE1, "$samtoolsBin view $in|head -1|" ) or die $!;
	my $firstLine = <BAMLINE1>;
	my @headers = split( ":", ( split( "\t", $firstLine ) )[0] );
	my ( $instrument, $runNumber, $flowcell, $lane ) =
	  ( $headers[0], $headers[1], $headers[2], $headers[3] );
	close(BAMLINE1);
	open( BAM, "$samtoolsBin view $in|" ) or die $!;
	my $flag_read_unmapped           = 0x0004;
	my $flag_read_mapped_proper_pair = 0x0002;

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
			and (  $flag & $flag_read_mapped_proper_pair ) )
		{
			$insertflag = 1;
		}
		my $chr = $line[2];
		my $pos = $line[3];

		if ( $flag & $flag_read_unmapped ) {

			#read unmapped
			$ummapped++;
		}
		else {
			$total++;
			if ($MQflag) {
				&storeData( $MQ, $totalMQ, $caculateMethod );
				$totalMQflag++;
			}
			if ($insertflag) {
				&storeData( $insert, $totalinsert, $caculateMethod );
				$totalInsertflag++;
			}

			if ( vec( $regionDatabaseRef->{$chr}, $pos, 2 ) == $inBedSign )
			{    #in bed file
				$ontarget++;
				if ($MQflag) {
					&storeData( $MQ, $ontargetMQ, $caculateMethod );
					$ontargetMQflag++;
				}
				if ($insertflag) {
					&storeData( $insert, $ontargetinsert, $caculateMethod );
					$ontargetInsertflag++;
				}

			}
			else {
				$offtarget++;
				if ($MQflag) {
					&storeData( $MQ, $offtargetMQ, $caculateMethod );
					$offtargetMQflag++;
				}
				if ($insertflag) {
					&storeData( $insert, $offtargetinsert, $caculateMethod );
					$offtargetInsertflag++;
				}

				if ( vec( $regionDatabaseRef->{$chr}, $pos, 2 ) == 3 )
				{    #in intron
					$offtargetintron++;
					if ($MQflag) {
						&storeData( $MQ, $offtargetintronMQ, $caculateMethod );
						$offtargetintronMQflag++;
					}
					if ($insertflag) {
						&storeData( $insert, $offtargetintroninsert,
							$caculateMethod );
						$offtargetintronInsertflag++;
					}
				}
				elsif ( !vec( $regionDatabaseRef->{$chr}, $pos, 2 ) )
				{    #==0, then not in exon and intron, must be intergenic
					$offtargetintergenic++;
					if ($MQflag) {
						&storeData( $MQ, $offtargetintergenicMQ,
							$caculateMethod );
						$offtargetintergenicMQflag++;
					}
					if ($insertflag) {
						&storeData( $insert, $offtargetintergenicinsert,
							$caculateMethod );
						$offtargetintergenicInsertflag++;
					}
				}
				if ( $chr =~ /M/i ) {
					$offtargetmito++;
					if ($MQflag) {
						&storeData( $MQ, $offtargetmitoMQ, $caculateMethod );
						$offtargetmitoMQflag++;
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
	$totalinsert =
	  findMedianMean( $totalinsert, $totalInsertflag, $caculateMethod );

	$ontargetMQ =
	  findMedianMean( $ontargetMQ, $ontargetMQflag, $caculateMethod );
	$ontargetinsert =
	  findMedianMean( $ontargetinsert, $ontargetInsertflag, $caculateMethod );

	$offtargetMQ =
	  findMedianMean( $offtargetMQ, $offtargetMQflag, $caculateMethod );
	$offtargetinsert =
	  findMedianMean( $offtargetinsert, $offtargetInsertflag, $caculateMethod );

	$offtargetintronMQ =
	  findMedianMean( $offtargetintronMQ, $offtargetintronMQflag,
		$caculateMethod );
	$offtargetintroninsert =
	  findMedianMean( $offtargetintroninsert, $offtargetintronInsertflag,
		$caculateMethod );

	$offtargetintergenicMQ =
	  findMedianMean( $offtargetintergenicMQ, $offtargetintergenicMQflag,
		$caculateMethod );
	$offtargetintergenicinsert =
	  findMedianMean( $offtargetintergenicinsert,
		$offtargetintergenicInsertflag,
		$caculateMethod );

	$offtargetmitoMQ =
	  findMedianMean( $offtargetmitoMQ, $offtargetmitoMQflag, $caculateMethod );
	$offtargetmitoinsert =
	  findMedianMean( $offtargetmitoinsert, $offtargetmitoInsertflag,
		$caculateMethod );

	# if need to calculate the depth
	if ($isdepth) {

		# use samtools depth to get the depth file
		#		pInfo("Getting depth of '$in'",$logFile);
		open( DEPTH, "$samtoolsBin depth $in|" ) or die $!;
		my ( $totaldepth, $ontargetdepth, $offtargetdepth,
			$offtargetintrondepth, $offtargetintergenicdepth,
			$offtargetmitodepth );
		my (
			$totalnum,               $ontargetnum,
			$ontargetD10num,         $ontargetD30num,
			$offtargetnum,           $offtargetintronnum,
			$offtargetintergenicnum, $offtargetmitonum
		) = ( 0, 0, 0, 0, 0, 0, 0, 0 );    # used as denominator
		while (<DEPTH>) {
			s/\r|\n//g;
			my ( $chr, $pos, $depth ) = split "\t";
			if ( !$depth ) { next; }
			$totalnum++;
			&storeData( $depth, $totaldepth, $caculateMethod );

			if ( vec( $regionDatabaseRef->{$chr}, $pos, 2 ) == $inBedSign ) {
				$ontargetnum++;
				&storeData( $depth, $ontargetdepth, $caculateMethod );

				if ( $depth > 10 ) { $ontargetD10num++; }
				if ( $depth > 30 ) { $ontargetD30num++; }
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

		my $exonRatio1 = $ontargetnum / $totalExonLength;
		my $exonRatio2 = $ontargetD10num / $totalExonLength;
		my $exonRatio3 = $ontargetD30num / $totalExonLength;

		return (
			$in,                        $instrument,
			$runNumber,                 $flowcell,
			$lane,                      $total,
			$ontarget,                  $offtarget,
			$ummapped,                  $offtargetintron,
			$offtargetintergenic,       $offtargetmito,
			$totalMQ,                   $ontargetMQ,
			$offtargetMQ,               $offtargetintronMQ,
			$offtargetintergenicMQ,     $offtargetmitoMQ,
			$totalinsert,               $ontargetinsert,
			$offtargetinsert,           $offtargetintroninsert,
			$offtargetintergenicinsert, $offtargetmitoinsert,
			$totaldepth,                $ontargetdepth,
			$offtargetdepth,            $offtargetintrondepth,
			$offtargetintergenicdepth,  $offtargetmitodepth,
			$exonRatio1,                $exonRatio2,
			$exonRatio3
		);
	}
	else {
		return (
			$in,                        $instrument,
			$runNumber,                 $flowcell,
			$lane,                      $total,
			$ontarget,                  $offtarget,
			$ummapped,                  $offtargetintron,
			$offtargetintergenic,       $offtargetmito,
			$totalMQ,                   $ontargetMQ,
			$offtargetMQ,               $offtargetintronMQ,
			$offtargetintergenicMQ,     $offtargetmitoMQ,
			$totalinsert,               $ontargetinsert,
			$offtargetinsert,           $offtargetintroninsert,
			$offtargetintergenicinsert, $offtargetmitoinsert
		);
	}
}

sub initializeDatabase {
	my ( $in, $databaseRef, $samtoolsBin ) = @_;
	open( BAMHEAD, "$samtoolsBin view -H $in|" ) or die $!;
	while (<BAMHEAD>) {
		chomp;
		if (/^\@SQ\tSN:(\w+)\tLN:(\d+)/) {
			my $chr = $1;
			$databaseRef->{$chr} = "";
			my $chrLength = $2;
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
		s/\r|\n//g;
		my ( $chr, $feature, $start, $end, $attributes ) =
		  ( split( /\t/, $_ ) )[ ( 0, 2, 3, 4, 8 ) ];

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

sub pInfo {
	my $s       = shift;
	my $logFile = shift;
	print "[", scalar(localtime), "] $s\n";
	print $logFile "[", scalar(localtime), "] $s\n";
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
