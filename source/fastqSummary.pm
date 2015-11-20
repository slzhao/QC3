package source::fastqSummary;

use strict;
use warnings;
use Exporter;
use File::Basename;
use threads;
use threads::shared;

our @ISA    = qw(Exporter);
our @EXPORT = qw(fastqSummary);

my $rSourceLocation = '/source/fastq_plot_byPerl.R';
my $current : shared;
my @resultOut : shared;
my @result1 : shared;
my @result2 : shared;
my @result3 : shared;
my @result4 : shared;

my $logRef = \@main::log;

1;

sub fastqSummary {

	# each line is a file name in $filelist
	my ( $filelist, $config ) = @_;
	my $resultDir  = $config->{'resultDir'};
	my $singleEnd  = $config->{'singleEnd'};
	my $RBin       = $config->{'RBin'};
	my $maxThreads = $config->{'maxThreads'};

	#	$logFile    = $config->{'log'};
	my $rePlot   = $config->{'rePlot'};
	my $showHelp = $config->{'showHelp'};

	my $usage =
"Program: qc3.pl (a quality control tool for DNA sequencing data in raw data, alignment, and variant calling stages)
Version: $main::version

Module: fastq QC
Usage:   perl qc3.pl -m f -i fastqFileList -o outputDirectory [-t threads] [other options]

Options:

	-i	input filelist      Required. A list file for fastq files to be analyzed. QC3 also support .gz files.
	-o	output directory    Required. Output directory for QC result. If the directory doesn't exist, it would be created.
	-t	threads             Optional. Threads used in analysis. The default value is 4.
	-se	                    Optional. Indicating the fastq files were single-end data. QC3 takes fastq files as pair-end data by default.
	
	-h	                    Optional. Show this information.
	
For more information, please refer to the readme file in QC3 directory. Or visit the QC3 website at https://github.com/slzhao/QC3

";
	die "\n$usage"
	  if ($showHelp);

	if ( !( -e $resultDir . '/fastqResult/' ) ) {
		mkdir $resultDir . '/fastqResult/';
	}
	my $outputFile = $resultDir . '/fastqResult/fastqSummary.txt';

	$| = 1;
	if ($rePlot) {    #re plot by R only, comment below codes
		pInfo(
"Parameter -rp was used. The program will use the result files in output directory to re-generate the report",
			$logRef
		);
	}
	else {

		#test R, comment below code

		open( IN, $filelist ) or die $!;
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
		close IN;

		my @threads;     #open threads and get results
		$current = 0;
		foreach my $x ( 1 .. $maxThreads ) {
			push @threads,
			  threads->new( \&fastqProcess, $x, \@fileList, $logRef );
		}

		foreach my $thread (@threads) {
			my $threadNum = $thread->join;
			pInfo( "Thread $threadNum finished", $logRef );
		}

		#		while (my $thread=shift @threads) {
		#			if ($thread->is_joinable()) {
		#				my $threadNum = $thread->join;
		#				pInfo( "Thread $threadNum finished", $logFile );
		#			} else {
		#				push @threads,$thread;
		#			}
		#		}

		#write results to file
		open OUT, ">$outputFile" or die $!;
		print OUT join "\t",
		  (
			"#Sample", "Instrument", "RunNumber", "Flowcell",
			"Lane",    "TotalReads", "Reads(Y)",  "Reads(N)",
			"BQ",      "BQ(Y)",      "BQ(N)",     "GC",
			"GC(Y)",   "GC(N)\n"
		  );
		open( INFORMATION1, ">$outputFile.score.txt" )  or die $!;
		open( INFORMATION2, ">$outputFile.nuc.txt" )    or die $!;
		open( INFORMATION3, ">$outputFile.scoreN.txt" ) or die $!;
		open( INFORMATION4, ">$outputFile.nucN.txt" )   or die $!;
		print INFORMATION1 "File\n";
		print INFORMATION2 "File\tA1\tT1\tC1\tG1\n";
		print INFORMATION3 "File\n";
		print INFORMATION4 "File\tA1\tT1\tC1\tG1\n";

		foreach my $result (@resultOut) {
			print OUT "$result\n";
		}
		foreach my $result (@result1) {
			print INFORMATION1 "$result\n";
		}
		foreach my $result (@result2) {
			print INFORMATION2 "$result\n";
		}
		foreach my $result (@result3) {
			print INFORMATION3 "$result\n";
		}
		foreach my $result (@result4) {
			print INFORMATION4 "$result\n";
		}
		close OUT;
		close INFORMATION1;
		close INFORMATION2;
		close INFORMATION3;
		close INFORMATION4;

		#test R, end comment
	}

	#plot by R
	my $Rsource =
	  dirname($0) . "/source/rFunctions.R " . dirname($0) . $rSourceLocation;
	my $rResult = system(
"cat $Rsource | $RBin --vanilla --slave --args $resultDir $singleEnd > $resultDir/fastqResult/fastqSummary.rLog"
	);
	pInfo( "Finish fastq summary!", $logRef );
	return ($rResult);
}

sub fastqProcess {
	my ( $threadNum, $fileListRef, $logRef ) = @_;
	pInfo( "Thread $threadNum stared", $logRef );
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
		my ( $metric, $info1, $info2, $info3, $info4 ) = &getmetric($file);
		{
			lock @resultOut;
			$resultOut[$current_temp] = join "\t", ( @{$metric} );
		}
		{
			lock @result1;
			$result1[$current_temp] = join "\t", ( @{$info1} );
		}
		{
			lock @result2;
			$result2[$current_temp] = join "\t", ( @{$info2} );
		}
		{
			lock @result3;
			$result3[$current_temp] = join "\t", ( @{$info3} );
		}
		{
			lock @result4;
			$result4[$current_temp] = join "\t", ( @{$info4} );
		}
	}
}

sub getmetric {
	my ($in) = @_;
	my @fileLable = ( split /\t/, $in );
	$in = $fileLable[0];
	my $label = $in;
	if ( defined $fileLable[1] ) {
		$label = $fileLable[1];
	}
	my @r1         = ($label);    #return values
	my @r2         = ($label);    #return values
	my @r3         = ($label);    #return values
	my @r4         = ($label);    #return values
	my @r5         = ($label);    #return values
	my $instrument = "";
	my $run        = "";
	my $flowcell   = "";
	my $lane       = "";
	my $failed     = "";

	my $first = 1;

	my $totalreads  = 0;
	my $totalreadsy = 0;
	my $totalreadsn = 0;

	my $bq  = 0;
	my $bqy = 0;
	my $bqn = 0;

	my $gc  = 0;
	my $gcy = 0;
	my $gcn = 0;

	my $totalnuclear  = 0;
	my $totalnucleary = 0;
	my $totalnuclearn = 0;

	my $offset = 33;    #default;

	#will use the top 40 reads to guess the offset
	#	my $guessoffset = 0;
	my @scores;

	#store some information needed
	my $readlen;
	my @scoreSumBase  = ();
	my @nucSumBase    = ();
	my @scoreSumBaseN = ();
	my @nucSumBaseN   = ();

#find information in title, such as machine, lane. Guess offset with top 40 sequences
	my $countTemp = 0;
	if ( $in =~ /\.gz$/ ) {
		open( TITLE, "zcat $in|" ) or die $!;
	}
	else {
		open( TITLE, $in ) or die $!;
	}
	while ( my $line1 = <TITLE> ) {
		<TITLE>;
		<TITLE>;
		my $line4 = <TITLE>;
		if ($first) {
			if ( $line1 =~
/@(.*?):(.*?):(.*?):(.*?):(.*?):(.*?):(.*?)\s+(.*?):(.*?):(.*?):(.*?)/
			  )
			{    #casava 1.8 format
				( $instrument, $run, $flowcell, $lane ) = ( $1, $2, $3, $4 );
				$first = 0;
			}    #not casava 1.8 format, but Illumina format
			elsif ( $line1 =~ /@(.*?):(.*?):(.*?):(.*?):(.*?)/ ) {
				pInfo(
"$in was Not casava 1.8 format, Can't find information for run ID and flowcell ID, use None instead",
					$logRef
				);
				( $instrument, $lane, $run, $flowcell ) = ( $1, $2, $3, $4 );
				$run      = 'None';
				$flowcell = 'None';
				$first    = 0;
			}
			else {
				die
"Sequence identifiers were Not Illumina format, Can't find information for instrument, run, flowcell and lane";
			}
		}
		my @tmpscores = map( ord, ( split //, $line4 ) );
		push @scores, @tmpscores;

		$countTemp++;
		if ( $countTemp >= 40 ) { last; }
	}
	close(TITLE);
	if ( $countTemp < 40 ) {
		print STDERR "Please note $in has less than 40 reads";
	}
	$offset = guessoffset(@scores);

	if ( $in =~ /\.gz$/ ) {
		open( IIN, "zcat $in|" ) or die $!;
	}
	else {
		open( IIN, $in ) or die $!;
	}
	while ( my $line1 = <IIN> ) {
		my $line2 = <IIN>;
		chomp $line2;
		<IIN>;
		my $line4 = <IIN>;
		$failed = ( split /:/, $line1 )[7];
		if ( !defined $failed ) {
			$failed = "N";    #in case Not casava 1.8 format
		}

		$totalreads++;

		$readlen = length($line2);

		$totalnuclear += $readlen;

		my $gcnum = $line2 =~ tr/[GCgc]//;
		$gc += $gcnum;

		my @tmpscores;
		foreach my $key ( 0 .. ( $readlen - 1 ) ) {

			#score for each base
			push @tmpscores, ord( substr $line4, $key, 1 );
			$scoreSumBase[$key] += $tmpscores[$key];

			#nuc for each base
			$nucSumBase[$key]{ substr $line2, $key, 1 }++;
		}

		my $tmpbq = mymean(@tmpscores);
		$bq += $tmpbq;

		if ( $failed eq "Y" ) {
			$bqy += $tmpbq;
			$gcy += $gcnum;
			$totalreadsy++;
			$totalnucleary += $readlen;
		}
		else {
			$bqn += $tmpbq;
			$gcn += $gcnum;
			$totalreadsn++;
			$totalnuclearn += $readlen;

			#no Y scores and nuc
			foreach my $key ( 0 .. ( $readlen - 1 ) ) {

				#score for each base
				$scoreSumBaseN[$key] += $tmpscores[$key];

				#nuc for each base
				$nucSumBaseN[$key]{ substr $line2, $key, 1 }++;
			}
		}
	}
	close IIN;

	#print to information file
	foreach my $key ( 0 .. ( $readlen - 1 ) ) {
		push @r2, $scoreSumBase[$key] / $totalreads - $offset;
		push @r4, $scoreSumBaseN[$key] / $totalreadsn - $offset;
		my $totalNuc =
		  $nucSumBase[$key]{"A"} +
		  $nucSumBase[$key]{"T"} +
		  $nucSumBase[$key]{"C"} +
		  $nucSumBase[$key]{"G"};    #percenty for each nuc, what about U?
		foreach my $nuc ( "A", "T", "C", "G" ) {
			push @r3, &myDivide( $nucSumBase[$key]{$nuc}, $totalNuc );
		}
		$totalNuc =
		  $nucSumBaseN[$key]{"A"} +
		  $nucSumBaseN[$key]{"T"} +
		  $nucSumBaseN[$key]{"C"} +
		  $nucSumBaseN[$key]{"G"};    #percenty for each nuc, what about U?
		foreach my $nuc ( "A", "T", "C", "G" ) {
			push @r5, &myDivide( $nucSumBaseN[$key]{$nuc}, $totalNuc );
		}
	}

	push @r1,
	  (
		$instrument, $run, $flowcell, $lane, $totalreads, $totalreadsy,
		$totalreadsn
	  );

	#add bq
	my $tempY = 0;
	my $tempN = 0;
	if ( $totalreadsy == 0 ) {
		$tempY = 'NA';
	}
	else {
		$tempY = &myDivide( $bqy, $totalreadsy ) - $offset;
	}
	if ( $totalreadsn == 0 ) {
		$tempN = 'NA';
	}
	else {
		$tempN = &myDivide( $bqn, $totalreadsn ) - $offset;
	}
	push @r1, ( &myDivide( $bq, $totalreads ) - $offset, $tempY, $tempN, );

	#add gc
	push @r1,
	  (
		&myDivide( $gc,  $totalnuclear ),
		&myDivide( $gcy, $totalnucleary ),
		&myDivide( $gcn, $totalnuclearn )
	  );

	return ( \@r1, \@r2, \@r3, \@r4, \@r5 );
}

sub mymean {
	my @a = @_;
	unless ( scalar(@a) ) {
		return 0;
	}
	my $t = 0;
	foreach my $tt (@a) {
		$t += $tt;
	}
	return $t / scalar(@a);
}

sub guessoffset {
	my (@a) = @_;
	@a = sort { $a <=> $b } @a;
	my $offset = 33;
	my $lowest = $a[0];
	if ( $lowest < 59 ) {
		$offset = 33;
	}
	elsif ( $lowest < 64 ) {
		$offset = 59;
	}
	else {
		$offset = 64;
	}
	return $offset;
}

sub pInfo {
	my $s = $_[0];
	print "[", scalar(localtime), "] $s\n";

	lock $_[1];
	{
		lock @{ $_[1] };
		push @{ $_[1] }, "[" . scalar(localtime) . "] $s\n";
	}
}

sub myDivide {
	my ( $a, $b ) = @_;
	if ( !( defined $a ) or !( defined $b ) or $b == 0 ) {
		return 'NA';
	}
	else {
		return ( $a / $b );
	}
}
