package source::bamSummary;

use strict;
use warnings;

#use forks;
#use threads;
#use threads::shared;
use Exporter;
use File::Basename;

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

	#	my $maxThreads       = $config->{'maxThreads'};
	my $logFile        = $config->{'log'};
	my $caculateMethod = $config->{'caculateMethod'};  #1 as mean 2 as median
	
	die $!
	  if ( !defined($filelist)
		or !-e $filelist
		or !defined($resultDir)
		or ( defined($targetregionfile)  and !-e $targetregionfile )
		or ( defined($gtffile)           and !-e $gtffile )
		or ( !defined($targetregionfile) and !defined($gtffile) ) );
	if ( !( -e $resultDir . '/bamResult/' ) ) {
		mkdir $resultDir . '/bamResult/';
	}
	my $outputFile = $resultDir . '/bamResult/bamSummary.txt';
	my $methodText;
	if ($caculateMethod==1) {$methodText="Mean";} else {
		$methodText="Median";
	}
	
	#test R, comment below codes
	if ( $isdepth == 1 ) {
		pInfo(
"Will calculate the depth in on-/off-target regions. It will take a long time. You can turn it off without -d in your command line",
			$logFile
		);
	}
	else {
		$isdepth = 0;
	}

	# 1) Load targetregionfile
	my %targetregion;
	if ( defined($targetregionfile) ) {
		pInfo( "Load target region file '$targetregionfile'", $logFile );
		loadbed( $targetregionfile, \%targetregion );
	} else {
					pInfo(
"No targetregion file, will use exon regions in gtf file instead",
				$logFile
			);
	}

	# 2) Load gtf file
	my %gtf;
	if ( defined($gtffile) ) {
		pInfo( "Load gtf file '$gtffile'", $logFile );
		loadgtf( $gtffile, \%gtf );
	} else {
					pInfo(
"No gtf file",
				$logFile
			);
	}

	# 3) read bam file
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
			"Total Reads(Depth)",           "On-target(Depth)",
			"Off-target(Depth)",            "Off-target-intron(Depth)",
			"Off-target-intergenic(Depth)", "Off-target-mito(Depth)\n"
		  );
	}
	else {
		print OUT "\n";
	}
	while ( my $f = <IN> ) {
		$f =~ s/\r|\n//g;
		pInfo( "Processing $f ", $logFile );
		my @metric;
		if ( !defined($targetregionfile) ) {
			@metric = &getbammetric( $f, \%{ $gtf{'exon'} },
				\%gtf, $isdepth, $samtoolsBin, $caculateMethod );
		}
		else {
			@metric =
			  &getbammetric( $f, \%targetregion, \%gtf, $isdepth, $samtoolsBin,
				$caculateMethod );
		}
		print OUT join "\t", (@metric);
		print OUT "\n";

		#		if ( scalar( threads->list() ) < $maxThreads ) {
		#			pInfo("Processing $f ",$logFile);
		#			my ($t) =
		#			  threads->new( \&getbammetric, $f, \%targetregion, \%gtf, $isdepth,
		#				$samtoolsBin );
		#		}
		#		else {
		#			foreach my $thread ( threads->list() ) {
		#				my @metric = $thread->join;
		#				print OUT join "\t", (@metric);
		#				print OUT "\n";
		#				pInfo("Processing $f ",$logFile);
		#				my ($t) =
		#				  threads->new( \&getbammetric, $f, \%targetregion, \%gtf,
		#					$isdepth, $samtoolsBin );
		#				last;
		#			}
		#		}
	}

	#	#join all left threads
	#	foreach my $thread ( threads->list() ) {
	#		my @metric = $thread->join;
	#		print OUT join "\t", (@metric);
	#		print OUT "\n";
	#	}
	close OUT;

	#end comment by test R

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
	my ( $in, $targetref, $gtfref, $isdepth, $samtoolsBin, $caculateMethod ) =
	  @_;
	my (
		$total,           $ontarget,
		$offtarget,       $ummapped,
		$offtargetintron, $offtargetintergenic,
		$offtargetmito
	) = ( 0, 0, 0, 0, 0, 0, 0 );
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
	my $flag_read_unmapped = 0x0004;

	if ( !%{$targetref} ) {
		$targetref = $gtfref->{'exon'};
	}

	while (<BAM>) {
		my @line   = split "\t", $_;
		my $flag   = $line[1];
		my $MQ     = $line[4];
		my $MQflag = 0;
		if ( $MQ =~ /^\d+$/ ) {
			$MQflag = 1;
		}
		my $insert     = $line[8];
		my $insertflag = 0;
		if ( $insert =~ /^\d+$/ ) {
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
			&storeData( $MQ, $totalMQ, $caculateMethod ) if ($MQflag);
			&storeData( $insert, $totalinsert, $caculateMethod )
			  if ($insertflag);
			if ( exists( $targetref->{$chr}->{$pos} ) ) {
				$ontarget++;
				&storeData( $MQ, $ontargetMQ, $caculateMethod ) if ($MQflag);
				&storeData( $insert, $ontargetinsert, $caculateMethod )
				  if ($insertflag);
			}
			else {
				$offtarget++;
				&storeData( $MQ, $offtargetMQ, $caculateMethod ) if ($MQflag);
				&storeData( $insert, $offtargetinsert, $caculateMethod )
				  if ($insertflag);
				if ( exists( $gtfref->{'intron'}->{$chr}->{$pos} ) ) {
					$offtargetintron++;
					&storeData( $MQ, $offtargetintronMQ, $caculateMethod )
					  if ($MQflag);
					&storeData( $insert, $offtargetintroninsert,
						$caculateMethod )
					  if ($insertflag);
				}
				elsif ( !exists( $gtfref->{'exon'}->{$chr}->{$pos} ) ) {
					$offtargetintergenic++;
					&storeData( $MQ, $offtargetintergenicMQ, $caculateMethod )
					  if ($MQflag);
					&storeData( $insert, $offtargetintergenicinsert,
						$caculateMethod )
					  if ($insertflag);
				}
				if ( $chr =~ /M/i ) {
					$offtargetmito++;
					&storeData( $MQ, $offtargetmitoMQ, $caculateMethod )
					  if ($MQflag);
					&storeData( $insert, $offtargetmitoinsert, $caculateMethod )
					  if ($insertflag);
				}
			}
		}
	}
	close BAM;

	$totalMQ     = findMedianMean( $totalMQ,     $total, $caculateMethod );
	$totalinsert = findMedianMean( $totalinsert, $total, $caculateMethod );

	$ontargetMQ = findMedianMean( $ontargetMQ, $ontarget, $caculateMethod );
	$ontargetinsert =
	  findMedianMean( $ontargetinsert, $ontarget, $caculateMethod );

	$offtargetMQ = findMedianMean( $offtargetMQ, $offtarget, $caculateMethod );
	$offtargetinsert =
	  findMedianMean( $offtargetinsert, $offtarget, $caculateMethod );

	$offtargetintronMQ =
	  findMedianMean( $offtargetintronMQ, $offtargetintron, $caculateMethod );
	$offtargetintroninsert =
	  findMedianMean( $offtargetintroninsert, $offtargetintron,
		$caculateMethod );

	$offtargetintergenicMQ =
	  findMedianMean( $offtargetintergenicMQ, $offtargetintergenic,
		$caculateMethod );
	$offtargetintergenicinsert =
	  findMedianMean( $offtargetintergenicinsert, $offtargetintergenic,
		$caculateMethod );

	$offtargetmitoMQ =
	  findMedianMean( $offtargetmitoMQ, $offtargetmito, $caculateMethod );
	$offtargetmitoinsert =
	  findMedianMean( $offtargetmitoinsert, $offtargetmito, $caculateMethod );

	# if need to calculate the depth
	if ($isdepth) {

		# use samtools depth to get the depth file
		#		pInfo("Getting depth of '$in'",$logFile);
		open( DEPTH, "$samtoolsBin depth $in|" ) or die $!;
		my ( $totaldepth, $ontargetdepth, $offtargetdepth,
			$offtargetintrondepth, $offtargetintergenicdepth,
			$offtargetmitodepth );
		my ( $totalnum, $ontargetnum, $offtargetnum, $offtargetintronnum,
			$offtargetintergenicnum, $offtargetmitonum )
		  = ( 0, 0, 0, 0, 0, 0 );    # used as denominator
		while (<DEPTH>) {
			s/\r|\n//g;
			my ( $chr, $pos, $depth ) = split "\t";
			$totalnum++;
			&storeData( $depth, $totaldepth, $caculateMethod );
#			$totaldepth += $depth;
			if ( exists( $targetref->{$chr}->{$pos} ) ) {
				$ontargetnum++;
				&storeData( $depth, $ontargetdepth, $caculateMethod );
#				$ontargetdepth += $depth;
			}
			else {
				$offtargetnum++;
				&storeData( $depth, $offtargetdepth, $caculateMethod );
#				$offtargetdepth += $depth;
				if ( exists( $gtfref->{'intron'}->{$chr}->{$pos} ) ) {
					$offtargetintronnum++;
					&storeData( $depth, $offtargetintrondepth, $caculateMethod );
#					$offtargetintrondepth += $depth;
				}
				elsif ( !exists( $gtfref->{'exon'}->{$chr}->{$pos} ) ) {
					$offtargetintergenicnum++;
					&storeData( $depth, $offtargetintergenicdepth, $caculateMethod );
#					$offtargetintergenicdepth += $depth;
				}
				if ( $chr =~ /M/i ) {
					$offtargetmitonum++;
					&storeData( $depth, $offtargetmitodepth, $caculateMethod );
#					$offtargetmitodepth += $depth;
				}
			}
		}
		close DEPTH;

		$totaldepth     = findMedianMean( $totaldepth, $totalnum, $caculateMethod );
		$ontargetdepth  = findMedianMean( $ontargetdepth,  $ontargetnum, $caculateMethod );
		$offtargetdepth = findMedianMean( $offtargetdepth, $offtargetnum, $caculateMethod );
		$offtargetintrondepth =
		  findMedianMean( $offtargetintrondepth, $offtargetintronnum, $caculateMethod );
		$offtargetintergenicdepth =
		  findMedianMean( $offtargetintergenicdepth, $offtargetintergenicnum , $caculateMethod);
		$offtargetmitodepth =
		  findMedianMean( $offtargetmitodepth, $offtargetmitonum, $caculateMethod );

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
			$offtargetintergenicdepth,  $offtargetmitodepth
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

#sub mydevide {
#	my ( $a, $b ) = @_;
#	if ( $b == 0 ) {
#		return 'NA';
#	}
#	else {
#		return $a / $b;
#	}
#}

sub loadgtf {
	my ( $in, $ref ) = @_;
	open( IIN, $in ) or die $!;
	my %transcripts;
	while (<IIN>) {
		s/\r|\n//g;

		#		my (
		#			$chr,   $source, $feature, $start, $end,
		#			$score, $strand, $frame,   $attributes
		#		) = split "\t";
		my ( $chr, $feature, $start, $end, $attributes ) =
		  ( split( /\t/, $_ ) )[ ( 0, 2, 3, 4, 8 ) ];

#gene_id "ENSG00000237375"; transcript_id "ENST00000327822"; exon_number "1"; gene_name "BX072566.1"; transcript_name "BX072566.1
		if ( $feature eq "exon" ) {
			foreach my $temp ( $start .. $end ) {
				$ref->{'exon'}->{$chr}->{$temp} = undef;
			}

			#@{ $ref->{'exon'}->{$chr} }{ ( $start .. $end ) } =undef;

			#			my ( $gene_id, $transcript_id, $exon_number ) =
			#			  ( split( / \"|\"; /, $attributes ) )[ ( 1, 3, 5 ) ];
			my $transcript_id = ( split( / \"|\"; /, $attributes ) )[3];
			$transcripts{$transcript_id}->{'chr'} = $chr;
			if ( !exists( $transcripts{$transcript_id}->{'start'} )
				or $start < $transcripts{$transcript_id}->{'start'} )
			{
				$transcripts{$transcript_id}->{'start'} = $start;
			}
			if ( !exists( $transcripts{$transcript_id}->{'end'} )
				or $end > $transcripts{$transcript_id}->{'end'} )
			{
				$transcripts{$transcript_id}->{'end'} = $end;
			}

			#			$transcripts{$transcript_id}->{'chr'}    = $chr;
			#			$transcripts{$transcript_id}->{'strand'} = $strand;
			#			$transcripts{$transcript_id}->{'exons'}->{$exon_number} =
			#			  [ $start, $end ];
		}
	}
	close IIN;

	#	#loop transcripts to get introns
	#	foreach my $t ( keys %transcripts ) {
	#		my @exons = sort { $a <=> $b } keys %{ $transcripts{$t}->{'exons'} };
	#		next if ( scalar(@exons) == 1 );    #only one exons
	#		my $chr = $transcripts{$t}->{'chr'};
	#		if ( $transcripts{$t}->{'strand'} eq "+" ) {
	#			for ( my $i = 1 ; $i < @exons ; $i++ ) {
	#				@{ $ref->{'intron'}->{$chr} }{
	#					(
	#						(
	#							$transcripts{$t}->{'exons'}->{ $exons[ $i - 1 ] }
	#							  ->[1] + 1
	#						) .. (
	#							$transcripts{$t}->{'exons'}->{ $exons[$i] }->[0] - 1
	#						)
	#					)
	#				  }
	#				  = undef;
	#			}
	#		}
	#		else {
	#			for ( my $i = 1 ; $i < @exons ; $i++ ) {
	#				@{ $ref->{'intron'}->{$chr} }{
	#					(
	#						(
	#							$transcripts{$t}->{'exons'}->{ $exons[$i] }->[1] + 1
	#						) .. (
	#							$transcripts{$t}->{'exons'}->{ $exons[ $i - 1 ] }
	#							  ->[0] - 1
	#						)
	#					)
	#				  }
	#				  = undef;
	#			}
	#		}
	#	}
	foreach my $transcript_id ( keys %transcripts ) {
		my $start = $transcripts{$transcript_id}->{'start'};
		my $end   = $transcripts{$transcript_id}->{'end'};
		my $chr   = $transcripts{$transcript_id}->{'chr'};
		foreach my $temp ( $start .. $end ) {
			$ref->{'intron'}->{$chr}->{$temp} = undef;
		}
	}
	foreach my $chr ( sort keys %{ $ref->{'intron'} } )
	{    #delete all exon regions in intron
		delete @{ $ref->{'intron'}->{$chr} }
		  { ( keys %{ $ref->{'exon'}->{$chr} } ) };
	}
	undef(%transcripts);
}

# Memory intensive
sub loadbed {
	my ( $in, $ref ) = @_;
	open( IIN, $in ) or die $!;
	while (<IIN>) {
		s/\r|\n//g;
		my ( $chr, $start, $end ) = split "\t";
		foreach ( my $i = $start + 1 ; $i <= $end ; $i++ ) {
			$ref->{$chr}->{$i} = undef;
		}
	}
	close IIN;
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
		${ $_[1] }[$input]++;
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
		if ( $arrayLength % 2 == 0 ) {
			my $medianPre;
			for ( my $count = 0 ; $count < scalar( @{ $_[0] } ) ; $count++ ) {
				if ( !defined( ${ $_[0] }[$count] ) ) { next; }
				$totalCount += ${ $_[0] }[$count];
				if ( $totalCount == $arrayLength / 2 ) {
					$medianPre = $count;
				}
				elsif ( $totalCount > $arrayLength / 2 ) {
					if ( defined $medianPre ) {
						return ( ( $medianPre + $count ) / 2 );
					}
					else {
						return ($count);
					}
				}
			}
		}
		else {
			for ( my $count = 0 ; $count < scalar( @{ $_[0] } ) ; $count++ ) {
				if ( !defined( ${ $_[0] }[$count] ) ) { next; }
				$totalCount += ${ $_[0] }[$count];
				if ( $totalCount > $arrayLength / 2 ) {
					return ($count);
				}
			}
		}
		return('NA1');
	}
}
