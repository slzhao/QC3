package source::fastqSummary;

use strict;
use warnings;
use Exporter;
use File::Basename;
use threads;

our @ISA    = qw(Exporter);
our @EXPORT = qw(fasqSummary);

1;

my $rSourceLocation = '/source/fastq_plot_byPerl.R';

sub fasqSummary {

	# each line is a file name in $filelist
	my ( $filelist, $config ) = @_;
	my $resultDir  = $config->{'resultDir'};
	my $pairEnd    = $config->{'pairEnd'};
	my $RBin       = $config->{'RBin'};
	my $maxThreads = $config->{'maxThreads'};
	my $logFile = $config->{'log'};
	my $usage =
"Please use the follow command:\n main.pl -m f -i inputFastqList -o outputDir [-p]\nFor more information, plase read the readme file\n";
	die "$!\n$usage"
	  if ( !defined($filelist) or !-e $filelist or !defined($resultDir) );

	if ( !( -e $resultDir . '/fastqResult/' ) ) {
		mkdir $resultDir . '/fastqResult/';
	}
	my $outputFile = $resultDir . '/fastqResult/fastqSummary.txt';

	#test R, comment below code
	$| = 1;
	open OUT, ">$outputFile" or die $!;
	print OUT join "\t",
	  (
		"#Sample", "Instrument", "RunNumber", "Flowcell",
		"Lane",    "TotalReads", "Reads(Y)",  "Reads(N)",
		"BQ",      "BQ(Y)",      "BQ(N)",     "GC",
		"GC(Y)",   "GC(N)\n"
	  );
	open( IN,           $filelist )                or die $!;
	open( INFORMATION1, ">$outputFile.score.txt" ) or die $!;
	open( INFORMATION2, ">$outputFile.nuc.txt" )   or die $!;
	open( INFORMATION3, ">$outputFile.scoreN.txt" ) or die $!;
	open( INFORMATION4, ">$outputFile.nucN.txt" )   or die $!;
	print INFORMATION1 "File\n";
	print INFORMATION2 "File\tA1\tT1\tC1\tG1\n";
	print INFORMATION3 "File\n";
	print INFORMATION4 "File\tA1\tT1\tC1\tG1\n";

	while ( my $f = <IN> ) {
		$f =~ s/\r|\n//g;
		if ( scalar( threads->list() ) < $maxThreads ) {
			pInfo("Processing $f ",$logFile);
			my ($t) = threads->new( \&getmetric, $f );
		}
		else {
			foreach my $thread ( threads->list() ) {
				my ( $metric, $info1, $info2, $info3, $info4 ) = $thread->join;
				print OUT join "\t", ( @{$metric} );
				print OUT "\n";
				print INFORMATION1 join "\t", ( @{$info1} );
				print INFORMATION1 "\n";
				print INFORMATION2 join "\t", ( @{$info2} );
				print INFORMATION2 "\n";
				print INFORMATION3 join "\t", ( @{$info3} );
				print INFORMATION3 "\n";
				print INFORMATION4 join "\t", ( @{$info4} );
				print INFORMATION4 "\n";
				pInfo("Processing $f ",$logFile);
				my ($t) = threads->new( \&getmetric, $f );
				last;
			}
		}
	}
	close IN;

	#join all left threads
	foreach my $thread ( threads->list() ) {
		my ( $metric, $info1, $info2, $info3, $info4 ) = $thread->join;
		print OUT join "\t", ( @{$metric} );
		print OUT "\n";
		print INFORMATION1 join "\t", ( @{$info1} );
		print INFORMATION1 "\n";
		print INFORMATION2 join "\t", ( @{$info2} );
		print INFORMATION2 "\n";
		print INFORMATION3 join "\t", ( @{$info3} );
		print INFORMATION3 "\n";
		print INFORMATION4 join "\t", ( @{$info4} );
		print INFORMATION4 "\n";
	}
	close OUT;
	close INFORMATION1;
	close INFORMATION2;
	close INFORMATION3;
	close INFORMATION4;

	#plot by R
	my $Rsource =
	  dirname($0) . "/source/rFunctions.R " . dirname($0) . $rSourceLocation;
	my $rResult = system(
"cat $Rsource | $RBin --vanilla --slave --args $resultDir $pairEnd > $resultDir/fastqResult/fastqSummary.rLog"
	);
	pInfo("Finish fastq summary!",$logFile);
	return ($rResult);
}

sub getmetric {
	my ($in) = @_;
	if ( $in =~ /\.gz$/ ) {
		open( IIN, "zcat $in|" ) or die $!;
	}
	else {
		open( IIN, $in ) or die $!;
	}
	my @r1         = ($in);    #return values
	my @r2         = ($in);    #return values
	my @r3         = ($in);    #return values
	my @r4         = ($in);    #return values
	my @r5         = ($in);    #return values
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
	my $guessoffset = 0;
	my @scores;

	#store some information needed
	my $readlen;
	my %scoreSumBase = ();
	my %nucSumBase   = ();
	my %scoreSumBaseN = ();
	my %nucSumBaseN   = ();

	while ( my $line1 = <IIN> ) {
		my $line2 = <IIN>;
		my $line3 = <IIN>;
		my $line4 = <IIN>;
		if ( $first == 1 ) {
			if ( $line1 =~
/@(.*?):(.*?):(.*?):(.*?):(.*?):(.*?):(.*?)\s+(.*?):(.*?):(.*?):(.*?)/
			  )
			{
				( $instrument, $run, $flowcell, $lane ) = ( $1, $2, $3, $4 );
				$first = 0;
			}
			else {
				print STDERR "Not in casava 1.8 format\n";
			}
		}
		$failed = ( split /:/, $line1 )[7];

#        if($line1=~/@(.*?):(.*?):(.*?):(.*?):(.*?):(.*?):(.*?)\s+(.*?):(.*?):(.*?):(.*?)/){
#            my($instrument1,$runid,$flowcellid,$lane1,$tilenumber,$x,$y,$pair,$failed,$n,$index)= ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11);
#            if($first){
#                $instrument=$instrument1;
#                $run=$runid;
#                $flowcell=$flowcellid;
#                $lane=$lane1;
#                $first=0;
#            }
		$totalreads++;

		$readlen = length($line2) - 1;    #note the "/n" at the end

		$totalnuclear += $readlen;

		#            my $gcnum=0;
		#            foreach my $s(split //,$line2){
		#                $gcnum++ if($s=~/c|g/i);
		#            }
		my $gcnum = () = $line2 =~ /[GC]/gi;
		$gc += $gcnum;

		#		my @tmpscores = ();
		#		foreach my $s ( split //, $line4 ) {
		#			push @tmpscores, ord($s);
		#		}
		my @tmpscores = map( ord, ( split //, $line4 ) );

		my @temp = ( split //, $line2 );
		foreach my $key ( 0 .. ( $readlen - 1 ) ) {

			#score for each base
			$scoreSumBase{$key} += $tmpscores[$key];

			#nuc for each base
			$nucSumBase{$key}{ $temp[$key] }++;
		}

		my $tmpbq = mymean(@tmpscores);
		$bq += $tmpbq;

		if ( $totalreads < 10 ) {
			push @scores, @tmpscores;
		}
		elsif ( !$guessoffset ) {

			# guess the offset
			$offset      = guessoffset(@scores);
			$guessoffset = 1;
		}

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
				$scoreSumBaseN{$key} += $tmpscores[$key];
				#nuc for each base
				$nucSumBaseN{$key}{ $temp[$key] }++;
			}
		}

		#        }else{
		#            print STDERR "Not in casava 1.8 format\n";
		#        }
	}
	close IIN;

	if ( !$guessoffset ) {
		$offset = guessoffset(@scores);
	}

	#print to information file
	foreach my $key ( 0 .. ( $readlen - 1 ) ) {
		push @r2, $scoreSumBase{$key} / $totalreads - $offset;
		push @r4, $scoreSumBaseN{$key} / $totalreadsn - $offset;
		my $totalNuc =
		  $nucSumBase{$key}{"A"} +
		  $nucSumBase{$key}{"T"} +
		  $nucSumBase{$key}{"C"} +
		  $nucSumBase{$key}{"G"};    #percenty for each nuc, what about U?
		foreach my $nuc ( "A", "T", "C", "G" ) {
			push @r3, $nucSumBase{$key}{$nuc} / $totalNuc;
		}
		$totalNuc =
		  $nucSumBaseN{$key}{"A"} +
		  $nucSumBaseN{$key}{"T"} +
		  $nucSumBaseN{$key}{"C"} +
		  $nucSumBaseN{$key}{"G"};    #percenty for each nuc, what about U?
		foreach my $nuc ( "A", "T", "C", "G" ) {
			push @r5, $nucSumBaseN{$key}{$nuc} / $totalNuc;
		}
	}

	push @r1,
	  (
		$instrument, $run, $flowcell, $lane, $totalreads, $totalreadsy,
		$totalreadsn
	  );

	#add bq
	push @r1,
	  (
		$bq / $totalreads - $offset,
		$bqy / $totalreadsy - $offset,
		$bqn / $totalreadsn - $offset
	  );

	#add gc
	push @r1,
	  ( $gc / $totalnuclear, $gcy / $totalnucleary, $gcn / $totalnuclearn );

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
	my $s = shift;
	my $logFile = shift;
	print "[", scalar(localtime), "] $s\n";
	print $logFile "[", scalar(localtime), "] $s\n";
}
