package source::vcfSummary;

use strict;
use warnings;
use Exporter;
use File::Basename;

#use List::Compare;

our @ISA    = qw(Exporter);
our @EXPORT = qw(vcfSummary);
use Env '@PATH';

1;

my $rSourceLocation = '/source/vcf_plot_byPerl.R';
my $logRef          = \@main::log;

sub vcfSummary {

	my ( $VCF, $config ) = @_;
	my $vcfCfgFile = $config->{'vcfCfgFile'};
	my $method     = $config->{'method'};
	my $resultDir  = $config->{'resultDir'};
	my $usePASS    = $config->{'usePASS'};

	#	my $logFile   = $config->{'log'};
	my $rePlot = $config->{'rePlot'};

	my $RBin           = $config->{'RBin'};
	my $annovarBin     = $config->{'annovarBin'};
	my $annovarConvert = $config->{'annovarConvert'};
	my $annovarOption  = $config->{'annovarOption'};
	my $annovarDb      = $config->{'annovarDb'};

	my $showHelp = $config->{'showHelp'};

	my $usage =
"Program: qc3.pl (a quality control tool for DNA sequencing data in raw data, alignment, and variant calling stages)
Version: $main::version

Module: vcf QC
Usage:   perl qc3.pl -m v -i vcfFile -o outputDirectory [other options]

Options:

	-i	vcfFile                 Required. A vcf file to be analyzed.
	-o	output directory        Required. Output directory for QC result. If the directory doesn't exist, it would be created.
	-s	method                  Optional. Method used in consistence calculation, should be 1 or 2.
	-c	config file             Optional. A file indicating the filter arguments for vcf files. If not specified, the default file 'GATK.cfg' in QC3 directory with GATK best practices recommended arguments will be used.
	-up	use PASS only           Optional. Only use PASS variants in vcf. QC3 will use all variants in vcf file except LowQual by default.
	-a	annotation directory	Optional. Directory of annovar database.

	-h	                        Optional. Show this information.
	
For more information, please refer to the readme file in QC3 directory. Or visit the QC3 website at https://github.com/slzhao/QC3

";

	die "\n$usage" if ($showHelp);
	die "Can't find the vcf QC config file\n$usage" if ( !( -e $vcfCfgFile ) );
	die "Method(-s) should be 1 or 2\n$usage"
	  if ( $method != 1 and $method != 2 );

	$| = 1;
	my $doAnnovar = 1;
	if ( ( grep -x "$_/$annovarBin", @PATH ) ) {
	}
	elsif ( -e $annovarBin ) {
		$annovarBin = "perl $annovarBin";
	}
	else {
		$doAnnovar = 0;
		pInfo(
"Can't find ANNOVAR bin. Will not perform ANNOVAR annotation but other functions will work well. You can read readme file to find how to download ANNOVAR",
			$logRef
		);
	}
	if ( !( defined $annovarDb ) or !( -e $annovarDb ) ) {
		$doAnnovar = 0;
		pInfo(
"Can't find ANNOVAR database. Will not perform ANNOVAR annotation but other functions will work well. You can read readme file to find how to download ANNOVAR",
			$logRef
		);
	}

	my %formatPos;
	my %IDList;
	my $done_filter;
	my @titles;
	my $sampleSize;
	my %Sample2NumberAll;
	my %Sample2NumberFilter;

	my %total;
	my %selected;
	my %total2;
	my %selected2;
	my %changNumber;
	my %changes;
	my %cfgFilter;
	my %snpCount;
	my %genderCount;
	my %snpNucCount;

	open CFG, "<$vcfCfgFile" or die "Can't read $vcfCfgFile\n$!";
	pInfo( "VariantFiltration based on $vcfCfgFile", $logRef );
	while (<CFG>) {
		chomp;
		if (/^#/) {
			next;
		}
		if (/[><=]+/) {
			$cfgFilter{$`} = $& . $';

			#			print "$` $& $'\n";
		}
	}

	my $filename         = basename($VCF);
	my $annovarResultDir = "$resultDir/vcfAnnovarResult/";

	if ($rePlot) {    #re plot by R only, comment below codes
		pInfo(
"Parameter -rp was used. The program will use the result files in output directory to re-generate the report",
			$logRef
		);
	}
	else {

		#test R, comment below code
		open READ, "<$VCF" or die "can't read $VCF\n";
		if ( !( -e ("$resultDir/vcfResult/") ) ) {
			mkdir "$resultDir/vcfResult/";
		}
		if ( !( -e $annovarResultDir ) ) { mkdir $annovarResultDir; }
		open RESULT1, ">$resultDir/vcfResult/$filename.IDlistAll.txt" or die $!;
		open RESULT2, ">$resultDir/vcfResult/$filename.IDlistFilter.txt"
		  or die $!;
		open RESULT3, ">$resultDir/vcfResult/$filename.SampleNumber.txt"
		  or die $!;

		my $resultFile = "$resultDir/vcfResult/$filename.Method$method.txt";
		open RESULT4, ">$resultFile" or die "can't write $resultFile\n";
		open RESULT5, ">$resultDir/vcfResult/$filename.snpCount.txt" or die $!;
		open RESULT6, ">$annovarResultDir$filename.pass"
		  or die $!;    #new vcf for annovar
		open RESULT7, ">$resultDir/vcfResult/$filename.sexCheck.txt" or die $!;
		open RESULT8, ">$resultDir/vcfResult/$filename.snpNucCount.txt"
		  or die $!;

		while (<READ>) {    #read title and ID list
			chomp;
			s/\r//g;
			if (/^##INFO=<ID=(\w+),Number=1/) {    #all IDs
				$IDList{$1} = "";
			}

  #			elsif (/^##FORMAT=<ID=(\w+),Number=[1G\.]/) { #all formats, find GT and AD
  #				print "$1\n";
  #				$formatPos{$1} = (scalar keys %formatPos);
  #			}
			elsif (/^##/) {
				next;
			}
			elsif (/^#/) {                         #title
				@titles = ( split /\t/, $_ );
				$sampleSize = scalar(@titles) - 9;
				pInfo( "The VCF file has $sampleSize samples", $logRef );
			}
			else {
				my @formats = ( split /:/, ( split /\t/, $_ )[8] );
				foreach my $format (@formats) {
					$formatPos{$format} = ( scalar keys %formatPos );
				}
				last;
			}
		}
		close(READ);
		if ( !exists( $formatPos{"GT"} ) or !exists( $formatPos{"AD"} ) ) {
			die "Can't find GT or AD in vcf file!\n";
		}
		if ( exists( $formatPos{"RD"} ) ) {
			pInfo(
"Find both RD and AD in VCF file, will use RD as ref and AD as alt",
				$logRef
			);
			if ( $method == 2 ) {
				$method = 1;
				pInfo(
"Method will be changed to 1 as this vcf file did NOT contain reads for all genetypes",
					$logRef
				);
			}
		}

		print RESULT1 "CHROM\tPOS\tREF\tALT";
		print RESULT2 "CHROM\tPOS\tREF\tALT";
		foreach my $ID ( sort keys %IDList ) {
			print RESULT1 "\t$ID";
			print RESULT2 "\t$ID";
		}
		print RESULT1 "\n";
		print RESULT2 "\n";

		open READ, "<$VCF" or die "can't read $VCF\n";
		while (<READ>) {    #read content
			chomp;
			s/\r//g;
			if ( /^#/ or $_ eq "" ) {
				print RESULT6 "$_\n";
				next;
			}
			my @lines = ( split /\t/, $_ );
			if (   length( $lines[3] ) != length( $lines[4] )
				or $lines[3] eq '.'
				or $lines[4] eq '.' )
			{               #insertion or deletion
				next;
			}

			my $REF  = $lines[3];
			my @ALTs = split( /,|\//, $lines[4] );
			my @SNPs = ( $REF, @ALTs );
			if (   ( $usePASS == 1 and $lines[6] ne "PASS" )
				or ( $usePASS != 1 and $lines[6] eq "LowQual" )
				or &filter( $lines[7], \%cfgFilter ) == 0 )
			{
				$done_filter = 1;
			}
			else {
				$done_filter = 0;
				print RESULT2 "$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]";
				print RESULT6 join( "\t", @lines ) . "\n";
			}
			print RESULT1 "$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]";

			my @ID1 = ( split /;/, $lines[7] );
			my %IDNumber;
			foreach my $ID (@ID1) {
				my @ID2 = ( split /=|,/, $ID );
				$IDNumber{ $ID2[0] } = $ID2[1];
			}
			if ( $done_filter == 0 ) {
				foreach my $ID ( sort keys %IDList ) {
					if ( defined $IDNumber{$ID} ) {
						print RESULT1 "\t$IDNumber{$ID}";
						print RESULT2 "\t$IDNumber{$ID}";
					}
					else {
						print RESULT1 "\t";
						print RESULT2 "\t";
					}
				}
				print RESULT1 "\n";
				print RESULT2 "\n";
			}
			else {
				foreach my $ID ( sort keys %IDList ) {
					if ( defined $IDNumber{$ID} ) {
						print RESULT1 "\t$IDNumber{$ID}";
					}
					else {
						print RESULT1 "\t";
					}
				}
				print RESULT1 "\n";
			}

			for ( my $x = 9 ; $x < ( 9 + $sampleSize ) ; $x++ ) {
				if ( !(defined $lines[$x]) or $lines[$x] =~ /^\./ ) { next; }
				$lines[$x] =~ /^(\d)\/(\d)/;
				my $allele1 = $1;
				my $allele2 = $2;
				if ( $allele1 ne $allele2 ) {
					if (   ( $SNPs[$allele1] . $SNPs[$allele2] eq "AG" )
						or ( $SNPs[$allele1] . $SNPs[$allele2] eq "GA" )
						or ( $SNPs[$allele1] . $SNPs[$allele2] eq "CT" )
						or ( $SNPs[$allele1] . $SNPs[$allele2] eq "TC" ) )
					{
						$Sample2NumberAll{ $titles[$x] }{"Transitions"}++;
						if ( $done_filter == 0 ) {
							$Sample2NumberFilter{ $titles[$x] }
							  {"Transitions"}++;
						}
					}
					else {
						$Sample2NumberAll{ $titles[$x] }{"Transversions"}++;
						if ( $done_filter == 0 ) {
							$Sample2NumberFilter{ $titles[$x] }
							  {"Transversions"}++;
						}
					}
					$Sample2NumberAll{ $titles[$x] }{"01Number"}++;
					if ( $done_filter == 0 ) {
						$Sample2NumberFilter{ $titles[$x] }{"01Number"}++;
						if ( $lines[0] =~ /[xX]/ ) {
							$genderCount{ $lines[0] }{"Heterozygous SNP"}
							  { $titles[$x] }++;
						}
					}
				}
				elsif ( ( $allele1 eq $allele2 )
					and ( $allele1 > 0 ) )
				{
					$Sample2NumberAll{ $titles[$x] }{"11Number"}++;
					if ( $done_filter == 0 ) {
						$Sample2NumberFilter{ $titles[$x] }{"11Number"}++;
						if ( $lines[0] =~ /[xX]/ ) {
							$genderCount{ $lines[0] }
							  {"Non-reference Homozygous SNP"}{ $titles[$x] }++;
						}
					}
				}

				if ( $done_filter == 0 ) {    #this line was kept
					                          #store X Y counts for sex check
					if ( $lines[0] =~ /[yY]/ ) {
						$genderCount{ $lines[0] }{"Total Reads"}
						  { $titles[$x] } += $allele1 + $allele2;
					}
					if ( ( $allele1 + $allele2 ) > 0 )
					{    #count SNP for each sample in each chromosome
						$snpCount{ $lines[0] }{ $titles[$x] }++;
					}
					if ( ( $allele1 + $allele2 ) > 0 )
					{    #count SNP NUC for each sample
						$snpNucCount{ $SNPs[$allele2] }{ $titles[$x] }++;
					}
					for ( my $y = ( $x + 1 ) ; $y < ( 9 + $sampleSize ) ; $y++ )
					{
						if ( !(defined $lines[$y]) or $lines[$y] =~ /^\./ ) { next; }
						else {
							my @result = &caculate_ratio(
								$lines[$x], $lines[$y],
								$method,    \%formatPos
							);
							$selected{ $titles[$x] }{ $titles[$y] } +=
							  $result[0];
							$total{ $titles[$x] }{ $titles[$y] } += $result[1];
							$selected{ $titles[$y] }{ $titles[$x] } +=
							  $result[2];
							$total{ $titles[$y] }{ $titles[$x] } += $result[3];
							$selected2{ $titles[$x] }{ $titles[$y] } +=
							  $result[4];
							$total2{ $titles[$x] }{ $titles[$y] } += $result[5];
						}
					}
				}
			}
		}
		close(READ);

		print RESULT3
"Sample\tTransitions:Transversions (Based on all SNPs)\tHeterozygous:Non-reference homozygous (Based on all SNPs)\tTransitions:Transversions (After filter)\tHeterozygous:Non-reference homozygous (After filter)\n";
		foreach my $sample ( sort keys %Sample2NumberAll ) {
			my $TTRatioAll = &myDivide(
				$Sample2NumberAll{$sample}{"Transitions"},
				$Sample2NumberAll{$sample}{"Transversions"}
			);
			my $TTRatioFilter = &myDivide(
				$Sample2NumberFilter{$sample}{"Transitions"},
				$Sample2NumberFilter{$sample}{"Transversions"}
			);
			my $Ratio0111All = &myDivide(
				$Sample2NumberAll{$sample}{"01Number"},
				$Sample2NumberAll{$sample}{"11Number"}
			);
			my $Ratio0111Filter = &myDivide(
				$Sample2NumberFilter{$sample}{"01Number"},
				$Sample2NumberFilter{$sample}{"11Number"}
			);

			printf RESULT3
"$sample\t$TTRatioAll\t$Ratio0111All\t$TTRatioFilter\t$Ratio0111Filter\n";
		}

		print RESULT4
"FileTitle\tSampleA\tSampleB\tCountB2A\tCountA\tHeterozygous Consistency (CountB2A:CountA)\tCountA2B\tCountB\tHeterozygous Consistency (CountA2B:CountB)\tOverall Consistent Genotypes\tOverall Overlapped Genotypes\tOverall Consistency\n";

		my $fileName = basename($VCF);
		foreach my $sample1 ( sort keys %selected2 ) {
			foreach my $sample2 ( sort keys %{ $selected2{$sample1} } ) {
				my ( $ratio1, $ratio2, $ratio3 ) = ( 0, 0, 0 );
				$ratio1 = myDivide( $selected{$sample1}{$sample2},
					$total{$sample1}{$sample2} );
				$ratio2 = myDivide( $selected{$sample2}{$sample1},
					$total{$sample2}{$sample1} );
				$ratio3 = myDivide(
					$selected2{$sample1}{$sample2},
					$total2{$sample1}{$sample2}
				);

				print RESULT4
"$fileName\t$sample1\t$sample2\t$selected{$sample1}{$sample2}\t$total{$sample1}{$sample2}\t$ratio1\t$selected{$sample2}{$sample1}\t$total{$sample2}{$sample1}\t$ratio2\t$selected2{$sample1}{$sample2}\t$total2{$sample1}{$sample2}\t$ratio3";

				print RESULT4 "\n";
			}
		}

		#export snp count
		my $temp = ( sort keys %snpCount )[0];
		print RESULT5 "Chromesome";
		foreach my $sample ( sort keys %{ $snpCount{$temp} } ) {
			print RESULT5 "\t$sample";
		}
		print RESULT5 "\n";
		foreach my $chrom ( sort keys %snpCount ) {
			print RESULT5 "$chrom";
			foreach my $sample ( sort keys %{ $snpCount{$chrom} } ) {
				if ( exists $snpCount{$chrom}{$sample} ) {
					print RESULT5 "\t$snpCount{$chrom}{$sample}";
				}
				else {
					print RESULT5 "\t0";
				}
			}
			print RESULT5 "\n";
		}

		#export information for sex check
		#		$genderCount{ "chro" }{"11Number"}{ "sample" }++;
		$temp = ( sort keys %snpCount )[0];
		print RESULT7 "Chromosome";
		foreach my $sample ( sort keys %{ $snpCount{$temp} } ) {
			print RESULT7 "\t$sample";
		}
		print RESULT7 "\n";
		foreach my $chrom ( sort keys %genderCount ) {
			foreach my $type ( sort keys %{ $genderCount{$chrom} } ) {
				print RESULT7 "$chrom: $type";
				foreach my $sample ( sort keys %{ $snpCount{$temp} } ) {
					if ( exists $genderCount{$chrom}{$type}{$sample} ) {
						print RESULT7 "\t$genderCount{$chrom}{$type}{$sample}";
					}
					else {
						print RESULT7 "\t0";
					}
				}
				print RESULT7 "\n";
			}
		}
		close(RESULT7);

		$temp = ( sort keys %snpCount )[0];
		print RESULT8 "Chromosome";
		foreach my $sample ( sort keys %{ $snpCount{$temp} } ) {
			print RESULT8 "\t$sample";
		}
		print RESULT8 "\n";
		foreach my $nuc ( sort keys %snpNucCount ) {
			print RESULT8 "$nuc";
			foreach my $sample ( sort keys %{ $snpCount{$temp} } ) {
				if ( exists $snpNucCount{$nuc}{$sample} ) {
					print RESULT8 "\t$snpNucCount{$nuc}{$sample}";
				}
				else {
					print RESULT8 "\t0";
				}
			}
			print RESULT8 "\n";
		}
		close(RESULT8);

		#do annovar
		if ($doAnnovar) {
			pInfo( "do ANNOVAR annotation", $logRef );
			system(
"$annovarConvert -format vcf4 $annovarResultDir$filename.pass -includeinfo > $annovarResultDir$filename.pass.avinput"
			);
			system(
"$annovarBin $annovarResultDir$filename.pass.avinput $annovarDb $annovarOption  --outfile $annovarResultDir$filename.pass.avinput.annovar"
			);
		}

		#end comment
	}

	#plot by R
	my $annovarBuildver = "";
	if ( $annovarOption =~ /-buildver (\S+) / ) {
		$annovarBuildver = $1;
	}
	my $Rsource = dirname($0) . $rSourceLocation;
	my $rResult = system(
"cat $Rsource | $RBin --vanilla --slave --args $resultDir/$filename $annovarBuildver > $resultDir/vcfResult/vcfSummary.rLog"
	);
	pInfo( "Finish vcf summary!", $logRef );
	return ($rResult);
}

sub filter {
	my $input     = $_[0];
	my %cfgFilter = %{ $_[1] };

	my %filter;

	my @filter1 = ( split /;/, $input );
	foreach my $filter (@filter1) {
		if ( $filter =~ /=/ ) {
			my @filter2 = ( split /=/, $filter );
			$filter{ $filter2[0] } = $filter2[1];
		}
	}

	#if keep this line:
	foreach my $key ( sort keys %cfgFilter ) {
		if ( exists( $filter{$key} ) ) {
			my @filterValues = ( split /,/, $filter{$key} );
			foreach my $filterValue (@filterValues) {
				if ( eval( $filterValue . $cfgFilter{$key} ) ) {
					return (0);
				}
			}
		}
	}
	return (1);
}

sub caculate_ratio {
	my $sample1   = $_[0];
	my $sample2   = $_[1];
	my $method    = $_[2];
	my %formatPos = %{ $_[3] };

	my (
		$a2bSelected, $b2aSelected, $abSelected,
		$a2bTotal,    $b2aTotal,    $abTotal
	) = ( 0, 0, 0, 0, 0, 0 );
	my $fenzi = 0;

	my @lines1 = ( split /:/, $sample1 );
	my @lines2 = ( split /:/, $sample2 );
	my $GTpos  = $formatPos{"GT"};
	my ( $sample1GTRef, $sample1GTAlt ) = ( split /\//, $lines1[$GTpos] );
	my ( $sample2GTRef, $sample2GTAlt ) = ( split /\//, $lines2[$GTpos] );

	my $Sample1RDreads = 0;
	my $Sample1ADreads = 0;
	my $Sample2RDreads = 0;
	my $Sample2ADreads = 0;
	if ( exists( $formatPos{"RD"} ) ) {  #RD and AD reads as RDreads and ADreads
		my $ADpos = $formatPos{"AD"};
		my $RDpos = $formatPos{"RD"};
		$Sample1RDreads = $lines1[$RDpos];
		$Sample1ADreads = $lines1[$ADpos];
		$Sample2RDreads = $lines2[$RDpos];
		$Sample2ADreads = $lines2[$ADpos];
	}
	else {    #RD and AD reads as 0Reads,1Reads,2Reads
		my $ADpos  = $formatPos{"AD"};
		my @reads1 = ( split /,/, $lines1[$ADpos] );
		my @reads2 = ( split /,/, $lines2[$ADpos] );
		$Sample1RDreads = $reads1[$sample1GTRef];
		$Sample1ADreads = $reads1[$sample1GTAlt];
		$Sample2RDreads = $reads2[$sample2GTRef];
		$Sample2ADreads = $reads2[$sample2GTAlt];
	}

	#deepth <10
	my $sum1 = 0;
	my $sum2 = 0;
	if ( $Sample1RDreads eq '.' or $Sample2RDreads eq '.' ) {
		return ( 0, 0, 0, 0, 0, 0 );
	}
	if ( $sample1GTRef eq $sample1GTAlt ) {
		$sum1 = $Sample1RDreads;
	}
	else {
		$sum1 = $Sample1RDreads + $Sample1ADreads;
	}
	if ( $sample2GTRef eq $sample2GTAlt ) {
		$sum2 = $Sample2RDreads;
	}
	else {
		$sum2 = $Sample2RDreads + $Sample2ADreads;
	}
	if ( $sum1 < 10 or $sum2 < 10 ) {
		(
			$a2bSelected, $a2bTotal,   $b2aSelected,
			$b2aTotal,    $abSelected, $abTotal
		) = ( 0, 0, 0, 0, 0, 0 );
	}
	else {
		$abTotal = 1;
		if ( $method == 1 ) {
			if (    ( $sample1GTRef eq $sample2GTRef )
				and ( $sample1GTAlt eq $sample2GTAlt ) )
			{
				$fenzi = 1;
			}
		}
		elsif ( $method == 2 ) {
			my $ADpos  = $formatPos{"AD"};
			my @reads1 = ( split /,/, $lines1[$ADpos] );
			my @reads2 = ( split /,/, $lines2[$ADpos] );
			if (
				$sample1GTRef eq $sample2GTRef
				and (
					( $sample1GTAlt eq $sample2GTAlt )
					or (

						#						(
						#							$lines1[ 2 + $lines1[1] ] *
						#							$lines2[ 2 + $lines1[1] ]
						#						) > 0
						#						and ( $lines1[ 2 + $lines2[1] ] *
						#							$lines2[ 2 + $lines2[1] ] ) > 0
						( $reads1[$sample1GTAlt] * $reads2[$sample1GTAlt] ) > 0
						and ( $reads1[$sample2GTAlt] * $reads2[$sample2GTAlt] )
						> 0
					)
				)
			  )
			{
				$fenzi = 1;
			}
		}
		if ( $fenzi == 1 ) {
			$abSelected = 1;
		}
		else { $abSelected = 0; }

		#Samples A
		if ( $sample1GTRef ne $sample1GTAlt ) {
			$a2bTotal = 1;
			if ( $fenzi == 1 ) {
				$a2bSelected = 1;
			}
			else {
				$a2bSelected = 0;
			}
		}
		else {
			$a2bTotal    = 0;
			$a2bSelected = 0;
		}

		#Samples B
		if ( $sample2GTRef ne $sample2GTAlt ) {
			$b2aTotal = 1;
			if ( $fenzi == 1 ) {
				$b2aSelected = 1;
			}
			else {
				$b2aSelected = 0;
			}
		}
		else {
			$b2aTotal    = 0;
			$b2aSelected = 0;
		}
	}

	return (
		$a2bSelected, $a2bTotal,   $b2aSelected,
		$b2aTotal,    $abSelected, $abTotal
	);
}

#sub pInfo {
#	my $s       = shift;
#	my $logFile = shift;
#	print "[", scalar(localtime), "] $s\n";
#	print $logFile "[", scalar(localtime), "] $s\n";
#}
sub pInfo {
	my $s = shift;
	print "[", scalar(localtime), "] $s\n";
	push @{ $_[1] }, "[" . scalar(localtime) . "] $s\n";
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
