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
my $logRef=\@main::log;

sub vcfSummary {

	my ( $VCF, $config ) = @_;
	my $vcfCfgFile   = $config->{'vcfCfgFile'};
	my $method    = $config->{'method'};
	my $resultDir = $config->{'resultDir'};
#	my $logFile   = $config->{'log'};
	my $rePlot    = $config->{'rePlot'};

	my $RBin           = $config->{'RBin'};
	my $annovarBin     = $config->{'annovarBin'};
	my $annovarConvert = $config->{'annovarConvert'};
	my $annovarOption  = $config->{'annovarOption'};
	my $annovarDb      = $config->{'annovarDb'};
	
	my $showHelp     = $config->{'showHelp'};

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
	-a	annotation directory	Optional. Directory of annovar database.

	-h	                        Optional. Show this information.
	
For more information, please refer to the readme file in QC3 directory. Or visit the QC3 website at https://github.com/slzhao/QC3

";

	die "\n$usage" if ($showHelp);
	die "Can't find the vcf QC config file\n$usage" if (!( -e $vcfCfgFile ));
	die "Method(-s) should be 1 or 2\n$usage" if ($method != 1 and $method != 2);

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

		while (<READ>) {    #read title and ID list
			chomp;
			s/\r//g;
			if (/^##INFO=<ID=(\w+),/) {    #all IDs
				$IDList{$1} = "";
			}
			elsif (/^##/) {
				next;
			}
			elsif (/^#/) {                 #title
				@titles = ( split /\t/, $_ );
				$sampleSize = scalar(@titles) - 9;
				pInfo( "The VCF file has $sampleSize samples", $logRef );
				last;
			}
		}
		close(READ);

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
			if (   $lines[6] eq "LowQual"
				or &filter( $lines[7], \%cfgFilter ) == 0 )
			{
				$done_filter = 1;
			}
			else {
				$done_filter = 0;
				print RESULT2 "$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]";
				print RESULT6 "$_\n";
			}
			print RESULT1 "$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]";

			my @ID1 = ( split /;/, $lines[7] );
			my %IDNumber;
			foreach my $ID (@ID1) {
				my @ID2 = ( split /=/, $ID );
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
				if ( $lines[$x] eq './.' ) { next; }
				elsif ( $lines[$x] =~ '0/1' ) {
					if (   ( $lines[3] . $lines[4] eq "AG" )
						or ( $lines[3] . $lines[4] eq "GA" )
						or ( $lines[3] . $lines[4] eq "CT" )
						or ( $lines[3] . $lines[4] eq "TC" ) )
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
					}
				}
				elsif ( $lines[$x] =~ '1/1' ) {
					$Sample2NumberAll{ $titles[$x] }{"11Number"}++;
					if ( $done_filter == 0 ) {
						$Sample2NumberFilter{ $titles[$x] }{"11Number"}++;
					}
				}
				if ( $done_filter == 0 ) {    #this line was kept
					if ( $lines[$x] =~ '0/1' or $lines[$x] =~ '1/1' )
					{    #count SNP for each sample in each chromosome
						$snpCount{ $lines[0] }{ $titles[$x] }++;
					}
					for ( my $y = ( $x + 1 ) ; $y < ( 9 + $sampleSize ) ; $y++ )
					{
						if ( $lines[$y] eq './.' ) { next; }
						else {
							my @result =
							  &caculate_ratio( $lines[$x], $lines[$y],
								$method );
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
#			printf RESULT3 (
#				"%s\t%.3f\t%.3f\t%.3f\t%.3f\n",
#				$sample, $TTRatioAll, $Ratio0111All, $TTRatioFilter,
#				$Ratio0111Filter
#			);
		}

		print RESULT4
"FileTitle\tSampleA\tSampleB\tCountB2A\tCountA\tHeterozygous Consistency (CountB2A:CountA)\tCountA2B\tCountB\tHeterozygous Consistency (CountA2B:CountB)\tOverall Consistent Genotypes\tOverall Overlapped Genotypes\tOverall Consistency\n";
#		foreach my $change ( sort keys %changes ) {
#			print RESULT4 "\t$change";
#		}
#		print RESULT4 "\n";

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

#				printf RESULT4 (
#					"%s\t%s\t%s\t%d\t%d\t%.3f\t%d\t%d\t%.3f\t%d\t%d\t%.3f",
#					$fileName,
#					$sample1,
#					$sample2,
#					$selected{$sample1}{$sample2},
#					$total{$sample1}{$sample2},
#					$ratio1,
#					$selected{$sample2}{$sample1},
#					$total{$sample2}{$sample1},
#					$ratio2,
#					$selected2{$sample1}{$sample2},
#					$total2{$sample1}{$sample2},
#					$ratio3
#				);
				print RESULT4 
					"$fileName\t$sample1\t$sample2\t$selected{$sample1}{$sample2}\t$total{$sample1}{$sample2}\t$ratio1\t$selected{$sample2}{$sample1}\t$total{$sample2}{$sample1}\t$ratio2\t$selected2{$sample1}{$sample2}\t$total2{$sample1}{$sample2}\t$ratio3";

#				foreach my $change ( sort keys %changes ) {
#					if ( exists( $changNumber{$sample1}{$sample2}{$change} ) ) {
#						print RESULT4
#						  "\t$changNumber{$sample1}{$sample2}{$change}";
#					}
#					else {
#						print RESULT4 "\t0";
#					}
#				}
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
				if (exists $snpCount{$chrom}{$sample}) {
					print RESULT5 "\t$snpCount{$chrom}{$sample}";
				} else {
					print RESULT5 "\t0";
				}
			}
			print RESULT5 "\n";
		}

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
		my @filter2 = ( split /=/, $filter );
		$filter{ $filter2[0] } = $filter2[1];
	}

	#if keep this line:
	foreach my $key ( sort keys %cfgFilter ) {
		if ( exists( $filter{$key} )
			and eval( $filter{$key} . $cfgFilter{$key} ) )
		{

			#			print "$filter{$key} $cfgFilter{$key}\n";
			return (0);
		}
	}
	return (1);
}

sub caculate_ratio {
	my $sample1 = $_[0];
	my $sample2 = $_[1];
	my $method  = $_[2];

	my (
		$a2bSelected, $b2aSelected, $abSelected,
		$a2bTotal,    $b2aTotal,    $abTotal
	) = ( 0, 0, 0, 0, 0, 0 );
	my $fenzi = 0;

	my @lines1 = ( split /;|:|,|\//, $sample1 );
	my @lines2 = ( split /;|:|,|\//, $sample2 );

	#deepth <10
	if ( ( $lines1[2] + $lines1[3] ) < 10 or ( $lines2[2] + $lines2[3] ) < 10 )
	{
		(
			$a2bSelected, $a2bTotal,   $b2aSelected,
			$b2aTotal,    $abSelected, $abTotal
		) = ( 0, 0, 0, 0, 0, 0 );
	}
	else {
		$abTotal = 1;
		if ( $method == 1 ) {
			if ( ( $lines1[0] == $lines2[0] ) and ( $lines1[1] == $lines2[1] ) )
			{
				$fenzi = 1;
			}
		}
		elsif ( $method == 2 ) {
			if (
				$lines1[0] == $lines2[0]
				and (  ( $lines1[1] == $lines2[1] )
					or ( $lines1[3] * $lines2[3] > 0 ) )
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
		if ( ( $lines1[0] + $lines1[1] ) == 1 ) {
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
		if ( ( $lines2[0] + $lines2[1] ) == 1 ) {
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
	push @{$_[1]}, "[".scalar(localtime)."] $s\n";
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
