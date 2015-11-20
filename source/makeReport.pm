package source::makeReport;

use strict;
use warnings;
use Exporter;
use Scalar::Util qw(looks_like_number);

our @ISA    = qw(Exporter);
our @EXPORT = qw(build_template file2table file2text dir2list cstmrow_filter);

1;

sub build_template {
	my $templateName = $_[0];
	my %hash         = %{ $_[1] };

	my $headEnd    = 0;
	my $summaryTop = '<div class="summary">
    <h2>Table of Contents</h2>
        <ul>';
	my $summaryBottom = '        </ul>
</div>';
	my @head;
	my @content;
	my @summary;
	open TEMP, "<$templateName" or die "can't open $templateName!";

	while (<TEMP>) {
		chomp;
		if ( $headEnd == 0 ) {
			push @head, $_;
		}
		else {
			push @content, $_;
			if (/TMPL_LOOP NAME=(MAKETABLE\d+)/) {
				if ( exists $hash{$1} ) {
					foreach my $tableTitle ( @{ ${ $hash{$1} }[0] } )
					{
						push @content,
						    '      <td align="center"><!-- TMPL_VAR NAME=\''
						  . $tableTitle
						  . '\' --></td>';
					}
					shift @{ $hash{$1} };
				}
			}
		}
		if (/<\/h1>/) {
			$headEnd = 1;
		}
		if (/<div><h2 id="(\w+)">([\s\w]+)/) {
			push @summary, '<li><a href=#' . $1 . '>' . $2 . '</a></li>';
		}
	}
	open TEMPRESULT, ">$templateName.temp" or die "can't write $!";
	print TEMPRESULT ( join( "\n", @head ) );
	print TEMPRESULT $summaryTop . "\n";
	print TEMPRESULT ( join( "\n", @summary ) );
	print TEMPRESULT $summaryBottom . "\n";
	print TEMPRESULT ( join( "\n", @content ) );
	close(TEMPRESULT);

	my $template = HTML::Template->new(
		filename          => "$templateName.temp",
		loop_context_vars => 1,
		die_on_bad_params => 0,
		filter            => \&cstmrow_filter
	);
	foreach my $key ( keys %hash ) {
		$template->param( $key => $hash{$key} );
	}
	unlink "$templateName.temp";
	return ($template);
}

sub file2text {
	my $file        = $_[0];
	my $result="";
	open READ, "<$file" or die "can't find $file\n";
	while (<READ>) {
		chomp;
		$result="$result$_\n";
	}
	return($result);
}

sub file2table {
	my $file        = $_[0];
	my $recordTitle = $_[2];

	my $firstLine = 0;
	my $rows;
	my @title;

	open READ, "<$file" or die "can't find $file\n";
	while (<READ>) {
		chomp;
		my @content = ( split /\t/, $_ );
		if ( $firstLine == 0 ) {
			if ( defined( $_[1] ) and $_[1] ne "" ) {
				@title = @content[ @{ $_[1] } ];
			}
			else {
				@title = @content;
			}

			#			print (join "\n",@title);
			if (defined($recordTitle) and $recordTitle==1) {
				push @{$rows}, \@title;
			}
			$firstLine++;
		}
		my $temp = {};
		if ( defined( $_[1] ) and $_[1] ne "" ) {
			my @temp=&arry2Number(@content[ @{ $_[1] } ]);
			@$temp{@title} = @temp;
		}
		else {
			my @temp=&arry2Number(@content);
			@$temp{@title} = @temp;
		}
		push @{$rows}, $temp;
	}
	close(READ);
	return ($rows);
}

sub dir2list {
	my $dir1     = $_[0];    #$_[0]
	my $dir2     = $_[1];    #"/scorefigure/"
	my $pattern  = $_[2];    #".png"
	my $tmplName = $_[3];    #"FIGURE2"
	my $herf     = $_[4];    #"FIGURE2"

	my $list;
	opendir( DIR, "$dir1$dir2" ) or die $!;
	while ( my $file = readdir(DIR) ) {

		#		print "$dir1$dir2/$file\n";
		if ( $file =~ /$pattern/ ) {
			my $temp = {};
			${$temp}{$tmplName} = ".$dir2/$file";
			if ($herf) {
				${$temp}{ $tmplName . "_NAME" } = "$file";
			}
			push @{$list}, $temp;
		}
	}
	close(DIR);
	return ($list);
}

sub cstmrow_filter {
	my $text_ref = shift;

	#no first, end with no space, don't match if with first
	$$text_ref =~ s/<CSTM_ROW\s+EVEN=(.+)\s+ODD=(\S+)\s*>
                 /<TMPL_IF NAME=__odd__>
                    <tr class="$2">
                  <TMPL_ELSE>
                    <tr class="$1">
                  <\/TMPL_IF>
                 /gx;

	#with first
	$$text_ref =~ s/<CSTM_ROW\s+EVEN=(.+)\s+ODD=(.+)\s+FIRST=(.+)\s*>
                 /<TMPL_IF NAME=__first__>
                 <tr class="$3">
                  <TMPL_ELSE>
                 <TMPL_IF NAME=__odd__>
                    <tr class="$1">
                  <TMPL_ELSE>
                    <tr class="$2">
                  <\/TMPL_IF>
                  <\/TMPL_IF>
                 /gx;
}

sub arry2Number {
	my @data     = @_;
	my @result;
	foreach my $data (@data) {
		if (looks_like_number($data)) {
			my $format;
			if ($data=~/\./) {
				$format="%.2f";
			} else {
				$format="%d";
			}
			my $temp=sprintf($format,$data);
			push @result,$temp;
		} else {
			push @result,$data;
		}
	}
	return(@result);
}
