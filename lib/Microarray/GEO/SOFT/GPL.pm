package Microarray::GEO::SOFT::GPL;


use List::Vectorize;
use strict;

require Microarray::GEO::SOFT;
our @ISA = ("Microarray::GEO::SOFT");

$| = 1;


1;


sub new {

	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = { "file" => "",
	             @_ };
	bless($self, $class);
	
	$self->set_class("GPL");
	
	return $self;
	
}

# 在soft文件中解析platform部分
sub parse {

	my $self = shift;
	
	my $fh;
	if(ref($self->{file}) ne "GLOB") {
	
		open F, $self->{file} or die "cannot open $self->{file}.\n";
		$fh = \*F;
		
		$self->{file} = $fh;
	}
	
	$self->parse_platform($self->{file});
	
	return 1;
}

sub parse_platform {

	my $self = shift;

	my $fh = shift;
	
	my $accession;
	my $title;
	my $field_explain;
	my $table_colnames = [];
	my $table_rownames = [];
	my $table_matrix;
	
	while(my $line = <$fh>) {
		
		chomp $line;
		if($line =~/^!Platform_geo_accession = (GPL\d+)$/) {
			$accession = $1;
		}
		
		if($line =~/^!Platform_title = (.*?)$/) {
			$title = $1;
		}
		
		if($line =~/^#(.*?) = (.*?)$/) {
			$field_explain->{$1} = $2;
		}
		
		if($line =~/^!platform_table_begin$/) {
			
			$line = <$fh>;
			chomp $line;
			
			@$table_colnames = split "\t", $line;
			shift(@$table_colnames);
			
			while($line = <$fh>) {
			
				if($line =~/^!platform_table_end$/) {
					last;
				}
			
				chomp $line;
				my @tmp = split "\t", $line;
				
				my $uid = shift(@tmp);
				
				push(@$table_rownames, $uid);
				push(@$table_matrix, [@tmp]);
				
			}
			
			
		}
		if($line =~/^!platform_table_end$/) {
			last;
		}
		
	}
	
	my $n_row = len($table_rownames);
	my $n_col = len($table_colnames);
	
	my $platform = $accession;
	
	print "Platform info:\n";
	print "Accession: $accession\n";
	print "Platform: $platform\n";
	print "Title: $title\n";
	print "Rows: $n_row\n";
	print "Columns: $n_col\n";
	print "Sorting UIDs...\n";
	print "\n";
	
	# ID列从大到小排序
	my $table_rownames_sorted = sort_array($table_rownames, sub {$_[0] cmp $_[1]});
	my $table_rownames_sorted_index = order($table_rownames, sub {$_[0] cmp $_[1]});
	my $table_matrix_sorted = subset($table_matrix, $table_rownames_sorted_index);
	
	$self->set_meta($accession, $title, $platform, $field_explain);
	$self->set_table($table_rownames_sorted, $table_colnames, $table_matrix_sorted);
	
	
	return $self;
}

sub mapping {

	my $self = shift;
	my $to_id = shift;
	
	my $mapping;
	
	my $to_index;
	my $colnames = $self->colnames;
	for(my $i = 0; $i < len($colnames); $i ++) {
		if($colnames->[$i] eq $to_id) {
			$to_index = $i;
			last;
		}
	}
	
	my $mat = $self->matrix;
	for(my $i = 0; $i < len($mat); $i ++) {
		push(@$mapping, $mat->[$i]->[$to_index]);
	}

	return $mapping;
	
}


__END__

=pod

=head1 NAME

Microarray::GEO::SOFT::GPL - GEO platform data class

=head1 SYNOPSIS

  use Microarray::GEO::SOFT:
  my $soft = Microarray::GEO::SOFT->new("file" => "GPL15181.soft");
  
  # or you can download from GEO website
  my $soft = Microarray::GEO::SOFT->new;
  $soft->download("GPL15181");
  
  # $gpl is a Microarray::GEO::SOFT::GPL class object
  my $gpl = $soft->parse;
  
  # the meta information
  $gpl->meta;
  $gpl->platform;
  $gpl->title;
  $gpl->field;
  $gpl->accession;
  
  # the platform data is a matrix
  $gpl->matrix;
  # the names for each column
  $gpl->colnames;
  $ the names for each row, it is the primary id for rows
  $gpl->rownames;
  
  # we want to get other ID for microarray data
  my $other_id = $gpl->mapping("miRNA_ID");
  # or
  my $other_id = $gpl->mapping($gpl->colnames[1]);

=head1 DESCRIPTION

This module retrieves platform information from microarray data. It is not usually
to get platform information from GPL record alone since a GPL record is such a huge
file. It is common to get platform information from GSE record.

=head2 Subroutines

=over 4

=item C<new("file" = $file)>

Initial a GPL class object. The only argument is the microarray data in SOFT format
or a file handle that has been openned.

=item C<$gpl-E<gt>parse>

Retrieve platform information from microarray data. The platform data in SOFT format
is alawys a table

=item C<$gpl-E<gt>meta>

Get meta information

=item C<$gpl-E<gt>platform>

Get accession number of the platform

=item C<$gpl-E<gt>title>

Title of the platform record

=item C<$gpl-E<gt>field>

Description of each field in the data matrix

=item C<$gpl-E<gt>accession>

Accession number for the platform

=item C<$gpl-E<gt>mapping>

get ID mappings

=back

=head1 AUTHOR

Zuguang Gu E<lt>jokergoo@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2012 by Zuguang Gu

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.1 or,
at your option, any later version of Perl 5 you may have available.

=head1 SEE ALSO

Microarray::GEO::SOFT

=cut
