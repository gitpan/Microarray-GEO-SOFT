package Microarray::GEO::SOFT::GSM;

# 解析SOFT文件
# 读取GSM部分，只读取文件中从读取位置开始的第一个GSM记录

# parse SOFT file
# get the GSM part, only the first record of GSM from the current position to read

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
	
	$self->set_class("GSM");
	
	return $self;
	
}

sub parse {

	my $self = shift;
	
	my $fh;
	if(! is_glob_ref($self->{file})) {
	
		open F, $self->{file} or die "cannot open $self->{file}.\n";
		$fh = \*F;
		
		$self->{file} = $fh;
	}
	
	$self->parse_sample($self->{file});
	
	return 1;
}

sub parse_sample {

	my $self = shift;

	my $fh = shift;
	
	my $accession;
	my $title;
	my $platform;
	my $field_explain;
	my $table_colnames = [];
	my $table_rownames = [];
	my $table_matrix;
	
	my $_TMP_SOFT_DIR = $self->soft_dir;
	
	while(my $line = <$fh>) {
		
		chomp $line;
		if($line =~/^!Sample_geo_accession = (GSM\d+)$/) {
			$accession = $1;
			;
		}
		
		if($line =~/^!Sample_title = (.*?)$/) {
			$title = $1;
		}
		
		if($line =~/^!Sample_platform_id = (GPL\d+)$/) {
			$platform = $1;
		}
		
		if($line =~/^#(.*?) = (.*?)$/) {
			$field_explain->{$1} = $2;
		}
		
		if($line =~/^!sample_table_begin$/) {
			
			$line = <$fh>;
			chomp $line;
			
			@$table_colnames = split "\t", $line;
			
			while($line = <$fh>) {
			
				if($line =~/^!sample_table_end$/) {
					last;
				}
			
				chomp $line;
				my @tmp = split "\t", $line;
				
				my $uid = shift(@tmp);
				
				push(@$table_rownames, $uid);
				# 第一列为VALUE
				push(@$table_matrix, [$tmp[0]]);
			}
			
			
		}
		if($line =~/^!sample_table_end$/) {
			last;
		}
		
	}
	
	my $n_row = len($table_rownames);
	my $n_col = len($table_colnames);
	
	print "Sample info:\n";
	print "Accession: $accession\n";
	print "Platform: $platform\n";
	print "Title: $title\n";
	print "Rows: $n_row\n";
	print "Columns: $n_col\n";
	print "Sorting UIDs...\n";
	print "\n";
	
	my $table_rownames_sorted = sort_array($table_rownames, sub {$_[0] cmp $_[1]});
	my $table_rownames_sorted_index = order($table_rownames, sub {$_[0] cmp $_[1]});
	my $table_matrix_sorted = subset($table_matrix, $table_rownames_sorted_index);
		
	open OUT, ">$_TMP_SOFT_DIR/$accession.tab";
	for(my $i = 0; $i < len($table_matrix_sorted); $i ++) {
		print OUT join "\t", @{$table_matrix_sorted->[$i]};
		print OUT "\n";
	}
	close OUT;
	
	$self->set_meta($accession, $title, $platform, $field_explain);
	$self->set_table($table_rownames_sorted, $table_colnames, undef);
	
	
	return $self;

}


__END__

=pod

=head1 NAME

Microarray::GEO::SOFT::GSM - GEO sample data class

=head1 SYNOPSIS

  use Microarray::GEO::SOFT:
  my $soft = Microarray::GEO::SOFT->new;
  $soft->download("GSE35505");
  
  my $gse = $soft->parse;
  my $gsm = $gse->list("GSM")->[0];
  
  # the meta information
  $gsm->meta;
  $gsm->platform;
  $gsm->title;
  $gsm->field;
  $gsm->accession;
  
  # the sample data is a matrix (in fact it is a vector)
  $gsm->matrix;
  # the names for each column
  $gsm->colnames;
  $ the names for each row, it is the primary id for rows
  $gsm->rownames;

=head1 DESCRIPTION

This module retrieves sample information from microarray data. The module is used
within parsing GSE data.

=head2 Subroutines

=over 4

=item C<new("file" = $file)>

Initial a GSM class object. The only argument is the path of the microarray data in SOFT format
or a file handle that has been openned.

=item C<$gsm-E<gt>parse>

Retrieve sample information from microarray data. The sample data in SOFT format
is alawys a table

=item C<$gsm-E<gt>meta>

Get meta information

=item C<$gsm-E<gt>platform>

Get accession number of the platform

=item C<$gsm-E<gt>title>

Title of the platform record

=item C<$gsm-E<gt>field>

Description of each field in the data matrix

=item C<$gsm-E<gt>accession>

Accession number for the sample

=back

=head1 AUTHOR

Zuguang Gu E<lt>jokergoo@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2012 by Zuguang Gu

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.1 or,
at your option, any later version of Perl 5 you may have available.

=head1 SEE ALSO

L<Microarray::GEO::SOFT>

=cut

