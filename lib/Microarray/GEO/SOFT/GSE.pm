package Microarray::GEO::SOFT::GSE;

# 解析SOFT文件

use List::Vectorize;
use Microarray::GEO::SOFT::GPL;
use Microarray::GEO::SOFT::GSM;
use Microarray::GEO::SOFT::GDS;
use strict;

require Microarray::GEO::SOFT;
our @ISA = ("Microarray::GEO::SOFT");


$| = 1;

# 全局的文件句柄
our $fh;

1;

sub new {

	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = { "file" => "",
	             @_ };
	bless($self, $class);
	
	$self->set_class("GSE");
	
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
	
	$self->parse_series($self->{file});
	
	return 1;
}

# 解析series
sub parse_series {

	my $self = shift;

	my $fh = shift;
	my $using_sqlite = shift;
	
	my $accession;
	my $title;
	my $platform;
	
	my $series;
	my $field_explain;
	
	my $gpl_list;
	my $gsm_list;
	
	while(my $line = <$fh>) {
	
		chomp $line;
		
		if($line =~/^\^SERIES = (GSE\d+)$/) {
			$accession = $1;
		}
		
		if($line =~/^!Series_title = (.*?)$/) {
			$title = $1;
		}
		
		if($line =~/^!Series_platform_id = (GPL\d+)$/) {
			push(@$platform, $1);
		}
		
		# platfrom 信息
		elsif($line =~/^\^PLATFORM = (GPL\d+)$/) {
		
			$fh = back_to_last_line($fh, length($line));
			
			# 返回GEO::GPL对象
			my $gpl = Microarray::GEO::SOFT::GPL->new(file => $fh);
			$gpl->parse;
			
			push(@$gpl_list, $gpl);
		}
		
		elsif($line =~/^\^SAMPLE = (GSM\d+)$/) {
		
			$fh = back_to_last_line($fh, length($line));
			
			# 返回SOFT::GSM对象
			my $gsm = Microarray::GEO::SOFT::GSM->new(file => $fh);
			$gsm->parse;
			
			push(@$gsm_list, $gsm);
			
		}
	}
	
	my $n_platform = len($gpl_list);
	my $n_sample = len($gsm_list);
	
	print "Series info:\n";
	print "Accession: $accession\n";
	print "Title:$title\n";
	print "Platforms: $n_platform\n";
	print "Samples: $n_sample\n";
	print "\n";
	
	$self->set_meta($accession, $title, $platform, $field_explain);
	$self->set_list($gpl_list, "GPL");
	$self->set_list($gsm_list, "GSM");
	
	return $self;
}

# 读取series信息
sub read_series {
	
	my $fh = shift;
	
	my $accession;
	my $sample_id = [];
	my $platform_id = [];
	
	while(my $line = <$fh>) {
	
		chomp $line;
		
		# 当series记录结束
		if($line =~/^\^/) {
		
			$fh = back_to_last_line($fh, length($line));
			last;
		}
		
		if($line =~/^!Series_geo_accession = (GSE\d+)$/) {
			$accession = $1;
		}
		elsif($line =~/^!Series_sample_id = (GSM\d+)$/) {
			push(@$sample_id, $1);
		}
		elsif($line =~/^!Series_platform_id = (GPL\d+)$/) {
			push(@$platform_id, $1);
		}
	}
	
	my $res;
	$res->{accession} = $accession;
	$res->{platform_id} = $platform_id;
	$res->{sample_id} = $sample_id;
	
	return $res;
}



# 回到上一行的开头
sub back_to_last_line {
	
	my $fh = shift;
	my $current_line_length = shift;
	
	my $position = tell($fh);
	seek($fh, $position - $current_line_length - 2, 0);
	
	return $fh;
}

# 把相同platform合并为matrix
# 返回一个GDS的数组
sub merge {
		
	my $self = shift;
	
	my $gse_platform = $self->platform;
	my $gds;
	
	for(my $i = 0; $i < len($gse_platform); $i ++) {
	
		my $sample_list = $self->list("GSM");
		
		# $s 是
		my $s = subset($sample_list, sub {$_[0]->platform eq $gse_platform->[$i]} );
		
		my $g = Microarray::GEO::SOFT::GDS->new;
		$g->merge_gsm($s);

		push(@$gds, $g);
	}
	
	return $gds;
	
}


# 把若干个相同platform的GSM合并为一个GDS
# 虽然不是GEO上定义的GDS (manually assembled)
# 但是数据格式相同
sub merge_gsm {

	my $self = shift;
	
	my $gsm_list = shift;
	
	# 检查gpl号是否都一样
	my $gpl_list = sapply($gsm_list, sub {$_[0]->platform});
	if(len(unique($gpl_list)) != 1) {
		die "Platform should be same\n";
	}
	
	my $r = sample(c(["0".."9"], ["A".."Z"]), 20, "replace" => 1);
	my $accession = "GDS_".(join "", @$r);
	my $title = "gsm_merged";
	my $platform = $gpl_list->[0];
	my $field_explain;
	my $table_colnames;
	
	for(my $i = 0; $i < len($gsm_list); $i ++) {
	
		$field_explain->{$gsm_list->[$i]->accession} = $gsm_list->[$i]->title;
		$table_colnames->[$i] = $gsm_list->[$i]->accession;
		
	}

	my $table_rownames = $gsm_list->[0]->rownames;
	
	# 把数据存入文本文件中
	my $_TMP_SOFT_DIR = $self->soft_dir;
	open OUT, ">$_TMP_SOFT_DIR/$accession.tab" or die "cannot open $_TMP_SOFT_DIR/$accession.tab\n";
	
	# 首先初始化文件句柄
	my $fh;
	for(my $i = 0; $i < len($gsm_list); $i ++) {
		local *F;
		$fh->[$i] = *F;
		
		my $acc = $gsm_list->[$i]->accession;
		open F, "$_TMP_SOFT_DIR/$acc.tab" or die "cannot open $_TMP_SOFT_DIR/$acc.tab\n";
		
	}
	
	# 写入数据
	my $flag = 0;
	while(1) {
		for(my $i = 0; $i < len($gsm_list); $i ++) {
			
			my $handle = $fh->[$i];
			my $line = <$handle>;
			
			if($line) {
				chomp $line;
				my @tmp = split "\t", $line;
				
				$i == len($gsm_list) - 1 ? print OUT "$tmp[0]\n"
										 : print OUT "$tmp[0]\t";
			}
			else {
				$flag = 1;
			}
		}
		
		if($flag) {
			last;
		}
	
	}
	
	# 关闭文件句柄
	for(my $i = 0; $i < len($gsm_list); $i ++) {
		my $handle = $fh->[$i];
		close $handle;
	}
	
	close OUT;
	
	
	my $n_row = len($table_rownames);
	my $n_col = len($table_colnames);
	
	print "Merge GSM into GDS:\n";
	print "Accession: $accession\n";
	print "Platform: $platform\n";
	print "Title: $title\n";
	print "Rows: $n_row\n";
	print "Columns: $n_col\n";
	print "\n";
	
	$self->set_meta($accession, $title, $platform, $field_explain);
	$self->set_table($table_rownames, $table_colnames, undef);
	
	return $self;
}


__END__

=pod

=head1 NAME

Microarray::GEO::SOFT::GSE - GEO series data class

=head1 SYNOPSIS

  use Microarray::GEO::SOFT:
  my $soft = Microarray::GEO::SOFT->new("file" => "GSE35505.soft");
  
  # or you can download from GEO website
  my $soft = Microarray::GEO::SOFT->new;
  $soft->download("GSE35505");
  
  # $gse is a Microarray::GEO::SOFT::GSE class object
  my $gse = $soft->parse;
  
  # the meta information
  $gse->meta;
  $gse->platform;
  $gse->title;
  $gse->field;
  $gse->accession;
	
  # since a GSE can contain more than one GSM and GPL, so the GPL and GSM stored
  # in GSE is a list or array
  my $samples = $gse->list("GSM");
  my $platforms = $gse->list("GPL");
  
  # data in single GSM can be merged as matrix by platforms
  # it is a GDS class object
  my $g = $gse->merge->[0];

=head1 DESCRIPTION

This module retrieves data storing as GEO series format.

=head2 Subroutines

=over 4

=item C<new("file" = $file)>

Initial a GSE class object. The only argument is the microarray data in SOFT format
or a file handle that has been openned.

=item C<$gse-E<gt>parse>

Retrieve series information from microarray data.

=item C<$gse-E<gt>meta>

Get meta information

=item C<$gse-E<gt>platform>

Get accession number of the platform

=item C<$gse-E<gt>title>

Title of the platform record

=item C<$gse-E<gt>field>

Description of each field in the data matrix

=item C<$gse-E<gt>accession>

Accession number for the platform

=item C<$gse-E<gt>list("GSM" | "GPL")>

Since a series can contain more than one samples and platforms. This method can
get GSM list or GPL list that belong to the GSE record.

=item C<$gse-E<gt>merge>

merge single GSMs into a expression value matrix. The merging process is by platforms.
Each matrix is a GDS class object.

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

