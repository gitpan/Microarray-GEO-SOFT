package Microarray::GEO::SOFT;

# 解析SOFT文件

use List::Vectorize;
require Microarray::GEO::SOFT::GPL;
require Microarray::GEO::SOFT::GSM;
require Microarray::GEO::SOFT::GDS;
require Microarray::GEO::SOFT::GSE;
use Microarray::ExprSet;


use File::Basename;
use LWP::UserAgent;
use Thread;
use threads::shared;
use Time::HiRes qw(usleep);

use strict;

our @ISA = qw();

our $VERSION = "0.10";

$| = 1;

# download geo files
our $ua;
our $response;


# 创建临时文件夹
our $_TMP_SOFT_DIR = ".tmp_soft";
opendir DIR, $_TMP_SOFT_DIR and closedir DIR
	or mkdir $_TMP_SOFT_DIR;

1;

#END {
#	rmtree($_TMP_SOFT_DIR);
#}
	
sub new {

	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = { "file" => "",
	             @_ };
	bless($self, $class);
	return $self;
	
}

sub soft_dir {
	return $_TMP_SOFT_DIR;
}

# 返回相应的子类
sub parse {

	my $self = shift;
	
	my $type = check_type($self->{file});
	
	my $obj;
	if($type eq "SERIES") {
	
		$obj = Microarray::GEO::SOFT::GSE->new(file => $self->{file});
		$obj->parse;
		
	}
	elsif($type eq "DATASET") {
		
		$obj = Microarray::GEO::SOFT::GDS->new(file => $self->{file});
		$obj->parse;
		
	}
	elsif($type eq "PLATFORM") {
		
		$obj = Microarray::GEO::SOFT::GPL->new(file => $self->{file});
		$obj->parse;
	}
	else {
	
		die "Format not supported.\n";
		
	}
	
	return $obj;
}

# 读取soft文件的前几行来判断是属于哪一种id
# 从而返回相应的类
sub check_type {

	my $file = shift;
		
	open F, $file or die "Cannot open $file.\n";
	
	while(my $line = <F>) {
	
		if($line =~/^\^SERIES /) {
			return "SERIES";
		}
		
		elsif($line =~/^\^DATASET /) {
			return "DATASET";
		}
		
		elsif($line =~/^\^PLATFORM /) {
			return "PLATFORM";
		}
	}
	
	return undef;
}

# 获得meta数据
sub meta {
	
	my $self = shift;
	
	return $self->{meta};

}

# 设置meta
sub set_meta {

	my $self = shift;
	
	my $accession = shift;
	my $title = shift;
	my $platform = shift;
	my $field = shift;
	
	$self->{meta}->{accession} = $accession;
	$self->{meta}->{title} = $title;
	$self->{meta}->{platform} = $platform;
	$self->{meta}->{field} = $field;

	return 1;
}

sub class {

	my $self = shift;
	
	return $self->{class};
}

sub set_class {

	my $self = shift;
	
	my $class = shift;
	
	$self->{class} = $class;
	
	return 1;
	
}

# GSE的list
# 包括GPL_list和GSM_list
sub list {

	my $self = shift;
	my $type = shift;
	
	if($self->class() ne "GSE") {
		die "only GSE is allowed.\n";
	}
	
	if($type ne "GPL" and $type ne "GSM") {
		die "wrong type in GSE list.\n";
	}
	
	return $self->{$type."_list"};

}

sub set_list {

	my $self = shift;
	
	my $list = shift;
	my $type = shift;
	
	if($self->class() ne "GSE") {
		die "only GSE is allowed.\n";
	}
	
	if($type ne "GPL" and $type ne "GSM") {
		die "wrong type in GSE list (GPL|GSM).\n";
	}
	
	$self->{$type."_list"} = $list;
	
	return 1;
}

sub rownames {
	
	my $self = shift;
	
	my $class = $self->class();
	
	if($class eq "GPL" or $class eq "GSM" or $class eq "GDS") {
		return $self->{table}->{table_rownames};
	}
	else {
		die "only GPL|GSM is allowed.\n";
	}
}

sub set_rownames {

	my $self = shift;
	my $new_rownames = shift;
	
	if(len($self->{table}->{table_rownames}) != len($new_rownames)) {
		die "length of new rownames should be equal to the old rownames.\n";
	}
	
	$self->{table}->{table_rownames} = $new_rownames;
	
	return 1;
}

sub colnames {
	
	my $self = shift;
	
	my $class = $self->class();
	
	if($class eq "GPL" or $class eq "GSM" or $class eq "GDS") {
		return $self->{table}->{table_colnames};
	}
	else {
		die "only GPL|GSM is allowed.\n";
	}
}

sub set_colnames {

	my $self = shift;
	my $new_colnames = shift;
	
	if(len($self->{table}->{table_colnames}) != len($new_colnames)) {
		die "length of new colnames should be equal to the old colnames.\n";
	}
	
	$self->{table}->{table_colnames} = $new_colnames;
	
	return 1;
}


sub matrix {
	
	my $self = shift;
	
	my $class = $self->class;
	if($class eq "GPL") {
		return $self->{table}->{table_matrix};
	}
	elsif($class eq "GSM" or $class eq "GDS") {
		my $accession = $self->accession;
		open F, "$_TMP_SOFT_DIR/$accession.tab" or die "cannot open $_TMP_SOFT_DIR/$accession.tab\n";
		
		my $matrix;
		while(my $line = <F>) {
			chomp $line;
			my @tmp = split "\t", $line;
			push(@$matrix, [@tmp]);
		}
		
		return $matrix;
	}
	else {
		die "wrong parameter\n";
	}
}

sub accession {

	my $self = shift;
	
	return $self->{meta}->{accession};

}

sub title {

	my $self = shift;
	
	return $self->{meta}->{title};

}

sub field {

	my $self = shift;
	
	return $self->{meta}->{field};

}

sub platform {

	my $self = shift;
	
	return $self->{meta}->{platform};

}

# 通过GPL类来修改rownames
# 允许修改的convert的类为GSM和GDS
sub id_convert {

	my $self = shift;
	my $gpl = shift;
	my $to_id = shift;
	
	my $class = $self->class;
	if($class ne "GSM" and $class ne "GDS") {
		die "id convert only permit GSM and GDS\n";
	}
	
	my $platform_id = $gpl->accession;
	my $available_field = $gpl->field;
	
	if(! $available_field->{$to_id}) {
		die "wrong ID ($to_id) in $platform_id\n";
	}
	
	my $new_rownames = $gpl->mapping($to_id);
	
	$self->set_rownames($new_rownames);
	
	return $self;

}

# 将soft类转变为exprset类
# 可以转的类为GSM和GDS
sub soft2exprset {

	my $self = shift;
	
	my $class = $self->class;
	if($class ne "GSM" and $class ne "GDS") {
		die "convert only permit GSM and GDS\n";
	}
	
	my $eset = Microarray::ExprSet->new;
	$eset->set_feature($self->rownames);
	
	# phenotype 应该是每个GSM的title
	if($class eq "GSM") {
		$eset->set_phenotype([$self->title]);
	}
	elsif($class eq "GDS") {
		$eset->set_phenotype(sapply($self->colnames, sub{$self->field->{$_[0]}}));
	}
	
	$eset->set_matrix($self->matrix);
	
	return $eset;
	
}







# 完整的从geo上下载数据
# 返回下载到本地的文件列表(数组索引)
sub download {

	my $self = shift;
	
	my $id = shift;
	
	my %option = ( "proxy" => "",        # proxy setting, only http, should be like "http://username:password@127.0.0.1:808/"
				   "timeout" => 30,
	               @_ );
	
	my $remote_file_list;
	my $remote_file_name;
	my $remote_file_size;
	my $local_file_name;
	
	$ua = LWP::UserAgent->new;
	$ua->timeout($option{timeout});
	
	if($option{proxy}) {
		$ua->proxy(["http"], $option{proxy});
	}
	
	my $url;
	
	# 不同类型数据的url格式
	my $url_format = { "gse" => "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series",
					   "gpl" => "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_platform",
					   "gds" => "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/GDS" };
	
	my $chip_file_format = lc($option{"format"});
	
	# 根据不同的芯片格式选择url
	if($id =~/^gse\d+$/i) {
		$url = "$url_format->{gse}/$id";
	}
	elsif($id =~/^gpl\d+$/i) {
		$url = "$url_format->{gpl}/$id";
	}
	elsif($id =~/^gds\d+$/i) {
		$url = "$url_format->{gds}/$id";
	}
	
	# 如果是GSE或者GPL
	if($id =~/^gse\d+$/i or $id =~/gpl\d+$/i) {
		
		# 先获得目录下的文件，因为某些GSE或者GPL文件下会有多个文件
		print "Reading dir from GEO FTP site:\n";
		print "  $url\n\n";
		$response = $ua->get($url);
		
		unless($response->is_success) {
			die $response->status_line;
		}
		
		my $content = $response->content;
		@$remote_file_list = split "\n", $content;
		
		print "found ", scalar(@$remote_file_list), " file.\n";
		
		for(my $i = 0; $i < scalar(@$remote_file_list); $i ++) {
			my @tmp = split " ", $remote_file_list->[$i];
			push(@$remote_file_name, $tmp[$#tmp]);
			push(@$remote_file_size, $tmp[4]);
			push(@$local_file_name, "$_TMP_SOFT_DIR/$tmp[$#tmp]");
		}
		
		foreach (@$remote_file_name) {
			print "  $_\n";
		}
		print "\n";
	} # 如果是GDS数据
	elsif($id =~/gds\d+$/i) {
		print "validating link from GEO FTP site:\n";
		print "  $url/$id.soft.gz\n\n";
		$response = $ua->head("$url/$id.soft.gz");
		
		unless($response->is_success) {
			die $response->status_line;
		}
		
		print "found $id.soft.gz on the server.\n\n";
		$remote_file_name->[0] = "$id.soft.gz";
		$remote_file_size->[0] = $response->header("content-length");
		$local_file_name->[0] = "$_TMP_SOFT_DIR/$id.soft.gz";
		
	}
	
	for(my $i = 0; $i < scalar(@$remote_file_name); $i ++) {
		my $local_file = "$_TMP_SOFT_DIR/$remote_file_name->[$i]";
		
		# 检查当地文件夹中是否有同名文件
		if(-e $local_file) {
			my $r = int(rand(10000000000));
			if($remote_file_name->[$i] =~/^(.*?)\.(\S+)$/) {
				my $base = $1;
				my $ext = $2;
				$local_file = "$_TMP_SOFT_DIR/$base.$r.$ext";
			}
			else {
				$local_file = "$_TMP_SOFT_DIR/$remote_file_name->[$i].$r";
			}
			
			$local_file_name->[$i] = $local_file;
		}
		
		$url = "$url/$remote_file_name->[$i]";
		print "downloading $url\n";
		print "file size: $remote_file_size->[$i] byte.\n";
		print "local file: $local_file\n\n";
		
		# 开始下载，并显示下载进度
		$response = undef;
		my $response : shared;
		my $f1 = Thread->new(\&_download, $url, $local_file);
		my $f2 = Thread->new(\&_progress, $local_file, $remote_file_size->[$i]);
		
		$f1->join;
		$f2->join;
		
	}
	
	my $local_uncompressed_file_name;
	foreach my $f (@$local_file_name) {
		my $fn = _decompress($f);
		push(@$local_uncompressed_file_name, $fn);
	}
	
	# 原先设计还可以洗在series matrix (目录下可能会有超过一个的文件)
	$self->{file} = $local_uncompressed_file_name->[0];
	return 1;
}

# 下载
sub _download {
	my $url = shift;
	my $local_file = shift;
		
	$response = $ua->get($url, ":content_file" => $local_file);	

}

# 显示下载进度
sub _progress {
	my $local_file = shift;
	my $remote_file_size = shift;
	
	my $s_sleep = 1000000;
	my $s_sleep_ms = int($s_sleep / 1000);
	
	# 一开始文件可能正在连接中，还未开始下载
	while(! -e $local_file) {
		#print "$local_file does not exist, sleep $s_sleep_ms ms.\n";
		usleep($s_sleep);
		
		if($response) {
			print "\n\n";
			last;
		}
	}
									 
	my $recieved_file_size = -s "$local_file";
	my $i = 0;
	my $bar = ["|", "\\", "-", "/"];
	while($recieved_file_size != $remote_file_size) {
		$recieved_file_size = -s "$local_file";
		my $percentage = $recieved_file_size / $remote_file_size;
		
		$percentage = sprintf("%.1f", $percentage * 100);
		print "\b" x 100;
		
		
		print "[", $bar->[$i % scalar(@$bar)], "]";
		print " Recieving $recieved_file_size byte. $percentage\%";
		$i ++;
		
		usleep(1000);
		
		# 如果下载完毕
		if($response) {
			last;
		}
	}
	
	print "\n\n";
}


# 解压缩
sub _decompress {
	
	# 压缩文件
	my $compressed_file = shift;
	
	my $version = `gzip --version`;
	unless($version =~/gzip/im) {
		die "You need to install gzip and put the path of gzip into you PATH envirionment variable.\n";
	}
	
	# 压缩文件的文件名
	my $basename = basename($compressed_file);
	
	print "decompress $compressed_file.\n";
	my $command;
	
	# 获得解压缩文件的文件名
	$command = "gzip -l \"$compressed_file\"";
	my $status = `$command`;
	my @foo = split "\n", $status;
	@foo = split " ", $foo[1];
	my $uncompressed_file = $foo[$#foo];
	
	# 解压缩
	$command = "gzip -cd \"$compressed_file\" > \"$uncompressed_file\"";
	
	$status = `$command`;
	if($status) {
		die "$status\n";
	}
	
	# 返回解压后的文件名
	return "$uncompressed_file";
}



__END__

=pod

=head1 NAME

Microarray::GEO::SOFT - Reading microarray data in SOFT format from GEO database.

=head1 SYNOPSIS

  use Microarray::GEO::SOFT;
  use strict;
  
  # initialize
  my $soft = Microarray::GEO::SOFT->new; 
  
  # download
  $soft->download("GSE19513");
  $soft->download("GPL6793");
  $soft->download("GDS3718");
  
  # or else you can read local data
  $soft = Microarray::GEO::SOFT->new(file => "GSE19513.soft");
  
  # parse
  # $data would be a object of Microarray:GEO::SOFT::GSE, Microarray::GEO::SOFT::GPL
  # or Microarray::GEO::SOFT::GDS class
  my $data = $soft->parse;
  
  # meta info
  $data->meta;
  $data->title;
  $data->platform;
  $data->field;
  
  # GPL belongs to GSE
  my $gpl = $data->list("GPL")->[0];
  
  # merge GSMs belonging to a same GPL into a whole
  my $g = $data->merge->[0];
  
  # transform the uid from probe id to gene symbol
  $g->id_convert($gpl, "Gene Symbol");
  
  # transform into Microarray::ExprSet class object
  my $e = $g->soft2exprset;
  
  # eliminate the blank lines
  $e->remove_empty_feature;
  
  # make all symbols unique
  $e->unique_feature;
  
  # obtain the expression matrix
  $e->matrix;	

=head1 DESCRIPTION

GEO (Gene Expression Omnibus) is the biggest database providing gene expression
profile data. This module provides method to download and parse files in GEO database
and transform them into format for common usage.

There are always four type of data in GEO which are GSE, GPL, GSM and GDS.

GPL: Platform of the microarray, like Affymetrix U133A

GSM: A single microarray

GSE: A complete microarray experiment, always contains multi GSMs and multi GPLS

GDS: manually collected data sets from GSE, only 1 platform

Data stored in GEO database has several formats. We provide method to parse the most
used format: SOFT formatted family files. The origin data is downloaded from GEO ftp site.

=head2 Subroutines

=over 4

=item C<new("file" = $file)>

Initial a Microarray::GEO::SOFT class object. The only argument is file path for 
the microarray data in SOFT format or a file handle that has been openned.

=item C<$soft-E<gt>download(ACC, %options>

Download GEO record from NCBI website. The first argument is the accession number
such as (GSExxx, GPLxxx or GDSxxx). Your can set the timeout and proxy via C<%options>.
the proxy should be set as http://username:password@server-addr:port.

=item C<$soft-E<gt>parse>

Proper parsing method is selected according to the accession number of GEO record.
E.g. if a GSExxx record is required, then the parsing function would choose method
to parse GSExxx part and return a C<Microarray::GEO::SOFT::GSE> class object.

=item C<$data-E<gt>meta>

Get meta information, more detailed meta information can be get via C<platform>, 
C<title>, C<field>, C<accession>.

=item C<$data-E<gt>platform>

Get accession number of the platform. If a record has multiple platforms, the function
return a reference of array.

=item C<$data-E<gt>title>

Title of the record

=item C<$data-E<gt>field>

Description of each field in the data matrix

=item C<$data-E<gt>accession>

Accession number for the record

=item C<$gds-E<gt>id_convert($gpl, id)>

Change the primary id for genes which always the rownames by default. Mapping information
is provided in GPL record. The first argument is the GPL record corresponding to 
the GDS record, the id argument is from colnames in the GPL record. Use 
C<$gpl->>field> or C<$gpl->>colnames> to find the ID names to convert.

=item C<$gds-E<gt>soft2exprset>

Transform C<Microarray::GEO::SOFT> class object to C<Microarray::ExprSet> class object.

=back

=head1 AUTHOR

Zuguang Gu E<lt>jokergoo@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2012 by Zuguang Gu

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.1 or,
at your option, any later version of Perl 5 you may have available.

=head1 SEE ALSO

L<Microarray::ExprSets>

=cut
