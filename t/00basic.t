use strict;
use warnings;

use Test::More tests => 12;
BEGIN {
	use_ok('List::Vectorize') or BAIL_OUT "Unable to load List::Vectorize";
	use_ok('Microarray::ExprSet') or BAIL_OUT "Unable to load Microarray::ExprSet";
	use_ok('File::Basename') or BAIL_OUT "Unable to load File::Basename";
	use_ok('LWP::UserAgent') or BAIL_OUT "Unable to load LWP::UserAgent";
	use_ok('Thread') or BAIL_OUT "Unable to load Thread";
	use_ok('threads::shared') or BAIL_OUT "Unable to load threads::shared";
	use_ok('Time::HiRes') or BAIL_OUT "Unable to load Time::HiRes";
	use_ok('Microarray::GEO::SOFT') or BAIL_OUT "Unable to load Microarray::GEO::SOFT";
	use_ok('Microarray::GEO::SOFT::GSE') or BAIL_OUT "Unable to load Microarray::GEO::SOFT::GSE";
	use_ok('Microarray::GEO::SOFT::GSM') or BAIL_OUT "Unable to load Microarray::GEO::SOFT::GSM";
	use_ok('Microarray::GEO::SOFT::GPL') or BAIL_OUT "Unable to load Microarray::GEO::SOFT::GPL";
	use_ok('Microarray::GEO::SOFT::GDS') or BAIL_OUT "Unable to load Microarray::GEO::SOFT::GDS";
}
