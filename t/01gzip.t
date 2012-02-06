use strict;
use warnings;

use Test::More tests => 1;

my $version = `gzip --version`;
my $has_gzip = ($version =~/gzip/im) + 0;

is($has_gzip, 1);
