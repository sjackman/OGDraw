# -*- perl -*-

# t/007_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'GeneMap::Linear' ); }

my $object = GeneMap::Linear->new (file => 't/test.gb');
isa_ok ($object, 'GeneMap::Linear', 'instantiate from file');

