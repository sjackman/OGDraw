# -*- perl -*-

# t/002_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'GeneMap::Plastome' , 'load class'); }

BEGIN {
	my $object = GeneMap::Plastome->new (file => 't/test.gb');
	isa_ok ($object, 'GeneMap::Plastome', 'instantiate from file');
}

