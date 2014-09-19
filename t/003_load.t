# -*- perl -*-

# t/003_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'GeneMap::Plastome::BoundingBox' ); }

my $object = GeneMap::Plastome::BoundingBox->new ();
isa_ok ($object, 'GeneMap::Plastome::BoundingBox');


