##########################################################
#                   Default.pm                           #
##########################################################
# Default set  of drawable features                      #
#                                                        #
#                                                        #
# Authors:                                               #
# Marc Lohse and Oliver Drechsel                         #
# Max-Planck-Institute of molecular                      #
# plant physiology                                       #
# Am Muehlenberg 1                                       #
# 14476 Potsdam-Golm                                     #
# lohse@mpimp-golm.mpg.de                                #
# drechsel@mpimp-golm.mpg.de                             #
#                                                        #
##########################################################

package GeneMap::Default;

#
# include modules
#

use strict;

#
# constants
#

use constant BLACK  => 0,   0,   0;
use constant WHITE  => 255, 255, 255;
use constant RED    => 255, 0,   0;
use constant GREEN  => 0,   255, 0;
use constant BLUE   => 0,   0,   255;
use constant ATP    => 151, 190, 13;
use constant PSA    => 0,   102, 44;
use constant PSB    => 50,  137, 37;
use constant RBCL   => 31,  161, 45;
use constant PET    => 121, 156, 19;
use constant TRN    => 22,  41,  131;
use constant ORF    => 87,  185, 168;
use constant NDH    => 255, 236, 0;
use constant CLP    => 233, 93,  15;
use constant RRN    => 226, 0,   26;
use constant RPO    => 189, 18,  32;
use constant RPS    => 219, 170, 115;
use constant RPL    => 158, 119, 66;
use constant YCF    => 255, 250, 208;
use constant ORI    => 255, 128, 128;
use constant VIOLET => 171, 37,  157;
use constant DARK_RED => 127,  0,  0;

use constant TRUE 		        => 1;    # quite "unperlish" but improves the readability (at least for me)
use constant FALSE              => undef();# was 0

# common template for the drawable features list for mitochondria
our @drawable_features = 
(
{type => 'gene', 		pattern => '^rpo.*',			color => [RPO],		fullname => 'RNA polymerase', 			drawflag => TRUE},
{type => 'gene|CDS', 	pattern => '^rps.*',			color => [RPS],		fullname => 'ribosomal proteins (SSU)', drawflag => TRUE},
{type => 'gene', 		pattern => '^rpl.*',			color => [RPL],		fullname => 'ribosomal proteins (LSU)', drawflag => TRUE},
{type => 'tRNA', 		pattern => '.*', 				color => [TRN],		fullname => 'transfer RNAs', 			drawflag => TRUE},
{type => 'rRNA', 		pattern => '.*', 				color => [RRN],		fullname => 'ribosomal RNAs', 			drawflag => TRUE},
{type => 'rep_origin', 	pattern => '^ori.*', 			color => [ORI],		fullname => 'origin of replication', 	drawflag => TRUE},
{type => 'intron', 		pattern => '.*',				color => [WHITE],	fullname => 'introns', 					drawflag => TRUE},
{type => 'other', 		pattern => '.*',				color => [VIOLET],	fullname => 'other genes', 				drawflag => TRUE},
{type => 'CDS|gene', 	pattern => '^orf.*',			color => [ORF],		fullname => 'ORFs', 					drawflag => TRUE},
{type => 'gene', 		pattern => '(^clp.*)|(^mat.*)',	color => [CLP],		fullname => 'clpP, matK', 				drawflag => TRUE},
{type => '_operon_',	pattern => '.*',				color => [RED],		fullname => 'polycistronic transcripts',drawflag => TRUE}
);
1;
