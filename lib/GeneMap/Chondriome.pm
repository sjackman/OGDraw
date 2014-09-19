##########################################################
#                   Chondriome.pm                        #
##########################################################
# Just the modified list of drawable features for        #
# mitochondria                                           #
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

package GeneMap::Chondriome;

#
# include modules
#

use strict;

use vars qw/@ISA/;

our @ISA = qw/GeneMap::Plastome/;

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
use constant SDH    =>  52, 211,  77;
use constant COB    => 200, 250,  40;
use constant COX    => 255, 180, 255;

use constant VIOLET => 171, 37,  157;
use constant DARK_RED => 127,  0,  0;

use constant TRUE 		        => 1;    # quite "unperlish" but improves the readability (at least for me)
use constant FALSE              => undef();# was 0

# common template for the drawable features list for mitochondria
our @drawable_features = 
(

{type => 'gene', 		pattern => '(^nad.*)|(^nd.*)', 	color => [NDH],		fullname => 'complex I (NADH dehydrogenase)', 		drawflag => TRUE},
{type => 'gene', 		pattern => '^sdh.*',			color => [SDH],		fullname => 'complex II (succinate dehydrogenase)', 	drawflag => TRUE},
{type => 'gene', 		pattern => '^cob.*',			color => [COB],		fullname => 'complex III (ubichinol cytochrome c reductase)', 	drawflag => TRUE},
{type => 'gene', 		pattern => '^cox.*',			color => [COX],		fullname => 'complex IV (cytochrome c oxidase)', 	drawflag => TRUE},
{type => 'gene', 		pattern => '^atp.*', 			color => [ATP],		fullname => 'ATP synthase', 			drawflag => TRUE},
{type => 'gene', 		pattern => '^ccb.*',			color => [PSB],		fullname => 'cytochrome c biogenesis', 	drawflag => TRUE},
{type => 'gene', 		pattern => '^rpo.*',			color => [RPO],		fullname => 'RNA polymerase', 			drawflag => TRUE},
{type => 'gene', 	pattern => '^rps.*',			color => [RPS],		fullname => 'ribosomal proteins (SSU)', drawflag => TRUE},
{type => 'gene', 		pattern => '^rpl.*',			color => [RPL],		fullname => 'ribosomal proteins (LSU)', drawflag => TRUE},
{type => 'gene', 		pattern => '(^clp.*)|(^mat.*)',	color => [CLP],		fullname => 'maturases', 				drawflag => TRUE},
{type => 'other', 		pattern => '.*',				color => [VIOLET],	fullname => 'other genes', 				drawflag => TRUE},
{type => 'gene', 	pattern => '^orf.*',			color => [ORF],		fullname => 'ORFs', 					drawflag => TRUE},
{type => 'tRNA', 		pattern => '.*', 				color => [TRN],		fullname => 'transfer RNAs', 			drawflag => TRUE},
{type => 'rRNA', 		pattern => '.*', 				color => [RRN],		fullname => 'ribosomal RNAs', 			drawflag => TRUE},
{type => 'rep_origin', 	pattern => '^ori.*', 			color => [ORI],		fullname => 'origin of replication', 	drawflag => TRUE},
{type => 'intron', 		pattern => '.*',				color => [WHITE],	fullname => 'introns', 					drawflag => TRUE},
{type => '_operon_',	pattern => '.*',				color => [RED],		fullname => 'polycistronic transcripts',drawflag => FALSE}
);



# sub _preprocess_feature_list
# {
#     my $self = shift;
#     my @seen_angles;
#     
#     # here we have to check whether the feature is a gene that is not covered by the 
#     # drawable features list and then define type as other. problem is that all tRNAs
#     # and rRNAs seem to meet these conditions because they are also annotated as genes
#     # and we only filter for type 'tRNA' resp. type 'rRNA' entries. We need a means to
#     # see whether a feature is already on the map....
#     
#     # save the label angles of all drawable features that will appear on the map
#     for (my $c = 0; $c <= $#{$self->{_featureList}}; $c++)
#     {
#         push @seen_angles, ${$self->{_featureList}}[$c]->{labelAngle} 
#             if (_is_drawable_feature($self, ${$self->{_featureList}}[$c], FALSE));
#     }
#     
#     for (my $c = 0; $c <= $#{$self->{_featureList}}; $c++)
#     {
#         # is the feature a gene that is not drawable and not
#         # at the same angle position as one of the drawable features?
#         # if yes tag it 'other' << problem the ones that are already tagged other
#         # will be drawable features ...... but it seems to work, anyway.
#         if (!_is_default_drawable_feature($self, ${$self->{_featureList}}[$c], FALSE) 
#                 && 
#             (${$self->{_featureList}}[$c]->{type} eq 'gene')
#             	&& 
#             (${$self->{_featureList}}[$c]->{name} !~ /^rn.*/i)
#                 &&
#              (!_is_element_of(${$self->{_featureList}}[$c]->{labelAngle}, @seen_angles )) )
#             
#         {
#              print "wir sind anders:".${$self->{_featureList}}[$c]->{name}."<br>";
#              ${$self->{_featureList}}[$c]->{type} = 'other';
#         }
#     }
#     
# } # end _preprocess_feature_list

sub _is_drawable_feature
{
	my ($self, $test, $p_color) = @_;
	for (my $d = 0; $d <= $#{$self->{'_drawableFeatures'}}; $d++)
			{
				if (($test->{type} =~ /${$self->{_drawableFeatures}}[$d]->{type}/i) and 
				($test->{name} =~ /${$self->{_drawableFeatures}}[$d]->{pattern}/i) and 
				(${$self->{_drawableFeatures}}[$d]->{drawflag})) 
				{
				${$p_color} = ${$self->{_drawableFeatures}}[$d]->{color} unless (!$p_color);
				return TRUE;
				}		 
			}
	 return FALSE;
}

sub _is_default_drawable_feature
{
	my ($self, $test, $p_color) = @_;
	
	
	
	for (my $d = 0; $d <= $#{$self->{'_drawableFeatures'}};$d++)
			{
				if (($test->{type} =~ /$GeneMap::Chondriome::drawable_features[$d]->{type}/i) and 
				($test->{name} =~ /$GeneMap::Chondriome::drawable_features[$d]->{pattern}/i) 
# 				and (@drawable_features[$d]->{drawflag})
				) 
				{
# 				print"CHONDRIOME";
				${$p_color} = $GeneMap::Chondriome::drawable_features[$d]->{color} unless (!$p_color);
				return TRUE;
				}		 
			}
	 return FALSE;
}

1;