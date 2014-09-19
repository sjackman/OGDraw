######################################
# This Module defines a simple class #
# to store the bounding box border   #
# coordinates of a text label. It    #
# implements a simple 2D collision   #
# detection mechanism                #
# part of the GeneMap package        #
######################################
# author : Marc Lohse                #
######################################


package GeneMap::Plastome::BoundingBox;



use strict;
use warnings;
use Carp;

# for standard helper subs 

use GeneMap::Plastome::ToolBox;

use constant TRUE   => 1;
use constant FALSE  => undef();


sub new
{
    my ($self, %args) = @_;
    _check_args(\%args, 'x', 'y', 'height', 'width');
    #
    # here some checking of the values - do they describe a rectangle?
    #
    return bless    {
                        x    => $args{x},
                        y    => $args{y},
                        height  => $args{height},
                        width   => $args{width}
                    }, $self;
} # end newCOMP_LABEL_POS	

sub collides_with
{
    # detect 2D spatial overlap between the instance on which
    # the method was called and a second instance of BoundingBox
    my ($self, $other_box) = @_;
    my 	(
    	$left1, $left2,
    	$right1, $right2,
    	$top1, $top2,
    	$bottom1, $bottom2
    	);
    
    $left1 = $self->getX();
    $left2 = $other_box->getX();
    $right1 = $self->getX() + $self->width();
    $right2 = $other_box->getX() + $other_box->width();
    $top1 = $self->getY();
    $top2 = $other_box->getY();
    $bottom1 = $self->getY() + $self->height();
    $bottom2 = $other_box->getY() + $other_box->height();
    
    return FALSE if ($bottom1 < $top2) ;
    return FALSE if ($top1 > $bottom2) ;
  
    return FALSE if ($right1 < $left2) ;
    return FALSE if ($left1 > $right2) ;

    return TRUE;
}

#####################
# accessory methods #
#####################


sub getX
{
    my ($self) = shift;
    return $self->{x};
}

sub setX
{
    my ($self, $new_x) = @_;
    $self->{x} = $new_x;
}

sub getY
{
    my ($self) = shift;
    return $self->{y};
}

sub setY
{
    my ($self, $new_x) = @_;
    $self->{y} = $new_x;
}

sub width
{
    my ($self) = shift;
    return $self->{width};
}

sub height
{
    my ($self) = shift;
    return $self->{height};
}

sub setWidth
{
	my ($self) = shift;
	my $new_width = shift;
	$self->{width} = $new_width;	
}

sub setHeight
{
	my ($self) = shift;
	my $new_height = shift;
	$self->{height} = $new_height;	
}

sub swapDimensions
{
	my ($self) = shift;
	my $width = $self->width();
	$self->setWidth($self->height);
	$self->setHeight($width);	
}


1;
