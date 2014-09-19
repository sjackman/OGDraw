package GeneMap;

# BIG LICENCE TEXT HERE ?

=head1 NAME

GeneMap - A set of Perl modules designed to create high quality 
custom graphical maps of small genomes.

=head1 VERSION

Version 0.01

=cut


=head1 SYNOPSIS

The GeneMap package reads annotation information from the common
GeneBank sequence file format and generates high-resolution graphical 
maps in a range of image formats (TIFF, JPEG, PNG, GIF and PS).
It has been specially designed an optimized for the display of
small organelle genomes but it may as well be used on any annotated
DNA sequence that is available as a GeneBank file or accession
number. The package contains different modules providing the API
for circular (GeneMap::Plastome, GeneMap::Chondriome and GeneMap::Default)
or linear (GeneMap::Linear) maps. Please see the code examples below
as a starter. The detailed documentation of the specific classes 
has been included in POD format. So

	perldoc GeneMap::Plastome ...
	
and

	perldoc GeneMap::Linear
	
will give you all the dirty details.

=head2 Examples

 use GeneMap::Plastome;
 
 # Create a new plastome object using the GenBank entry of the 
 # tobacco chloroplast genome; include a %GC graph in the map
 my $map = GeneMap::Plastome->new( accnum => "Z00044",
                                   gc_cont => 1);
                                   
 # ... or from a local GenBank file ...
 my $nobelPrize = GeneMap::Plastome->new( file => "/home/userX/myBreakthrough.gb",
                                          gc_cont => 0);

 # Automatically detect the inverted repeat borders and include
 # them in the map - this doesn't work if the irscan tool is 
 # not properly installed
 $map->findIRborders();
 
 # ...or specify them manually (values are not correct...)
 $map->setIRA(1, 100000); 
 $map->setIRB(100001, 150000);
 
 # Read in a custom configuration file. Please read the FAQ
 # section on our website L<http://ogdraw.mpimp-golm.mpg.de/ssi_faqs.shtml>
 # to learn more about configuration files
 $map->readDrawableFeatures(file => 'myConfig.xml');
 
 # Write the map to a TIFF image
 # other possible formats are 'jpg', 'png', 'gif' and 'ps'
 # the map will be circular
 $map->createMap(outputfile => "plastome.tif", 
                 type => 'tif',
                 density => '300x300');
                 
                 
 use GeneMap::Linear;
 
 # Create a linear map of the same chloroplast genome
 my $linMap = GeneMap::Linear->new( accnum => "Z00044",
                                   gc_cont => 1);
 
 # Zoom into a smaller region 
 $linMap->setZoomRange(10000, 40000);
 
 # write the map to a PostScript(R) file which can subsequently
 # be edited manually using editors like Inkscape
 $linMap->createMap(outputfile => "linear_region.ps", 
                    type => 'ps');  


=cut


=head1 AUTHOR

 Marc Lohse and Oliver Drechsel, C<< <lohse<at>mpimp-golm.mpg.de, drechsel<at>mpimp-golm.mpg.de> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-genemap at rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=GeneMap>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc GeneMap

You can also look for information at:

=over 4

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/GeneMap>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/GeneMap>

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=GeneMap>

=item * Search CPAN

L<http://search.cpan.org/dist/GeneMap>

=back

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2008 Marc Lohse and Oliver Drechsel, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of GeneMap
