

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'GeneMap',
    VERSION_FROM => 'lib/GeneMap/Plastome.pm', 
    AUTHOR       => 'Marc Lohse and Oliver Drechsel (lohse@mpimp-golm.mpg.de, drechsel@mpimp-golm.mpg.de)',
    ABSTRACT     => 'Create high quality graphical genome maps',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                     'Bio::Root::Root' => 0.0,
                     'Bio::SeqIO' => 0.0,
                     'Bio::SeqFeature::Generic' => 0.0,
                     'Bio::Restriction::Analysis' => 0.0,
                     'Bio::Perl' => 0.0,
                     'PostScript::Simple' => 0.0,
                     'Image::Magick' => 0.0,
                     'XML::Generator' => 0.0,
                     'XML::Simple' => 0.0                     
                    },
	EXE_FILES	=> ['bin/drawgenemap'],
);
