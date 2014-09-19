############################
#  scribbled together by Marc Lohse
#  general purpose subroutine library
############################

package GeneMap::Plastome::ToolBox;

use strict;
use warnings;
use Bio::Root::Root;
use Bio::Seq;
use Carp;

use constant TRUE 		=>  1;
use constant FALSE      => undef();

#
# Inline C routines
#



# use Inline C => <<'END_OF_C_CODE';
# 
# #include <stdio.h>
# #include <string.h>
# 
# void greet(SV* name1, ...) 
# {
# 	Inline_Stack_Vars;
# 	int i;
# 	for (i = 0; i < Inline_Stack_Items; i++) 
# 		printf("Hello %s!\n", SvPV(Inline_Stack_Item(i), PL_na));
# 	Inline_Stack_Void;
# }
# 
# char* test_add (int a, int b)
# {
# 	char *s;
# 	s = "TESTEST\0";
# 	printf("in C sub: %s\n", s);
# 	return s;
# }
# 
# /*
# char* convert(const char *A) 
# {
# 	int i;
# 	int len = strlen(A);
# 	char *B;
# 	
# 	for (i= 0; i < len; i++)
# 	{                        
# 		if      (A[i] == 'A') B[len-1-i] = 'T';
# 		else if (A[i] == 'T') B[len-1-i] = 'A';
#   		else if (A[i] == 'G') B[len-1-i] = 'C';
#   		else if (A[i] == 'C') B[len-1-i] = 'G';	
#   	}
#   	return B;
# }
# */
# 
# END_OF_C_CODE
# #
#
#


require Exporter;
use vars       qw($VERSION @ISA @EXPORT);

$VERSION = 0.15;

@ISA         = qw(Exporter);

@EXPORT      = 	qw(	
					&format_seq 
					&is_element_of 
					&break_sequence
					&unique
					&complementary
					&rev_com
					&insert_element
					&calcGCcontent
					&count_digits
                    &min
                    &max
                    &_check_args
                    &print_in_columns
                    &getStandardEnzymes
                    &print_in_html_columns
                    &convert
                    &greet
                    &test_add
                    &trunc_with_features
                    &is_same_feature
                    &string_width_of
                    &US_format_number
                    &isVNTIfile
                  );
                  
my %standardEnzymes = 
( 	'AatII' => 'GACGTC 5', 'AccI' => 'GTMKAC 2', 'AclI' => 'AACGTT 2', 'AcyI' => 'GRCGYC 2',
	'AflII' => 'CTTAAG 1', 'AflIII' => 'ACRYGT 1', 'AgeI' => 'ACCGGT 1', 'AhaIII' => 'TTTAAA 3',
	'AhdI' => 'GACNNNNNGTC 6', 'AluI' => 'AGCT 2', 'AlwNI' => 'CAGNNNCTG 6', 'ApaBI' => 'GCANNNNNTGC 8',
	'ApaI' => 'GGGCCC 5', 'ApaLI' => 'GTGCAC 1', 'ApoI' => 'RAATTY 1', 'AscI' => 'GGCGCGCC 2', 
	'AsuI' => 'GGNCC 1', 'AsuII' => 'TTCGAA 2', 'AvaI' => 'CYCGRG 1', 'AvaII' => 'GGWCC 1', 
	'AvrII' => 'CCTAGG 1', 'BalI' => 'TGGCCA 3', 'BamHI' => 'GGATCC 1', 'BclI' => 'TGATCA 1', 
	'BetI' => 'WCCGGW 1', 'BglI' => 'GCCNNNNNGGC 7', 'BglII' => 'AGATCT 1', 'BsaAI' => 'YACGTR 3', 
	'BsaBI' => 'GATNNNNATC 5', 'BsePI' => 'GCGCGC 1', 'BsiYI' => 'CCNNNNNNNGG 7', 'Bsp1407I'=> 'TGTACA 1', 
	'BspHI' => 'TCATGA 1', 'BspLU11I'=> 'ACATGT 1', 'BspMII' => 'TCCGGA 1', 'BstEII' => 'GGTNACC 1', 
	'BstXI' => 'CCANNNNNNTGG 8', 'Cac8I' => 'GCNNGC 3', 'CauII' => 'CCSGG 2', 'Cfr10I' => 'RCCGGY 1', 
	'CfrI' => 'YGGCCR 1', 'ClaI' => 'ATCGAT 2', 'CviJI' => 'RGCY 2', 'CviRI' => 'TGCA 2', 'DdeI' => 'CTNAG 1', 
	'DpnI' => 'GATC 2',	 'DraI' => 'TTTAAA 3', 'DraII' => 'RGGNCCY 2', 'DraIII' => 'CACNNNGTG 6', 
	'DrdI' => 'GACNNNNNNGTC 7', 'DsaI' => 'CCRYGG 1', 'Eam1105I'=> 'GACNNNNNGTC 6', 'Eco47III'=> 'AGCGCT 3', 
	'EcoNI' => 'CCTNNNNNAGG 5', 'EcoRI' => 'GAATTC 1', 'EcoRII' => 'CCWGG 0', 'EcoRV' => 'GATATC 3', 
	'EspI' => 'GCTNAGC 2', 'Fnu4HI' => 'GCNGC 2', 'FnuDII' => 'CGCG 2', 'FseI' => 'GGCCGGCC 6', 
	'HaeI' => 'WGGCCW 3', 'HaeII' => 'RGCGCY 5', 'HaeIII' => 'GGCC 2', 'HgiAI' => 'GWGCWC 5', 
	'HgiCI' => 'GGYRCC 1', 'HgiJII' => 'GRGCYC 5', 'HhaI' => 'GCGC 3',	 'HincII' => 'GTYRAC 3', 
	'HindII' => 'GTYRAC 3', 'HindIII' => 'AAGCTT 1', 'HinfI' => 'GANTC 1', 'HpaI' => 'GTTAAC 3', 
	'HpaII' => 'CCGG 1', 'KpnI' => 'GGTACC 5', 'MaeI' => 'CTAG 1', 'MaeII' => 'ACGT 1', 'MaeIII' => 'GTNAC 0', 
	'MboI' => 'GATC 0', 'McrI' => 'CGRYCG 4', 'MfeI' => 'CAATTG 1', 'MluI' => 'ACGCGT 1', 'MseI' => 'TTAA 1', 
	'MslI' => 'CAYNNNNRTG 5', 'MstI' => 'TGCGCA 3', 'MwoI' => 'GCNNNNNNNGC 7', 'NaeI' => 'GCCGGC 3', 
	'NarI' => 'GGCGCC 2', 'NcoI' => 'CCATGG 1', 'NdeI' => 'CATATG 2', 'NheI' => 'GCTAGC 1', 'NlaIII' => 'CATG 4', 
	'NlaIV' => 'GGNNCC 3', 'NotI' => 'GCGGCCGC 2', 'NruI' => 'TCGCGA 3', 'NspBII' => 'CMGCKG 3', 'NspI' => 'RCATGY 5', 
	'PacI' => 'TTAATTAA 5', 'PflMI' => 'CCANNNNNTGG 7', 'PmaCI' => 'CACGTG 3', 'PmeI' => 'GTTTAAAC 4', 
	'PpuMI' => 'RGGWCCY 2', 'PshAI' => 'GACNNNNGTC 5', 'PstI' => 'CTGCAG 5', 'PvuI' => 'CGATCG 4', 
	'PvuII' => 'CAGCTG 3', 'RsaI' => 'GTAC 2', 'RsrII' => 'CGGWCCG 2', 'SacI' => 'GAGCTC 5', 'SacII' => 'CCGCGG 4',	 
	'Sau96I' => 'GGNCC 1', 'SalI' => 'GTCGAC 1', 'SanDI' => 'GGGWCCC 2', 'SauI' => 'CCTNAGG 2',	 'SbfI' => 'CCTGCAGG 6', 
	'ScaI' => 'AGTACT 3', 'ScrFI' => 'CCNGG 2', 'SduI' => 'GDGCHC 5', 'SecI' => 'CCNNGG 1', 'SexAI' => 'ACCWGGT 1', 
	'SfeI' => 'CTRYAG 1', 'SfiI' => 'GGCCNNNNNGGCC 8', 'SgfI' => 'GCGATCGC 5', 'SgrAI' => 'CRCCGGYG 2', 
	'SmaI' => 'CCCGGG 3', 'SmlI' => 'CTYRAG 1', 'SnaBI' => 'TACGTA 3', 'SpeI' => 'ACTAGT 1', 'SphI' => 'GCATGC 5', 
	'SplI' => 'CGTACG 1', 'SrfI' => 'GCCCGGGC 4', 'Sse8387I'=> 'CCTGCAGG 6', 'Sse8647I'=> 'AGGWCCT 2', 
	'SspI' => 'AATATT 3', 'StuI' => 'AGGCCT 3', 'StyI' => 'CCWWGG 1', 'SwaI' => 'ATTTAAAT 4', 'TaqI' => 'TCGA 1', 
	'TatI' => 'WGTACW 1', 'TfiI' => 'GAWTC 1', 'TseI' => 'GCWGC 1', 'Tsp45I' => 'GTSAC 0', 'Tsp4CI' => 'ACNGT 3', 
	'TspEI' => 'AATT 0', 'TspRI' => 'CASTGNN 7', 'Tth111I' => 'GACNNNGTC 4', 'VspI' => 'ATTAAT 2', 
	'XbaI' => 'TCTAGA 1', 'XcmI' => 'CCANNNNNNNNNTGG 8', 'XhoI' => 'CTCGAG 1', 'XhoII' => 'RGATCY 1', 
	'XmaIII' => 'CGGCCG 1', 'XmnI' => 'GAANNNNTTC 5', 'BssHII' => 'GCGCGC', 'SapI' => 'not_important' 
);

#TODO add more special characters
my %charwidths = 
(
	'a' => 6.627, 'b' => 6.627, 'c' => 6.025, 'd' => 6.627, 'e' => 6.627, 'f' => 3.313, 'g' => 6.627, 'h' => 6.627,
	'i' => 2.711, 'j' => 2.711, 'k' => 6.025, 'l' => 2.711, 'm' => 9.941, 'n' => 6.627, 'o' => 6.627, 'p' => 6.627,
	'q' => 6.627, 'r' => 3.916, 's' => 6.025, 't' => 3.313, 'u' => 6.627, 'v' => 6.025, 'w' => 8.736, 'x' => 6.025,
	'y' => 6.025, 'z' => 6.025, 'A' => 8.133,'B' => 8.133, 'C' => 8.736, 'D' => 8.736, 'E' => 8.133, 'F' => 7.23,
	'G' => 9.338, 'H' => 8.736, 'I' => 3.313, 'J' => 6.025, 'K' => 8.133, 'L' => 6.627, 'M' => 9.941, 'N' => 8.736,
	'O' => 9.338, 'P' => 8.133, 'Q' => 9.338, 'R' => 8.736, 'S' => 8.133, 'T' => 7.23, 'U' => 8.736, 'V' => 8.133,
	'W' => 11.447, 'X' => 8.133, 'Y' => 8.133, 'Z' => 7.23, ' ' => 3.313, '.' => 3.313, ',' => 3.313, '-' => 3.916,
	'_' => 6.627, '(' => 3.916, ')' => 3.916, '1' => 6.627, '2' => 6.627, '3' => 6.627, '4' => 6.627, '5' => 6.627,
	'6' => 6.627, '7' => 6.627, '8' => 6.627, '9' => 6.627, '0' => 6.627, '!' => 3.313, '\'' => 3.313, '/' => 3.916,
	'\\' => 3.916, '&' => 6.627
);
sub getStandardEnzymes
{
    return %standardEnzymes;
}

sub calcGCcontent
{
    my $sequenz = shift;
    my @hits    = $sequenz =~ /(g|c)/gi;
    return ((scalar @hits / length($sequenz)) * 100);
}

sub count_digits
{
	my $number = shift;
	my @digits = $number =~ /\d/g;
	return scalar(@digits);

}

sub format_seq
{
	my %args = @_;	
	_check_args(\%args, qw/-linewidth -blocklength -seq -linenumbers -frame -strand -linespacing/);

    my ($output, @outputs);
	my $leading_pos = 0;				  
	my $seq = $args{-seq};
    $seq =~ s/ //gi;
    $seq =~ s/\n//gi;
    
    #insert 5'->3'translations
    if (defined($args{-frame}))
    {
        if (is_element_of(1,@{$args{-frame}}))
        {
            push @outputs, _translate(1,'',$seq);
            $leading_pos++;
        }
        if (is_element_of(2,@{$args{-frame}}))
        {
            push @outputs, _translate(1,' ',substr($seq,1));
            $leading_pos++;
        }
        if (is_element_of(3,@{$args{-frame}}))
        {
            push @outputs, _translate(1,'  ',substr($seq,2));
            $leading_pos++;
        }        
    }
	push @outputs, $seq;
	
    #insert complementary strand
    if ((defined($args{-strand})) and  ($args{-strand} eq 'both'))
    {
        my $compSeq = complementary($seq);
        push @outputs, $compSeq;
    }
    $args{-blocklength} = $args{-linewidth} 
		if (!defined($args{-blocklength}) and $args{-linewidth});

    #insert 3'->5' translations
    if (defined($args{-frame}))
    {
    	my @revstuffer;
    	my $l = (length($seq) % 3);
    	#print "------->$l\n";
    	if ($l == 0) 
    	{
    		$revstuffer[0] = '';
    		$revstuffer[1] = '  ';
    		$revstuffer[2] = ' ';
    		
    	} elsif ($l == 1) {
    		$revstuffer[0] = ' ';
    		$revstuffer[1] = '';
    		$revstuffer[2] = '  ';
    		
    	} elsif ($l == 2) {
    		$revstuffer[0] = '  ';
    		$revstuffer[1] = ' ';
    		$revstuffer[2] = '';
    	}
    	
        if (is_element_of(-1,@{$args{-frame}}))
        {
            my $s =  reverse _translate(-1,$revstuffer[0],rev_com($seq));
            push @outputs, $s;
        }
        if (is_element_of(-2,@{$args{-frame}}))
        {
            my $s = reverse _translate(-1,$revstuffer[1],substr(rev_com($seq),1));
            push @outputs, $s;
        }
        if (is_element_of(-3,@{$args{-frame}}))
        {
            my $s = reverse _translate(-1,$revstuffer[2],substr(rev_com($seq),2));
            push @outputs, $s;
        }
    }
    
    #insert line breaks
    for (my $c= 0;$c <= $#outputs;$c++)
    {
        $outputs[$c] = break_sequence(-seq => $outputs[$c], -length => $args{-linewidth});
    }

    #combine all lineArrays into one output string
    my @lineArrays;
    for (my $c= 0;$c <= $#outputs;$c++)
    {
        
        my @tmp = split /\n/, $outputs[$c];
        push @lineArrays, \@tmp;
    }
    
    #insert block spacing
    unless ($args{-blocklength} == $args{-linewidth})
	{
	    for (my $q=0; $q <= $#lineArrays; $q++)
	    {
	    	$lineArrays[$q] = _split_blocks($args{-blocklength}, $lineArrays[$q]);
	    }
	}

    # insert extra newline in last output line to improve readability
    # optional
    if ((defined($args{-linespacing})) and ($args{-linespacing} == 1))
    {
        for (my $a =0; $a <= $#{$lineArrays[-1]}; $a++)
        {
            ${$lineArrays[-1]}[$a].="\n";
        }
    }
    
    
    #insert line numbers before leading strand
    my $pad ='   ';
    for (1..count_digits(length($seq))) {$pad .= ' '}    
    if (defined($args{-linenumbers}) and ($args{-linenumbers} == 1))
        {
            for (my $c=0;$c <= scalar($#{$lineArrays[$leading_pos]}); $c++)
            {
                ${$lineArrays[$leading_pos]}[$c] = 
                    sprintf("%".(count_digits(length($seq))+1)."d: ",($c*$args{-linewidth}+1)).${$lineArrays[$leading_pos]}[$c];

                for (my $f = 0; $f <= scalar($#lineArrays); $f++)
                {
                    next if ($f == $leading_pos);
                    ${$lineArrays[$f]}[$c] = $pad.${$lineArrays[$f]}[$c];
                }
            }
        }
    
    #combining continued
    for (my $d=0; $d <= (scalar(@{$lineArrays[0]})-1); $d++)
    {
        for (my $e=0; $e <= $#lineArrays; $e++)
        {
            $output .= ${$lineArrays[$e]}[$d]."\n";
        }
    }

	return $output;
}

sub _check_args
{
    my $args    = shift;
    my @allowed = @_;
    foreach my $arg (keys(%$args))
    {
        croak("Unknown argument: $arg\n")
          unless (&is_element_of($arg, @allowed));
    }
}           #end sub _check_args

sub is_element_of
{
    my $test  = shift;
    my @array = @_;
    foreach my $e (@array)
    {
        return TRUE if ($e eq $test);
    }
    return FALSE;
}    #end sub

sub is_same_feature
{
#	print "ARGS_TO_TESTROUTINE: @_\n";   
    my $testfeat  = shift;
    my @featarray = @_;
    
#    print Dumper ($testfeat);
#    print Dumper (@featarray);
#    exit(1);
    foreach my $e (@featarray)
    {
        return TRUE if ( ($e->start == $testfeat->start) &&
        				 ($e->end == $testfeat->end) );
    }
    return FALSE;
}    #end sub

sub break_sequence
{
	my %args = @_;
	_check_args(\%args, qw/-seq -length/);
	my ($length, $sequence) = ($args{-length}, $args{-seq});
	my @seqArray = split //, $sequence;
	my $chunks = int(length($sequence) / $length);
	
	for (my $c = 1; $c <= $chunks; $c++)
	{
		$seqArray[$c*$length-1] .= "\n";		
	}
	return join("", @seqArray);
}

sub unique #http://www.rocketaware.com/perl/perlfaq4/How_can_I_extract_just_the_uniqu.htm
{
    my @orig_array = @_;
    my @out;
    my %saw;
    undef %saw;
    @out = grep(!$saw{$_}++, @orig_array);
    return @out;
}           #end sub _u

sub complementary
{
	my $seq = shift;
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	return $seq;
}

sub rev_com
{
	my $seq = shift;
	$seq = reverse (complementary($seq));
	return $seq;
}

sub insert_element
{
	my $element = shift;
	my $position = shift;
	my @array1=@_;
	my @array2 = splice (@array1, $position);
	push @array1, $element;
	push @array1, @array2;
	return @array1;
}

sub _split_blocks
{
	my ($blocklength, $p_lineArray) = @_;
	my @lineArray = @{$p_lineArray};
	my $lc = 0;
	my @outputArray;
	foreach my $line (@lineArray)
	{
		my @chars = split //, $line;
		
		for (my $c = 1; ; ++$c)
		{
			next if ($chars[$c] eq "\n");
			
			$chars[$c-1] .= ' ' if ($c % $blocklength == 0);
			last if ($c >= $#chars);
		}
		my $newline = join("", @chars);
		push @outputArray, $newline;
		$lc++;
	}
	return \@outputArray;
}

sub _translate
{
    my ($ori, $stuffer, $seq) = @_;
    my $seqObj = Bio::Seq->new(-seq => $seq);
    my @prot = split //, $seqObj->translate()->seq();
    for (my $c=0; $c <= $#prot; $c++) 
    {
        $prot[$c] =" $prot[$c] ";
    }
    if ($ori == 1)
    {
      $prot[0]=$stuffer.$prot[0];
    } elsif ($ori == -1)
    {
      $prot[-1]=$prot[-1].$stuffer;
    }
    
    return join("", @prot);
}

sub min
{
    my @values = @_;
    my $min = $values[0];
    foreach my $value (@values)
    {
        $min = $value if ($value < $min);
    }
    return $min;
}

sub max
{
    my @values = @_;
    my $max = $values[0];
    foreach my $value (@values)
    {
        $max = $value if ($value > $max);
    }
    return $max;
}

sub print_in_columns
{
    my $no_columns = shift;
    my @fields = @_;
    my $c = 1;
    my $output = '';
    foreach my $text (@fields)
{
        $output .= $text;
        if ($c % $no_columns == 0)
{
            $output .= "\n";
} else {
            $output .= "\t";
}
        $c++;
}
    return $output;
}

sub print_in_html_columns
{
    my $no_columns = shift;
    my @fields = @_;
    my $c = 1;
    my $output = '';
    foreach my $text (@fields)
{
        $output .= $text;
        if ($c % $no_columns == 0)
{
            $output .= "<br>";
} else {
            $output .= ", <TAB>";
}
        $c++;
}
    return $output;
}

sub convert
{
	my $dna = @_;
	my $rev_com_dna = convert($dna);
	return $rev_com_dna;
}

sub trunc_with_features
{
	my ($seqobj, $start, $end, $VERBOSE) = @_;
	my $sub_seq = $seqobj->trunc($start, $end);
	$start--;
	#
	# get all features from the original sequence object
	# identify the ones that lie within the truncated
	# seq region and add them to the subsequence object
	# - difficulty: cut lies within a feature - we need
	# to generate a partial feature and mark in it the display!
	#
	my @all_features = $seqobj->get_all_SeqFeatures;
	foreach my $feature (@all_features)
	{
		next if ($feature->primary_tag =~ /source/i);
		if ($feature->location->location_type() eq "IN-BETWEEN") {
#			print "IN_BETWEEN FEATURE: ".$feature->display_name()."\n";
#			print "startpos:".$feature->start."\n";
#			print "temporarily switching to EXACT...\n";
			$feature->location->location_type("EXACT");
		}
		if (($feature->start() >= $start) && ($feature->end() <= $end))
		{
			print "Feature lies completely within trunc region: ".$feature->primary_tag()."\t".$feature->start()."..".$feature->end()."\n" 
				if ($VERBOSE);
			if ($feature->location->isa('Bio::Location::SplitLocationI'))
			{
				print "-> Split location!\n" if ($VERBOSE);
				my @sub_locs = $feature->location->sub_Location;
				
				# create a new split location
				# and add modified copies of 
				# the original sublocations
				# to it.
				my $new_split_loc = new Bio::Location::Split();
				
				foreach my $sub_loc (@sub_locs)
                {
                        next if ($sub_loc->end <= $start); # skip sublocation of it's out of range
                        printf "\toriginal location: %d..%d\n", $sub_loc->start, $sub_loc->end if ($VERBOSE);
                        printf "\tnow set to: %d..%d\n", ($sub_loc->start - $start), ($sub_loc->end - $start) if ($VERBOSE);
                        $new_split_loc->add_sub_Location(new Bio::Location::Simple( -start=>($sub_loc->start - $start),
                                                                                    -end=>($sub_loc->end - $start),
                                                                                    -strand=>$sub_loc->strand));
                        
                    
                }# foreach loc
                
                # set the feature location
                # to the new split location
                # and add the feature to the
                # sequence object
                $feature->location($new_split_loc);
				$sub_seq->add_SeqFeature($feature);
				next; # do something
			}
			$feature->start($feature->start - $start );
			$feature->end($feature->end - $start);
			$sub_seq->add_SeqFeature($feature);
			
		} elsif ((($feature->start() >= $start) && ($feature->start() <= $end)) && ($feature->end() >= $end)) {
		
			print "Feature overlaps the end of the trunc region: ".$feature->primary_tag()."\t".$feature->start()."..".$feature->end()."\n" 
				if ($VERBOSE);
			if ($feature->location->isa('Bio::Location::SplitLocationI'))
			{
				print "-> Split location!\n" if ($VERBOSE);
				my @sub_locs = $feature->location->sub_Location;
				
				# create a new split location
				# and add modified copies of 
				# the original sublocations
				# to it.
				my $new_split_loc = new Bio::Location::Split();
				
				foreach my $sub_loc (@sub_locs)
                {
                        next if ($sub_loc->end <= $start); # skip sublocation of it's out of range
                        printf "\toriginal location: %d..%d\n", $sub_loc->start, $sub_loc->end if ($VERBOSE);
                        printf "\tnow set to: %d..%d\n", ($sub_loc->start - $start), ($sub_loc->end - $start) if ($VERBOSE);
                        $new_split_loc->add_sub_Location(new Bio::Location::Simple( -start=>($sub_loc->start - $start),
                                                                                    -end=>($sub_loc->end > $end) ? ($sub_seq->length-1) : ($sub_loc->end - $start),
                                                                                    -strand=>$sub_loc->strand));
                        
                    
                }# foreach loc
                
                # set the feature location
                # to the new split location
                # and add the feature to the
                # sequence object
                $feature->location($new_split_loc);
				$sub_seq->add_SeqFeature($feature);
				next; # do something
			}
			# create truncated, marked copy of feature!
			$feature->end($sub_seq->length-1);
			$feature->start($feature->start - $start);
			$sub_seq->add_SeqFeature($feature);
			
		} elsif ((($feature->end() >= $start) && ($feature->end() <= $end)) && ($feature->start() <= $start)) {
		
			print "Feature overlaps the start of the trunc region: ".$feature->primary_tag()."\t".$feature->start()."..".$feature->end()."\n" 
				if ($VERBOSE);
			if ($feature->location->isa('Bio::Location::SplitLocationI'))
			{
				print "-> Split location!\n" if ($VERBOSE);
				my @sub_locs = $feature->location->sub_Location;
				
				# create a new split location
				# and add modified copies of 
				# the original sublocations
				# to it.
				my $new_split_loc = new Bio::Location::Split();
				
				foreach my $sub_loc (@sub_locs)
                {
                        next if ($sub_loc->end <= $start); # skip sublocation of it's out of range
                        printf "\toriginal location: %d..%d\n", $sub_loc->start, $sub_loc->end if ($VERBOSE);
                        printf "\tnow set to: %d..%d\n", ($sub_loc->start - $start), ($sub_loc->end - $start) if ($VERBOSE);
                        $new_split_loc->add_sub_Location(new Bio::Location::Simple( -start=>($sub_loc->start < $start) ? 1 :($sub_loc->start - $start),
                                                                                    -end=>($sub_loc->end - $start),
                                                                                    -strand=>$sub_loc->strand));
                        
                    
                }# foreach loc
                
                # set the feature location
                # to the new split location
                # and add the feature to the
                # sequence object
                $feature->location($new_split_loc);
				$sub_seq->add_SeqFeature($feature);
				next; # do something
			}
			# create truncated, marked copy of feature!
			$feature->start(1);
			$feature->end($feature->end - $start);
			$sub_seq->add_SeqFeature($feature);
		}
	}
	$sub_seq->desc($sub_seq->desc()." Truncated region ".$start++."..$end");
	$sub_seq->species($seqobj->species());
	return $sub_seq;
}

sub string_width_of
{
	my $string = shift;
	my $scalingfactor = shift;
	my $stringlength;
	foreach my $char (split //, $string) {
		$stringlength += ($charwidths{$char} * $scalingfactor);
	}
	return $stringlength;
	
}

sub US_format_number {
	my $number = shift;
	my @seqlen_arr = reverse (split(//, $number));
    my ($seqlength, @arr);
    my $c = 0;
    foreach my $s (@seqlen_arr)#
    {
        push @arr, "," if (($c > 0) && ($c % 3 == 0));
        push @arr, $s;
        $c++;
    }
    $seqlength = reverse (join("",@arr));
    return $seqlength;
}

sub isVNTIfile {
	my $file = shift;
	open FILE,  $file or die "Could not open file: $!\n";
	while (<FILE>) {
		chomp;
		if ($_ =~ /This file is created by Vector NTI/i) {
			return TRUE;
		}		
	}
	close FILE;
	return FALSE;
	
}

1;

