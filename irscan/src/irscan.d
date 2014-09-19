/*
    irscreen
    ---------------------
    Program to find inverted 
    repeat regions
    in chloroplast genomes (and
    other DNA). Since it is 
    optimized for chloroplast
    genomes it will stop
    after finding the first 
    inverted reapeat pair.
    Part of the GeneMap::Plastome
    package (but written in D for
    better performance)
*/

import std.stdio;
import std.file;
import std.string;
import std.conv;
import getopt;
import bio.base;
import bio.sequence;
import bio.genbank;
import bio.fasta;

/*
	global constants
	and variables
*/

long WORD_SIZE = 2000; // min length of IR
long PRECISION = 1;	// not used anymore
bit  VERBOSE = false;
char[] help_msg =
"
Usage: irscan [-f|--file=<file>] {options}
 
Options:
        -v|--verbose                     Set program output to verbose level
        -w|--wordsize=<WORDSIZE>         Define minimal size for inverted repeat. Defaults to 2000
        -p|--precision=<PRECISION>       Set precision of algorithm. Default value is 1
        
          ";

struct Repeat
{
	long start;
	long end;
}

/*
    function takes a 
    char[] DNA sequence
    as argument and returns
    char[] the reverse 
    complement of the 
    sequence
*/

char[] rev_com (char[] seq)
{
	char[] revcom = "";
	foreach (int index, char i ; seq)
	{
		switch (i)
		{
			case 'A' :
				revcom ~= "T";
				break;
			case 'T' :
				revcom ~= "A";
				break;
			case 'C' :
				revcom ~= "G";
				break;
			case 'G' :
				revcom ~= "C";
				break;
			default :
				writefln("found illegal character ",i," at position ",index);
				return "";
		}
	}			
	return revcom.reverse;	 
}

/*
	Function takes a DNA string
	and calulates the GC content
*/

float gc_cont (char[] sequence)
{
	float gc_cont;
	long gc_num = 0;
	for (long i = 1; i < sequence.length; i++)
	{
		if (cast(char)sequence[i] == 'G' || cast(char)sequence[i] == 'C') 
		{
			gc_num++;
		}
	}
	gc_cont = (cast(float)gc_num / sequence.length) * 100;
	//writefln("gc_num: ", gc_num, " len: ", sequence.length, " gccont: ", (gc_num / sequence.length));
	return gc_cont;
}

/*
	Function takes a
	filename and loads the
	DNA sequence from the
	file. Only raw text 
	format supported so far
*/

char[] load_seq(char[] file)
{
	char[] raw_sequence = "";
	
	try
	{
		raw_sequence = cast(char[])read(file);
	} catch (FileException xy) {
		writefln();
		writefln("A file exception occured: " ~ xy.toString());
	}
	catch (Exception xx)
	{
		writefln();
		writefln("An exception occured: " ~ xx.toString());
	}
	
	/* 
		strip newlines
		convert to uppercase
		and return
	*/
	raw_sequence = std.string.toupper(raw_sequence);
	return std.regexp.sub(raw_sequence, "\n", "", "g");
	
}

/*
    function to find 
    the borders
*/

long[] screen_from_to_pos(char[] seq, long pos, long endpos, int step)
{
    long start, end;
    char[] test_word;
    
    while (true)
    {
        if (pos + WORD_SIZE >= endpos) 
        {
            return [-1];
        }
        
        test_word = rev_com(seq[pos..(pos + WORD_SIZE)]);
        
        end = std.string.find(seq, test_word);		
        if (end != -1)
        {
            end += WORD_SIZE;
            start = pos + 1;
            break;
        }
        pos += step;
    }
    return [start, end];
}

/*
	Function determines
	IR borders with a
	precision of +/- 500 bp
*/

long[] rough_scan (char[] sequence)
{
	char[] test_word = "";
	long position = 0;
	long ira_start,
		 ira_end,
		 irb_start,
		 irb_end;
		 
	while (true)
	{
		if ((position + WORD_SIZE >= sequence.length) && (! VERBOSE))
		{
			writefln("false");
			return [-1];
		} else if (position + WORD_SIZE >= sequence.length) {
			writefln("This sequence does not contain\ninverted repeats longer than ", WORD_SIZE, " bp");
			return [-1];
		}
		
		test_word = rev_com(sequence[position..(position + WORD_SIZE)]);
		
		irb_end = std.string.find(sequence, test_word);		
		if (irb_end != -1)
		{
			irb_end += WORD_SIZE;
			ira_start = position;
			break;
		}
		position += 500;
	}
	
	/*
		now we try to find 
		the other two border
		positions. To save
		time we start between
		them.
	*/
	
	long ssc_center = (ira_start + irb_end) / 2;
	position = ssc_center;
	long[] inner_borders = screen_from_to_pos(sequence, position,sequence.length, 500);
	irb_start = inner_borders[0];
	ira_end = inner_borders[1];
	return [ira_start, ira_end, irb_start, irb_end];
}

void full_scan (char [] sequence)
{
	for (long i = 0; i <= sequence.length; i++)
	{
		writefln("Doing complete scan of sequence...");
	}
}

/***************************
	Main function
	-------------
	reads in first 
	command line argument
	as filename and prints
	the exact border positions 
	to seperated by ";" to STDOUT
	if no IRs are found the text
	"false" is printed.
	If option -v is passed the
	output will be more
	informative
****************************/

int main(char[][] args)
{
	args = args[1..args.length]; // remove first arg (filename) 
	if (!args.length > 0) 
	{
		writefln(help_msg);
		return 1;
	}
    getopt_t[] usedoptions;

	char [] filename, sequence, seq_format;
	Repeat IRA, IRB;
	long[] pre_borders, IR_borders;
    char[] shortoptions = "hvp:w:f:";
    char[][6]longoptions;
	int idxop = -1;
	
	longoptions[++idxop] = "help";
    longoptions[++idxop] = "verbose";
	longoptions[++idxop] = "precision=";
	longoptions[++idxop] = "wordsize=";
	longoptions[++idxop] = "file=";
	longoptions[++idxop] = "full";


    try
    {
        usedoptions = getopt.getopt(args, shortoptions, longoptions);
    } catch(GetOptException e) {
        writefln("Error in passed options: ", e.msg);
        return 1;
    }

    foreach(getopt_t item; usedoptions)
    {
        
			if (item.option == "-h" || item.option == "--help") 
			{
				writefln(help_msg);
				return 0;
			} else if (item.option == "-v" || item.option == "--verbose") {
				VERBOSE = true;
				
			} else if (item.option == "-p" || item.option == "--precision") {
				PRECISION = std.conv.toLong(item.value);
				
			} else if (item.option == "-w" || item.option == "--wordsize") {
				WORD_SIZE = std.conv.toLong(item.value);
				
			} else if (item.option == "-f" || item.option == "--file") {
				filename = item.value;
							
			} else if (item.option == "--full") {
				full_scan(sequence);
				return 0;
			}
	}	
	
	seq_format = bio.base.guess_format(filename);
	debug writefln("Sequence format is %s", seq_format); 
	if (seq_format == "fasta")
	{
		Fasta tmp_seq = new Fasta(filename);
		sequence = tmp_seq.seq();
		delete tmp_seq;
		
	} else if (seq_format == "genbank") {
		GenBank tmp_seq = new GenBank(filename);
		sequence = tmp_seq.seq;
		delete tmp_seq;
		
	} else {
		if (VERBOSE) writefln("Unknown sequence file format. Known are [fasta|genbank]\nAssuming sequence is raw text/DNA\n");
		sequence = load_seq(filename);
	}
	
	if (VERBOSE) 
	{
		writefln("--------------------------");
		writefln("loading sequence from file...OK\n","length:\t\t", sequence.length, " bp");
		writefln("word size:\t", WORD_SIZE);
		writefln("precision:\t", PRECISION);
	}
	
	pre_borders = rough_scan(sequence);
	
	//no IRs? -> exit program
	if (pre_borders[0] == -1) return 1;
	
	if (VERBOSE) writefln ("rough scan\t",pre_borders);
	
	IR_borders = screen_from_to_pos(sequence, pre_borders[0]-500, sequence.length, 1);
	IRA.start = IR_borders[0];
	IRB.end = IR_borders[1];
	
	IR_borders = screen_from_to_pos(sequence, pre_borders[2]-500, sequence.length, 1);
	IRA.end = IR_borders[1];
	IRB.start = IR_borders[0];
	
	if (! VERBOSE) 
	{
		writefln (IRA.start,";",IRA.end,";", IRB.start,";",IRB.end, ";", 
					gc_cont(sequence[IRA.start..IRA.end]), ";", //GC IRs
					gc_cont(sequence[IRA.end+1..IRB.start-1]), ";", //GC SSC
					gc_cont(sequence[1..IRA.start-1]) //GC LSC
					);
	}
	else {
		writefln("\nInverted repeat A\t\t%6d..%6d\tlength:%6d bp\tGC: %3.3f", IRA.start, IRA.end, (IRA.end - IRA.start), gc_cont(sequence[IRA.start..IRA.end]));
		writefln("Inverted repeat B\t\t%6d..%6d\tlength:%6d bp\tGC: %3.3f", IRB.start, IRB.end, (IRB.end - IRB.start), gc_cont(sequence[IRB.start..IRB.end]));
		writefln("Small single copy region\t%6d..%6d\tlength:%6d bp\tGC: %3.3f", (IRA.end+1), (IRB.start-1), ((IRB.start-1)-(IRA.end+1)) , gc_cont(sequence[IRA.end+1..IRB.start-1]));
		writefln("Large single copy region\t%6d..%6d\tlength:%6d bp\tGC: %3.3f", 1, (IRA.start-1), (IRA.start-1), gc_cont(sequence[1..IRA.start-1]));	
		writefln("--------------------------");
	}
		
	return 0;
}