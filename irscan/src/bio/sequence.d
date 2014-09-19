
/***************************************************************************
 *    Copyright (C) 2007 by Marc Lohse,   <lohse@mpimp-golm.mpg.de>        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

module bio.sequence;

private 
{
    import std.stdio;
    import std.string;
    import std.regexp;
    import std.conv;
    import std.file;
    import bio.base;
}


/*
	globals & constants
*/


const char [] AminoAcids = "ACDEFGHIKLMNPQRSTVWY*";
const char [] Nucleotides = "ACGT";
const char [] RNA_Nucleotides = "ACGU";


interface SequenceI
{
	void setSeq(char[] newseq);
	char[] seq ();
	char[] subSeq(long start, long end);
	bit is_dna();
	bit is_rna();
	char[] type();
	char[] formatted (int linewidth);
	char[] composition();
	Sequence translate(byte frame);
	Sequence revcom();
}

class Sequence : SequenceI
{
	private char[] alphabet;
	private char[] seq_string;
	
	this (char[] initseq)
	{
		initseq = std.string.toupper(initseq);
		initseq = std.regexp.sub(initseq, "\n", "", "g");
		_check_sequence(initseq);
		this.alphabet = _guess_alphabet(initseq);
		this.seq_string = initseq;
		
	}
	
	void setSeq(char[] newseq)
	{
		newseq = std.string.toupper(newseq);
		newseq = std.regexp.sub(newseq, "\n", "", "g");
		_check_sequence(newseq);
		this.alphabet = _guess_alphabet(newseq);
		this.seq_string = newseq;
	}
	
	Sequence revcom()
	{
		return new Sequence(bio.base.revcom(this.seq()));
	}
	
	char[] seq ()
	{
		return this.seq_string;
	}
	
	char[] subSeq(long start, long end)
	{
		if ((start < 1) || (end < 1) || (start > this.seq_string.length) || (end > this.seq_string.length) || (end < start))
		{
			throw new bioException("Subsequence boundaries out of sequencence range or endindex < startindex.");
		}
		return this.seq_string[start..end];
	}
	
	/*
		guess type via
		complexity of 
		alphabet ... plus
		checking of occurring
		characters. This
		should give a fairly
		confident guess
	*/
	private char [] _guess_alphabet(char[] seq)
	{
		char[int] composition;
		foreach (char character; seq)
		{
			composition[character]++;
		}
		debug writefln("---COMPOSITION: ", composition.length, " ", cast(char[])composition.keys );
		
		if (composition.length > 4)
		{
			return "protein";
		} else {
			foreach (char c; cast(dchar[])composition.keys )
			{
				debug writefln(".", c);
				if ((std.string.find(Nucleotides,c) == -1)
					&& (std.string.find(RNA_Nucleotides,c) == -1))
				{
					debug writefln("--------------_>positive");
					return "protein";
				}
				if (std.string.find(RNA_Nucleotides,c) >= 0)
				{
					return "RNA";
				}	
			}
			return "DNA";
		}
	}
	
	private bit _check_sequence(char[] seq)
	{
		foreach (uint index, char n; cast(char[]) seq)
		{
            if ((std.string.find(AminoAcids,n) == -1)
					&& (std.string.find(RNA_Nucleotides,n) == -1))
            {
            	writefln("Trying to set sequence to [char<"~n~">, pos "~std.string.toString(index)~"] which is neither amino acid nor nucleotide.");
            	throw new bioException("Illegal seq character");
            }
		}
		return true;
	}
	
	bit is_dna()
	{
		if (this.alphabet == "DNA")
		{
			return true;
		} else {
			return false;
		}
	}
	
	bit is_rna()
	{
		if (this.alphabet == "RNA")
		{
			return true;
		} else {
			return false;
		}
	}
	
	char[] type()
	{
		return this.alphabet;
	}
	
	char[] composition()
	{
		char[] report;
		/*
			... generate composition report
		*/
		return report;
	}
	
	char[] formatted (int linewidth)
	{
		char [] formatseq;
		formatseq = break_line(this.seq_string, linewidth);
		return formatseq;
	}
	
	Sequence translate(byte frame) // returns a new bio.sequence object
	{
		Sequence protein;
		
		if ((frame ==0) || (frame < -3) || (frame > 3))
			throw new bioException("unknown frame. Possible translation frames are -3|-2|-1|1|2|3");
		
		
		if (this.alphabet == "protein") 
		{
			// will be changed in the future to 
			// provide backtranslation functionality
			throw new bioException("can't translate a protein sequence to a protein sequence");
		}
		
		if (frame > 0 )
		{
			protein = new Sequence(DNA_to_prot(this.seq_string[frame-1..this.seq_string.length]));
		} else {
			frame *= -1;
			protein = new Sequence(DNA_to_prot(bio.base.revcom(this.seq_string)[frame-1..this.seq_string.length]));
		}
		return protein;
	}
	
	
}

unittest
{
	writefln("UNITTEST bio.sequence\n------------------");
	Sequence seq = new Sequence("AGCATAGACAGACGACTACAGCATCAGCATC");
	writefln("construction of Sequence %s ...OK", seq.seq());
	writefln(">>UNITTEST bio.sequence successful\n");
	// to be continued...
	
}