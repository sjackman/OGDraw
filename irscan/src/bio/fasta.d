
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

module bio.fasta;

private
{
	import std.stdio;
    import std.string;
    import std.regexp;
    import std.conv;
    import std.file;
    import bio.sequence;
	import bio.base;
}


class Fasta : Sequence
{
	private char[] description;
	//private bio.sequence sequence;
	
	/*
		constructor for inline
		Fasta entry definition
	*/
	this(char[] desc, char[] seq)
	{
		super(seq);
		this.description = desc;
	}
	
	/*
		dummy constructor
	*/
	this()
	{
		super("");
	}
	
	/*
		overloaded constructor
		if sequence is read 
		from file
	*/
	this(char[] filename)
	{	
		super("");
		this.read(filename);
	}
	
	/*
	void retrieve (char[] accno)
	{
		if (this.seq.length > 0) 
			throw new bioException("Trying to overwrite existing data in object!\nPlease instantiate an empty new object for DB retrieval.");
		_parse(bio.base.http_get("www.ncbi.nlm.nih.gov",
						  "/entrez/viewer.fcgi?db=nucleotide&qty=1&c_start=1&list_uids="~accno~"&uids=&dopt=fasta&dispmax=5&sendto=t&from=begin&to=end"));
	}
	*/
	
	
	void read (char[] filename)
	{
		debug writefln("reading in fasta file: ", filename);
		char[] raw_sequence = "";
	
		try
		{
			raw_sequence = cast(char[])std.file.read(filename);
		}
		catch (FileException xy)
		{
			writefln();
			writefln("A file exception occured: " ~ xy.toString());
			throw new bioException("");
		}
		catch (Exception xx)
		{
			writefln();
			writefln("An exception occured: " ~ xx.toString());
			throw new bioException("");
		}
		_parse(raw_sequence);
		
	}
	
	int write(char[] path)
	{
		try
        {
        std.file.write(path, ">" ~ this.description ~ "\n");
        std.file.append(path, break_line(this.seq(),80));
        }
        catch (FileException xy)
        {
        	throw new bioException("A file exception occured: " ~ xy.toString());
        }
        catch (Exception xx)
        {
        	throw new bioException("A file exception occured: " ~ xx.toString());
        }
		return 0;
	}
	
	Sequence seqobj ()
	{
		return new Sequence(this.seq());
	}
	
// 	Fasta revcom ()
// 	{
// 		return new Fasta("revcom of "~this.desc, bio.base.revcom(this.seq));
// 	}
	
	char[] desc()
	{
		return this.description;
	}
	
	private void _parse(char[] raw_sequence)
	{
		int pos = std.regexp.find(raw_sequence, "^>.*", "");
		debug writefln("MATCH start: ", pos);
		while (cast(char)raw_sequence[pos] != '\n')
		{
			this.description ~= raw_sequence[pos];
			pos++;
		}
		debug writefln("MATCH end: ", pos);
		
		raw_sequence = std.regexp.sub(raw_sequence, "^>.*", "", "g");
		this.description = std.regexp.sub(this.description, ">", "", "g");
		
		/* 
			strip newlines
			convert to uppercase
			and return
		*/
		raw_sequence = std.string.toupper(raw_sequence);
		raw_sequence = std.regexp.sub(raw_sequence, "\n", "", "g");
		this.setSeq(raw_sequence);
	}
	
}
