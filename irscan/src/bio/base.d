
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

module bio.base;

private 
{
    import std.stdio;
    import std.string;
    import std.regexp;
    import std.conv;
    import std.file;
    import std.socket;
	import std.socketstream;
}

const int WEBPORT=80;

class bioException : Exception
{
	this(char[] message)
	{
		super(message);
	}
}


char[] guess_format(char[] file)
{
	char[] raw_data = "";
	
	try
	{
		raw_data = cast(char[])std.file.read(file);
	}
	catch (FileException xy)
	{
		writefln();
		writefln("BIO.BASE.GUESS_FORMAT: A file exception occurred: " ~ xy.toString());
		throw new bioException("");
	}
	catch (Exception xx)
	{
		writefln();
		writefln("BIO.BASE.GUESS_FORMAT: An exception occurred: " ~ xx.toString());
		throw new bioException("");
	}
	
	/*
		examine file content
		... not very thorough
	*/
	
	char[][] lines = std.string.split(raw_data, "\n");
	debug writefln("line[0]:%s\nmatch pos:%d",lines[0], std.regexp.find(lines[0], "^>.+", ""));
	
	if (std.regexp.find(lines[0], "^>.+", "") >= 0)
	{
		return "fasta";
	} else if (std.regexp.find(lines[0], "^LOCUS.+", "i") >= 0) {
		return "genbank";
	} else {
		return "unknown";	
	}
}

char[] revcom (char[] seq)
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
				throw new bioException("found illegal character "~cast(char)i~" at position "~cast(char)index~"\n");
				
		}
	}			
	return revcom.reverse;	 
}

char [] break_line (char[] seq, int length)
{
	char[] broken_seq = "";
	long index = 0;
	foreach (char c; seq)
	{
		if ((index > 0) && (index % length == 0))
		{
			broken_seq ~= '\n';
		}
		broken_seq ~= c;
		index++;
	}
	return broken_seq;
}


/***
	This is not a real "bio" function
	as it just executes http GET requests
	but it is nevertheless handy
	args: 	host address (char[])
			path (char[])
	returns:http response (wchar[]) 

char[] http_get(char[] host, char[] path)
{
	// open a socket:
  	InternetHost ih = new InternetHost;
  	ih.getHostByName(host);
  	
  	debug(http_get) printf("addrList.length = %d\n", ih.addrList.length);
  
  	InternetAddress ia = new InternetAddress(ih.addrList[0], WEBPORT);
  	
    debug(http_get) 
    {
        printf("IP address = %.*s\nname = %.*s\n", ia.toAddrString(), ih.name);
        foreach(int i, char[] s; ih.aliases)
        {
            printf("aliases[%d] = %.*s\n", i, s);
        }
        
        printf("---\n");
    
        assert(ih.getHostByAddr(ih.addrList[0]));
        printf("name = %.*s\n", ih.name);
        foreach(int i, char[] s; ih.aliases)
        {
            printf("aliases[%d] = %.*s\n", i, s);
        }
	}
	
   	TcpSocket sock = new TcpSocket();
   	sock.connect(ia);
   
   	// send the HTTP request
   	sock.send("GET " ~ path ~ " HTTP/1.0\n\n");
         
	// read  the result:
  	char[] response ="";
  	char[] line;
  	char c;
  	SocketStream stream = new SocketStream(sock);
  	while (! stream.eof()) 
  	{
        if (c == '\0') break;
        line = stream.readLine();
        response ~= line ~ "\n"; 
  	}
  	// strip server response and empty lines from sequence
	response = std.regexp.sub(response, "^HTTP[0-9a-z/:,. \t-]+\n", "", "gi");
	response = std.regexp.sub(response, "^Date[0-9a-z/:,. \t-]+\n", "", "gi");
	response = std.regexp.sub(response, "^Server[0-9a-z/:,. \t-]+\n", "", "gi");
	response = std.regexp.sub(response, "^Content[0-9a-z/:,. \t-]+\n", "", "gi");
	response = std.regexp.sub(response, "^Connection[0-9a-z/:,. \t-]+\n", "", "gi");
	response = std.regexp.sub(response, "^\n", "", "gm");
	
	debug writefln("retrieved data: ", response.length * response.sizeof, " bytes");
	
	return response;
}

*/

/***
	function takes gi_ID or accession number
	and retrieves the corresponding GenBank
	entry in the specified format. The entry
	will be saved directly to the disk.
	args:	gi_ID|accession , filename, format (char[])
	returns: true on success.
	
	To retrieve entries directly into Fasta or
	GenBank objects use the respective
	XXX.retrieve(ACCNO) methods



bit get_sequence(char[] gi_ID, char[] filename = "sequence", char[] format = "fasta")
{
	
	assert(	(format == "genbank") || 
			(format == "fasta") ||
			(format == "gbwithparts")  );
			
	if ((format == "genbank") || (format == "gbwithparts")) 
	{
		filename ~= ".gb";
	} else {
		filename ~= ".fasta";
	}
	char[] host = "www.ncbi.nlm.nih.gov";
	char[] path = "/entrez/viewer.fcgi?db=nucleotide&qty=1&c_start=1&list_uids="~gi_ID~"&uids=&dopt="~format~"&dispmax=5&sendto=t&from=begin&to=end";
	char[] sequence;
	try 
	{
		sequence = http_get(host, path);
	} catch (Exception e) {
		throw new bioException("Could not retrieve sequence"~e.toString);
	}
	
	// analyse server response
	if (std.string.find(sequence, "is not found") >= 0)
	{
		//writefln("Accession number %s is not found", gi_ID);
		return false;
	}
	
	// return false if response is empty
	if (std.regexp.find(sequence, "[0-9a-z/:,. \t-]+", "gi") == -1) return false;
	
	/*
		write sequence to file
	
	try
    {
    std.file.write(filename, sequence);
    } catch (FileException xy) {
        throw new bioException("A file exception occured: " ~ xy.toString());
    }
    catch (Exception xx) {
        throw new bioException("An exception occured: " ~ xx.toString());
    }
    
	//writefln(sequence);
	return true;
}
*/

char [] DNA_to_prot (char [] dna)
{
	/+debug if ((dna.length % 3) != 0)
		throw new bioException("sequence is not a valid reading frame");+/
	
	// Standard genetic code
	char[char[]] TranslationTable;
	
	TranslationTable["GGG"] = 'G'; 
	TranslationTable["GGA"] = 'G'; 
	TranslationTable["GGT"] = 'G'; 
	TranslationTable["GGC"] = 'G'; 
	 
	TranslationTable["GAG"] = 'E'; 
	TranslationTable["GAA"] = 'E'; 
	TranslationTable["GAT"] = 'D'; 
	TranslationTable["GAC"] = 'D'; 
	 
	TranslationTable["GTG"] = 'V'; 
	TranslationTable["GTA"] = 'V'; 
	TranslationTable["GTT"] = 'V'; 
	TranslationTable["GTC"] = 'V'; 
	 
	TranslationTable["GCG"] = 'A'; 
	TranslationTable["GCA"] = 'A'; 
	TranslationTable["GCT"] = 'A'; 
	TranslationTable["GCC"] = 'A'; 
	 
	TranslationTable["AGG"] = 'R'; 
	TranslationTable["AGA"] = 'R'; 
	TranslationTable["AGT"] = 'S'; 
	TranslationTable["AGC"] = 'S'; 
	 
	TranslationTable["AAG"] = 'K'; 
	TranslationTable["AAA"] = 'K'; 
	TranslationTable["AAT"] = 'N'; 
	TranslationTable["AAC"] = 'N'; 
	 
	TranslationTable["ATG"] = 'M'; 
	TranslationTable["ATA"] = 'I'; 
	TranslationTable["ATT"] = 'I'; 
	TranslationTable["ATC"] = 'I'; 
	 
	TranslationTable["ACG"] = 'T'; 
	TranslationTable["ACA"] = 'T'; 
	TranslationTable["ACT"] = 'T'; 
	TranslationTable["ACC"] = 'T'; 
	 
	TranslationTable["TGG"] = 'W'; 
	TranslationTable["TGA"] = '*'; 
	TranslationTable["TGT"] = 'C'; 
	TranslationTable["TGC"] = 'C'; 
	 
	TranslationTable["TAG"] = '*'; 
	TranslationTable["TAA"] = '*'; 
	TranslationTable["TAT"] = 'Y'; 
	TranslationTable["TAC"] = 'Y'; 
	 
	TranslationTable["TTG"] = 'L'; 
	TranslationTable["TTA"] = 'L'; 
	TranslationTable["TTT"] = 'F'; 
	TranslationTable["TTC"] = 'F'; 
	 
	TranslationTable["TCG"] = 'S'; 
	TranslationTable["TCA"] = 'S'; 
	TranslationTable["TCT"] = 'S'; 
	TranslationTable["TCC"] = 'S'; 
	 
	TranslationTable["CGG"] = 'R'; 
	TranslationTable["CGA"] = 'R'; 
	TranslationTable["CGT"] = 'R'; 
	TranslationTable["CGC"] = 'R'; 
	 
	TranslationTable["CAG"] = 'Q'; 
	TranslationTable["CAA"] = 'Q'; 
	TranslationTable["CAT"] = 'H'; 
	TranslationTable["CAC"] = 'H'; 
	 
	TranslationTable["CTG"] = 'L'; 
	TranslationTable["CTA"] = 'L'; 
	TranslationTable["CTT"] = 'L'; 
	TranslationTable["CTC"] = 'L'; 
	 
	TranslationTable["CCG"] = 'P';
	TranslationTable["CCA"] = 'P';
	TranslationTable["CCT"] = 'P';
	TranslationTable["CCC"] = 'P';	
	char [] translated, codon = "";
	for (long i = 0; i <= (dna.length -3); i+=3)
	{
		codon = dna[i..i+3];
		translated ~= TranslationTable[codon];
		//writefln(codon, "->",TranslationTable[codon]);
	}
	return translated;
}
