
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

module bio.genbank;

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


class GenBank : Sequence
{
	// this structure will hold all info on one feature
	struct feature
	{
		char[] 	ID;
		long 	start;
		long 	end;
		byte	strand;
		char[]	type;
		char[char[]] tags;
	}
	
	/*
		properties - only accessible via corr. methods
	*/
	
	private 
	{
		int _feature_counter = -1;
		char[] description;
        char[] _accnum;
        char[] _type; // _type of whole sequence (e.g. mRNA, DNA, cDNA, EST etc.)
        char[] _organism; 
        char[][] _keywords;
        char[] _date;
        bit _circular;
        feature[500] features;  // max 500 features....might be to little
	}
	/*..... to be continued*/
	
	/*
		constructor for inline
		entry definition
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
		char[] seq = this.read(filename);
		super(seq);
	}
	
	bit addFeature (GenBank.feature feat)
	{
		// add a feature to the list
		return true;
	}
	
	/***
	* Accessory methods
	*/
	
	char[] desc()
	{
		return this.description;
	}
	
	char[] type()
	{
		return this._type;
	}
	
	char[] accnum()
	{
		return this._accnum;
	}
	
	char[] date()
	{
		return this._date;
	}
	
	char[] organism()
	{
		return this._organism;
	}
	
	char[][] keywords()
	{
		return this._keywords;
	}
	
	bit circular()
	{
		return this._circular;
	}
	
	
	feature getFeature(int index)
	{
		if (index > this._feature_counter)
		{
			throw new bioException("Your're trying to access feature no."~std.conv.toString(index)~
				". However there are only "~std.conv.toString(this._feature_counter+1)~
				" features. Hence the last index position is "~std.conv.toString(this._feature_counter));
		}
		return this.features[index];
	}
	/**
	* getAllFeatures()
	* DESC: get all features stored in the bio.genbank object
	* ARGS: none
	* RETURNS: an array of bio.genbank.feature objects
	*/
	
	feature[] getAllFeatures()
	{
		return this.features[0..this._feature_counter];
	}
	
	/***
	* read()
	* DESC: Reads a genbank file into
	* a bio.genbank.GenBank object
	* ARGS: char[]
	* RETURNS: 	raw sequence.
	*/
	
	char[] read (char[] path)
	{
		char[] 	raw_sequence, 
				stripped_sequence, 
				header_segment, 
				feature_segment, 
				intermediate,
				seq_segment;
		
		try
		{
			raw_sequence = cast(char[])std.file.read(path);
		}
		catch (FileException xy)
		{
			writefln();
			writefln("A file exception occurred: " ~ xy.toString());
			throw new bioException("");
		}
		catch (Exception xx)
		{
			writefln();
			writefln("An exception occurred: " ~ xx.toString());
			throw new bioException("");
		}
	
		return _parse(raw_sequence);
	}
	
	/***
	* retrieve()
	* DESC: Takes an accession number or genInfo 
	* identofier as argument and directly reads
	* the corresponding entry into a bio.genbank.GenBank
	* object.
	* ARGS: char[]	
	* RETURNS: 	true on success.
	*/
	/*
	bit retrieve(char[] accno)
	{
		
		if (this.seq().length > 0) 
			throw new bioException("Trying to overwrite existing data in object!\nPlease instantiate an empty new object for DB retrieval.");
			
		this.setSeq(_parse(bio.base.http_get("www.ncbi.nlm.nih.gov",
						  "/entrez/viewer.fcgi?db=nucleotide&qty=1&c_start=1&list_uids="~accno~"&uids=&dopt=genbank&dispmax=5&sendto=t&from=begin&to=end")));
		
		return true;
	}
	*/
	
	/*---INTERNAL FUNCTIONS---------*/
	
	private char[] _parse(char[] raw_sequence)
	{
		char[] 	stripped_sequence, 
				header_segment, 
				feature_segment, 
				intermediate,
				seq_segment;
				
		auto segments = std.string.split(raw_sequence, "FEATURES");
		header_segment = segments[0];
		auto segments2 = std.string.split(segments[1], "ORIGIN");
		feature_segment = segments2[0];
		seq_segment = segments2[1];
		
		//stripped_sequence =	sub(seq_segment, "[0-9\n\t/ ]", "", "g");
		stripped_sequence = std.string.toupper(seq_segment);
		
		/*The regex above performs SEVERAL HUNDRED FOLD SLOWER than this ugly crap*/
		char[] clean_seq = "";
		foreach (char N; stripped_sequence)
		{
			if ((N == '0') || (N == '1') || (N == '2') || (N == '3') ||
				(N == '4') || (N == '5') || (N == '6') || (N == '7') ||
				(N == '8') || (N == '9') || (N == ' ') || (N == '\n') ||
				(N == '\t') || (N == '/')) continue;
			clean_seq ~= N;
		}
		

		this.setSeq(clean_seq);
		debug 
		{
            writefln("HEADER-------------------\n", header_segment);
            writefln("FEATURES-------------------\n", feature_segment);
            writefln("SEQUENCE-------------------\n", clean_seq);
        }
		
		_parse_header(header_segment);
		//_parse_features(feature_segment); // temporarily deactivated... bugs
		
		return clean_seq;
	}
	

	private void _parse_header (char[] header)
	{
		char[][] lines = std.regexp.split(header, "\n", "");
		foreach (int index, char[] line; lines)
		{
			if (auto m = std.regexp.search(line, "locus[ \t]*", "i"))
			{
				char[][] words = std.regexp.split(m.post, "[ \t]+", "");
				this._accnum = words[0];
				this._type = words[3];
				this._date = words[6];
				if (words[4] == "linear")
				{
					this._circular = false;
				} else {
					this._circular = true;	
				}
			}
			
			if (auto m = std.regexp.search(line, "organism[ \t]*", "i"))
			{
				this._organism = m.post;
			}
			
			if (auto m = std.regexp.search(line, "^keywords", "i"))
			{
				auto m2 = std.regexp.search(m.post, "SOURCE", "");
				//char[] keywords = sub(m2.pre, "\n", "", "g");
				//this.keywords = std.string.split(keywords, ";");
			}
			
			// this doesn't work for multiline DEFINITIONs !! --- needs refinement
			if (auto m = std.regexp.search(line, "definition[ \t]*", "i"))
			{
				this.description = m.post;
			}
		}
	}
	
	private void _parse_features (char[] features)
	{
		/*
			Parse the features and
			extract info on each feature
			into a feature struct 
		*/
		
		//split features into lines
		char[][] lines = std.regexp.split(features, "\n", "");
		long i, chunks, pos_idx;
		long[500] positions;
		
		pos_idx = -1;
		i = 0;
		chunks = 0;
		
		while (i <= (lines.length -1))
		{
			char[] line = lines[i];
			// replace leading whitespaces
			line = std.regexp.sub(line, "^[\t ]+", "", "");
			
			char[][] words = std.regexp.split(line, "[\t ]+", "g");
			//writefln(words);
			if (is_feature_key(words[0]))
			{
				//store positions of feature starts
				try
				{
					positions[++pos_idx] = i;
				} catch (Exception e) {
					throw new bioException("To many features (>500): " ~ e.toString());
				}
			}
			i++;
		}
		for (long c = 0; c <= positions.length; c++)
		{
			long endpos = positions[c+1];
			if (positions[c] > positions[c+1])
			{
				endpos = lines.length;
			}
			// slice line array into feature chunks
			char[][] raw_feature = lines[positions[c]..endpos];
			debug(featparse)writefln("featurechunk: -----------------------\n", raw_feature, "\n-----------------------------------\n");
			chunks++;
			this._feature_counter++;
			/*
				process feature chunk
			*/
			_process_feature_chunk(raw_feature);
			if (positions[c] > positions[c+1])
				break;
		}
		debug(featparse) writefln("No of chunks: ", chunks);
		this._feature_counter++;
	}
	
	private void _process_feature_chunk(char[][] raw_feat)
	{
		feature feat;
		int feat_idx = -1;
		char[] line1 = std.regexp.sub(raw_feat[0], "^[\t ]+", "", "");
		char[][] firstline = std.regexp.split(line1, "[\t ]+", "g");
		//writefln(raw_feat, "----------------------------"); 
		//writefln(firstline);
		feat.type = firstline[0];
		
		if (auto m = std.regexp.search(firstline[1], "join|complement", "ig"))
		{
			writefln("found complementary and/or composite feature: %s <%s> %s", m.pre, m.match(0), m.post);
			writefln("<<<<<<<<<<<<<<<<<<<<<<< %s|%s",m.match(0), m.match(1));
			
			if ( (std.string.find(firstline[1], "complement") != -1) && 
				 (std.string.find(firstline[1], "join") != -1) )
			{
				// process complement composite
				
				auto m1 = std.regexp.sub(firstline[1], "^complement.", "", "");
				auto m2 = std.regexp.sub(m1, ".$", "", "");
				auto m3 = std.regexp.sub(m2, "<|>", "", "g");
				//writefln("----------------------", m3);
				
				
				feat.strand = -1;
				goto TAG_PROCESSING;
				
			} else if (m.match(0) == "complement") {
				// process complement
				auto m3 = std.regexp.sub(firstline[1], "<|>", "", "g");
				auto positions = std.regexp.search(m3, "([0-9]+..[0-9]+)", "g");
				auto range = std.string.split(positions.match(0), "..");
				writefln("Cleaned line:%s\tMatch: %s\tRANGE:", m3, positions.match(0),range);
				try {
                    feat.start = std.conv.toLong(range[0]);
                    feat.end = std.conv.toLong(range[1]);
                } catch (Exception e) {
                	throw new bioException("Feature range could not be extracted:"~e.toString);
                }
				feat.strand = -1;
				goto TAG_PROCESSING;
				
			} else if (m.match(0) == "join") {
				//process composite
				feat.strand = 1;
				goto TAG_PROCESSING;
				
			} else {
				feat.strand = 1;
			}
		}
		
		feat.strand = 1;
		char[][] range = std.string.split(firstline[1], "..");
		
		try 
		{
			feat.start = std.conv.toLong(range[0]);
			feat.end = std.conv.toLong(range[1]);
		} catch (Exception e) {
			throw new bioException("Error parsing feature range: "~firstline[1]~" : "~e.toString());
		}
		
		TAG_PROCESSING:
		// process the tags
		foreach (char[] line; raw_feat[1..raw_feat.length])
		{
			char[] cleanline = std.regexp.sub(line, "^[\t ]+", "", "");
			char[][] words = std.regexp.split(cleanline, "[\t ]+", "g");
			
		}
		
		this.features[this._feature_counter] = feat;
	}
}


/*
	non member functions
*/

bit is_feature_key (char[] test)
{
	int[char[]] featureKeys;
	featureKeys["source" ] = 1;
	featureKeys["misc_feature" ] = 1;
    featureKeys["misc_difference"] = 1;
    featureKeys["conflict" ] = 1;
    featureKeys["unsure" ] = 1;
    featureKeys["old_sequence"] = 1;
    featureKeys["variation"  ] = 1;
    featureKeys["modified_base"] = 1;
    featureKeys["gene" ] = 1;
    featureKeys["misc_signal" ] = 1;
    featureKeys["promoter"] = 1;
    featureKeys["CAAT_signal" ] = 1;
    featureKeys["TATA_signal" ] = 1;
    featureKeys["-35_signal" ] = 1;
    featureKeys["-10_signal" ] = 1;
    featureKeys["GC_signal"] = 1;
    featureKeys["RBS"] = 1;
    featureKeys["polyA_signal" ] = 1;
    featureKeys["enhancer" ] = 1;
    featureKeys["attenuator" ] = 1;
    featureKeys["terminator" ] = 1;
    featureKeys["rep_origin"] = 1;
    featureKeys["oriT"] = 1;
    featureKeys["misc_RNA" ] = 1;
    featureKeys["prim_transcript" ] = 1;
    featureKeys["precursor_RNA" ] = 1;
    featureKeys["mRNA" ] = 1;
    featureKeys["5'clip" ] = 1;
    featureKeys["3'clip" ] = 1;
    featureKeys["5'UTR" ] = 1;
    featureKeys["3'UTR" ] = 1;
    featureKeys["exon" ] = 1;
    featureKeys["CDS" ] = 1;
    featureKeys["sig_peptide" ] = 1;
    featureKeys["transit_peptide" ] = 1;
    featureKeys["mat_peptide" ] = 1;
    featureKeys["intron" ] = 1;
    featureKeys["polyA_site" ] = 1;
    featureKeys["rRNA" ] = 1;
    featureKeys["tRNA" ] = 1;
    featureKeys["scRNA" ] = 1;
    featureKeys["snRNA" ] = 1;
    featureKeys["snoRNA"] = 1;
    featureKeys["C_region" ] = 1;
    featureKeys["D_segment" ] = 1;
    featureKeys["J_segment" ] = 1;
    featureKeys["N_region" ] = 1;
    featureKeys["S_region" ] = 1;
    featureKeys["V_region" ] = 1;
    featureKeys["V_segment"] = 1;
    featureKeys["repeat_region" ] = 1;
    featureKeys["repeat_unit" ] = 1;
    featureKeys["LTR" ] = 1;
    featureKeys["satellite"] = 1;
    featureKeys["misc_binding" ] = 1;
    featureKeys["primer_bind" ] = 1;
    featureKeys["protein_bind"] = 1;
    featureKeys["misc_recomb" ] = 1;
    featureKeys["iDNA"] = 1;
    featureKeys["misc_structure" ] = 1;
    featureKeys["stem_loop" ] = 1;
    featureKeys["D-loop"] = 1;
    featureKeys["gap"] = 1;
    featureKeys["operon" ] = 1;
    
    if (test in featureKeys)
    {
    	return true;
    } else {
    	return false;
    }
}

bit is_gb_keyword (char[] test)
{
	int[char[]] gb_keywords;
	gb_keywords["LOCUS"] = 1;
    gb_keywords["DEFINITION"] = 1;
    gb_keywords["ACCESSION"] = 1;
    gb_keywords["VERSION"] = 1;
    gb_keywords["KEYWORDS"] = 1;
    gb_keywords["SOURCE"] = 1;
    gb_keywords["ORGANISM"] = 1;
    gb_keywords["REFERENCE"] = 1;
    gb_keywords["AUTHORS"] = 1;
    gb_keywords["TITLE"] = 1;
    gb_keywords["JOURNAL"] = 1;
    gb_keywords["MEDLINE"] = 1;
    gb_keywords["REFERENCE"] = 1;
    gb_keywords["AUTHORS"] = 1;
    gb_keywords["TITLE"] = 1;
    gb_keywords["JOURNAL"] = 1;
    gb_keywords["FEATURES"] = 1;
    gb_keywords["ORIGIN"] = 1;
    
    if (test in gb_keywords)
    {
    	return true;
    } else {
    	return false;
    }
}
