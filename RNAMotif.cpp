// ==========================================================================
//                                  RNAMotif
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Benjamin Strauch
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/graph_types.h>

#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C"{
	#include  <ViennaRNA/data_structures.h>
	#include  <ViennaRNA/mfe.h>
	#include  <ViennaRNA/params.h>
	#include  <ViennaRNA/utils.h>
	#include  <ViennaRNA/eval.h>
	#include  <ViennaRNA/fold.h>
	#include  <ViennaRNA/part_func.h>
}

// ==========================================================================
// Classes
// ==========================================================================

//types used in the program
typedef unsigned TPosition;
typedef float TScoreValue;
typedef seqan::CharString TString;
typedef float TBioval;
typedef std::map<TPosition, TScoreValue> TMap;
typedef seqan::String<TMap > TMapLine;

typedef float TCargo;
typedef seqan::Graph<seqan::Directed<TCargo> > TDgraph;
typedef seqan::VertexDescriptor<TDgraph>::Type TDVertexDescriptor;
typedef seqan::EdgeDescriptor<TDgraph>::Type TDEdgeDescriptor;
typedef seqan::Graph<seqan::Undirected<TCargo> > TUgraph;
typedef seqan::VertexDescriptor<TUgraph>::Type TUVertexDescriptor;
typedef seqan::EdgeDescriptor<TUgraph>::Type TUEdgeDescriptor;

typedef seqan::Map<TPosition, TDVertexDescriptor> TMapDGraph;
typedef seqan::String<TMapDGraph > TMapDGraphStr;
typedef seqan::Map<TPosition, TUVertexDescriptor> TMapUGraph;
typedef seqan::String<TMapUGraph > TMapUGraphStr;

typedef seqan::RnaString TSequence;
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;      // align type

struct vectGraphElement {
	std::vector<TUVertexDescriptor > uVertexVect;
	TUgraph interGraph; // this graph represents all the computed interaction edges
	std::vector<TDVertexDescriptor > dVertexVect;
	TDgraph interGraphUpdated;
};

template <typename TString, typename TPosition>
struct fixedStructElement {
	TString method; // place the method and parameters used to compute the structure
//	seqan::String<unsigned> structure;
	seqan::String<TPosition> seqPos;
	seqan::String<TPosition> interPos;
};
template <typename TString, typename TBioval>
struct bioValStructElement {
	TString method; // place the method and parameters used to compute the structure
//	seqan::String<TBioval> val;
	seqan::String<TBioval> val;
};

template <typename TSequence, typename TString, typename TPosition, typename TBioval, typename TMapLine>
struct RnaStructSeq {
	TSequence seq;  // from fasta input
	TString qual;  // from fasta input
	TString id;  // from fasta input
	TString info; // raw info from bpseq or ebpseq input
	seqan::String<fixedStructElement<TString, TPosition> > structPairMate; //TODO use this field to collect all the structural information of this sequence
	seqan::String<bioValStructElement<TString, TBioval> > structBioVal;
	vectGraphElement bpProb;
//	seqan::String<std::map<TPosition, TValue> > *bpProb;
//	TMapLine *bpProb;
};

typedef RnaStructSeq<TSequence, TString, TPosition, TBioval, TMapLine> TRnaStruct;
// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The first (and only) argument of the program is stored here.
    seqan::CharString rna_file;

    AppOptions() :
        verbosity(1)
    {}
};

template <typename TAlphabet>
struct StockholmRecord {
	// header (GF tags -> tag value)
	std::map<std ::string, std::string > header;
	// seqence names -> sequence maps
	std::map<std::string, std::string > seqences;
	// per column (GC) -> annotation string
	std::map<std::string, std::string > seqence_information;

	//TODO: Maybe support per-sequence (GS) and per-residue (GR) annotation
};

// ==========================================================================
// Functions
// ==========================================================================

void read_Stockholm_file(char * file, StockholmRecord<seqan::Rna>& record) {
	// read Stockholm format 
	std::ifstream inStream(file, std::ios::in);

	// go to block with metadata (#=GF blocks, no other right now)
	std::string line;

	do { // skip newlines to the first #= block
		std::getline(inStream, line);
	} while (line.substr(0, 2) != "#=");

	std::stringstream test;

	int line_num = 2;

	// read the file and save the lines depending on their flags
	do {
		++line_num;
		std::string tag, feature, value;
		std::istringstream iss(line);
		
		// skip if line is a newline
		if (line.find_first_not_of("\t\r\n ") == std::string::npos){
			continue;
		}

		// check if the line contains metadata
		if (line[0] == '#'){
			iss >> tag; 						// get the tag
			tag = tag.substr(2, tag.size());	// remove #= from tag
			iss >> feature; 					// get the feature description
			std::getline(iss, value);		 	// the rest of the line is the value
			// strip space before the line and newline at the end
			value.erase(0, value.find_first_not_of(" \t"));
			value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
			value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());

			// store tag in the respective map
			if (tag == "GF"){
				if (record.header.find(feature) == record.header.end())
					record.header[feature] = value;
				// append the tag contents if they were spread over multiple lines
				else
					record.header[feature].append(" " + value);
			}

			if (tag == "GC"){
				record.seqence_information[feature] = value;
			}
		}
		// we have a sequence record
		else {
			std::string name, sequence;
			iss >> name >> sequence;
			record.seqences[name] = sequence;

		}
	}while (std::getline(inStream, line));

	/*for (auto elem : record.header)
	{
		std::cout << elem.first << " " << elem.second << "\n";
	}*/

//	for (auto elem : record.seqences)
//	{
//		std::cout << elem.first << " " << elem.second << "\n";
//	}

	/*for (auto elem : record.seqence_information)
	{
		std::cout << elem.first << " " << elem.second << "\n";
	}*/
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("RNAMotif");
    // Set short description, version, and date.
    setShortDescription(parser, "RNA motif generator");
    setVersion(parser, "0.1");
    setDate(parser, __DATE__);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIMULTIPLE ALIGNMENT\\fP\"");
    addDescription(parser, "Generate a searchable RNA motif from a multiple structural RNA alignment.");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "INPUT FILE"));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBRNAMotif\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;
    seqan::getArgumentValue(options.rna_file, parser, 0);

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    AppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
	
    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cout << "RNA motif generator\n"
              << "===============\n\n";

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << '\n'
                  << "RNA      \t" << options.rna_file << "\n\n";
    }

    // read the stockholm alignment
	StockholmRecord<seqan::Rna> record;
	read_Stockholm_file(seqan::toCString(options.rna_file), record);

	std::cout << "Structure prediction\n";

	std::vector<TRnaStruct> RNASeqs;


	for (auto elem : record.seqences)
	{
		size_t seq_size = elem.second.length();

		// create RNA data structure
		TRnaStruct rna;
		rna.id = seqan::CharString(elem.first);		// sequence name/id
		rna.seq = seqan::RnaString(elem.second);	// sequence

		// save sequence structure
		seqan::appendValue(rna.structPairMate, fixedStructElement<TString, TPosition>());
		//rna.structPairMate[0].seqPos


		// add new vertex for each base and save the vertex references
		rna.bpProb.uVertexVect.resize(seq_size);
		for (unsigned i=0; i < seq_size; ++i){
			rna.bpProb.uVertexVect[i] = seqan::addVertex(rna.bpProb.interGraph);
		}


		// create RNAlib data structures
		const char * seq = elem.second.c_str();
		char  *structure = (char*)vrna_alloc(sizeof(char) * (strlen(seq) + 1));
		vrna_fold_compound_t *vc = vrna_fold_compound(seq, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);

		// predict secondary structure (will create base pair probs in compound)
		double gibbs = (double)vrna_pf(vc, structure);
		printf("%s %zu\n%s (%6.2f)\n", seq, seq_size, structure, gibbs);

		// extract base pair probabilities from fold compound
		vrna_plist_t *pl1, *ptr;
		pl1 = vrna_plist_from_probs(vc, 0.05);

		// get size of probability list (a full list, if not filtered with a threshold,
		// contains the upper half of the n x n probability matrix (size (n x n)/2 - n)
		unsigned size;
		for(size = 0, ptr = pl1; ptr->i; size++, ptr++);

		// store base pair prob. in the interaction graph
		for(unsigned i=0; i<size;++i)
		{
			//std::cout << i << "_"<< pl1[i].i <<":"<< pl1[i].j <<"|"<< pl1[i].p <<"|"<< pl1[i].type << "\n";
			seqan::addEdge(rna.bpProb.interGraph, rna.bpProb.uVertexVect[pl1[i].i-1], rna.bpProb.uVertexVect[pl1[i].j-1], pl1[i].p);
		}

		RNASeqs.push_back(rna);
	}

    return 0;
}
