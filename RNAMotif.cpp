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

// SeqAn headers
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/graph_types.h>

// C++ headers
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// App headers
#include "RNAlib_utils.h"
#include "IPknot_utils.h"
#include "motif.h"

// ==========================================================================
// Classes
// ==========================================================================

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
    bool constrain;
    bool pseudoknot;

    // The first (and only) argument of the program is stored here.
    seqan::CharString rna_file;

    AppOptions() :
        verbosity(1),
		constrain(0)
    {}
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
		
		// skip if line is an empty line, a newline or the alignment end \\ line
		if (line.find_first_not_of("\t\r\n/ ") == std::string::npos){
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
			record.seqNames.push_back(name);
			record.seqs.push_back(sequence);
		}
	}while (std::getline(inStream, line));

	seqan::resize(seqan::rows(record.alignment), record.seqences.size());

	int i = 0;
	for (auto elem : record.seqences)
	{
		//TODO: remove gaps from sequences. Those will be conserved in the align-object.

		// erase all gaps from the string and insert it into the alignment
		std::string tmp = elem.second;
		tmp.erase(std::remove(tmp.begin(), tmp.end(), '-'), tmp.end());
		seqan::assignSource(seqan::row(record.alignment, i), seqan::RnaString(tmp));
	    TRow & row = seqan::row(record.alignment, i++);

		// find all gap positions in the sequence and insert into alignment
	    int offset = 0;
		size_t pos = elem.second.find("-", 0);
		while (pos != std::string::npos){
			// find how long the gap is
			int len = 1;
			while (elem.second[pos+len] == '-') ++len;
			//std::cout << pos << " " << len << " " << offset << "\n";
		    seqan::insertGaps(row, pos, len);
		    pos = elem.second.find("-", pos+len);
		    offset += len;
		}

	}
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
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fISEED ALIGNMENT\\fP>");
    addDescription(parser, "Generate a searchable RNA motif from a seed alignment.");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "INPUT FILE"));

    addOption(parser, seqan::ArgParseOption("ps", "pseudoknot", "Predict structure with IPknot to include pseuoknots."));
    addOption(parser, seqan::ArgParseOption("co", "constrain", "Constrain individual structures with the seed consensus structure."));
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

    options.constrain  = isSet(parser, "constrain");
    options.pseudoknot = isSet(parser, "pseudoknot");

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
                  << "CONSTRAINT\t" << options.constrain << '\n'
				  << "PSEUDOKNOTS\t" << options.pseudoknot << '\n'
                  << "RNA      \t" << options.rna_file << "\n\n";
    }

    // read the stockholm alignment
	StockholmRecord<seqan::Rna> record;
	read_Stockholm_file(seqan::toCString(options.rna_file), record);

	std::cout << record.alignment << "\n";

	// convert Rfam WUSS structure to normal brackets to get a constraint
	char *constraint_bracket = NULL;
	if (options.constrain){
		constraint_bracket = new char[record.seqence_information["SS_cons"].length() + 1];
		WUSStoPseudoBracket(record.seqence_information["SS_cons"], constraint_bracket);
	}

   	int i = 0;

   	Motif rna_motif;
   	rna_motif.seedAlignment = record.alignment;
   	rna_motif.interactionGraphs.resize(record.seqences.size());
   	rna_motif.interactionPairs.resize(record.seqences.size());

	// build interaction graphs for each sequence
	for (auto elem : record.seqences)
	{
		createInteractions(rna_motif.interactionGraphs[i],
						   rna_motif.interactionPairs[i],
						   elem.second, constraint_bracket);
		i++;
	}

	// create structure for the whole multiple alignment
	std::cout << "Rfam:   " << record.seqence_information["SS_cons"] << std::endl;
	//getConsensusStructure((const char**)seqs, rna_motif.consensusStructure, bracket, IPknotFold());
	getConsensusStructure(record, rna_motif.consensusStructure, constraint_bracket, IPknotFold());

	structurePartition(rna_motif);

	std::cout << std::endl;

	free(constraint_bracket);

    return 0;
}
