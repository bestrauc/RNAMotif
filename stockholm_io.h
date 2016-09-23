// ==========================================================================
//                               Stockholm_io.h
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef APPS_RNAMOTIF_STOCKHOLM_IO_H_
#define APPS_RNAMOTIF_STOCKHOLM_IO_H_

#include "stockholm_file.h"
#include <seqan/stream.h>

namespace seqan{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

void read_Stockholm_file(char * file, std::vector<StockholmRecord<seqan::Rna> >& records) {
	// read Stockholm format
	std::ifstream inStream(file, std::ios::in);

	// go to block with metadata (#=GF blocks, no other right now)
	std::string line;


	while (std::getline(inStream, line)){
		StockholmRecord<seqan::Rna> record;

		// Skip the #Stockholm and newlines at the start
		bool new_record;
		while ( (new_record = static_cast<bool>(std::getline(inStream, line))) && line.substr(0, 2) != "#=");

		if (!new_record)
			break;

		do {
			std::string tag, feature, value;
			std::istringstream iss(line);

			// skip if line is an empty line or a newline (i.e. if there is no non-whitespace)
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
					record.seqence_information[feature].append(value);
				}
			}
			// we have a sequence record
			else {
				std::string name, sequence;
				iss >> name >> sequence;
				// if sequence is new - add to dictionary
				if (record.seqences.find(name) == record.seqences.end())
					record.sequence_names.push_back(name);

				// append to empty string if new, else append to existing string
				record.seqences[name].append(sequence);
			}

		} while (std::getline(inStream, line) && line.substr(0,2) != "//");

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

		records.push_back(record);
	}

	std::cout << records.size() << " Stockholm records read\n";
}

template <typename TRecord, typename TContext, typename TForwardIter>
void readRecord(TRecord& record, TContext & context, TForwardIter & iter, Stockholm const & /*tag*/)
{
	seqan::CharString & buffer = context.buffer;

	// go to the start: "# STOCKHOLM 1.0"
	//skipUntil(iter, EqualsChar<'#'>());
	readLine(buffer, iter);
	std::cout << "TEST: " << buffer << "\n";

	readLine(buffer, iter);
	std::cout << "TEST: " << buffer << "\n";

	readLine(buffer, iter);
	std::cout << "TEST: " << buffer << "\n";

	readLine(buffer, iter);
	std::cout << "TEST: " << buffer << "\n";

	std::cout << record.header.size() << "\n";

}

}

#endif  // #ifndef APPS_RNAMOTIF_STOCKHOLM_IO_H_
