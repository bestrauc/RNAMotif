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

template <typename TRecord, typename TContext, typename TForwardIter>
void readRecord(TRecord& record, TContext & context, TForwardIter & iter, Stockholm const & /*tag*/)
{
	seqan::CharString & buffer = context.buffer;

	// go to the start: "# STOCKHOLM 1.0"
	skipUntil(iter, EqualsChar<'#'>());
	readLine(buffer, iter);
	if (buffer != "# STOCKHOLM 1.0"){
		ParseError("'# STOCKHOLM 1.0' header not found. Valid Stockholm format?");
	}

	skipUntil(iter, NotFunctor<IsWhitespace>());

	do{
		// both if-branches end with readUntil(..,IsNewLine) -> read one line
		if (value(iter) == '#'){
			std::string tag, feature, val;

			// skip '#='
			skipUntil(iter, IsAlpha());
			readUntil(tag, iter, IsWhitespace());
			skipUntil(iter, NotFunctor<IsWhitespace>());
			readUntil(feature, iter, IsWhitespace());
			skipUntil(iter, NotFunctor<IsWhitespace>());
			readUntil(val, iter, IsNewline());

			//std::cout << tag << "," << feature << "," << val << "\n";

			// store tag in the respective map
			if (tag == "GF"){
				if (record.header.find(feature) == record.header.end())
					record.header[feature] = val;
				// append the tag contents if they were spread over multiple lines
				else
					record.header[feature].append(" " + val);
			}

			if (tag == "GC"){
				record.seqence_information[feature].append(val);
			}
		}
		else{
			std::string name, sequence;

			readUntil(name, iter, IsWhitespace());
			skipUntil(iter, NotFunctor<IsWhitespace>());
			readUntil(sequence, iter, IsNewline());

			//std::cout << name << "," << sequence << "\n";

			// if sequence is new - add to dictionary
			if (record.seqences.find(name) == record.seqences.end())
				record.sequence_names.push_back(name);

			// append to empty string if new, else append to existing string
			record.seqences[name].append(sequence);
		}

		// skip possibly empty lines
		skipUntil(iter, NotFunctor<IsWhitespace>());

	// end the record at the '//' line
	} while (value(iter) != '/');
	//TODO: end of record is actually '//', but we just check for '/' ?

	// skip to start of next record (if possible, else we get to the end of the file)
	skipUntil(iter, IsWhitespace());
	skipUntil(iter, NotFunctor<IsWhitespace>());

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

}

#endif  // #ifndef APPS_RNAMOTIF_STOCKHOLM_IO_H_
