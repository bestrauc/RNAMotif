// ==========================================================================
//                                  motif.h
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

#ifndef APPS_RNAMOTIF_MOTIF_H_
#define APPS_RNAMOTIF_MOTIF_H_

#include "motif_structures.h"
#include <stack>

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

bool isOpen(char c){
	return (c == '(') || (c == '[') || (c == '{') || (c == '<');
}

const std::string parentheses("{}[]()<>");

// helper function to generate interaction vectors via bracket notation
void bracketToInteractions(const char* structure, TInteractionPairs& interaction){
	// create stacks for different types of bracket matches
	std::unordered_map<char, std::stack<int>* > match_stacks;

	// curly braces
	std::stack<int> round_stack;
	std::stack<int> curly_stack;
	std::stack<int> angle_stack;
	std::stack<int> edge_stack;
	match_stacks['('] = &round_stack;
	match_stacks[')'] = &round_stack;
	match_stacks['{'] = &curly_stack;
	match_stacks['}'] = &curly_stack;
	match_stacks['['] = &angle_stack;
	match_stacks[']'] = &angle_stack;
	match_stacks['<'] = &edge_stack;
	match_stacks['>'] = &edge_stack;

	unsigned n = strlen(structure);
	interaction = TInteractionPairs(n, std::make_pair(ROUND, -1));

	for (unsigned i = 0; i < n; ++i){
		char c = structure[i];

		// if we see an opening bracket (== not in match_table), put on corresponding stack
		if (isOpen(c))
			match_stacks[c]->push(i);
		// skip non-brackets
		else if (parentheses.find(c) == std::string::npos)
			continue;
		// else we have a closing bracket and record a match
		else if (!match_stacks[c]->empty()){
			int left = match_stacks[c]->top();
			BracketType btype = bracket_to_type[c];
			interaction[left] = std::make_pair(btype, i);
			interaction[i]    = std::make_pair(btype, left);
			match_stacks[c]->pop();
		}
		// we had a closing bracket, but no corresponding opening bracket (error)
		else
			std::cerr << "Missing opening bracket for " << c << " at position " << i << "\n";
	}
}

RNAProfileString addRNAProfile(StructureElement &structureElement, unsigned start, unsigned end, TAlign &align){
	typedef seqan::Row<TAlign>::Type TRow;

	RNAProfileString profileString;
	seqan::resize(profileString, end-start+1);
	StructureStatistics stats;

	// min and max length initialized with their most extreme possible values
	stats.min_length = end - start + 1;
	stats.max_length = 0;
	stats.mean_length = 0;

	// store the profile of the alignment in [start,end]
	for (unsigned row=0; row < length(rows(align)); ++row){
		TRow & i_row = seqan::row(align,row);
		unsigned source_start = seqan::toSourcePosition(i_row, start);
		unsigned source_end = seqan::toSourcePosition(i_row, end);
		unsigned seqLength = source_end -  source_start+1;

		// check if we only had gaps in the region (end-start+1 doesn't work then)
		if (source_start == source_end && seqan::isGap(i_row, end))
			seqLength = 0;

		// set the statistics (min, max, average lengths)
		stats.mean_length += seqLength;

		if (seqLength < stats.min_length){
			stats.min_length = seqLength;
		}

		if (seqLength > stats.max_length)
			stats.max_length = seqLength;

		// create profile of bases in this column
		for (unsigned i=0; i < seqan::length(profileString); ++i){
			int index = i+start;
			profileString[i].count[seqan::ordValue(seqan::row(align, row)[index])] += 1;
		}
	}

	stats.mean_length = stats.mean_length / length(rows(align));

	if (structureElement.type == StructureType::HAIRPIN && stats.min_length < 3)
		stats.min_length = 3;

	structureElement.components.push_back(profileString);
	structureElement.statistics.push_back(stats);

	return profileString;
}

TSequenceRegions findStemLoops(TInteractionPairs &consensus){
	TSequenceRegions stemLoops;

	std::vector<std::stack<int> > pairStacks;
	std::vector<std::pair<int, int> > lastHairpins;
	for (size_t i = 0; i < MAX; ++i){
		pairStacks.push_back(std::stack<int>());
		lastHairpins.push_back(std::make_pair(0,0));
	}

	// locate the stem loops
	for (size_t i = 0; i < consensus.size(); ++i){
		std::pair<BracketType, int> b = consensus[i];
		BracketType btype = b.first;
		int bracket = b.second;

		// skip unmatched regions
		if (bracket == -1)
			continue;

		// opening bracket, save previous stem-loop if there was one
		if (bracket > i){
			// if we found a open/close match before, save that and reset
			if (lastHairpins[btype].second > 0){
				stemLoops.push_back(std::make_pair(btype, lastHairpins[btype]));

				// clear stack and reset last hairpin found
				std::stack<int>().swap(pairStacks[btype]);
				lastHairpins[btype] = std::pair<int, int>();
			}

			// save the opening bracket on the stack
			pairStacks[btype].push(i);
		}
		// closing bracket
		else {
			if (!pairStacks[btype].empty() && bracket == pairStacks[btype].top()){
				lastHairpins[btype].first = bracket;
				lastHairpins[btype].second = i;
				pairStacks[btype].pop();
			}
		}
	}

	// save last stem-loop if there is one
	for (unsigned i=0; i< lastHairpins.size(); ++i)
		if (lastHairpins[i].second > 0)
			stemLoops.push_back(std::make_pair((BracketType)i, lastHairpins[i]));

	return stemLoops;
}

// * a run of opening brackets is the left side of a stem
// * any unpaired bases in between are interior loops
// 		* if no corresponding unpaired bases: bulge
// * innermost unpaired bases are the hairpin
void partitionStemLoop(Motif &motif, BracketType btype, std::pair<int, int > stemLoopRegion){
	//TInteractionPairs &consensus = motif.consensusStructure;
	TInteractionPairs consensus(motif.consensusStructure);

	// consider only the brackets of type btype (TODO: ignore more elegantly without replacing anything?)
	for (size_t i=0; i < consensus.size(); ++i){
		if (consensus[i].first != btype)
			consensus[i] = std::make_pair(MAX, -1);
	}

	TStructure stemStructure;

	size_t i = stemLoopRegion.first;
	do {
		int pos = i;
		int right = consensus[pos].second;

		// count the row of opening brackets
		if (right > pos){					// an opening bracket
			while (right > pos){			// while it's a series of opening brackets
				// check if the corresponding closing bracket follows
				// or if there is a bulge on the right side
				if ((consensus[pos+1].second > pos+1) && consensus[right-1].second == -1){
					++pos;

					// add the stem up to the right bulge
					StructureElement stem;

					DEBUG_MSG("Stem: [" << i << "," << pos-1 << " " << pos-i << "] ; [" << consensus[pos-1].second << "," << consensus[i].second << " " << consensus[i].second - consensus[pos-1].second+1 << "]");

					stem.type = STEM;
					RNAProfileString leftProfile  = addRNAProfile(stem, i, pos-1, motif.seedAlignment);
					RNAProfileString rightProfile = addRNAProfile(stem, consensus[pos-1].second, consensus[i].second, motif.seedAlignment);

					// find the extension of the right bulge and add it
					// there is one run of unpaired bases from right+1
					int unpaired = right-1;
					while (consensus[unpaired].second ==-1) --unpaired;

					DEBUG_MSG("Right bulge in [" << unpaired+1 << "," << right-1 << " " << right - unpaired-1 << "]");

					StructureElement bulge;

					bulge.type = LOOP;
					RNAProfileString bulgeProfile = addRNAProfile(bulge, unpaired+1, right-1, motif.seedAlignment);
					stemStructure.push_back(bulge);

					i = pos;
				}

				++pos;
				right = consensus[pos].second;
			}

			// add the uninterrupted stem we found so far
			StructureElement stem;

			DEBUG_MSG("Stem: [" << i << "," << pos-1 << " " << pos-i << "] ; [" << consensus[pos-1].second << "," << consensus[i].second << " " << consensus[i].second - consensus[pos-1].second+1 << "]");

			stem.type = STEM;
			RNAProfileString leftProfile  = addRNAProfile(stem, i, pos-1, motif.seedAlignment);
			RNAProfileString rightProfile = addRNAProfile(stem, consensus[pos-1].second, consensus[i].second, motif.seedAlignment);

			stemStructure.push_back(stem);
		}
		// if unpaired, count the length of the loop and check for a bulge
		else if (right == -1){
			StructureElement structure;

			// get right border bracket of other half of loop (->(...(..)...)<=)
			int run = pos;
			int rb = consensus[pos-1].second;
			while (consensus[run].second == -1) ++run;

			// get partner of end bracket ((...->(..)<=...))
			int lb = consensus[run].second;

			// Left bulge
			if (rb - lb == 1){
				DEBUG_MSG("Left bulge in [" << pos << "," << run-1 << " " << run-pos << "]");

				structure.type = LOOP;
				RNAProfileString bulgeProfile = addRNAProfile(structure, pos, run-1, motif.seedAlignment);
			}
			// Hairpin (stop the outer loop here since all structures found)
			else if (rb == run){
				DEBUG_MSG("Hairpin in [" << pos << "," << run-1 << " " << run-pos << "]");

				structure.type = HAIRPIN;
				RNAProfileString hairpinProfile = addRNAProfile(structure, pos, run-1, motif.seedAlignment);
				break;
			}
			// Interior loop with left and right side
			else{
				DEBUG_MSG("Left loop: [" << pos << "," << run-1 << "]" << " " << run-pos << " ; " << "Right loop: [" << lb+1 << "," << rb-1 << " " << rb-1-lb << "]");

				structure.type = LOOP;
				RNAProfileString leftProfile  = addRNAProfile(structure, pos, run-1, motif.seedAlignment);
				RNAProfileString rightProfile = addRNAProfile(structure, lb+1, rb-1, motif.seedAlignment);
			}

			stemStructure.push_back(structure);

			pos = run;
			right = consensus[pos].second;
		}
		else
			pos = pos + 1;

		i = pos;

	} while (i <= stemLoopRegion.second);

	motif.profile.push_back(stemStructure);
}

// Example: (((((((.......((((((((..(((((..(((((.....((((((....((((...))))))))))...((((....((.....))....))))))))).)))))....))).))))))))))))........((((((.....)))))).................
// take a structure table and determine the structural elements (stem, bulge, internal loop, hairpin)
void structurePartition(Motif &motif){
	TSequenceRegions stemLoops = findStemLoops(motif.consensusStructure);

	int pos = -8;
	for (auto pair : stemLoops){
		BracketType btype = pair.first;
		TRegion region = pair.second;
		//std::cout << pair.first << " " << region.first << " " << region.second << "\n";

		// visualize stem loops for debugging
		//std::cout << std::string(region.first-pos, ' ');
		//std::cout << std::string(region.second-region.first+1, '+');
		//pos += (region.first-pos) + (region.second-region.first+1);
	}

	//std::cout << "\n";

	// after locating stem loops, separate structural elements
	for (auto pair : stemLoops){
		BracketType btype = pair.first;
		TRegion region = pair.second;

		// find structural elements
		partitionStemLoop(motif, btype, region);
	}

	return;
}

#endif  // #ifndef APPS_RNAMOTIF_MOTIF_H_
