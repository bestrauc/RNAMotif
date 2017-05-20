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
#include "motif_search.h"

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

template <typename TAlphabet>
double characterEntropy(seqan::ProfileChar<TAlphabet> &profChar){
	double h = 0;
	int size = seqan::ValueSize<TAlphabet>::VALUE-1;

	double norm_sum = 0;
	for (int i=0; i < size; ++i){
		norm_sum += profChar.count[i];
	}

	for (int i=0; i < size; ++i){
		double freq = profChar.count[i] / norm_sum;

		if (freq == 0)
			continue;

		h += freq * std::log(freq);
	}

	return -h;
}

template <typename TAlphabet>
double profileEntropy(seqan::String<seqan::ProfileChar<TAlphabet> > &profileString){
	double h = 0;
	int n = seqan::length(profileString);

	for (int i = 0; i < n; ++i){
		double ch = characterEntropy(profileString[i]);
		h += ch;

		//std::cout << ch << " ";
	}

	//std::cout << "\n";

	return h/n;
}

double loopEntropy(std::vector<TLoopProfileString> &profileStringComponents){
	//std::cout << profileStringComponents[0] << " " << profileStringComponents[1] << "\n";

	int n = seqan::length(profileStringComponents);

	double H1 = profileEntropy(profileStringComponents[0]);
	double H2 = n == 2 ? profileEntropy(profileStringComponents[1]) : 0;

	return (H1+H2)/n;
}

double stemEntropy(TStemProfileString &stemProfile){
	//std::cout << profileStringComponents[0] << " " << profileStringComponents[1] << "\n";

	return profileEntropy(stemProfile);
}

// helper function to generate interaction vectors via bracket notation
void bracketToInteractions(const char* structure, TConsensusStructure& interaction){
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
	interaction = TConsensusStructure(n, std::make_pair(ROUND, -1));

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

TLoopProfileString addProfile(StructureElement &structureElement, unsigned start, unsigned end, TAlign &align, std::set<int> &excludeSet, bool reverse = false){
	typedef seqan::Row<TAlign>::Type TRow;
	typedef typename seqan::Value<TLoopProfileString>::Type TProfileChar;

	TLoopProfileString profileString;
	std::vector<std::map<int, int> > gapString;
	int n = end-start+1;
	seqan::resize(profileString, n);
	seqan::resize(gapString, n);
	StructureStatistics stats;

	// min and max length initialized with their most extreme possible values
	stats.min_length = end - start + 1;
	stats.max_length = 0;
	stats.mean_length = 0;

	// store the profile of the alignment in [start,end]
	for (unsigned row=0; row < length(rows(align)); ++row){
		//if (row != 45)
		//	continue;

		if (excludeSet.find(row) != excludeSet.end()){
			continue;
		}

		TRow & i_row = seqan::row(align,row);
		//std::cout << seqan::infix(i_row, start, end+1) << "\n";

		unsigned source_start = seqan::toSourcePosition(i_row, start);
		unsigned source_end = seqan::toSourcePosition(i_row, end);
		// get length of sequence without gaps (don't add +1 if the last position is a gap)
		unsigned seqLength = source_end -  source_start + (1-seqan::isGap(i_row, end));

		// check if we only had gaps in the region (end-start+1 doesn't work then)
		// remove the whole alignment row in case the hairpin is just a gap
		if (source_start == source_end && seqan::isGap(i_row, end)){
			if (structureElement.type == StructureType::HAIRPIN){
				std::cout << "SKIPPING\n";
				excludeSet.insert(row);
				continue;
			}

			seqLength = 0;
		}

		//if (seqLength == 0)
		//	continue;

		// set the statistics (min, max, average lengths)
		stats.mean_length += seqLength;

		if (seqLength < stats.min_length){
			stats.min_length = seqLength;
		}

		if (seqLength > stats.max_length)
			stats.max_length = seqLength;

		//std::cout << seqan::infix(i_row, start, end+1) << " " << start << " " << end << " " << source_start << " " << source_end << "\n";

		unsigned n = seqan::length(profileString);

		int gap_run = 0;
		int run_start = -1;

		// create profile of bases in this column
		for (unsigned i=0; i < n; ++i){
			int index = i+start;

			// store length of gaps at each position
			if (seqan::isGap(i_row, index)){
				//std::cout << index << " ";
				++gap_run;
				run_start = run_start == -1 ? i : run_start;
			}
			else if (gap_run > 0){
				//std::cout << "Run: " << run_start << " " << gap_run << "\n";
				gapString[run_start][gap_run]++;
				gap_run = 0;
				run_start = -1;
			}

			//std::cout << "(" << seqan::ordValue(seqan::row(align, row)[index]) << "," << seqan::row(align, row)[index] << ")" << " ";
			unsigned ord_val = seqan::ordValue(seqan::row(align, row)[index]);

			// save gaps as last character in alphabet
			if (ord_val > seqan::ValueSize<TProfileChar>::VALUE){
				continue;
				ord_val = AlphabetSize-1;
				//std::cout << "Gappé! " << ord_val << "\n";
			}

			if (reverse)
				profileString[n-i-1].count[ord_val] += 1;
			else
				profileString[i].count[ord_val] += 1;
		}

		if (gap_run > 0){
			//std::cout << "Run: " << run_start << " " << gap_run << "\n";
			gapString[run_start][gap_run]++;
		}
		//std::cout << std::endl;
	}

	//std::cout << "\n";

	stats.mean_length = stats.mean_length / length(rows(align));

	if (structureElement.type == StructureType::HAIRPIN && stats.min_length < 3)
		stats.min_length = 3;

	structureElement.loopComponents.push_back(profileString);
	structureElement.statistics.push_back(stats);
	structureElement.gap_lengths = gapString;

	return profileString;
}

TStemProfileString addProfile(StructureElement &structureElement, unsigned start1, unsigned end1, unsigned start2, unsigned end2, TAlign &align, std::set<int> &excludeSet){

	typedef typename seqan::Value<TStemProfileString>::Type TProfileChar;

	TStemProfileString profileString;
	std::vector<std::map<int, int> > gapString;
	int n = end1-start1+1;
	seqan::resize(profileString, n);
	seqan::resize(gapString, n);

	// store the profile of the alignment in [start,end]
	for (unsigned row=0; row < length(rows(align)); ++row){
		//if (row != 45)
		//	continue;

		// skip row if we excluded it because of gaps that are too large
		if (excludeSet.find(row) != excludeSet.end()){
			continue;
		}

		// create binucleotide profile of bases in the two columns
		TRow & i_row = seqan::row(align,row);

		//std::cout << seqan::infix(i_row,  << "\n";

		int gap_run = 0;
		int run_start = -1;

		for (unsigned i=0; i < n; ++i){
			int l_index = i+start1; // left stem
			int r_index = end2-i;	// right stem

			// store length of gaps at each position
			if (seqan::isGap(i_row, l_index) && seqan::isGap(i_row, r_index)){
				++gap_run;
				run_start = run_start == -1 ? i : run_start;
			}
			else if (gap_run > 0){
				gapString[run_start][gap_run]++;
				gap_run = 0;
				run_start = -1;
			}

			unsigned l_ord_val = seqan::ordValue(i_row[l_index]);
			unsigned r_ord_val = seqan::ordValue(i_row[r_index]);

			//std::cout << l_ord_val << " " << r_ord_val << "\n";

			if (l_ord_val == AlphabetSize-1 || r_ord_val == AlphabetSize-1 || l_ord_val == 4 || r_ord_val == 4){
				throw(std::runtime_error("Whaaaaa?"));
			}

			// if neither character is a gap
			if (l_ord_val < AlphabetSize && r_ord_val < AlphabetSize){
				unsigned pair_val = (l_ord_val*AlphabetSize) + r_ord_val;
				assert((pair_val / AlphabetSize) == l_ord_val);
				assert((pair_val % AlphabetSize) == r_ord_val);

				//assert(false);

				profileString[n-i-1].count[pair_val] += 1;
			}
			/*
			// if both characters are gaps
			else if (l_ord_val > AlphabetSize && r_ord_val > AlphabetSize){
				l_ord_val = AlphabetSize-1;
				r_ord_val = AlphabetSize-1;
				unsigned pair_val = (l_ord_val*AlphabetSize) + r_ord_val;
				//std::cout << "Found a stem gap! " << pair_val << "\n";
				profileString[n-i-1].count[pair_val] += 1;
			}
			*/
			// ignore cases where only one of the sides is a gap (can't deal with those)
		}

		if (gap_run > 0){
			gapString[run_start][gap_run]++;
		}
	}

	structureElement.stemProfile = profileString;
	structureElement.gap_lengths = gapString;

	return profileString;
}

std::vector<TLoopProfileString> getExternal(TConsensusStructure &consensus, TAlign &align){
	TLoopProfileString profile;
	typedef typename seqan::Value<TLoopProfileString>::Type TProfileChar;
	int aln_len = seqan::length(seqan::row(align,0));

	unsigned i=0;

	while (i < aln_len){
		std::pair<BracketType, int> b = consensus[i];
 		int partner = b.second;

		if (partner == -1){
			TProfileChar profChar;
			// store the profile of the alignment in [start,end]
			for (unsigned row=0; row < length(rows(align)); ++row){
				// create profile of bases in this column
				unsigned ord_val = seqan::ordValue(seqan::row(align, row)[i]);

				// save gaps as last character in alphabet
				if (ord_val > seqan::ValueSize<TProfileChar>::VALUE)
					ord_val = AlphabetSize-1;

				profChar.count[ord_val] += 1;
			}

			seqan::appendValue(profile, profChar);
		}
		else{
			i = partner;
		}

		++i;
	}

	std::vector<TLoopProfileString> ret = {profile};

	return ret;
}

TStemLoopProfile findStemLoops(TConsensusStructure const &consensus){
	TStemLoopProfile stemLoops;

	std::vector<std::stack<int> > pairStacks;
	std::vector<std::pair<int, int> > lastHairpins;

	std::vector<TInteractions> filteredInteractions;

	for (size_t i = 0; i < MAX; ++i){
		pairStacks.push_back(std::stack<int>());
		lastHairpins.push_back(std::make_pair(0,0));

		// create all filtered consensus structures
		filteredInteractions.push_back(TInteractions(consensus.size()));

		for (size_t k=0; k < consensus.size(); ++k){
			if (consensus[k].first != (BracketType)i || consensus[k].second == -1){
				filteredInteractions[i][k] = -1;
				//std::cout << ".";
			}
			else{
				filteredInteractions[i][k] = consensus[k].second;
				//std::cout << ((k < filteredInteractions[i][k]) ? "(" : ")");
			}
		}

		//std::cout << "\n";
	}

	// locate the stem loops
	for (size_t i = 0; i < consensus.size(); ++i){
		std::pair<BracketType, int> b = consensus[i];
		BracketType btype = b.first;
		int partner = b.second;

		// skip unmatched regions
		if (partner == -1)
			continue;

		// opening bracket, save previous stem-loop if there was one
		if (partner > i){
			// if we found a open/close match before, save that and reset
			if (lastHairpins[btype].second > 0){
				// save filtered consensus structure in the bracket type
				stemLoops.push_back(TStructure(btype, lastHairpins[btype],
											  filteredInteractions[btype]));

				// clear stack and reset last hairpin found
				std::stack<int>().swap(pairStacks[btype]);
				lastHairpins[btype] = std::pair<int, int>();
			}

			// save the opening bracket on the stack
			pairStacks[btype].push(i);
		}
		// closing bracket
		else {
			if (pairStacks[btype].empty())
				continue;

			lastHairpins[btype].first = partner;
			lastHairpins[btype].second = i;
			pairStacks[btype].pop();

			if (pairStacks[btype].empty()){
				stemLoops.push_back(TStructure(btype, lastHairpins[btype],
											  filteredInteractions[btype]));

				// clear stack and reset last hairpin found
				lastHairpins[btype] = std::pair<int, int>();
			};
		}
	}
	return stemLoops;
}

// * a run of opening brackets is the left side of a stem
// * any unpaired bases in between are interior loops
// 		* if no corresponding unpaired bases: bulge
// * innermost unpaired bases are the hairpin
void partitionStemLoop(TAlign &seedAlignment, TStructure &stemStructure){
	TInteractions &consensus = stemStructure.interactions;
	//TInteractionPairs consensus(motif.consensusStructure);

	{
		std::set<int> tmpSet;
		StructureElement tmp;
		tmp.type == HAIRPIN;
		addProfile(tmp, stemStructure.pos.first, stemStructure.pos.second, seedAlignment, tmpSet);
		std::cout << "VARIATION: " << tmp.statistics[0].min_length << " Max: " << tmp.statistics[0].max_length << "\n";
	}

	std::set<int> excludeSet;
	size_t i = stemStructure.pos.first;
	do {
		int pos = i;
		int right = consensus[pos];

		// count the row of opening brackets
		if (right > pos){					// an opening bracket
			while (right > pos){			// while it's a series of opening brackets
				// check if the corresponding closing bracket follows
				// or if there is a bulge on the right side
				if ((consensus[pos+1] > pos+1) && consensus[right-1] == -1){
					++pos;

					// add the stem up to the right bulge
					StructureElement stem;

					DEBUG_MSG("Stem: [" << i << "," << pos-1 << " " << pos-i << "] ; [" << consensus[pos-1] << "," << consensus[i] << " " << consensus[i] - consensus[pos-1]+1 << "]");

					stem.type = STEM;
					stem.location = i;
					TLoopProfileString leftProfile  = addProfile(stem, i, pos-1, seedAlignment, excludeSet, true);
					TLoopProfileString rightProfile = addProfile(stem, consensus[pos-1], consensus[i], seedAlignment, excludeSet);
					TStemProfileString stemProfile  = addProfile(stem, i, pos-1, consensus[pos-1], consensus[i], seedAlignment, excludeSet);

					stemStructure.elements.push_back(stem);

					// find the extension of the right bulge and add it
					// there is one run of unpaired bases from right+1
					int unpaired = right-1;
					while (consensus[unpaired] ==-1) --unpaired;

					DEBUG_MSG("Right bulge in [" << unpaired+1 << "," << right-1 << " " << right - unpaired-1 << "]");

					StructureElement bulge;

					bulge.type = RBULGE;
					bulge.location = unpaired+1;
					TLoopProfileString bulgeProfile = addProfile(bulge, unpaired+1, right-1, seedAlignment, excludeSet);
					bulge.loopLeft = false;
					stemStructure.elements.push_back(bulge);

					i = pos;
				}

				++pos;
				right = consensus[pos];
			}

			// add the uninterrupted stem we found so far
			StructureElement stem;

			DEBUG_MSG("Stem: [" << i << "," << pos-1 << " " << pos-i << "] ; [" << consensus[pos-1] << "," << consensus[i] << " " << consensus[i] - consensus[pos-1]+1 << "]");

			stem.type = STEM;
			stem.location = i;
			TLoopProfileString leftProfile  = addProfile(stem, i, pos-1, seedAlignment, excludeSet, true);
			TLoopProfileString rightProfile = addProfile(stem, consensus[pos-1], consensus[i], seedAlignment, excludeSet);
			TStemProfileString stemProfile  = addProfile(stem, i, pos-1, consensus[pos-1], consensus[i], seedAlignment, excludeSet);

			stemStructure.elements.push_back(stem);
		}
		// if unpaired, count the length of the loop and check for a bulge
		else if (right == -1){
			StructureElement structure;
			StructureElement structure2;

			// get right border bracket of other half of loop (->(...(..)...)<=)
			int run = pos;
			int rb = consensus[pos-1];
			while (consensus[run] == -1) ++run;

			// get partner of end bracket ((...->(..)<=...))
			int lb = consensus[run];

			// Left bulge
			if (rb - lb == 1){
				DEBUG_MSG("Left bulge in [" << pos << "," << run-1 << " " << run-pos << "]");

				structure.type = LBULGE;
				TLoopProfileString bulgeProfile = addProfile(structure, pos, run-1, seedAlignment, excludeSet, true);
				structure.loopLeft = true;
			}
			// Hairpin (stop the outer loop here since all structures found)
			else if (rb == run){
				DEBUG_MSG("Hairpin in [" << pos << "," << run-1 << " " << run-pos << "]");

				structure.type = HAIRPIN;
				TLoopProfileString hairpinProfile = addProfile(structure, pos, run-1, seedAlignment, excludeSet);
				std::cout << "HAIRPIN: " << structure.statistics[0].min_length << " Max: " << structure.statistics[0].max_length << "\n";
				structure.loopLeft = false;
			}
			// Interior loop with left and right side
			else{
				DEBUG_MSG("Left loop: [" << pos << "," << run-1 << "]" << " " << run-pos << " ; " << "Right loop: [" << lb+1 << "," << rb-1 << " " << rb-1-lb << "]");

				structure.type = LOOP;
				TLoopProfileString leftProfile  = addProfile(structure, pos, run-1, seedAlignment, excludeSet, true);
				TStemProfileString loopProfile  = addProfile(structure, pos, run-1, lb+1, rb-1, seedAlignment, excludeSet);
				structure.loopLeft = true;
				TLoopProfileString rightProfile = addProfile(structure2, lb+1, rb-1, seedAlignment, excludeSet);
				TStemProfileString loopProfile2  = addProfile(structure2, pos, run-1, lb+1, rb-1, seedAlignment, excludeSet);
				structure2.loopLeft = false;

				//TLoopProfileString leftProfile  = addProfile(structure, pos, run-1, seedAlignment, true);
				//TLoopProfileString rightProfile = addProfile(structure, lb+1, rb-1, seedAlignment);
				//TStemProfileString loopProfile  = addProfile(structure, pos, run-1, lb+1, rb-1, seedAlignment);
			}

			structure.location = pos;
			stemStructure.elements.push_back(structure);

			if (structure.type == HAIRPIN)
				break;

			pos = run;
			right = consensus[pos];
		}
		else
			pos = pos + 1;

		i = pos;

	} while (i <= stemStructure.pos.second);
}

#endif  // #ifndef APPS_RNAMOTIF_MOTIF_H_
