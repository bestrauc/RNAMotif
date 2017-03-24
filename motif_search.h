// ==========================================================================
//                               motif_search.h
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

#ifndef APPS_RNAMOTIF_MOTIF_SEARCH_H_
#define APPS_RNAMOTIF_MOTIF_SEARCH_H_

#include "motif_structures.h"
#include "motif.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/* ------------------------------------------------------- */
/*!
 * @class ProfileCharIter
 *
 * @brief Abstract base class of the ProfileChar iterator.
 *
 *	Allows us to use different subclass iterators for different ProfileChar
 *	Alphabet types while searching. (Mainly for Nucleotides and BiNucleotides)
 */

class ProfileCharIter{
public:
	bool inner  = false;
	bool gapped = false;

	virtual int getNextChar() = 0;
	virtual bool atEnd() = 0;
	virtual ~ProfileCharIter() {};
};

/*!
 * @class ProfileCharIterImpl
 *
 * @brief Implement the ProfileCharIter for all types of ProfileChars.
 *
 *	Templated implementation for ProfileChar types of different alphabets.
 *	Return the next characters as for all as int for consistency.
 */

template <typename TProfileChar>
class ProfileCharIterImpl : public ProfileCharIter{
	typedef typename seqan::SourceValue<TProfileChar>::Type TProfileAlphabet;
	TProfileChar c;
	const int char_size = seqan::ValueSize<TProfileChar>::VALUE;

	int state = 0;
	std::vector<int> idx;
	int total;

public:
	ProfileCharIterImpl (TProfileChar c) : c(c), idx(char_size), total(seqan::totalCount(c)){
		// fill index vector with indies [0,1,..,N-1]
		std::iota(idx.begin(), idx.end(), 0);

		// sort index by comparing the entries in tmp to get order
		std::sort(idx.begin(), idx.end(),
			[&] (int i1, int i2) {
				return (c.count[i1] > c.count[i2]);
			}
		);
	}

	int getNextChar(){
		if (atEnd())
			return -1;
		return TProfileAlphabet(idx[state++]);
	}

	// end if all chars that occurred were returned
	bool atEnd(){
		return ((state == char_size) || (c.count[idx[state]] == 0));
	}
};

/* ------------------------------------------------------- */

template <typename TBidirectionalIndex>
class MotifIterator{
	typedef typename seqan::Iterator<TBidirectionalIndex, seqan::TopDown<seqan::ParentLinks<> > >::Type TIterator;
	typedef typename seqan::SAValue<TBidirectionalIndex>::Type THitPair;
	typedef seqan::String< THitPair > TOccurenceString;

	typedef std::shared_ptr<ProfileCharIter> ProfilePointer;

	// Search state to keep track of the left/right borders of the extension.
	// In the case of stems, only one object, the Binucleotide iterator.
	int id;
	std::stack<ProfilePointer> state;
	//int l, r;
	int pos;
	TStructure structure;
	TIterator it;
	int active_element;
	bool cont = true;

	// structure for reference


	// threshold: stop expanding when below this likelihood for the sequence
	unsigned min_match;

	bool setEnd(){
		return (cont = false);
	}

	bool stepIterator(int next_char, StructureType stype){
		bool success = false;

		if (stype == STEM){
			// partition stem char into components (using / and % instead
			// of bit shifting since AlphabetSize might not be power of 2)
			int lchar = next_char / AlphabetSize;
			int rchar = next_char % AlphabetSize;

			//std::cout << TBaseAlphabet(lchar) << " " << TBaseAlphabet(rchar) << " " << stem_char << "\n";

			// needs to extend into both directions:
			bool wentLeft  = seqan::goDown(it, lchar, seqan::Fwd());
			bool wentRight = seqan::goDown(it, rchar, seqan::Rev());

			success = wentLeft && wentRight;

			if (wentLeft ^ wentRight)
				seqan::goUp(it);

		}
		else{
			success = goDown(it, next_char, seqan::Rev());
		}

		return success;
	}

public:
	MotifIterator(TStructure &structure, TBidirectionalIndex &index, double min_match, unsigned id)
		: id(id), structure(structure), it(index),active_element(structure.elements.size()-1), min_match(min_match){
		// push center elements of hairpin onto the stack
		StructureElement &hairpinElement = structure.elements.back();
		auto &hairpin = hairpinElement.loopComponents[0];
		pos = 0;
		ProfilePointer posPointer(new ProfileCharIterImpl<TAlphabetProfile>(hairpin[pos]));
		state.push(posPointer);
	}

	// next() returns true as long as the motif is not exhausted.
	// Only 'valid' matches are iterated: exclude those who do not match
	// or do not represent the family well (prob. below threshold)
	bool next(){
		if (!cont)
			return false;

		// restore state from last call of next
		StructureType stype = structure.elements[active_element].type;
		ProfilePointer statePointer = state.top();

		// the iterator still points to the previous match, backtrack from that
		if (!statePointer->gapped){
			seqan::goUp(it);
			// for stem pairs, need to backtrack right too
			if (stype == STEM)	seqan::goUp(it);
		}

		// loop until a new match is found and save the state
		// In each iteration, extend one step.
		// When no further extension possible, check if long enough.

		do {
			stype = structure.elements[active_element].type;
			int n = seqan::length(structure.elements[active_element].loopComponents[0]);

			// extend to the right from pos 0
			int next_char = statePointer->getNextChar();

			// skip gap characters (i.e. don't consider them while searching)
			if (next_char == (AlphabetSize-1))
				next_char = statePointer->getNextChar();

			//std::cout << next_char << " " << TBaseAlphabet(next_char) << "\n";

			//statePointer->gapped = (next_char == AlphabetSize-1);
			//std::cout << "(" << next_char << "," << pos << ") ";

			// if chars exhausted, backtrack
			if (next_char == -1){
				--pos;
				state.pop();
				// if we have backtracked to the start, no sequences left
				if (state.empty()){
					return setEnd();
				}

				// reset iterator (if we didn't come from a gap) and get previous state
				if (statePointer->inner){
					statePointer = state.top();

					// only non-gap characters advanced the iterator
					// only goUp if the character wasn't a gap
					//if (!statePointer->gapped)
					seqan::goUp(it);
					continue;
				}
				else
					break;
			}

			// advance to the next character if we have a gap or match
			//if (statePointer->gapped || goDown(it, next_char, seqan::Rev())){
			if (stepIterator(next_char, stype)){
				statePointer->inner = true;
				// end of this descriptor region, swap to next one
				if (pos == n-1){
					break;
					--active_element;

					// swap differently if we have a stem following (binucleotide alphabet)
					pos = 0; //seqan::length(structure.elements[active_element].loopComponents[0])-1;
					auto next_char = structure.elements[active_element].stemProfile[pos];
					ProfilePointer posPointer(new ProfileCharIterImpl<TBiAlphabetProfile>(next_char));
					state.push(posPointer);
					statePointer = posPointer;
				}
				// else just increment position
				else{
					++pos;
					auto next_char = structure.elements[active_element].loopComponents[0][pos];
					ProfilePointer posPointer(new ProfileCharIterImpl<TAlphabetProfile>(next_char));
					state.push(posPointer);
					statePointer = posPointer;
				}
			}

			/*
			if (stype == HAIRPIN){
				//std::cout << this->id << ": " << seqan::representative(it) << "\n";

				// extend to the right from pos 0
				int next_char = statePointer->getNextChar();

				// skip gap characters (i.e. don't consider them while searching)
				if (next_char == (AlphabetSize-1))
					next_char = statePointer->getNextChar();

				//std::cout << next_char << " " << TBaseAlphabet(next_char) << "\n";

				statePointer->gapped = (next_char == AlphabetSize-1);
				//std::cout << "(" << next_char << "," << pos << ") ";

				// if chars exhausted, backtrack
				if (next_char == -1){
					--pos;
					state.pop();
					// if we have backtracked to the start, no sequences left
					if (state.empty()){
						return setEnd();
					}

					// reset iterator (if we didn't come from a gap) and get previous state
					if (statePointer->inner){
						statePointer = state.top();

						// only non-gap characters advanced the iterator
						// only goUp if the character wasn't a gap
						if (!statePointer->gapped)
							seqan::goUp(it);

						continue;
					}
					else
						break;
				}

				// advance to the next character if we have a gap or match
				if (statePointer->gapped || goDown(it, next_char, seqan::Rev())){
					statePointer->inner = true;

					// if we exit the hairpin, go to stem
					// (currently disabled)
					if (pos == n-1){
						//break;
						--active_element;
						pos = 0; //seqan::length(structure.elements[active_element].loopComponents[0])-1;
						auto next_char = structure.elements[active_element].stemProfile[pos];
						ProfilePointer posPointer(new ProfileCharIterImpl<TBiAlphabetProfile>(next_char));
						state.push(posPointer);
						statePointer = posPointer;
					}
					else{
						++pos;
						auto next_char = structure.elements[active_element].loopComponents[0][pos];
						ProfilePointer posPointer(new ProfileCharIterImpl<TAlphabetProfile>(next_char));
						state.push(posPointer);
						statePointer = posPointer;
					}
				}
			}
			else if (stype == STEM){
				//std::cout << this->id << ": " << seqan::representative(it) << "\n";

				int stem_char = statePointer->getNextChar();

				// backtrack if this stem char is exhausted
				// (backtrack: go to previous char, adjust index and rewind iterator)
				if (stem_char == -1){
					--pos;
					state.pop();

					// if backtracked back into the hairpin, change active element
					// and report the match ending at the last hairpin element
					if (pos == -1){
						++active_element;
						pos = seqan::length(structure.elements[active_element].loopComponents[0])-1;
						break;
					}

					// if it's an inner node, simply backtrack
					// (a match was found in the children)
					if (statePointer->inner){
						seqan::goUp(it);
						seqan::goUp(it);
						statePointer = state.top();
						continue;
					} // else break to report a match up to this position
					else{
						break;
					}
				}

				// partition stem char into components (using / and % instead
				// of bit shifting since AlphabetSize might not be power of 2)
				int lchar = stem_char / AlphabetSize;
				int rchar = stem_char % AlphabetSize;

				//std::cout << TBaseAlphabet(lchar) << " " << TBaseAlphabet(rchar) << " " << stem_char << "\n";

				// needs to extend into both directions:
				bool wentLeft  = seqan::goDown(it, lchar, seqan::Fwd());
				bool wentRight = seqan::goDown(it, rchar, seqan::Rev());

				// if we went both left and right, go to next stem pair
				if (wentLeft && wentRight){
					// end of first stem, stop searching
					statePointer->inner = true;
					if (pos == n-1){
						break;
					}
					else{
						++pos;
						auto next_char = structure.elements[active_element].stemProfile[pos];
						ProfilePointer next_ptr(new ProfileCharIterImpl<TBiAlphabetProfile>(next_char));
						state.push(next_ptr);
						statePointer = next_ptr;
					}
				}
				// if we only successfully went into one direction, reset this one
				else if (wentLeft || wentRight){
					seqan::goUp(it);
				}
			}
			*/

		} while(true);

		// if match too short, start next again
		if (seqan::repLength(it) < min_match)
			return this->next();

		//std::cout << "EndpointH: " << pos << " " << seqan::representative(it) << " " <<  seqan::repLength(it) << "\n";


		if (active_element == seqan::length(structure.elements)-1)
			std::cout << "EndpointH: " << pos << " " << seqan::representative(it) << " " <<  seqan::repLength(it) << "\n";
		else{
			int ll = seqan::length(structure.elements[active_element].loopComponents[0]);
			std::cout << "EndpointS: " << seqan::length(structure.elements.back().loopComponents[0])-1 + (ll-pos) << " " << seqan::representative(it) << " " << seqan::repLength(it) << "\n";
		}

		return true;
	}

	TOccurenceString getOccurrences(){
		// FIXME: getting occurrences directly via
		// occs = seqan::getOccurrences(it) doesn't work for some reason
		TOccurenceString occs;
		for (THitPair test: seqan::getOccurrences(it))
			seqan::appendValue(occs, test);

		return occs;
	}

	unsigned countOccurrences(){
		return seqan::countOccurrences(it);
	}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TBidirectionalIndex>
std::vector<TProfileInterval> getStemloopPositions(TBidirectionalIndex &index, Motif &motif){
	//typedef typename seqan::SAValue<TBidirectionalIndex>::Type THitPair;
	typedef seqan::String< TIndexPosType > TOccurenceString;

	//std::vector<TOccurenceString> result(profile.size());

	// create one interval tree for each contig in the reference genome
	std::vector<TProfileInterval> intervals(seqan::countSequences(index));

	int id = 0;
	int n = seqan::length(motif.seedAlignment);
	int stems = motif.profile.size();

	for (TStructure &structure : motif.profile){
		// start of the hairpin in the whole sequence
		std::cout <<  structure.pos.first << " " << structure.pos.second << ": " << motif.profile.size() << "\n";

		int loc = motif.profile[id].elements.back().location;

		MotifIterator<TBidirectionalIndex> iter(structure, index, 11, id);

		std::cout << "Starting iterator\n";

		while (iter.next()){
			TOccurenceString occs = iter.getOccurrences();

			for (TIndexPosType pos : occs){
				TProfileInterval &interval = intervals[pos.i1];

				seqan::String<TProfileCargo> hits;

				// check if the stem occurs in an already existing match region
				findIntervals(hits, interval, pos.i2);

				// if this stem isn't located in an existing match region
				if (seqan::length(hits) == 0){
					//std::cout << loc << " " << n << " " << pos.i2 << " " << pos.i2-loc << " " << pos.i2 + (n-loc) << "\n";

					std::shared_ptr<std::vector<bool> > stemSet(new std::vector<bool>(stems));
					(*stemSet)[id] = true;
					// add an interval around the matched stem
					seqan::addInterval(interval, pos.i2-loc, pos.i2+(n-loc), stemSet);
				}
				// else add to the list of hits
				else{
					for (unsigned i=0; i < seqan::length(hits); ++i){
						(*hits[i].cargo)[id] = true;
					}
				}
			}
			//std::cout << "\n";

		}

		for (unsigned i=0; i < seqan::countSequences(index); ++i)
			std::cout << "After " << id << "," << i << ": " << intervals[i].interval_counter << "\n";

		//std::cout << seqan::length(result[id]) << "\n";
		//seqan::sort(result[id]);
		++id;
	}

	return intervals;
}

void countHits(TProfileInterval positions, int window_size){
	seqan::String<TProfileCargo> hits;
	seqan::getAllIntervals(hits, positions);

	for (unsigned i=0; i < seqan::length(hits); ++i)
		// only report hits where all stems occurred in the region
		if (std::all_of(hits[i].cargo->begin(), hits[i].cargo->end(), [](bool v) { return v; }))
			std::cout << hits[i].i1 << " " << hits[i].i2 << "\n";
}


template <typename TStringType>
std::vector<seqan::Tuple<int, 3> > findFamilyMatches(seqan::StringSet<TStringType> &seqs, std::vector<Motif> &motifs){
	std::vector<seqan::Tuple<int, 3> > results;

	TBidirectionalIndex index(seqs);

    for (Motif &motif : motifs){
    	std::cout << motif.header.at("ID") << "\n";
    	// find the locations of the motif matches
    	//std::cout << motif.header.at("AC") << "\n";
    	std::vector<TProfileInterval> result = getStemloopPositions(index, motif);

		// cluster results into areas (i.e. where hairpins of a given type cluster together)
    	std::cout << result[1].interval_counter << "\n";
    	//countHits(result[1], motif.consensusStructure.size());
    	//for (TProfileInterval intervals : result)
    		//countHits(intervals, motif.consensusStructure.size());
    }

	return results;
}

#endif  // #ifndef APPS_RNAMOTIF_MOTIF_SEARCH_H_
