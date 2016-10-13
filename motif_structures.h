// ==========================================================================
//                             motif_structures.h
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

#ifndef APPS_RNAMOTIF_MOTIF_STRUCTURES_H_
#define APPS_RNAMOTIF_MOTIF_STRUCTURES_H_

// SeqAn headers
#include <seqan/align.h>
#include <seqan/index.h>
#include "stockholm_file.h"

// C++ headers
#include <unordered_map>
#include <stack>
#include <numeric>
#include <memory>

#ifndef NDEBUG
#define DEBUG_MSG(str) do { std::cout << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/* Various types related to secondary structure */

typedef enum {ROUND = 0, SQUARE, CURLY, ANGLE, MAX} BracketType;

typedef enum {HAIRPIN = 1, STEM = 2, LBULGE = 4, RBULGE = 8, LOOP = 12} StructureType;

const std::string parentheses("{}[]()<>");

std::unordered_map<char, BracketType> bracket_to_type =
	{
	 {'(', ROUND},  {')', ROUND},
	 {'[', SQUARE}, {']', SQUARE},
	 {'{', CURLY},  {'}', CURLY},
	 {'<', ANGLE},  {'>', ANGLE}
	};

std::unordered_map<char, char> match_table =
	{
	 {')', '('}, {'>', '<'},
	 {']', '['}, {'}', '{'}
	};

/* Various functions related to brackets and bases */

bool isOpen(char c){
	return (c == '(') || (c == '[') || (c == '{') || (c == '<');
}


// Types for the alignment
typedef seqan::String<seqan::Rna> TSequence;
typedef seqan::StringSet<TSequence> TStringSet;
typedef seqan::StringSet<TSequence, seqan::Dependent<> > TDepStringSet;
typedef seqan::Graph<seqan::Alignment<TDepStringSet, void, seqan::WithoutEdgeId> > TAlignGraph;

typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;      // align type
typedef seqan::Row<TAlign>::Type TRow;
typedef seqan::Iterator<TRow>::Type TRowIterator;


// Types for the interaction graphs for each sequence
typedef float TCargo;
typedef seqan::Graph<seqan::Undirected<TCargo> > TUgraph;
typedef seqan::VertexDescriptor<TUgraph>::Type TUVertexDescriptor;

struct vectGraphElement {
	std::vector<TUVertexDescriptor > uVertexVect;
	TUgraph interGraph; // this graph represents all the computed interaction edges
};

// tag for the type of folding algorithm
struct RNALibFold__;
typedef seqan::Tag<RNALibFold__> RNALibFold;

struct IPknotFold__;
typedef seqan::Tag<IPknotFold__> IPknotFold;

struct InteractionGraph {
	std::vector<TUVertexDescriptor > vertices;
	TUgraph graph; // this graph represents all the computed interaction edges
};

// From the Alphabet used (Dna, Rna), define a Binucleotide Alphabet and Profile strings.
// Binucleotide alphabet: // AA, AC, AG, AT
typedef seqan::Rna TAlphabet;

const seqan::ValueSize<TAlphabet>::Type AlphabetSize = seqan::ValueSize<TAlphabet>::VALUE;
typedef seqan::SimpleType<unsigned char, seqan::Finite<AlphabetSize*AlphabetSize> > TBiAlphabet;

typedef seqan::ProfileChar<TAlphabet> TAlphabetProfile;
typedef seqan::ProfileChar<TBiAlphabet> TBiAlphabetProfile;

typedef seqan::String<seqan::ProfileChar<TAlphabet> > TLoopProfileString;
typedef seqan::String<seqan::ProfileChar<TBiAlphabet> > TStemProfileString;

struct StructureStatistics{
	unsigned min_length;
	unsigned max_length;
	double  mean_length;
};


// Describes a type of secondary structure
// Stem, Internal Loop, Hairpin Loop
struct StructureElement{
	StructureType type;

	// for HAIRPIN: just one string
	// for STEM   : two strings (left and right side of the stem)
	// for LOOP   : two strings (left and right side of the loop)
	//			    if one side is empty, loop is a bulge
	std::vector<TLoopProfileString> loopComponents;
	TStemProfileString stemProfile;

	// length of the sequence that makes up the structure
	// the variation in length comes due to gaps in the alignment
	std::vector<StructureStatistics> statistics;
};

typedef std::vector<std::pair<BracketType, int> > TInteractionPairs;
typedef std::pair<int, int> TRegion;

typedef std::vector<std::pair<BracketType, TRegion> > TSequenceRegions;
typedef std::vector<StructureElement> TStructure;
typedef std::vector<TStructure> TStemLoopProfile;

struct Motif{
	// header data, same as the Stockholm header.
	std::unordered_map<std ::string, std::string > header;

	// stores interaction information for the N sequences in a N-vector
	std::vector<InteractionGraph> interactionGraphs; // a graph of interaction probabilities
	std::vector<TInteractionPairs> interactionPairs; // fixed structure predictions

	// the consensus interactions for the whole alignment
	TInteractionPairs consensusStructure;

	// the alignment structure of the seed RNA family
	TAlign seedAlignment;

	// the profiles for the individual stem loops
	TStemLoopProfile profile;
};

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

	bool atEnd(){
		return (state == char_size);
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
	std::stack<ProfilePointer> state;
	int l, r;
	int pos;
	TStructure structure;
	TIterator it;
	int active_element;
	bool cont = true;
	bool search_limit = false;

	// structure for reference


	// threshold: stop expanding when below this likelihood for the sequence
	int min_match;

	bool setEnd(){
		return (cont = false);
	}

public:
	MotifIterator(TStructure &structure, TBidirectionalIndex &index, double min_match = 11)
		: structure(structure), it(index),active_element(structure.size()-1), min_match(min_match){
		// push center elements of hairpin onto the stack
		StructureElement &hairpinElement = structure.back();
		auto &hairpin = hairpinElement.loopComponents[0];
		pos = 0;
		ProfilePointer posPointer(new ProfileCharIterImpl<TAlphabetProfile>(pos));
		state.push(posPointer);

		//l = n/2-1;
		//r = n/2;

		//auto left_char  = hairpin[n/2-1];
		//auto right_char = hairpin[n/2+1];

		//ProfilePointer left(new ProfileCharIterImpl<TAlphabetProfile>(left_char));
		//ProfilePointer right(new ProfileCharIterImpl<TAlphabetProfile>(right_char));

		//state.push(left);
		//rightState.push(right);
	}

	// next() returns true as long as the motif is not exhausted.
	// Only 'valid' matches are iterated: exclude those who do not match
	// or do not represent the family well (prob. below threshold)
	bool next(){
		if (!cont)
			return false;

		// restore state from last call of next
		ProfilePointer statePointer = state.top();

		bool loopCondition = true;

		// loop until a new match is found and save the state
		// In each iteration, extend one step.
		// When no further extension possible, check if long enough.

		do {
			StructureType stype = structure[active_element].type;
			int n = seqan::length(structure[active_element].loopComponents[0]);

			if (stype == HAIRPIN){
				// extend to the right from pos 0
				int next_char = statePointer->getNextChar();

				// if chars exhausted, backtrack
				if (next_char == -1){
					--pos;
					state.pop();
					// if we have backtracked to the start, no sequences left
					if (state.empty()){
						return setEnd();
					}

					// reset iterator and get previous state
					seqan::goUp(it);
					statePointer = state.top();

					continue;
				}

				// advance to the next character
				if ( (loopCondition = goDown(it, next_char, seqan::Fwd() )) ){
					++pos;
					// if we exit the hairpin, go to stem
					int next_char;
					if (pos == n){
						--active_element;
						pos = seqan::length(structure[active_element])-1;
						next_char = structure[active_element].stemProfile[pos];
					}
					else
						next_char = structure[active_element].loopComponents[0][pos];

					ProfilePointer posPointer(new ProfileCharIterImpl<TAlphabetProfile>(next_char));
					state.push(posPointer);
					statePointer = posPointer;
				}


			}
			else if (stype & STEM){
				// if we are at the end of the stem, set iterator back left and right
				if (search_limit){
					seqan::goUp(it);
					seqan::goUp(it);
					search_limit = false;
				}

				int stem_char = statePointer->getNextChar();

				// backtrack if this stem char is exhausted
				// (backtrack: go to previous char, adjust index and rewind iterator)
				if (stem_char == -1){
					++pos;
					state.pop();
					statePointer = state.top();

					// if backtracked back into the hairpin, change active element
					if (pos == n){
						++active_element;
						pos = n-1;
					}

					continue;
				}

				// partition stem char into its components
				int lchar = stem_char >> AlphabetSize;
				int rchar = stem_char & AlphabetSize;

				// needs to extend into both directions:
				bool wentLeft  = seqan::goDown(it, lchar, seqan::Rev());
				bool wentRight = seqan::goDown(it, rchar, seqan::Fwd());

				// if we went both left and right
				if ( (loopCondition = (wentLeft && wentRight)) ){
					if (pos == 0){
						loopCondition = false;
						search_limit = true;
					}
					else{
						--pos;
						auto next_char = structure[active_element].stemProfile[pos];
						ProfilePointer next_ptr(new ProfileCharIterImpl<TBiAlphabetProfile>(next_char));
						state.push(next_ptr);
						statePointer = next_ptr;
					}
				}
				//
				else if (wentLeft || wentRight){
					seqan::goUp(it);
				}
			}
		} while(loopCondition);

		// if match too short, start next again
		if (seqan::repLength(it) < min_match)
			return this->next();

		if (active_element == seqan::length(structure)-1)
			std::cout << "Endpoint: " << pos << " " << seqan::repLength(it) << "\n";
		else
			std::cout << "Endpoint: " << seqan::length(structure.back().loopComponents[0])-1 + pos << " " << seqan::repLength(it) << "\n";

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

#endif  // #ifndef APPS_RNAMOTIF_MOTIF_STRUCTURES_H_

/*
if (!lEnd){
					// what is next char to the left?
					int lchar = left->getNextChar();

					// if all characters exhausted, backtrack and try again
					// (backtrack: go to previous char, adjust index and rewind iterator)
					if (lchar == -1){
						++l;
						leftState.pop();
						if (leftState.empty())
							return setEnd();

						left = leftState.top();
						continue;
					}

					// if match extended, move to next profile column
					if (loopCondition |= seqan::goDown(it, lchar, seqan::Rev())){
						lEnd = (l == 0);
						if (!lEnd){
							auto next_char = structure[active_element].loopComponents[0][--l];
							ProfilePointer next_ptr(new ProfileCharIterImpl<TAlphabetProfile>(next_char));
							leftState.push(next_ptr);
							left = next_ptr;
						}
					}
				}

				if (!rEnd){
					// what is next char to the right?
					int rchar = right->getNextChar();

					// if all characters exhausted, backtrack and try again
					// (backtrack: go to previous char, adjust index and rewind iterator)
					if (rchar == -1){
						--r;
						rightState.pop();
						if (rightState.empty())
							return setEnd();

						right = rightState.top();
						continue;
					}

					// if match extended, move to next profile column
					if (loopCondition |= seqan::goDown(it, rchar, seqan::Fwd())){
						rEnd = (r == n-1);
						if (!rEnd){
							auto next_char = structure[active_element].loopComponents[0][++r];
							ProfilePointer next_ptr(new ProfileCharIterImpl<TAlphabetProfile>(next_char));
							rightState.push(next_ptr);
							right = next_ptr;
						}
					};
				}

				// if the hairpin reached the end, move state to stem and set indices
				if (lEnd && rEnd){
					--active_element;
					l = seqan::length(structure[active_element].stemProfile)-1;
					r = l;
					auto next_char = structure[active_element].stemProfile[l];
					ProfilePointer next_ptr(new ProfileCharIterImpl<TBiAlphabetProfile>(next_char));
					leftState.push(next_ptr);
					left = next_ptr;
				}*/
