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
#include <seqan/misc/interval_tree.h>

#include "stockholm_file.h"

// C++ headers
#include <unordered_map>
#include <unordered_set>
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

namespace seqan {

// specialize seqan's Interval tree 'findIntervals' function to enable it to
// return IntervalAndCargo to include boundaries of the overlapping intervals
template <typename TValue, typename TCargo, typename TValue2>
inline void
findIntervals(
        String<IntervalAndCargo<TValue, TCargo> > & result,
        IntervalTree<TValue, TCargo> const & it,
        TValue2 query)
{
    findIntervals(result, it.g, it.pm, query);
}

template <typename TSpec, typename TPropertyMap, typename TValue, typename TCargo, typename TValue2>
inline void
findIntervals(
        String< IntervalAndCargo<TValue, TCargo> > & result,
        Graph<TSpec> const & g,
        TPropertyMap const & pm,
		TValue2 query)
{
    typedef Graph<TSpec> const TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Value<TPropertyMap>::Type TProperty;
    typedef typename Value<TProperty>::Type TPropertyValue;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

    resize(result, 0);
    if (empty(g))
        return;

    // start at root
    TVertexDescriptor act_knot = 0;
    TProperty act_prop = property(pm, act_knot);
    TProperty next_prop;

    while (true)
    {
        //typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type it7;
        //Iter<Graph<TSpec>, GraphIterator<InternalOutEdgeIterator<OutEdgeIterator> > > it5(g, act_knot);
        //TOutEdgeIterator it4;
        TOutEdgeIterator it(g, act_knot);
        act_prop = property(pm, act_knot);
        if (act_prop.center < (TPropertyValue)query) // look in current node and right subtree
        {
            unsigned int i = 0;
            while (i<length(act_prop.list2) && rightBoundary(value(act_prop.list2, i))>(TPropertyValue) query)
            {
            	auto &val = value(act_prop.list2, i);
            	IntervalAndCargo<TValue, TCargo> tmp(leftBoundary(val), rightBoundary(val), cargo(val));
                appendValue(result, tmp, Generous());
                ++i;
            }
            if (atEnd(it))
                break;

            next_prop = property(pm, targetVertex(it));
            if (next_prop.center <= act_prop.center)
            {
                goNext(it);
                if (atEnd(it))
                    break;
            }
            act_knot = targetVertex(it);
        }
        else
        {
            if ((TPropertyValue)query < act_prop.center) // look in current node and left subtree
            {
                unsigned int i = 0;
                while (i < length(act_prop.list1) && leftBoundary(value(act_prop.list1, i)) <= (TPropertyValue)query)
                {
                	auto &val = value(act_prop.list1, i);
					IntervalAndCargo<TValue, TCargo> tmp(leftBoundary(val), rightBoundary(val), cargo(val));
					appendValue(result, tmp, Generous());
                    ++i;
                }
                if (atEnd(it))
                    break;

                next_prop = property(pm, targetVertex(it));
                if (next_prop.center >= act_prop.center)
                {
                    goNext(it);
                    if (atEnd(it))
                        break;
                }
                act_knot = targetVertex(it);
            }
            else  // look in current node only, as query is center
            {
                for (unsigned int i = 0; i < length(act_prop.list1); ++i){
                	auto &val = value(act_prop.list1, i);
					IntervalAndCargo<TValue, TCargo> tmp(leftBoundary(val), rightBoundary(val), cargo(val));
					appendValue(result, tmp, Generous());
                }
                break;
            }
        }
    }
}

// do an inorder traversal of the tree and report all intervals
template <typename TValue, typename TCargo>
inline void
getAllIntervals(
        String<IntervalAndCargo<TValue, TCargo> > & result,
		IntervalTree<TValue, TCargo, StoreIntervals> const & it)
{
	getAllIntervals(result, it.g, it.pm);
}

template <typename TSpec, typename TPropertyMap, typename TValue, typename TCargo>
inline void
getAllIntervals(
        String<IntervalAndCargo<TValue, TCargo> > & result,
        Graph<TSpec> const & g,
        TPropertyMap const & pm)
{
    typedef Graph<TSpec> const TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	resize(result, 0);
	// start at root
	TVertexDescriptor act_knot = 0;
	getAllIntervals(result, g, pm, act_knot);
}

template <typename TSpec, typename TPropertyMap, typename TValue, typename TCargo, typename TVertexDescriptor>
inline void
getAllIntervals(
        String< IntervalAndCargo<TValue, TCargo> > & result,
        Graph<TSpec> const & g,
        TPropertyMap const & pm,
		TVertexDescriptor & act_knot)
{
	typedef Graph<TSpec> const TGraph;
	typedef typename Value<TPropertyMap>::Type TProperty;
	//typedef typename Value<TProperty>::Type TPropertyValue;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	if (empty(g))
		return;

	TProperty act_prop = property(pm, act_knot);
	TProperty next_prop;

	TOutEdgeIterator it(g, act_knot);

	// go left
	if (!atEnd(it)){
		TVertexDescriptor next_knot = targetVertex(it);
		getAllIntervals(result, g, pm, next_knot);
		goNext(it);
	}

	// append center list
	for (unsigned int i = 0; i < length(act_prop.list1); ++i){
		auto &val = value(act_prop.list1, i);
		IntervalAndCargo<TValue, TCargo> tmp(leftBoundary(val), rightBoundary(val), cargo(val));
		appendValue(result, tmp, Generous());
	}

	// go right
	if (!atEnd(it))
	{
		TVertexDescriptor next_knot = targetVertex(it);
		getAllIntervals(result, g, pm, next_knot);
	}
}

}

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

typedef seqan::Rna5 TBaseAlphabet;

// Types for the alignment
typedef seqan::String<TBaseAlphabet> TSequence;
typedef seqan::StringSet<TSequence> TStringSet;
typedef seqan::StringSet<TSequence, seqan::Dependent<> > TDepStringSet;

typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign; // align type
typedef seqan::Row<TAlign>::Type TRow;
typedef seqan::Iterator<TRow>::Type TRowIterator;

// tag for the type of folding algorithm
struct RNALibFold__;
typedef seqan::Tag<RNALibFold__> RNALibFold;

struct IPknotFold__;
typedef seqan::Tag<IPknotFold__> IPknotFold;

// From the Alphabet used (Dna, Rna), define a Binucleotide Alphabet and Profile strings.
// Binucleotide alphabet: // AA, AC, AG, AT
// Add one field to store gap characters.
typedef seqan::Index<seqan::StringSet<seqan::String<TBaseAlphabet> >, seqan::BidirectionalIndex<seqan::FMIndex< > > >  TBidirectionalIndex;
typedef typename seqan::SAValue<TBidirectionalIndex>::Type TIndexPosType;

// Base Alphabet + gap character for profiles
//typedef seqan::SimpleType<unsigned char, seqan::Finite<seqan::ValueSize<TBaseAlphabet>::VALUE+1> > TAlphabet;
typedef seqan::SimpleType<unsigned char, seqan::Finite<seqan::ValueSize<TBaseAlphabet>::VALUE> > TAlphabet;

const seqan::ValueSize<TAlphabet>::Type AlphabetSize = seqan::ValueSize<TAlphabet>::VALUE;
const seqan::BitsPerValue<TAlphabet>::Type AlphabetBitSize = seqan::BitsPerValue<TAlphabet>::VALUE;
typedef seqan::SimpleType<unsigned char, seqan::Finite<AlphabetSize*AlphabetSize> > TBiAlphabet;

typedef seqan::ProfileChar<TAlphabet> TAlphabetProfile;
typedef seqan::ProfileChar<TBiAlphabet> TBiAlphabetProfile;

typedef seqan::String<seqan::ProfileChar<TAlphabet> > TLoopProfileString;
typedef seqan::String<seqan::ProfileChar<TBiAlphabet> > TStemProfileString;

typedef seqan::IntervalAndCargo<long unsigned int, std::shared_ptr<std::vector<bool> > > TProfileCargo;
typedef seqan::IntervalTree<long unsigned int, std::shared_ptr<std::vector<bool> >, seqan::StoreIntervals> TProfileInterval;

typedef uint32_t THashType;
const size_t HashTabLength = 1*1024*1024;

struct StructureStatistics{
	unsigned min_length;
	unsigned max_length;
	double  mean_length;
};


// Describes a type of secondary structure
// Stem, Internal Loop, Hairpin Loop
struct StructureElement{
	StructureType type;
	bool loopLeft;
	int location;

	// for HAIRPIN: just one string
	// for STEM   : two strings (left and right side of the stem)
	// for LOOP   : two strings (left and right side of the loop)
	//			    if one side is empty, loop is a bulge
	TLoopProfileString loopComponents;
	// gaps that can occur at each position of the descriptor
	std::vector<std::map<int, int> > gap_lengths;
	TStemProfileString stemProfile;

	// length of the sequence that makes up the structure
	// the variation in length comes due to gaps in the alignment
	StructureStatistics statistics;
};

typedef std::vector<std::pair<BracketType, int> > TConsensusStructure;
typedef std::vector<int> TInteractions;

// legacy
typedef std::vector<std::pair<BracketType, std::pair<int, int> > > TSequenceRegions;

typedef struct ProfileStructure{
	BracketType btype;
	std::pair<int, int> pos;

	TInteractions interactions;

	double prob = 0;

	// hairpin, loop, etc. elements of the descriptor
	std::vector<StructureElement> elements;

	ProfileStructure(BracketType btype, std::pair<int, int> pos)
		: btype(btype), pos(pos) {};

	ProfileStructure(BracketType btype, std::pair<int, int> pos, TInteractions interactions)
		: btype(btype), pos(pos), interactions(interactions) {};
} TStructure;

typedef std::vector<TStructure> TStemLoopProfile;

struct Motif{
	// header data, same as the Stockholm header.
	std::unordered_map<std ::string, std::string > header;
	std::unordered_map<std::string, std::string > seqence_information;

	// the consensus interactions for the whole alignment
	//TInteractionPairs consensusStructure;

	// the alignment structure of the seed RNA family
	TAlign seedAlignment;

	// the profiles for the individual stem loops
	TStemLoopProfile profile;

	std::vector<TLoopProfileString> externalBases;

	double mcc = 0;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef APPS_RNAMOTIF_MOTIF_STRUCTURES_H_
