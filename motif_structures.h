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
#include <unordered_map>
#include "stockholm_file.h"

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

typedef seqan::String<seqan::ProfileChar<seqan::Rna> > RNAProfileString;
typedef enum {HAIRPIN = 1, STEM = 2, LBULGE = 4, RBULGE = 8, LOOP = 12} StructureType;

struct StructureStatistics{
	unsigned min_length;
	unsigned max_length;
	double  mean_length;
};

// Describes a type of secondary structure
// Stem, Internal Loop, Hairpin Loop
struct StructureElement{
	StructureType type;
	typedef RNAProfileString TProfileString;

	// for HAIRPIN: just one string
	// for STEM   : two strings (left and right side of the stem)
	// for LOOP   : two strings (left and right side of the loop)
	//			    if one side is empty, loop is a bulge
	std::vector<TProfileString> components;

	// length of the sequence that makes up the structure
	// the variation in length comes due to gaps in the alignment
	std::vector<StructureStatistics> statistics;
};

std::unordered_map<char, char> match_table =
	{
	 {')', '('}, {'>', '<'},
	 {']', '['}, {'}', '{'}
	};

typedef enum {ROUND = 0, SQUARE, CURLY, ANGLE, MAX} BracketType;

std::unordered_map<char, BracketType> bracket_to_type =
	{
	 {'(', ROUND},  {')', ROUND},
	 {'[', SQUARE}, {']', SQUARE},
	 {'{', CURLY},  {'}', CURLY},
	 {'<', ANGLE},  {'>', ANGLE}
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

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef APPS_RNAMOTIF_MOTIF_STRUCTURES_H_
