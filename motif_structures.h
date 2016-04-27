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

typedef std::vector<int> TInteractionPairs;
typedef std::vector<std::pair<int, int > > TStructureRegions;

// Alphabets usually found in Stockholm format files are Rna or AminoAcid
template <typename TAlphabet>
struct StockholmRecord {
	// SeqAn specific data	==========================
	// store the alignment given in the Stockholm file
	typedef seqan::String<TAlphabet> TSequence;
	typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;

	TAlign alignment;

	// Raw string data		===========================
	// header (GF tags -> tag value)
	std::unordered_map<std ::string, std::string > header;
	// seqence names -> sequence maps
	std::unordered_map<std::string, std::string > seqences;
	// per column (GC) -> annotation string
	std::unordered_map<std::string, std::string > seqence_information;

	//TODO: Maybe support per-sequence (GS) and per-residue (GR) annotation
};

struct InteractionGraph {
	std::vector<TUVertexDescriptor > vertices;
	TUgraph graph; // this graph represents all the computed interaction edges
};

// Describes a type of secondary structure and its statistics (varying length, base probabilities, etc.)
// Can be: Stem, Bulge, Internal Loop, Hairpin Loop
struct StructureElement{

};

struct Motif{
	// stores interaction information for the N sequences in a N-vector
	// TODO: graph
	std::vector<InteractionGraph> interactionGraphs; // a graph of interaction probabilities
	std::vector<TInteractionPairs> interactionPairs; // fixed structure predictions
	TInteractionPairs consensusStructure;

	// the alignment structure of the seed RNA family
	TAlign seedAlignment;

	// the regions identified as hairpins. Nested hairpin structures are stored as successive pairs.
	TStructureRegions hairpinLoops;

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef APPS_RNAMOTIF_MOTIF_STRUCTURES_H_
