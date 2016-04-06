// ==========================================================================
//                               RNAlib_utils.h
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

#ifndef APPS_RNAMOTIF_RNALIB_UTILS_H_
#define APPS_RNAMOTIF_RNALIB_UTILS_H_

#include "motif_structures.h"

// Import RNAlib. 'extern "C"', since it's a C library
extern "C"{
	#include  <ViennaRNA/data_structures.h>
	#include  <ViennaRNA/mfe.h>
	#include  <ViennaRNA/params.h>
	#include  <ViennaRNA/utils.h>
	#include  <ViennaRNA/eval.h>
	#include  <ViennaRNA/fold.h>
	#include  <ViennaRNA/part_func.h>
	#include  <ViennaRNA/alifold.h>
	#include  <ViennaRNA/constraints.h>
}

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// simplify WUSS notation for use in RNAlib
std::unordered_map<char, char> WUSS_Table =
	{
	 {'(', '('}, {')', ')'}, // (
	 {'<', '('}, {'>', ')'},
	 {'[', '('}, {']', ')'},
	 {'{', '('}, {'}', ')'},
	 {',', '.'}, {'_', '.'},
	 {'-', '.'}, {':', '.'},
	 {'~', '.'}
	};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// convert Rfam's WUSS structure notation to RNAlib's pseudo-bracket notation
// NOTE: RNAlib's pseudo-bracket notation cannot handle pseudoknots
void WUSStoPseudoBracket(std::string& structure, char* pseudoBracketString){
	int i=0;

	for (char& c: structure){
		auto iter = WUSS_Table.find(c);
		char replace = (iter == WUSS_Table.end()) ? '.' : iter->second;
		pseudoBracketString[i++] = replace;
	}

	// terminate pseudo bracket string with 0
	pseudoBracketString[i] = 0;
}

void createInteractions(InteractionGraph &interGraph, TInteractionPairs& interPairs, std::string& seq_str, const char * constraint = NULL){
	size_t seq_size = seq_str.length();

	// add new vertex for each base and save the vertex references
	interGraph.vertices.resize(seq_size);
	interPairs.resize(seq_size);

	for (unsigned i=0; i < seq_size; ++i){
		// add interaction edge in the interaction graph
		interGraph.vertices[i] = seqan::addVertex(interGraph.graph);
	}

	// create RNAlib data structures
	const char * seq = seq_str.c_str();
	char  *structure 		= (char*)vrna_alloc(sizeof(char) * (strlen(seq) + 1));
	char  *weird_structure 	= (char*)vrna_alloc(sizeof(char) * (strlen(seq) + 1));
	vrna_fold_compound_t *vc = vrna_fold_compound(seq, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);

	// add constraints if available
	if (constraint)
		vrna_constraints_add(vc, constraint, VRNA_CONSTRAINT_DB | VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_RND_BRACK);

	// predict secondary structure (will create base pair probs in compound)
	double mfe 	 = (double)vrna_mfe(vc, structure);
	double gibbs = (double)vrna_pf(vc, weird_structure);
	//printf("%s %zu\n%s\n%s (%6.2f)\n", seq, seq_size, structure, weird_structure, gibbs);
	printf("%s %zu\n%s (%6.2f)\n", seq, seq_size, weird_structure, gibbs);

	// extract base pair probabilities from fold compound
	vrna_plist_t *pl1, *ptr;
	pl1 = vrna_plist_from_probs(vc, 0.05);

	// get size of probability list (a full list, if not filtered with a threshold,
	// contains the upper half of the n x n probability matrix (size (n x n)/2 - n)
	unsigned size;
	for(size = 0, ptr = pl1; ptr->i; size++, ptr++);

	// store base pair prob. in the interaction graph
	for(unsigned i=0; i<size;++i){
		seqan::addEdge(interGraph.graph, interGraph.vertices[pl1[i].i-1], interGraph.vertices[pl1[i].j-1], pl1[i].p);
	}

	// save sequence structure
	short * struc_table = vrna_ptable(structure);
	for (size_t i=0; i < seq_size; i++){
		interPairs[i] = struc_table[i];
	}

	free(structure);
	free(weird_structure);
	vrna_fold_compound_free(vc);
}

void getConsensusStructure(const char** seqs, char* structure, const char* constraint = NULL){
	vrna_fold_compound_t *vc = vrna_fold_compound_comparative(seqs, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);

	// add constraints if available
	if (constraint){
		vrna_constraints_add(vc, constraint, VRNA_CONSTRAINT_DB | VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_RND_BRACK);
		std::cout << "Rfam converted:\n" << constraint << std::endl;
	}

	vrna_mfe(vc, structure);
	vrna_fold_compound_free(vc);

	std::cout << "Vienna consensus:\n" << structure << std::endl;

	vc = vrna_fold_compound_comparative(seqs, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
	if (constraint)
		vrna_constraints_add(vc, constraint, VRNA_CONSTRAINT_DB | VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_RND_BRACK);

	vrna_pf(vc, structure);
	std::cout << "Weird consensus:\n" << structure << std::endl;

	vrna_fold_compound_free(vc);
}

#endif  // #ifndef APPS_RNAMOTIF_RNALIB_UTILS_H_
