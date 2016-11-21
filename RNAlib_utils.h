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
#include "motif.h"
#include "stockholm_file.h"
#include <functional>

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
	#include  <ViennaRNA/PS_dot.h>
	#include  <ViennaRNA/aln_util.h>
	#include  <ViennaRNA/Lfold.h>
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
void WUSStoPseudoBracket(std::string const & structure, char* pseudoBracketString){
	int i=0;

	for (char const & c: structure){
		auto iter = WUSS_Table.find(c);
		char replace = (iter == WUSS_Table.end()) ? '.' : iter->second;
		pseudoBracketString[i++] = replace;
	}

	// terminate pseudo bracket string with 0
	pseudoBracketString[i] = 0;
}

void structureToInteractions(const char * const structure, TInteractionPairs &interactions){
	short * struc_table = vrna_ptable(structure);

	interactions.resize(struc_table[0]);
	int j = 0;
	for (size_t i=1; i <= (size_t)struc_table[0]; i++){
		if (struc_table[i] == 0)
			interactions[i-1] = std::make_pair(MAX,-1);
		else
			interactions[i-1] = std::make_pair(ROUND, struc_table[i]-1);
		j++;
	}

	free(struc_table);
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
	//printf("%s %zu\n%s (%6.2f)\n", seq, seq_size, weird_structure, gibbs);

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
	structureToInteractions(structure, interPairs);

	free(pl1);
	free(structure);
	free(weird_structure);
	vrna_fold_compound_free(vc);
}

void getConsensusStructure(Motif &motif, seqan::StockholmRecord<TBaseAlphabet> const & record, const char* constraint, RNALibFold const &){
//void getConsensusStructure(Motif &motif, seqan::StockholmRecord<TBaseAlphabet> const & record, TInteractionPairs &consensusStructure, const char* constraint, RNALibFold const &){
	TInteractionPairs &consensusStructure = motif.consensusStructure;
	int n_seq = record.seqences.size();
   	char** seqs = new char*[n_seq+1];
   	seqs[n_seq] = 0;

   	{
		int i = 0;
		for (auto elem : record.seqences)
		{
			// TODO: Get row out of the alignment object? (not always Stockholm)
			seqs[i] = new char[elem.second.size()+1];
			std::strcpy(seqs[i], elem.second.c_str());
			++i;
		}
   	}

	vrna_init_rand();
	int length = strlen(seqs[0]);

	char *structure  	 = (char*)vrna_alloc(sizeof(char) * (length + 1));
	char *prob_structure = (char*)vrna_alloc(sizeof(char) * (length + 1));

	vrna_md_t md;
	vrna_md_set_default(&md);
	md.uniq_ML = 1;
	vrna_fold_compound_t *vc 	= vrna_fold_compound_comparative((const char**)seqs, &md, VRNA_OPTION_MFE | VRNA_OPTION_PF);

	// add constraints if available
	if (constraint){
		vrna_constraints_add(vc, constraint, VRNA_CONSTRAINT_DB | VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_RND_BRACK);
		std::cout << "bracket:" << constraint << std::endl;
	}

	double min_en = vrna_mfe(vc, structure);

	/* rescale parameters for Boltzmann factors */
	vrna_exp_params_rescale(vc, &min_en);

	float energy, kT;
	float betaScale = 1.;

	kT = (betaScale*((temperature+K0)*GASCONST))/1000.; /* in Kcal */

    energy = vrna_pf(vc, prob_structure);

    std::unordered_map<unsigned, bool> seenStructures;

	std::cout << "Vienna: " << structure << " Freq.: " << std::exp((energy-min_en)/kT) << " " << vrna_mean_bp_distance(vc) << "\n";
	std::cout << "        " << consens_mis((const char**)seqs) << "\n";
//	std::cout << "Weird:  " << prob_structure << "\n";

	std::cout << "\n";

	// sample sequences
	//std::cout << "Alternatives: \n";

	double sum = 0;
	for (int i=0; i<5; i++) {
	//while (sum < 0.9){
		char *s;
		double prob=1.;
		s = vrna_pbacktrack(vc);

		double e  = (double)vrna_eval_structure(vc, s);
		e -= (double)vrna_eval_covar_structure(vc, s);
		prob = exp((energy - e)/kT);

		//double index = std::round(15*e)/15;
		unsigned h = std::hash<std::string>()(std::string(s));
		if (!seenStructures[h]){
			sum += prob;
			printf("        %s ", s);
			printf("%6g %.2f  ",prob, -1*(kT*log(prob)-energy));
			printf("\n");

	//		vrna_eval_loop_pt(vc, 15, struc_table);

			seenStructures[h] = true;
		}
		else
			--i;

		free(s);
    }

	std::cout << "SUM: " << sum << "\n";


	structureToInteractions(structure, consensusStructure);
	TStemLoopProfile stemLoops = findStemLoops(consensusStructure);

	for (auto pair : stemLoops){
		std::cout << pair.pos.first << " " << pair.pos.second << "\n";
	}

	FLT_OR_DBL* probs = vc->exp_matrices->probs;
	int *iindex = vc->iindx;
	int n = vc->length;

	std::unordered_map<int, std::tuple<int,int,int> > hairpins;
	std::vector<FLT_OR_DBL> diagonals(2*n, -1);

	FLT_OR_DBL threshold = 0.2;

	// locate potential hairpins
	// iterate along the matrix diagonals
	for (int k=1; k<length; k++){
		for (int offset=0; offset <= 1; ++offset){
			int i = k-1+offset;
			int j = k+1;

			bool inHairpin = false;

			std::tuple<int,int,int> hairpin;

			while (i > 0 && j <= length){
				FLT_OR_DBL prob = (float)probs[iindex[i] - j];

				if (prob > threshold && !inHairpin){
					inHairpin = true;
					hairpin = std::make_tuple(i,j,1);
				}

				if (prob < threshold && inHairpin){
					std::get<2>(hairpin) = j-std::get<1>(hairpin);
					break;
				}

				i -= 1;
				j += 1;
			}

			if (inHairpin){
				int k = std::get<0>(hairpin) + std::get<1>(hairpin);
				hairpins[k] = hairpin;
			}
		}
	}

	TStemLoopProfile result_regions;

	bool first = true;

	int count = 0;
	int max_count = hairpins.size();

	// iteratively enforce most likely bases from hairpins that haven't been encountered yet
	do {
		// create new structure
		if (!first){
			//int k = unused_hairpins.front();
			int k = hairpins.begin()->first;

			int i,j,l;
			std::tie(i,j,l) = hairpins[k];

			// enforce base pair with highest probability in hairpin
			int maxi, maxj;
			float maxp = -1;

			for (int off=0; off <= l; ++off){
				FLT_OR_DBL prob = (float)probs[iindex[i-off] - (j+off)];
				if (prob > maxp){
					maxp = prob;
					maxi = i-off;
					maxj = j+off;
				}
			}

			vrna_hc_init(vc);
			//std::cout << maxi << " " << maxj << "\n";
			vrna_hc_add_bp(vc, maxi, maxj, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);

			vrna_mfe(vc, structure);
			double e  = (double)vrna_eval_structure(vc, structure);
			e -= (double)vrna_eval_covar_structure(vc, structure);

			std::cout << "Vienna: " << structure << " Freq.: " << std::exp((energy-e)/kT) << " " << vrna_mean_bp_distance(vc) << "\n";

			structureToInteractions(structure, consensusStructure);
			stemLoops = findStemLoops(consensusStructure);

			for (auto pair : stemLoops){
				std::cout << pair.pos.first << " " << pair.pos.second << "\n";
			}
		}

		first = false;
		std::vector<TInteractionPairs> structureVariants;

		// check with which hairpins our stemLoops overlap
		for (auto pair : stemLoops){
			bool stemAdded = false;
			std::unordered_map<int, std::tuple<int,int,int> > tmpHairpin = hairpins;
			for (auto hairpin : tmpHairpin){
				int k = hairpin.first;
				int i,j,l;
				std::tie(i,j,l) = hairpin.second;

				int posi = i-l;
				int posj = j+l;

				// check if the hairpin region overlaps with any stemloop
				for (int off=0; off <= l; ++off){
					if ((pair.pos.first < posi+off) && (posj-off < pair.pos.second)){
						std::cout << pair.pos.first << " " << pair.pos.second << "\n";
						partitionStemLoop(motif.seedAlignment, consensusStructure, pair);

						if (!stemAdded){
							stemAdded = true;
							result_regions.push_back(pair);
						}

						std::cout << "Hairpin from: (" << i-l << "," << j+l << ") to (" << i << "," << j << ")\n";
						hairpins.erase(k);

						break;
					}
				}
			}
		}
	} while (!hairpins.empty() && ++count < max_count);

	std::cout << "Regions\n";
	for (auto pair : result_regions){
		std::cout << pair.pos.first << " " << pair.pos.second << "\n";
	}

	// get dot plot structures
	vrna_plist_t *pl1, *pl2;
	pl1 = vrna_plist_from_probs(vc, 0.005);
	pl2 = vrna_plist(structure, 0.95*0.95);

	//DEBUG_MSG("Vienna: " << structure);
	//DEBUG_MSG("        " << consens_mis((const char**)seqs));

	//structureToInteractions(structure, consensusStructure);

	// write dot-plot
	//	Function used to plot the dot_plot graph
	//(void) PS_dot_plot_list((char*)seqs[0], (record.header.at("AC") + std::string(".ps")).c_str(), pl1, pl1, "");
	char *tmp = (char*)(record.header.at("AC") + std::string(".ps")).c_str();
	(void) PS_dot_plot_list((char*)consens_mis((const char**)seqs), tmp, pl1, pl2, structure);

	// free all used RNAlib data structures
	free(structure);
	free(prob_structure);
	vrna_fold_compound_free(vc);

	// free the sequence array (terminated with a 0 at the end)
	for (size_t k = 0; seqs[k] != 0; ++k)
		free(seqs[k]);

	free(seqs);}

#endif  // #ifndef APPS_RNAMOTIF_RNALIB_UTILS_H_
