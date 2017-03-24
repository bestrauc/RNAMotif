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
#include <stack>

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

void structureToInteractions(const char * const structure, TConsensusStructure &interactions){
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

double structureDefinedness(TInteractions& interactions, FLT_OR_DBL* probs, int *iindex){
	double S = 0;
	double prob;
	for (int k=0; k < interactions.size()-1; ++k){
		int i = k+1;
		int j = interactions[i]+1;
		if (j != 0 && i > j)
			continue;

		// TODO: why index error?
		//std::cout << k << "\n";

		//std::cout << i << " " << j << " " << (float)probs[iindex[i] - j] << "\n";

		if (j == 0)
			prob = (1-probs[iindex[i] - j]);
		else
			prob = probs[iindex[i] - j];

		if (prob > 0.01)
			S -= prob*std::log2(prob);
	}

	return S/interactions.size();
}

double getInteractionProbability(TInteractions& interactions, FLT_OR_DBL* probs, int *iindex, int start, int end){
	double prob = 1;
	std::cout << "Prob from " << start << " to " << end << "\n";
	//std::cout << interactions.size() << " size\n";
	for (int k=start; k <= end; ++k){
		// index +1 to correspond to ViennaRNA's numbering
		int i = k+1;
		int j = interactions[k]+1;

		// if we look 'back', we already processed this pair
		if (j != 0 && i > j)
			continue;

		//std::cout << i << " " << j << " " << (float)probs[iindex[i] - j] << "\n";

		if (j == 0)
			prob *= 1;//(1-probs[iindex[i] - j]);
		else
			prob *= probs[iindex[i] - j];
	}

	return prob;
}

char* interactionsToStructure(TInteractions& interactions, int start, int end){
	char *structure  	 = (char*)vrna_alloc(sizeof(char) * (interactions.size() + 1));
	structure[interactions.size()] = 0;

	for (int i = 0; i < interactions.size(); ++i)
		structure[i] = '.';

	for (int i=start; i <= end; ++i){
		int j = interactions[i];

		if (j == -1)
			continue;
		else if (i < j)
			structure[i] = '(';
		else
			structure[i] = ')';
	}

	return structure;
}


double getStructureProbability(TStructure& structure, FLT_OR_DBL* probs, int *iindex){
	return getInteractionProbability(structure.interactions, probs, iindex, structure.pos.first, structure.pos.second);
}

bool checkMatch(TStructure& region, short *strucTab){
//	for (int i=0; i < region.interactions.size(); ++i)
//		std::cout << region.interactions[i] << " ";
//	std::cout << "\n\n";
//
//	for (int i=1; i <= strucTab[0]; ++i)
//		std::cout << strucTab[i] << " ";
//	std::cout << "\n\n";
//
//	for (int i=0; i < region.interactions.size(); ++i)
//			std::cout << region.interactions[i] << "," << strucTab[i+1]-1 << " ";
//		std::cout << "\n\n";

	for (int i=region.pos.first; i <= region.pos.second; ++i){
		int j = region.interactions[i];

		//if (j == -1)
		//	continue;

		if ((strucTab[i+1]-1) != j)
			return false;
	}

	return true;
}

void getHairpinCandidates(std::unordered_map<int, std::tuple<int,int,int> > &hairpins, std::vector<int> &hairpinKeys, vrna_fold_compound_t *vc, double threshold){
	double* probs = vc->exp_matrices->probs;
	int *iindex = vc->iindx;
	unsigned length = vc->length;

	// locate potential hairpins
	// iterate along the matrix diagonals
	// TODO: i = 0 happens at the start?
	for (unsigned k=1; k<length; k++){
		for (int offset=0; offset <= 1; ++offset){
			int i = k-1+offset;
			int j = k+1;

			bool inHairpin = false;

			std::tuple<int,int,int> hairpin;

			// scan the diagonal upwards
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
				hairpinKeys.push_back(k);
			}
		}
	}
}

TStemLoopProfile enforceHairpins(Motif &motif, std::unordered_map<int, std::tuple<int,int,int> > &hairpins, std::vector<int> &hairpinKeys, vrna_fold_compound_t *vc){
	bool skipEnforce = true;

	double* probs = vc->exp_matrices->probs;
	int *iindex = vc->iindx;
	unsigned length = vc->length;

	char *structure = (char*)vrna_alloc(sizeof(char) * (length + 1));
	vrna_mfe(vc, structure);

	TStemLoopProfile result_regions;
	TConsensusStructure consensusStructure;
	structureToInteractions(structure, consensusStructure);
	TStemLoopProfile stemLoops = findStemLoops(consensusStructure);

	if (skipEnforce){
		for (auto stemLoop : stemLoops){
			result_regions.push_back(stemLoop);
		}

		return result_regions;
	}

	vrna_plist_t *pl1, *pl2;
	bool first = true;

	// iteratively enforce most likely bases from hairpins that haven't been encountered yet
	do {
		// create new structure
		if (!first){
			int k = hairpinKeys.back();
			hairpinKeys.pop_back();

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
			//std::cout << "Enforcing " << maxi << " " << maxj << "\n";
			vrna_hc_add_bp(vc, maxi, maxj, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_ENFORCE);
			vrna_mfe(vc, structure);

			std::cout << "Vienna: " << structure << "\n";

			structureToInteractions(structure, consensusStructure);
			stemLoops = findStemLoops(consensusStructure);

			for (auto pair : stemLoops){
				std::cout << pair.pos.first << " " << pair.pos.second << "\n";
			}
		}

		first = false;
		std::vector<TConsensusStructure> structureVariants;

		// check with which hairpins our stemLoops overlap
		for (auto stemLoop : stemLoops){
			partitionStemLoop(motif.seedAlignment, stemLoop);

			StructureElement &hairpinStem = stemLoop.elements[stemLoop.elements.size()-2];
			int hairpinStart = hairpinStem.location;
			int hairpinEnd = hairpinStart + seqan::length(hairpinStem.loopComponents[0])
										  + seqan::length(stemLoop.elements.back().loopComponents[0])
										  + seqan::length(hairpinStem.loopComponents[1]);

			std::cout << hairpinStart << " -- " << hairpinEnd << "\n";

			bool stemAdded = false;
			std::unordered_map<int, std::tuple<int,int,int> > tmpHairpin = hairpins;
			for (auto hairpin : tmpHairpin){
				int k = hairpin.first;
				int i,j,l;
				std::tie(i,j,l) = hairpin.second;

				int posi = i-l;
				int posj = j+l;

				std::cout << "Testing " << i << " " << j << "\n";

				// check if the hairpin region overlaps with any stemloop
				for (int off=0; off <= l; ++off){
					// if overlap -> partition and save among the motifs
					if ((stemLoop.pos.first < posi+off) && (posj-off < stemLoop.pos.second)){
						//std::cout << stemLoop.pos.first << "-" << stemLoop.pos.second << " overlaps " << (posi+off) << "-" << (posj-off) << "\n";
						//std::cout << hairpinStart << " " << hairpinEnd << "-" << (posi+off) << " " << (posj-off) << "\n";
					}
					if ((hairpinStart <= (posi+off) && (posi+off) < hairpinEnd && hairpinEnd <= (posj-off)) ||
						((posi+off) <=hairpinStart && hairpinStart < (posj-off) && (posj-off) <= hairpinEnd)){
						//std::cout << hairpinStart << " " << hairpinEnd << " overlaps" << (posi+off) << "-" << (posj-off) << "\n";

						if (!stemAdded){
							stemAdded = true;
							std::cout << "Added loop " << stemLoop.pos.first << " " << stemLoop.pos.second << "\n";
							result_regions.push_back(stemLoop);

							//pl1 = vrna_plist_from_probs(vc, 0.005);
							//pl2 = vrna_plist(structure, 0.95*0.95);

							//char *tmp = (char*)("hairpin" + std::to_string(k) + std::string(".ps")).c_str();
							//(void) PS_dot_plot_list((char*)consens_mis((const char**)vc->sequences), tmp, pl1, pl2, structure);

						}

						//std::cout << "Hairpin " << k << " from: (" << i-l << "," << j+l << ") to (" << i << "," << j << ")\n";

						hairpins.erase(k);
						hairpinKeys.erase(std::remove(hairpinKeys.begin(), hairpinKeys.end(), k), hairpinKeys.end());

						break;
					}
				}
			}
		}

		// TODO: Actually enforce hairpins.
		//break;

	} while (!hairpins.empty() && !hairpinKeys.empty());

	free(structure);

	return result_regions;
}

void getConsensusStructure(Motif &motif, seqan::StockholmRecord<TBaseAlphabet> const & record, const char* constraint, RNALibFold const &){
//void getConsensusStructure(Motif &motif, seqan::StockholmRecord<TBaseAlphabet> const & record, TInteractionPairs &consensusStructure, const char* constraint, RNALibFold const &){
	TConsensusStructure consensusStructure;
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
	md.ribo = 1;
	vrna_fold_compound_t *vc 	= vrna_fold_compound_comparative((const char**)seqs, &md, VRNA_OPTION_MFE | VRNA_OPTION_PF);

	// add constraints if available
	if (constraint){
		vrna_constraints_add(vc, constraint, VRNA_CONSTRAINT_DB | VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_RND_BRACK);
		std::cout << "bracket:" << constraint << std::endl;
	}

	double min_en = vrna_mfe(vc, structure);

	/* rescale parameters for Boltzmann factors */
	//vrna_exp_params_rescale(vc, &min_en);

	float energy, kT, sfact;
	float betaScale = 1.;

	// scaling parameters
	kT = (betaScale*((temperature+K0)*GASCONST))/1000.; /* in Kcal */
	vrna_exp_param_t  *pf_parameters;
	sfact         = 1.07;

	pf_scale = std::exp(-(sfact*min_en)/kT/length);

    /* rescale energy parameters according to above calculated pf_scale */
    pf_parameters = vrna_exp_params_comparative(n_seq, &md);
    pf_parameters->pf_scale = pf_scale;

    /* change energy parameters in vc */
    vrna_exp_params_subst(vc, pf_parameters);

    energy = vrna_pf(vc, prob_structure);

	structureToInteractions(structure, consensusStructure);
	TStemLoopProfile stemLoops = findStemLoops(consensusStructure);

	for (auto pair : stemLoops){
		std::cout << pair.pos.first << " " << pair.pos.second << "\n";
	}

	int n = vc->length;

    std::unordered_map<unsigned, bool> seenStructures;

    if (stemLoops.empty()){
    	std::cout << "No stem loops found.\n";
    	return;
    }

	std::cout << "Vienna: " << structure << " Freq.: " << std::exp((energy-min_en)/kT) << " " << min_en << " " << energy << "\n";
	//std::cout << "Cent. : " << cent_structure << " Dist.: " << cdist << "\n";
	std::cout << "        " << consens_mis((const char**)seqs) << "\n";
//	std::cout << "Weird:  " << prob_structure << "\n";

	std::cout << "\n";

	std::unordered_map<int, std::tuple<int,int,int> > hairpins;
	std::vector<int> hairpinKeys;

	std::vector<FLT_OR_DBL> diagonals(2*n, -1);

	double threshold = 0.1;
	getHairpinCandidates(hairpins, hairpinKeys, vc, threshold);

	/*
	std::cout << "Candidates:\n";

	for (auto key : hairpinKeys){
		int i,j,l;
		std::tie(i,j,l) = hairpins[key];
		std::cout << i << " " << j << " - " << l << "\n";
	}
	*/

	TStemLoopProfile result_regions = enforceHairpins(motif, hairpins, hairpinKeys, vc);

	vrna_hc_init(vc);
	vrna_pf(vc, NULL);

	// assign probabilities to the hairpins we found
	// get dot plot structures
	vrna_plist_t *pl1, *pl2;
	//pl1 = vrna_plist_from_probs(vc, 0.005);
	//pl2 = vrna_plist(structure, 0.95*0.95);

	std::cout << "Regions (" << result_regions.size() << ")\n";
	for (auto& pair : result_regions){
		partitionStemLoop(motif.seedAlignment, pair);

		//pair.countProb /= boltzmann_samples;
		vrna_fold_compound_t *vc2 	= vrna_fold_compound_comparative((const char**)seqs, &md, VRNA_OPTION_MFE | VRNA_OPTION_PF);

		char* stemLoopStruc = interactionsToStructure(pair.interactions, pair.pos.first, pair.pos.second);

		std::cout << "       " << stemLoopStruc << "\n";

		// scaling parameters
		kT = (betaScale*((temperature+K0)*GASCONST))/1000.; /* in Kcal */
		vrna_exp_param_t  *pf_parameters;
		sfact         = 1.07;

		pf_scale = std::exp(-(sfact*min_en)/kT/length);

	    /* rescale energy parameters according to above calculated pf_scale */
	    pf_parameters = vrna_exp_params_comparative(n_seq, &md);
	    pf_parameters->pf_scale = pf_scale;

	    /* change energy parameters in vc */
	    vrna_exp_params_subst(vc2, pf_parameters);

		vrna_constraints_add(vc2, stemLoopStruc, VRNA_CONSTRAINT_DB | VRNA_CONSTRAINT_DB_DOT | VRNA_CONSTRAINT_DB_RND_BRACK | VRNA_CONSTRAINT_DB_ENFORCE_BP);

		double sub_energy = vrna_pf(vc2, prob_structure);
		std::cout << "       " << prob_structure << " " << sub_energy << "\n";

		pair.prob = std::exp((energy-sub_energy)/kT);

		std::cout << pair.pos.first << " " << pair.pos.second <<  " " << pair.prob << "\n";

		vrna_fold_compound_free(vc2);

		std::cout << "------ \n";
	}

	motif.profile = result_regions;

	//DEBUG_MSG("Vienna: " << structure);
	//DEBUG_MSG("        " << consens_mis((const char**)seqs));

	// write dot-plot
	//	Function used to plot the dot_plot graph
	//(void) PS_dot_plot_list((char*)seqs[0], (record.header.at("AC") + std::string(".ps")).c_str(), pl1, pl1, "");
	//std::cout << "Writing dot plot, no?\n";
	//char *tmp = (char*)(record.header.at("AC") + std::string(".ps")).c_str();
	//(void) PS_dot_plot_list((char*)consens_mis((const char**)seqs), tmp, pl1, pl2, structure);

	// free all used RNAlib data structures
	free(structure);
	free(prob_structure);
	vrna_fold_compound_free(vc);

	// free the sequence array (terminated with a 0 at the end)
	for (size_t k = 0; seqs[k] != 0; ++k)
		free(seqs[k]);

	free(seqs);}

#endif  // #ifndef APPS_RNAMOTIF_RNALIB_UTILS_H_
