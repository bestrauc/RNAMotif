/*
 * vienna_rna_funct.h
 *
 *  Created on: Feb 18, 2016
 *      Author: vitrusky8
 */

#ifndef INCLUDE_VIENNA_RNA_FUNCT_H_
#define INCLUDE_VIENNA_RNA_FUNCT_H_

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <string.h>

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function computeBppMatrix()
// ----------------------------------------------------------------------------
template <typename TOption, typename TSeq, typename TRnaStruct>
void computeBppMatrix(TOption const & options, TSeq const & seq, bool const & write_dot_plot, TRnaStruct & rnaSeq)
{
	char *structure = new char[seqan::length(seq) + 1];
	vrna_md_t md_p;

//  apply default model details
	set_model_details(&md_p);

//	get a vrna_fold_compound with MFE and PF DP matrices and default model details
	vrna_fold_compound_t *vc = vrna_fold_compound(seq, &md_p, VRNA_OPTION_MFE | VRNA_OPTION_PF);
	double gibbs = (double)vrna_pf(vc, structure);

	vrna_plist_t *pl1, *ptr;
	pl1 = vrna_plist_from_probs(vc, options.cmd.bppmThreshold);
	if(write_dot_plot)
	{
		vrna_plist_t *pl2;
		pl2= vrna_plist(structure, 0.95*0.95);
//	Function used to plot the dot_plot graph
		(void) PS_dot_plot_list(seq, "prova_dot_plot", pl1, pl2, "");
	}

//	get size of pl1
	unsigned size;
	for(size = 0, ptr = pl1; ptr->i; size++, ptr++);

	std::cout << "BPPM2 => " << size  << std::endl;
	std::cout << "BPPM2 => " << std::endl;

	std::cout << "size seq = " << length(seq) << std::endl;
//	std::cout << "size graph = " << length(interGraph) << std::endl;
//	std::cout << interGraph << std::endl;
	for(unsigned i=0; i<size;++i)
	{
		std::cout << i << "_"<< pl1[i].i <<":"<< pl1[i].j <<"|"<< pl1[i].p <<"|"<< pl1[i].type << "\t";
		addEdge(rnaSeq.bpProb.interGraph, rnaSeq.bpProb.uVertexVect[pl1[i].i-1], rnaSeq.bpProb.uVertexVect[pl1[i].j-1], pl1[i].p);
		addEdge(rnaSeq.bpProb.interGraphUpdated, rnaSeq.bpProb.dVertexVect[pl1[i].i-1], rnaSeq.bpProb.dVertexVect[pl1[i].j-1], pl1[i].p);
		addEdge(rnaSeq.bpProb.interGraphUpdated, rnaSeq.bpProb.dVertexVect[pl1[i].j-1], rnaSeq.bpProb.dVertexVect[pl1[i].i-1], pl1[i].p);
	}
	std::cout << std::endl;


	std::cout << seq << std::endl;
	std::cout << structure << "\tgibbs = " << gibbs << std::endl;

//	free memory occupied by vrna_fold_compound
	vrna_fold_compound_free(vc);
//	clean up
	free(structure);
	//	return 0;
}

#endif /* INCLUDE_VIENNA_RNA_FUNCT_H_ */
