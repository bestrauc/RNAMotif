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

#ifndef APPS_RNAMOTIF_IPKNOT_UTILS_H_
#define APPS_RNAMOTIF_IPKNOT_UTILS_H_

#include "motif_structures.h"

// Import IPknot libraries
#include "ipknot/ip.h"
#include "ipknot/fa.h"
#include "ipknot/aln.h"
#include "ipknot/fold.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// reads the multiple alignment sequences into the IPknot Aln object
void getConsensusStructure(StockholmRecord<seqan::Rna> & record, TInteractionPairs &consensusStructure, const char* constraint, IPknotFold const & tag){
	std::list<std::string> names;
	std::list<std::string> seqs;
	for (auto elem : record.seqences){
		names.push_back(elem.first);
		seqs.push_back(elem.second);
	}

	Aln aln(names, seqs);
	std::vector<float> bp;
	std::vector<int> offset;
	std::vector<int> bpseq;
	std::vector<int> plevel;

	BPEngineAln* mix_en=NULL;
	std::vector<BPEngineSeq*> en_s;
	std::vector<BPEngineAln*> en_a;
	const char* param=NULL;

	BPEngineSeq* e = new CONTRAfoldModel();
	en_s.push_back(e);
	en_a.push_back(new AveragedModel(e));

	//BPEngineAln* en= mix_en ? mix_en : en_a[0];
	BPEngineAln* en = new MixtureModel(en_a);
	en->calculate_posterior(aln.seq(), bp, offset);
//	if (max_pmcc)
//	  ipknot.solve(aln->size(), bp, offset, ep, bpseq, plevel);
//	else
//	  ipknot.solve(aln->size(), bp, offset, t, bpseq, plevel);
//	for (int i=0; i!=n_refinement; ++i)
//	{
//	  update_bpm(pk_level, aln->seq(), *en, bpseq, plevel, bp, offset);
//	  if (max_pmcc)
//		ipknot.solve(aln->size(), bp, offset, ep, bpseq, plevel);
//	  else
//		ipknot.solve(aln->size(), bp, offset, t, bpseq, plevel);
//	}
//	if (os_bpseq)
//	  output_bpseq(*os_bpseq, aln->name().front(), aln->consensus(), bpseq, plevel);
//	if (os_mfa)
//	  output_mfa(*os_mfa, *aln, bpseq, plevel);
//	if (os_bpseq!=&std::cout && os_mfa!=&std::cout)
//	  output_fa(std::cout, aln->name().front(), aln->consensus(), bpseq, plevel);
}

#endif  // #ifndef APPS_RNAMOTIF_IPKNOT_UTILS_H_
