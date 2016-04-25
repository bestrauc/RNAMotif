// ==========================================================================
//                                  motif.h
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

#ifndef APPS_RNAMOTIF_MOTIF_H_
#define APPS_RNAMOTIF_MOTIF_H_

#include "motif_structures.h"
#include <stack>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// take a structure table and determine:
// Stem Loop -> Stem, Hairpin
void structurePartition(Motif &motif){
	TInteractionPairs & consensus = motif.consensusStructure;
	std::pair<int, int > lastHairpin;
	int closing = 0;
	int open = 0;

	std::stack<int> pairStack;

	for (size_t i = 0; i < motif.consensusStructure.size(); ++i){
		int bracket = consensus[i];
		if (bracket == -1)
			continue;

		// opening bracket, save previous stem-loop if there was one
		if (bracket > i){
			// if we found a open/close match before, save that and reset
			if (lastHairpin.second > 0){
				motif.hairpinLoops.push_back(lastHairpin);

				// clear stack and reset last hairpin found
				std::stack<int>().swap(pairStack);
				lastHairpin = std::pair<int, int>();
			}

			pairStack.push(i);
		}
		// closing bracket
		else {
			if (!pairStack.empty() && bracket == pairStack.top()){
				lastHairpin.first = bracket;
				lastHairpin.second = i;
				pairStack.pop();
			}
		}
	}

	// save last stem-loop if there is one
	if (lastHairpin.second > 0)
		motif.hairpinLoops.push_back(lastHairpin);

	return;
}

#endif  // #ifndef APPS_RNAMOTIF_MOTIF_H_
