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

// Example: (((((((.......((((((((..(((((..(((((.....((((((....((((...))))))))))...((((....((.....))....))))))))).)))))....))).))))))))))))........((((((.....)))))).................
// take a structure table and determine the structural elements (stem, bulge, internal loop, hairpin)
void structurePartition(Motif &motif){
	TStemLoopRegions stemLoops;

	TInteractionPairs & consensus = motif.consensusStructure;
	std::pair<int, int > lastHairpin;
	std::stack<int> pairStack;

	// locate the stem loops
	for (size_t i = 0; i < motif.consensusStructure.size(); ++i){
		int bracket = consensus[i];

		// skip unmatched regions
		if (bracket == -1)
			continue;

		// opening bracket, save previous stem-loop if there was one
		if (bracket > i){
			// if we found a open/close match before, save that and reset
			if (lastHairpin.second > 0){
				stemLoops.push_back(lastHairpin);

				// clear stack and reset last hairpin found
				std::stack<int>().swap(pairStack);
				lastHairpin = std::pair<int, int>();
			}

			// save the opening bracket on the stack
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
		stemLoops.push_back(lastHairpin);

	int pos = -8;
	for (auto pair : stemLoops){
		// visualize stem loops for debugging
		std::cout << std::string(pair.first-pos, ' ');
		std::cout << std::string(pair.second-pair.first+1, '+');
		pos += (pair.first-pos) + (pair.second-pair.first+1);
	}

	std::cout << "\n";

	// after locating stem loops, separate structural elements:
	// * a run of opening brackets is the left side of a stem
	// * any unpaired bases in between are interior loops
	// 		* if no corresponding unpaired bases: bulge
	// * innermost unpaired bases are the hairpin
	for (auto pair : stemLoops){
		size_t i = pair.first;
		do {
			int pos = i;
			int right = consensus[pos];

			// count the row of opening brackets
			if (right > pos){					// an opening bracket
				while (right > pos){			// while it's a series of opening brackets
					// check if the corresponding closing bracket follows
					// or if there is a bulge on the right side
					if ((consensus[pos+1] > pos+1) && consensus[right-1] == -1){
						std::cout << pos << " " << consensus[pos+1] << "\n";
						// there is one run of unpaired bases from right+1
						int unpaired = right+1;
						while (consensus[unpaired]==-1) ++unpaired;
						std::cout << "Left bulge in [" << right+1 << "," << unpaired-1 << "]\n";
					}

					++pos;
					right = consensus[pos];
				}

				std::cout << "Stem: [" << i << "," << pos-1 << "] ; [" << consensus[pos-1] << "," << consensus[i] << "]\n";
			}
			// if unpaired, count the length of the loop and check for a bulge
			else if (right == -1){
				// get right border bracket of other half of loop (->(...(..)...)<=)
				int run = pos;
				int rb = consensus[pos-1];
				while (consensus[run] == -1) ++run;

				// get partner of end bracket ((...->(..)<=...))
				int lb = consensus[run];

				if (rb - lb == 1)
					std::cout << "Left bulge in [" << pos << "," << run-1 << "]\n";
				else if (rb == run){
					std::cout << "Hairpin in [" << pos << "," << run-1 << "]\n";
					break;
				}
				else
					std::cout << "Left loop: " << "[" << pos << "," << run-1 << "]" << " ; " << "Right loop: " << lb+1 << "," << rb-1 << "\n";

				pos = run;
				right = consensus[pos];
			}
			else
				pos = pos + 1;

			i = pos;

		} while (i <= pair.second);

	}

	return;
}

#endif  // #ifndef APPS_RNAMOTIF_MOTIF_H_
