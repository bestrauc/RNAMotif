// ==========================================================================
//                               motif_search.h
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

#ifndef APPS_RNAMOTIF_MOTIF_SEARCH_H_
#define APPS_RNAMOTIF_MOTIF_SEARCH_H_

#include "motif_structures.h"
#include "motif.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/* ------------------------------------------------------- */
/*!
 * @class ProfileCharIter
 *
 * @brief Abstract base class of the ProfileChar iterator.
 *
 *	Allows us to use different subclass iterators for different ProfileChar
 *	Alphabet types while searching. (Mainly for Nucleotides and BiNucleotides)
 */

class ProfileCharIter{
public:
	bool gapped = false;
	bool left = false;

	int charNum = 0;
	int charPos = 0;
	int gapsLeft = 0;

	virtual int getNextChar() = 0;
	virtual int lastChar() = 0;
	virtual bool atGap() = 0;
	virtual bool atEnd() = 0;
	virtual void setEnd() {};
	virtual ~ProfileCharIter() {};
};

class ProfileIterEmpty : public ProfileCharIter{
public:
	ProfileIterEmpty (){};

	int getNextChar(){
		return -1;
	}

	int lastChar(){
		return 0;
	}

	bool atGap(){
		return false;
	}

	// end if all chars that occurred were returned
	bool atEnd(){
		return true;
	}
};

/*!
 * @class ProfileCharIterImpl
 *
 * @brief Implement the ProfileCharIter for all types of ProfileChars.
 *
 *	Templated implementation for ProfileChar types of different alphabets.
 *	Return the next characters as for all as int for consistency.
 */

template <typename TProfileChar>
class ProfileCharIterImpl : public ProfileCharIter{
	typedef typename seqan::SourceValue<TProfileChar>::Type TProfileAlphabet;
	TProfileChar c;
	const uint64_t char_size = seqan::ValueSize<TProfileChar>::VALUE;

	int state = 0;
	std::vector<int> idx;
	int total;
	int threshold;

public:
	ProfileCharIterImpl (TProfileChar c, int charNum, int charPos, int gapsLeft, bool left, int threshold=1) : c(c), idx(char_size), total(seqan::totalCount(c)){
		this->charNum = charNum;
		this->charPos = charPos;
		this->left = left;
		this->gapsLeft = gapsLeft;
		this->threshold = threshold;

		// if no gaps are left, clear gap character
		if (gapsLeft == 0){
			c.count[char_size-2] = 0;
		}

		// fill index vector with indies [0,1,..,N-1]
		std::iota(idx.begin(), idx.end(), 0);

		// sort index by comparing the entries in tmp to get order
		std::sort(idx.begin(), idx.end(),
			[&] (int i1, int i2) {
				return (c.count[i1] > c.count[i2]);
			}
		);
	}

	int getNextChar(){
		if (atEnd())
			return -1;

		return TProfileAlphabet(idx[state++]);
	}

	int lastChar(){
		return TProfileAlphabet(idx[state-1]);
	}

	bool atGap(){
		//std::cout << idx[state] << " " <<  (char_size-2) << "\n";
		return !atEnd() && (idx[state] == (char_size-2));
	}

	// end if all chars that occurred were returned
	bool atEnd(){
		return ((state == char_size) || (c.count[idx[state]] <= threshold));
	}

	void setEnd(){
		state = char_size;
	}
};

class StructureIterator{
public:
	typedef ProfileCharIterImpl<TAlphabetProfile> TSinglePointer;
	typedef ProfileCharIterImpl<TBiAlphabetProfile> TPairPointer;
	typedef std::shared_ptr<ProfileCharIter> ProfilePointer;

	std::vector<StructureElement> structure_elements;
	int element;
	int elem_length;
	int pos;
	int intital_pos;
	ProfilePointer prof_ptr;

	//bool full_pattern = false;
	uint64_t single = 0, paired = 0;
	uint64_t count;

	std::stack<ProfilePointer> state;
	std::tuple<int, int, int> end;

	StructureIterator(std::vector<StructureElement> &structure_elements)
		: structure_elements(structure_elements), end(-1,-1,-1) {
		uint64_t sum = 1;
		if (true){
			for (StructureElement elem : structure_elements){
				for (int i=0; i < seqan::length(elem.loopComponents[0]); ++i){
					int nonzero = 0;

					if (elem.type != StructureType::STEM){
						TAlphabetProfile &pchar = elem.loopComponents[0][i];
						int count = std::count_if(std::begin(pchar.count), std::end(pchar.count), [&] (int x) {return (x > 0);});
						nonzero = count;
						std::cout << "No stem: " << count << "\n";
					}
					else {
						TBiAlphabetProfile &pchar = elem.stemProfile[i];
						int count = std::count_if(std::begin(pchar.count), std::end(pchar.count), [&] (int x) {return (x > 0);});
						nonzero = count;
						std::cout << "Stem   : " << count << "\n";
					}

					sum *= nonzero;
					//std::cout << sum <<"\n";
				}
			}

			std::cout << "Number of sequences: " << sum << "\n";
		}

		//this->structure_elements = structure_elements;
		this->element = structure_elements.size()-1;
		this->elem_length = seqan::length(structure_elements[element].loopComponents[0]);
		this->pos = 0;
		this->count = 0;

		auto &hairpin = structure_elements[element].loopComponents[0];
		this->intital_pos = structure_elements[element].location;

		prof_ptr = ProfilePointer(new TSinglePointer(hairpin[pos], 1, this->intital_pos,
								  elem_length - structure_elements[element].statistics[0].min_length, false));

		//std::cout << "INIT gaps " << elem_length << " " << structure_elements[element].statistics[0].min_length << " " << structure_elements[element].statistics[0].max_length << "\n";
	}

	std::tuple<int, int, int> get_next_char(){
		std::tuple<int, int, int> ret;

		int backtracked = 0;

		// get previous state (backtrack to valid state if necessary)
		while (prof_ptr->atEnd()) {
			// if no more characters can be generated
			if (state.size() == 0){
				//std::cout << "ENDING\n";
				return this->end;
			}

			//full_pattern = false;

			//std::cout << "Popping " << (typeid(*prof_ptr) == typeid(TSinglePointer)) << " from stack\n";
			if(typeid(*prof_ptr) == typeid(ProfileIterEmpty)){
				prof_ptr = state.top();
				state.pop();
				continue;
			}

			--pos;

			// check if the substructure has to be changed again
			if (pos < 0) {
				//std::cout << "Going up\n";
				++element;
				elem_length = seqan::length(structure_elements[element].loopComponents[0]);
				pos = elem_length - 1;
			}

			prof_ptr = state.top();
			state.pop();

			// backtrack 1 character for unidirectional char, 2 for bidirectional
			// do not backtrack if we popped a gap character (since the iterator will have ignored it)
			if (typeid(*prof_ptr) == typeid(TSinglePointer)){
				backtracked += 1;
			}
			else{
				backtracked += 2;
			}
		}

		// adjust pattern length and position temporarily if we have a gap here
		// ********************************************************************
		int gap_correct = (1 + (structure_elements[element].type == StructureType::STEM));
		// if gap characters, do not account for them in the size
		if (prof_ptr->atGap()){
			//std::cout << "Gap correct\n";
			//std::cout << prof_ptr->charNum << "\n";
			prof_ptr->gapped = true;
			prof_ptr->charNum -= gap_correct;
		}
		// remove gap adjustment if there currently is one
		else if (prof_ptr->gapped){
			//std::cout << "Undo gap correct\n";
			prof_ptr->gapped = false;
			prof_ptr->charNum += gap_correct;
		}
		// ********************************************************************

		// subtract a gap if the current character is a gap
		int new_gaps = prof_ptr->gapsLeft - prof_ptr->gapped;

		// advance to the next character in the active state
		int next_char_val = prof_ptr->getNextChar();

		// return one or two characters to search, depending on the substructure
		if (typeid(*prof_ptr) == typeid(TSinglePointer)){
			//std::cout << "Single\n";
			++single;

			if (structure_elements[element].loopLeft == false){
				ret = std::make_tuple(backtracked, -1, next_char_val);
			}
			else {
				ret = std::make_tuple(backtracked, next_char_val, -1);
			}
		}
		else{
			++paired;
			//std::cout << "Double\n";
			int lchar = next_char_val / AlphabetSize;
			int rchar = next_char_val % AlphabetSize;

			ret = std::make_tuple(backtracked, lchar, rchar);
		}

		state.push(prof_ptr); // save last character state

		// advance in the structure elements if possible
		if (pos < (elem_length-1)) {
			++pos;
		}
		// go to next substructure if necessary
		else if ((pos == (elem_length-1)) && (element > 0)){
			//std::cout << "Changing structures\n";
			--element;
			elem_length = seqan::length(structure_elements[element].loopComponents[0]);
			pos = 0;

			// set allowed number of gaps for this substructures
			new_gaps = elem_length - structure_elements[element].statistics[0].min_length;
			//std::cout << "Gaps " << new_gaps << "\n";
		}
		// element == 0, searched everything
		else{
			//++count;
			prof_ptr = ProfilePointer(new ProfileIterEmpty());
			return ret;
		/*
			if (count % 1000000 == 0){
				//std::cout << state.size() << " - " << count << " / " << sum << "\n";
				std::stack<ProfilePointer> stateCopy = state;

				while (!stateCopy.empty()){
					std::cout << (int)stateCopy.top()->atEnd();
					stateCopy.pop();
				}
				std::cout << "\n";
		*/
		}

		//state.push(prof_ptr); // save last character state

		// Update character pointer.
		// Special cases:
		// - keep track of leftward start position in the alignment
		if (structure_elements[element].type == StructureType::STEM) {

			prof_ptr = ProfilePointer(new TPairPointer(structure_elements[element].stemProfile[pos],
									  prof_ptr->charNum+2,
									  prof_ptr->charPos-1,
									  new_gaps,
									  false));
		}
		else {
			int newpos = prof_ptr->charPos - structure_elements[element].loopLeft;

			prof_ptr = ProfilePointer(new TSinglePointer(structure_elements[element].loopComponents[0][pos],
									  prof_ptr->charNum+1,
									  newpos,
									  new_gaps,
									  structure_elements[element].loopLeft));
		}

		return ret;
	}

	bool at_char_end(){
		return prof_ptr->atEnd();
	}

	int patLen(){
		//if (full_pattern)			return prof_ptr->charNum;

		return state.empty() ? 0 : state.top()->charNum;
	}

	int patPos(){
		//if (full_pattern)		return prof_ptr->charPos;

		return state.empty() ? this->intital_pos : state.top()->charPos;
	}

	// print the current pattern
	std::string printPattern(){
		std::stack<ProfilePointer> stateCopy = state;
		std::stringstream ss;

		ProfilePointer tmpptr;

		std::stack<char> rightPrint;

		while(!stateCopy.empty()){
			tmpptr = stateCopy.top();
			stateCopy.pop();

			int cur_char_val = tmpptr->lastChar();

			if (typeid(*tmpptr) == typeid(TSinglePointer)){
				char cur_print_char;

				if (cur_char_val > AlphabetSize -1){
					std::cout << "\nWrong: " << cur_char_val << "\n";
					throw(std::runtime_error("Can't be AlphabetSize?!"));
				}

				if (cur_char_val == AlphabetSize-1){
					continue;
					cur_print_char = '-';
					//std::cout << "Waaaa\n";
				}
				else{
					cur_print_char = TBaseAlphabet(cur_char_val);
				}

				if (!tmpptr->left)
					rightPrint.push(cur_print_char);
				else
					ss << cur_print_char;

				//std::cout << "'" << cur_char_val << "' ";
			}
			else{
				int lchar = cur_char_val / AlphabetSize;
				int rchar = cur_char_val % AlphabetSize;

				char lprint_char;
				char rprint_char;

				if (lchar > AlphabetSize-1 || rchar > AlphabetSize-1){
					std::cout << "\nWrong: " << cur_char_val << " " << lchar << " " << rchar << "\n";
					throw(std::runtime_error("Can't be AlphabetSize?!"));
				}

				if (lchar == AlphabetSize-1 || rchar == AlphabetSize-1){
					continue;
					lprint_char = '-';
					rprint_char = '-';
					//std::cout << "Wuuuu\n";
				}
				else{
					lprint_char = TBaseAlphabet(lchar);
					rprint_char = TBaseAlphabet(rchar);
				}

				ss << lprint_char;
				//std::cout << "'(" << lchar << "','" << rchar << "') ";
				rightPrint.push(rprint_char);
			}
		};

		while (!rightPrint.empty()){
			ss << rightPrint.top();
			rightPrint.pop();
		}

		//std::cout << "\n";

		return ss.str();
	}

	// do not extend the current word further
	void skip_char(){
		count++;

		if (state.size() == 0){
			//std::cout << "FULL PATTERN\n";
			return;
		}

		if(typeid(*prof_ptr) == typeid(ProfileIterEmpty)){
			prof_ptr = state.top();
			state.pop();
			return;
		}

		--pos;

		if (pos < 0) {
			++element;
			elem_length = seqan::length(structure_elements[element].loopComponents[0]);
			pos = elem_length - 1;
		}

		prof_ptr = state.top();
		state.pop();
	}
};

/* ------------------------------------------------------- */

template <typename TBidirectionalIndex>
class MotifIterator{
	typedef typename seqan::Iterator<TBidirectionalIndex, seqan::TopDown<seqan::ParentLinks<> > >::Type TIterator;
	typedef typename seqan::SAValue<TBidirectionalIndex>::Type THitPair;
	typedef seqan::String< THitPair > TOccurenceString;

	StructureIterator structure_iter;
	TIterator it;
	bool cont = true;

	// threshold: stop expanding when below this likelihood for the sequence
	unsigned min_match;

	bool setEnd(){
		return (cont = false);
	}

public:
	MotifIterator(TStructure &structure, TBidirectionalIndex &index, double min_match)
		: structure_iter(structure.elements), it(index), min_match(min_match){
	}

	int patternPos(){
		return this->structure_iter.patPos();
	}

	// next() returns true as long as the motif is not exhausted.
	// Only 'valid' matches are iterated: exclude those who do not match
	// or do not represent the family well (prob. below threshold)
	bool next(){
		if (!cont){
			return false;
		}

		bool extension = true;
		bool paired;
		int backtracked, lchar, rchar;

		// search until the pattern cannot be extended further
		// or the required seed length has been reached.
		while (extension && seqan::repLength(it) < min_match){
			// get the next characters to search for (either one or two, depending on the search direction)
			std::tuple<int, int, int> n_char = structure_iter.get_next_char();

			// check if no next pattern existed, iterator exhausted
			if (n_char == structure_iter.end){
				return setEnd();
			}

			// unpack next pattern and backtracking information
			std::tie(backtracked, lchar, rchar) = n_char;

			// check if the next character was the result
			// of backtracking into a different pattern
			// reset search iterator accordingly
			while (backtracked > 0){
				//std::cout << seqan::representative(it) << "\n";
				bool wentUp = goUp(it);

				//std::cout << seqan::representative(it) << "\n";

				// Assertion: We have to be able to backtrack, else
				// structure iterator and search iterator are not in sync.
				if (!wentUp){
					throw(std::runtime_error("Can't go the way we came?"));
				}
				--backtracked;
			}

			//std::cout << "-------------------------\n";

			// remember if we tried to extend both sides, which we later have to undo
			paired = (lchar != -1) && (rchar != -1);

			// search for the pattern provided
			// one-directional extension rightwards
			if (lchar == -1) {
				extension = seqan::goDown(it, rchar, seqan::Rev());
			}
			// one-directional extension leftwards
			else if (rchar == -1) {
				extension = seqan::goDown(it, lchar, seqan::Fwd());
			}
			// bi-directional extension
			else {
				bool wentLeft  = seqan::goDown(it, lchar, seqan::Fwd());
				bool wentRight = seqan::goDown(it, rchar, seqan::Rev());

				// successful if we went both left and right
				extension = wentLeft && wentRight;

				// if we just followed one direction but not the other - go back
				if (wentLeft ^ wentRight){
					seqan::goUp(it);
				}
			}

			//std::cout << seqan::representative(it) << "\n";
		}

		// check if the minimum length was satisfied
		// if yes, return true as a hit
		// else, try to extend a different pattern

		/*
		int bla = structure_iter.patLen();

		std::stack<StructureIterator::ProfilePointer> stateCopy = structure_iter.state;

		while (!stateCopy.empty()){
			std::cout << (int)stateCopy.top()->atEnd();
			stateCopy.pop();
		}

		std::cout << "\n";
		*/

		// reset the current search pattern being generated
		// to prepare for the next invocation of next()
		structure_iter.skip_char();

		// if we stopped because no next character was found,
		// we do not have to reset the search iterator
		if (!extension){
			return this->next();
			if (false){
				std::cout	<< "No pattern    ";
				if (lchar != -1)
					std::cout << TBaseAlphabet(lchar);

				std::cout << seqan::representative(it);
				if (rchar != -1)
					std::cout << TBaseAlphabet(rchar);
				std::cout << "\n";
			}

		}
		// else we stopped because the minimum length was reached
		else {
			//std::cout << "Found pattern (pos " << structure_iter.patPos() << ")" << seqan::representative(it) << " " << seqan::repLength(it) << "\n";

			goUp(it);
			if (paired){
				goUp(it);
			}
		}

		return true;
	}

	TOccurenceString getOccurrences(){
		// FIXME: getting occurrences directly via
		// occs = seqan::getOccurrences(it) doesn't work for some reason
		TOccurenceString occs;
		for (THitPair test: seqan::getOccurrences(it))
			seqan::appendValue(occs, test);

		return occs;
	}

	unsigned countOccurrences(){
		return seqan::countOccurrences(it);
	}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

//template <typename TBidirectionalIndex>
//void searchProfile(TBidirectionalIndex &index, TStructure &profile){
void searchProfile(seqan::StringSet<seqan::String<TBaseAlphabet> > seqs, TStructure profile){
	typedef seqan::String< TIndexPosType > TOccurenceString;

	TBidirectionalIndex index(seqs);
	MotifIterator<TBidirectionalIndex> iter(profile, index, 11);

	std::cout << "Searching for " << profile.pos.first << " " << profile.pos.second << "\n";
	while (iter.next()){
		TOccurenceString occs = iter.getOccurrences();
	//	//std::cout << iter.countOccurrences() << "\n";
	}
}

template <typename TBidirectionalIndex>
//std::vector<TProfileInterval> getStemloopPositions(TBidirectionalIndex &index, Motif &motif){
std::vector<TProfileInterval> getStemloopPositions(TBidirectionalIndex &index, Motif *motif){
	typedef typename seqan::SAValue<TBidirectionalIndex>::Type THitPair;
	typedef seqan::String< TIndexPosType > TOccurenceString;

	//std::vector<TOccurenceString> result(profile.size());

	// create one interval tree for each contig in the reference genome
	std::vector<TProfileInterval> intervals(seqan::countSequences(index));
	//std::vector<TProfileInterval> intervals;

	std::cout << "INIT\n";

	int id = 0;
	//int n = seqan::length(motif->seedAlignment);
	int n = seqan::length(seqan::row(motif->seedAlignment,0));
	std::cout << "ALN LEN: " << n << "\n";
	int stems = motif->profile.size();

	// tolerance for the seach window to account for novel inserts
	int const eps = 10;

	for (TStructure &structure : motif->profile){
		// start of the hairpin in the whole sequence
		std::cout <<  structure.pos.first << " " << structure.pos.second << ": " << motif->profile.size() << "\n";

		std::cout << motif->header.at("AC") << " is being searched..\n";

		//std::cout << "ITER\n";

		MotifIterator<TBidirectionalIndex> iter(structure, index, 15);

		unsigned inters = 0;
		unsigned more1  = 0;

		unsigned occ_sum   = 0;
		unsigned pat_count = 0;

		while (iter.next()){
			++pat_count;
			TOccurenceString occs = iter.getOccurrences();
			unsigned occ_count = iter.countOccurrences();
			occ_sum += occ_count;

			int pattern_pos = iter.patternPos();

			for (TIndexPosType pos : occs){
				TProfileInterval &interval = intervals[pos.i1];

				seqan::String<TProfileCargo> hits;

				// check if the stem occurs in an already existing match region
				findIntervals(hits, interval, pos.i2);
				inters += seqan::length(hits);

				// if this stem isn't located in an existing match region
				// create a new match region and store the stem loop there
				if (seqan::length(hits) == 0){
					//std::cout << loc << " " << n << " " << pos.i2 << " " << pos.i2-loc << " " << pos.i2 + (n-loc) << "\n";

					std::shared_ptr<std::vector<bool> > stemSet(new std::vector<bool>(stems));
					(*stemSet)[id] = true;
					// add an interval around the matched stem
					//std::cout << "Adding interval from " << (pos.i2 - (  pattern_pos) - eps) << " to " << (pos.i2 + (n-pattern_pos) + eps) << "\n";
					//std::cout << "(Match: " << pos.i2 << " " << pattern_pos << ")" << "\t" << (pos.i2 - pattern_pos) << " " << (pos.i2 + (n-pattern_pos)) << "\n";
					seqan::addInterval(interval, pos.i2 - (  pattern_pos) - eps,
												 pos.i2 + (n-pattern_pos) + eps, stemSet);
				}
				// else add to the list of hits
				else{
					if (seqan::length(hits) > 1) ++more1;

					for (unsigned i=0; i < seqan::length(hits); ++i){
						(*hits[i].cargo)[id] = true;
					}
				}
			}
		}

		std::cout << "Intervals: " << inters << " " << more1 << "\n";
		std::cout << "Hits: " << occ_sum << " in " << pat_count << " patterns.\n";

		for (unsigned i=0; i < seqan::countSequences(index); ++i)
			std::cout << "After " << id << "," << i << ": " << intervals[i].interval_counter << "\n";

		//std::cout << seqan::length(result[id]) << "\n";
		//seqan::sort(result[id]);
		++id;
	}

	return intervals;
}

void countHits(TProfileInterval positions){
	seqan::String<TProfileCargo> hits;
	seqan::getAllIntervals(hits, positions);

	for (unsigned i=0; i < seqan::length(hits); ++i) {
		// only report hits where all stems occurred in the region
		auto it = hits[i].cargo->begin();
		int sum = 0;

		while (it != hits[i].cargo->end()){
			sum += *it;
			++it;
		}

		//if (std::all_of(hits[i].cargo->begin(), hits[i].cargo->end(), [](bool v) { return v; }))
		if (sum > (hits[i].cargo->size()/2)){
			std::cout << hits[i].i1 << "    " << hits[i].i2 << "\n";
		}
	}
}


template <typename TStringType>
std::vector<seqan::Tuple<int, 3> > findFamilyMatches(seqan::StringSet<TStringType> &seqs, std::vector<Motif*> &motifs){
	std::vector<seqan::Tuple<int, 3> > results;

	TBidirectionalIndex index(seqs);

    for (Motif *motif : motifs){
    	std::cout << motif->header.at("ID") << "\n";
    	// find the locations of the motif matches
    	//std::cout << motif.header.at("AC") << "\n";
    	//std::vector<TProfileInterval> result = getStemloopPositions(index, *motif);
    	std::vector<TProfileInterval> result = getStemloopPositions(index, motif);

		// cluster results into areas (i.e. where hairpins of a given type cluster together)
    	//std::cout << result[1].interval_counter << "\n";
    	//countHits(result[1], motif.consensusStructure.size());
    	int id = 1;
    	for (TProfileInterval intervals : result){
    		std::cout << "rmark" << id++ << "\n";
    		countHits(intervals);
    	}
    }

	return results;
}

#endif  // #ifndef APPS_RNAMOTIF_MOTIF_SEARCH_H_
