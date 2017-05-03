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
	bool left = false;

	int charNum = 0;
	int charPos = 0;
	int gapsLeft = 0;
	int offset = 0;
	int gapsize = 0;

	virtual std::pair<int,int> getNextChar() = 0;
	virtual int nextLength() = 0;
	virtual int lastChar() = 0;
	virtual bool atGap() = 0;
	virtual bool atEnd() = 0;
	virtual void setEnd() {};
	virtual ~ProfileCharIter() {};
};

class ProfileIterEmpty : public ProfileCharIter{
public:
	ProfileIterEmpty (int offset, int charNum) {
		this->offset = offset;
		this->charNum = charNum;
	};

	std::pair<int,int> getNextChar(){
		return std::make_pair(-1,-1);
	}

	int nextLength(){
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
	const uint64_t char_size = seqan::ValueSize<TProfileChar>::VALUE;

	TProfileChar c;
	std::vector<std::pair<int, int>> gapVals;
	int gap_pos = 0;

	int state = 0;
	std::vector<int> idx;
	int total;
	int threshold;
	bool at_gap = false;


public:
	ProfileCharIterImpl (TProfileChar c, std::map<int, int> &gapMap, int charNum, int charPos, int gapsLeft, int offset, bool left, int threshold=0)
							: c(c), idx(char_size), total(seqan::totalCount(c)){
		this->charNum = charNum;
		this->charPos = charPos;
		this->left = left;
		this->gapsLeft = gapsLeft;
		this->offset = offset;
		this->threshold = threshold;

		// fill index vector with indies [0,1,..,N-1]
		std::iota(idx.begin(), idx.end(), 0);

		// sort index by comparing the entries in tmp to get order
		std::sort(idx.begin(), idx.end(),
			[&] (int i1, int i2) {
				return (c.count[i1] > c.count[i2]);
			}
		);

		int ind = 0;
		while (c.count[idx[ind]] > 0){
			//std::cout << c.count[idx[ind]] << " char seen " << idx[ind] << "\n";
			++ind;
		}

		//std::cout << "Gaps: " << gapsLeft << "\n";


		// insert gaps into gapVal vector and sort by value
		for (auto itr = gapMap.begin(); itr != gapMap.end(); ++itr){
			// do not include gaps that don't occur often enough or which would exceed
			// the maximum number of gaps for this section
			if (itr->second < threshold || itr->first > gapsLeft)
				continue;

			gapVals.push_back(*itr);
		}

		std::sort(gapVals.begin(), gapVals.end(),
			[&](std::pair<int, int>& a, std::pair<int, int>& b) {
				return a.second > b.second;
			}
		);
	}

	std::pair<int,int> getNextChar(){
		if (atEnd())
			return std::make_pair(-1,-1);

		at_gap = false;
		gapsize = 0;

		// return a gap
		if (state == char_size || (gap_pos < gapVals.size() && c.count[idx[state]] < gapVals[gap_pos].second)){
			at_gap = true;
			gapsize = gapVals[gap_pos].first;
			return std::make_pair(gapVals[gap_pos++].first, -1);
		}

		// else return an alphabet character
		return std::make_pair(0, idx[state++]);
	}

	// after this profile has been processed, the pattern length is increased
	// by 1 or, in the case of a bidirectional pattern, by 2
	int nextLength(){
		if (at_gap)
			return this->charNum;

		return (this->charNum + 1 + std::is_same<TProfileChar, TBiAlphabetProfile>::value);
	}

	int lastChar(){
		return TProfileAlphabet(idx[state-1]);
	}

	bool atGap(){
		return this->at_gap;
		//std::cout << idx[state] << " " <<  (char_size-2) << "\n";
		//return !atEnd() && (idx[state] == (char_size-2));
	}

	// end if all chars that occurred were returned
	bool atEnd(){
		return (state == char_size || (c.count[idx[state]] <= threshold)) && (gap_pos >= gapVals.size());
	}

	void setEnd(){
		state = char_size;
	}
};

class StructureIterator{
private:
	typedef ProfileCharIterImpl<TAlphabetProfile> TSinglePointer;
	typedef ProfileCharIterImpl<TBiAlphabetProfile> TPairPointer;
	typedef std::shared_ptr<ProfileCharIter> ProfilePointer;

	int total_length = 0;
	std::vector<int> cumulative_lens;
	std::vector<StructureElement> structure_elements;

	int element;
	int elem_length;
	int pos;
	int intital_pos;
	ProfilePointer prof_ptr;

	// return the profile pointer of the profile
	// that corresponds to position pos
	ProfilePointer next_profile(int offset){
		int next_pos  = prof_ptr->charPos;
		int last_gaps = prof_ptr->gapsLeft;

		bool changed = false;

		//std::cout << element << " " << pos << "\n";

		//std::cout << "Before " << last_gaps << " " << offset << "\n";

		for (int i=0; i < offset; ++i){
			// advance in the structure elements if possible
			if (pos < (elem_length-1)){
				++pos;
			}
			else if ((pos == (elem_length-1)) && (element > 0)){
				--element;
				elem_length = seqan::length(structure_elements[element].loopComponents[0]);
				last_gaps  = elem_length - structure_elements[element].statistics[0].min_length;
				pos = 0;
				//std::cout << "Changed!\n";
				changed = true;
			}
			// else return empty profile pointer as an end marker
			else {
				//std::cout << "AT ENEDEDEDD\n";
				return ProfilePointer(new ProfileIterEmpty(i, prof_ptr->nextLength()));
			}

			// the next position in the pattern goes to the left for stems or left-sided loops
			// else it does not change (-0)
			next_pos = next_pos - (structure_elements[element].type == StructureType::STEM || structure_elements[element].loopLeft);
		}

		if (!changed && prof_ptr->atGap()){
			last_gaps -= offset;
		}

		//std::cout << "After " << last_gaps << "\n";

		ProfilePointer next_ptr;

		if (structure_elements[element].type == StructureType::STEM) {

			next_ptr = ProfilePointer(new TPairPointer(structure_elements[element].stemProfile[pos],
									  structure_elements[element].gap_lengths[pos],
									  prof_ptr->nextLength(),
									  next_pos,
									  last_gaps,
									  offset,
									  false));
		}
		else {
			next_ptr = ProfilePointer(new TSinglePointer(structure_elements[element].loopComponents[0][pos],
									  structure_elements[element].gap_lengths[pos],
									  prof_ptr->nextLength(),
									  next_pos,
									  last_gaps,
									  offset,
									  structure_elements[element].loopLeft));
		}

		//std::cout << "new: " << next_ptr->charNum << "\n";

		//std::cout << element << " " << pos << "\n";

		return next_ptr;
	}

	// get next profile pointer position
	ProfilePointer next_profile(){
		return this->next_profile(1);
	}

	bool is_initial_state(){
		return state.size() == 0;
	}

	bool atEnd(){
		return (typeid(*prof_ptr) == typeid(ProfileIterEmpty));
	}

	std::pair<int, ProfilePointer> prev_profile(){
		//std::cout << "Removing " << (typeid(*prof_ptr) == typeid(TSinglePointer) ? "Single" : "Double") << " - " << prof_ptr->offset << " " << prof_ptr->charNum << "\n";
		std::pair<int, ProfilePointer> ret = {0, state.top()};
		//std::cout << "Next " << (typeid(*next_ptr) == typeid(TSinglePointer) ? "Single" : "Double") << " - " << next_ptr->offset << " " << next_ptr->charNum << "\n";

		//std::cout << "OFFSET: " << prof_ptr->offset << "\n";

		// decrement position in the motif
		//std::cout << "Resetting offset " << prof_ptr->offset << "\n";
		for (int i=0; i < prof_ptr->offset; ++i){
		//for (int i=0; i < ret.second->offset; ++i){
			--pos;

			if (pos < 0) {
				++element;
				elem_length = seqan::length(structure_elements[element].loopComponents[0]);
				pos = elem_length - 1;
			}
		}

		//std::cout << element << " " << pos << "\n";

		// backtrack 1 character for unidirectional char, 2 for bidirectional
		// do not backtrack if we popped a gap character (since the iterator will have ignored it)

		if (!ret.second->atGap()){
			ret.first += 1 + (typeid(*ret.second) == typeid(TPairPointer));
		}

		state.pop();

		return ret;
	}

public:
	//bool full_pattern = false;
	uint64_t single = 0, paired = 0;
	uint64_t count;

	std::stack<ProfilePointer> state;
	std::tuple<int, int, int> end;

	StructureIterator(std::vector<StructureElement> &structure_elements)
		: structure_elements(structure_elements), end(-1,-1,-1) {
		//for (StructureElement elem : structure_elements){
		//	cumulative_lens.push_back(total_length);
		//	total_length += seqan::length(elem.loopComponents[0]);
		//}

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

		prof_ptr = ProfilePointer(new TSinglePointer(hairpin[pos], structure_elements[element].gap_lengths[pos], 0,
								  this->intital_pos, this->elem_length-structure_elements[element].statistics[0].min_length,
								  1, false));

		//std::cout << "INIT gaps " << elem_length << " " << structure_elements[element].statistics[0].min_length << " " << structure_elements[element].statistics[0].max_length << "\n";
	}

	std::tuple<int, int, int> get_next_char(int backtrack=0){
		std::tuple<int, int, int> ret;

		int backtracked = backtrack;

		//std::cout << prof_ptr->atEnd() << "\n";

		// get previous state (backtrack to valid state if necessary)
		while (prof_ptr->atEnd()) {
			//std::cout << "At End\n";

//			// if no more characters can be generated
			if (is_initial_state()){
				return this->end;
			}

			//std::cout << "Backtracking\n";

			int tmp_back;
			std::tie(tmp_back, prof_ptr) = this->prev_profile();
			backtracked += tmp_back;

			//std::cout << backtracked << "\n";
		}

		//std::cout << prof_ptr->charNum << " num -- \n";

		// advance to the next character in the active state
		int skip_count, next_char_val;
		std::tie(skip_count, next_char_val) = prof_ptr->getNextChar();

		// save state of advanced pointer
		state.push(prof_ptr); // save last character state

		// if we are at a gap, skip the gap to the next character
		if (skip_count > 0){
			//std::cout << "Skipping gap of length " << skip_count << " " << prof_ptr->charNum << " " << prof_ptr->nextLength() << " " << state.size() << "\n";

			prof_ptr = this->next_profile(skip_count);
			return get_next_char(backtracked);
		}


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
			unsigned lchar = next_char_val / AlphabetSize;
			unsigned rchar = next_char_val % AlphabetSize;

			ret = std::make_tuple(backtracked, lchar, rchar);
		}

		prof_ptr = this->next_profile();

		return ret;
	}

	bool at_char_end(){
		return prof_ptr->atEnd();
	}

	int patLen(){
		//if (full_pattern)			return prof_ptr->charNum;

		//return state.empty() ? 0 : state.top()->charNum;
		return prof_ptr->charNum;
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

			//if (tmpptr->atGap())
			//	continue;

			int cur_char_val = tmpptr->lastChar();

			if (typeid(*tmpptr) == typeid(TSinglePointer)){
				char cur_print_char;

				if (cur_char_val > AlphabetSize -1){
					std::cout << "\nWrong: " << cur_char_val << "\n";
					throw(std::runtime_error("Can't be AlphabetSize?!"));
				}

				//if (cur_char_val == AlphabetSize-1){
				if (tmpptr->atGap()){
					//continue;
					for (int i=0; i < tmpptr->gapsize; ++i){
						//std::cout << "Offset: " << tmpptr->offset << "\n";
						if (tmpptr->left)
							ss << '-';
						else
							rightPrint.push('-');
					}
						//cur_print_char = '-';
					continue;
					//std::cout << "Waaaa\n";
				}
				//else{
				cur_print_char = TBaseAlphabet(cur_char_val);
				//}

				if (!tmpptr->left)
					rightPrint.push(cur_print_char);
				else
					ss << cur_print_char;

				//std::cout << "'" << cur_char_val << "' ";
			}
			else{
				//std::cout << cur_char_val << "--\n";
				unsigned lchar = cur_char_val / AlphabetSize;
				unsigned rchar = cur_char_val % AlphabetSize;

				char lprint_char;
				char rprint_char;

				if (lchar > AlphabetSize-1 || rchar > AlphabetSize-1){
					std::cout << "\nWrong: " << cur_char_val << " " << lchar << " " << rchar << "\n";
					throw(std::runtime_error("Can't be AlphabetSize?!"));
				}

				//if (lchar == AlphabetSize-1 || rchar == AlphabetSize-1){
				if (tmpptr->atGap()){
					//continue;
					for (int i=0; i < tmpptr->offset; ++i){
						ss << '-';
						rightPrint.push('-');
					}
					//lprint_char = '-';
					//rprint_char = '-';
					//std::cout << "Wuuuu\n";
					continue;
				}

				//else{
				lprint_char = TBaseAlphabet(lchar);
				rprint_char = TBaseAlphabet(rchar);
				//}

				//std::cout << lchar << "\n";

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
	int skip_char(){
		if (is_initial_state()){
			return -1;
		}

		//std::cout << "Skipping\n";

		count++;

		int backtracked;
		std::tie(backtracked, prof_ptr) = this->prev_profile();

		return backtracked;
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
	unsigned count = 0;

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
		int backtracked, lchar, rchar;

		do {
			// search until the pattern cannot be extended further
			// or the required seed length has been reached.
			while (extension && seqan::repLength(it) < min_match){
				// get the next characters to search for (either one or two, depending on the search direction)
				std::tuple<int, int, int> n_char = structure_iter.get_next_char();

				//std::string strString = structure_iter.printPattern();
				//std::cout << "Target: " << strString << " " << structure_iter.patPos() << " " << structure_iter.patLen() << "\n";
				//std::cout << "State:  " << seqan::representative(it) << "\n";

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
					//std::cout << "Backtracking.. " << seqan::representative(it) << " to ";
					//std::cout << seqan::representative(it) << "\n";
					bool wentUp = goUp(it);
					//bool wentUp = true;
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
				//paired = (lchar != -1) && (rchar != -1);

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

				if (structure_iter.patLen() >= min_match){
					extension = true;
					break;
				}

				//std::cout << " Lala\n";

				//std::cout << "After:  " << seqan::representative(it) << "\n";
			}

			// reset the current search pattern being generated
			// to prepare for the next invocation of next()
			/*
			std::string strString = structure_iter.printPattern();
			if (!extension){
				std::cout << strString << " not found. " << structure_iter.patPos() << " " << structure_iter.patLen() << "\n";
			}
			else{
				std::cout << strString << "     found. " << structure_iter.patPos() << " " << structure_iter.patLen() << "\n";
			}
			*/

		// repeat until a pattern of sucfficient length was searched for
		} while (!extension);

		count++;

		backtracked = structure_iter.skip_char();
		while(backtracked > 0){
			goUp(it);
			backtracked--;
		}

		//std::string strString = structure_iter.printPattern();
		//std::cout << "Searching for " << strString << " " << structure_iter.patPos() << " " << structure_iter.patLen() << "\n";
		//std::cout << "Iterator      " << seqan::representative(it) << "\n";

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

//		if (extension){
//			std::string strString = structure_iter.printPattern();
//			std::cout << "Searching for " << strString << " " << structure_iter.patPos() << " " << structure_iter.patLen() << "\n";
//		}
//
//		// reset the current search pattern being generated
//		// to prepare for the next invocation of next()
//		backtracked = structure_iter.skip_char();
//
//		// if we stopped because no next character was found,
//		// we do not have to reset the search iterator
//		if (!extension){
//			std::cout << "Recursive \n";
//			return this->next();
//			if (false){
//				std::cout	<< "No pattern    ";
//				if (lchar != -1)
//					std::cout << TBaseAlphabet(lchar);
//
//				std::cout << seqan::representative(it);
//				if (rchar != -1)
//					std::cout << TBaseAlphabet(rchar);
//				std::cout << "\n";
//			}
//
//		}
//		// else we stopped because the minimum length was reached
//		else {
//			count++;
//
//			while(backtracked > 0){
//				goUp(it);
//				backtracked--;
//			}
//			//std::cout << "Found pattern (pos " << structure_iter.patPos() << ")" << seqan::representative(it) << " " << seqan::repLength(it) << "\n";
//
//			//goUp(it);
//			//if (paired){
//			//	goUp(it);
//			//}
//		}

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
void searchProfile(seqan::StringSet<seqan::String<TBaseAlphabet> > seqs, TStructure profile, int threshold){
	typedef seqan::String< TIndexPosType > TOccurenceString;

	TBidirectionalIndex index(seqs);
	MotifIterator<TBidirectionalIndex> iter(profile, index, threshold);

	std::cout << "Searching for " << profile.pos.first << " " << profile.pos.second << "\n";

	if ((profile.pos.second - profile.pos.first + 1) < threshold){
		std::cout << "Skipping, length " << (profile.pos.second - profile.pos.first + 1) << " shorter than " << threshold << "\n";
		return;
	}

	while (iter.next()){
		TOccurenceString occs = iter.getOccurrences();
		std::cout << iter.countOccurrences() << "\n";
	}

	std::cout << iter.count << " counted \n";
}

template <typename TBidirectionalIndex>
//std::vector<TProfileInterval> getStemloopPositions(TBidirectionalIndex &index, Motif &motif){
std::vector<TProfileInterval> getStemloopPositions(TBidirectionalIndex &index, Motif *motif, int threshold){
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

		if ((structure.pos.second - structure.pos.first + 1) < threshold){
			std::cout << "Skipping, length " << (structure.pos.second - structure.pos.first + 1) << " shorter than " << threshold << "\n";
			continue;
		}

		//std::cout << "ITER\n";

		MotifIterator<TBidirectionalIndex> iter(structure, index, threshold);

		unsigned inters = 0;
		unsigned more1  = 0;

		unsigned occ_sum   = 0;
		unsigned pat_count = 0;

		while (iter.next()){
			//if (pat_count % 1000 == 0){
			//	std::cout << occ_sum << " " << pat_count << "\n";
			//}

			++pat_count;

			unsigned occ_count = iter.countOccurrences();
			continue;
			TOccurenceString occs = iter.getOccurrences();
			occ_sum += occ_count;

			int pattern_pos = iter.patternPos();

			std::cout << occ_count << "\n";

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

		std::cout << "Count " << pat_count << "\n";

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
std::vector<seqan::Tuple<int, 3> > findFamilyMatches(seqan::StringSet<TStringType> &seqs, std::vector<Motif*> &motifs, int threshold){
	std::vector<seqan::Tuple<int, 3> > results;

	TBidirectionalIndex index(seqs);

    for (Motif *motif : motifs){
    	std::cout << motif->header.at("ID") << "\n";
    	// find the locations of the motif matches
    	//std::cout << motif.header.at("AC") << "\n";
    	//std::vector<TProfileInterval> result = getStemloopPositions(index, *motif);
    	std::vector<TProfileInterval> result = getStemloopPositions(index, motif, threshold);

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
