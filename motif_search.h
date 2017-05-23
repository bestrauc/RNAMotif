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
	int endPos = 0;
	int gapsLeft = 0;
	int offset = 0;
	int gapsize = 0;
	THashType seqHash = 0;

	virtual std::pair<int,int> getNextChar() = 0;
	virtual int nextLength() = 0;
	virtual int lastChar() = 0;
	virtual bool atGap() = 0;
	virtual bool atEnd() = 0;
	virtual THashType nextHash() = 0;
	virtual void setEnd() {};
	virtual ~ProfileCharIter() {};
};

class ProfileIterEmpty : public ProfileCharIter{
public:
	ProfileIterEmpty (int offset, int charNum, THashType empty_hash) {
		this->offset = offset;
		this->charNum = charNum;
		this->seqHash = empty_hash;
	};

	std::pair<int,int> getNextChar(){
		return std::make_pair(-1,-1);
	}

	int nextLength(){
		return -1;
	}

	THashType nextHash(){
		return this->seqHash;
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
	int threshold;
	bool at_gap = false;


public:
	ProfileCharIterImpl (TProfileChar c, std::map<int, int> &gapMap, int charNum, int charPos, int gapsLeft, int offset, THashType prev_hash, bool left, double threshold_freq=0.05)
							: c(c), idx(char_size) {
		this->charNum = charNum;
		this->charPos = charPos;
		this->left = left;
		this->gapsLeft = gapsLeft;
		this->offset = offset;
		this->seqHash = prev_hash;

		int gap_count = std::accumulate(gapMap.begin(), gapMap.end(), 0,
										[](const int p, const std::pair<int,int>& n) {
											return p+n.second;
										});
		//this->threshold = (int)(seqan::totalCount(c) + gap_count)*threshold_freq;
		this->threshold = (int)seqan::totalCount(c)*threshold_freq;

		// fill index vector with indies [0,1,..,N-1]
		std::iota(idx.begin(), idx.end(), 0);

		// sort index by comparing the entries in tmp to get order
		std::sort(idx.begin(), idx.end(),
			[&] (int i1, int i2) {
				return (c.count[i1] > c.count[i2]);
			}
		);

		// insert gaps into gapVal vector and sort by value
		for (auto itr = gapMap.begin(); itr != gapMap.end(); ++itr){
			// do not include gaps that don't occur often enough or which would exceed
			// the maximum number of gaps for this section
			//if (itr->second < threshold || itr->first > gapsLeft)
			//	continue;

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

	THashType nextHash(){
		if (this->at_gap){
			return this->seqHash;
		}

		THashType next_hash = this->seqHash;

		int last_char = lastChar();
		int lchar, rchar;
		if (std::is_same<TProfileChar, TBiAlphabetProfile>::value){
			//std::cout << "Double\t";
			lchar = last_char / AlphabetSize;
			rchar = last_char % AlphabetSize;
		}
		else{
			//std::cout << "Single\t";
			lchar = this->left ? last_char : -1;
			rchar = this->left ? -1 : last_char;
		}

		if (lchar != -1){
			next_hash = (lchar << (2*this->charNum)) + next_hash;
		}

		if (rchar != -1){
			next_hash = next_hash*4 + rchar;
		}

		//std::cout << this->seqHash << " - " << lchar << " " << rchar << " - " << next_hash << "\n";

		return next_hash;
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

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };
}

class StructureIterator{
public:
	typedef ProfileCharIterImpl<TAlphabetProfile> TSinglePointer;
	typedef ProfileCharIterImpl<TBiAlphabetProfile> TPairPointer;
	typedef std::shared_ptr<ProfileCharIter> ProfilePointer;

	int total_length = 0;
	typedef std::pair<uint8_t, THashType> TPosHashPair;
	typedef std::vector<TPosHashPair> TPairArray;
	std::vector<std::vector<TPairArray> > prefix_states;
	std::unordered_set<TPosHashPair> end_state;

	std::vector<StructureElement> structure_elements;

	int element;
	int max_length;
	int elem_length;
	int pos;
	int intital_pos;
	int end_count = 0;
	double threshold_freq;
	ProfilePointer prof_ptr;
	bool hashLast;

	// return the profile pointer of the profile
	// that corresponds to position pos
	ProfilePointer next_profile(int offset, bool &duplicate){
		int next_pos  = prof_ptr->charPos;
		int last_gaps = prof_ptr->gapsLeft;

		bool changed = false;

		for (int i=0; i < offset; ++i){
			// advance in the structure elements if possible
			if (pos < (elem_length-1)){
				++pos;
			}
			else if ((pos == (elem_length-1)) && (element > 0)){
				--element;
				elem_length = seqan::length(structure_elements[element].loopComponents);
				last_gaps  = elem_length - structure_elements[element].statistics.min_length;
				pos = 0;
				//std::cout << "Changed! Now at " << element << " " << pos << "\n";
				changed = true;
			}
			// else return empty profile pointer as an end marker
			else {
				//std::cout << "AT ENEDEDEDD\n";
				duplicate = false;
				//std::cout << "End with " << prof_ptr->seqHash << " " << prof_ptr->nextHash() << "\n";
				return ProfilePointer(new ProfileIterEmpty(i, prof_ptr->nextLength(), prof_ptr->nextHash()));
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
									  prof_ptr->nextHash(),
									  false, threshold_freq));
		}
		else {
			next_ptr = ProfilePointer(new TSinglePointer(structure_elements[element].loopComponents[pos],
									  structure_elements[element].gap_lengths[pos],
									  prof_ptr->nextLength(),
									  next_pos,
									  last_gaps,
									  offset,
									  prof_ptr->nextHash(),
									  structure_elements[element].loopLeft, threshold_freq));
		}

		TPosHashPair hash_pair = std::make_pair(next_ptr->charNum, next_ptr->seqHash);

		size_t prefix_key = std::hash<TPosHashPair>()(hash_pair);
		if (prefix_states[element][pos][(prefix_key % HashTabLength)] != hash_pair){
			prefix_states[element][pos][prefix_key % HashTabLength] = hash_pair;
			duplicate = false;
		}
		else{
			duplicate = true;
			//std::cout << this->printPattern() << "\t   duplicate at " << element << " " << pos << " " << prefix_key << "\n";
		}


//		duplicate = false;

		return next_ptr;
	}

	// get next profile pointer position
	ProfilePointer next_profile(bool &duplicate){
		return this->next_profile(1, duplicate);
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

		// decrement position in the motif
		//std::cout << "Resetting offset " << prof_ptr->offset << "\n";
		for (int i=0; i < prof_ptr->offset; ++i){
		//for (int i=0; i < ret.second->offset; ++i){
			--pos;

			if (pos < 0) {
				++element;
				elem_length = seqan::length(structure_elements[element].loopComponents);
				pos = elem_length - 1;
			}
		}

		// backtrack 1 character for unidirectional char, 2 for bidirectional
		// do not backtrack if we popped a gap character (since the iterator will have ignored it)

		if (!ret.second->atGap()){
			ret.first += 1 + (typeid(*ret.second) == typeid(TPairPointer));
		}

		state.pop();

		return ret;
	}

//public:
	//bool full_pattern = false;
	uint64_t single = 0, paired = 0;
	uint64_t count;

	std::stack<ProfilePointer> state;
	std::tuple<int, int, int> end;

	StructureIterator(std::vector<StructureElement> &structure_elements, int length, bool hashLast, double threshold_freq=0.05)
		: structure_elements(structure_elements), max_length(length), threshold_freq(threshold_freq), hashLast(hashLast), end(-1,-1,-1) {

		uint64_t sum = 1;
		if (true){
			for (StructureElement elem : structure_elements){
				int elem_len = seqan::length(elem.loopComponents);
				//prefix_states.push_back(std::vector<std::unordered_set<std::pair<uint8_t, THashType> > >(elem_len));
				prefix_states.push_back(std::vector<TPairArray>(elem_len, TPairArray(HashTabLength, TPosHashPair(255,0))));

				for (int i=0; i < elem_len; ++i){
					//std::cout << i << " " << elem_len << " " << elem.type << "\n";
					int nonzero = 0;

					if (elem.type != StructureType::STEM){
						TAlphabetProfile &pchar = elem.loopComponents[i];
						//std::cout << seqan::totalCount(pchar) << " " << threshold_freq << " - " << seqan::totalCount(pchar)*threshold_freq << "\n";
						int count = std::count_if(std::begin(pchar.count), std::end(pchar.count), [&] (int x) {return (x > seqan::totalCount(pchar)*threshold_freq);});
						nonzero = count;
						std::cout << "No stem: " << count << "\n";
					}
					else {
						TBiAlphabetProfile &pchar = elem.stemProfile[i];
						//std::cout << seqan::totalCount(pchar) << " " << threshold_freq << " - " << seqan::totalCount(pchar)*threshold_freq << "\n";
						int count = std::count_if(std::begin(pchar.count), std::end(pchar.count), [&] (int x) {return (x > seqan::totalCount(pchar)*threshold_freq);});
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
		this->elem_length = seqan::length(structure_elements[element].loopComponents);
		this->pos = 0;
		this->count = 0;

		auto &hairpin = structure_elements[element].loopComponents;
		this->intital_pos = structure_elements[element].location;

		prof_ptr = ProfilePointer(new TSinglePointer(hairpin[pos], structure_elements[element].gap_lengths[pos], 0,
								  this->intital_pos, this->elem_length-structure_elements[element].statistics.min_length,
								  1, 0, false, threshold_freq));

		//std::cout << "INIT gaps " << elem_length << " " << structure_elements[element].statistics[0].min_length << " " << structure_elements[element].statistics[0].max_length << "\n";
	}

	std::tuple<int, int, int> make_return_char(int next_char_val, int backtracked, bool single_type){
		std::tuple<int, int, int> ret;

		// return one or two characters to search, depending on the substructure
		if (single_type){
			//std::cout << "Single\n";
			++single;

			if (structure_elements[element].loopLeft == false){
			//if (prof_ptr->left){
				ret = std::make_tuple(backtracked, -1, next_char_val);
			}
			else {
				ret = std::make_tuple(backtracked, next_char_val, -1);
			}
		}
		else{
			//std::cout << "Double\n";
			++paired;
			unsigned lchar = next_char_val / AlphabetSize;
			unsigned rchar = next_char_val % AlphabetSize;

			//std::cout << next_char_val << " ? " << lchar << " -- " << rchar << "\n";

			ret = std::make_tuple(backtracked, lchar, rchar);
		}

		return ret;
	}

	std::tuple<int, int, int> get_next_char(){
		std::tuple<int, int, int> ret;

		int next_char_val, backtracked = 0;
		bool single_type;
		bool duplicate = true;

		if (this->patLen() >= this->max_length){
			skip_char();
		}

		//std::cout << prof_ptr->atEnd() << "\n";

		while (duplicate){
			// get previous state (backtrack to valid state if necessary)
			while (prof_ptr->atEnd()) {
				//std::cout << "At End\n";

				// if no more characters can be generated
				if (is_initial_state()){
					return this->end;
				}

				//std::cout << "Backtracking\n";

				int tmp_back;
				std::tie(tmp_back, prof_ptr) = this->prev_profile();
				backtracked += tmp_back;

				//std::cout << backtracked << "\n";
			}

			// advance to the next character in the active state
			int skip_count;
			std::tie(skip_count, next_char_val) = prof_ptr->getNextChar();
			single_type = typeid(*prof_ptr) == typeid(TSinglePointer);

			//std::cout << skip_count << " skip, " << next_char_val << " - next char\n";

			// save state of advanced pointer
			state.push(prof_ptr); // save last character state

			// if we are at a gap, skip the gap to the next character
			if (skip_count > 0){
				//std::cout << "Skipping " << skip_count << " times at pos " << this->pos << "\n";
				prof_ptr = this->next_profile(skip_count, duplicate);

				// can't return the current character - it was a gap
				// get the actual next character after we skipped the gap
				duplicate = true;
				continue;
			}
			else {
				prof_ptr = this->next_profile(duplicate);
			}

			//duplicate = false;

			// if the prefix pattern was already seen before, abort and skip
			// the current extension. (Start the loop again top)
			if (duplicate){
				skip_char();
				continue;
			}

			// from here: only executed if !duplicate

			// check if the target length has been reached
			// if so, exclude duplicates again
			if (hashLast)
				if (this->patLen() >= this->max_length){
					TPosHashPair pair_hash = std::make_pair(this->patLen(), prof_ptr->seqHash);

					// skip
					if (end_state.find(pair_hash) != end_state.end()){
						skip_char();
						duplicate = true;
					}
					else{
						//end_state[pair_hash] = true;
						end_state.insert(pair_hash);
					}
				}
		}

		return make_return_char(next_char_val, backtracked, single_type);
	}

	bool at_char_end(){
		return prof_ptr->atEnd();
	}

	int patLen(){
		//if (full_pattern)			return prof_ptr->charNum;

		//return state.empty() ? 0 : state.top()->charNum;
		return prof_ptr->charNum;
	}

	uint64_t patHash(){
		//return prof_ptr->nextHash();
		return prof_ptr->seqHash;
	}

	std::string prevHash(){
		//if (full_pattern)			return prof_ptr->charNum;

		std::stringstream ss;
		auto stackCopy = state;
		while (!stackCopy.empty()){
			auto top = stackCopy.top();
			stackCopy.pop();
			ss << top->seqHash << " ";
		}

		return ss.str();
	}

	int patPos(){
		//if (full_pattern)		return prof_ptr->charPos;

		return state.empty() ? this->intital_pos : state.top()->charPos;
	}

	// print the current pattern
	std::string printPattern(bool printGaps = true){
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
					if (!printGaps)
						continue;

					if (tmpptr->left){
						ss << '|';
					}
					else{
						rightPrint.push('|');
					}

					for (int i=0; i < tmpptr->gapsize; ++i){
						//std::cout << "Offset: " << tmpptr->offset << "\n";
						if (tmpptr->left){
							ss << '-';
						}
						else{
							rightPrint.push('-');
						}
					}

					if (tmpptr->left){
						ss << '|';
					}
					else{
						rightPrint.push('|');
					}

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
					if (!printGaps)
						continue;

					ss << '|';
					rightPrint.push('|');

					for (int i=0; i < tmpptr->offset; ++i){
						ss << '-';
						rightPrint.push('-');
					}

					ss << '|';
					rightPrint.push('|');

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
	bool first = true;

	// threshold: stop expanding when below this likelihood for the sequence
	unsigned min_match;

	bool setEnd(){
		return (cont = false);
	}

public:
	unsigned count = 0;

	MotifIterator(TStructure &structure, TBidirectionalIndex &index, double min_match)
		: structure_iter(structure.elements, min_match, true), it(index), min_match(min_match){
	}

	int patternPos(){
		return this->structure_iter.patPos();
	}

	auto printRep(){
		return seqan::representative(it);
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

		if (!first){
			backtracked = structure_iter.skip_char();
			while(backtracked > 0){
				goUp(it);
				backtracked--;
			}
		}

		first = false;

		do {
			extension = true;

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

				//std::string strString = structure_iter.printPattern();
				//std::cout << "Target: " << strString << " " << structure_iter.patPos() << " " << structure_iter.patLen() << "\n";
				//std::cout << "State:  " << seqan::representative(it) << " | " << lchar << " " << rchar << " | " << backtracked << "\n";
				//std::cout << lchar << " " << rchar << backtracked << "\n";

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

				//if (structure_iter.patLen() >= min_match){
				////	extension = true;
				//	break;
				//}

				//std::cout << " Lala\n";

				//std::cout << "After:  " << seqan::representative(it) << "\n";
			}

			// reset the current search pattern being generated
			// to prepare for the next invocation of next()

			std::string strString = structure_iter.printPattern();
			if (!extension){
				//std::cout << strString << " not found. " << structure_iter.patPos() << " " << structure_iter.patLen() << "\n";
				structure_iter.skip_char();
			}
			else{
				//std::cout << strString << "     found. " << structure_iter.patPos() << " " << structure_iter.patLen() << "\n";
			}

		// repeat until a pattern of sucfficient length was searched for
		} while (!extension);

		count++;

		return true;
	}

	auto getOccurrences(){
		// FIXME: getting occurrences directly via
		// occs = seqan::getOccurrences(it) doesn't work for some reason
		//TOccurenceString occs = seqan::getOccurrences(it);
		//auto tests = seqan::getOccurrences(it);
		//for (THitPair test: seqan::getOccurrences(it))
			//seqan::appendValue(occs, test);

		return seqan::getOccurrences(it);
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
		auto occs = iter.getOccurrences();
		//std::cout << iter.countOccurrences() << "\n";
		//#pragma omp parallel for
		for (int i=0; i < seqan::length(occs); ++i){
			const TIndexPosType &pos = seqan::value(occs, i);
		}
	}

	std::cout << iter.count << " counted \n";
}

typedef std::vector<std::vector<std::vector<bool> > > TBoolVec;

template <typename TBidirectionalIndex>
//std::vector<TProfileInterval> getStemloopPositions(TBidirectionalIndex &index, Motif *motif, int threshold){
std::pair<TBoolVec,TBoolVec> getStemloopPositions(TBidirectionalIndex &index, Motif *motif, int threshold){
	typedef typename seqan::SAValue<TBidirectionalIndex>::Type THitPair;
	typedef seqan::String< TIndexPosType > TOccurenceString;

	int n_refs   = seqan::length(seqan::indexText(index));
	int n_strucs = motif->profile.size();

	//std::vector<TOccurenceString> result(profile.size());

	// create one interval tree for each contig in the reference genome
	std::vector<TProfileInterval> intervals(seqan::countSequences(index));
	TBoolVec pos_hits;
	TBoolVec pos_ends;
	//std::vector<TProfileInterval> intervals;

	std::cout << "INIT " << " " << n_refs << "\n";
	auto indText = seqan::indexText(index);


	for (int i=0; i < seqan::length(indText); ++i){
		int textLen = seqan::length(seqan::value(indText, i));
		pos_hits.push_back(std::vector<std::vector<bool> >(textLen, std::vector<bool>(n_strucs)));
	}

	int id = 0;
	//int n = seqan::length(motif->seedAlignment);
	int n = seqan::length(seqan::row(motif->seedAlignment,0));
	std::cout << "ALN LEN: " << n << "\n";
	int stems = motif->profile.size();

	// tolerance for the seach window to account for novel inserts
	int const eps = 5;

	for (unsigned i=0; i < motif->profile.size(); ++i){
		TStructure &structure = motif->profile[i];

		// start of the hairpin in the whole sequence
		std::cout <<  structure.pos.first << " " << structure.pos.second << ": " << motif->profile.size() << "\n";

		std::cout << motif->header.at("AC") << " is being searched..\n";

		int struclen = structure.pos.second - structure.pos.first + 1;

		/*
		if (struclen < threshold){
			std::cout << "Skipping, length " << (structure.pos.second - structure.pos.first + 1) << " shorter than " << threshold << "\n";
			++id;
			continue;
		}*/

		//std::cout << "ITER\n";

		MotifIterator<TBidirectionalIndex> iter(structure, index, std::min(struclen, threshold));
		//MotifIterator<TBidirectionalIndex> iter(structure, index, threshold);

		unsigned occ_sum   = 0;
		unsigned pat_count = 0;

		while (iter.next()){
			//if (pat_count % 1000 == 0){
			//	std::cout << occ_sum << " " << pat_count << "\n";
			//}

			++pat_count;

			unsigned occ_count = iter.countOccurrences();
			auto occs = iter.getOccurrences();
			occ_sum += occ_count;

			int pattern_pos = iter.patternPos();

			for (int i=0; i < seqan::length(occs); ++i){
				const TIndexPosType &pos = seqan::value(occs, i);

				// count a match at this position
				int stem_loop_pos = pattern_pos - structure.pos.first;
				int index = (pos.i2 >= stem_loop_pos) ? (pos.i2 - stem_loop_pos) : 0;
				pos_hits[pos.i1][index][id] = true;
				//int index = (pos.i2 >= pattern_pos) ? (pos.i2 - pattern_pos) : 0;

				/*
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
					const long unsigned int left  = (pos.i2 < (pattern_pos + eps)) ? 0 : (pos.i2 - (  pattern_pos) - eps);
					const long unsigned int right = pos.i2 + (n-pattern_pos) + eps;

					//std::cout << left << " (" << pos.i2 << " " << pos.i2 + n-pattern_pos << " )" << " " << right << "\n";

					//seqan::addInterval(interval, pos.i2 - (  pattern_pos) - eps, pos.i2 + (n-pattern_pos) + eps, stemSet);
					seqan::addInterval(interval, left, right, stemSet);
				}
				// else add to the list of hits
				else{
					if (seqan::length(hits) > 1) ++more1;

					for (unsigned i=0; i < seqan::length(hits); ++i){
						(*hits[i].cargo)[id] = true;
					}
				}
				*/
			}
		}

		//std::cout << "Count " << pat_count << " - " << occ_sum << " / " << (double)occ_sum / (double)pat_count << "\n";
		/*

		std::cout << "Intervals: " << inters << " " << more1 << "\n";
		std::cout << "Hits: " << occ_sum << " in " << pat_count << " patterns.\n";

		for (unsigned i=0; i < seqan::countSequences(index); ++i)
			std::cout << "After " << id << "," << i << ": " << intervals[i].interval_counter << "\n";
		*/

		std::cout << occ_sum << " found for " << id << "\n";

		//for (unsigned i=0; i < seqan::countSequences(index); ++i)
		//	std::cout << "After " << id << "," << i << ": " << std::accumulate(pos_hits[i].begin(), pos_hits[i].end(), 0) << "\n";
		++id;
	}

	return std::make_pair(pos_hits, pos_ends);
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

template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
    os<<"{ ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, " "));
    os<<"}";
    return os;
}

// scan the hit count vector with a sliding window

void countHits(Motif *motif, std::vector<std::vector<bool> > hits, int aln_len){
	int hitsize = hits[0].size();
	int n = hits.size();
	std::vector<int> state(hitsize);

	std::cout << hits.size() << " " << (aln_len+1) << " " << ((int)hits.size() - aln_len +1) << "--- \n";

	int hitThreshold = hitsize*0.8;

	std::cout << hitThreshold << " threshold\n";

	int window_start = -1;

	std::cout << "Counting hits?\n";

	// scan over the hits with a window of length aln_len
	for (int i=0; i < n; ++i){
		std::cout << hits[i] << "\n";
		std::cout << state << " " << i << "\n";

		// remove the contribution of the position that dropped from the window
		if (i >= aln_len){
			int hitsum = std::accumulate(state.begin(), state.end(), 0, [](int a, int b){return a + (b > 0);});

			// first check whether the current window had sufficient matches
			//if (hitsum > hitsize/2){
			if (hitsum > hitThreshold){
				//std::cout << state << " - " << i << "\n";
				if (window_start == -1){
					std::cout << "Window started at " << i - aln_len << "\n";
					window_start = i-aln_len;
				}

			}
			else if (window_start != -1){
				std::cout << "Windows ended at " << (i-1) << "\n";
				std::cout << window_start << " " << (i-1) << "\n";
				window_start = -1;
			}

			for (int j=0; j < hitsize; ++j){
				state[j] -= hits[i-aln_len][j];
			}
		}

		// update matches in the window in a sliding fashion
		for (int j=0; j < hitsize; ++j){
			state[j] += hits[i][j];
		}
	}

	int hitsum = std::accumulate(state.begin(), state.end(), 0, [](int a, int b){return a + (b > 0);});
	//std::cout << hitsum << "\n";
	// first check whether the current window had sufficient matches
	if (hitsum > hitThreshold){
		std::cout << std::max(0, n-aln_len) << " " << n-1 << "\n";
	}
}

void countHits2(Motif *motif, std::vector<std::vector<bool> > hits, int aln_len){
	int hitsize = hits[0].size();
	int n = hits.size();
	std::vector<int> state(hitsize);

	std::cout << hits.size() << " " << (aln_len+1) << " " << ((int)hits.size() - aln_len +1) << "--- \n";

	int hitThreshold = hitsize*0.8;
	std::cout << hitThreshold << " threshold\n";

	int window_start = -1;
	int hit_state = 0;

	std::cout << "Counting hits?\n";

	// scan over the hits with a window of length aln_len
	// as soon as the stem loops were found in order, report the window
	/*
	for (int i=0; i < n; ++i){
		if (hits[i][hit_state]){
			if (hit_state == 0){
				std::cout << "Start at " << motif-> << "\n";
			}
		}
	}
	*/

	//int hitsum = std::accumulate(state.begin(), state.end(), 0, [](int a, int b){return a + (b > 0);});
	//std::cout << hitsum << "\n";
	// first check whether the current window had sufficient matches
	//if (hitsum > hitThreshold){
	//	std::cout << std::max(0, n-aln_len) << " " << n-1 << "\n";
	//}
}

template <typename TStringType>
std::vector<seqan::Tuple<int, 3> > findFamilyMatches(seqan::StringSet<TStringType> &seqs, std::vector<Motif*> &motifs, int threshold){
	std::vector<seqan::Tuple<int, 3> > results;

	TBidirectionalIndex index(seqs);

	#pragma omp parallel for schedule(dynamic)
    for (unsigned i=0; i < motifs.size(); ++i){
    	Motif *motif = motifs[i];

    	if (motif == 0)
    		continue;

    	std::cout << motif->header.at("ID") << "\n";

    	//std::cout << motif->seedAlignment << "\n";

    	// find the locations of the motif matches
    	//std::cout << motif.header.at("AC") << "\n";
    	//std::vector<TProfileInterval> result = getStemloopPositions(index, motif, threshold);
    	int aln_len = seqan::length(seqan::row(motif->seedAlignment,0));
    	std::pair<TBoolVec, TBoolVec> result = getStemloopPositions(index, motif, threshold);

		// cluster results into areas (i.e. where hairpins of a given type cluster together)
    	//std::cout << result[1].interval_counter << "\n";
    	//countHits(result[1], motif.consensusStructure.size());
    	int id = 1;
    	for (auto region : result.first){
    		std::cout << "rmark" << id++ << "\n";
    		countHits(motif, region, aln_len);
    	}
    }

	return results;
}

#endif  // #ifndef APPS_RNAMOTIF_MOTIF_SEARCH_H_
