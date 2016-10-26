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

// Header to access the IPknot RNA structure prediction code.
// Code adapted from ipknot's main.cpp to work with seqan alignments.

#ifndef APPS_RNAMOTIF_IPKNOT_UTILS_H_
#define APPS_RNAMOTIF_IPKNOT_UTILS_H_

#include "motif_structures.h"

// Import IPknot libraries
#include "ipknot/ip.h"
#include "ipknot/aln.h"
#include "ipknot/fold.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/* IPknot functions */

typedef unsigned int uint;
typedef std::vector<float> VF;
typedef std::vector<VF> VVF;
typedef std::vector<int> VI;
typedef std::vector<VI> VVI;
typedef std::vector<VVI> VVVI;

class IPknot
{
public:
  template < class T > class EnumParam;

public:
  IPknot(uint pk_level, const float* alpha,
         bool levelwise, bool stacking_constraints, int n_th)
    : pk_level_(pk_level),
      alpha_(alpha, alpha+pk_level_),
      levelwise_(levelwise),
      stacking_constraints_(stacking_constraints),
      n_th_(n_th)
  {
  }

  void solve(uint L, const std::vector<float>& bp, const std::vector<int>& offset,
             const std::vector<float>& th, std::vector<int>& bpseq, std::vector<int>& plevel) const
  {
    IP ip(IP::MAX, n_th_);
    VVVI v(pk_level_, VVI(L, VI(L, -1)));
    VVVI w(pk_level_, VVI(L));

    // make objective variables with their weights
    for (uint j=1; j!=L; ++j)
    {
      for (uint i=j-1; i!=-1u; --i)
      {
        const float& p=bp[offset[i+1]+(j+1)];
        for (uint lv=0; lv!=pk_level_; ++lv)
          if (p>th[lv])
          {
            v[lv][i][j] = ip.make_variable(p*alpha_[lv]);
            w[lv][i].push_back(j);
          }
      }
    }

    ip.update();

    // constraint 1: each s_i is paired with at most one base
    for (uint i=0; i!=L; ++i)
    {
      int row = ip.make_constraint(IP::UP, 0, 1);
      for (uint lv=0; lv!=pk_level_; ++lv)
      {
        for (uint j=0; j<i; ++j)
          if (v[lv][j][i]>=0)
            ip.add_constraint(row, v[lv][j][i], 1);
        for (uint j=i+1; j<L; ++j)
          if (v[lv][i][j]>=0)
            ip.add_constraint(row, v[lv][i][j], 1);
      }
    }

    if (levelwise_)
    {
      // constraint 2: disallow pseudoknots in x[lv]
      for (uint lv=0; lv!=pk_level_; ++lv)
        for (uint i=0; i<w[lv].size(); ++i)
          for (uint p=0; p<w[lv][i].size(); ++p)
          {
            uint j=w[lv][i][p];
            for (uint k=i+1; k<j; ++k)
              for (uint q=0; q<w[lv][k].size(); ++q)
              {
                uint l=w[lv][k][q];
                if (j<l)
                {
                  int row = ip.make_constraint(IP::UP, 0, 1);
                  ip.add_constraint(row, v[lv][i][j], 1);
                  ip.add_constraint(row, v[lv][k][l], 1);
                }
              }
          }

      // constraint 3: any x[t]_kl must be pseudoknotted with x[u]_ij for t>u
      for (uint lv=1; lv!=pk_level_; ++lv)
        for (uint k=0; k<w[lv].size(); ++k)
          for (uint q=0; q<w[lv][k].size(); ++q)
          {
            uint l=w[lv][k][q];
            for (uint plv=0; plv!=lv; ++plv)
            {
              int row = ip.make_constraint(IP::LO, 0, 0);
              ip.add_constraint(row, v[lv][k][l], -1);
              for (uint i=0; i<k; ++i)
                for (uint p=0; p<w[plv][i].size(); ++p)
                {
                  uint j=w[plv][i][p];
                  if (k<j && j<l)
                    ip.add_constraint(row, v[plv][i][j], 1);
                }
              for (uint i=k+1; i<l; ++i)
                for (uint p=0; p<w[plv][i].size(); ++p)
                {
                  uint j=w[plv][i][p];
                  if (l<j)
                    ip.add_constraint(row, v[plv][i][j], 1);
                }
            }
          }
    }

    if (stacking_constraints_)
    {
      for (uint lv=0; lv!=pk_level_; ++lv)
      {
        // upstream
        for (uint i=0; i<L; ++i)
        {
          int row = ip.make_constraint(IP::LO, 0, 0);
          for (uint j=0; j<i; ++j)
            if (v[lv][j][i]>=0)
              ip.add_constraint(row, v[lv][j][i], -1);
          if (i>0)
            for (uint j=0; j<i-1; ++j)
              if (v[lv][j][i-1]>=0)
                ip.add_constraint(row, v[lv][j][i-1], 1);
          if (i+1<L)
            for (uint j=0; j<i+1; ++j)
              if (v[lv][j][i+1]>=0)
                ip.add_constraint(row, v[lv][j][i+1], 1);
        }

        // downstream
        for (uint i=0; i<L; ++i)
        {
          int row = ip.make_constraint(IP::LO, 0, 0);
          for (uint j=i+1; j<L; ++j)
            if (v[lv][i][j]>=0)
              ip.add_constraint(row, v[lv][i][j], -1);
          if (i>0)
            for (uint j=i; j<L; ++j)
              if (v[lv][i-1][j]>=0)
                ip.add_constraint(row, v[lv][i-1][j], 1);
          if (i+1<L)
            for (uint j=i+2; j<L; ++j)
              if (v[lv][i+1][j]>=0)
                ip.add_constraint(row, v[lv][i+1][j], 1);
        }
      }
    }

    // execute optimization
    ip.solve();

    // build the result
    bpseq.resize(L);
    std::fill(bpseq.begin(), bpseq.end(), -1);
    plevel.resize(L);
    std::fill(plevel.begin(), plevel.end(), -1);
    for (uint lv=0; lv!=pk_level_; ++lv)
    {
      for (uint i=0; i<L; ++i)
        for (uint j=i+1; j<L; ++j)
          if (v[lv][i][j]>=0 && ip.get_value(v[lv][i][j])>0.5)
          {
            bpseq[i]=j; bpseq[j]=i;
            plevel[i]=plevel[j]=lv;
          }
    }

    if (!levelwise_)
      decompose_plevel(bpseq, plevel);
  }

public:
  template < class T >
  class EnumParam
  {
  public:
    EnumParam(const std::vector<std::vector<T> >& p)
      : p_(p), m_(p.size()), v_(p.size(), 0)
    {
      for (uint i=0; i!=p.size(); ++i)
        m_[i] = p[i].size();
    }

    uint size() const { return m_.size(); }

    void get(std::vector<T>& q) const
    {
      for (uint i=0; i!=v_.size(); ++i)
        q[i] = p_[i][v_[i]];
    }

    bool succ()
    {
      return succ(m_.size(), &m_[0], &v_[0]);
    }

  private:
    static bool succ(int n, const int* m, int* v)
    {
      if (n==0) return true;
      if (++(*v)==*m)
      {
        *v=0;
        return succ(n-1, ++m, ++v);
      }
      return false;
    }

  private:
    const std::vector<std::vector<T> >& p_;
    std::vector<int> m_;
    std::vector<int> v_;
  };

private:
  struct cmp_by_degree : public std::less<int>
  {
    cmp_by_degree(const std::vector< std::vector<int> >& g) : g_(g) {}
    bool operator()(int x, int y) const { return g_[y].size()<g_[x].size(); }
    const std::vector< std::vector<int> >& g_;
  };

  struct cmp_by_count : public std::less<int>
  {
    cmp_by_count(const std::vector<int>& count) : count_(count) { }
    bool operator()(int x, int y) const { return count_[y]<count_[x]; }
    const std::vector<int>& count_;
  };

  static void
  decompose_plevel(const std::vector<int>& bpseq, std::vector<int>& plevel)
  {
    // resolve the symbol of parenthsis by the graph coloring problem
    uint L=bpseq.size();

    // make an adjacent graph, in which pseudoknotted base-pairs are connected.
    std::vector< std::vector<int> > g(L);
    for (uint i=0; i!=L; ++i)
    {
      if (bpseq[i]<0 || bpseq[i]<=(int)i) continue;
      uint j=bpseq[i];
      for (uint k=i+1; k!=L; ++k)
      {
        uint l=bpseq[k];
        if (bpseq[k]<0 || bpseq[k]<=(int)k) continue;
        if (k<j && j<l)
        {
          g[i].push_back(k);
          g[k].push_back(i);
        }
      }
    }
    // vertices are indexed by the position of the left base
    std::vector<int> v;
    for (uint i=0; i!=bpseq.size(); ++i)
      if (bpseq[i]>=0 && (int)i<bpseq[i])
        v.push_back(i);
    // sort vertices by degree
    std::sort(v.begin(), v.end(), cmp_by_degree(g));

    // determine colors
    std::vector<int> c(L, -1);
    int max_color=0;
    for (uint i=0; i!=v.size(); ++i)
    {
      // find the smallest color that is unused
      std::vector<int> used;
      for (uint j=0; j!=g[v[i]].size(); ++j)
        if (c[g[v[i]][j]]>=0) used.push_back(c[g[v[i]][j]]);
      std::sort(used.begin(), used.end());
      used.erase(std::unique(used.begin(), used.end()), used.end());
      int j=0;
      for (j=0; j!=(int)used.size(); ++j)
        if (used[j]!=j) break;
      c[v[i]]=j;
      max_color=std::max(max_color, j);
    }

    // renumber colors in decentant order by the number of base-pairs for each color
    std::vector<int> count(max_color+1, 0);
    for (uint i=0; i!=c.size(); ++i)
      if (c[i]>=0) count[c[i]]++;
    std::vector<int> idx(count.size());
    for (uint i=0; i!=idx.size(); ++i) idx[i]=i;
    sort(idx.begin(), idx.end(), cmp_by_count(count));
    std::vector<int> rev(idx.size());
    for (uint i=0; i!=rev.size(); ++i) rev[idx[i]]=i;
    plevel.resize(L);
    for (uint i=0; i!=c.size(); ++i)
      plevel[i]= c[i]>=0 ? rev[c[i]] : -1;
  }

  static void
  compute_expected_accuracy(const std::vector<int>& bpseq,
                            const std::vector<float>& bp, const std::vector<int>& offset,
                            float& sen, float& ppv, float& mcc)
  {
    int L  = bpseq.size();
    int L2 = L*(L-1)/2;
    int N = 0;

    float sump = 0.0;
    float etp  = 0.0;

    for (uint i=0; i!=bp.size(); ++i) sump += bp[i];

    for (uint i=0; i!=bpseq.size(); ++i)
    {
      if (bpseq[i]!=-1 && bpseq[i]>(int)i)
      {
        etp += bp[offset[i+1]+bpseq[i]+1];
        N++;
      }
    }

    float etn = L2 - N - sump + etp;
    float efp = N - etp;
    float efn = sump - etp;

    sen = ppv = mcc = 0;
    if (etp+efn!=0) sen = etp / (etp + efn);
    if (etp+efp!=0) ppv = etp / (etp + efp);
    if (etp+efp!=0 && etp+efn!=0 && etn+efp!=0 && etn+efn!=0)
      mcc = (etp*etn-efp*efn) / std::sqrt((etp+efp)*(etp+efn)*(etn+efp)*(etn+efn));
  }

  static uint length(const std::string& seq) { return seq.size(); }
  static uint length(const std::list<std::string>& aln) { return aln.front().size(); }

private:
  // options
  uint pk_level_;
  std::vector<float> alpha_;
  bool levelwise_;
  bool stacking_constraints_;
  int n_th_;
};

std::string
make_parenthesis(const std::vector<int>& bpseq, const std::vector<int>& plevel)
{
  const int n_support_parens=4;
  const char* left_paren="([{<";
  const char* right_paren=")]}>";

  std::string r(bpseq.size(), '.');
  for (int i=0; i!=(int)bpseq.size(); ++i)
  {
    if (bpseq[i]>=0 && i<bpseq[i])
    {
      int j=bpseq[i];
      if (plevel[i]<n_support_parens)
      {
        r[i]=left_paren[plevel[i]];
        r[j]=right_paren[plevel[i]];
      }
      else if (plevel[i]<n_support_parens+'Z'-'A'+1)
      {
        r[i]='A'+plevel[i]-n_support_parens;
        r[j]='a'+plevel[i]-n_support_parens;
      }
    }
  }
  return r;
}

void
make_interaction_pairs(const std::vector<int>& bpseq, const std::vector<int>& plevel, TInteractionPairs & interactions)
{
  const int n_support_parens = BracketType::MAX;
  interactions.resize(bpseq.size());

  const char* left_paren="([{<";
  const char* right_paren=")]}>";

  std::string r(bpseq.size(), '.');
  for (int i=0; i!=(int)bpseq.size(); ++i)
  {
    if (bpseq[i]>=0 && i<bpseq[i])
    {
    	int j=bpseq[i];
    	BracketType b_type = static_cast<BracketType>(plevel[i]);
    	interactions[i] = std::make_pair(b_type, j);
    	interactions[j] = std::make_pair(b_type, i);
    }
    else if (bpseq[i] == -1) {
    	interactions[i] = std::make_pair(MAX, -1);
    }
  }
}

template < class SEQ, class EN >
void
update_bpm(uint pk_level, const SEQ& seq, EN& en,
           const std::vector<int>& bpseq, const std::vector<int>& plevel,
           std::vector<float>& bp, std::vector<int>& offset)
{
  // update the base-pairing probability matrix by the previous result
  uint L=bpseq.size();
  //bp.resize(L);
  std::fill(bp.begin(), bp.end(), 0.0);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;

  std::vector<float> bpl;
  std::vector<int> offsetl;
  for (uint l=0; l!=pk_level; ++l)
  {
    // make the constraint string
    std::string str(L, '?');
    for (uint i=0; i!=bpseq.size(); ++i)
      if (bpseq[i]>=0 && (int)i<bpseq[i])
      {
        if ((int)l==plevel[i])
        {
          str[i]='('; str[bpseq[i]]=')';
        }
        else
        {
          str[i]=str[bpseq[i]]='.';
        }
      }

    // re-folding the seq with the constraint
    std::fill(bpl.begin(), bpl.end(), 0.0);
    en.calculate_posterior(seq, str, bpl, offsetl);
    assert(bp.size()==bpl.size());
    // update the base-pairing probability matrix
    for (uint k=0; k!=bp.size(); ++k) bp[k]+=bpl[k];
  }
#ifndef NDEBUG
  for (uint k=0; k!=bp.size(); ++k) assert(bp[k]<=1.0);
#endif
}

static
void
output_fa(std::ostream& os,
          const std::string& desc, const std::string& seq,
          const std::vector<int>& bpseq, const std::vector<int>& plevel)
{
  os << ">" << desc << std::endl
     << seq << std::endl
     << make_parenthesis(bpseq, plevel) << std::endl;
}

/* End of IPknot functions */

// reads the multiple alignment sequences into the IPknot Aln object
void getConsensusStructure(seqan::StockholmRecord<TBaseAlphabet> const & record, TInteractionPairs &consensusStructure, const char* constraint, IPknotFold const &){
	std::list<std::string> names;
	std::list<std::string> seqs;
	for (std::string name : record.sequence_names){
		names.push_back(name);
		seqs.push_back(record.seqences.at(name));
	}

	//for (size_t i = 0; i < record.seqences.size(); ++i){
		//names.push_back(record.seqNames[i]);
		//seqs.push_back(record.seqs[i]);
	//}

	uint pk_level=0;
	bool isolated_bp=false;
	int n_th=1;
	int n_refinement=0;
	std::vector< std::vector<float> > th;
	std::vector<float> alpha;
	bool levelwise=true;
	//bool max_pmcc=false; 	// TODO: Deleted the pseudo expected MMC solver? Don't use max_pmcc from ipknot.cpp

	if (th.empty())
	{
		th.resize(2);
		if (n_refinement==0)
		{
			th[0].resize(1, 1/(2.0+1)); // -g 2
			th[1].resize(1, 1/(4.0+1)); // -g 4
		}
		else
		{
			th[0].resize(1, 1/(1.0+1)); // -g 1
			th[1].resize(1, 1/(1.0+1)); // -g 1
		}
	}

	if (alpha.empty())
	{
		alpha.resize(th.size());
		for (uint i=0; i!=alpha.size(); ++i)
			alpha[i]=1.0/alpha.size();
	}
	pk_level=alpha.size();

    IPknot ipknot(pk_level, &alpha[0], levelwise, !isolated_bp, n_th);
    std::vector<float> bp;
    std::vector<int> offset;
    std::vector<int> plevel;

    IPknot::EnumParam<float> ep(th);
    std::vector<float> t(th.size());
    ep.get(t);

	Aln aln(names, seqs);

	BPEngineAln* mix_en=NULL;
	std::vector<BPEngineSeq*> en_s;
	std::vector<BPEngineAln*> en_a;
	const char* param=NULL;
//

    BPEngineSeq* e = new CONTRAfoldModel();
    en_s.push_back(e);
    en_a.push_back(new AveragedModel(e));

	//mix_en = new MixtureModel(en_a);

	BPEngineAln* en= mix_en ? mix_en : en_a[0];
	en->calculate_posterior(aln.seq(), bp, offset);

	std::vector<int> bpseq;

	ipknot.solve(aln.size(), bp, offset, t, bpseq, plevel);

	for (int i=0; i!=n_refinement; ++i)
	{
		update_bpm(pk_level, aln.seq(), *en, bpseq, plevel, bp, offset);
		ipknot.solve(aln.size(), bp, offset, t, bpseq, plevel);
	}

	make_interaction_pairs(bpseq, plevel, consensusStructure);

	std::cout << "IPknot: " << make_parenthesis(bpseq, plevel) << "\n";

	//output_fa(std::cout, aln.name().front(), aln.consensus(), consensusStructure, plevel);

    if (mix_en) delete mix_en;
    for (uint i=0; i!=en_s.size(); ++i) delete en_s[i];
    for (uint i=0; i!=en_a.size(); ++i) delete en_a[i];
}

#endif  // #ifndef APPS_RNAMOTIF_IPKNOT_UTILS_H_
