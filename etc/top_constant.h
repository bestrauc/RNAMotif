/*
 * top_constant.h
 *
 *  Created on: Dec 15, 2015
 *      Author: vitrusky8
 */

#ifndef TOP_CONSTANT_H_
#define TOP_CONSTANT_H_

#define I 7

// ============================================================================
// Prerequisites
// ============================================================================

#define FASTA 0
#define EXTENDED_FASTA 1
#define DOTPLOT 2
#define UNKNOWN 3
#define NOT_FOUND 4
#define FASTQ 5
#define BPSEQ 6
#define EBPSEQ 7

#define VERBOSE_OUTPUT 1

// define different schemata (./lisa/incl/LISA/base/common.h)
#define LOGARITHMIC 0
#define SCALE       1
#define ORIGINAL    2
#define RIBOSUM     3

#define UNUSED 0
#define USED 1


//types used in the program
typedef unsigned TPosition;
typedef float TScoreValue;
typedef seqan::CharString TString;
typedef float TBioval;
typedef std::map<TPosition, TScoreValue> TMap;
typedef seqan::String<TMap > TMapLine;

typedef float TCargo;
typedef Graph<Directed<TCargo> > TDgraph;
typedef VertexDescriptor<TDgraph>::Type TDVertexDescriptor;
typedef EdgeDescriptor<TDgraph>::Type TDEdgeDescriptor;
typedef Graph<Undirected<TCargo> > TUgraph;
typedef VertexDescriptor<TUgraph>::Type TUVertexDescriptor;
typedef EdgeDescriptor<TUgraph>::Type TUEdgeDescriptor;

typedef Map<TPosition, TDVertexDescriptor> TMapDGraph;
typedef seqan::String<TMapDGraph > TMapDGraphStr;
typedef Map<TPosition, TUVertexDescriptor> TMapUGraph;
typedef seqan::String<TMapUGraph > TMapUGraphStr;


//typedef seqan::String<std::map<TPosition, TScoreValue> > TMapLine;
typedef seqan::RnaString TSequence;
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;      // align type

template <typename TString, typename TPosition>
struct fixedStructElement {
	TString method; // place the method and parameters used to compute the structure
//	seqan::String<unsigned> structure;
	seqan::String<TPosition> seqPos;
	seqan::String<TPosition> interPos;
};
template <typename TString, typename TBioval>
struct bioValStructElement {
	TString method; // place the method and parameters used to compute the structure
//	seqan::String<TBioval> val;
	seqan::String<TBioval> val;
};

struct vectGraphElement {
	std::vector<TUVertexDescriptor > uVertexVect;
	TUgraph interGraph; // this graph represents all the computed interaction edges
	std::vector<TDVertexDescriptor > dVertexVect;
	TDgraph interGraphUpdated;
};


template <typename TSequence, typename TString, typename TPosition, typename TBioval, typename TMapLine>
struct RnaStructSeq {
	TSequence seq;  // from fasta input
	TString qual;  // from fasta input
	TString id;  // from fasta input
	TString info; // raw info from bpseq or ebpseq input
	seqan::String<fixedStructElement<TString, TPosition> > structPairMate; //TODO use this field to collect all the structural information of this sequence
	seqan::String<bioValStructElement<TString, TBioval> > structBioVal;
	vectGraphElement bpProb;
//	seqan::String<std::map<TPosition, TValue> > *bpProb;
//	TMapLine *bpProb;
};

typedef RnaStructSeq<TSequence, TString, TPosition, TBioval, TMapLine> TRnaStruct;
typedef std::vector<TRnaStruct > TRnaVect;

//	Set up the generic scoring function
typedef Score<TScoreValue, ScoreMatrix<Value<TSequence>::Type, Ribosum65> > TScoringSchemeRib;
typedef Score<TScoreValue, RnaStructureScore<TScoringSchemeRib> > TScoringSchemeRibStruct;

struct alignPos
{
	unsigned seqPos1;
	unsigned seqPos2;
	unsigned interMate1;
	unsigned interMate2;
	float interScore;
};

template <typename TScoreValue, typename TRnaStruct, typename TPosition,
typename TMapLine, typename TMapDGraphStr, typename TMapUGraphStr,
typename TAlign, typename TScoringSchemeRibStruct>
struct RnaStructAlign
{
//public:
	TRnaStruct rna1;
	TRnaStruct rna2;
	TMapLine mapline; // string with size seq1 storing all the aligned lines
	TMapDGraphStr lambdaGraph;
	TMapDGraphStr upboundGraph;
	TAlign align;
	TScoreValue alignscore;
	TScoringSchemeRibStruct riboStructScore; // it is necessary to embed the score system in the alignment structure
	TAlign bestalign;
	TScoreValue bestalignscore;
	std::vector<alignPos> xlines;
	unsigned xlinesiterator;
};// rnaStructAlign;

typedef RnaStructAlign<TScoreValue, TRnaStruct, TPosition, TMapLine, TMapDGraphStr, TMapUGraphStr, TAlign, TScoringSchemeRibStruct> TRnaAlign;

typedef std::vector<TRnaAlign> TRnaAlignVect;


template <typename TScoreValue, typename TSequence, typename TPosition, typename TMapLine, typename TAlign>
struct RnaStructAlignOld
{
//public:
	TSequence seq1;
	TSequence seq2;
	TMapLine mapline; // string with size seq1 storing all the aligned lines
	TAlign align;
	TScoreValue alignscore;
	TAlign bestalign;
	TScoreValue bestalignscore;
};// rnaStructAlign;

typedef RnaStructAlignOld<TScoreValue, TSequence, TPosition, TMapLine, TAlign> TRnaAlignOld;

#endif /* TOP_CONSTANT_H_ */
