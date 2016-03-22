template <typename TVect, typename TFileList>
void readBpseqUpdateRnas(TVect & rnaSeqs, TFileList const & fileList)
{
	CharString fileBp;
	CharString tmpStr;
	bool flagSeqUsed=UNUSED;
	fixedStructElement<TString, TPosition> tmpStructPairMate;
	RnaStructSeq<TSequence, TString, TPosition, TBioval, TMapLine> tmpStructSeq;
	for(unsigned i=0;i<fileList.size();++i)
	{
		std::cout << fileList[i] << std::endl;
// Get path to example file.
		fileBp = getAbsolutePath(toCString(fileList[i]));
// Open Bpseq input file.
	    BpseqFileIn bpseqIn(toCString(fileBp));
// Get file extension
	    std::vector<std::string> formats = getFileExtensions(bpseqIn);
	    std::cout << "file format in bpseq = " << formats[0] << std::endl;
// Attach Bpseq to standard output.
	    BpseqFileOut bpseqOut(bpseqIn);
	    open(bpseqOut, std::cout, Bpseq());
//Open an output file
	    BpseqFileOut bpseqOutB(bpseqIn);
	    tmpStr = "/demo/"+std::to_string(i)+"-seq.bpseq";
	    open(bpseqOutB, toCString(getAbsolutePath(toCString(tmpStr)))); //TODO for now files are called with numbers but they can be called with the original file name
// Check the output file format
	    formats = getFileExtensions(bpseqOut);
	    std::cout << "file format out bpseq = " << formats[0] << std::endl;
	    std::cout << std::endl;
// Write all the output on the std::out and on the files
	    BpseqHeader headerBp;
	    BpseqRecord recordBp;
	    while (!atEnd(bpseqIn))
	    {
// Copy over headerBp.
			readHeader(headerBp, bpseqIn);
//			writeHeader(bpseqOut, headerBp);  // Print on video the file contents
			writeHeader(bpseqOutB, headerBp);  // Print on an output file the file contents
// Copy the Bpseq file record by record.
	        readRecord(recordBp, bpseqIn);
//	        writeRecord(bpseqOut, recordBp); // Print on video the file contents
	        writeRecord(bpseqOutB, recordBp); // Print on an output file the file contents
	        flagSeqUsed=UNUSED;
		    for(unsigned i=0; i<length(rnaSeqs);++i)
		    {
		    	if(!seqCompare(rnaSeqs[i].seq, recordBp.seq))
		    	{
//		    		std::cout << "seq recordBp.seq" << " = " << recordBp.seq << std::endl;
//		    		std::cout << "seq rnaSeqs" << i << " = " << rnaSeqs[i].seq << std::endl;
//		    		std::cout << "header " << headerBp[0].key << tmpStr << std::endl;
//		    		tmpStructPairMate.method = fileList[i];
		    		tmpStructPairMate.method = headerBp[0].value; // TODO in this position must be added the tool name and the parameters and not the sequence name
		    		tmpStructPairMate.seqPos = recordBp.seqPos;
		    		tmpStructPairMate.interPos = recordBp.interPos;
		    		appendValue(rnaSeqs[i].structPairMate, tmpStructPairMate);
		    		flagSeqUsed=USED;
		    		elementClear(tmpStructPairMate);
		    	}
		    }
		    if(flagSeqUsed==UNUSED)
		    {
// create the new item in the sequence structure
	    		tmpStructPairMate.method = headerBp[0].value;  // TODO in this position must be added the tool name and the parameters and not the sequence name
//		    	tmpStructPairMate.method = fileList[i];
		    	tmpStructPairMate.seqPos = recordBp.seqPos;
	    		tmpStructPairMate.interPos = recordBp.interPos;
// Add the new reported sequence in the rnaSeqs vector
	    		tmpStructSeq.seq = recordBp.seq;
	    		tmpStructSeq.id = headerBp[0].value;
				appendValue(tmpStructSeq.structPairMate, tmpStructPairMate);
				appendValue(rnaSeqs,tmpStructSeq);
//		    	std::cout << "seq recordBp.seq has not been used" << " = " << recordBp.seq << std::endl;
		    	flagSeqUsed=USED;
		    	elementClear(tmpStructPairMate);
		    }
	    }
// Print the last cached sequences from .bpseq file
//		for (unsigned i = 0; i < length(recordBp.seqPos); ++i)
//		{
//			std::cout << recordBp.seqPos[i] << "\t";
//		}
//		std::cout << std::endl;
//		for (unsigned i = 0; i < length(recordBp.seq); ++i)
//		{
//			std::cout << recordBp.seq[i] << "\t";
//		}
//		std::cout << std::endl;
//		for (unsigned i = 0; i < length(recordBp.interPos); ++i)
//		{
//			std::cout << recordBp.interPos[i] << "\t";
//		}
//		std::cout << std::endl;
//	    std::cout << recordBp.seq << std::endl;

// Print all the content of rnaSeqs vector

// Fill the RNA fields
//	    appendValue(header, record);
// Check the rnaSeqs[i].seq vector and in case you find the right sequence append the structure info
	}
}
//template <typename TVect>
//void createRnaBppGraph(TVect & rnaSeqs_i, unsigned const & seqSize)
//{
//	TUVertexDescriptor uVertex[seqSize];
//	TUEdgeDescriptor uEdge[seqSize];
//	rnaSeqs_i.seq;
//}

template <typename TOption, typename TVect>
void bppInteractionGraphBuild(TOption const & options, TVect  & rnaSeqs)
{
//variable used to plot or not the energy matrix (for now is disabled but the function is tested)
	bool write_dot_plot = 0; //TODO decide if provide this option in the cmd options
// Create a c-style string object for str:
	String<char, CStyle> cStr;
// Now use cStyle as char array:
//	Call RNAfold in order to fix the energy score matrix
//	TDgraph interGraph;
//	TDVertexDescriptor dVertex;
//	TDEdgeDescriptor dEdge;
//	std::vector<TDVertexDescriptor > dVertexVect;
//	std::vector<TDEdgeDescriptor > dEdgeVect;
	for(unsigned i=0;i<length(rnaSeqs);++i)  // TODO Execute this part in PARALLEL
	{
		cStr = rnaSeqs[i].seq;
//		CharString curSeq = rnaSeqs[0].seq;
		std::cout << toCString(cStr) << std::endl;
		for(unsigned j=0;j<length(rnaSeqs[i].seq);++j)
		{
//			dVertex = addVertex(interGraph);
			rnaSeqs[i].bpProb.uVertexVect.push_back(addVertex(rnaSeqs[i].bpProb.interGraph));
			rnaSeqs[i].bpProb.dVertexVect.push_back(addVertex(rnaSeqs[i].bpProb.interGraphUpdated));
		}
//TODO once is given support to introduce several BPP from the extended bpseq file or from the dot plot file in this position must be placed a condition capable to discriminate which bpp matrix to use
		if(1) // if(dotplot or extended bpseq data are not present)
		{
//			computeBppMatrix(options, toCString(cStr), write_dot_plot, rnaSeqs[i].bpProb.interGraph, rnaSeqs[i].bpProb.uVertexVect);
			computeBppMatrix(options, toCString(cStr), write_dot_plot, rnaSeqs[i]);
	//		createRnaBppGraph(rnaSeqs[i], length(rnaSeqs[i].seq));
		}else{
// TODO read the BPP matrix from the files
		}
// TODO add the filtering step that involve the biological input and the majority voter of predicted structures

		std::cout << "Input Degree nodo 0 = " << inDegree(rnaSeqs[i].bpProb.interGraph, rnaSeqs[i].bpProb.uVertexVect[0]) << std::endl;
		std::cout << "Output Degree nodo 0 = " << outDegree(rnaSeqs[i].bpProb.interGraph, rnaSeqs[i].bpProb.uVertexVect[0]) << std::endl;
//		std::cout << rnaSeqs[i].bpProb.interGraph << std::endl;

	    // Output distances of shortest paths
//	    Iterator<TUgraph, VertexIterator>::Type it(rnaSeqs[i].bpProb.interGraph);
		typedef Iterator<TUgraph, OutEdgeIterator>::Type TOutEdgeIterator;
//TODO place the printing part in an external function printWeightedGraph()
		for(unsigned j=0;j<length(rnaSeqs[i].bpProb.uVertexVect);++j)
		{
			TOutEdgeIterator it(rnaSeqs[i].bpProb.interGraph, rnaSeqs[i].bpProb.uVertexVect[j]);
			for(;!atEnd(it);goNext(it)) {
	//		    std::cout << "Edge id = " << it;
// build the directed graph that will be updated over the iterations
//				addEdge(rnaSeqs[i].bpProb.interGraphUpdated, source(*it), target(*it), cargo(*it));
//				addEdge(rnaSeqs[i].bpProb.interGraphUpdated, target(*it), source(*it), cargo(*it));
			    std::cout << "source = " << source(*it) << "\ttarget = " << target(*it) << "\tcargo = " << cargo(*it) << std::endl;
			}
		}

//	    while (!atEnd(it))
//	    {
//	        std::cout << "Distance from 0 to " << getValue(it) << ": ";
//	        std::cout << getProperty(rnaSeqs[i].bpProb.uVertexVect, getValue(it)) << std::endl;
//	        goNext(it);
//	    }

//		dEdge = addEdge(g,v0,v1);

	}

}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
	if (argc==1)
		std::cout << "type ./lara_gu --help to get the parameters table (-i option is mandatory)" << std::endl;

    seqan::ArgumentParser parser;

    Options options;

    setupArgumentParserCmd(parser, options.cmd);
    ArgumentParser::ParseResult res;
	res = parseCmd(options.cmd, parser, argc, argv); //riempio le due strutture che contengono le opzioni

	if (res != ArgumentParser::PARSE_OK)  // controllo argomenti
	   return res == ArgumentParser::PARSE_ERROR;

	loadParameterFile(options.cmd, options.parameters);

	TRnaVect rnaSeqs;
	ReadFastaFastqFormat(options, rnaSeqs);

	for (unsigned i = 0; i < length(rnaSeqs); ++i)
	{
	   std::cout << rnaSeqs[i].id << '\t' << rnaSeqs[i].seq << '\t' << rnaSeqs[i].qual << '\n';
	}

//	Read file list produced by Ipknot, dotknot and other tools
	std::vector<std::string> fileList;
	splitBpseqFilenames(options, fileList);

//	Read .bpseq file list and update the data structure
	readBpseqUpdateRnas(rnaSeqs, fileList);

//	Add the weight interaction edges vector map in the data structure
	bppInteractionGraphBuild(options, rnaSeqs);

//	printAllRnaSeqs(rnaSeqs);


	return 0;
}

