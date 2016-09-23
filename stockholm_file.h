// ==========================================================================
//                              Stockholm_file.h
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

#ifndef APPS_RNAMOTIF_STOCKHOLM_FILE_H_
#define APPS_RNAMOTIF_STOCKHOLM_FILE_H_

// ============================================================================
// Forwards
// ============================================================================

namespace seqan{

struct Stockholm_;
typedef Tag<Stockholm_> Stockholm;

// ============================================================================
// Typedefs
// ============================================================================

// ----------------------------------------------------------------------------
// Type StockholmFileIn
// ----------------------------------------------------------------------------

/*!
 * @class StockholmFileIn
 * @signature typedef FormattedFile<Stockholm, Input> StockholmFileIn;
 * @extends FormattedFileIn
 * @brief Class for reading Stockholm files.
 */

typedef FormattedFile<Stockholm, Input>	StockholmFileIn;

// ----------------------------------------------------------------------------
// Typedef StockholmFileOut
// ----------------------------------------------------------------------------

/*!
 * @class StockholmFileOut
 * @signature typedef FormattedFile<Stockholm, Output> StockholmFileOut;
 * @extends FormattedFileOut
 * @brief Class for writing Stockholm files.
 *
 * @see StockholmHeader
 * @see StockholmRecord
 */

typedef FormattedFile<Stockholm, Output>  StockholmFileOut;


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class StockholmContext
 * @headerfile <stockholm_io.h>
 * @brief Context to use for Stockholm file I/O.
 *
 * @signature struct StockholmContext;
 */
struct StockholmContext
{
    /*!
     * @var String StockholmContext::buffer
     * @brief Buffer used during I/O.
     */
    seqan::CharString buffer;
};

// Alphabets usually found in Stockholm format files are Rna or AminoAcid
template <typename TAlphabet>
struct StockholmRecord {
	// SeqAn specific data	==========================
	// store the alignment given in the Stockholm file
	typedef seqan::String<TAlphabet> TSequence;
	typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;

	TAlign alignment;

	// Raw string data		===========================
	// header (GF tags -> tag value)
	std::unordered_map<std ::string, std::string > header;

	// seqence names -> sequence maps
	std::unordered_map<std::string, std::string > seqences;
	std::vector<std::string > sequence_names;

	// per column (GC) -> annotation string
	std::unordered_map<std::string, std::string > seqence_information;

	//TODO: Maybe support per-sequence (GS) and per-residue (GR) annotation
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

//template <typename T>
//struct MagicHeader<Stockholm, T> :
//    public MagicHeader<Nothing, T> {};

template <typename T>
struct MagicHeader<Stockholm, T>
{
    static unsigned char const VALUE[11];
};

template <typename T>
unsigned char const MagicHeader<Stockholm, T>::VALUE[11] =
{
    '#', ' ', 'S', 'T', 'O', 'C', 'K', 'H', 'O', 'L', 'M'  // Stockholm format's magic header
};

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Stockholm, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Stockholm, T>::VALUE[1] =
{
    ".msa"     // default output extension
};

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Stockholm, TDirection, TSpec>, TStorageSpec>
{
    typedef StockholmContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<Stockholm, TDirection, TSpec> >
{
    typedef Stockholm Type;
};


// ============================================================================
// Functions
// ============================================================================

// convient Stockholm variant
template <typename TAlphabet, typename TSpec>
inline void
readRecord(StockholmRecord<TAlphabet> & record, FormattedFile<Stockholm, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

}

#endif  // #ifndef APPS_RNAMOTIF_STOCKHOLM_FILE_H_
