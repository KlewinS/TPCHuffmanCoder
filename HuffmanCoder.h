#ifndef _HUFFMANCODER_H
#define _HUFFMANCODER_H

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <algorithm>

namespace TPC {

/*
 * @class HuffmanCoder
 * Compression class for truncated Huffman encoding and decoding.
 *
 */

class HuffmanCoder 
{
	public:
		/*
		 * constructor for already existing Huffman table
		 */
		HuffmanCoder(
			const char *fileName,
			int debugLevel = 3,
			unsigned int huffmanRange = 2048,
			unsigned int maxCodeLength = 12
			);
	
		/*
		 * constructor without loading Huffman table
		 * TODO: adding of training values has to be implemented
		 */
		HuffmanCoder(
			int debugLevel = 3,
			unsigned int huffmanRange = 2048,
			unsigned int maxCodeLength = 12
			);

		~HuffmanCoder();

		/*
		 * Encodes the given value, encoded bit string is written in 'code'.
		 * Return value indicades wether value was found in table or not. If not,
		 * a marker is written to 'code'.
		 */
		bool EncodeValue(int value, std::string *code);

		/*
		 * Decodes the first value in the given data stream and writes it to 'value'.
		 * If no code could be found in the whole stream, the return value is false.
		 */
		bool DecodeFirstValue(std::string *stream, unsigned int* value);

		/*
		 * Generates the Huffman table for the length-limited Huffman according to the 
		 * givene parameters from the existing Huffman table.
		 */
		bool GenerateLengthLimitedHuffman(unsigned int maxCodeLength, unsigned int numberOfWords);

		/*
		 * Returns max code length of length-limited Huffman
		 */
		unsigned int GetLLMaxCodeLength() { return mLLMaxCodeLength;};

		/*
		 * Returns number of words in length-limited Huffman table
		 */
		unsigned int GetLLNumberOfWords() { return mLLNumberOfWords;};

		/*
		 * struct for Huffman code
		 */
		struct HuffmanCode
		{
			unsigned int	value;			// encoded value
			float 			weight;			// weight of value
			unsigned int	length;			// length of Huffman code
			std::string 	code;			// Huffman code
			HuffmanCode(unsigned int Value = 0, int Length = 0, float Weight=0, std::string Code = "0")
				: value(Value)
				, length(Length)
				, weight(Weight)
				, code(Code)
			{
//				if (length != code.size()) if (mDebug > kWarning) std::cout << "Warning: Given length does not match the code size" << std::endl;
			};
			HuffmanCode(const HuffmanCode& toCopy)
				: value(toCopy.value)
				, length(toCopy.length)
				, weight(toCopy.weight)
				, code(toCopy.code)
			{};
		};
		struct HuffmanCodeCompareCodelength
		{
			bool operator() (const HuffmanCode& lhs, const HuffmanCode& rhs) const { return lhs.length < rhs.length;};
		};

		/*
		 * Converts an integer value into a string with the binary representative with the given length
		 */
		std::string IntToBinaryString(int value, int length = 10);

		/*
		 * Converts a binary string into an integer value
		 */
		unsigned int BinaryStringToInt(std::string value, int length = 10);


	private:
		/*
		 * Loads a Huffman table from the specified file. The format has to be the following one:
		 *
		 * ########################
		 * ...
		 * Huffman codes:
		 * value=XX		weight=XX	length=XX	code=XX
		 * value=XX		weight=XX	length=XX	code=XX
		 * ...
		 *
		 * #######################
		 *
		 * everything above the line 'Huffman codes:' will be ignored. The value of 'value' is
		 * casted as integer, of 'weight' as float, of 'length' as integer and of 'code' as string
		 * with length of 'length'. The separators of the columns doesn't matter because the '=' 
		 * sign is looked for.
		 */
		 /*
		  * TODO: correct implementation of generation of raw data marker
		  */
		bool LoadHuffmanTableFromFile(const char * fileName);

		// Function which processes the package-merge algorithm
		std::multimap<float, std::vector<int> > PackageMerge(std::multimap<float, std::vector<int> >& map, unsigned int iterations);

		// Merges two multimaps
		std::multimap<float, std::vector<int> > Merge(std::multimap<float, std::vector<int> >& mapOne, std::multimap<float, std::vector<int> >& mapTwo);

		// Packs one map
		std::multimap<float, std::vector<int> > Pack(std::multimap<float, std::vector<int> >& map);

		// debug level
		int mDebug;
		enum DebugLevel {
			kInfo = 2,
			kWarning = 1,
			kError = 0
		};

		// Indicates wether a Huffman table is already generated/loaded
		bool mHuffmanTableExists;

		// Range of values to encode with Huffman
		unsigned int mHuffmanRange;

		// Indicades wether Huffman codes are length limited or not
		bool mLengthLimitedHuffman;

		// Max length of length-limited Huffman code
		unsigned int mLLMaxCodeLength;

		// Number of words contained in length-limited Huffman code table
		unsigned int mLLNumberOfWords;

		// Max length of truncated Huffman
		unsigned int mMaxCodeLength;

		// marker for raw data (not found in truncated Huffman table)
		std::string mRawDataMarker;
		unsigned int mRawDataMarkerSize;

		// marker for raw data (not found in length-limited Huffman table)
		std::string mLLRawDataMarker;
		unsigned int mLLRawDataMarkerSize;

		// Table with Huffman Codes
		std::map<unsigned int, HuffmanCode> mHuffmanTable;
		std::map<unsigned int, HuffmanCode> mTruncatedHuffmanTable;
		std::map<unsigned int, HuffmanCode> mLengthLimitedHuffmanTable;
};

}
#endif
