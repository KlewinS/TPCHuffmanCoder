#include "HuffmanCoder.h"

namespace TPC {

HuffmanCoder::HuffmanCoder(
	const char *fileName,
	int debugLevel,
	unsigned int huffmanRange,
	unsigned int maxCodeLength
	)
	: HuffmanCoder(debugLevel, huffmanRange, maxCodeLength)
{
	LoadHuffmanTableFromFile(fileName);
}

HuffmanCoder::HuffmanCoder(
	int debugLevel,
	unsigned int huffmanRange,
	unsigned int maxCodeLength
	)
	: mDebug(debugLevel)
	, mHuffmanTableExists(false)
	, mHuffmanRange(huffmanRange)
	, mLengthLimitedHuffman(false)
	, mLLMaxCodeLength(0)
	, mLLNumberOfWords(0)
	, mMaxCodeLength(maxCodeLength)
	, mRawDataMarker("")
	, mRawDataMarkerSize(0)
	, mLLRawDataMarker("")
	, mLLRawDataMarkerSize(0)
	, mHuffmanTable()
	, mTruncatedHuffmanTable()
	, mLengthLimitedHuffmanTable()
{
}

HuffmanCoder::~HuffmanCoder()
{
}

bool HuffmanCoder::LoadHuffmanTableFromFile(const char* fileName)
{
	if (mDebug > kInfo) std::cout << "Huffman table will be loaded from file " << fileName << std::endl;
	if (mHuffmanTableExists) {
		if (mDebug > kWarning) std::cout << "WARNING: Huffman table is already available" << std::endl;
		return true;
	}

	std::ifstream tableStream(fileName,std::ios::in);
	if (tableStream.is_open()) {
		std::string temp;
		getline(tableStream,temp,'\n');

		// Skip first lines of file
		while (temp != "Huffman codes:") getline(tableStream,temp,'\n');

		int value = 0;
		int length = 0;
		float weight = 0;
		std::string code = "";

		// read file and parse the lines
		do {
			// read line
			getline(tableStream,temp,'\n');

			// parse values
			value = atoi(temp.substr(temp.find("=")+1,temp.find("weight")).c_str());
			temp = temp.substr(temp.find("weight"));

			weight = atof(temp.substr(temp.find("=")+1,temp.find("length")).c_str());
			temp = temp.substr(temp.find("length"));

			length = atoi(temp.substr(temp.find("=")+1,temp.find("code")).c_str());
			temp = temp.substr(temp.find("code"));

			code = temp.substr(temp.find("=")+1);
			code = code.substr(code.size()-length);

			// Look for raw data marker in file
			if (value == 0 && weight != 0) {
				mRawDataMarker = code;
				mRawDataMarkerSize = code.size();
			}

			// Store Huffman code in map
			HuffmanCode hCode(value, length, weight, code);
			if (mHuffmanTable.find(value) == mHuffmanTable.end()) {
				mHuffmanTable[value] = hCode;
			} else {
				if (mDebug > kWarning) std::cout << "WARNING: Huffman code for value " << value << " exists already in Huffman table" << std::endl;
			}
			
			// if Huffman code is short enough, it is stored in the truncated Huffman table
			if ((unsigned int) length <= mMaxCodeLength) {
				if (mTruncatedHuffmanTable.find(value) == mTruncatedHuffmanTable.end()) {
					mTruncatedHuffmanTable[value] = hCode;
				} else {
					if (mDebug > kWarning) std::cout << "WARNING: Huffman code for value " << value << " exists already in truncated Huffman table" << std::endl;
				}
			}
		} while ((unsigned int) value < (mHuffmanRange-1));

		// done, closing file, creating raw data marker if not found and return
		if (mDebug > kInfo) std::cout << "Reading of Huffman table done." << std::endl;
		tableStream.close();

		if (mRawDataMarkerSize == 0) {
			// if no raw data marker was found, it is set to a serie of 1s, one more than the 
			// max length of Huffman code in order to not influence the correct decoding.
			for (unsigned int i = 0; i <= mMaxCodeLength; i++) mRawDataMarker += "1";
			mRawDataMarkerSize = mRawDataMarker.size();
			if (mDebug > kWarning) std::cout << "WARNING: No raw data marker was found, it is set to " << mRawDataMarker << std::endl;
		}

		mHuffmanTableExists = true;

		return true;
	} else {
		if (mDebug > kError) std::cout << "ERROR: Could not open file " << fileName << " to read the Huffman Table, no table is loaded!" << std::endl;
		return false;
	}
}

bool HuffmanCoder::EncodeValue(int value, std::string *code) {
	// if length limited Huffman is available and can be used, use it, else the standard one
	std::map<unsigned int, HuffmanCode>* map = NULL;
	std::string marker;
	if (mLengthLimitedHuffman) {
		map = &mLengthLimitedHuffmanTable;
		marker = mLLRawDataMarker;
	} else {
		map = &mTruncatedHuffmanTable;
		marker = mRawDataMarker;
	}

	// if the value is not in the corresponding map, the code is set to the raw data marker
	if (map->find(value) == map->end()) {
		*code = marker;
		return false;
	} else {
		*code = (*map)[value].code;
		return true;
	}
}

std::string HuffmanCoder::IntToBinaryString(int value, int length) {
	std::string returnString = ""; 
	if (value > pow(2,length)-1) {   
		if (mDebug > kWarning) std::cout << "WARNING: Value (" << value << ") bigger than a " << length << "bit number" << std::endl;
	} else {
		for (int i = length-1; i >= 0; i--) {   
			returnString += std::to_string((value >> i) & 0x1);
		}
	}
	return returnString;
}

unsigned int HuffmanCoder::BinaryStringToInt(std::string value, int length) {
	unsigned int returnValue = 0;
	std::string valueCopy = value;
	for (int i = 0; i < length; i++) {   
		returnValue += (std::stoi(valueCopy.substr(0,1)) << ((length-1)-i));
		valueCopy = valueCopy.erase(0,1);
		if (valueCopy.size() == 0) break;
	}   
	return returnValue;
}

bool HuffmanCoder::DecodeFirstValue(std::string *stream, unsigned int* value) {

	// if length limited Huffman is available and can be used, use it, else the standard one
	std::map<unsigned int, HuffmanCode>* map = NULL;
	unsigned int maxCodeLength;
	std::string marker;
	if (mLengthLimitedHuffman) {
		map = &mLengthLimitedHuffmanTable;
		maxCodeLength = mLLMaxCodeLength;
		marker = mLLRawDataMarker;
	} else {
		map = &mTruncatedHuffmanTable;
		maxCodeLength = mMaxCodeLength + 1; // +1 because the raw data marker is 1 bit longer than other codes
		marker = mRawDataMarker;
	}


	unsigned int l = (maxCodeLength<stream->size()) ? maxCodeLength : stream->size();
	// For each possible code lengths i, starting from 1 to the max code length,
	// the first i characters are compared to the Huffman codes.
	for (unsigned int i = 1; i <= l; i++) {
		// should speed up, since integer comparison should be faster than string comparison
		if (i == marker.size()) { 
			if(stream->substr(0,i) == marker) {
				stream->erase(0,i);
				*value = BinaryStringToInt(stream->substr(0,10));
				stream->erase(0,10);
				return true;
			}
		}

		for (std::map<unsigned int, HuffmanCode>::iterator h_code = (*map).begin(); h_code != (*map).end(); h_code++) {
			// should speed up, since integer comparison should be faster than string comparison
			if (i != h_code->second.length) continue; 

			if (stream->substr(0,i) == h_code->second.code) {
				*value = h_code->second.value;
				stream->erase(0,h_code->second.length);
				return true;
			}
		}
	}

	*value = -1;
	return false;
}

bool HuffmanCoder::GenerateLengthLimitedHuffman(unsigned int maxCodeLength, unsigned int numberOfWords) {
	if (mDebug > kInfo) std::cout << "INFO: Generating length-limited Huffman codes." << std::endl;
	
	if (!mHuffmanTableExists) {
		if (mDebug > kError) std::cout << "ERROR: Can't generate length-limited Huffman table because no Huffman table exists." << std::endl;
		return false;
	}

	mLLMaxCodeLength = maxCodeLength;
	mLLNumberOfWords = numberOfWords;

	/*
	 * To generate the length-limited Huffman codes, the so called package-merge algorithm is used.
	 */


	// 1) fill map with (one-element) vectors containing the value of the Huffman codes, 
	// sorted by their weight
	std::multimap<float, std::vector<int> > LengthLimitedHuffmanMap;
	for (std::map<unsigned int, HuffmanCode>::iterator it = mHuffmanTable.begin(); it != mHuffmanTable.end(); ++it) {
		std::vector<int> vec(1,it->second.value);
		LengthLimitedHuffmanMap.insert(std::pair<float,std::vector<int> >(it->second.weight, vec));
	}

	// Delete elements with the lowest weight from the map again until the the length of the map is 
	// at the given limit. The weights are added for the marker.
	float weightOfRawDataMarker = 0;
	while (LengthLimitedHuffmanMap.size() > mLLNumberOfWords) {
		weightOfRawDataMarker += LengthLimitedHuffmanMap.begin()->first;
		LengthLimitedHuffmanMap.erase(LengthLimitedHuffmanMap.begin());
	}
	std::vector<int> vec(1,-1);
	LengthLimitedHuffmanMap.insert(std::pair<float,std::vector<int> >(weightOfRawDataMarker, vec));
	int SizeOfSourceAlphabet = LengthLimitedHuffmanMap.size();

	// check wether the given alphabet fits into the given max. code range
	if (mLLMaxCodeLength < log(SizeOfSourceAlphabet)/log(2.)) {
		if (mDebug > kWarning) std::cout << "WARNING: The alphabet doesn't fit into the given maximum code length  of " << mLLMaxCodeLength << ". The max. code length has to be greater than / equal to " << log(SizeOfSourceAlphabet)/log(2.) << ". No length-limited Huffman table was generated." << std::endl;
		return false;
	}

	// 2) proceed with the actual package-merge algorithm
	std::multimap<float, std::vector<int> > mergedMap = PackageMerge(LengthLimitedHuffmanMap, mLLMaxCodeLength- 1);

	// calculate code lengths out of the occurence of the symbols according to the package-merge algorithm
	std::multiset<HuffmanCode,HuffmanCodeCompareCodelength> codeSet;
	for (int i = -1; i < (int) mHuffmanRange; i++) {
		int occurence = 0;
		std::multimap<float, std::vector<int> >::iterator it = mergedMap.begin();
		for (int set = 0; set < (2*SizeOfSourceAlphabet - 2); set++) {
			occurence += std::count(it->second.begin(), it->second.end(), i);
			++it;
		}
		if (occurence == 0) continue;

		if (i < 0) {
			// raw data marker
			mLLRawDataMarkerSize = occurence;
		} else {
			// rest
			HuffmanCode hCode(i, occurence, mHuffmanTable[(unsigned int) i].weight, "");
			codeSet.insert(hCode);
		}
	}

	// 3) generate actual Huffman codes out of the lengths of the codes
	int previousValue = -1;
	int previousCodeLength = -1;
	bool markerIncluded = false;
	for (std::multiset<HuffmanCode,HuffmanCodeCompareCodelength>::reverse_iterator rit = codeSet.rbegin(); rit != codeSet.rend(); ++rit) {
		int value;
		int length;
		if (rit->length <= mLLRawDataMarkerSize && !markerIncluded) {
			if (previousCodeLength == -1) previousCodeLength = mLLRawDataMarkerSize;
			rit--;
			value = -1;
			length = mLLRawDataMarkerSize;
			markerIncluded = true;
		} else {
			if (previousCodeLength == -1) previousCodeLength = rit->length;
			value = rit->value;
			length = rit->length;
		}
		std::string code = IntToBinaryString((previousValue+1) >> (previousCodeLength - length), length);
		if (value == -1) {
			mLLRawDataMarker = code;
		} else {
			HuffmanCode h_code = (*rit);
			h_code.code = code; 
			if (mLengthLimitedHuffmanTable.find(h_code.value) == mLengthLimitedHuffmanTable.end()) {
				mLengthLimitedHuffmanTable[h_code.value] = h_code;
			} else {
				if (mDebug > kWarning) std::cout << "WARNING: Huffman code for value " << h_code.value << " exists already in length-limited Huffman table" << std::endl;
			}
		}

		previousValue = (previousValue + 1) >> (previousCodeLength - length);
		previousCodeLength = length;
	}

	mLengthLimitedHuffman = true;

	return true;
}

std::multimap<float, std::vector<int> > HuffmanCoder::PackageMerge(std::multimap<float, std::vector<int> >& map, unsigned int iterations) {
	std::multimap<float, std::vector<int> > M = map;
	for (unsigned int i = 0; i < iterations; i++) {
		// Pack the map M
		std::multimap<float, std::vector<int> > P = Pack(M);
		// Merge it with the original one
		M = Merge(P,map);
	}

	return M;
}

std::multimap<float, std::vector<int> > HuffmanCoder::Merge(std::multimap<float, std::vector<int> >& mapOne, std::multimap<float, std::vector<int> >& mapTwo) {
	std::multimap<float, std::vector<int> > newMap(mapOne);
	for (std::multimap<float, std::vector<int> >::iterator it = mapTwo.begin(); it != mapTwo.end(); it++) {
		newMap.insert(*it);
	}
	return newMap;
}

std::multimap<float, std::vector<int> > HuffmanCoder::Pack(std::multimap<float, std::vector<int> >& map) {
	std::multimap<float, std::vector<int> > newMap;

	std::multimap<float, std::vector<int> >::iterator it = map.begin();
	for (unsigned int i = 0; i < map.size()/2; i++) {
		std::vector<int> newSymbolVector;
		float sumOfWeights = 0;

		for (int j = 0; j < 2; j++) {
			for (std::vector<int>::const_iterator itt = it->second.begin(); itt != it->second.end(); itt++) {
				newSymbolVector.push_back(*itt);
			}
			sumOfWeights += it->first;
			++it;
		}

		newMap.insert(std::pair<float,std::vector<int> >(sumOfWeights, newSymbolVector));
	}
	return newMap;
}
}
