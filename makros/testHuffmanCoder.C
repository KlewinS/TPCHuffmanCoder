#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TList.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPad.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPaletteAxis.h"
#include "../../Huffman/Huffman.h"
#include "../HuffmanCoder.h"
#include "AliHLTHuffman.h"
#include "../../generator/DataGenerator.h"
#include "TStyle.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>

void printHist(TH1* hist, TString baseName);
void printHist(TH2* hist, TString baseName);
void printProf(TProfile* prof, TString baseName,TString axis);

void testHuffmanCoder(TString CurrentMacroName, float rateForTable)
{
	////////////////////////////////////////////////////////////////////////////////
    // input / output files
//	const char* dataFiles = "datafiles_firstHalf.txt";
//	const char* dataFiles = "datafiles_onlyfirst.txt";
//	const char* dataFiles = "datafiles_firstTen.txt";
//	const char* dataFiles = "datafiles_secondHalf.txt";
	const char* dataFiles = "datafiles_all.txt";
	TString HuffmanBaseFileName("../../Huffman/TPCRawSignalDifference_HuffmanTable_");
	TString MappingFileName("../../generator/mapping.dat");
	TString PedestalFileName("../../generator/pedestal-statistics.txt");

	TString baseName(CurrentMacroName);
//		baseName.Remove(baseName.Last('.'));
		baseName += "_";
		if (rateForTable < 10) baseName += "0";
		baseName += rateForTable;
		baseName += "tableRate_";
		baseName.Remove(baseName.Last('.'),2);
	TString plotBaseName(baseName);
    	plotBaseName.Insert(plotBaseName.Last('/'),"/Plots");
    TString folderName(plotBaseName);
    	folderName.Remove(folderName.Last('/'));
    	folderName = "mkdir -p " + folderName;
    	system(folderName.Data());
	TString targetFileName(baseName);
		targetFileName += "summary.root";

	////////////////////////////////////////////////////////////////////////////////
	// configuration

	const int nFrames = -1;//10; 
	int TimeFrameNumber = 0;
	// DDL range
    const int DDLmin = 0;
    const int DDLmax = 71;

	// Padrow range
	const int PadrowMin = -1;
	const int PadrowMax = -1;

	// signal bit length
	const int signalBitLength=10;
	const int signalRange=0x1<<signalBitLength;

	// the length of one channel
	const int maxChannelLength=1024;

	// zero suppression
	const int ZSbaseline = 5;
	const int ZSthreshold = 2;

	// occupancy
	const float occupancyMax = 1.;//0.4;

	// Data format
	const int oldHeaderSize = 50; // bit
	const int newHeaderSize = 50; // bit
	const bool enableHeader = false;

	// data generator
	const int numberOfGenerators = 7;
	const int mode = 3;
	const int nCol = 5;
	const float rate = 5.0;

	////////////////////////////////////////////////////////////////////////////////
	// Huffman
	HuffmanBaseFileName += (mode==0||mode==2)?nCol:rateForTable;//nCol;
	HuffmanBaseFileName += "mergedCollisions";
    TString HuffmanTableNameRoot(HuffmanBaseFileName);
    HuffmanTableNameRoot += ".root";
    TString HuffmanTableNameTxt(HuffmanBaseFileName);
    HuffmanTableNameTxt += ".txt";
    AliHLTHuffman* hltHuffman=NULL;

	TFile* htf = TFile::Open(HuffmanTableNameRoot);
	if (!htf || htf->IsZombie()) {
		std::cerr << "ERROR: Can't open file " << HuffmanTableNameRoot << std::endl;
		return;
	}   
	TObject* obj = NULL;
	htf->GetObject("TPCRawSignalDifference",obj);
	if (obj==NULL) {
		std::cerr << "ERROR: Can't load Huffman decoder object " << "TPCRawSignalDifference" << "from file " <<     HuffmanTableNameRoot << std::endl;
		return;
	}
	hltHuffman = (AliHLTHuffman*)obj;

    const unsigned int numDiffBits = 2;
    unsigned int Bits[numDiffBits] = {10,12}; 
    const unsigned int numDiffWords = 4; 
    unsigned int Words[numDiffWords] = {50,70,90,110};
    unsigned int numberOfHuffmans = numDiffBits*numDiffWords;
    HuffmanCoder* huffman[numberOfHuffmans];
    TPC::HuffmanCoder* newHuffman[numberOfHuffmans];
    unsigned int count = 0;
    unsigned int countNew = 0;
	for (unsigned int i = 0; i < numDiffBits; i++) {
		for (unsigned int j = 0; j < numDiffWords; j++) {
			std::cout << "Generating Huffman encoder for " << Bits[i] << " bits and " << Words[j] << " words" << std::endl << "\t";
				clock_t t1 = clock();
			huffman[count] = new HuffmanCoder(HuffmanTableNameTxt.Data());
				clock_t t2 = clock();
			newHuffman[count] = new TPC::HuffmanCoder(HuffmanTableNameTxt.Data());
				clock_t t3 = clock();
				double durationOld = double(t2 - t1) / CLOCKS_PER_SEC;
				double durationNew = double(t3 - t2) / CLOCKS_PER_SEC;
				std::cout << "time old: " << durationOld << "\ttime new: " << durationNew << std::endl;

			TString codeTableName("VerilogTruncatedHuffmanCodes.v");
			TString decoderCodeTableName("VerilogTruncatedHuffmanDecoderCodes.v");
			TString lengthTableName("VerilogTruncatedHuffmanLengths.v");
			newHuffman[count]->WriteVerilogEncoderTable(codeTableName.Data(), lengthTableName.Data());
			newHuffman[count]->WriteVerilogDecoderTable(decoderCodeTableName.Data());

			if (huffman[count] && newHuffman[count]) {
				bool result1 = huffman[count]->GenerateLLHuffmanCode(Bits[i],Words[j]);
				newHuffman[count]->SetLLRawDataMarkerSize(count+1);
				bool result2 = newHuffman[count]->GenerateLengthLimitedHuffman(Bits[i],Words[j]);
				if (result1 && result2) {
					TString codeTableNameLL("VerilogLengthLimitedHuffmanCodes_");
					codeTableNameLL += count;
					codeTableNameLL += ".v";
					TString decoderCodeTableNameLL("VerilogLengthLimitedHuffmanDecoderCodes_");
					decoderCodeTableNameLL += count;
					decoderCodeTableNameLL += ".v";
					TString lengthTableNameLL("VerilogLengthLimitedHuffmanLengths_");
					lengthTableNameLL += count;
					lengthTableNameLL += ".v";
					newHuffman[count]->WriteVerilogEncoderTable(codeTableNameLL.Data(), lengthTableNameLL.Data());
					newHuffman[count]->WriteVerilogDecoderTable(decoderCodeTableNameLL.Data());

					++count;
				}
				else { 
					delete huffman[count];
					delete newHuffman[count];
					std::cout << "WARNING: deleted last Huffman encoder" << std::endl;
				}
			}       
		}       
	}       
	numberOfHuffmans = count;

	////////////////////////////////////////////////////////////////////////////////

	for (unsigned int i = 0; i < numberOfHuffmans; i++) {
		std::cout << i << std::endl;
		std::string temp = "";
		std::string temp2 = "";
		unsigned int value;
		for (unsigned int j = 0; j < 2047; j++) {
			huffman[i]->findHuffmanCode(j,&temp);
			if (!newHuffman[i]->EncodeValue(j,&temp2)) {
				temp2 += newHuffman[i]->IntToBinaryString(682);
			}
			//if (temp != temp2) std::cout << "\t" << j << "\t" << temp << "\t" << temp2 << std::endl;
			std::cout << "\t" << j << "\t" << temp << "\t" << temp2;
			newHuffman[i]->DecodeFirstValue(&temp2,&value);
			std::cout << "\t" << value << "\t" << temp2.size() << std::endl;
		}
	}


	////////////////////////////////////////////////////////////////////////////////
	// Histograms

	////////////////////////////////////////////////////////////////////////////////
	// Data generator
/*
	TPC::DataGenerator* dg;
		dg = new TPC::DataGenerator(mode);
		dg->SetMappingFileName(MappingFileName);
		dg->SetPedestalFileName(PedestalFileName);
		dg->SetDDLRange(DDLmin,DDLmax);
		dg->SetPadrowRange(PadrowMin,PadrowMax);
		dg->InitZeroSuppression(ZSthreshold, ZSbaseline);
		dg->SetApplyZeroSuppression(false);
		dg->SetNormalizeChannels(false);
		dg->SetApplyCommonModeEffect(false);
		dg->Init( (mode==0||mode==2)?nCol:rate, dataFiles);


	////////////////////////////////////////////////////////////////////////////////
	// start

	int result = 0;
	while ((TimeFrameNumber++ < nFrames || nFrames < 0) && result >= 0) {
		result = dg->SimulateFrame();
    
		if (result < 0) {
	  		std::cout << "Simulation of timeframe failed with error code " << result << "; aborting at timeframe no " << TimeFrameNumber << std::endl;
	  		continue;
		} 

		std::vector<unsigned int> chindices=dg->GetChannelIndices();
		for (std::vector<unsigned int>::iterator index = chindices.begin(); index != chindices.end(); ++index) {
			TPC::DataGenerator::ChannelDesc_t desc = dg->GetChannelDescriptor(*index);
			if (desc.ptr == NULL) continue;
			int DDLnumber = ((*index)&0xffff0000)>>16;
			int padrow = desc.padrow;
			if (DDLnumber > 71) padrow += 63;
	//			std::cout << desc.occupancy << std::endl;
			if (desc.occupancy > occupancyMax) continue;
			int signal = 0;
			int lastSignal = -1;
			int diff = 0;
			std::string temp = "";
			std::string temp2 = "";
			double cf = 0;
			AliHLTUInt64_t codeLength;
			AliHLTUInt64_t v;
			int oldLength = 0;
			int newLength[numberOfHuffmans+1];
			for (int i = 0; i < numberOfHuffmans+1; i++) {
				newLength[i] = 0;
			}
			for (int i = desc.size-1; i >= 0; i--) {
				signal = desc.ptr[i];
				if (signal < 0 || signal >= signalRange) continue;

//				if (i == desc.size-1 || desc.ptr[i+1] >= signalRange) // fist signal
				if (lastSignal < 0) // fist signal
					diff = signal;
				else
					diff = (signal - lastSignal);

				lastSignal = signal;
				diff += signalRange;
				
				v = diff;
				hltHuffman->Encode(v,codeLength);
				if (codeLength > 12) codeLength = 22;
				newLength[numberOfHuffmans] += codeLength;

				for (int j = 0; j < numberOfHuffmans; j++) {
					if (!huffman[j]->findLLHuffmanCode(diff,&temp)) {
						newLength[j] += 10;
					}
					newLength[j] += temp.size();

//					newHuffman[j]->EncodeValue(diff,&temp2);
//					if (temp != temp) std::cout << diff << "\t" << j << "\t" << temp << "\t" << temp2 << std::endl;
				}
				oldLength += 10;
			}

			while (oldLength % 40 != 0) oldLength++;
			if (enableHeader) oldLength += oldHeaderSize;
			for (int j = 0; j < numberOfHuffmans; j++) {
				while (newLength[j] % 10 != 0) newLength[j]++;
				if (enableHeader) newLength[j] += newHeaderSize;
				cf = (double)oldLength / (double)newLength[j];
//				hCompFactorLLHuffman[j]->Fill(desc.occupancy,cf);
			}

			while (newLength[numberOfHuffmans] % 10 != 0) newLength[numberOfHuffmans]++;
			if (enableHeader) newLength[numberOfHuffmans] += oldHeaderSize;
			cf = (double)oldLength / (double)newLength[numberOfHuffmans];
//			hCompFactorStdHuffman->Fill(desc.occupancy,cf);
		}
	}

	if (dg) delete dg;
*/
	////////////////////////////////////////////////////////////////////////////////
	// print and delete Histograms

	TFile* of = TFile::Open(targetFileName, "RECREATE");
	if (!of || of->IsZombie())
	{
		std::cerr << "can not open file " << targetFileName << std::endl;
		return;
	}

	////////////////////////////////////////////////////////////////////////////////
	// clean-up

	if (hltHuffman) {
		delete hltHuffman;
	}

	for (int i = 0; i < numberOfHuffmans; i++) { 
		if (huffman[i]) {
			delete huffman[i];
		}
		if (newHuffman[i]) {
			delete newHuffman[i];
		}
	}


}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// Additional functions
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

void printHist(TH1* hist, TString baseName) {
    TCanvas* cnv2 = new TCanvas("cnv2", "cnv2",1000,1000);
    cnv2->SetLeftMargin(0.11);
    cnv2->SetRightMargin(0.10);
    
    cnv2->cd();
    hist->Draw();
	cnv2->SetLogy();
	cnv2->SetLogx();
    cnv2->SetGrid(1,1);
    cnv2->Update();
    TString printName(baseName);
    printName += hist->GetName();
    printName += ".pdf";
    cnv2->Print(printName);
    delete cnv2;
}

void printProf(TProfile* prof, TString baseName, TString axis) {
    TCanvas* cnv3 = new TCanvas("cnv3", "cnv3",1000,1000);
    cnv3->SetLeftMargin(0.11);
    cnv3->SetRightMargin(0.10);
    
    cnv3->cd();
	if (axis == "X") prof->GetYaxis()->SetRangeUser(0,10);
	if (axis == "Y") prof->GetYaxis()->SetRangeUser(0,1);
    prof->Draw();
    cnv3->SetGrid(1,1);
    cnv3->Update();
    TString printName(baseName);
    printName += prof->GetName();
    printName += ".pdf";
    cnv3->Print(printName);
    delete cnv3;
}

void printHist(TH2* hist, TString baseName) {
    TCanvas* cnv1 = new TCanvas("cnv1", "cnv1",1000,1000);
    cnv1->SetLeftMargin(0.11);
    cnv1->SetRightMargin(0.10);
    
    cnv1->cd();
    hist->Draw("COLZ1");
//    hist->Draw("lego2z");
//    cnv1->SetLogx();
    cnv1->SetLogz();
//    cnv1->SetPhi(70);
    cnv1->Update();
//    TPaletteAxis *palette = NULL;
//    palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
//    palette->SetY1NDC(0.4);
//    cnv1->Update();
    TString printName(baseName);
    printName += hist->GetName();
    printName += ".pdf";
    cnv1->Print(printName);
    delete cnv1;
}

