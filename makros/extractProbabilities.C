#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TList.h"
#include "TPad.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector.h"
#include "TCanvas.h"
#include "TPaletteAxis.h"
#include "../HuffmanCoder.h"
#include "AliHLTHuffman.h"
#include "../../generator/DataGenerator.h"

#include <iostream>
#include <iomanip>
#include <string>

void printHist(TH1* hist, TString baseName);
void printHist(TH2* hist, TString baseName);

void extractProbabilities(TString CurrentMacroName, int numberOfBits)
{
	////////////////////////////////////////////////////////////////////////////////
    // input / output files
	const char* dataFiles = "datafiles_firstHalf.txt";
//	const char* dataFiles = "datafiles_all.txt";
//	const char* dataFiles = "datafiles_onlyfirst.txt";
	TString MappingFileName("../../generator/mapping.dat");
	TString PedestalFileName("../../generator/pedestal-statistics.txt");
	TString ofilename("probabilities_");
	ofilename += numberOfBits;
	ofilename += ".txt";
	

	TString baseName(CurrentMacroName);
		baseName += "_";
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

	const int nFrames = 10; 
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

	// data generator
	int mode = 3;
	int nCol = 5;
	const float rate = 5.0;

	////////////////////////////////////////////////////////////////////////////////
	// Huffman
	TString huffmanTableName("TPCRawSignalDifference_HuffmanTable_");
	huffmanTableName += (mode==0||mode==2)?nCol:rate;
	huffmanTableName += "mergedCollisions";
	TString HuffmanTableNameRoot(huffmanTableName);
	HuffmanTableNameRoot += ".root";
	TString HuffmanTableNameTxt(huffmanTableName);
	HuffmanTableNameTxt += ".txt";
	const char* huffmanDecoderName="TPCRawSignalDifference";

	AliHLTHuffman* hltHuffman = NULL;
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

	// generating LL Huffman coders
    unsigned int minBits = numberOfBits;//5;
    unsigned int maxBits = numberOfBits;//15; 
    const unsigned int numDiffWords = 4; 
    unsigned int Words[numDiffWords] = {50,70,90,110};
    unsigned int numberOfHuffmans = ((maxBits-minBits)+1) * numDiffWords;
    TPC::HuffmanCoder* huffman[numberOfHuffmans];
    unsigned int count = 0;
    for (unsigned int i = minBits; i <= maxBits; i++) {
        for (unsigned int j = 0; j < numDiffWords; j++) {
            std::cout << "Generating Huffman encoder for " << i << " bits and " << Words[j] << " words" << std::endl << "\t";
            huffman[count] = new TPC::HuffmanCoder(HuffmanTableNameTxt.Data());
            if (huffman[count]) {
                if (huffman[count]->GenerateLengthLimitedHuffman(i,Words[j])) ++count;
                else {
                    delete huffman[count];
                    std::cout << "WARNING: deleted last Huffman encoder" << std::endl;
                }   
            }   
        }   
    }
	numberOfHuffmans = count;
    // done

	////////////////////////////////////////////////////////////////////////////////
	// Histograms
	TH1* hConLongWords[numberOfHuffmans+1];
	std::cout << numberOfHuffmans << std::endl;
	hConLongWords[numberOfHuffmans] = new TH1D("hConLongWordsStdHuff", "Connsecutive occurrence of long words (> 10 bit) for standard Huffman", 50, 0, 50);
	hConLongWords[numberOfHuffmans]->GetXaxis()->SetTitle("number of connsecutive words");
	hConLongWords[numberOfHuffmans]->GetYaxis()->SetTitle("counts");
	for (int i = 0; i < numberOfHuffmans; i++) {
		TString shortName("hConLongWords");
		shortName += huffman[i]->GetLLMaxCodeLength();
		shortName += "b";
		shortName += huffman[i]->GetLLNumberOfWords();
		shortName += "w";
		TString longName("Connsecutive occurrence of long words (> 10 bit) for ");
		longName += huffman[i]->GetLLMaxCodeLength();
		longName += " bit and ";
		longName += huffman[i]->GetLLNumberOfWords();
		longName += " words Huffman";
		hConLongWords[i] = new TH1D(shortName, longName, 50, 0, 50);
		hConLongWords[i]->GetXaxis()->SetTitle("number of connsecutive words");
		hConLongWords[i]->GetYaxis()->SetTitle("counts");
	}

	////////////////////////////////////////////////////////////////////////////////
	// Data generator

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
	Int_t sampleMaxSignal = 0;
	double longAfterShort[numberOfHuffmans+1];
	double shortAfterLong[numberOfHuffmans+1];
	double longAfterLong[numberOfHuffmans+1];
	double shortAfterShort[numberOfHuffmans];
	for (int i = 0; i < numberOfHuffmans+1; i++) {
		longAfterShort[i] = 0;
		shortAfterLong[i] = 0;
		longAfterLong[i] = 0;
		shortAfterShort[i] = 0;
	}
	double someSignal = 0;


	while (TimeFrameNumber++ < nFrames || nFrames < 0) {
		int result = dg->SimulateFrame();
    
		if (result < 0) {
	  		std::cout << "Simulation of timeframe failed with error code " << result << "; aborting at timeframe no " << TimeFrameNumber << std::endl;
	  		break;
		} 

		std::vector<unsigned int> chindices=dg->GetChannelIndices();
		for (std::vector<unsigned int>::iterator index = chindices.begin(); index != chindices.end(); ++index) {
			TPC::DataGenerator::ChannelDesc_t desc = dg->GetChannelDescriptor(*index);
			if (desc.ptr == NULL) continue;
			int DDLnumber = ((*index)&0xffff0000)>>16;
			int padrow = desc.padrow;
			if (DDLnumber > 71) padrow += 63;
			int size = desc.size;

                            
			bool lastWordWasLong[numberOfHuffmans+1];
			int conOcc[numberOfHuffmans+1];
			for (int i = 0; i < numberOfHuffmans+1; i++) {
				lastWordWasLong[i] = false;
				conOcc[i] = 0;
			} 
			unsigned int diff;
			int signal = 0;
			int lastSignal = -1;
			AliHLTUInt64_t codeLength;
			AliHLTUInt64_t v;
			std::string temp = "";
			int cl = 0;
			bool firstSignal = false;
			for (int i = (size - 1); i >= 0; i--) {
				signal = desc.ptr[i];
				if (signal >= signalRange || signal < 0) {
//					if (desc.ptr[i] != 65535) std::cout << "WARNING: value out of range (" << desc.ptr[i] << "), current TF: " << TimeFrameNumber << ", current index: " << *index << ", current iterator: " << i << std::endl;
					continue;
				}

//				if (desc.ptr[i] <= signalRange && desc.ptr[i+1] > signalRange && !firstSignal) { // first signal
				if (lastSignal < 0) { // first signal
					diff = signal;
				} else
					diff = signal - lastSignal;

				lastSignal = signal;
				diff += signalRange;

				if (diff < 0 || diff >= 2*signalRange) {
//					std::cout << i << " " << size << " " << desc.ptr[i] << " " << desc.ptr[i+1] << " " << signalRange << " " << diff << std::endl;

					std::cerr << "WARNING: value to encode out of range (" << diff << ")" << std::endl;
					continue;
				}
				v = diff;
				hltHuffman->Encode(v,codeLength);
				if (codeLength > 12) codeLength = 22;
				if (codeLength > 10) {
					if (lastWordWasLong[numberOfHuffmans]) {
						longAfterLong[numberOfHuffmans]++;
					} else {
						longAfterShort[numberOfHuffmans]++;
					}
					conOcc[numberOfHuffmans]++;
					lastWordWasLong[numberOfHuffmans] = true;
				} else {
					if (lastWordWasLong[numberOfHuffmans]) {
						shortAfterLong[numberOfHuffmans]++;
						hConLongWords[numberOfHuffmans]->Fill(conOcc[numberOfHuffmans]);
					} else {
						shortAfterShort[numberOfHuffmans]++;
					}
					lastWordWasLong[numberOfHuffmans] = false;
					conOcc[numberOfHuffmans] = 0;
				}

				for (int k = 0; k < numberOfHuffmans; k++) {
					if (!huffman[k]->EncodeValue(diff,&temp)) {
						cl = temp.size() + 10;
					} else {
						cl = temp.size();
					}
					if (cl > 10) {
						if (lastWordWasLong[k]) {
							longAfterLong[k]++;
						} else {
							longAfterShort[k]++;
						}
						conOcc[k]++;
						lastWordWasLong[k] = true;
					} else {
						if (lastWordWasLong[k]) {
							shortAfterLong[k]++;
							hConLongWords[k]->Fill(conOcc[k]);
						} else {
							shortAfterShort[k]++;
						}
						lastWordWasLong[k] = false;
						conOcc[k] = 0;
					}
				}
				someSignal++;
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	// clean-up

	if (dg) {
		delete dg;
	}

	if (hltHuffman) {
		delete hltHuffman;
	}

	////////////////////////////////////////////////////////////////////////////////
	// print and delete Histograms

	TFile* of = TFile::Open(targetFileName, "RECREATE");
	if (!of || of->IsZombie())
	{
		std::cerr << "can not open file " << targetFileName << std::endl;
		return;
	}

	std::ofstream output(ofilename.Data());
	if (!output.good()) {
		std::cerr << "can not open file '" << ofilename << "' for writing" << std::endl;
	}

	for (int i = 0; i < numberOfHuffmans+1; i++) {
		if (output.good()) {
			if (i == numberOfHuffmans) {
				output 
					<< 0 << "\t"
					<< 0 << "\t";
			} else {
				output 
					<< huffman[i]->GetLLMaxCodeLength() << "\t"
					<< huffman[i]->GetLLNumberOfWords() << "\t";
				delete huffman[i];
			}
			output 
				<< (double)longAfterLong[i]/(double)someSignal << "\t"
				<< (double)longAfterShort[i]/(double)someSignal << "\t"
				<< (double)shortAfterLong[i]/(double)someSignal << "\t"
				<< (double)shortAfterShort[i]/(double)someSignal << "\t"
				<< std::endl;
		}

		if (hConLongWords){
			of->cd();
			hConLongWords[i]->Print();
			hConLongWords[i]->Write();
			printHist(hConLongWords[i],plotBaseName);
			delete hConLongWords[i];
		}
	}
	if (output.good()) output.close();

	of->Close();
}

void printHist(TH1* hist, TString baseName) {
    TCanvas* cnv2 = new TCanvas("cnv2", "cnv2",1000,1000);
    cnv2->SetLeftMargin(0.11);
    cnv2->SetRightMargin(0.10);
    
    cnv2->cd();
    hist->Draw();
	cnv2->SetLogy();
//	cnv2->SetLogx();
    cnv2->SetGrid(1,1);
    cnv2->Update();
    TString printName(baseName);
    printName += hist->GetName();
    printName += ".pdf";
    cnv2->Print(printName);
    delete cnv2;
}

void printHist(TH2* hist, TString baseName) {
    TCanvas* cnv1 = new TCanvas("cnv1", "cnv1",1000,1000);
    cnv1->SetLeftMargin(0.11);
    cnv1->SetRightMargin(0.10);
    
    cnv1->cd();
    hist->Draw("COLZ1");
//    hist->Draw("lego2z");
//    cnv1->SetLogx();
//    cnv1->SetLogz();
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

