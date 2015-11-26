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
#include "../HuffmanCoder.h"
#include "AliHLTHuffman.h"
#include "../../generator/DataGenerator.h"
#include "TStyle.h"

#include <iostream>
#include <iomanip>
#include <string>

void printHist(TH1* hist, TString baseName);
void printHist(TH2* hist, TString baseName);
void printProf(TProfile* prof, TString baseName,TString axis);

void extractCompressionFactors_newTest(TString CurrentMacroName, float rateForTable, const int numberOfUncompressedWords)
{
	////////////////////////////////////////////////////////////////////////////////
    // input / output files
//	const char* dataFiles = "datafiles_firstHalf.txt";
//	const char* dataFiles = "datafiles_onlyfirst.txt";
	const char* dataFiles = "datafiles_firstTen.txt";
//	const char* dataFiles = "datafiles_secondHalf.txt";
//	const char* dataFiles = "datafiles_all.txt";
	TString HuffmanBaseFileName("TPCRawSignalDifference_HuffmanTable_");
	TString MappingFileName("../../generator/mapping.dat");
	TString PedestalFileName("../../generator/pedestal-statistics.txt");

	TString baseName(CurrentMacroName);
//		baseName.Remove(baseName.Last('.'));
		baseName += "_";
		if (rateForTable < 10) baseName += "0";
		baseName += rateForTable;
		baseName += "tableRate_";
		if (numberOfUncompressedWords < 10) baseName += "0";
		baseName += numberOfUncompressedWords;
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
	const int headerSize = 50; // bit
	const bool enableHeader = false;

	// data generator
	const int numberOfGenerators = 3;//7;
	const int mode = 3;
	const int nCol = 5;
	const float rate = 5.0;


//	const int numberOfUncompressedWords = 0;

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

    const unsigned int numDiffBits = 1;//2;
    unsigned int Bits[numDiffBits] = {10};//{10,12}; 
    const unsigned int numDiffWords = 1;//5; 
    unsigned int Words[numDiffWords] = {30};//{30,50,70,90,110};
    unsigned int numberOfHuffmans = numDiffBits*numDiffWords;
    TPC::HuffmanCoder* huffman[numberOfHuffmans];
    unsigned int count = 0;
    for (unsigned int i = 0; i < numDiffBits; i++) {
        for (unsigned int j = 0; j < numDiffWords; j++) {
            std::cout << "Generating Huffman encoder for " << Bits[i] << " bits and " << Words[j] << " words" << std::endl << "\t";
            huffman[count] = new TPC::HuffmanCoder(HuffmanTableNameTxt.Data());
            if (huffman[count]) {
				huffman[count]->SetLLRawDataMarkerMaxSize(10);
                if (huffman[count]->GenerateLengthLimitedHuffman(Bits[i],Words[j])) ++count;
                else { 
                    delete huffman[count];
                    std::cout << "WARNING: deleted last Huffman encoder" << std::endl;
                }       
            }       
        }       
    }       
    numberOfHuffmans = count;

	////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////
	// Histograms
	TH2* hCompFactorStdHuffman;
	hCompFactorStdHuffman = new TH2D("hCompFactorStdHuffman", "Compression Factor of Standard Huffman", occupancyMax*200, 0, occupancyMax, 100, 0, 10);
		hCompFactorStdHuffman->GetXaxis()->SetTitle("occupancy");
		hCompFactorStdHuffman->GetYaxis()->SetTitle("compression factor");
		hCompFactorStdHuffman->SetContour(99);
		hCompFactorStdHuffman->SetStats(0);

	TH2* hCompFactorLLHuffman[numberOfHuffmans];
	for (int i = 0; i < numberOfHuffmans; i++) {
		int bits = huffman[i]->GetLLMaxCodeLength();
		int words = huffman[i]->GetLLNumberOfWords();
		TString name("hCompFactor");
		name += bits;
		name += "b";
		name += words;
		name += "w";
		TString title("Compression Factor of ");
		title += bits;
		title += " bits and ";
		title += words;
		title += " words LL Huffman";
		hCompFactorLLHuffman[i] = new TH2D(name, title, occupancyMax*200, 0, occupancyMax, 100, 0, 10);
		hCompFactorLLHuffman[i]->GetXaxis()->SetTitle("occupancy");
		hCompFactorLLHuffman[i]->GetYaxis()->SetTitle("compression factor");
		hCompFactorLLHuffman[i]->SetContour(99);
		hCompFactorLLHuffman[i]->SetStats(0);
	}

	////////////////////////////////////////////////////////////////////////////////
	// Data generator

	TPC::DataGenerator* dg[numberOfGenerators];
	for (int i = 0; i < numberOfGenerators; i++) {
		dg[i] = new TPC::DataGenerator(mode);
		dg[i]->SetMappingFileName(MappingFileName);
		dg[i]->SetPedestalFileName(PedestalFileName);
		dg[i]->SetDDLRange(DDLmin,DDLmax);
		dg[i]->SetPadrowRange(PadrowMin,PadrowMax);
		dg[i]->InitZeroSuppression(ZSthreshold, ZSbaseline);
		dg[i]->SetApplyZeroSuppression(false);
		dg[i]->SetNormalizeChannels(false);
		dg[i]->SetApplyCommonModeEffect(false);
		if (numberOfGenerators > 1) dg[i]->Init( (mode==0||mode==2)?nCol:((i*2)+1), dataFiles);
		else dg[i]->Init( (mode==0||mode==2)?nCol:rate, dataFiles);
	}


	////////////////////////////////////////////////////////////////////////////////
	// start

	for (int k = 0; k < numberOfGenerators; k++) {
		int result = 0;
		while ((TimeFrameNumber++ < nFrames || nFrames < 0) && result >= 0) {
			result = dg[k]->SimulateFrame();
	    
			if (result < 0) {
		  		std::cout << "Simulation of timeframe failed with error code " << result << "; aborting at timeframe no " << TimeFrameNumber << std::endl;
		  		continue;
			} 
	
			std::vector<unsigned int> chindices=dg[k]->GetChannelIndices();
			for (std::vector<unsigned int>::iterator index = chindices.begin(); index != chindices.end(); ++index) {
				TPC::DataGenerator::ChannelDesc_t desc = dg[k]->GetChannelDescriptor(*index);
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
				double cf = 0;
				AliHLTUInt64_t codeLength;
				AliHLTUInt64_t v;
				int oldLength = 0;
				int newLength[numberOfHuffmans+1];
				for (int i = 0; i < numberOfHuffmans+1; i++) {
					newLength[i] = 0;
				}
				int sendUncompressed[numberOfHuffmans];
				for (int i = 0; i < numberOfHuffmans; i++) {
					sendUncompressed[i] = 0;
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
						if (sendUncompressed[j] == 0) {
							if (!huffman[j]->EncodeValue(diff,&temp)) {
								newLength[j] += 10;
								sendUncompressed[j] = numberOfUncompressedWords;
							}
							newLength[j] += temp.size();
						} else {
							newLength[j] += 10;
							sendUncompressed[j]--;
						}
					}
					oldLength += 10;
				}
	
				if (enableHeader) oldLength += headerSize;
				for (int j = 0; j < numberOfHuffmans; j++) {
					while (newLength[j] % 10 != 0) newLength[j]++;
					if (enableHeader) newLength[j] += headerSize;
					cf = (double)oldLength / (double)newLength[j];
					hCompFactorLLHuffman[j]->Fill(desc.occupancy,cf);
				}
	
				while (newLength[numberOfHuffmans] % 10 != 0) newLength[numberOfHuffmans]++;
				if (enableHeader) newLength[numberOfHuffmans] += headerSize;
				cf = (double)oldLength / (double)newLength[numberOfHuffmans];
				hCompFactorStdHuffman->Fill(desc.occupancy,cf);
			}
		}
		delete dg[k];
	}

	////////////////////////////////////////////////////////////////////////////////
	// print and delete Histograms

	TFile* of = TFile::Open(targetFileName, "RECREATE");
	if (!of || of->IsZombie())
	{
		std::cerr << "can not open file " << targetFileName << std::endl;
		return;
	}

	int markerPos = 0;
	int color = 1;
	int marker[9] = {2,5,24,25,26,27,28,30,32};
	TProfile* hProfXStdHuffman;
	if (hCompFactorStdHuffman) {
		TString nameX(hCompFactorStdHuffman->GetName());
		nameX += "_pfx";
		hProfXStdHuffman = hCompFactorStdHuffman->ProfileX(nameX);

		TString nameY(hCompFactorStdHuffman->GetName());
		nameY += "_pfy";
		TProfile* hProfY;
		hProfY = hCompFactorStdHuffman->ProfileY(nameY);
		if (of) {
			of->cd();
			hCompFactorStdHuffman->Write();
			hProfXStdHuffman->Write();
			hProfY->Write();
		}
		printHist(hCompFactorStdHuffman,plotBaseName);
		printProf(hProfXStdHuffman,plotBaseName,"X");
		printProf(hProfY,plotBaseName,"Y");
		delete hProfY;
		delete hCompFactorStdHuffman;
	}
	TProfile* hProfXLLHuffman[numberOfHuffmans];
	for (int i = 0; i < numberOfHuffmans; i++) {
		if (hCompFactorLLHuffman[i]) {
			TString nameX(hCompFactorLLHuffman[i]->GetName());
			nameX += "_pfx";
			hProfXLLHuffman[i] = hCompFactorLLHuffman[i]->ProfileX(nameX);
	
			TString nameY(hCompFactorLLHuffman[i]->GetName());
			nameY += "_pfy";
			TProfile* hProfY;
			hProfY = hCompFactorLLHuffman[i]->ProfileY(nameY);
			if (of) {
				of->cd();
				hCompFactorLLHuffman[i]->Write();
				hProfXLLHuffman[i]->Write();
				hProfY->Write();
			}
			printHist(hCompFactorLLHuffman[i],plotBaseName);
			printProf(hProfXLLHuffman[i],plotBaseName,"X");
			printProf(hProfY,plotBaseName,"Y");
			delete hProfY;
			delete hCompFactorLLHuffman[i];
		}
	}
/*
	gStyle->SetOptTitle(0);
    TCanvas* cnv4 = new TCanvas("cnv4", "cnv4",1000,1000);
    cnv4->SetLeftMargin(0.10);
    cnv4->SetRightMargin(0.10);
    cnv4->SetTopMargin(0.10);
    cnv4->SetBottomMargin(0.10);
	TPad* upperPad = new TPad("upperPad", "upperPad", .005, .35, .995, .995);
	TPad* lowerPad = new TPad("lowerPad", "lowerPad", .005, .1, .995, .35);
	upperPad->SetBottomMargin(0);
	upperPad->SetTopMargin(0.1);
	upperPad->Draw();
    
    upperPad->cd();
	TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
	hProfXStdHuffman->GetXaxis()->SetLabelFont(63);
	hProfXStdHuffman->GetXaxis()->SetLabelSize(14);
	hProfXStdHuffman->GetXaxis()->SetTitleFont(63);
	hProfXStdHuffman->GetXaxis()->SetTitleSize(20);
	hProfXStdHuffman->GetXaxis()->SetTitleOffset(1.3);
	hProfXStdHuffman->GetXaxis()->SetDecimals();
	hProfXStdHuffman->GetXaxis()->SetTitle("");
	hProfXStdHuffman->GetYaxis()->SetRangeUser(1.9,4.5);
	hProfXStdHuffman->GetYaxis()->SetTitle("Compression Factor");
	hProfXStdHuffman->GetYaxis()->SetLabelFont(63);
	hProfXStdHuffman->GetYaxis()->SetLabelSize(14);
	hProfXStdHuffman->GetYaxis()->SetTitleFont(63);
	hProfXStdHuffman->GetYaxis()->SetTitleSize(20);
	hProfXStdHuffman->GetYaxis()->SetTitleOffset(1.3);
	hProfXStdHuffman->GetYaxis()->SetDecimals();
	hProfXStdHuffman->SetMarkerStyle(marker[markerPos++]);
	hProfXStdHuffman->SetMarkerSize(1);
	hProfXStdHuffman->SetLineColor(color);
	hProfXStdHuffman->SetMarkerColor(color++);
	hProfXStdHuffman->SetStats(0);
    hProfXStdHuffman->Draw("E1");
	leg->AddEntry(hProfXStdHuffman,"std Huffman","p");
    upperPad->SetGrid(1,1);
    upperPad->Update();
	for (int i = 0; i < numberOfHuffmans; i++) {
		hProfXLLHuffman[i]->SetMarkerStyle(marker[markerPos++]);
		hProfXLLHuffman[i]->SetMarkerSize(1);
		hProfXLLHuffman[i]->SetLineColor(color);
		hProfXLLHuffman[i]->SetMarkerColor(color++);
		hProfXLLHuffman[i]->SetStats(0);
		hProfXLLHuffman[i]->Draw("E1SAME");
		TString legendEntryName("");
		legendEntryName += huffman[i]->GetLLMaxCodeLength();
		legendEntryName += " bits and ";
		legendEntryName += huffman[i]->GetLLNumberOfWords();
		legendEntryName += " words LL Huffman ";
		leg->AddEntry(hProfXLLHuffman[i],legendEntryName,"p");
	}
	upperPad->Update();
	leg->Draw("SAME");
	upperPad->Update();

	TLine* bandWithLimit = new TLine(0,2.5,occupancyMax,2.5);
	bandWithLimit->SetLineWidth(1);
	bandWithLimit->SetLineColor(kRed);
	bandWithLimit->SetLineStyle(1);
	bandWithLimit->Draw("SAME");

	upperPad->Update();
	
	cnv4->cd();
	lowerPad->SetTopMargin(0);
	lowerPad->SetBottomMargin(0.1);
	lowerPad->Draw();
	lowerPad->cd();
	TH1D* h1 = hProfXStdHuffman->ProjectionX();
	TH1D* h2 = hProfXLLHuffman[0]->ProjectionX();
	h2->Divide(h1);
	h2->GetXaxis()->SetLabelFont(63);
	h2->GetXaxis()->SetLabelSize(14);
	h2->GetXaxis()->SetTitleFont(63);
	h2->GetXaxis()->SetTitleSize(20);
	std::cout << h2->GetXaxis()->GetTitleOffset() << std::endl;
	h2->GetXaxis()->SetTitleOffset(1.3);
	std::cout << h2->GetXaxis()->GetTitleOffset() << std::endl;
	h2->GetXaxis()->SetDecimals();
	h2->GetXaxis()->SetTitle("occupancy");
	h2->GetYaxis()->SetTitle("LL Huffman / std Huffman");
	h2->GetYaxis()->SetLabelFont(63);
	h2->GetYaxis()->SetLabelSize(14);
	h2->GetYaxis()->SetTitleSize(20);
	h2->GetYaxis()->SetTitleOffset(1.3);
	h2->GetYaxis()->SetTitleFont(63);
	h2->GetYaxis()->SetNdivisions(505);
	h2->GetYaxis()->SetRangeUser(0.9,1.11);
	h2->GetYaxis()->SetDecimals();
	h2->SetMarkerStyle(hProfXLLHuffman[0]->GetMarkerStyle());
	h2->SetMarkerSize(hProfXLLHuffman[0]->GetMarkerSize());
	h2->SetLineColor(hProfXLLHuffman[0]->GetLineColor());
	h2->SetMarkerColor(hProfXLLHuffman[0]->GetMarkerColor());
	h2->SetStats(0);
	h2->Draw("E1");
    lowerPad->SetGrid(1,1);
	lowerPad->Update();
	TLine* zeroLine = new TLine(0,1,occupancyMax,1);
	zeroLine->SetLineWidth(1);
	zeroLine->SetLineColor(kRed);
	zeroLine->SetLineStyle(1);
	zeroLine->Draw("SAME");

	TH1D* h3 = NULL;
	for (int i = 1; i < numberOfHuffmans; i++) {
		h3 = hProfXLLHuffman[i]->ProjectionX();
		h3->Divide(h1);
		h3->GetYaxis()->SetTitle("std Huffman / LL Huffman");
		h3->SetMarkerStyle(hProfXLLHuffman[i]->GetMarkerStyle());
		h3->SetMarkerSize(hProfXLLHuffman[i]->GetMarkerSize());
		h3->SetLineColor(hProfXLLHuffman[i]->GetLineColor());
		h3->SetMarkerColor(hProfXLLHuffman[i]->GetMarkerColor());
		h3->SetStats(0);
		h3->Draw("E1SAME");
	}
	lowerPad->Update();
	cnv4->Update();
    TString printName(plotBaseName);
    printName += "all_pfx.pdf";
    cnv4->Print(printName);
	if (of) {
		of->cd();
		cnv4->Write();
	}

	if (occupancyMax != 0.4) {
		hProfXStdHuffman->GetXaxis()->SetRangeUser(0,0.4);
		h2->GetXaxis()->SetRangeUser(0,0.4);
		bandWithLimit->SetX2(0.4);
		zeroLine->SetX2(0.4);
		cnv4->Update();
		printName.Insert(printName.Last('_'),"_shortRange");
	    cnv4->Print(printName);
	}

    delete cnv4;
*/	delete hProfXStdHuffman;
	for (int i = 0; i < numberOfHuffmans; i++) {
		delete hProfXLLHuffman[i];
	}
//	delete bandWithLimit;
//	delete zeroLine;
//	delete leg;
	of->cd();
	of->Close();

	////////////////////////////////////////////////////////////////////////////////
	// clean-up

	if (hltHuffman) {
		delete hltHuffman;
	}

	for (int i = 0; i < numberOfHuffmans; i++) { 
		if (huffman[i]) {
			delete huffman[i];
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

