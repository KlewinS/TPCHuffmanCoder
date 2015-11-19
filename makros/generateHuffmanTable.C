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

void generateHuffmanTable(TString CurrentMacroName, float rate)
{
	////////////////////////////////////////////////////////////////////////////////
    // input / output files
//	const char* dataFiles = "datafiles_firstHalf.txt";
	const char* dataFiles = "datafiles_all.txt";
//	const char* dataFiles = "datafiles_onlyfirst.txt";
	TString MappingFileName("../../generator/mapping.dat");
	TString PedestalFileName("../../generator/pedestal-statistics.txt");
	TString VerilogLLHuffmanDecoderTable("VerilogHuffmanDecoderTable.v");
	TString VerilogLLHuffmanCodeTable("VerilogHuffmanEnocderTable.v");
	TString VerilogLLHuffmanLengthTable("VerilogHuffmanEncoderLengthTable.v");

	TString baseName(CurrentMacroName);
		baseName += "_";
		baseName.Remove(baseName.Last('.'),2);
		baseName += rate;
		baseName += "rate_";
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
    const int DDLmin = -1;
    const int DDLmax = -1;

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
	const bool enableHeader = true;

	// data generator
	int mode = 3;
	int nCol = 5;
//	const float rate = 5.0;

	////////////////////////////////////////////////////////////////////////////////
	// Huffman
	TString huffmanTableName("TPCRawSignalDifference_HuffmanTable_");
	huffmanTableName += (mode==0||mode==2)?nCol:rate;
	huffmanTableName += "mergedCollisions.root";
	const char* huffmanDecoderName="TPCRawSignalDifference";

	AliHLTHuffman* hltHuffman = NULL;
	hltHuffman = new AliHLTHuffman(huffmanDecoderName, signalBitLength+1);

	////////////////////////////////////////////////////////////////////////////////
	// Histograms
	TH1* hDiffSignals = NULL;
		hDiffSignals = new TH1D("hDiffSignals", "Raw data signal differences", 2*signalRange, 0, 2*signalRange);
		hDiffSignals->GetXaxis()->SetTitle("difference");
		hDiffSignals->GetYaxis()->SetTitle("counts");

	TH1* hSignals = NULL;
		hSignals = new TH1D("hSignals", "Raw data signal", signalRange, 0, signalRange);
		hSignals->GetXaxis()->SetTitle("signal");
		hSignals->GetYaxis()->SetTitle("counts");
		
	
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
			unsigned int diff;
			int signal = 0;
			int lastSignal = -1;
			for (int i = size-1; i >= 0; i--) {
				signal = desc.ptr[i];
				if (signal >= signalRange || signal < 0) continue;

//				if (i == size - 1) // first signal
				if (lastSignal < 0) // first signal
					diff = signal;
				else
					diff = signal - lastSignal;

				lastSignal = signal;
				diff += signalRange;

				hSignals->Fill(signal);
				hDiffSignals->Fill(diff);
				hltHuffman->AddTrainingValue(diff);
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	// clean-up

	if (dg) {
		delete dg;
	}

	////////////////////////////////////////////////////////////////////////////////
	// print and delete Histograms

	TFile* of = TFile::Open(targetFileName, "RECREATE");
	if (!of || of->IsZombie()) {
		std::cerr << "Can't open file " << targetFileName << std::endl;
	}

	if (hltHuffman) {
		hltHuffman->GenerateHuffmanTree();
		TFile* htf=TFile::Open(huffmanTableName, "RECREATE");
		if (!htf || htf->IsZombie()) {
			std::cerr << "can not open file " << htf << endl;
			return;
		}
		htf->cd();
		hltHuffman->Write();
		htf->Close();
	

		huffmanTableName.Remove(huffmanTableName.Last('.'));
		huffmanTableName += ".txt";
		std::ofstream out(huffmanTableName);
		std::streambuf *coutbuf = std::cout.rdbuf();
		std::cout.rdbuf(out.rdbuf());
		hltHuffman->Print("full");
		std::cout.rdbuf(coutbuf);
		delete hltHuffman;
	}

	TPC::HuffmanCoder* huffman = NULL;
	huffman = new TPC::HuffmanCoder(huffmanTableName.Data());
	huffman->GenerateLengthLimitedHuffman(12,30);
	huffman->WriteVerilogEncoderTable(VerilogLLHuffmanCodeTable.Data(),VerilogLLHuffmanLengthTable.Data());
	huffman->WriteVerilogDecoderTable(VerilogLLHuffmanDecoderTable.Data());
	delete huffman;

	

	if (hSignals) {
		if (of) {
			of->cd();
			hSignals->Write();
		}
		printHist(hSignals,plotBaseName);
		delete hSignals;
	}

	if (hDiffSignals) {
		if (of) {
			of->cd();
			hDiffSignals->Write();
		}
		printHist(hDiffSignals,plotBaseName);
		delete hDiffSignals;
	}

	if (of) {
		of->Close();
	}
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

