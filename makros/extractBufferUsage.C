#if defined(__CINT__) && !defined(__MAKECINT__)
{
	gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I.");
	TString macroname=gInterpreter->GetCurrentMacroName();
	macroname+="++";
	gROOT->LoadMacro("../HuffmanCoder.cpp++");
	gROOT->LoadMacro(macroname);
	
	extractBufferUsage();
}
#else

#include "AliRawReader.h"
#include "AliAltroRawStreamV3.h"
#include "TInterpreter.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TList.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TF2.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "AliHLTHuffman.h"

#include <iostream>
#include <iomanip>
#include <string>
#include "../HuffmanCoder.h"

void printHist(TH2* hist, TString baseName) {
	TCanvas* cnv1 = new TCanvas("cnv1", "cnv1",1000,1000);
	cnv1->SetLeftMargin(0.11);
	cnv1->SetRightMargin(0.10);

	cnv1->cd();
	hist->Draw("lego2z");
	cnv1->SetPhi(70);
	cnv1->Update();
	TPaletteAxis *palette = NULL;
	palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
	palette->SetY1NDC(0.4);
	cnv1->Update();
	TString printName(baseName);
	printName += hist->GetName();
	printName += "_";
	printName += hist->GetXaxis()->GetXmin()+0.5;
	printName += "-";
	printName += hist->GetXaxis()->GetXmax()-0.5;
	printName += "bits_";
	printName += hist->GetYaxis()->GetXmin()+5;
	printName += "-";
	printName += hist->GetYaxis()->GetXmax()-5;
	printName += "words.pdf";
	cnv1->Print(printName);
	delete cnv1;
}

void printHist(TH1* hist, TString baseName) {
	TCanvas* cnv2 = new TCanvas("cnv2", "cnv2",1000,1000);
	cnv2->SetLeftMargin(0.11);
	cnv2->SetRightMargin(0.10);

	cnv2->cd();
	hist->Draw();
	cnv2->SetGrid(1,1);
	cnv2->Update();
	TString printName(baseName);
	printName += hist->GetName();
	printName += "_";
	printName += hist->GetXaxis()->GetXmin()+5;
	printName += "-";
	printName += hist->GetXaxis()->GetXmax()-5;
	printName += "words.pdf";
	cnv2->Print(printName);
	delete cnv2;
}

void extractBufferUsage()
{

	TString baseName(gInterpreter->GetCurrentMacroName());
	baseName.Remove(baseName.Last('.'));
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

	int huffmanRange = 2048;
	int signalRange = huffmanRange/2;
	int binMargin = 50;
	int nBins = 2*(signalRange+binMargin)+1;

	TString HuffmanTableName("TPCRawSignalDifference_HuffmanTable_5mergedCollisions");
	TString HuffmanTableNameRoot(HuffmanTableName);
	HuffmanTableNameRoot += ".root";
	TString HuffmanTableNameTxt(HuffmanTableName);
	HuffmanTableNameTxt += ".txt";

	// loading singal spectrum from file
	TH1* hSignalDiffs = NULL;
	TString signalDiffFile("generateHuffmanTable_5rate_summary.root");
	TFile* signalDiffs = TFile::Open(signalDiffFile);
	if (!signalDiffs || signalDiffs->IsZombie()) {
		std::cerr << "ERROR: Cant' open file " << signalDiffFile << std::endl;
		return;
	}
	TObject* sigObject = NULL;
	TString histDiffName("hDiffSignals");
	signalDiffs->GetObject(histDiffName,sigObject);
	if (sigObject==NULL) {
		std::cerr << "ERROR: Can't load signal Differences object " << histDiffName << " from file " << signalDiffFile << std::endl;
		return;
	}
	hSignalDiffs = (TH1D*)sigObject;
	// done


	// loading standard Huffman from file
	AliHLTHuffman* hltHuffman=NULL;
	TFile* htf = TFile::Open(HuffmanTableNameRoot);
	if (!htf || htf->IsZombie()) {
		std::cerr << "ERROR: Can't open file " << HuffmanTableNameRoot << std::endl;
		return;
	}
	TObject* obj = NULL;
	htf->GetObject("TPCRawSignalDifference",obj);
	if (obj==NULL) {
		std::cerr << "ERROR: Can't load Huffman decoder object " << "TPCRawSignalDifference" << "from file " << HuffmanTableNameRoot << std::endl;
		return;
	}
	hltHuffman = (AliHLTHuffman*)obj;
	// done

	// generating LL Huffman coders
	const unsigned int numDiffBits = 2;
	unsigned int Bits[numDiffBits] = {10,12};
	const unsigned int numDiffWords = 4;
	unsigned int Words[numDiffWords] = {50,70,90,110};
	unsigned int numberOfHuffmans = numDiffBits * numDiffWords;
	TPC::HuffmanCoder* huffman[numberOfHuffmans];
	unsigned int count = 0;
	for (unsigned int i = 0; i < numDiffBits; i++) {
		for (unsigned int j = 0; j < numDiffWords; j++) {
			std::cout << "Generating Huffman encoder for " << Bits[i] << " bits and " << Words[j] << " words" << std::endl << "\t";
			huffman[count] = new TPC::HuffmanCoder(HuffmanTableNameTxt.Data());
			if (huffman[count]) {
				if (huffman[count]->GenerateLengthLimitedHuffman(Bits[i],Words[j]))	++count;
				else {
					delete huffman[count];
					std::cout << "WARNING: deleted last Huffman encoder" << std::endl;
				}
			}
		}
	}
	numberOfHuffmans = count;
	// done

	const char* probFilename = "probabilities.txt";
	std::ifstream input(probFilename);
	if (!input.good()) return;
	std::cout << "Reading probabilities from file " << probFilename << std::endl;
	int b = -1;
	int w = -1;
	double lAl = -1;
	double lAs = -1;
	double sAl = -1;
	double sAs = -1;
	int i = 0;
	double longAfterLong[numberOfHuffmans+1];
	double shortAfterLong[numberOfHuffmans+1];
	double longAfterShort[numberOfHuffmans+1];
	double shortAfterShort[numberOfHuffmans+1];
	do {
		input >> b;
		input >> w;
		input >> lAl;
		input >> lAs;
		input >> sAl;
		input >> sAs;
		if (b != 0) { 
			if (b != huffman[i]->GetLLMaxCodeLength() || w != huffman[i]->GetLLNumberOfWords()) {
				std::cerr << "ERROR: probability file (" << probFilename << ") does not fit to current configuration." << std::endl;
				return;
			}
//			std::cout << i << " " << b << " " << w << " " << lAl << " " << lAs << " " << sAl << " " << sAs << std::endl;
		}
//		else {
//			std::cout << i << " " << "std Huffman " << b << " " << w << " " << lAl << " " << lAs << " " << sAl << " " << sAs << std::endl;
//		}
		longAfterLong[i] = lAl;
		shortAfterLong[i] = sAl;
		longAfterShort[i] = lAs;
		shortAfterShort[i] = sAs;
/*		if (b != 0)
			std::cout << i << " " << huffman[i]->GetLLHuffmanMaxLength() << " " << huffman[i]->GetLLHuffmanMaxWords() << " " << longAfterShort[i] << " " << longAfterLong[i] << std::endl;
		else
			std::cout << i << " " << "std Huffman"  << " " << longAfterShort[i] << " " << longAfterLong[i] << std::endl;
*/
		i++;
	} while (input.good());
//	std::cout << numberOfHuffmans << std::endl;


	TFile* of = TFile::Open(targetFileName, "RECREATE");
	if (!of || of->IsZombie()) {
		std::cerr << "Can't open file " << targetFileName << std::endl;
		return;
	}
	


	// calculating probabilities of signals
	double maxProb = pow(10,-17);
	double entries = 0;
	double data = 0;
	double occurence[huffmanRange];
	double factor = 0; 
	for (int i = 0; i < huffmanRange; i++) {
//		std::cout << hSignalDiffs->GetXaxis()->FindBin(i-1024) << "\t" << i << "\t" << hSignalDiffs->GetBinContent(hSignalDiffs->GetXaxis()->FindBin(i-1024)) << std::endl;
		occurence[i] = hSignalDiffs->GetBinContent(hSignalDiffs->GetXaxis()->FindBin(i));
		entries += occurence[i];
//		std::cout << i << " " << occurence[i] << " " << entries << std::endl;
//		std::cout << entries << std::endl;
	}
	double probabilities[huffmanRange];
	for (int i = 0; i < huffmanRange; i++) {
		probabilities[i] = (double)occurence[i]/(double)entries;
	}

/*	
	double probTot = 0;
	for (int i = 0; i < huffmanRange; i++) {
		std::cout << probabilities[i] << std::endl;
		probTot += probabilities[i];
	}
	std::cout << "Total probability = " << probTot << std::endl;
*/
	// done

//	AliHLTUInt64_t codeLength;
	double probBigger10Standard = longAfterLong[numberOfHuffmans]+longAfterShort[numberOfHuffmans];
//	double probSmaller10Standard = 0;
//	for (int i = 0; i < huffmanRange; i++) {
//		AliHLTUInt64_t v = i;
//		hltHuffman->Encode(v,codeLength);
//		if (codeLength > 12) codeLength = 12 + 10;
//		if (codeLength > 10) probBigger10Standard += probabilities[i];
//		else probSmaller10Standard += probabilities[i];
//	}
//	std::cout << probBigger10Standard*100. << " %" << std::endl;
//	std::cout << probSmaller10Standard*100. << " %" << std::endl;
//	std::cout << longAfterShort[numberOfHuffmans] << " " << 0.00334189 << std::endl;
//	std::cout << longAfterLong[numberOfHuffmans] << " " << 0.00454116 << std::endl;
	double prob = longAfterShort[numberOfHuffmans];//probBigger10Standard;
	int countsStd = 1;
//	while (prob * (shortAfterLong[numberOfHuffmans]/ (longAfterLong[numberOfHuffmans] + shortAfterLong[numberOfHuffmans])) > maxProb) {
	while (prob > maxProb) {
		prob *= longAfterLong[numberOfHuffmans] / (longAfterLong[numberOfHuffmans] + shortAfterLong[numberOfHuffmans]);//probBigger10Standard;
		++countsStd;
		std::cout << countsStd << " " << prob << " " << maxProb << std::endl;
	}

	int bufferSizeStd = 9;
	for (int i = 0 ; i < countsStd-1; i++) {
		bufferSizeStd += 22;
		bufferSizeStd -= 10;
	}
	bufferSizeStd += 22;
	std::cout << bufferSizeStd << std::endl;

	TH2* hAllHuffProb;
	hAllHuffProb = new TH2D(
			"hAllHuffProb", 
			"Probability of code longer then 10 bits of LL Huffman",
			Bits[numDiffBits-1]-Bits[0]+1,
			Bits[0]-0.5,
			Bits[numDiffBits-1]+0.5,
			(Words[numDiffWords-1]-Words[0])/10 +1, 
			Words[0]-5, 
			Words[numDiffWords-1]+5);
		hAllHuffProb->GetXaxis()->SetTitle("code length [bits]");
		hAllHuffProb->GetXaxis()->SetTitleOffset(2.0);
		hAllHuffProb->GetYaxis()->SetTitle("encoded words");
		hAllHuffProb->GetYaxis()->SetTitleOffset(1.4);
		hAllHuffProb->GetZaxis()->SetTitle("prob LL Huffman / prob std Huffman");
		hAllHuffProb->GetZaxis()->SetTitleOffset(1.6);
		hAllHuffProb->GetZaxis()->SetRangeUser(0,20./18.);

	TH1* h12bHuffProb;
	h12bHuffProb = new TH1D(
			"h12bHuffProb", 
			"Probability of code longer then 10 bits of 12bit LL Huffman",
			(Words[numDiffWords-1]-Words[0])/10 +1, 
			Words[0]-5, 
			Words[numDiffWords-1]+5);
		h12bHuffProb->GetXaxis()->SetTitle("encoded words");
		h12bHuffProb->GetXaxis()->SetTitleOffset(1.4);
		h12bHuffProb->GetYaxis()->SetTitle("prob LL Huffman / prob std Huffman");
		h12bHuffProb->GetYaxis()->SetTitleOffset(1.6);
		h12bHuffProb->GetYaxis()->SetRangeUser(0,4.2);

	TH1* h12bHuffProb_becOfShortCode;
	h12bHuffProb_becOfShortCode = new TH1D(
			"h12bHuffProb_becOfShortCode", 
			"Probability because of short code",
			(Words[numDiffWords-1]-Words[0])/10 +1, 
			Words[0]-5, 
			Words[numDiffWords-1]+5);
		h12bHuffProb_becOfShortCode->GetXaxis()->SetTitle("encoded words");
		h12bHuffProb_becOfShortCode->GetXaxis()->SetTitleOffset(1.4);
		h12bHuffProb_becOfShortCode->GetYaxis()->SetTitle("prob LL Huffman / prob std Huffman");
		h12bHuffProb_becOfShortCode->GetYaxis()->SetTitleOffset(1.6);
		h12bHuffProb_becOfShortCode->GetYaxis()->SetRangeUser(0,3);


	TH1* h12bHuffProb_becOfLongCode;
	h12bHuffProb_becOfLongCode = new TH1D(
			"h12bHuffProb_becOfLongCode", 
			"Probability because of long code",
			(Words[numDiffWords-1]-Words[0])/10 +1, 
			Words[0]-5, 
			Words[numDiffWords-1]+5);
		h12bHuffProb_becOfLongCode->GetXaxis()->SetTitle("encoded words");
		h12bHuffProb_becOfLongCode->GetXaxis()->SetTitleOffset(1.4);
		h12bHuffProb_becOfLongCode->GetYaxis()->SetTitle("prob LL Huffman / prob std Huffman");
		h12bHuffProb_becOfLongCode->GetYaxis()->SetTitleOffset(1.6);
		h12bHuffProb_becOfLongCode->GetYaxis()->SetRangeUser(0,3);


	TH2* hAllHuffMarkLen;
	hAllHuffMarkLen = new TH2D(
			"hAllHuffMarkLen", 
			"Length of marker for orig. data LL Huffman",
			Bits[numDiffBits-1]-Bits[0]+1,
			Bits[0]-0.5,
			Bits[numDiffBits-1]+0.5,
			(Words[numDiffWords-1]-Words[0])/10 +1, 
			Words[0]-5, 
			Words[numDiffWords-1]+5);
		hAllHuffMarkLen->GetXaxis()->SetTitle("code length [bits]");
		hAllHuffMarkLen->GetXaxis()->SetTitleOffset(2.0);
		hAllHuffMarkLen->GetYaxis()->SetTitle("encoded words");
		hAllHuffMarkLen->GetYaxis()->SetTitleOffset(1.4);
		hAllHuffMarkLen->GetZaxis()->SetTitle("marker length [bits]");
		hAllHuffMarkLen->GetZaxis()->SetTitleOffset(1.6);
//		hAllHuffMarkLen->GetZaxis()->SetRangeUser(0,20./18.);

	TH1* h12bHuffMarkLen;
	h12bHuffMarkLen = new TH1D(
			"h12bHuffMarkLen", 
			"Length of marker for orig. data of 12bit LL Huffman",
			(Words[numDiffWords-1]-Words[0])/10 +1, 
			Words[0]-5, 
			Words[numDiffWords-1]+5);
		h12bHuffMarkLen->GetXaxis()->SetTitle("encoded words");
		h12bHuffMarkLen->GetXaxis()->SetTitleOffset(1.4);
		h12bHuffMarkLen->GetYaxis()->SetTitle("marker length [bits]");
		h12bHuffMarkLen->GetYaxis()->SetTitleOffset(1.6);
//		h12bHuffMarkLen->GetYaxis()->SetRangeUser(0.2,1.8);

	TH2* hAllHuffConOcc;
	hAllHuffConOcc = new TH2D(
			"hAllHuffConOcc", 
			"Allowed consecutive occurrence of max bit length of LL Huffman",
			Bits[numDiffBits-1]-Bits[0]+1,
			Bits[0]-0.5,
			Bits[numDiffBits-1]+0.5,
			(Words[numDiffWords-1]-Words[0])/10 +1, 
			Words[0]-5, 
			Words[numDiffWords-1]+5);
		hAllHuffConOcc->GetXaxis()->SetTitle("code length [bits]");
		hAllHuffConOcc->GetXaxis()->SetTitleOffset(2.0);
		hAllHuffConOcc->GetYaxis()->SetTitle("encoded words");
		hAllHuffConOcc->GetYaxis()->SetTitleOffset(1.4);
		hAllHuffConOcc->GetZaxis()->SetTitle("occ LL Huffman / occ std Huffman");
		hAllHuffConOcc->GetZaxis()->SetTitleOffset(1.6);
		hAllHuffConOcc->GetZaxis()->SetRangeUser(0,20./18.);

	TH1* h12bHuffConOcc;
	h12bHuffConOcc = new TH1D(
			"h12bHuffConOcc", 
			"Allowed consecutive occurrence of max bit length of 12bit LL Huffman",
			(Words[numDiffWords-1]-Words[0])/10 +1, 
			Words[0]-5, 
			Words[numDiffWords-1]+5);
		h12bHuffConOcc->GetXaxis()->SetTitle("encoded words");
		h12bHuffConOcc->GetXaxis()->SetTitleOffset(1.4);
		h12bHuffConOcc->GetYaxis()->SetTitle("occ LL Huffman / occ std Huffman");
		h12bHuffConOcc->GetYaxis()->SetTitleOffset(1.6);
		h12bHuffConOcc->GetYaxis()->SetRangeUser(0.8,1.2);


	TH2* hAllHuffBuffer;
	hAllHuffBuffer = new TH2D(
			"hAllHuffBuffer", 
			"Buffer size of LL Huffman",
			Bits[numDiffBits-1]-Bits[0]+1,
			Bits[0]-0.5,
			Bits[numDiffBits-1]+0.5,
			(Words[numDiffWords-1]-Words[0])/10 +1, 
			Words[0]-5, 
			Words[numDiffWords-1]+5);
		hAllHuffBuffer->GetXaxis()->SetTitle("code length [bits]");
		hAllHuffBuffer->GetXaxis()->SetTitleOffset(2.0);
		hAllHuffBuffer->GetYaxis()->SetTitle("encoded words");
		hAllHuffBuffer->GetYaxis()->SetTitleOffset(1.4);
		hAllHuffBuffer->GetZaxis()->SetTitle("buffer size LL Huffman / buffer size std Huffman");
		hAllHuffBuffer->GetZaxis()->SetTitleOffset(1.6);
		hAllHuffBuffer->GetZaxis()->SetRangeUser(0,20./18.);

	TH1* h12bHuffBuffer;
	h12bHuffBuffer = new TH1D(
			"h12bHuffBuffer", 
			"Buffer size of LL Huffman",
			(Words[numDiffWords-1]-Words[0])/10 +1, 
			Words[0]-5, 
			Words[numDiffWords-1]+5);
		h12bHuffBuffer->GetXaxis()->SetTitle("encoded words");
		h12bHuffBuffer->GetXaxis()->SetTitleOffset(1.4);
		h12bHuffBuffer->GetYaxis()->SetTitle("buffer size LL Huffman / buffer size std Huffman");
		h12bHuffBuffer->GetYaxis()->SetTitleOffset(1.6);
		h12bHuffBuffer->GetYaxis()->SetRangeUser(0.,1.5);

	TH1* h10bHuffBuffer;
	h10bHuffBuffer = new TH1D(
			"h10bHuffBuffer", 
			"Buffer size of LL Huffman",
			(Words[numDiffWords-1]-Words[0])/10 +1, 
			Words[0]-5, 
			Words[numDiffWords-1]+5);
		h10bHuffBuffer->GetXaxis()->SetTitle("encoded words");
		h10bHuffBuffer->GetXaxis()->SetTitleOffset(1.4);
		h10bHuffBuffer->GetYaxis()->SetTitle("buffer size LL Huffman / buffer size std Huffman");
		h10bHuffBuffer->GetYaxis()->SetTitleOffset(1.6);
		h10bHuffBuffer->GetYaxis()->SetRangeUser(0.,1.5);




	for (unsigned int i = 0; i < numberOfHuffmans; i++) {
		int bits = huffman[i]->GetLLMaxCodeLength();
		int words = huffman[i]->GetLLNumberOfWords();

		hAllHuffMarkLen->SetBinContent(
				hAllHuffMarkLen->GetXaxis()->FindBin(bits),
				hAllHuffMarkLen->GetYaxis()->FindBin(words),
				huffman[i]->GetLLMarkerLength());
		if (bits==12) h12bHuffMarkLen->SetBinContent(
				h12bHuffMarkLen->GetXaxis()->FindBin(words),
				huffman[i]->GetLLMarkerLength());

		double probBigger10 = longAfterLong[i]+longAfterShort[i];
//		double probBigger10_becOfShortCode = 0;
//		double probBigger10_becOfLongCode = 0;
//		for (int j = 0; j < huffmanRange; j++) {
//			std::string temp = "";
//			int codeLength = 0;
//			if (!huffman[i]->findLLHuffmanCode(j,&temp)){
//				codeLength = temp.size() + 10;
//				if (codeLength > 10) probBigger10_becOfLongCode += probabilities[j];
//			} else {
//				codeLength = temp.size();
//				if (codeLength > 10) probBigger10_becOfShortCode += probabilities[j];
//			}
//			if (codeLength > 10) probBigger10 += probabilities[j];
//		}
//
//		std::cout << bits 
//			<< "\t"	<< words
//			<< "\t" << probBigger10_1
//			<< "\t"	<< probBigger10_2
//			<< std::endl;

		hAllHuffProb->SetBinContent(hAllHuffProb->GetXaxis()->FindBin(bits),hAllHuffProb->GetYaxis()->FindBin(words),probBigger10/probBigger10Standard);
		if(bits==12) {
			h12bHuffProb->SetBinContent(h12bHuffProb->GetXaxis()->FindBin(words),probBigger10/probBigger10Standard);
//			h12bHuffProb_becOfShortCode->SetBinContent(h12bHuffProb_becOfShortCode->GetXaxis()->FindBin(words),probBigger10_becOfShortCode/probBigger10Standard);
//			h12bHuffProb_becOfLongCode->SetBinContent(h12bHuffProb_becOfLongCode->GetXaxis()->FindBin(words),probBigger10_becOfLongCode/probBigger10Standard);
		}

		double prob = longAfterShort[i];
		int counts = 1;
//		std::cout << probBigger10*100. << " %" << std::endl;
//		while (prob * (shortAfterLong[i] / (longAfterLong[i] + shortAfterLong[i])) > maxProb) {
		while (prob > maxProb) {
			prob *= longAfterLong[i] / (longAfterLong[i] + shortAfterLong[i]);
			++counts;
		}
//		std::cout << bits << "\t" << words << "\t" << counts << "\t" << prob << "\t" << maxProb << std::endl;
//		std::cout << std::endl;
		hAllHuffConOcc->SetBinContent(hAllHuffConOcc->GetXaxis()->FindBin(bits),hAllHuffProb->GetYaxis()->FindBin(words),(double)counts/(double)countsStd);
		if(bits==12) h12bHuffConOcc->SetBinContent(h12bHuffConOcc->GetXaxis()->FindBin(words),(double)counts/(double)countsStd);

		int bufferSize = 9;
		for (int j = 0 ; j < counts-1; j++) {
			bufferSize += 10 + huffman[i]->GetLLMarkerLength();
			bufferSize -= 10;
		}
		bufferSize += 10 + huffman[i]->GetLLMarkerLength();
		std::cout << bits << " " << words << "\t" << bufferSize << std::endl;
		hAllHuffBuffer->SetBinContent(hAllHuffBuffer->GetXaxis()->FindBin(bits),hAllHuffBuffer->GetYaxis()->FindBin(words),(double)bufferSize/(double)bufferSizeStd);
		if(bits==12) h12bHuffBuffer->SetBinContent(h12bHuffBuffer->GetXaxis()->FindBin(words),(double)bufferSize/(double)bufferSizeStd);
		if(bits==10) h10bHuffBuffer->SetBinContent(h10bHuffBuffer->GetXaxis()->FindBin(words),(double)bufferSize/(double)bufferSizeStd);
	}

	of->cd();
/*	if (hHltHuffCL) {
		hHltHuffCL->Print();
		hHltHuffCL->Write();
		delete hHltHuffCL;
	}
*/
	gStyle->SetOptStat(0);
	if (hAllHuffProb) {
		if (of) {
			of->cd();
			hAllHuffProb->Write();
		}
		printHist(hAllHuffProb,plotBaseName);
		delete hAllHuffProb;
	}

	if (hAllHuffConOcc) {
		if (of) {
			of->cd();
			hAllHuffConOcc->Write();
		}
		printHist(hAllHuffConOcc,plotBaseName);
		delete hAllHuffConOcc;
	}

	if (hAllHuffBuffer) {
		if (of) {
			of->cd();
			hAllHuffBuffer->Write();
		}
		printHist(hAllHuffBuffer,plotBaseName);
		delete hAllHuffBuffer;
	}

	if (hAllHuffMarkLen) {
		if (of) {
			of->cd();
			hAllHuffMarkLen->Write();
		}
		printHist(hAllHuffMarkLen,plotBaseName);
		delete hAllHuffMarkLen;
	}

	if (h12bHuffProb) {
		if (of) {
			of->cd();
			h12bHuffProb->Write();
		}
		printHist(h12bHuffProb,plotBaseName);
		delete h12bHuffProb;
	}

	if (h12bHuffConOcc) {
		if (of) {
			of->cd();
			h12bHuffConOcc->Write();
		}
		printHist(h12bHuffConOcc,plotBaseName);
		delete h12bHuffConOcc;
	}

	if (h12bHuffBuffer) {
		if (of) {
			of->cd();
			h12bHuffBuffer->Write();
		}
		printHist(h12bHuffBuffer, plotBaseName);
		delete h12bHuffBuffer;
	}

	if (h10bHuffBuffer) {
		if (of) {
			of->cd();
			h10bHuffBuffer->Write();
		}
		printHist(h10bHuffBuffer, plotBaseName);
		delete h10bHuffBuffer;
	}

	if (h12bHuffMarkLen) {
		if (of) {
			of->cd();
			h12bHuffMarkLen->Write();
		}
		printHist(h12bHuffMarkLen,plotBaseName);
		delete h12bHuffMarkLen;
	}

	if (h12bHuffProb_becOfShortCode && h12bHuffProb_becOfLongCode) {
		THStack *hs = new THStack("hs","Probability of code longer then 10 bits of 12bit LL Huffman");
		h12bHuffProb_becOfShortCode->SetFillColor(kBlue);
		h12bHuffProb_becOfLongCode->SetFillColor(kGreen);
		hs->Add(h12bHuffProb_becOfLongCode);
		hs->Add(h12bHuffProb_becOfShortCode);
		if (of) {
			of->cd();
			hs->Write();
			h12bHuffProb_becOfShortCode->Write();
			h12bHuffProb_becOfLongCode->Write();
		}
		TCanvas* cnv3 = new TCanvas("cnv3", "cnv1",1000,1000);
		cnv3->SetLeftMargin(0.11);
		cnv3->SetRightMargin(0.10);
		cnv3->cd();
		hs->Draw();
		hs->SetMinimum(0);
		hs->SetMaximum(4);
		hs->Draw();
		cnv3->SetGrid(1,1);
		cnv3->Update();
		cnv3->BuildLegend();
		cnv3->Update();
		TString printName(plotBaseName);
		printName += "h12bHuffProb_ShortLong_";
		printName += h12bHuffProb_becOfShortCode->GetXaxis()->GetXmin()+5;
		printName += "-";
		printName += h12bHuffProb_becOfShortCode->GetXaxis()->GetXmax()-5;
		printName += "words.pdf";
		cnv3->Print(printName);
		delete cnv3;
		delete hs;
		delete h12bHuffProb_becOfShortCode;
		delete h12bHuffProb_becOfLongCode;
	}



	for (unsigned int i = 0; i < numberOfHuffmans; i++) {
		if (huffman[i]) delete huffman[i];
		huffman[i] = NULL;
	}

	of->Close();
}

#endif
