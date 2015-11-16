void generateHuffmanTable_run(float rate = 5) {
	gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I.");
	gSystem->Load("../../inst/lib/libGenerator.so");
	gROOT->LoadMacro("../HuffmanCoder.cpp+");
	gROOT->LoadMacro("generateHuffmanTable.C+");
	TString CurrentMacroName(gInterpreter->GetCurrentMacroName());
	CurrentMacroName.Remove(CurrentMacroName.Last('_'));
	generateHuffmanTable(CurrentMacroName, rate);
}
