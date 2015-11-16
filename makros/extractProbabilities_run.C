void extractProbabilities_run(int numberOfBits = 10) {
	gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I.");
	gSystem->Load("../../inst/lib/libGenerator.so");
	gROOT->LoadMacro("../HuffmanCoder.cpp+");
	gROOT->LoadMacro("extractProbabilities.C+");
	TString CurrentMacroName(gInterpreter->GetCurrentMacroName());
	CurrentMacroName.Remove(CurrentMacroName.Last('_'));
	extractProbabilities(CurrentMacroName,numberOfBits);
}
