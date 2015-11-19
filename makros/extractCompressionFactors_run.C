void extractCompressionFactors_run(float rate = 5) {
	gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I.");
	gSystem->Load("../../inst/lib/libGenerator.so");
	gROOT->LoadMacro("../HuffmanCoder.cpp+");
	gROOT->LoadMacro("extractCompressionFactors.C+");
	TString CurrentMacroName(gInterpreter->GetCurrentMacroName());
	CurrentMacroName.Remove(CurrentMacroName.Last('_'));
	extractCompressionFactors(CurrentMacroName, rate);
}