void extractCompressionFactors_newTest_run(int number = 0, float rate = 5) {
	gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I.");
	gSystem->Load("../../inst/lib/libGenerator.so");
	gROOT->LoadMacro("../HuffmanCoder.cpp+");
	gROOT->LoadMacro("extractCompressionFactors_newTest.C+");
	TString CurrentMacroName(gInterpreter->GetCurrentMacroName());
	CurrentMacroName.Remove(CurrentMacroName.Last('_'));
	extractCompressionFactors_newTest(CurrentMacroName, rate, number);
}
