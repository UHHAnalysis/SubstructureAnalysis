// Dear emacs, this is -*- c++ -*-
// $Id$
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;

// Add the declarations of your cycles, and any other classes for which you
// want to generate a dictionary, here. The usual format is:
//
// #pragma link C++ class MySuperClass+;
#pragma link C++ class GenJetProps+;
#pragma link C++ class SubstructureCycle+;
#pragma link C++ class SubstructurePreSelectionCycle+;
#pragma link C++ class SubstructureJetHTPreSelectionCycle+;
#pragma link C++ class SubstructureSelectionCycle+;
//#pragma link C++ class PileupCycle+;
#pragma link C++ class WSelectionCycle+;
#pragma link C++ class WPostSelectionCycle+;
		       //#pragma link C++ class WCycle+;
#pragma link C++ class WHists+;
//#pragma link C++ class testCycle+;
#pragma link C++ class SubstructurePostSelectionCycle+;
		       //#pragma link C++ class GenCycle+;
		       //#pragma link C++ class GenPreSelectionCycle+;
		       //#pragma link C++ class GenPreSelectionTightCycle+;
		       //#pragma link C++ class GenSelectionCycle+;
		       //#pragma link C++ class JetmassCycle+;
		       //#pragma link C++ class ControlHists+;
#pragma link C++ class TruthHists+;
		       //#pragma link C++ class PileupHists+;
#pragma link C++ class ComparisonHists+;
#pragma link C++ class ComparisonHistsReco+;
//#pragma link C++ class SoftDrop+;


#endif // __CINT__
