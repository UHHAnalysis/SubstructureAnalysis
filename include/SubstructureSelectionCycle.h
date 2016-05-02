// Dear emacs, this is -*- c++ -*-
#ifndef SubstructureSelectionCycle_H
#define SubstructureSelectionCycle_H

// SFrame include(s):
#include "SFrameAnalysis/include/AnalysisCycle.h"
#include "SFrameAnalysis/include/Cleaner.h"
#include "SFrameTools/include/HypothesisDiscriminator.h"
#include "SFrameAnalysis/include/HypothesisHists.h"
#include "SFrameAnalysis/include/SelectionModules.h"
#include "SFrameTools/include/HypothesisStatistics.h"

#include "SFrameAnalysis/include/EventHists.h"
#include "SFrameAnalysis/include/JetHists.h"
#include "SFrameAnalysis/include/ElectronHists.h"
#include "SFrameAnalysis/include/MuonHists.h"
#include "SFrameAnalysis/include/TauHists.h"
#include "SFrameAnalysis/include/TopJetHists.h"

/**
 *  @short Selection cycle to perform 
 *         full selection for Z'->ttbar analysis
 *  @author Thomas Peiffer
 */

class SubstructureSelectionCycle : public AnalysisCycle {

public:
  /// Default constructor
  SubstructureSelectionCycle();
  /// Default destructor
  ~SubstructureSelectionCycle();

  /// Function called at the beginning of the cycle
  void BeginCycle() throw( SError );
  /// Function called at the end of the cycle
  void EndCycle() throw( SError );

  /// Function called at the beginning of a new input data
  void BeginInputData( const SInputData& ) throw( SError );
  /// Function called after finishing to process an input data
  void EndInputData  ( const SInputData& ) throw( SError );

  /// Function called after opening each new input file
  void BeginInputFile( const SInputData& ) throw( SError );

  /// Function called for every event
  void ExecuteEvent( const SInputData&, Double_t ) throw( SError );

  /// Fill control histograms
  void FillControlHists(TString postfix="");

private:
  //
  // Put all your private variables here
  //
  
  std::string m_Electron_Or_Muon_Selection;

  // flag to cut on generated mttbar
  bool m_mttgencut;

  // Flg use to reverse electron selection
  bool m_reversed_electron_selection;
 
  int m_Nbtags_max;
  int m_Nbtags_min;  

  bool m_veto_electron_trigger;
  bool m_useORTriggerWithPFJet320;

  Cleaner* m_cleaner;
  Chi2Discriminator* m_chi2discr;
  BestPossibleDiscriminator* m_bpdiscr;
  SumDeltaRDiscriminator* m_sumdrdiscr;
  CorrectMatchDiscriminator* m_cmdiscr;

  HypothesisStatistics* m_bp_chi2;
  HypothesisStatistics* m_bp_sumdr;
  HypothesisStatistics* m_cm_chi2;
  HypothesisStatistics* m_cm_sumdr; 
  HypothesisStatistics* m_cm_bp;


  TTree* m_intree; 
  bool m_rec_presel;
  bool m_gen_presel;
  bool m_gen_presel_hadron;

  /////test///
  TTree* m_tree; 
  bool m_recsel_basic;
  bool m_gensel_basic;
  ///////////////

  // Macro adding the functions for dictionary generation
  ClassDef( SubstructureSelectionCycle, 0 );

}; // class ZprimeSelectionCycle

#endif // ZprimeSelectionCycle_H

