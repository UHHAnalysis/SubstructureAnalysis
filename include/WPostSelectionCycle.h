// Dear emacs, this is -*- c++ -*-
// $Id: WPostSelectionCycle.h,v 1.4 2013/06/12 12:40:14 peiffer Exp $
#ifndef WPostSelectionCycle_H
#define WPostSelectionCycle_H

#include "SFrameAnalysis/include/AnalysisCycle.h"
#include "SFrameAnalysis/include/Cleaner.h"
#include "SFrameTools/include/HypothesisDiscriminator.h"
#include "SFrameTools/include/HypothesisStatistics.h"


#include "SFrameAnalysis/include/EventHists.h"
#include "SFrameAnalysis/include/JetHists.h"
#include "SFrameAnalysis/include/ElectronHists.h"
#include "SFrameAnalysis/include/MuonHists.h"
#include "SFrameAnalysis/include/TopJetHists.h"
/**
 *   @short W of an analysis cycle
 *
 *          This is an example of an analysis cycle. It can be used
 *          as a template for writing your own analysis. Also should
 *          be used for quick cross checks of the system setup.
 *
 *  @author Roman Kogler
 *  @version $Revision: 1.4 $
 */

class WPostSelectionCycle : public AnalysisCycle {

public:
  /// Default constructor
  WPostSelectionCycle();
  /// Default destructor
  ~WPostSelectionCycle();

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
 bool m_reversed_electron_selection;

  std::string m_Electron_Or_Muon_Selection;

 bool m_mttgencut;

  bool m_monitoring;

  int m_Nbtags_max;
  int m_Nbtags_min;

   std::string m_dobsf;
  E_BtagType m_btagtype;
  E_BtagType x_btagtype;
  BTaggingScaleFactors* m_bsf;
 

  bool m_veto_electron_trigger;
  bool m_useORTriggerWithPFJet320;

  std::string m_matching; 

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
  

  Selection* BSelT, * NoBSel, *TopSel, *TopSel2, *chi2_selection,* BSelM, * BSelL, *BSelT2;

  // Macro adding the functions for dictionary generation
  ClassDef( WPostSelectionCycle, 0 );

}; // class WPostSelectionCycle

#endif // WPostSelectionCycle_H

