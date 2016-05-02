// Dear emacs, this is -*- c++ -*-
// $Id: SubstructurePostSelectionCycle.h,v 1.4 2013/06/12 12:40:14 peiffer Exp $
#ifndef SubstructurePostSelectionCycle_H
#define SubstructurePostSelectionCycle_H

#include "SFrameAnalysis/include/AnalysisCycle.h"
#include "SFrameAnalysis/include/Cleaner.h"

#include "include/PileupHists.h"
#include "include/TruthHists.h"
#include "include/GroomingHists.h"
#include "include/WHists.h"
#include "include/ComparisonHists.h"
#include "include/ComparisonHistsReco.h"
#include "include/JetMassTools.h"

#include "SFrameAnalysis/include/SelectionModules.h"
#include "SFrameAnalysis/include/HypothesisHists.h"
#include "SFrameAnalysis/include/TopJetHists.h"
#include "SFrameAnalysis/include/JetHists.h"
#include "SFrameAnalysis/include/EventHists.h"
#include "SFrameTools/include/TopFitCalc.h"
#include "include/SubstructureSelectionModules.h"
#include "NtupleWriter/include/JetProps.h"
#include "fastjet/tools/Filter.hh"



class SubstructurePostSelectionCycle : public AnalysisCycle {

public:
  /// Default constructor
  SubstructurePostSelectionCycle();
  /// Default destructor
  ~SubstructurePostSelectionCycle();

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

  void SetDefaults();
  void SetOutputTreeBranches();


private:
 
 bool m_reversed_electron_selection;

  std::string m_Electron_Or_Muon_Selection;

  bool m_mttgencut;

  bool m_monitoring;
  bool m_groom_recojets;
  bool m_groom_genjets;

  std::string m_matching_gen;
  std::string m_matching_reco;

  int m_Nbtags_max;
  int m_Nbtags_min;

   std::string m_dobsf;
  E_BtagType m_btagtype;
  E_BtagType x_btagtype;
  BTaggingScaleFactors* m_bsf;
 

  bool m_veto_electron_trigger;
  bool m_useORTriggerWithPFJet320;

  bool m_applyEleORJetTriggerSF;

  Cleaner* m_cleaner;
  
  // out put tree variables
  double m_weight;
  double m_lumiweight;
  double m_recweight;
  double m_genweight;
  double m_jetmass_reco;
  double m_jetmass_reco2;

  double m_jetmass_gen;
  double m_jetmass_gen2;

  Bool_t m_Nsubjets;
  int m_NPV;

  double m_genjet_pt1;
  double m_genjet_pt2;
  double m_genjet_pt3;

  double m_recojet_pt1;
  double m_recojet_pt2;
  double m_recojet_pt3;
  double m_top_pt;
  double m_hadr_top_pt;
  double m_antitop_pt;

 
  Bool_t m_recosel_300_100;
  Bool_t m_recosel_300_150;
  Bool_t m_recosel_400_100;
  Bool_t m_recosel_400_150;

  Bool_t m_gensel_300_100;
  Bool_t m_gensel_400_100;
  Bool_t m_gensel_500_100;
  Bool_t m_gensel_300_150;
  Bool_t m_gensel_400_150;
  
  Bool_t m_gensel;
  Bool_t m_gensel_v100;
  Bool_t m_gensel_v150;

  Bool_t m_gen_mass_sel;
  Bool_t m_gen_mass_sel2;
  Bool_t m_gen_mass_sel3;
  Bool_t m_gen_boosted;

  Bool_t m_reco_mass_sel;
  Bool_t m_reco_mass_sel2;
  Bool_t m_reco_boosted;

  Bool_t m_recosel;

  Bool_t m_reco_jet_veto;
  Bool_t m_gen_jet_veto;
  
  TTree *m_tree;

  //input selection tree variables
  TTree* m_intree; 
  bool m_recsel_basic;
  bool m_gensel_basic;
  bool m_rec_presel;
  bool m_gen_presel;
  bool m_gen_presel_hadron;


  //==================================selections============================== 
  Selection* BSelT, *RecoLeptonSelection, *TopJetSelection_400_150, *TopJetSelection_300_150;
  Selection *TopJetSelection_400_100, *TopJetSelection_500_150, *TopJetSelection_300_100, *WSel, *RecoJetVeto;
  Selection *GenJetSelection_500_100, *GenJetSelection_400_100 , *GenJetSelection_400_300, *GenJetSelection_400_400, 
    *GenJetSelection_300_100, * GenJetSelection_400_150, *GenJetSelection_300_150 ;
  Selection *GenJetVeto_30, *GenJetVeto_100, *GenJetVeto_150;
  Selection *mttbar_gen_selection, *HTlep_selection;

  Selection *BoostedGenSelection, *GenMassSelection, *GenMassSelectionLepton, *GenMassSelectionLeptonAndNeutrino;
 Selection *BoostedRecoSelection, *RecoMassSelection, *RecoMassSelectionLepton;

  Selection *GenMatchedSelection, *GenMatchedSelectionHadron, *RecoMatchedSelection; 

  //==============================base hists===================================
  //full hadronic gen
  BaseHists *Hists_hadronic_gen, *Hists_hadronic_gen_mass;

  //gen level historgams
  BaseHists *Hists_gen_presel;
  BaseHists *Hists_gen_noVeto_no2nd, *Hists_gen_noVeto_no2nd_300, *Hists_gen_noVeto_no2nd_500;
  BaseHists *Hists_gen_noVeto, *Hists_gen, *Hists_gen_dr, *Hists_gen_mass, *Hists_gen_mass2,
    *Hists_gen_mass3; 
  BaseHists *Hists_gen_recgen, *Hists_gen_recgen_NOdr; 
  BaseHists *Hists_gen_rec_sel, *Hists_gen_rec_notgen_sel, *Hists_reco_rec_notgen_sel_dilept, *Hists_gen_gen_notrec_sel;
  BaseHists *Hists_gen_v50, *Hists_gen_v100;

  //reco level histograms
  BaseHists *Hists_reco_basic, *Hists_reco_noVeto_no2nd, *Hists_reco_noVeto;
  BaseHists *Hists_reco, *Hists_reco_mass, *Hists_reco_mass2; 
  BaseHists *Hists_reco_mass_highPU, *Hists_reco_mass_lowPU;
  BaseHists *Hists_reco_mass_300to400,  *Hists_reco_mass_400to500, *Hists_reco_mass_500toInf;
  BaseHists *Hists_reco_recgen, *Hists_reco_recgen_NOdr;

  BaseHists *Hists_signal, *Hists_SideBand1, *Hists_SideBand2, *Hists_SideBand3, *Hists_SideBand4, *Hists_SideBand5;

  //side studeis
  BaseHists *Hists_grooming, *Hists_W, *Hists_truth;

  //event historgrams
  BaseHists *Hists_event, *Hists_jet, *Hists_event_basic, *Hists_jet_basic;

  // Macro adding the functions for dictionary generation
  ClassDef( SubstructurePostSelectionCycle, 0 );

}; // class SubstructurePostSelectionCycle

#endif // SubstructurePostSelectionCycle_H

