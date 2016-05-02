// $Id: WCycle.cxx,v 1.10 2012/12/07 14:21:51 peiffer Exp $

#include <iostream>

using namespace std;

#include "include/WPostSelectionCycle.h"

#include "SFrameAnalysis/include/SelectionModules.h"
#include "SFrameAnalysis/include/HypothesisHists.h"
#include "SFrameAnalysis/include/TopJetHists.h"
#include "SFrameAnalysis/include/JetHists.h"
#include "SubstructureAnalysis/include/WHists.h"
#include "SFrameTools/include/TopFitCalc.h"
#include "include/SubstructureSelectionModules.h"
#include "NtupleWriter/include/JetProps.h"


ClassImp( WPostSelectionCycle );

WPostSelectionCycle::WPostSelectionCycle()
  : AnalysisCycle() {

  // constructor, declare additional variables that should be 
  // obtained from the steering-xml file
  
  // set the integrated luminosity per bin for the lumi-yield control plots
  SetIntLumiPerBin(500.);
  m_sys_var = e_Default;
  m_sys_unc = e_None;

  m_mttgencut = false;
  DeclareProperty( "ApplyMttbarGenCut", m_mttgencut );
  DeclareProperty( "Electron_Or_Muon_Selection", m_Electron_Or_Muon_Selection );

  //default: no btagging cuts applied, other cuts can be defined in config file
  m_Nbtags_min=0;
  m_Nbtags_max=int_infinity();
  DeclareProperty( "Nbtags_min", m_Nbtags_min);
  DeclareProperty( "Nbtags_max", m_Nbtags_max);

  // steering property for data-driven qcd in electron channel
  m_reversed_electron_selection = false;
  DeclareProperty( "ReversedElectronSelection", m_reversed_electron_selection);

  // veto electron trigger when running on the JetHT data stream
  m_veto_electron_trigger = false;
  DeclareProperty( "vetoElectronTrigger", m_veto_electron_trigger);

  // put the selected trigger in OR with HLT_PFJet320_v* (electron channel)
  m_useORTriggerWithPFJet320 = false;
  DeclareProperty( "useORTriggerWithPFJet320", m_useORTriggerWithPFJet320);

  m_dobsf = "None";
  DeclareProperty( "BTaggingScaleFactors", m_dobsf );

  x_btagtype = e_CSVM;

  //
  m_matching = "None"; 
  DeclareProperty( "Matching", m_matching); 

}

WPostSelectionCycle::~WPostSelectionCycle() 
{
  // destructor
}

void WPostSelectionCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void WPostSelectionCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void WPostSelectionCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{
  // declaration of histograms and selections.
  // AnalysisCyle expects Selections and HistCollections to be registered here.
  // Their memory will be released in AnalysisCycle::EndInputData.

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );

 
  // -------------------- set up the selections ---------------------------

  bool doEle=false;
  bool doMu=false;

  if(m_Electron_Or_Muon_Selection=="Electrons" || m_Electron_Or_Muon_Selection=="Electron" || m_Electron_Or_Muon_Selection=="Ele" || m_Electron_Or_Muon_Selection=="ELE") {
    doEle=true;
  } else if(m_Electron_Or_Muon_Selection=="Muon" || m_Electron_Or_Muon_Selection=="Muons" || m_Electron_Or_Muon_Selection=="Mu" || m_Electron_Or_Muon_Selection=="MU") {
    doMu=true;
  } else {
    m_logger << ERROR << "Electron_Or_Muon_Selection is not defined in your xml config file --- should be either `ELE` or `MU`" << SLogger::endmsg;
  }
 
  /////////////////////////////////////////////////////////////////////////
  //reco postselections
  ///////////////////////////////////////////////////////////////////////// 

  Selection* LeptonicHemisphere_Selection = new Selection("LeptonicHemisphereSelection");
  LeptonicHemisphere_Selection->addSelectionModule(new LeptonicHemisphereSelection());
  RegisterSelection(LeptonicHemisphere_Selection);

  Selection* HadronicHemisphere_Selection = new Selection("HadronicHemisphereSelection");
  HadronicHemisphere_Selection->addSelectionModule(new HadronicHemisphereSelection());
  RegisterSelection(HadronicHemisphere_Selection);

  // ---------------- set up the histogram collections --------------------
 
  RegisterHistCollection( new WHists("WLepHists") );
  RegisterHistCollection( new WHists("WFullHists") );
  RegisterHistCollection( new WHists("W_PTCutHists") );
  RegisterHistCollection( new WHists("W_DRHists") );
  RegisterHistCollection( new WHists("WLowMassHists") );


  RegisterHistCollection( new EventHists("Event_PreSel") );
  RegisterHistCollection( new JetHists("Jets_PreSel") );
  RegisterHistCollection( new ElectronHists("Electron_PreSel") );
  RegisterHistCollection( new MuonHists("Muon_PreSel") );
  RegisterHistCollection( new TopJetHists("TopJets_PreSel") );

  RegisterHistCollection( new EventHists("Event_LeptonicSel") );
  RegisterHistCollection( new JetHists("Jets_LeptonicSel") );
  RegisterHistCollection( new ElectronHists("Electron_LeptonicSel") );
  RegisterHistCollection( new MuonHists("Muon_LeptonicSel") );
  RegisterHistCollection( new TopJetHists("TopJets_LeptonicSel") );

  RegisterHistCollection( new EventHists("Event_FullSel") );
  RegisterHistCollection( new JetHists("Jets_FullSel") );
  RegisterHistCollection( new ElectronHists("Electron_FullSel") );
  RegisterHistCollection( new MuonHists("Muon_FullSel") );
  RegisterHistCollection( new TopJetHists("TopJets_FullSel") );
 
  // important: initialise histogram collections after their definition
  InitHistos();


  // Data-MC b-tagging reweighting
  m_bsf = NULL;
  std::transform(m_dobsf.begin(), m_dobsf.end(), m_dobsf.begin(), ::tolower);
  if(m_dobsf != "none") {
    E_SystShift sys_bjets = e_Default;
    E_SystShift sys_ljets = e_Default;
    if (m_dobsf == "default") {
      m_logger << INFO << "Applying btagging scale factor" << SLogger::endmsg;
    } else if (m_dobsf == "up-bjets") {
      m_logger << INFO << "Applying btagging up scale factor for b-jets" << SLogger::endmsg;
      sys_bjets = e_Up;
    } else if (m_dobsf == "down-bjets") {
      m_logger << INFO << "Applying btagging down scale factor for b-jets" << SLogger::endmsg;
      sys_bjets = e_Down;
    } else if (m_dobsf == "up-ljets") {
      m_logger << INFO << "Applying btagging up scale factor for l-jets" << SLogger::endmsg;
      sys_ljets = e_Up;
    } else if (m_dobsf == "down-ljets") {
      m_logger << INFO << "Applying btagging down scale factor for l-jets" << SLogger::endmsg;
      sys_ljets = e_Down;


    }
    else
      m_logger << ERROR << "Unknown BTaggingScaleFactors option, default option is applied --- should be either `Default`, `Up-bjets`, `Down-bjets`, `Up-ljets`, or `Down-ljets`" << SLogger::endmsg;
    if(doEle)
      m_bsf = new BTaggingScaleFactors(x_btagtype, e_Electron, sys_bjets, sys_ljets);
    else if(doMu)
      m_bsf = new BTaggingScaleFactors(x_btagtype, e_Muon, sys_bjets, sys_ljets);
  }

}

void WPostSelectionCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );
}

void WPostSelectionCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer

  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );
}

void WPostSelectionCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError ) 
{
 
  // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results

  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset
  AnalysisCycle::ExecuteEvent( id, weight );

  // get the histogram collections. NOTE: this could be done more performant by making
  // all thse BaseHists* vairables private member variables of SubstructureCycle and
  // setting them in BeginInputData. Then, there is no need here to call GetHistColletion ...


  BaseHists* Hists_WLep = GetHistCollection("WLepHists");
  BaseHists* Hists_WFull = GetHistCollection("WFullHists");
  BaseHists* Hists_W_PTCut = GetHistCollection("W_PTCutHists");
  BaseHists* Hists_W_DRCut = GetHistCollection("W_DRHists");
  BaseHists* Hists_W_LowMass = GetHistCollection("WLowMassHists");

  static Selection* LeptonicHemisphereSelection = GetSelection("LeptonicHemisphereSelection");
  static Selection* HadronicHemisphereSelection = GetSelection("HadronicHemisphereSelection"); 
 
  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();


  //if(!mttbar_gen_selection->passSelection())  throw SError( SError::SkipEvent );

 
  // b tagging scale factor
  if(m_bsf && m_addGenInfo) {
    calc->ProduceWeight(m_bsf->GetWeight());
    calc->ProduceRecWeight(m_bsf->GetWeight());
  }
  



  //matching by hand, should be cahnged!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TTbarHemisphereReconstruction* tt = new TTbarHemisphereReconstruction(bcc);

  tt->ReconstructLeptonicHemisphere();
  tt->ReconstructHadronicHemisphere();

  if(m_matching == "Matched" || m_matching == "MisMatched"){
    TTbarGen *ttgen = new TTbarGen(bcc); 

    double dr1 = 10;
    double dr2 = 10;
    double drqq = 10;

    if( ttgen->DecayChannel() == TTbarGen::e_ehad || ttgen->DecayChannel() == TTbarGen::e_muhad ){
      dr1 = deltaR(tt->WJet().v4(), ttgen->Q1().v4());
      dr2 = deltaR(tt->WJet().v4(), ttgen->Q2().v4());
      drqq = deltaR(ttgen->Q1().v4(), ttgen->Q2().v4());
    }

    if(m_matching == "Matched" && (dr1 > 1.0 || dr2 > 1.0 || drqq > 1.2) ) 
      throw SError( SError::SkipEvent );
    else if( m_matching == "MisMatched" && dr1 < 1.0 && dr2 < 1.0 && drqq < 1.2 ) 
      throw SError( SError::SkipEvent );

    delete ttgen;
  }

  //Fill Histogram collections

  FillControlHists("_PreSel");

  if(!LeptonicHemisphereSelection->passSelection()) throw SError( SError::SkipEvent );

  FillControlHists("_LeptonicSel");
   Hists_WLep->Fill();

  if(!HadronicHemisphereSelection->passSelection()) throw SError( SError::SkipEvent );
  
  FillControlHists("_FullSel");
  Hists_WFull->Fill();

  if(tt->WJet().v4().mass() < 65) Hists_W_LowMass->Fill();

  if(deltaR(tt->WJet().v4(), tt->bJetHad().v4()) < 2.1){
    Hists_W_DRCut->Fill();
    if(tt->WJet().pt() < 260)  Hists_W_PTCut->Fill();
  }

  delete tt;
  return;
  
}

void WPostSelectionCycle::FillControlHists(TString postfix)
{
  // fill some control histograms, need to be defined in BeginInputData

  BaseHists* eventhists = GetHistCollection((std::string)("Event"+postfix));
  BaseHists* jethists = GetHistCollection((std::string)("Jets"+postfix));
  BaseHists* elehists = GetHistCollection((std::string)("Electron"+postfix));
  BaseHists* muonhists = GetHistCollection((std::string)("Muon"+postfix));
  //  BaseHists* tauhists = GetHistCollection((std::string)("Tau"+postfix));
  BaseHists* topjethists = GetHistCollection((std::string)("TopJets"+postfix));

  eventhists->Fill();
  jethists->Fill();
  elehists->Fill();
  muonhists->Fill();
  // tauhists->Fill();
  topjethists->Fill();

}


