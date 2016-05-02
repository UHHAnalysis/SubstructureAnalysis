// $Id: WSelectionCycle.cxx,v 1.10 2012/12/07 14:21:51 peiffer Exp $

#include <iostream>

using namespace std;

#include "include/WSelectionCycle.h"

#include "SFrameAnalysis/include/SelectionModules.h"
//#include "SFrameAnalysis/include/HypothesisHists.h"
#include "SFrameAnalysis/include/TopJetHists.h"
#include "SFrameAnalysis/include/JetHists.h"
#include "SFrameTools/include/TopFitCalc.h"
#include "include/SubstructureSelectionModules.h"
#include "NtupleWriter/include/JetProps.h"


ClassImp( WSelectionCycle );

WSelectionCycle::WSelectionCycle()
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
}

WSelectionCycle::~WSelectionCycle() 
{
  // destructor
}

void WSelectionCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void WSelectionCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void WSelectionCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{
  // declaration of histograms and selections.
  // AnalysisCyle expects Selections and HistCollections to be registered here.
  // Their memory will be released in AnalysisCycle::EndInputData.

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );

 
  // -------------------- set up the selections ---------------------------


  //////////////////////////////////////////////////////////////////////////
  // ZprimePreSelection (applied here again to be able to run over SubstructurPreselection)
  //////////////////////////////////////////////////////////////////////////
  Selection* preselection = new Selection("preselection");

  if(m_Electron_Or_Muon_Selection=="Electrons" || m_Electron_Or_Muon_Selection=="Electron" || m_Electron_Or_Muon_Selection=="Ele" || m_Electron_Or_Muon_Selection=="ELE") {
    preselection->addSelectionModule(new NElectronSelection(1,int_infinity()));//at least one electron
    preselection->addSelectionModule(new NMuonSelection(0,0));//no muons
  } else if(m_Electron_Or_Muon_Selection=="Muon" || m_Electron_Or_Muon_Selection=="Muons" || m_Electron_Or_Muon_Selection=="Mu" || m_Electron_Or_Muon_Selection=="MU") {
    preselection->addSelectionModule(new NElectronSelection(0,0));//no electron
    preselection->addSelectionModule(new NMuonSelection(1,int_infinity()));//at least one muon
  } else {
    m_logger << ERROR << "Electron_Or_Muon_Selection is not defined in your xml config file --- should be either `ELE` or `MU`" << SLogger::endmsg;
  }

  preselection->addSelectionModule(new NJetSelection(2));//at least two jets

  RegisterSelection(preselection);
  ///////////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////////
  //ZprimeSelection
  ///////////////////////////////////////////////////////////////////////////////
 
  static Selection* mttbar_gen_selection = new Selection("Mttbar_Gen_Selection");
  if ( m_mttgencut && ((id.GetVersion() == "TTbar_0to700") || (id.GetVersion() == "TTbar") )  ) {
    m_logger << INFO << "Applying mttbar generator cut from 0 to 700 GeV." << SLogger::endmsg;
    mttbar_gen_selection->addSelectionModule(new MttbarGenCut(0,700));
    mttbar_gen_selection->EnableSelection();
  } else {
    m_logger << INFO << "Disabling mttbar generator cut." << SLogger::endmsg;
    mttbar_gen_selection->DisableSelection();
  }

  bool doEle=false;
  bool doMu=false;

  if(m_reversed_electron_selection)
    m_logger << INFO << "Applying reversed electron selection (data-driven qcd) !!!!" << SLogger::endmsg;

  if(m_Electron_Or_Muon_Selection=="Electrons" || m_Electron_Or_Muon_Selection=="Electron" || m_Electron_Or_Muon_Selection=="Ele" || m_Electron_Or_Muon_Selection=="ELE") {
    doEle=true;
  } else if(m_Electron_Or_Muon_Selection=="Muon" || m_Electron_Or_Muon_Selection=="Muons" || m_Electron_Or_Muon_Selection=="Mu" || m_Electron_Or_Muon_Selection=="MU") {
    doMu=true;
  } else {
    m_logger << ERROR << "Electron_Or_Muon_Selection is not defined in your xml config file --- should be either `ELE` or `MU`" << SLogger::endmsg;
  }
    
  Selection* Ele30trig_selection= new Selection("Ele30trig_selection");
  Ele30trig_selection->addSelectionModule(new TriggerSelection("HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v"));

  Selection* PFJet320trig_selection = new Selection("PFJet320trig_selection");
  PFJet320trig_selection->addSelectionModule(new TriggerSelection("HLT_PFJet320_v"));

  Selection* trig_selection= new Selection("trig_selection");
  trig_selection->addSelectionModule(new TriggerSelection(m_lumi_trigger));

  Selection* first_selection= new Selection("first_selection");
  first_selection->addSelectionModule(new NPrimaryVertexSelection(1)); //at least one good PV
  first_selection->addSelectionModule(new NJetSelection(2,int_infinity(),30,2.4));//at least two jets

  if(doEle) {
    first_selection->addSelectionModule(new NElectronSelection(1,int_infinity()));//at least one electron
    first_selection->addSelectionModule(new NElectronSelection(1,1));//exactly one electron
    first_selection->addSelectionModule(new NMuonSelection(0,0));//no muons
  }
  if(doMu) {
    first_selection->addSelectionModule(new NMuonSelection(1,int_infinity()));//at least one muon
    first_selection->addSelectionModule(new NMuonSelection(1,1));//exactly one muon
    first_selection->addSelectionModule(new NElectronSelection(0,0));//no ided electrons
  }

  first_selection->addSelectionModule(new TwoDCut());

  Selection* second_selection= new Selection("second_selection");

  second_selection->addSelectionModule(new BJetCloseToLeptonSelection(30, 2.1, e_CSVM, TMath::Pi()/2. )); //!!!!Malte


  Selection* trangularcut_selection= new Selection("trangularcut_selection");
  if(doEle) trangularcut_selection->addSelectionModule(new TriangularCut());//triangular cuts

 
  RegisterSelection(mttbar_gen_selection);
  RegisterSelection(Ele30trig_selection);
  RegisterSelection(PFJet320trig_selection);
  RegisterSelection(trig_selection);
  RegisterSelection(first_selection);
  RegisterSelection(second_selection);
  RegisterSelection(trangularcut_selection);

  ///////////////////////////////////////////////////////////////////////////////

 
 
  // ---------------- set up the histogram collections --------------------
 
  // RegisterHistCollection( new ComparisonHists("comp_presel") );
 
  // important: initialise histogram collections after their definition
  // control histogras
  RegisterHistCollection( new EventHists("Event_Presel") );
  RegisterHistCollection( new JetHists("Jets_Presel") );
  RegisterHistCollection( new ElectronHists("Electron_Presel") );
  RegisterHistCollection( new MuonHists("Muon_Presel") );
  RegisterHistCollection( new TauHists("Tau_Presel") );
  RegisterHistCollection( new TopJetHists("TopJets_Presel") );

  RegisterHistCollection( new EventHists("Event_Cleaned") );
  RegisterHistCollection( new JetHists("Jets_Cleaned") );
  RegisterHistCollection( new ElectronHists("Electron_Cleaned") );
  RegisterHistCollection( new MuonHists("Muon_Cleaned") );
  RegisterHistCollection( new TauHists("Tau_Cleaned") );
  RegisterHistCollection( new TopJetHists("TopJets_Cleaned") );

  RegisterHistCollection( new EventHists("Event_Postsel") );
  RegisterHistCollection( new JetHists("Jets_Postsel") );
  RegisterHistCollection( new ElectronHists("Electron_Postsel") );
  RegisterHistCollection( new MuonHists("Muon_Postsel") );
  RegisterHistCollection( new TauHists("Tau_Postsel") );
  RegisterHistCollection( new TopJetHists("TopJets_Postsel") );
  // RegisterHistCollection( new NeutrinoHists("Neutrino_Postsel" , m_chi2discr) );

  InitHistos();

}

void WSelectionCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );
}

void WSelectionCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer

  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );
}

void WSelectionCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError ) 
{
 
  // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results

  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset
  AnalysisCycle::ExecuteEvent( id, weight );

  // get the histogram collections. NOTE: this could be done more performant by making
  // all thse BaseHists* vairables private member variables of SubstructureCycle and
  // setting them in BeginInputData. Then, there is no need here to call GetHistColletion ...

  static Selection* preselection = GetSelection("preselection");
  static Selection* Ele30trig_selection = GetSelection("Ele30trig_selection");
  static Selection* PFJet320trig_selection = GetSelection("PFJet320trig_selection");
  static Selection* trig_selection = GetSelection("trig_selection");
  static Selection* first_selection = GetSelection("first_selection");
  static Selection* second_selection = GetSelection("second_selection");
  static Selection* trangularcut_selection = GetSelection("trangularcut_selection");
 

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

 if(!preselection->passSelection()) throw SError( SError::SkipEvent );
 FillControlHists("_Presel"); 
 

  m_cleaner = new Cleaner();
  m_cleaner->SetJECUncertainty(m_jes_unc);

  // settings for jet correction uncertainties
  if (m_sys_unc==e_JEC){
    if (m_sys_var==e_Up) m_cleaner->ApplyJECVariationUp();
    if (m_sys_var==e_Down) m_cleaner->ApplyJECVariationDown();
  }
  if (m_sys_unc==e_JER){
    if (m_sys_var==e_Up){ m_cleaner->ApplyJERVariationUp();
       m_cleaner->ApplyfatJERVariationUp(); 
    }
    if (m_sys_var==e_Down){ m_cleaner->ApplyJERVariationDown();
       m_cleaner->ApplyfatJERVariationDown(); 
    }
  }

  
  if(bcc->pvs) m_cleaner->PrimaryVertexCleaner(4, 24., 2.);
  if(bcc->electrons) m_cleaner->ElectronCleaner_noIso(35,2.5, m_reversed_electron_selection,true);
  if(bcc->muons) m_cleaner->MuonCleaner_noIso(45,2.1);
  if(bcc->jets) m_cleaner->JetLeptonSubtractor(m_corrector,false);
  if(bcc->topjets) m_cleaner->TopJetLeptonSubtractor(m_corrector,true);


  if(!bcc->isRealData && bcc->jets) m_cleaner->JetEnergyResolutionShifter();
  if(!bcc->isRealData && bcc->topjets) m_cleaner->JetEnergyResolutionShifterFat();


  //apply loose jet cleaning for 2D cut
  if(bcc->jets) m_cleaner->JetCleaner(25,double_infinity(),true);

  // control histograms
  FillControlHists("_Cleaned");

  if(m_veto_electron_trigger && Ele30trig_selection->passSelection()){
    throw SError( SError::SkipEvent );
  }

  bool triggerbit(0);
  if(!m_useORTriggerWithPFJet320) triggerbit = trig_selection->passSelection();
  else triggerbit = trig_selection->passSelection() || PFJet320trig_selection->passSelection();

  if(!triggerbit) throw SError( SError::SkipEvent );


  if(!first_selection->passSelection()) throw SError( SError::SkipEvent );


    
  //apply tighter jet cleaning for further cuts and analysis steps
  if(bcc->jets) m_cleaner->JetCleaner(25,2.5,true);

  //remove all taus from collection for HTlep calculation
  if(bcc->taus) m_cleaner->TauCleaner(double_infinity(),0.0);
 
  if(!second_selection->passSelection()) throw SError( SError::SkipEvent );

  if(!m_reversed_electron_selection) {
    if(!trangularcut_selection->passSelection()) throw SError( SError::SkipEvent );
  } else {
    if(!trangularcut_selection->passInvertedSelection()) throw SError( SError::SkipEvent );
  }


  // control histograms 
  FillControlHists("_Postsel");

  //debugging!!!
  //cout << unshifted_pt << ",  " << shifted_pt << endl;
  //  TTbarHemisphereReconstruction* tt = new TTbarHemisphereReconstruction(bcc);

  //  tt->ReconstructLeptonicHemisphere();

  
 
  WriteOutputTree();


  delete m_cleaner;
 
  return;
  
}

void WSelectionCycle::FillControlHists(TString postfix)
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


