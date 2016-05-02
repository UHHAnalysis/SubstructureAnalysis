#include <iostream>

using namespace std;

#include "include/SubstructurePreSelectionCycle.h"
#include "include/JetMassTools.h"
#include "SFrameAnalysis/include/SelectionModules.h"
#include "include/SubstructureSelectionModules.h"
#include "SFrameTools/JetMETObjects/interface/JetCorrectorParameters.h"


ClassImp( SubstructurePreSelectionCycle );

SubstructurePreSelectionCycle::SubstructurePreSelectionCycle()
  : AnalysisCycle() {

  // constructor, declare additional variables that should be 
  // obtained from the steering-xml file
  
  // set the integrated luminosity per bin for the lumi-yield control plots
  SetIntLumiPerBin(25.);

  DeclareProperty( "Electron_Or_Muon_Selection", m_Electron_Or_Muon_Selection );
 
  // steering property for data-driven qcd in electron channel
  m_reversed_electron_selection = false;
  DeclareProperty( "ReversedElectronSelection", m_reversed_electron_selection);

  //apply tight selection on reco level (pt1>300GeV, pt2>100 GeV) for masspoints and systematics
  m_tight_reco_selection = false;
  DeclareProperty( "ApplyTightRecoSelection", m_tight_reco_selection);
 
}

SubstructurePreSelectionCycle::~SubstructurePreSelectionCycle() 
{
  // destructor
}

void SubstructurePreSelectionCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void SubstructurePreSelectionCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void SubstructurePreSelectionCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{
  // declaration of histograms and selections.
  // AnalysisCyle expects Selections and HistCollections to be registered here.
  // Their memory will be released in AnalysisCycle::EndInputData.

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );

  //  m_tree = GetOutputMetadataTree("PreSelectionTree");
  // m_tree->Branch("rec_presel" , &m_rec_presel , "rec_presel/O");
  // m_tree->Branch("gen_presel" , &m_gen_presel , "gen_presel/O");
 
  
  Book( TH1F("top_pt","p_{t, hadronic top}",80,0,800)); //to have this without any cuts 
  Book( TH1F("genweights","genweights",1,0,1));

  // -------------------- set up the selections ---------------------------

  
  preselection = new Selection("preselection");

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

  if(m_tight_reco_selection){  //default = false
    preselection->addSelectionModule(new NTopJetSelection(1,int_infinity(),150,2.5 ));
    preselection->addSelectionModule(new NTopJetSelection(2,int_infinity(),50,2.5 ));
  }

  RegisterSelection(preselection);
 
  
  genpreselection = new Selection("genpreselection");
  if(m_Electron_Or_Muon_Selection=="Electrons" || m_Electron_Or_Muon_Selection=="Electron" || m_Electron_Or_Muon_Selection=="Ele" || m_Electron_Or_Muon_Selection=="ELE") {
    genpreselection->addSelectionModule(new SubstructureGeneratorPreSelectionElectron(0));
  } else if(m_Electron_Or_Muon_Selection=="Muon" || m_Electron_Or_Muon_Selection=="Muons" || m_Electron_Or_Muon_Selection=="Mu" || m_Electron_Or_Muon_Selection=="MU") {
    genpreselection->addSelectionModule(new SubstructureGeneratorPreSelectionMuon(0));
  }
  genpreselection->addSelectionModule(new NCAGenJetSelection(1,int_infinity(),300,2.5));
  genpreselection->addSelectionModule(new NCAGenJetSelection(2,int_infinity(),100,2.5));
  RegisterSelection(genpreselection);

  genpreselection_hadron = new Selection("genpreselection_hadron");
  genpreselection_hadron->addSelectionModule(new SubstructureGeneratorPreSelectionHadron());
  genpreselection_hadron->addSelectionModule(new NCAGenJetSelection(2,int_infinity(), 300, 2.5));
  RegisterSelection(genpreselection_hadron);
  
  return;
}



void SubstructurePreSelectionCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );
}

void SubstructurePreSelectionCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer

  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );

  return;
}

void SubstructurePreSelectionCycle::ExecuteEvent( const SInputData& id, Double_t weight ) throw( SError ) 
{
  // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results

  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset
  AnalysisCycle::ExecuteEvent( id, weight);

  
  Cleaner cleaner;
  GenCleaner gen_cleaner;
  //static Selection* preselection = GetSelection("preselection");
  //static Selection* genpreselection = GetSelection("genpreselection");
  //static Selection* jetsel = GetSelection("jetsel");


  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();


  Hist("genweights")->Fill(0.5, calc->GetGenWeight());

  // save uncleaned jet collection and MET to be stored in output
  std::vector<Jet> uncleaned_jets;
  for(unsigned int i=0; i<bcc->jets->size(); ++i) {
    uncleaned_jets.push_back(bcc->jets->at(i));
  }
  MET uncleaned_met = *bcc->met;
 
  //clean collections here
  
  if(bcc->muons) cleaner.MuonCleaner_noIso(45,2.1);
  if(bcc->electrons) cleaner.ElectronCleaner_noIso(35,2.5,m_reversed_electron_selection,true);
  if(bcc->jets) cleaner.JetLeptonSubtractor(m_corrector,false);
  if(!bcc->isRealData && bcc->jets) cleaner.JetEnergyResolutionShifter();
  if(bcc->jets) cleaner.JetCleaner(30,2.5,true);
  //if(!bcc->isRealData && bcc->cagenjets) gen_cleaner.CAGenJetCleaner(50,2.5);

  // -------------------------------------------------------------------------
  //get the jet selections

  m_rec_presel = false;
  m_gen_presel = false;
  m_gen_presel_hadron = false;

  if(id.GetVersion().Contains("TT_powheg") || id.GetVersion().Contains("TT_Mtt"))
    m_gen_presel_hadron = genpreselection_hadron->passSelection();

  if( id.GetVersion().Contains("TT")) m_gen_presel = genpreselection->passSelection();

  m_rec_presel = preselection->passSelection();


 if(id.GetVersion().Contains("TT_powheg") || id.GetVersion().Contains("TT_Mtt")){ //hadron selection only for POWHEG ttbar
    if(!genpreselection->passSelection() && !genpreselection_hadron->passSelection() && !preselection->passSelection())throw SError( SError::SkipEvent );
  }else{
    if( id.GetVersion().Contains("TT")){
      if( !genpreselection->passSelection() && !preselection->passSelection()){ throw SError( SError::SkipEvent );}
    }else{
    if(!preselection->passSelection())  throw SError( SError::SkipEvent );
    }
  }
  


  //fill the uncleaned collections back to bcc to store them in output tree
  bcc->met->set_pt (uncleaned_met.pt());
  bcc->met->set_phi (uncleaned_met.phi());
  bcc->jets->clear();
  for(unsigned int i=0; i<uncleaned_jets.size(); ++i) {
    bcc->jets->push_back(uncleaned_jets.at(i));
  }

  
  WriteOutputTree();
   
  return;

}

