// $Id: SubstructureCycle.cxx,v 1.10 2012/12/07 14:21:51 peiffer Exp $

#include <iostream>

using namespace std;

#include "include/SubstructureCycle.h"
#include "include/SubstructureHists.h"
#include "include/GenHists.h"

#include "SFrameAnalysis/include/SelectionModules.h"
#include "SFrameAnalysis/include/HypothesisHists.h"
#include "SFrameAnalysis/include/TopJetHists.h"
#include "SFrameAnalysis/include/JetHists.h"
#include "SFrameTools/include/TopFitCalc.h"
#include "include/SubstructureSelectionModules.h"

ClassImp( SubstructureCycle );

SubstructureCycle::SubstructureCycle()
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

SubstructureCycle::~SubstructureCycle() 
{
  // destructor
}

void SubstructureCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void SubstructureCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void SubstructureCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{
  // declaration of histograms and selections.
  // AnalysisCyle expects Selections and HistCollections to be registered here.
  // Their memory will be released in AnalysisCycle::EndInputData.

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );


  Book( TH1F("eff","efficiancy after the cuts",14,0,14));
  // -------------------- set up the selections ---------------------------


  //////////////////////////////////////////////////////////////////////////
  //for efficiancy check of the ZprimePreSelection
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
  //for efficiancy check of the ZprimeSelection Muon
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

  //Set-Up Selection

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
  first_selection->addSelectionModule(new NJetSelection(2,int_infinity(),50,2.4));//at least two jets

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

  second_selection->addSelectionModule(new NJetSelection(1,int_infinity(),150,2.5)); //leading topjet with pt>150 GeV
  second_selection->addSelectionModule(new NBTagSelection(m_Nbtags_min,m_Nbtags_max)); //b tags from config file
  //  second_selection->addSelectionModule(new HTlepCut(150));
  second_selection->addSelectionModule(new METCut(20));
  // second_selection->addSelectionModule(new TriangularCut());

  Selection* HTLSel_50 = new Selection("HTLSel_50");
  HTLSel_50->addSelectionModule(new HTlepCut(50));
  RegisterSelection(HTLSel_50);

  Selection* HTLSel_100 = new Selection("HTLSel_100");
  HTLSel_100->addSelectionModule(new HTlepCut(100));
  RegisterSelection(HTLSel_100);

  Selection* HTLSel_150 = new Selection("HTLSel_150");
  HTLSel_150->addSelectionModule(new HTlepCut(150));
  RegisterSelection(HTLSel_150);

  Selection* trangularcut_selection= new Selection("trangularcut_selection");
  if(doEle) trangularcut_selection->addSelectionModule(new TriangularCut());//triangular cuts

  Selection* chi2_selection= new Selection("chi2_selection");
  m_chi2discr = new Chi2Discriminator();
  chi2_selection->addSelectionModule(new HypothesisDiscriminatorCut( m_chi2discr, -1*double_infinity(), 10));

  m_bpdiscr = new BestPossibleDiscriminator();
  m_sumdrdiscr = new SumDeltaRDiscriminator();
  m_cmdiscr = new CorrectMatchDiscriminator();

  Selection* matchable_selection = new Selection("matchable_selection");
  matchable_selection->addSelectionModule(new HypothesisDiscriminatorCut( m_cmdiscr, -1*double_infinity(), 999));

  Selection* TopTagSel = new Selection("TopTagSelection");
  TopTagSel->addSelectionModule(new NTopJetSelection(1,int_infinity(),350,2.5));// top jet
  TopTagSel->addSelectionModule(new NTopTagSelection(1,int_infinity()));
  TopTagSel->addSelectionModule(new TopTagOverlapSelection());

  RegisterSelection(mttbar_gen_selection);
  RegisterSelection(Ele30trig_selection);
  RegisterSelection(PFJet320trig_selection);
  RegisterSelection(trig_selection);
  RegisterSelection(first_selection);
  RegisterSelection(second_selection);
  RegisterSelection(trangularcut_selection);
  RegisterSelection(chi2_selection);
  RegisterSelection(matchable_selection);
  RegisterSelection(TopTagSel);
  ///////////////////////////////////////////////////////////////////////////////


  Selection* genpreselection = new Selection("genpreselection");
  genpreselection->addSelectionModule(new SubstructureGeneratorPreSelection(0));
  //genpreselection->addSelectionModule(new NCAGenJetSelection(100,100,2.5));
  genpreselection->addSelectionModule(new NCAGenJetSelection(2,int_infinity(),100,2.5));
  RegisterSelection(genpreselection);

  Selection* LeptonSelection = new Selection("LeptonSelection");
  LeptonSelection->addSelectionModule(new GenLeptonSelection(45,2.5));
  RegisterSelection(LeptonSelection);

 Selection* Jetsel_350 = new Selection("Jetsel_350");
 // Jetsel_350->addSelectionModule(new NCAGenJetSelection(400,150,2.5));
  Jetsel_350->addSelectionModule(new NCAGenJetSelection(1,int_infinity(),400,2.5));
  Jetsel_350->addSelectionModule(new NCAGenJetSelection(2,int_infinity(),100,2.5));
  RegisterSelection(Jetsel_350);

  Selection* Jetselveto_100 = new Selection("Jetselveto_100");
  Jetselveto_100->addSelectionModule(new NCAGenJetveto(150,2,2.5));
  RegisterSelection(Jetselveto_100);


  BSelT = new Selection( "BSelection_tight");
  // addSelectionModule transfers memory release responsibility to the Selection instance.
  BSelT->addSelectionModule(new NBTagSelection(1, int_infinity(), e_CSVT)); //at least one b tag

  BSelT2 = new Selection( "BSelT2");
  BSelT2->addSelectionModule(new NBTagSelection(2, int_infinity(), e_CSVT)); //at least two b tag


  BSel2M = new Selection( "BSel2M");
  BSel2M->addSelectionModule(new NBTagSelection(2, int_infinity(), e_CSVM)); 

  BSel2L = new Selection( "BSel2L");
  BSel2L->addSelectionModule(new NBTagSelection(2, int_infinity(), e_CSVL)); 

  BSelM = new Selection( "BSelection_medium");
  BSelM->addSelectionModule(new NBTagSelection(1, int_infinity(), e_CSVM));

  BSelL = new Selection( "BSelection_loose");
  BSelL->addSelectionModule(new NBTagSelection(1, int_infinity(), e_CSVL));

  NoBSel = new Selection( "NoBSelection");
  NoBSel->addSelectionModule(new NBTagSelection(0,0)); //no b tags

  chi2_selection= new Selection("chi2_selection");
  Chi2Discriminator* m_chi2discr = new Chi2Discriminator();
  chi2_selection->addSelectionModule(new HypothesisDiscriminatorCut( m_chi2discr, -1*double_infinity(), 10));
  //chi2_selection->addSelectionModule(new MttbarGenCut(0,700));
 
  TopSel = new Selection("TopSelection");
  //DO NOT use trigger selection in PROOF mode at the moment
  //TopSel->addSelectionModule(new TriggerSelection("HLT_PFJet320_v"));
  TopSel->addSelectionModule(new NTopJetSelection(1,int_infinity(),400,2.5));
  TopSel->addSelectionModule(new NTopJetSelection(2,int_infinity(),150,2.5));
  // TopSel->addSelectionModule(new NTopTagSelection(1,int_infinity()));

 


  // RegisterSelection transfers memory release responsibility the Selection instance to AnalysisCycle
  // (will be done in EndInputData)
  RegisterSelection(BSelT);
  RegisterSelection(BSelT2);
  RegisterSelection(BSelM);
  RegisterSelection(BSelL);
  RegisterSelection(BSel2M);
  RegisterSelection(BSel2L);
  RegisterSelection(NoBSel);
  RegisterSelection(TopSel);
  RegisterSelection(chi2_selection);

  // ---------------- set up the histogram collections --------------------

  RegisterHistCollection( new GenHists("350_150") );
  RegisterHistCollection( new GenHists("350_150e") );
  // histograms without any cuts
  RegisterHistCollection( new SubstructureHists("NoCutsNoChi2") );
  RegisterHistCollection( new SubstructureHists("NoCuts") );
  RegisterHistCollection( new HypothesisHists("Chi2_NoCuts", m_chi2discr ) );
  RegisterHistCollection( new TopJetHists("Top_NoCuts") );
  RegisterHistCollection( new JetHists("Jet_NoCuts") );
  RegisterHistCollection( new GenHists("350_150NoCuts") );
  RegisterHistCollection( new GenHists("350_150NoCutse") );
  //histograms with and without b tagging
  RegisterHistCollection( new SubstructureHists("BTagT") );
  RegisterHistCollection( new SubstructureHists("BTagT_2L") );
  RegisterHistCollection( new SubstructureHists("BTagM") );
  RegisterHistCollection( new SubstructureHists("BTagL") );
  RegisterHistCollection( new SubstructureHists("BTag2M") );
  RegisterHistCollection( new SubstructureHists("BTag2L") );
  RegisterHistCollection( new SubstructureHists("BTagTM") );
  RegisterHistCollection( new SubstructureHists("BTagTL") );
  RegisterHistCollection( new SubstructureHists("BTagML") );
  RegisterHistCollection( new GenHists("350_150BTagT") );
  RegisterHistCollection( new GenHists("350_150BTagM") );
  RegisterHistCollection( new GenHists("350_150BTagL") );

  RegisterHistCollection( new GenHists("350_150BTagML") );
  RegisterHistCollection( new GenHists("350_150BTag2M") );
  RegisterHistCollection( new GenHists("350_150BTag2L") );
  RegisterHistCollection( new GenHists("350_150BTagTM") );
  RegisterHistCollection( new GenHists("350_150BTagTL") );

  RegisterHistCollection( new TopJetHists("Top_BTagT") );
  RegisterHistCollection( new TopJetHists("Top_BTagM") );
  RegisterHistCollection( new TopJetHists("Top_BTagL") );
  RegisterHistCollection( new SubstructureHists("NoBTag") );
  RegisterHistCollection( new HypothesisHists("Chi2_BTag", m_chi2discr ) );
  RegisterHistCollection( new HypothesisHists("Chi2_NoBTag", m_chi2discr ) );

  // histograms after the top selection
  RegisterHistCollection( new SubstructureHists("TopSel") );
  RegisterHistCollection( new SubstructureHists("TopSelB") );
  RegisterHistCollection( new SubstructureHists("TopSelB_HTL100") );
  RegisterHistCollection( new SubstructureHists("TopSelB_noHTL") );
  RegisterHistCollection( new SubstructureHists("TopSelB_HTL50") );
  RegisterHistCollection( new SubstructureHists("TopSel2B") );
  RegisterHistCollection( new SubstructureHists("TopSelB_pt") );
  RegisterHistCollection( new SubstructureHists("TopSelB_m") );
  RegisterHistCollection( new GenHists("350_150TopSel") );
  RegisterHistCollection( new GenHists("350_150TopSelB") );
  RegisterHistCollection( new GenHists("350_150TopSelB_noHTL") );
  RegisterHistCollection( new GenHists("350_150TopSelB_HTL50") );
  RegisterHistCollection( new GenHists("350_150TopSelB_HTL100") );
  RegisterHistCollection( new GenHists("350_150TopSel2B") );
  RegisterHistCollection( new GenHists("350_150TopSelBe") );
  RegisterHistCollection( new TopJetHists("Top_TopSel") );
  RegisterHistCollection( new JetHists("Jet_TopSel") );
  RegisterHistCollection( new HypothesisHists("Chi2_TopSel", m_chi2discr ) );

  // important: initialise histogram collections after their definition
  InitHistos();

}

void SubstructureCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );
}

void SubstructureCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer

  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );
}

void SubstructureCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError ) 
{
  // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results

  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset
  AnalysisCycle::ExecuteEvent( id, weight );

  // get the histogram collections. NOTE: this could be done more performant by making
  // all thse BaseHists* vairables private member variables of SubstructureCycle and
  // setting them in BeginInputData. Then, there is no need here to call GetHistColletion ...
  // BaseHists* HistsNoCutsNoChi2 = GetHistCollection("NoCutsNoChi2");
  BaseHists* Gen_Hists = GetHistCollection("350_150");
  BaseHists* Gen_Hists_e = GetHistCollection("350_150e");
  BaseHists* Gen_HistsNoCuts = GetHistCollection("350_150NoCuts");
  BaseHists* Gen_HistsNoCuts_e = GetHistCollection("350_150NoCutse");
  BaseHists* HistsNoCuts = GetHistCollection("NoCuts");
  BaseHists* HistsBTagT = GetHistCollection("BTagT");
  BaseHists* HistsBTagT_2L = GetHistCollection("BTagT_2L");
  BaseHists* HistsBTagM = GetHistCollection("BTagM");
  BaseHists* HistsBTagL = GetHistCollection("BTagL");

  BaseHists* HistsBTag2M = GetHistCollection("BTag2M");
  BaseHists* HistsBTag2L = GetHistCollection("BTag2L");
  BaseHists* HistsBTagTM = GetHistCollection("BTagTM");
  BaseHists* HistsBTagTL = GetHistCollection("BTagTL");
  BaseHists* HistsBTagML = GetHistCollection("BTagML");

  BaseHists* Gen_HistsBTagT = GetHistCollection("350_150BTagT");
  BaseHists* Gen_HistsBTagM = GetHistCollection("350_150BTagM");
  BaseHists* Gen_HistsBTagL = GetHistCollection("350_150BTagL");

  BaseHists* Gen_HistsBTagML = GetHistCollection("350_150BTagML");
  BaseHists* Gen_HistsBTag2M = GetHistCollection("350_150BTag2M");
  BaseHists* Gen_HistsBTag2L = GetHistCollection("350_150BTag2L");
  BaseHists* Gen_HistsBTagTM = GetHistCollection("350_150BTagTM");
  BaseHists* Gen_HistsBTagTL = GetHistCollection("350_150BTagTL");

  BaseHists* HistsNoBTag = GetHistCollection("NoBTag");
  BaseHists* HistsTopSel = GetHistCollection("TopSel");
  BaseHists* HistsTopSelB = GetHistCollection("TopSelB");
  BaseHists* HistsTopSelB_noHTL = GetHistCollection("TopSelB_noHTL");
  BaseHists* HistsTopSelB_HTL50 = GetHistCollection("TopSelB_HTL50");
  BaseHists* HistsTopSelB_HTL100 = GetHistCollection("TopSelB_HTL100");
  BaseHists* HistsTopSelB_pt = GetHistCollection("TopSelB_pt");
  BaseHists* HistsTopSelB_m = GetHistCollection("TopSelB_m");
  BaseHists* HistsTopSel2B = GetHistCollection("TopSel2B");
  BaseHists* Gen_HistsTopSel = GetHistCollection("350_150TopSel");
  BaseHists* Gen_HistsTopSelB = GetHistCollection("350_150TopSelB");
  BaseHists* Gen_HistsTopSelB_noHTL = GetHistCollection("350_150TopSelB_noHTL");
  BaseHists* Gen_HistsTopSelB_HTL50 = GetHistCollection("350_150TopSelB_HTL50");
  BaseHists* Gen_HistsTopSelB_HTL100 = GetHistCollection("350_150TopSelB_HTL100");
  BaseHists* Gen_HistsTopSel2B = GetHistCollection("350_150TopSel2B");
  BaseHists* Gen_HistsTopSelB_e = GetHistCollection("350_150TopSelBe");

  BaseHists* Chi2_HistsNoCuts = GetHistCollection("Chi2_NoCuts"); 
  BaseHists* Top_HistsNoCuts = GetHistCollection("Top_NoCuts");
  BaseHists* Top_HistsBTagT = GetHistCollection("Top_BTagT");
  BaseHists* Top_HistsBTagM = GetHistCollection("Top_BTagM");
  BaseHists* Top_HistsBTagL = GetHistCollection("Top_BTagL");
  BaseHists* Jet_HistsNoCuts = GetHistCollection("Jet_NoCuts");
  BaseHists* Top_HistsTopSel = GetHistCollection("Top_TopSel");
  BaseHists* Jet_HistsTopSel = GetHistCollection("Jet_TopSel");
  BaseHists* Chi2_HistsBTag = GetHistCollection("Chi2_BTag");
  BaseHists* Chi2_HistsNoBTag = GetHistCollection("Chi2_NoBTag");
  BaseHists* Chi2_HistsTopSel = GetHistCollection("Chi2_TopSel");


  static Selection* Ele30trig_selection = GetSelection("Ele30trig_selection");
  static Selection* trangularcut_selection = GetSelection("trangularcut_selection");
  static Selection* preselection = GetSelection("preselection");
  static Selection* genpreselection = GetSelection("genpreselection");
  static Selection* trig_selection = GetSelection("trig_selection");
  static Selection* first_selection = GetSelection("first_selection");
  static Selection* second_selection = GetSelection("second_selection");
  static Selection* PFJet320trig_selection = GetSelection("PFJet320trig_selection");
  static Selection* mttbar_gen_selection = GetSelection("Mttbar_Gen_Selection");
  static Selection* chi2_selection = GetSelection("chi2_selection");
  static Selection* BSelT = GetSelection("BSelection_tight");
  static Selection* BSelT2 = GetSelection("BSelT2");
  static Selection* BSel2M = GetSelection("BSel2M");
  static Selection* BSel2L = GetSelection("BSel2L");
  static Selection* HTLSel_50 = GetSelection("HTLSel_50");
  static Selection* HTLSel_100 = GetSelection("HTLSel_100");
  static Selection* HTLSel_150 = GetSelection("HTLSel_150");


  static Selection* LeptonSelection = GetSelection("LeptonSelection");
  static Selection* Jetselveto_100 = GetSelection("Jetselveto_100");
  static Selection* Jetsel_350 = GetSelection("Jetsel_350");


  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();


  bool ele = false;
  if ( bcc->genparticles){
    std::vector<GenParticle>* genparts = bcc->genparticles;
    int Ngenparts = genparts->size();
 
    for(int l=0; l<Ngenparts ; l++){
      GenParticle genpart = genparts->at(l);
      if(genpart.status()<3){break;}
      if( fabs(genpart.pdgId()) == 6){

	const  GenParticle* daughter = genpart.daughter(genparts, 1)->daughter(genparts, 1);
	if(daughter){
	  if( fabs(daughter->pdgId())>10 && fabs(daughter->pdgId())<13) ele=true;
	}
      }
    }

  }




  bool gensel=false;
  if(genpreselection->passSelection() && LeptonSelection->passSelection()  &&  Jetsel_350->passSelection()  && Jetselveto_100->passSelection()){
    Gen_Hists->Fill();
    gensel=true;
    if(ele){Gen_Hists_e->Fill();}
  }

  ////////////////////////////////////////////////////////////////////////////
  //selections
  ///////////////////////////////////////////////////////////////////////////

  if(gensel && !ele) Hist("eff")->Fill(0.5,weight);

  if( !preselection->passSelection()) throw SError( SError::SkipEvent );
  if(gensel &&  !ele)	Hist("eff")->Fill(1.5,weight);

  if(!mttbar_gen_selection->passSelection())  throw SError( SError::SkipEvent );
  if(gensel && !ele )	  Hist("eff")->Fill(2.5,weight);

  if(m_veto_electron_trigger && Ele30trig_selection->passSelection()){
    throw SError( SError::SkipEvent );
  }
  if(gensel &&  !ele)	    Hist("eff")->Fill(3.5,weight);

  bool triggerbit(0);
  if(!m_useORTriggerWithPFJet320) triggerbit = trig_selection->passSelection();
  else triggerbit = trig_selection->passSelection() || PFJet320trig_selection->passSelection();

  if(!triggerbit) throw SError( SError::SkipEvent );
  if(gensel &&  !ele)     Hist("eff")->Fill(4.5,weight);

  if(!first_selection->passSelection()) throw SError( SError::SkipEvent );
  if(gensel &&  !ele)	  Hist("eff")->Fill(5.5,weight);

  if(!second_selection->passSelection()) throw SError( SError::SkipEvent );
  if(gensel &&  !ele)	  Hist("eff")->Fill(6.5,weight);

  if(!m_reversed_electron_selection) {
    if(!trangularcut_selection->passSelection())  throw SError( SError::SkipEvent );
  } else {
    if(!trangularcut_selection->passInvertedSelection())  throw SError( SError::SkipEvent );
  }

  if(gensel &&  !ele)Hist("eff")->Fill(7.5,weight);
  ////////////////////////////////////////////////////////////////////////////////////

  Chi2Discriminator* tagchi2discr;
  TopFitCalc* topfit = TopFitCalc::Instance();
   
   
  bcc->recoHyps->clear();
  topfit->Reset();
  topfit->CalculateSelection();
  // std::cout << bcc->recoHyps->size() << std::endl;
  if(bcc->recoHyps->size()<1) {
    topfit->FillHighMassTTbarHypotheses();}
  tagchi2discr = new Chi2Discriminator();
  tagchi2discr->FillDiscriminatorValues();

  ReconstructionHypothesis *discr = tagchi2discr->GetBestHypothesis();
  double discr_cut=discr->discriminator("Chi2_tlep");


  if(discr_cut>25)  throw SError( SError::SkipEvent );
 if(gensel)Hist("eff")->Fill(8.5,weight);
  // if(!chi2_selection->passSelection())  throw SError( SError::SkipEvent );

  LorentzVector top_lep = discr->toplep_v4();
  // for (unsigned int i =0; i<bcc->topjets->size(); ++i){
  TopJet topjet;
  if(bcc->topjets->size()>0)   topjet = bcc->topjets->at(0);
  bool jet_distance=true;
  bool selection_thad=true;
  if(sqrt(pow(topjet.phi()-top_lep.phi(),2)+pow(topjet.eta()-top_lep.eta(),2))<2.7 || sqrt(pow(topjet.phi()-top_lep.phi(),2)+pow(topjet.eta()-top_lep.eta(),2))>3.5) selection_thad = false;
  if(selection_thad/*sqrt(pow(topjet.phi()-top_lep.phi(),2)+pow(topjet.eta()-top_lep.eta(),2))<2.7 || sqrt(pow(topjet.phi()-top_lep.phi(),2)+pow(topjet.eta()-top_lep.eta(),2))>3.5*/) {
    for(unsigned t=0;t<bcc->jets->size();++t){
      Jet jet = bcc->jets->at(t);
      if((sqrt(pow(topjet.phi()-jet.phi(),2)+pow(topjet.eta()-jet.eta(),2))>0.8 &&sqrt(pow(topjet.phi()-jet.phi(),2)+pow(topjet.eta()-jet.eta(),2))<1.8) ) jet_distance=false;
      if((sqrt(pow(top_lep.phi()-jet.phi(),2)+pow(top_lep.eta()-jet.eta(),2))>1 &&sqrt(pow(top_lep.phi()-jet.phi(),2)+pow(top_lep.eta()-jet.eta(),2))<2.2)) jet_distance=false;
      if((sqrt(pow(top_lep.phi()-jet.phi(),2)+pow(top_lep.eta()-jet.eta(),2))>4)) jet_distance=false;
    }
  }

  if(!selection_thad && !jet_distance)  throw SError( SError::SkipEvent );
 if(gensel && !ele )Hist("eff")->Fill(9.5,weight);

  bool jetmass2 = true;
  bool ptratio = true;
  if(calc->GetCAJets()->size()>1){
    if(calc->GetCAJets()->at(1).v4().M()>150)  jetmass2=false;
    if(calc->GetCAJets()->at(1).pt()>0.75*topjet.pt()) ptratio=false;
      
    //     std::cout << "mass: " <<calc->GetCAJets()->at(1).v4().M()  << "   pt: " << calc->GetCAJets()->at(1).pt();
  }
  //   std::cout<<"  had: " << had << std::endl;

  if(calc->GetJets()->size()>=12){
    std::cout << "run: " << calc->GetRunNum() << "   lb: " << calc->GetLumiBlock() << "  event: " << calc->GetEventNum() << "   N(jets): " << calc->GetJets()->size() << std::endl;
  }

  // start the analysis
  //(mttbar_gen_selection->passSelection() && first_selection->passSelection() && triggerbit && second_selection->passSelection() && preselection->passSelection()){


  // HistsNoCutsNoChi2->Fill();

  // if(!chi2_selection->passSelection())  throw SError( SError::SkipEvent );


  if(HTLSel_150->passSelection()){
    HistsNoCuts->Fill();
    if(gensel){  Gen_HistsNoCuts->Fill();
      if(ele){Gen_HistsNoCuts_e->Fill();}
    }
    Top_HistsNoCuts->Fill();
    Jet_HistsNoCuts->Fill();
    Chi2_HistsNoCuts->Fill();

    if(BSelT->passSelection()){
      HistsBTagT->Fill();
      if(gensel){ Gen_HistsBTagT->Fill();}
      Top_HistsBTagT->Fill();
      Chi2_HistsBTag->Fill();
    }
    if(BSelM->passSelection()){
      Top_HistsBTagM->Fill();
      HistsBTagM->Fill();
      if(gensel){ Gen_HistsBTagM->Fill();}

    }
    if(BSelL->passSelection()){
      HistsBTagL->Fill();
      Top_HistsBTagL->Fill();
      if(gensel){  Gen_HistsBTagL->Fill();}

    }

    if(BSel2M->passSelection()){
      HistsBTag2M->Fill();
      if(gensel){ Gen_HistsBTag2M->Fill();}
      if(BSelT->passSelection()){
        HistsBTagTM->Fill();
	if(gensel){ Gen_HistsBTagTM->Fill();}
      }
    }

    if(BSel2L->passSelection()){
      HistsBTag2L->Fill();
      if(gensel){ Gen_HistsBTag2L->Fill();}
      if(BSelT->passSelection()){
        HistsBTagTL->Fill();
	if(gensel){ Gen_HistsBTagTL->Fill();}
      }
      if(BSelM->passSelection()){
        HistsBTagML->Fill();
	if(gensel){ Gen_HistsBTagML->Fill();}
      }
    }



    if(NoBSel->passSelection()){
      HistsNoBTag->Fill();  
      Chi2_HistsNoBTag->Fill();
    } 
  }
  
  if(!TopSel->passSelection())  throw SError( SError::SkipEvent );

  if(BSelT->passSelection() && ptratio){
    HistsTopSelB_noHTL->Fill();
    if(gensel){   Gen_HistsTopSelB_noHTL->Fill();}
  }

  if(HTLSel_50->passSelection()){
    if(BSelT->passSelection() && ptratio){
      HistsTopSelB_HTL50->Fill();
      if(gensel){   Gen_HistsTopSelB_HTL50->Fill();}
    }
  }


  if(HTLSel_100->passSelection()){
    if(BSelT->passSelection() && ptratio){
      HistsTopSelB_HTL100->Fill();
      if(gensel){   Gen_HistsTopSelB_HTL100->Fill();}
    }
  }

  if(HTLSel_150->passSelection()){
    if(gensel && !ele )Hist("eff")->Fill(10.5,weight);

    HistsTopSel->Fill();
    if(BSelT->passSelection()){
      HistsTopSelB->Fill();
      if(gensel  && !ele )Hist("eff")->Fill(11.5,weight);
      //    std::cout<<"FILL  ";
      if(ptratio){ 
	HistsTopSelB_pt->Fill();
	if(gensel  && !ele )Hist("eff")->Fill(12.5,weight);
	if (jetmass2){
	  if(gensel  && !ele )Hist("eff")->Fill(13.5,weight);
	  HistsTopSelB_m->Fill();
	  if(gensel){   Gen_HistsTopSelB->Fill();
	    if(!ele){Gen_HistsTopSelB_e->Fill();}
	  }
	}
      }
    }
    if(BSelT2->passSelection()){
      HistsTopSel2B->Fill();
      if(gensel ){   Gen_HistsTopSel2B->Fill();}
    }
    if(gensel){   Gen_HistsTopSel->Fill();}
    Top_HistsTopSel->Fill();
    Jet_HistsTopSel->Fill();
    Chi2_HistsTopSel->Fill();
  }
  return;
  
}


