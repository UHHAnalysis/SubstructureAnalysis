// $Id: SubstructurePostSelectionCycle.cxx,v 1.10 2012/12/07 14:21:51 peiffer Exp $

#include <iostream>

using namespace std;
#include "include/SubstructurePostSelectionCycle.h"


ClassImp( SubstructurePostSelectionCycle );


SubstructurePostSelectionCycle::SubstructurePostSelectionCycle()
  : AnalysisCycle() {

  // constructor, declare additional variables that should be 
  // obtained from the steering-xml file
  
  // set the integrated luminosity per bin for the lumi-yield control plots
  SetIntLumiPerBin(500.);
  m_sys_var = e_Default;
  m_sys_unc = e_None;

  m_mttgencut = false;
  m_dobsf = "None";
  DeclareProperty( "ApplyMttbarGenCut", m_mttgencut );
  DeclareProperty( "Electron_Or_Muon_Selection", m_Electron_Or_Muon_Selection );
  DeclareProperty( "BTaggingScaleFactors", m_dobsf );

  x_btagtype = e_CSVT;

  //default: no btagging cuts applied, other cuts can be defined in config file
  m_Nbtags_min=0;
  m_Nbtags_max=int_infinity();
  DeclareProperty( "Nbtags_min", m_Nbtags_min);
  DeclareProperty( "Nbtags_max", m_Nbtags_max);

  // apply SF for the "Ele30 OR PFJet320" trigger (electron channel)
  m_applyEleORJetTriggerSF = false;
  DeclareProperty( "applyEleORJetTriggerSF", m_applyEleORJetTriggerSF);

  //set grooming hists (should be renamed)
  m_monitoring = false;
  DeclareProperty( "DoMonitoring", m_monitoring);

  //grooming
  m_groom_recojets = false;
  DeclareProperty( "GroomRecoJets", m_groom_recojets);

  m_groom_genjets = false;
  DeclareProperty( "GroomGenJets", m_groom_genjets);

  //  m_groom_genjets = false; //not working at the moment

  //matching 
  m_matching_gen = "None";
  DeclareProperty("GenMatching", m_matching_gen);

  m_matching_reco = "None";
  DeclareProperty("RecoMatching", m_matching_reco);
}

SubstructurePostSelectionCycle::~SubstructurePostSelectionCycle() 
{
  // destructor
}

void SubstructurePostSelectionCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void SubstructurePostSelectionCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void SubstructurePostSelectionCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{
  // declaration of histograms and selections.
  // AnalysisCyle expects Selections and HistCollections to be registered here.
  // Their memory will be released in AnalysisCycle::EndInputData.

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );
 
  
  //decare a metadata output tree that is later used for the unfolding procedure
  m_tree = GetOutputMetadataTree("UnfoldingTree");
  SetOutputTreeBranches();

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

  //========================================================================
  //the mttbar_gen_selection is defined again to provide the possibility to 
  //apply it in the post selection only
  //========================================================================
 
  mttbar_gen_selection = new Selection("Mttbar_Gen_Selection");
  if ( m_mttgencut && ((id.GetVersion() == "TTbar_0to700") || (id.GetVersion() == "TTbar") )  ) {
    m_logger << INFO << "Applying mttbar generator cut from 0 to 700 GeV." << SLogger::endmsg;
    mttbar_gen_selection->addSelectionModule(new MttbarGenCut(0,700));
    mttbar_gen_selection->EnableSelection();
  } else {
    m_logger << INFO << "Disabling mttbar generator cut." << SLogger::endmsg;
    mttbar_gen_selection->DisableSelection();
  }
  RegisterSelection(mttbar_gen_selection);

  //========================================================================
  //The HTlep cut was excluded from the previous selection steps in order to   
  //be able to perform the JEC variation in the post selection
  //========================================================================
    
  HTlep_selection = new Selection("HTlep_selection");
  HTlep_selection->addSelectionModule(new HTlepCut(150));
 
  RegisterSelection(HTlep_selection);


  //===============================
  //post selections on gen level
  //===============================
 

  GenJetSelection_500_100 = new Selection("GenJetSelection_500_100");
  GenJetSelection_500_100->addSelectionModule(new NCAGenJetSelection(1,int_infinity(),500,2.5));
  GenJetSelection_500_100->addSelectionModule(new NCAGenJetSelection(2,int_infinity(),100,2.5));
  RegisterSelection(GenJetSelection_500_100);

  GenJetSelection_400_100 = new Selection("GenJetSelection_400_100");
  GenJetSelection_400_100->addSelectionModule(new NCAGenJetSelection(1,int_infinity(),400,2.5));
  GenJetSelection_400_100->addSelectionModule(new NCAGenJetSelection(2,int_infinity(),100,2.5));
  RegisterSelection(GenJetSelection_400_100);

  GenJetSelection_300_100 = new Selection("GenJetSelection_300_100");
  GenJetSelection_300_100->addSelectionModule(new NCAGenJetSelection(1,int_infinity(),300,2.5));
  GenJetSelection_300_100->addSelectionModule(new NCAGenJetSelection(2,int_infinity(),100,2.5));
  RegisterSelection(GenJetSelection_300_100);

  GenJetSelection_400_150 = new Selection("GenJetSelection_400_150");
  GenJetSelection_400_150->addSelectionModule(new NCAGenJetSelection(1,int_infinity(),400,2.5));
  GenJetSelection_400_150->addSelectionModule(new NCAGenJetSelection(2,int_infinity(),150,2.5));
  RegisterSelection(GenJetSelection_400_150);

  GenJetSelection_300_150 = new Selection("GenJetSelection_300_150");
  GenJetSelection_300_150->addSelectionModule(new NCAGenJetSelection(1,int_infinity(),300,2.5));
  GenJetSelection_300_150->addSelectionModule(new NCAGenJetSelection(2,int_infinity(),150,2.5));
  RegisterSelection(GenJetSelection_300_150);

  // for the hadron selection
  GenJetSelection_400_300 = new Selection("GenJetSelection_400_300");
  GenJetSelection_400_300->addSelectionModule(new NCAGenJetSelection(2,int_infinity(),300,2.5));
  GenJetSelection_400_300->addSelectionModule(new NCAGenJetSelection(1,int_infinity(),400,2.5));
  RegisterSelection(GenJetSelection_400_300);

  GenJetSelection_400_400 = new Selection("GenJetSelection_400_400");
  GenJetSelection_400_400->addSelectionModule(new NCAGenJetSelection(2,int_infinity(),400,2.5));
  GenJetSelection_400_400->addSelectionModule(new NCAGenJetSelection(1,int_infinity(),400,2.5));
  RegisterSelection(GenJetSelection_400_400);



  GenJetVeto_150 = new Selection("GenJetVeto_150");
  GenJetVeto_150->addSelectionModule(new NCAGenJetSelection(0,2,150,2.5));
  //  GenJetVeto_150->addSelectionModule(new NCAGenJetveto(150,2,2.5));
  RegisterSelection(GenJetVeto_150);

  GenJetVeto_100 = new Selection("GenJetVeto_100");
  GenJetVeto_100->addSelectionModule(new NCAGenJetveto(100,2,int_infinity()));
  RegisterSelection(GenJetVeto_100);

  GenJetVeto_30 = new Selection("GenJetVeto_30");
  GenJetVeto_30->addSelectionModule(new NCAGenJetveto(30,2,int_infinity()));
  RegisterSelection(GenJetVeto_30);




  BoostedGenSelection = new Selection("BoostedGenSelection");
  BoostedGenSelection->addSelectionModule(new DR2ndGenCAJetLeptonSelection(1.2));
  RegisterSelection(BoostedGenSelection);




  GenMassSelection = new Selection("GenMassSelection");
  GenMassSelection->addSelectionModule(new LeadingGenCAJetHigherMassSelection("Default"));
  RegisterSelection(GenMassSelection);

  GenMassSelectionLepton = new Selection("GenMassSelectionLepton");
  GenMassSelectionLepton->addSelectionModule(new LeadingGenCAJetHigherMassSelection("Lepton"));
  RegisterSelection(GenMassSelectionLepton);

  GenMassSelectionLeptonAndNeutrino = new Selection("GenMassSelectionLeptonAndNeutrino");
  GenMassSelectionLeptonAndNeutrino->addSelectionModule(new LeadingGenCAJetHigherMassSelection("LeptonAndNeutrino"));
  RegisterSelection(GenMassSelectionLeptonAndNeutrino);

  //matching 
  GenMatchedSelection = new Selection("GenMatchedSelection"); 
  GenMatchedSelection->addSelectionModule(new GenCAJetMatchedSelection(1.2, "LeptonPlusJets"));
  RegisterSelection(GenMatchedSelection);

  GenMatchedSelectionHadron = new Selection("GenMatchedSelectionHadron"); 
  GenMatchedSelectionHadron->addSelectionModule(new GenCAJetMatchedSelection(1.2, "AllHadronic"));
  RegisterSelection(GenMatchedSelectionHadron);


  //========================================
  //post selections on reconstruction level 
  //========================================

  BSelT = new Selection( "BSelection_tight");
  BSelT->addSelectionModule(new NBTagSelection(1, int_infinity(), e_CSVT)); //at least one b tag
  RegisterSelection(BSelT);

  RecoLeptonSelection = new Selection("RecoLeptonSelection");
  if(doEle) RecoLeptonSelection->addSelectionModule(new NElectronSelection(1,int_infinity(), 45 , 2.5));
  if(doMu) RecoLeptonSelection->addSelectionModule(new NMuonSelection(1,int_infinity(), 45 , 2.5));
  RegisterSelection(RecoLeptonSelection);

  TopJetSelection_400_150 = new Selection("TopJetSelection_400_150");
  TopJetSelection_400_150->addSelectionModule(new NTopJetSelection(1,int_infinity(),400,2.5));
  TopJetSelection_400_150->addSelectionModule(new NTopJetSelection(2,int_infinity(),150,2.5));
  RegisterSelection(TopJetSelection_400_150);

  TopJetSelection_500_150 = new Selection("TopJetSelection_500_150");
  TopJetSelection_500_150->addSelectionModule(new NTopJetSelection(1,int_infinity(),500,2.5));
  TopJetSelection_500_150->addSelectionModule(new NTopJetSelection(2,int_infinity(),150,2.5));
  RegisterSelection(TopJetSelection_500_150);

  TopJetSelection_400_100 = new Selection("TopJetSelection_400_100");
  TopJetSelection_400_100->addSelectionModule(new NTopJetSelection(1,int_infinity(),400,2.5));
  TopJetSelection_400_100->addSelectionModule(new NTopJetSelection(2,int_infinity(),100,2.5));
  RegisterSelection(TopJetSelection_400_100);

  TopJetSelection_300_100 = new Selection("TopJetSelection_300_100");
  TopJetSelection_300_100->addSelectionModule(new NTopJetSelection(1,int_infinity(),300,2.5));
  TopJetSelection_300_100->addSelectionModule(new NTopJetSelection(2,int_infinity(),100,2.5));
  RegisterSelection(TopJetSelection_300_100);

  /* WSelection = new Selection("WSelection");
     WSelection->addSelectionModule(new NTopJetSelection(3,int_infinity(),150,2.5));
     RegisterSelection(WSelection);
  */
 

  TopJetSelection_300_150 = new Selection("TopJetSelection_300_150");
  TopJetSelection_300_150->addSelectionModule(new NTopJetSelection(1,int_infinity(),300,2.5));
  TopJetSelection_300_150->addSelectionModule(new NTopJetSelection(2,int_infinity(),150,2.5));
  RegisterSelection(TopJetSelection_300_150);

  RecoJetVeto = new Selection("RecoJetVeto");
  RecoJetVeto->addSelectionModule(new NTopJetSelection(0,2,150,2.5));
  //RecoJetVeto->addSelectionModule(new NTopJetSelection(2,2,150,2.5));
  RegisterSelection(RecoJetVeto);


  BoostedRecoSelection = new Selection("BoostedRecoSelection");
  BoostedRecoSelection->addSelectionModule(new DR2ndTopJetLeptonSelection(1.2));
  RegisterSelection(BoostedRecoSelection);


  //mass selections 
  RecoMassSelection = new Selection("RecoMassSelection");
  RecoMassSelection->addSelectionModule(new LeadingTopJetHigherMassSelection("Default"));
  RegisterSelection(RecoMassSelection);

  RecoMassSelectionLepton = new Selection("RecoMassSelectionLepton");
  RecoMassSelectionLepton->addSelectionModule(new LeadingTopJetHigherMassSelection("Lepton"));
  RegisterSelection(RecoMassSelectionLepton);

  RecoMatchedSelection = new Selection("RecoMatchedSelection"); 
  RecoMatchedSelection->addSelectionModule(new TopJetMatchedSelection(1.2));
  RegisterSelection(RecoMatchedSelection);
 

  // ---------------- set up the histogram collections --------------------
  //hadronic selections
  RegisterHistCollection( Hists_hadronic_gen = new ComparisonHists("hadronic_gen") );
  RegisterHistCollection( Hists_hadronic_gen_mass = new ComparisonHists("hadronic_gen_mass") );
 
  //pt thresholds
  RegisterHistCollection( Hists_gen_presel = new ComparisonHists("gen_presel") );
  RegisterHistCollection( Hists_gen_noVeto_no2nd = new ComparisonHists("gen_noVeto_no2nd") );
  RegisterHistCollection( Hists_gen_noVeto_no2nd_300 = new ComparisonHists("gen_noVeto_no2nd_300") );
  RegisterHistCollection( Hists_gen_noVeto_no2nd_500 = new ComparisonHists("gen_noVeto_no2nd_500") );

  //genlecvel selections
  RegisterHistCollection( Hists_gen_noVeto = new ComparisonHists("gen_noVeto") );
  RegisterHistCollection( Hists_gen = new ComparisonHists("gen") );
  RegisterHistCollection( Hists_gen_dr = new ComparisonHists("gen_dr") );
  RegisterHistCollection( Hists_gen_mass = new ComparisonHists("gen_mass") );
  RegisterHistCollection( Hists_gen_mass2 = new ComparisonHists("gen_mass2") );
  RegisterHistCollection( Hists_gen_mass3 = new ComparisonHists("gen_mass3") );

  RegisterHistCollection( Hists_truth = new TruthHists("truth") );
  RegisterHistCollection( Hists_gen_v50 = new ComparisonHists("gen_v50") );
  RegisterHistCollection( Hists_gen_v100 = new ComparisonHists("gen_v100") );

  RegisterHistCollection( Hists_gen_recgen = new ComparisonHists("gen_recgen") );
  RegisterHistCollection( Hists_gen_recgen_NOdr = new ComparisonHists("gen_recgen") );
 
 
  RegisterHistCollection( Hists_gen_rec_sel = new ComparisonHists("gen_rec_sel") );
  RegisterHistCollection( Hists_gen_rec_notgen_sel = new ComparisonHists("gen_rec_notgen_sel") );
  RegisterHistCollection( Hists_gen_gen_notrec_sel = new ComparisonHists("gen_gen_notrec_sel") );

  //reco level selections
  RegisterHistCollection( Hists_reco_basic = new ComparisonHistsReco("reco_basic") );
  RegisterHistCollection( Hists_reco_noVeto_no2nd = new ComparisonHistsReco("reco_noVeto_no2nd") );
  RegisterHistCollection( Hists_reco_noVeto = new ComparisonHistsReco("reco_noVeto") );

  RegisterHistCollection( Hists_reco = new ComparisonHistsReco("reco") );
  RegisterHistCollection( Hists_reco_mass = new ComparisonHistsReco("reco_mass") );
  RegisterHistCollection( Hists_reco_mass_highPU = new ComparisonHistsReco("reco_mass_highPU") );
  RegisterHistCollection( Hists_reco_mass_lowPU = new ComparisonHistsReco("reco_mass_lowPU") );
  RegisterHistCollection( Hists_reco_mass_300to400 = new ComparisonHistsReco("reco_mass_300to400") );
  RegisterHistCollection( Hists_reco_mass_400to500 = new ComparisonHistsReco("reco_mass_400to500") );
  RegisterHistCollection( Hists_reco_mass_500toInf = new ComparisonHistsReco("reco_mass_500toInf") );
  RegisterHistCollection( Hists_reco_mass2 = new ComparisonHistsReco("reco_mass2") );


  RegisterHistCollection( Hists_reco_recgen = new ComparisonHistsReco("reco_recgen") );
  RegisterHistCollection( Hists_reco_recgen_NOdr = new ComparisonHistsReco("reco_recgen_NOdr") );
 
  RegisterHistCollection( Hists_reco_rec_notgen_sel_dilept = new ComparisonHistsReco("reco_rec_notgen_sel_dilept") );

  RegisterHistCollection( Hists_grooming = new GroomingHists("grooming") );
  RegisterHistCollection( Hists_W = new WHists("W") );
 
  //unfolding sidebands
  RegisterHistCollection( Hists_signal = new ComparisonHistsReco("signal") );
  RegisterHistCollection( Hists_SideBand1 = new ComparisonHistsReco("SideBand1") );
  RegisterHistCollection( Hists_SideBand2 = new ComparisonHistsReco("SideBand2") );
  RegisterHistCollection( Hists_SideBand3 = new ComparisonHistsReco("SideBand3") );
  RegisterHistCollection( Hists_SideBand4 = new ComparisonHistsReco("SideBand4") );
  RegisterHistCollection( Hists_SideBand5 = new ComparisonHistsReco("SideBand5") );

  //event hists
  RegisterHistCollection( Hists_event = new EventHists("event_fullselection") );
  RegisterHistCollection( Hists_jet = new JetHists("jet_fullselection") ) ;
  RegisterHistCollection( Hists_event_basic = new EventHists("event_basicselection") ); 
  RegisterHistCollection( Hists_jet_basic = new JetHists("jet_basicselection") ) ;

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

void SubstructurePostSelectionCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );
}

void SubstructurePostSelectionCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer

  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );

  
  //get the metadata input tree with the previous selection steps
  m_intree = GetInputMetadataTree("SelectionTree");
  m_intree->SetBranchAddress("recsel_basic" , &m_recsel_basic);
  m_intree->SetBranchAddress("gensel_basic" , &m_gensel_basic);
  m_intree->SetBranchAddress("rec_presel" , &m_rec_presel);
  m_intree->SetBranchAddress("gen_presel" , &m_gen_presel);
  m_intree->SetBranchAddress("gen_presel_hadron" , &m_gen_presel_hadron);

}

void SubstructurePostSelectionCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError ) 
{
 
  // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results

  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset
  AnalysisCycle::ExecuteEvent( id, weight );

 
  m_intree->GetEvent(id.GetEventTreeEntry());  //get input tree with basic selections
  SetDefaults();                                //set default values for all private members

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  if(calc->GetJets()->size()>=12){
    std::cout << "run: " << calc->GetRunNum() << "   lb: " << calc->GetLumiBlock() << "  event: " << calc->GetEventNum() << "   N(jets): " << calc->GetJets()->size() << std::endl;
  }


  //===========================================
  //scale factors and systematic variations
  //===========================================
  // if geninfo == true: store and plot events passing the genselection, particle level plots 
  bool geninfo = false; 

  if (id.GetVersion().Contains("TT"))geninfo = true; 

  if(!mttbar_gen_selection->passSelection())  throw SError( SError::SkipEvent );
 
  // b tagging scale factor
  if(m_bsf && m_addGenInfo) {
    calc->ProduceWeight(m_bsf->GetWeight());
    calc->ProduceRecWeight(m_bsf->GetWeight());
  }

  // Ele30_OR_PFJet320 trigger Scale Factor
  if(m_applyEleORJetTriggerSF && !calc->IsRealData() && calc->GetJets()->size()){
    calc->ProduceWeight( m_lsf->GetElectronORJetTrigWeight() );
    calc->ProduceRecWeight( m_lsf->GetElectronORJetTrigWeight() );
  }

  m_cleaner = new Cleaner();
  m_cleaner->SetJECUncertainty(m_jes_unc);

  // settings for jet correction uncertainties
  if (m_sys_unc==e_JEC){
    if (m_sys_var==e_Up) m_cleaner->ApplyJECVariationUp();
    if (m_sys_var==e_Down) m_cleaner->ApplyJECVariationDown();
  }
  if (m_sys_unc==e_JER){
    if (m_sys_var==e_Up){// m_cleaner->ApplyJERVariationUp();
      m_cleaner->ApplyfatJERVariationUp(); 
    }
    if (m_sys_var==e_Down){// m_cleaner->ApplyJERVariationDown();
      m_cleaner->ApplyfatJERVariationDown(); 
    }
  }

  //if(!bcc->isRealData && bcc->jets) m_cleaner->JetEnergyResolutionShifter();
  if(!bcc->isRealData && bcc->jets) m_cleaner->JetEnergyResolutionShifterFat();
  // delete m_cleaner;


  //std::sort(bcc->topjets->begin(), bcc->topjets->end(), HigherPt()); //sort after resolution shifts

  //===========================================
  //gen level selections
  //===========================================

  //============ hadronic selection =========== 
  if(id.GetVersion().Contains("TTbar"))
    {

      //matching 
      bool matched = true;
      if(m_matching_gen == "Matched" && !GenMatchedSelectionHadron->passSelection()) matched = false;
      if(m_matching_gen == "MisMatched" && !GenMatchedSelectionHadron->passInvertedSelection() ) matched = false;

      if(matched
	 && m_gen_presel_hadron
	 && GenJetSelection_400_300->passSelection()
	 && GenJetVeto_30->passSelection())
	{
	  Hists_hadronic_gen->Fill();
	  if(GenMassSelection->passSelection()) Hists_hadronic_gen_mass->Fill();
	}
    }
  
  //============ lepton + jets selection =======
  if(geninfo)
    {
      //  calc->PrintGenParticles("");
      //matching 
      bool matched = true;
     

      //plot presel
      if(m_gen_presel)
	{
	  if(m_matching_gen == "Matched" && !GenMatchedSelection->passSelection()) matched = false;
	  if(m_matching_gen == "MisMatched" && !GenMatchedSelection->passInvertedSelection() ) matched = false;
	  Hists_gen_presel->Fill();
	}
    
      // basic selections (all selections will be saved in the output tree)  
      if(m_gensel_basic && matched)
	{
	  /*
	  int NLep = 0;
	  int NLep2 = 0;
	  int Ntop = 0;
	  std::vector<GenParticle>* genparts = bcc->genparticles;
	  for(unsigned int j = 0; j < genparts->size(); j++){
	    if(genparts->at(j).status() < 3) break; //loop only over status three particles
	    if(fabs(genparts->at(j).pdgId()) == 11 || fabs(genparts->at(j).pdgId()) == 13 ) NLep ++;
	  }

	  for(unsigned int j = 0; j < genparts->size(); j++){
	    if(genparts->at(j).status() < 3) continue; //loop only over status three particles
	    if(fabs(genparts->at(j).pdgId()) == 6 ) Ntop ++;
	    if(fabs(genparts->at(j).pdgId()) == 11 || fabs(genparts->at(j).pdgId()) == 13 ) NLep2 ++;
	  }
	  if(NLep != 1) cout << "NLep = " << NLep << endl;
	  if(NLep2 != 1) cout << "NLep2 = " << NLep2 << endl;
	  if(Ntop != 2){ 
	    cout << "Ntop = " << Ntop << endl;
	  }
	  */

	  //	  if(NLep != 1 || NLep2 != 1) 
	  // calc->PrintGenParticles("");
	  //grooming 
	  if(m_groom_genjets) SoftDropNGenJetsWithParts(3);

	  if(GenJetSelection_300_100->passSelection()) m_gensel_300_100 = true;
	  if(GenJetSelection_400_100->passSelection()) m_gensel_400_100 = true;
	  if(GenJetSelection_500_100->passSelection()) m_gensel_500_100 = true;
	  if(GenJetSelection_300_150->passSelection()) m_gensel_300_150 = true;
	  //just tests
	  if(m_gensel_300_100) Hists_gen_noVeto_no2nd_300->Fill();
	  if(m_gensel_500_100) Hists_gen_noVeto_no2nd_500->Fill();

	  if(GenJetVeto_150->passSelection()) m_gen_jet_veto = true;
	  if(BoostedGenSelection->passSelection()) m_gen_boosted = true;
	  if(GenMassSelection->passSelection()) m_gen_mass_sel = true;
	  if(GenMassSelectionLepton->passSelection()) m_gen_mass_sel2 = true;
	  if(GenMassSelectionLeptonAndNeutrino->passSelection()) m_gen_mass_sel3 = true;
	}

      //main selection!!!
      if(m_gensel_400_100 && m_gensel_basic)
	{
	  Hists_gen_noVeto_no2nd->Fill();

	  if(GenJetSelection_400_150->passSelection())
	    {
	      m_gensel_400_150=true;
	      Hists_gen_noVeto->Fill();  

	      //different vetos
	      if(GenJetVeto_100->passSelection())
		{
		  m_gensel_v100=true;
		  Hists_gen_v100->Fill();
		}
	      if(GenJetVeto_150->passSelection())
		{
		  m_gensel_v150=true;
		  Hists_gen->Fill();
		  // Hists_truth->Fill();

		  if(m_gen_boosted) //DR(jet2, lepton) < 1.2
		    {   
		      Hists_gen_dr->Fill();

		      if(m_gen_mass_sel) //mjet1 > mjet2
			Hists_gen_mass->Fill(); 
		      if(m_gen_mass_sel2)//mjet1 > m(jet2 + lepton)
			{      
			  Hists_gen_mass2->Fill();
			  m_gensel = true;
			}
		      if(m_gen_mass_sel3)//mjet1 > m(jet2 + lepton + neutrino)
			Hists_gen_mass3->Fill(); 
    
		    }
		}
	    }
	}
  
    }//geninfo

  
  //===========================================
  //reco level selections
  //===========================================
 
  //matching 
  bool matched = true;
  if(m_matching_reco == "Matched" && !RecoMatchedSelection->passSelection()) matched = false;
  if(m_matching_reco == "MisMatched" && !RecoMatchedSelection->passInvertedSelection() ) matched = false;

 
  //======== basic selection ===================
  //now apply the HTlep cut
 
  if( HTlep_selection->passSelection() 
      && m_recsel_basic
      && BSelT->passSelection()
      && matched)
    {
      //grooming 
      if(m_groom_recojets) SoftDropNTopJets(3);

      Hists_jet_basic->Fill();
      Hists_event_basic->Fill();
      Hists_reco_basic->Fill();

      //===========================additional selections==========================
      //first selection step includes pt_leading > 300 GeV, pt_2nd > 100 GeV and no veto. All events passing this selection are stored in the output tree 
      if(TopJetSelection_300_100->passSelection() && RecoLeptonSelection->passSelection())
	{
	
	  m_recosel_300_100 = true;
	  if(TopJetSelection_400_100->passSelection()) m_recosel_400_100 = true; //pt_leading >  400 GeV
	  if(TopJetSelection_300_150->passSelection()) m_recosel_300_150 = true;                       //pt_leading >  300 GeV, pt_2nd > 150 GeV

	  //control distributions for the different pt bins
	  if(RecoJetVeto->passSelection() 
	     && BoostedRecoSelection->passSelection()
	     && RecoMassSelection->passSelection())
	    { 
	      if(m_recosel_300_150 && !TopJetSelection_400_150->passSelection()) 
		Hists_reco_mass_300to400->Fill();
	      if(TopJetSelection_400_150->passSelection() && !TopJetSelection_500_150->passSelection())
		Hists_reco_mass_400to500->Fill();
	      if(TopJetSelection_500_150->passSelection()) 
		Hists_reco_mass_500toInf->Fill();
	    }
		

	  //Set reco selections

	  if(BoostedRecoSelection->passSelection()) m_reco_boosted = true;
	  if(RecoMassSelection->passSelection()) m_reco_mass_sel = true;
	  if(RecoMassSelectionLepton->passSelection()) m_reco_mass_sel2 = true;
	  if(RecoJetVeto->passSelection()) m_reco_jet_veto = true;
	 
	  //=================main selection====================
	  if(TopJetSelection_400_150->passSelection() && RecoJetVeto->passSelection())  //pt_leading > 400 GeV, pt_2nd > 150 GeV, veto
	    {    
	      m_recosel_400_150=true;
	      Hists_reco->Fill();
 
	      if(m_reco_boosted && m_reco_mass_sel)                  //DR(jet2, lepton) > 1.2 
		{                             
		  m_recosel = true;                          
		  Hists_reco_mass->Fill();

		  if((calc->GetCAJets()->size()>2 && calc->GetCAJets()->at(2).pt()> 300) || 
		     (calc->GetCAJets()->size()>2 && calc->GetCAJets()->at(2).pt()>calc->GetCAJets()->at(1).pt())||
		     (calc->GetCAJets()->size()>1 && calc->GetCAJets()->at(1).pt()>calc->GetCAJets()->at(0).pt()) ){

		    cout << "run: " << calc->GetRunNum() << "   lb: " << calc->GetLumiBlock() << "  event: " << calc->GetEventNum()<<endl;
		    cout << "----------------------" <<endl;
		    cout << "1st jet    pt: " << calc->GetCAJets()->at(0).pt() << "eta" << calc->GetCAJets()->at(0).eta()<<endl;
		    cout << "2nd jet    pt: " << calc->GetCAJets()->at(1).pt() << "eta" << calc->GetCAJets()->at(1).eta()<<endl;
		    cout << "3rd jet    pt: " << calc->GetCAJets()->at(2).pt() << "eta" << calc->GetCAJets()->at(2).eta()<<endl;
		  }


		  if(calc->GetPrimaryVertices()->size() > 15) Hists_reco_mass_highPU->Fill();
		  if(calc->GetPrimaryVertices()->size() < 15) Hists_reco_mass_lowPU->Fill();
		  if(m_monitoring) Hists_grooming->Fill();
		  Hists_jet->Fill();
		  Hists_event->Fill();
		
		  if(geninfo){
		
		    if(m_gensel)
		      {
			Hists_gen_recgen->Fill();
			Hists_reco_recgen->Fill();
		      }
		    else Hists_gen_rec_notgen_sel->Fill();			  
 
		    //if(!ttbar->IsTopHadronicDecay() && !ttbar->IsTopHadronicDecay()) 
		    //Hists_reco_rec_notgen_sel_dilept->Fill();
		    //delete ttbar;	  
		  }
		  if(m_reco_mass_sel2) Hists_reco_mass2->Fill();
		}

	      if(geninfo 
		 && m_reco_mass_sel 
		 && m_gensel_v150 
		 && m_gen_mass_sel2)//----mjet1 > mjet2 (but no DR cut)!!!
		{				
		  Hists_gen_recgen_NOdr->Fill(); 
		  Hists_reco_recgen_NOdr->Fill();	    
		}
		
	    }

	  //=========fill sideband histos =========
	  Int_t NSel = 0;
	  if(m_recosel_400_100) NSel++; // pt1 > 400, pt2 > 100
	  if(m_recosel_300_150) NSel++;//pt1 > 300, pt2 > 150
	  if(m_reco_jet_veto) NSel++; 
	  if(m_reco_boosted) NSel++;
	  if(m_reco_mass_sel) NSel++;
	  
	  //signal region
	  if(NSel == 5) Hists_signal->Fill();
	  
	  //sidebands
	  if(NSel == 4 && !m_recosel_400_100) Hists_SideBand1->Fill();
	  if(NSel == 4 && !m_recosel_300_150) Hists_SideBand2->Fill();
	  if(NSel == 4 && !m_reco_jet_veto && calc->GetCAJets()->at(2).pt() < 200) Hists_SideBand3->Fill();
	  if(NSel == 4 && !m_reco_boosted) Hists_SideBand4->Fill();
	  if(NSel == 4 && !m_reco_mass_sel) Hists_SideBand5->Fill();



	  //===============================================
	  //declare output tree variables
	  //===============================================

	  if(calc->GetCAJets()->size()>1) 
	    m_recojet_pt1 = calc->GetCAJets()->at(0).pt();
	  if(calc->GetCAJets()->size()>1) 
	    m_recojet_pt2 = calc->GetCAJets()->at(1).pt();
	  if(calc->GetCAJets()->size()>2) 
	    m_recojet_pt3 = calc->GetCAJets()->at(2).pt();

	  m_jetmass_reco= calc->GetCAJets()->at(0).v4().mass();
	  m_jetmass_reco2= calc->GetCAJets()->at(1).v4().mass();
	
	  if (calc->GetCAJets()->size() && calc->GetCAJets()->at(0).subjets().size()>2) 
	    m_Nsubjets=true;
	  m_NPV = calc->GetPrimaryVertices()->size();
	  
	}
    }
 
  if(geninfo)
    {    
      //set pt
      if(bcc->cagenjets->size()>0) 
	m_genjet_pt1 = bcc->cagenjets->at(0).pt();
      if(bcc->cagenjets->size()>1) 
	m_genjet_pt2 = bcc->cagenjets->at(1).pt();
      if(bcc->cagenjets->size()>2)
	m_genjet_pt3 = bcc->cagenjets->at(2).pt();
   
      //set masses
      if( bcc->cagenjets->size()>0) 
	m_jetmass_gen = bcc->cagenjets->at(0).v4().mass();
      if( bcc->cagenjets->size()>1) 
	m_jetmass_gen2 = bcc->cagenjets->at(1).v4().mass();

      TTbarGen *ttbar = new TTbarGen(bcc);
      m_top_pt = ttbar->Top().pt();
      m_antitop_pt = ttbar->Antitop().pt();
      delete ttbar;

    }

  
  //set weights that will be written to the output tree and needed for the unfolding
  m_weight = calc->GetWeight();
  m_lumiweight = weight; //just a check
  m_recweight = calc->GetRecWeight();
  m_genweight = calc->GetGenWeight();
 



  //===========================================
  //Fill the output tree
  //===========================================
  if(m_recosel_300_100 || m_gensel_300_100)  m_tree->Fill();
 
  return;
}



void SubstructurePostSelectionCycle::SetDefaults(){
  //set default selections to false
  m_gensel = false;
  m_gensel_v100 = false;
  m_gensel_v150 = false;
  m_gensel_400_100 = false;
  m_gensel_300_100 = false;
  m_gensel_500_100 = false;
  m_gen_mass_sel = false;
  m_gen_mass_sel2 = false;
  m_gen_mass_sel3 = false;
  m_gen_boosted = false;
  m_gensel_400_150 = false;
  m_gensel_300_150 = false;
 
  m_Nsubjets = false;
  m_reco_jet_veto = false;
  m_gen_jet_veto= false;

  m_reco_mass_sel = false;
  m_reco_mass_sel2 = false;
  m_reco_boosted = false;
  m_recosel_400_100 = false;
  m_recosel_300_100 = false;
  m_recosel = false;
  m_recosel_300_150 = false;

  m_jetmass_reco = 0.;
  m_jetmass_reco2 = 0.;
  m_jetmass_gen = 0.;
  m_jetmass_gen2 = 0.;
  m_NPV = 0;

  m_genjet_pt1 = 0.;
  m_genjet_pt2 = 0.;
  m_genjet_pt3 = 0.;
  m_recojet_pt1 = 0.;
  m_recojet_pt2 = 0.;
  m_recojet_pt3 = 0.;
  m_top_pt = 0.;
  m_hadr_top_pt = 0.;
  m_antitop_pt = 0.;
}

void SubstructurePostSelectionCycle::SetOutputTreeBranches(){

  m_tree->Branch("weight" , &m_weight , "weight/D");
  m_tree->Branch("lumiweight" , &m_lumiweight , "lumiweight/D");
  m_tree->Branch("recweight" , &m_recweight , "recweight/D");
  m_tree->Branch("genweight" , &m_genweight , "genweight/D");

  m_tree->Branch("jetmass_reco" , &m_jetmass_reco , "jetmass_reco/D");
  m_tree->Branch("jetmass_reco2" , &m_jetmass_reco2 , "jetmass_reco2/D");
  m_tree->Branch("jetmass_gen" , &m_jetmass_gen , "jetmass_gen/D");
  m_tree->Branch("jetmass_gen2" , &m_jetmass_gen2 , "jetmass_gen2/D");

  m_tree->Branch("Nsubjets" , &m_Nsubjets , "Nsubjets/O");
  m_tree->Branch("NPV" , &m_NPV , "NPV/I");

  m_tree->Branch("genjet_pt1" , &m_genjet_pt1 , "genjet_pt1/D");
  m_tree->Branch("genjet_pt2" , &m_genjet_pt2 , "genjet_pt2/D");
  m_tree->Branch("genjet_pt3" , &m_genjet_pt3 , "genjet_pt3/D");
  m_tree->Branch("recojet_pt1" , &m_recojet_pt1 , "recojet_pt1/D");
  m_tree->Branch("recojet_pt2" , &m_recojet_pt2 , "recojet_pt2/D"); 
  m_tree->Branch("recojet_pt3" , &m_recojet_pt3 , "recojet_pt3/D");
 
  m_tree->Branch("hadr_top_pt" , &m_hadr_top_pt , "hadr_top_pt/D");
  m_tree->Branch("top_pt" , &m_top_pt , "top_pt/D");
  m_tree->Branch("antitop_pt" , &m_antitop_pt , "antitop_pt/D");

  m_tree->Branch("gensel_400_100" , &m_gensel_400_100, "gensel_400_100/O");
  m_tree->Branch("recosel_400_100" , &m_recosel_400_100, "recosel_400_100/O");
  m_tree->Branch("gensel_300_100" , &m_gensel_300_100, "gensel_300_100/O");
  m_tree->Branch("recosel_300_100" , &m_recosel_300_100, "recosel_300_100/O");
  m_tree->Branch("gensel" , &m_gensel , "gensel/O");
  m_tree->Branch("gensel_v100" , &m_gensel_v100 , "gensel_v100/O");
  m_tree->Branch("gensel_v150" , &m_gensel_v150 , "gensel_v150/O");

  m_tree->Branch("gensel_300_150" , &m_gensel_300_150 , "gensel_300_150/O");
  m_tree->Branch("gensel_400_150" , &m_gensel_400_150 , "gensel_400_150/O");
  m_tree->Branch("recosel" , &m_recosel , "recosel/O");
  m_tree->Branch("recosel_300_150" , &m_recosel_300_150 , "recosel_300_150/O");

  m_tree->Branch("reco_jet_veto" , &m_reco_jet_veto , "reco_jet_veto/O");
  m_tree->Branch("gen_jet_veto" , &m_gen_jet_veto , "gen_jet_veto/O");

  m_tree->Branch("gen_mass_sel" , &m_gen_mass_sel, "gen_mass_sel/O");
  m_tree->Branch("gen_mass_sel2" , &m_gen_mass_sel2 , "gen_mass_sel2/O");
  m_tree->Branch("gen_mass_sel3" , &m_gen_mass_sel3 , "gen_mass_sel3/O");
  m_tree->Branch("gen_boosted" , &m_gen_boosted , "gen_boosted/O");
  m_tree->Branch("reco_mass_sel" , &m_reco_mass_sel , "reco_mass_sel/O");
  m_tree->Branch("reco_mass_sel2" , &m_reco_mass_sel2 , "reco_mass_sel2/O");
  m_tree->Branch("reco_boosted" , &m_reco_boosted , "reco_boosted/O");
}
