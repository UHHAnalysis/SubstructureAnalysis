#include <iostream>

using namespace std;

// Local include(s):

#include "include/SubstructureSelectionCycle.h"
#include "include/JetMassTools.h"
#include "SFrameTools/include/TopFitCalc.h"
#include "ZprimeAnalysis/include/NeutrinoHists.h"
#include "SubstructureAnalysis/include/SubstructureSelectionModules.h"


ClassImp(SubstructureSelectionCycle );

SubstructureSelectionCycle::SubstructureSelectionCycle()
  : AnalysisCycle()
{

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

SubstructureSelectionCycle::~SubstructureSelectionCycle()
{
  // destructor
}

void SubstructureSelectionCycle::BeginCycle() throw( SError )
{
  // Start of the job, general set-up and definition of
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;
}

void SubstructureSelectionCycle::EndCycle() throw( SError )
{
  // clean-up, info messages and final calculations after the analysis


  // call the base cycle class for all standard methods
  // and a summary of the made selections
  AnalysisCycle::EndCycle();

  return;
}

void SubstructureSelectionCycle::BeginInputData( const SInputData& id ) throw( SError )
{
  // declaration of histograms and selections

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );


  m_tree = GetOutputMetadataTree("SelectionTree");
  m_tree->Branch("rec_presel" , &m_rec_presel , "rec_presel/O");
  m_tree->Branch("gen_presel" , &m_gen_presel , "gen_presel/O");
  m_tree->Branch("gen_presel_hadron" , &m_gen_presel_hadron , "gen_presel_hadron/O");
  m_tree->Branch("recsel_basic" , &m_recsel_basic , "recsel_basic/O");
  m_tree->Branch("gensel_basic" , &m_gensel_basic , "gensel_basic/O");


  if (id.GetVersion().Contains("TT")) m_logger << INFO << "Apply gen selections" << SLogger::endmsg; 


  // -------------------- set up the selections ---------------------------

  // cut out mttbar events for the 0_to_700 sample to not double-count them
  static Selection* mttbar_gen_selection = new Selection("Mttbar_Gen_Selection");
  if ( m_mttgencut && ((id.GetVersion() == "TTbar_0to700") || (id.GetVersion() == "TTbar") ) ) {
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

  second_selection->addSelectionModule(new NJetSelection(1,int_infinity(),150,2.4)); //leading jet with pt>150 GeV
  second_selection->addSelectionModule(new NBTagSelection(m_Nbtags_min,m_Nbtags_max)); //b tags from config file
  //second_selection->addSelectionModule(new HTlepCut(150));
  second_selection->addSelectionModule(new METCut(20));
  // second_selection->addSelectionModule(new TriangularCut());

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

  ////////////////////////////////////////
  Selection* genpreselection = new Selection("genpreselection");
  if(m_Electron_Or_Muon_Selection=="Electrons" || m_Electron_Or_Muon_Selection=="Electron" || m_Electron_Or_Muon_Selection=="Ele" || m_Electron_Or_Muon_Selection=="ELE") {
    genpreselection->addSelectionModule(new SubstructureGeneratorPreSelectionElectron(0));
  } else if(m_Electron_Or_Muon_Selection=="Muon" || m_Electron_Or_Muon_Selection=="Muons" || m_Electron_Or_Muon_Selection=="Mu" || m_Electron_Or_Muon_Selection=="MU") {
    genpreselection->addSelectionModule(new SubstructureGeneratorPreSelectionMuon(0));
  }
  genpreselection->addSelectionModule(new NCAGenJetSelection(1,int_infinity(),300,2.5));
  genpreselection->addSelectionModule(new NCAGenJetSelection(2,int_infinity(),100,2.5));
  RegisterSelection(genpreselection);

  Selection* genpreselection_hadron = new Selection("genpreselection_hadron");
  genpreselection_hadron->addSelectionModule(new SubstructureGeneratorPreSelectionHadron());
  genpreselection_hadron->addSelectionModule(new NCAGenJetSelection(2,int_infinity(), 300, 2.5));
  RegisterSelection(genpreselection_hadron);
 

  Selection* LeptonSelection = new Selection("LeptonSelection");
  LeptonSelection->addSelectionModule(new GenLeptonSelection(45,2.1));
  RegisterSelection(LeptonSelection);

  Selection* Jetsel_350 = new Selection("Jetsel_350");
  Jetsel_350->addSelectionModule(new NCAGenJetSelection(1,int_infinity(),350,2.5));
  Jetsel_350->addSelectionModule(new NCAGenJetSelection(2,int_infinity(),100,2.5));
  //  Jetsel_350->addSelectionModule(new NCAGenJetSelection(350,100,2.5));
  RegisterSelection(Jetsel_350);

  Selection* Jetselveto_150 = new Selection("Jetselveto_150");
  Jetselveto_150->addSelectionModule(new NCAGenJetveto(150,2,2.5));
  RegisterSelection(Jetselveto_150);



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

  ////////////////////////////////////////
    
  // ---------------- set up the histogram collections --------------------

  // histograms without any cuts
  RegisterHistCollection( new HypothesisHists("Chi2", m_chi2discr) );
  RegisterHistCollection( new HypothesisHists("BestPossible", m_bpdiscr) );
  RegisterHistCollection( new HypothesisHists("SumDeltaR", m_sumdrdiscr) );
  RegisterHistCollection( new HypothesisHists("CorrectMatch", m_cmdiscr) );

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
  RegisterHistCollection( new NeutrinoHists("Neutrino_Postsel" , m_chi2discr) );
  // important: initialise histogram collections after their definition
  InitHistos();

  m_bp_chi2 = new HypothesisStatistics("b.p. vs. Chi2");
  m_bp_sumdr = new HypothesisStatistics("b.p. vs. SumDR");
  m_cm_chi2 = new HypothesisStatistics("matched vs. Chi2");
  m_cm_sumdr = new HypothesisStatistics("matched vs. SumDR");
  m_cm_bp = new HypothesisStatistics("matched vs. b.p.");

  return;

}

void SubstructureSelectionCycle::EndInputData( const SInputData& id ) throw( SError )
{

  AnalysisCycle::EndInputData( id );

  m_bp_chi2->PrintStatistics();
  m_bp_sumdr->PrintStatistics();
  m_cm_chi2->PrintStatistics();
  m_cm_sumdr->PrintStatistics();
  m_cm_bp->PrintStatistics();

  delete m_bp_chi2;
  delete m_bp_sumdr;
  delete m_cm_chi2;
  delete m_cm_sumdr;
  delete m_cm_bp;

  return;
}

void SubstructureSelectionCycle::BeginInputFile( const SInputData& id ) throw( SError )
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer

  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );

  //read in the selections applyed in the preselection cycle (not working!!!!!!!)
  //m_intree = GetInputMetadataTree("PreSelectionTree");
  //m_intree->SetBranchAddress("recsel_basic" , &m_rec_presel);
  //m_intree->SetBranchAddress("gensel_basic" , &m_gen_presel);

  return;
}

void SubstructureSelectionCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError )
{
 
  // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results

  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset

  AnalysisCycle::ExecuteEvent( id, weight);

  m_recsel_basic = false; 
  m_gensel_basic= false; 
  m_rec_presel = false; 
  m_gen_presel= false; 
  m_gen_presel_hadron= false; 

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  TopFitCalc* topcalc = TopFitCalc::Instance();
  topcalc->Reset();

  static Selection* mttbar_gen_selection = GetSelection("Mttbar_Gen_Selection");
  static Selection* genpreselection = GetSelection("genpreselection");
  static Selection* genpreselection_hadron = GetSelection("genpreselection_hadron");
  static Selection* LeptonSelection = GetSelection("LeptonSelection");
  static Selection* Jetsel_350 = GetSelection("Jetsel_350");
  static Selection* Jetselveto_150 = GetSelection("Jetselveto_150");
  static Selection* preselection = GetSelection("preselection");
  bool gensel=false;


  if(!mttbar_gen_selection->passSelection()) throw SError( SError::SkipEvent );

  //hadron generator selection only for POWHEG ttbar 
  if (id.GetVersion().Contains("TTbar")){
    m_gen_presel_hadron = genpreselection_hadron->passSelection();
  }

  //other generator selections on all ttbar samples
  if (id.GetVersion().Contains("TT")){
    if(genpreselection->passSelection()) m_gen_presel = true;
    if(genpreselection->passSelection() && LeptonSelection->passSelection() ) gensel = true;
  }
  m_gensel_basic = gensel;

  if(preselection->passSelection()) m_rec_presel = true;
 
  if(!m_rec_presel && !gensel && !m_gen_presel_hadron) throw SError( SError::SkipEvent );

 
  // control histograms
  if(m_rec_presel) FillControlHists("_Presel");
 
  static Selection* Ele30trig_selection = GetSelection("Ele30trig_selection");
  static Selection* PFJet320trig_selection = GetSelection("PFJet320trig_selection");
  static Selection* trig_selection = GetSelection("trig_selection");
  static Selection* first_selection = GetSelection("first_selection");
  static Selection* second_selection = GetSelection("second_selection");
  static Selection* trangularcut_selection = GetSelection("trangularcut_selection");
  static Selection* chi2_selection = GetSelection("chi2_selection");
  static Selection* matchable_selection = GetSelection("matchable_selection");
  static Selection* TopTagSel = GetSelection("TopTagSelection");

  m_cleaner = new Cleaner();
  m_cleaner->SetJECUncertainty(m_jes_unc);

  // settings for jet correction uncertainties
  if (m_sys_unc==e_JEC){
    if (m_sys_var==e_Up) m_cleaner->ApplyJECVariationUp();
    if (m_sys_var==e_Down) m_cleaner->ApplyJECVariationDown();
  }
  if (m_sys_unc==e_JER){
    if (m_sys_var==e_Up){
       m_cleaner->ApplyJERVariationUp();
      // m_cleaner->ApplyfatJERVariationUp(); 
    }
    if (m_sys_var==e_Down){
      m_cleaner->ApplyJERVariationDown();
      // m_cleaner->ApplyfatJERVariationDown(); 
    }
  }

  GenCleaner *gen_cleaner = new GenCleaner();


  // if(!(m_Electron_Or_Muon_Selection=="Muon") && !(id.GetVersion() == "TTJets_SemiLeptMG") && !(id.GetVersion() == "TTJets_FullLeptMG") && !(id.GetVersion() == "TTJets_HadronicMG") ){
    if(gensel) gen_cleaner->GenCaJetLeptonSubstractor(true);
    if(!bcc->isRealData && bcc->cagenjets) gen_cleaner->CAGenJetCleaner(10,2.5);
    // }
  
  if(bcc->pvs) m_cleaner->PrimaryVertexCleaner(4, 24., 2.);
  if(bcc->electrons) m_cleaner->ElectronCleaner_noIso(35,2.5, m_reversed_electron_selection,true);
  if(bcc->muons) m_cleaner->MuonCleaner_noIso(45,2.1);
  if(bcc->jets) m_cleaner->JetLeptonSubtractor(m_corrector,false);
  if(bcc->topjets) m_cleaner->TopJetLeptonSubtractor(m_corrector,true);
  if(bcc->topjets) m_cleaner->TopJetCleaner(10,2.5,false);
  if(!bcc->isRealData && bcc->jets) m_cleaner->JetEnergyResolutionShifter();
  // if(!bcc->isRealData && bcc->topjets) m_cleaner->JetEnergyResolutionShifterFat();

  //apply loose jet cleaning for 2D cut
  if(bcc->jets) m_cleaner->JetCleaner(25,double_infinity(),true);

  // control histograms
  bool reco_presel = preselection->passSelection();
  // if(reco_presel) FillControlHists("_Cleaned");

  if(m_veto_electron_trigger && Ele30trig_selection->passSelection() && !gensel && !m_gen_presel_hadron){
    throw SError( SError::SkipEvent );
  }

  bool triggerbit(0);
  if(!m_useORTriggerWithPFJet320) triggerbit = trig_selection->passSelection();
  else triggerbit = trig_selection->passSelection() || PFJet320trig_selection->passSelection();

  if(!triggerbit && !gensel && !m_gen_presel_hadron) throw SError( SError::SkipEvent );

  if(triggerbit) FillControlHists("_Cleaned");

  bool reco_sel1 = first_selection->passSelection();
  if(!reco_sel1 && !gensel && !m_gen_presel_hadron) throw SError( SError::SkipEvent );


    
  //apply tighter jet cleaning for further cuts and analysis steps
  if(bcc->jets) m_cleaner->JetCleaner(50,2.5,true);

  //remove all taus from collection for HTlep calculation
  if(bcc->taus) m_cleaner->TauCleaner(double_infinity(),0.0);

  bool reco_sel2 = second_selection->passSelection();
  if(!reco_sel2 && !gensel && !m_gen_presel_hadron) throw SError( SError::SkipEvent );

  bool reco_sel3(0);
  if(!m_reversed_electron_selection) {
    if(trangularcut_selection->passSelection()) reco_sel3 = true;
  } else {
    if(trangularcut_selection->passInvertedSelection()) reco_sel3 = true;
  }
  if(!reco_sel3 && !gensel  && !m_gen_presel_hadron) throw SError( SError::SkipEvent );

  if(reco_presel && !(m_veto_electron_trigger && Ele30trig_selection->passSelection()) && triggerbit && reco_sel1 && reco_sel2 && reco_sel3){
    FillControlHists("_Postsel");
    m_recsel_basic = true;
  }


  //do reconstruction here
  //if(!bcc->recoHyps) cout<<"no Hyp list"<<endl;
  // if(first_selection->passSelection() && second_selection->passSelection()){
  

    //  if(TopTagSel->passSelection()){

    //   topcalc->CalculateTopTag();
      //topcalc->FillHighMassTTbarHypotheses();
 
  //  } else{
  //    topcalc->FillHighMassTTbarHypotheses();
  //  }

  //  m_chi2discr->FillDiscriminatorValues();

  //  m_bpdiscr->FillDiscriminatorValues();

  //  m_sumdrdiscr->FillDiscriminatorValues();
 
    //  m_cmdiscr->FillDiscriminatorValues();
 
    //if(!chi2_selection->passSelection()) throw SError( SError::SkipEvent );
    //if(!matchable_selection->passSelection()) throw SError( SError::SkipEvent );

  //  ReconstructionHypothesis *hyp = m_chi2discr->GetBestHypothesis();

    //double mttbar= (hyp->toplep_v4()+hyp->tophad_v4()).M();
    //double nu_pz = hyp->neutrino_v4().pz();
    //std::cout << "event: " << bcc->event << " mttbar : " << mttbar << " chi2 : " << hyp->discriminator("Chi2") << " chi2 tlep: " << hyp->discriminator("Chi2_tlep") <<" chi2 thad: " << hyp->discriminator("Chi2_thad") <<std::endl;

    //std::cout << " leptonic side: top pt : " << hyp->toplep_v4().Pt() << " top m : " << hyp->toplep_v4().mass() << " lepton pt : " << hyp->lepton().pt() << " nu pz : " << nu_pz << " jet pt : " << bcc->jets->at( hyp->blep_index() ).pt() <<" njets : " << hyp->toplep_jets_indices().size() << std::endl;
    //std::cout << " hadronic side: top pt : " << hyp->tophad_v4().Pt() << " top m : " << hyp->tophad_v4().mass() << " njets : " << hyp->tophad_jets_indices().size();
    // for(unsigned int i=0; i< hyp->tophad_jets_indices().size(); ++i){
    // std::cout << " jet " << i<< " pt : " << bcc->jets->at( hyp->tophad_jets_indices().at(i) ).pt() ;
    // }
    // std::cout << std::endl;

  //  ReconstructionHypothesis *bp_hyp = m_bpdiscr->GetBestHypothesis();
    //double bp_mttbar = (bp_hyp->toplep_v4()+bp_hyp->tophad_v4()).M();
    //std::cout << " bp mttbar : " << bp_mttbar << " bp deltaR : " << bp_hyp->discriminator("BestPossible") << std::endl;

  //  ReconstructionHypothesis *sdr_hyp = m_sumdrdiscr->GetBestHypothesis();
    // ReconstructionHypothesis *cm_hyp = m_cmdiscr->GetBestHypothesis();
    
  /* 
    // control histograms
    if(preselection->passSelection() && !(m_veto_electron_trigger && Ele30trig_selection->passSelection()) && triggerbit && first_selection->passSelection() && second_selection->passSelection() && trangularcut_selection->passSelection()){ 
      FillControlHists("_Postsel");
      m_recsel_basic = true;
    }
  */


    // neutrino hists
    //   BaseHists* nuhists = GetHistCollection("Neutrino_Postsel");
    //   nuhists->Fill();

    // get the histogram collections
    //  BaseHists* Chi2Hists = GetHistCollection("Chi2");
    //  BaseHists* BPHists = GetHistCollection("BestPossible");
    //  BaseHists* SumDRHists = GetHistCollection("SumDeltaR");
    // BaseHists* CorrectMatchHists = GetHistCollection("CorrectMatch");

    //  Chi2Hists->Fill();
    //  BPHists->Fill();
    //  SumDRHists->Fill();
    // if(cm_hyp) CorrectMatchHists->Fill();

    //  m_bp_chi2->FillHyps(bp_hyp,hyp);
    //    m_bp_sumdr->FillHyps(bp_hyp,sdr_hyp);

    /* if(cm_hyp) {
       m_cm_chi2->FillHyps(cm_hyp,hyp);
       m_cm_sumdr->FillHyps(cm_hyp,sdr_hyp);
       m_cm_bp->FillHyps(cm_hyp, bp_hyp);
       }*/
    //}
  //calc->PrintEventContent();

 
  m_tree->Fill();
  WriteOutputTree();


  delete m_cleaner;
  delete gen_cleaner;
 
  return;

}

void SubstructureSelectionCycle::FillControlHists(TString postfix)
{
  // fill some control histograms, need to be defined in BeginInputData

  BaseHists* eventhists = GetHistCollection((std::string)("Event"+postfix));
  BaseHists* jethists = GetHistCollection((std::string)("Jets"+postfix));
  BaseHists* elehists = GetHistCollection((std::string)("Electron"+postfix));
  BaseHists* muonhists = GetHistCollection((std::string)("Muon"+postfix));
  BaseHists* tauhists = GetHistCollection((std::string)("Tau"+postfix));
  BaseHists* topjethists = GetHistCollection((std::string)("TopJets"+postfix));

  eventhists->Fill();
  jethists->Fill();
  elehists->Fill();
  muonhists->Fill();
  tauhists->Fill();
  topjethists->Fill();

}
