
#include "include/GroomingHists.h"
#include "SFrameTools/include/EventCalc.h"
#include "NtupleWriter/include/JetProps.h"
#include "../include/JetMassTools.h"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include <iostream>


using namespace std;

GroomingHists::GroomingHists(const char* name) : BaseHists(name)
{
  // named default constructor
   
}

GroomingHists::~GroomingHists(){
}


void GroomingHists::Init(){


 Book( TH1F("r","",50,0,500));
 Book( TH1F("ungroomed_jet_mass","m_{jet}",50,0,500));
 Book( TH1F("trimmed_jet_mass_default","m_{jet}",50,0,500));
 Book( TH1F("trimmed_jet_mass_opti","m_{jet}",50,0,500));
 Book( TH1F("pruned_jet_mass_opti","m_{jet}",50,0,500));
 Book( TH1F("filtered_jet_mass_opti","m_{jet}",50,0,500));
 Book( TH1F("pull_ungroomed_jet_mass", "m_{jet}^{rec}-m_{jet}^{gen}/m_{jet}^{gen}",100,-3,3));
 Book( TH1F("pull_trimmed_jet_mass_default", "m_{jet}^{rec}-m_{jet}^{gen}/m_{jet}^{gen}",100,-3,3));
 Book( TH1F("pull_trimmed_jet_mass_opti", "m_{jet}^{rec}-m_{jet}^{gen}/m_{jet}^{gen}",100,-3,3));
 Book( TH1F("pull_pruned_jet_mass_opti", "m_{jet}^{rec}-m_{jet}^{gen}/m_{jet}^{gen}",100,-3,3));
 Book( TH1F("pull_filtered_jet_mass_opti", "m_{jet}^{rec}-m_{jet}^{gen}/m_{jet}^{gen}",100,-3,3));
}


void GroomingHists::Fill(){

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  bool IsRealData = calc->IsRealData();
  LuminosityHandler* lumih = calc->GetLumiHandler();

  // important: get the event weight
  double weight = calc->GetWeight();

  std::vector<TopJet> *recojets = calc->GetCAJets();
  sort(recojets->begin(), recojets->end(), HigherPt());



  //set up grooming methods
  JetProps jp(&recojets->at(0), calc->GetPFParticles() );

  std::vector<fastjet::PseudoJet> jets = jp.GetFastJet(2.0);


  //filtering
  fastjet::Filter filter3(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3)); //recommended


  //pruning
  double zcut = 0.1;
  double rcut_factor = recojets->at(0).v4().M()/recojets->at(0).pt();  //recommended
  fastjet::Pruner pruner(fastjet::cambridge_algorithm, zcut, rcut_factor); 

  //trimming
  fastjet::Filter trimm_2_03(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03));//recommended
  fastjet::Filter trimm_2_003(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.003));//optimized

  double jec = 1/bcc->topjets->at(0).JEC_factor_raw();



  if(jets.size()){
    fastjet::PseudoJet filtered_jet3 = filter3(jets[0]);
    fastjet::PseudoJet pruned_jet = pruner(jets[0]);
    fastjet::PseudoJet trimmed_jet = trimm_2_03(jets[0]);
    fastjet::PseudoJet trimmed03_jet = trimm_2_003(jets[0]);
    
    Hist("ungroomed_jet_mass")->Fill(recojets->at(0).v4().M(),weight);
    Hist("trimmed_jet_mass_default")->Fill(trimmed_jet.m()*jec,weight);
    Hist("trimmed_jet_mass_opti")->Fill(trimmed03_jet.m()*jec,weight);
    Hist("pruned_jet_mass_opti")->Fill( pruned_jet.m()*jec,weight);
    Hist("filtered_jet_mass_opti")->Fill(filtered_jet3.m()*jec,weight);
    

    //pull distributions
    if(bcc->cagenjets){
      Hist("pull_ungroomed_jet_mass")->Fill(recojets->at(0).v4().M()/bcc->cagenjets->at(0).v4().M()-1,weight);
      Hist("pull_trimmed_jet_mass_default")->Fill(trimmed_jet.m()*jec/bcc->cagenjets->at(0).v4().M()-1,weight);
      Hist("pull_trimmed_jet_mass_opti")->Fill(trimmed03_jet.m()*jec/bcc->cagenjets->at(0).v4().M()-1,weight);
      Hist("pull_pruned_jet_mass_opti")->Fill( pruned_jet.m()*jec/bcc->cagenjets->at(0).v4().M()-1,weight);
      Hist("pull_filtered_jet_mass_opti")->Fill(filtered_jet3.m()*jec/bcc->cagenjets->at(0).v4().M()-1,weight);

    }
  }
}

void GroomingHists::Finish()
{
  // final calculations, like division and addition of certain histograms
  EventCalc* calc = EventCalc::Instance();
  bool IsRealData = calc->IsRealData(); 
}
