#include "include/WHists.h"

#include <iostream>


using namespace std;

WHists::WHists(const char* name) : BaseHists(name)
{
  // named default constructor
   
}

WHists::~WHists(){
}


void WHists::Init()
{
  TString object[] = {"WJet", "CAbJet" , "WJet1" , "WJet2" , "bJetLep", "bJetHad" , "WRecoLep", "WRecoHad", "Lepton",  "Neutrino"};

  for(unsigned int i = 0; i < (sizeof(object)/sizeof(*object)); i++){
    Book( TH1F("pt_"+object[i] , "p_{T} "+object[i]+" [GeV]", 50, 0, 500));
    Book( TH1F("mass_"+object[i] , "mass "+object[i]+" [GeV]", 100, 0, 500));
    Book( TH1F("eta_"+object[i] , "#eta "+object[i], 60, -3, 3));
    Book( TH1F("phi_"+object[i] , "#phi "+object[i], 64, -3.2, 3.2));
    Book( TH1F("dR_TopLep_"+object[i] , "#DeltaR( t_{lep}, "+object[i]+")", 50, 0, 5));
    Book( TH1F("dPhi_TopLep_"+object[i] , "#Delta#phi( t_{lep}, "+object[i]+")", 63, 0, 6.3));
  }

  Book( TH1F("pt_TopHad" , "p_{T} TopHad [GeV]", 50, 0, 500));
  Book( TH1F("mass_TopHad" , "mass TopHad [GeV]", 50, 0, 500));
  Book( TH1F("eta_TopHad" , "#eta TopHad", 60, -3, 3));
  Book( TH1F("phi_TopHad" , "#phi TopHad", 64, -3.2, 3.2));
  Book( TH1F("dR_TopLep_TopHad" , "#DeltaR( t_{lep}, TopHad)", 50, 0, 5));
    Book( TH1F("dPhi_TopLep_TopHad" , "#Delta#phi( t_{lep}, TopHad)", 63, 0, 6.3));

  Book( TH1F("pt_TopLep" , "p_{T} TopLep [GeV]", 50, 0, 500));
  Book( TH1F("mass_TopLep" , "mass TopLep [GeV]", 50, 0, 500));
  Book( TH1F("eta_TopLep" , "#eta TopLep", 60, -3, 3));
  Book( TH1F("phi_TopLep" , "#phi TopLep", 64, -3.2, 3.2));
  Book( TH1F("dR_TopLep_TopLep" , "#DeltaR( t_{lep}, TopLep)", 50, 0, 5));
  Book( TH1F("dPhi_TopLep_TopLep" , "#Delta#phi( t_{lep}, TopLep)", 63, 0, 6.3));

  Book( TH1F("dR_Lepton_bJetLep", "#DeltaR( lepton, bJetLep)", 50, 0, 5));

  Book( TH1F("dR_WJet_bJetHad", "#DeltaR( Wjet, bJetHad)", 50, 0, 5));
  Book( TH1F("dR_WJet_bJetLep", "#DeltaR( Wjet, bJetLep)", 50, 0, 5));
  Book( TH1F("dR_WJet_W1", "#DeltaR( Wjet, WJet1)", 50, 0, 5));
  Book( TH1F("dR_WJet_W2", "#DeltaR( Wjet, WJet2)", 50, 0, 5));
  Book( TH1F("dR_WJet_CAbJet", "#DeltaR( Wjet, CAbJet)", 50, 0, 5));

  Book( TH1F("dR_bJetHad_W1", "#DeltaR( bJetHad, WJet1)", 50, 0, 5));
  Book( TH1F("dR_bJetHad_WRecoHad", "#DeltaR( bJetHad, WJet2)", 50, 0, 5));
  Book( TH1F("dR_bJetHad_CAbJet", "#DeltaR( bJetHad, CAbJet)", 50, 0, 5));
  Book( TH1F("dR_bJetHad_bJetLep", "#DeltaR( bJetHad, bJetLep)", 50, 0, 5));

  Book( TH1F("dR_W1_W2", "#DeltaR( WJet1, WJet2)", 50, 0, 5));
  Book( TH1F("dR_Lepton_Neutrino", "#DeltaR( Lepton, Neutrino)", 50, 0, 5));
  Book( TH1F("dPhi_Lepton_Neutrino", "#Delta#phi( Lepton, Neutrino)", 50, 0, 5));  
  Book( TH1F("dEta_Lepton_Neutrino", "#Delta#eta( Lepton, Neutrino)", 50, 0, 5));
}

void WHists::Fill()
{
  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  bool IsRealData = calc->IsRealData();
  LuminosityHandler* lumih = calc->GetLumiHandler();

  // important: get the event weight
  double weight = calc->GetWeight();

  m_tt = new TTbarHemisphereReconstruction(bcc);

  bool Lept_reco = m_tt->ReconstructLeptonicHemisphere();
  bool Had_reco = m_tt->ReconstructHadronicHemisphere(); 

  //  if(m_tt->LeptonicHemisphereFound() && m_tt->HadronicHemisphereFound()){

  if(Lept_reco){
    //pt eta phi mass
    FillPtEtaPhiMass(m_tt->bJetLep(), "bJetLep", weight);
    FillPtEtaPhiMass(m_tt->TopLep(), "TopLep", weight);
    FillPtEtaPhi(m_tt->Lepton(), "Lepton", weight);
    FillPtEtaPhiMass(m_tt->Neutrino(), "Neutrino", weight);

    //dRplots
    Hist("dR_Lepton_bJetLep")->Fill(deltaR(m_tt->Lepton().v4(), m_tt->bJetLep().v4()),weight);
    Hist("dR_Lepton_Neutrino")->Fill(deltaR(m_tt->Lepton().v4(), m_tt->Neutrino()),weight);
    Hist("dPhi_Lepton_Neutrino")->Fill(deltaPhiAbs(m_tt->Lepton().phi(), m_tt->Neutrino().Phi()),weight);
    Hist("dEta_Lepton_Neutrino")->Fill(fabs(m_tt->Lepton().eta()-m_tt->Neutrino().Eta()),weight);

    if(Had_reco){
      //pt eta phi mass
      FillPtEtaPhiMass(m_tt->WJet(), "WJet", weight);
      FillPtEtaPhiMass(m_tt->CAbJet(), "CAbJet", weight);
      FillPtEtaPhiMass(m_tt->WJet1(), "WJet1", weight);
      FillPtEtaPhiMass(m_tt->WJet2(), "WJet2", weight);
      FillPtEtaPhiMass(m_tt->bJetHad(), "bJetHad", weight);

      FillPtEtaPhiMass(m_tt->TopHad(), "TopHad", weight);
      FillPtEtaPhiMass(m_tt->WRecoLep(), "WRecoLep", weight);
      FillPtEtaPhiMass(m_tt->WRecoHad(), "WRecoHad", weight);
     
      //dRplots
      Hist("dR_bJetHad_W1")->Fill(deltaR(m_tt->bJetHad().v4(), m_tt->WJet1().v4()),weight);
      Hist("dR_bJetHad_WRecoHad")->Fill(deltaR(m_tt->bJetHad().v4(), m_tt->WRecoHad()),weight);
      Hist("dR_bJetHad_CAbJet")->Fill(deltaR(m_tt->bJetHad().v4(), m_tt->CAbJet().v4()),weight);
      Hist("dR_bJetHad_bJetLep")->Fill(deltaR(m_tt->bJetHad().v4(), m_tt->bJetLep().v4()),weight);
       
      Hist("dR_WJet_bJetHad")->Fill(deltaR(m_tt->bJetHad().v4(), m_tt->WJet().v4()),weight);
      Hist("dR_WJet_bJetLep")->Fill(deltaR(m_tt->bJetLep().v4(), m_tt->WJet().v4()),weight);
      Hist("dR_WJet_W1")->Fill(deltaR(m_tt->WJet().v4(), m_tt->WJet1().v4()),weight);
      Hist("dR_WJet_W2")->Fill(deltaR(m_tt->WJet().v4(), m_tt->WJet2().v4()),weight);
      Hist("dR_WJet_CAbJet")->Fill(deltaR(m_tt->WJet().v4(), m_tt->CAbJet().v4()),weight);

      Hist("dR_W1_W2")->Fill(deltaR(m_tt->WJet2().v4(), m_tt->WJet1().v4()),weight);
    
    }
  } 

  delete m_tt;
}

void WHists::FillPtEtaPhiMass(Jet jet, TString postfix, double weight){
  Hist("pt_"+postfix)->Fill(jet.pt(), weight);
  Hist("eta_"+postfix)->Fill(jet.eta(), weight);
  Hist("phi_"+postfix)->Fill(jet.phi(), weight);
  Hist("mass_"+postfix)->Fill(jet.v4().mass(), weight);
  Hist("dR_TopLep_"+postfix)->Fill(deltaR(m_tt->TopLep(),jet.v4() ),weight);
  Hist("dPhi_TopLep_"+postfix)->Fill(deltaPhiAbs(m_tt->TopLep().Phi(), jet.phi()),weight); 
}

void WHists::FillPtEtaPhiMass(TopJet jet, TString postfix, double weight){
  Hist("pt_"+postfix)->Fill(jet.pt(), weight);
  Hist("eta_"+postfix)->Fill(jet.eta(), weight);
  Hist("phi_"+postfix)->Fill(jet.phi(), weight);
  Hist("mass_"+postfix)->Fill(jet.v4().mass(), weight);
  Hist("dR_TopLep_"+postfix)->Fill(deltaR(m_tt->TopLep(),jet.v4() ),weight);
  Hist("dPhi_TopLep_"+postfix)->Fill(deltaPhiAbs(m_tt->TopLep().Phi(), jet.phi()),weight); 
}

void WHists::FillPtEtaPhiMass(LorentzVector vec, TString postfix, double weight){
  Hist("pt_"+postfix)->Fill(vec.Pt(), weight);
  Hist("eta_"+postfix)->Fill(vec.Eta(), weight);
  Hist("phi_"+postfix)->Fill(vec.Phi(), weight);
  Hist("mass_"+postfix)->Fill(vec.mass(), weight);
  Hist("dR_TopLep_"+postfix)->Fill(deltaR(m_tt->TopLep(),vec ),weight);
  Hist("dPhi_TopLep_"+postfix)->Fill(deltaPhiAbs(m_tt->TopLep().Phi(), vec.Phi()),weight); 
}

void WHists::FillPtEtaPhi(Particle particle, TString postfix, double weight){
  Hist("pt_"+postfix)->Fill(particle.pt(), weight);
  Hist("eta_"+postfix)->Fill(particle.eta(), weight);
  Hist("phi_"+postfix)->Fill(particle.phi(), weight);
  Hist("dR_TopLep_"+postfix)->Fill(deltaR(m_tt->TopLep(),particle.v4() ),weight);
  Hist("dPhi_TopLep_"+postfix)->Fill(deltaPhiAbs(m_tt->TopLep().Phi(), particle.phi()),weight); 
}






void WHists::Finish()
{
  // final calculations, like division and addition of certain histograms
  // EventCalc* calc = EventCalc::Instance();
  // bool IsRealData = calc->IsRealData();
 
}
