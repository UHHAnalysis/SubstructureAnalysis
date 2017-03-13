#include "include/ComparisonHistsReco.h"
#include "NtupleWriter/include/JetProps.h"
#include "SFrameTools/include/EventCalc.h"
#include "TH2F.h"
//#include "vector.h"
#include "../include/JetMassTools.h"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "NtupleWriter/include/JetProps.h"
#include "TLorentzVector.h"
#include <iostream>


using namespace std;

ComparisonHistsReco::ComparisonHistsReco(const char* name) : BaseHists(name)
{
  // named default constructor
   
}

ComparisonHistsReco::~ComparisonHistsReco(){
}

void ComparisonHistsReco::Init()
{
  Book( TH1F("r","",50,0,500));
  //masses
  Book( TH1F("reco_trimmed","m_{jet, trimmed} [GeV]",50,0,500));
  Book( TH1F("2nd_jet_trimmed","m_{2nd jet, trimmed} [GeV]",50,0,500));
  Book( TH1F("reco","m_{jet} [GeV]",50,0,500));
  Book( TH1F("2nd_jet_mass","m_{2nd jet} [GeV]",50,0,500));
  Book( TH1F("3rd_jet_mass","m_{3rd jet} [GeV]",50,0,500));
  Book( TH1F("4th_jet_mass","m_{4th jet} [GeV]",50,0,500));
  Book( TH1F("5th_jet_mass","m_{5th jet} [GeV]",50,0,500));
  Book( TH1F("mass_diff","m_{leading jet}-m_{2nd jet} [GeV]",100,-500,500));
  Book( TH1F("mtt", "m_{t#bar{t}} (NoNu) [GeV]", 100, 0, 2000));

  //momenta
  Book( TH1F("reco_pt1_trimmed","p_{t, leading jet trimmed} [GeV]",50,0,1000));
  Book( TH1F("reco_pt2_trimmed","p_{t, 2nd jet trimmed} [GeV]",50,0,500));
  Book( TH1F("reco_pt1","p_{t, leading jet} [GeV]",50,0,1000));
  Book( TH1F("reco_pt2","p_{t, 2nd jet} [GeV]",50,0,500));
  Book( TH1F("reco_pt3","p_{t, 3rd jet} [GeV]",50,0,500));
  Book( TH1F("reco_pt4","p_{t, 4th jet} [GeV]",50,0,500));
  Book( TH1F("reco_pt5","p_{t, 5th jet} [GeV]",50,0,500));
  Book( TH1F("pt_AddAK5","p_{t, additional ak5 jet} [GeV]",50,0,500));

  Book( TH1F("lepton_pt","p_{t,lepton} [GeV]",50,0,500));

  //distances
  Book( TH1F("drlep","#DeltaR_{leading jet,lepton}",20,0,10));
  Book( TH1F("drlep2","#DeltaR_{2nd jet,lepton}",20,0,10));
  Book( TH1F("dr_diff","#DeltaR_{leading jet,lepton}-#DeltaR_{2nd jet,lepton}",20,-10,10));

  Book( TH1F("dr_jets","#DeltaR_{leading jet, 2nd jet}",20,0,10));

  Book( TH1F("dr","#DeltaR_{jet,top}",20,0,10));
  Book( TH1F("dr_recgen","#DeltaR_{reco jet, gen jet}",20,0,10));
  Book( TH1F("drnewlep","#DeltaR_{jet,lepton}",20,0,10));


  

  //agles
  Book( TH1F("eta1","#eta_{leading jet}",60,-3,3));
  Book( TH1F("eta2","#eta_{2nd jet}",60,-3,3));
  Book( TH1F("eta3","#eta_{3rd jet}",120,-6,6));
  Book( TH1F("eta4","#eta_{4th jet}",120,-6,6));
  Book( TH1F("eta5","#eta_{5th jet}",120,-6,6));
  Book( TH1F("eta_lep","#eta_{lepton}",60,-3,3));

  Book( TH1F("phi1","#phi_{leading jet}",40,-3.5,3.5));
  Book( TH1F("phi2","#phi_{2nd jet}",40,-3.5,3.5));
  Book( TH1F("phi3","#phi_{3rd jet}",40,-3.5,3.5));
  Book( TH1F("phi4","#phi_{4th jet}",40,-3.5,3.5));
  Book( TH1F("phi5","#phi_{5th jet}",40,-3.5,3.5));
  Book( TH1F("phi_lep","#phi_{lepton}",40,-3.5,3.5));

  //numbers
  Book( TH1F("Njets","N_{jets}",10,0,10));
  Book( TH1F("Nsubjets","N_{subjets}",10,0,10));
  Book( TH1F("Nsubjets_2","N_{subjets, 2nd jet}",10,0,10));
  Book( TH1F("NCAjets","N_{C/A jets} (p_{t} > 150 GeV, #eta < 2.5)",10,0,10));

  Book( TH1F("mass/NPFC", "m_{jet}/Number of PF candidates" ,50, 0, 10)); 



  //Scatter Plots (versus gen masses)
  Book( TH2F("mass_vs_LeadingJetPt","mass_vs_LeadingJetPt",50,0,500,70,0,700));
  Book( TH2F("mass_vs_2ndJetPt","mass_vs_2ndJetPt",50,0,500,70,0,700));
  Book( TH2F("mass_vs_3rdJetPt","mass_vs_3rdJetPt",50,0,500,70,0,700));
  Book( TH2F("mass_vs_TopPt","mass_vs_TopPt",50,0,500,70,0,700));
  Book( TH2F("mass_vs_LeptonPt","mass_vs_LeptonPt",50,0,500,70,0,700));
  Book( TH2F("mass_vs_Nsubjets","mass_vs_Nsubjets",50,0,500,30,0,6));
  

  
  Book( TH1F("Wjet_mass","m_{jet, Wjet}",50,0,500));
  Book( TH1F("Wjet_mass_gen","m_{jet, Wjet}",50,0,500));
  Book( TH1F("Wjet_mass_matched","m_{jet, Wjet}",50,0,500));
  Book( TH1F("Wjet_mass_gen_matched","m_{jet, Wjet}",50,0,500));
  Book( TH2F("Wjet_mass_vs_Wjet_pt","m_{jet, Wjet}",50,0,500,50,0,500));

}

void ComparisonHistsReco::Fill()
{
  // fill the histograms


  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  bool IsRealData = calc->IsRealData();
  LuminosityHandler* lumih = calc->GetLumiHandler();
  Particle Lepton = *calc->GetPrimaryLepton();

  // important: get the event weight
  double weight = calc->GetWeight();

  int run = calc->GetRunNum();
  int lumiblock = calc->GetLumiBlock();
  int Npvs = calc->GetPrimaryVertices()->size();


 
  bool geninfo = true;
  if(IsRealData) geninfo = false; 

  GenParticle hadr_top;

  std::vector<TopJet> topjets = *calc->GetCAJets();

  std::vector<Jet> ak5jets = *calc->GetJets();
 

  if(topjets.size()>1){

    double dr_jet_lep = -1;
    double dr_2ndjet_lep = -1;
    double dr_diff;  
    dr_jet_lep = deltaR(topjets.at(0).v4(), Lepton.v4());
    dr_2ndjet_lep = deltaR(topjets.at(1).v4(), Lepton.v4());  
 
    dr_diff =  dr_jet_lep-dr_2ndjet_lep;
 

    //masses
    Hist("reco")->Fill( topjets.at(0).v4().mass(),weight);
    if(topjets.size()>1) Hist("2nd_jet_mass")->Fill(topjets.at(1).v4().mass(),weight); 
    if(topjets.size()>2) Hist("3rd_jet_mass")->Fill(topjets.at(2).v4().mass(),weight);
    if(topjets.size()>3) Hist("4th_jet_mass")->Fill(topjets.at(3).v4().mass(),weight);
    if(topjets.size()>4) Hist("5th_jet_mass")->Fill(topjets.at(4).v4().mass(),weight);
    Hist("mass_diff")->Fill(topjets.at(0).v4().mass()-topjets.at(1).v4().mass(),weight);

    std::vector<LorentzVector> nus = calc->NeutrinoReconstruction(Lepton.v4(), calc->GetMET()->v4());
    LorentzVector nu = nus.at(0);
    if(nus.size()>1){
      if( fabs(topjets.at(0).v4().mass()-(topjets.at(1).v4()+Lepton.v4()+nus.at(0)).mass()) > fabs(topjets.at(0).v4().mass()-(topjets.at(1).v4()+Lepton.v4()+nus.at(1)).mass()))
	nu = nus.at(1);
    }
    Hist("mtt")->Fill((topjets.at(0).v4()+(topjets.at(1).v4()+Lepton.v4()+nu)).mass(),weight);



    //momenta
    Hist("reco_pt1")->Fill(topjets.at(0).pt(),weight);
    if(topjets.size()>1) Hist("reco_pt2")->Fill(topjets.at(1).pt(),weight);
    if(topjets.size()>2)  Hist("reco_pt3")->Fill(topjets.at(2).pt(),weight);
    if(topjets.size()>3)  Hist("reco_pt4")->Fill(topjets.at(3).pt(),weight);
    if(topjets.size()>4)  Hist("reco_pt5")->Fill(topjets.at(4).pt(),weight);
    Hist("lepton_pt")->Fill(Lepton.pt(),weight);

 

    for(unsigned int i = 0; i < ak5jets.size(); i++){
      if(topjets.size()>1 && deltaR(ak5jets.at(i).v4(), topjets.at(0).v4())>1.7 && deltaR(ak5jets.at(i).v4(), topjets.at(1).v4())>1.7 ){
	Hist("pt_AddAK5")->Fill(ak5jets.at(i).pt(),weight);
	break;
      }
    }
 
    //distances
    //  Hist("dr")->Fill(deltaR(topjets.at(0).v4(), hadr_top.v4()),weight);
    Hist("dr_jets")->Fill(deltaR(topjets.at(0).v4(),topjets.at(1).v4()),weight);
    Hist("drlep")->Fill(dr_jet_lep,weight);
    Hist("drlep2")->Fill(dr_2ndjet_lep,weight);
    Hist("dr_diff")->Fill(dr_diff,weight);
 
    //angles
    Hist("eta1")->Fill(topjets.at(0).eta(),weight);
    if(topjets.size()>1) Hist("eta2")->Fill(topjets.at(1).eta(),weight);
    if(topjets.size()>2)  Hist("eta3")->Fill(topjets.at(2).eta(),weight);
    if(topjets.size()>3)  Hist("eta4")->Fill(topjets.at(3).eta(),weight);
    if(topjets.size()>4)  Hist("eta5")->Fill(topjets.at(4).eta(),weight);
    Hist("eta_lep")->Fill(Lepton.eta(),weight);

    Hist("phi1")->Fill(topjets.at(0).phi(),weight);
    if(topjets.size()>1) Hist("phi2")->Fill(topjets.at(1).phi(),weight);
    if(topjets.size()>2)  Hist("phi3")->Fill(topjets.at(2).phi(),weight);
    if(topjets.size()>3)  Hist("phi4")->Fill(topjets.at(3).phi(),weight);
    if(topjets.size()>4)  Hist("phi5")->Fill(topjets.at(4).phi(),weight);
    Hist("phi_lep")->Fill(Lepton.phi(),weight);


    //numbers
    // Hist("Njets")->Fill(,weight);
    Hist("NCAjets")->Fill(topjets.size(),weight);
    Hist("Nsubjets")->Fill(topjets.at(0).subjets().size(),weight);
    Hist("Nsubjets_2")->Fill(topjets.at(1).subjets().size(),weight);
    Hist("mass/NPFC")->Fill(topjets.at(0).v4().mass()/topjets.at(0).pfconstituents_indices().size(), weight);

 
    //scatter plots
    if(geninfo && bcc->cagenjets){
      if(bcc->cagenjets->size()>1){
	bool gen_newsel1 = false;
	bool gen_newsel2 = false;
	bool gen_newsel1_400 = false;
	bool gen_newsel2_400 = false;
      

	if(bcc->cagenjets->at(0).v4().mass() >  bcc->cagenjets->at(1).v4().mass()){ 
	  gen_newsel1 = true;
	  if(bcc->cagenjets->at(0).pt()>400){ 
	    gen_newsel1_400 = true;
	  }
	}
	if(bcc->cagenjets->at(0).v4().mass() <  bcc->cagenjets->at(1).v4().mass()){
	  gen_newsel2 = true;
	  if(bcc->cagenjets->at(1).pt()>400){
	    gen_newsel2_400 = true;
	  }
	}
  
    
	// for(int i = 0; i<2 ; i++){
	bool reco_sel1;
	bool reco_sel2;
	bool gen_sel1;
	bool gen_sel2;
	TString mode;

	mode ="";
	TopJet topjet;
	Particle genjet;
	topjet = topjets.at(0);

	genjet = bcc->cagenjets->at(0);
	  
	Hist("dr_recgen")->Fill(deltaR(topjet.v4(), genjet.v4()), weight);
	Hist("drnewlep")->Fill(deltaR(Lepton.v4(), topjet.v4()),weight);

	((TH2F*)Hist("mass_vs_LeadingJetPt"+mode))->Fill(genjet.v4().mass(),topjets.at(0).pt(),weight);
	if(topjets.size()>1) ((TH2F*)Hist("mass_vs_2ndJetPt"+mode))->Fill(genjet.v4().mass(),topjets.at(1).pt(),weight);
	if(topjets.size()>2) ((TH2F*)Hist("mass_vs_3rdJetPt"+mode))->Fill(genjet.v4().mass(),topjets.at(2).pt(),weight);

	((TH2F*)Hist("mass_vs_TopPt"+mode))->Fill(genjet.v4().mass(),hadr_top.pt(),weight);

	((TH2F*)Hist("mass_vs_LeptonPt"+mode))->Fill(genjet.v4().mass(),Lepton.pt(),weight);
	((TH2F*)Hist("mass_vs_Nsubjets"+mode))->Fill(genjet.v4().mass(),topjet.subjets().size(),weight); 
      }
    }
  }
  
 
 
}

void ComparisonHistsReco::Finish()
{
  // final calculations, like division and addition of certain histograms
  /* EventCalc* calc = EventCalc::Instance();
     bool IsRealData = calc->IsRealData();
     if (IsRealData){
    
     }
  */
}

