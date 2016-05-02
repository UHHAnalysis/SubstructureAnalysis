#include "include/TruthHists.h"
#include "NtupleWriter/include/JetProps.h"
#include "SFrameTools/include/EventCalc.h"
#include "NtupleWriter/include/GenTopJet.h"

#include "TH2F.h"
#include <iostream>


using namespace std;

TruthHists::TruthHists(const char* name) : BaseHists(name)
{
  // named default constructor
   
}

TruthHists::~TruthHists(){
}

void TruthHists::Init()
{
  // book all histograms here
  //  Book( TH1F("jetmass","jetmass",30,0,450));
  Book( TH1F("dr_hadr_W_hadr_b", "#DeltaR(hadronic W, hadronic b)", 20,0,10));
  Book( TH1F("dr_hadr_W_lept_b", "#DeltaR(hadronic W, leptonic b)", 20,0,10));
  Book( TH1F("dr_hadr_W_lepton", "#DeltaR(hadronic W, lepton)", 20,0,10));
  Book( TH1F("dr_lept_b_hadr_b", "#DeltaR(leptonic b, hadronic b)", 20,0,10));
  Book( TH1F("dr_lept_b_lepton", "#DeltaR(leptonic b, lepton)", 20,0,10));

  Book( TH1F("dphi_hadr_W_hadr_b", "#Delta#phi(hadronic W, hadronic b)", 20,0,10));
  Book( TH1F("dphi_hadr_W_lept_b", "#Delta#phi(hadronic W, leptonic b)", 20,0,10));
  Book( TH1F("dphi_hadr_W_lepton", "#Delta#phi(hadronic W, lepton)", 20,0,10));
  Book( TH1F("dphi_lept_b_hadr_b", "#Delta#phi(leptonic b, hadronic b)", 20,0,10));
  Book( TH1F("dphi_lept_b_lepton", "#Delta#phi(leptonic b, lepton)", 20,0,10));
 
  Book( TH1F("pt_hadr_b", "p_{t, hadronic b}",50 ,0,500));
  Book( TH1F("pt_leponic_b", "p_{t, leptonic b}",50 ,0,500));
}

void TruthHists::Fill()
{
  // fill the histograms


  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  bool IsRealData = calc->IsRealData();
  LuminosityHandler* lumih = calc->GetLumiHandler();

  // important: get the event weight
  double weight = calc->GetWeight();

  int run = calc->GetRunNum();
  int lumiblock = calc->GetLumiBlock();
  int Npvs = calc->GetPrimaryVertices()->size();

 

  //   Hist("pt_leading")->Fill(Top_Jets->at(0).pt(),weight);





  Bool_t isMuon = false;
  Bool_t isEle = false;

  const GenParticle* electron = NULL;
  const GenParticle* muon = NULL;
  const GenParticle* lepton = NULL;
  const GenParticle* neutrino = NULL;
  const GenParticle* hadr_W = NULL;
  const GenParticle* hadr_b = NULL;
  const GenParticle* lept_b = NULL;

  GenParticle hadr_top;

  if ( bcc->genparticles){
    std::vector<GenParticle>* genparts = bcc->genparticles;
    int Ngenparts = genparts->size();
 
    for(int l=0; l<Ngenparts ; l++){
      GenParticle genpart = genparts->at(l);
      if(genpart.status()<3){break;}
      if( fabs(genpart.pdgId()) == 6){
	//     if((fabs(genpart.daughter(genparts, 1)->pdgId())!=24)) daughter1 = genpart.daughter(genparts, 2);//W-bosons
	const  GenParticle* W; 
	const  GenParticle* b; 
	const  GenParticle* daughter1 = genpart.daughter(genparts, 1);
	const  GenParticle* daughter2 = genpart.daughter(genparts, 2);
	if((fabs(genpart.daughter(genparts, 1)->pdgId())!=24)){
	  W = daughter2;
	  b = daughter1;
	}
	else{
	  W = daughter1;
	  b = daughter2; 
	} 
	if(W){
	  // if( fabs(daughter->pdgId())<6) hadr_top=genpart;
	  if( fabs(W->daughter(genparts, 1)->pdgId())<6){
	    hadr_top=genpart;
	    hadr_W = W;
	    hadr_b = b;  
	  }
	  else lept_b = b; 
	  //	 if( fabs(daughter->pdgId())>12 && fabs(daughter->pdgId())<15){
	  if( fabs(W->daughter(genparts, 1)->pdgId()) == 13){
	    isMuon=true;
	    muon = W->daughter(genparts, 1);
	    neutrino = W->daughter(genparts, 2);
	  }
	  else if( fabs(W->daughter(genparts, 2)->pdgId()) == 13){
	    isMuon=true;
	    muon = W->daughter(genparts, 2);
	    neutrino = W->daughter(genparts, 1);
	  }
	  else if( fabs(W->daughter(genparts, 1)->pdgId()) == 11){
	    isEle=true;
	    electron = W->daughter(genparts, 1);
	    neutrino = W->daughter(genparts, 2);
	  }
	  else if( fabs(W->daughter(genparts, 2)->pdgId()) == 11){
	    isEle=true;
	    electron = W->daughter(genparts, 2);
	    neutrino = W->daughter(genparts, 1);
	  }
	}
      }
    }

  }

  if(isMuon) lepton = muon;
  else if(isEle) lepton = electron; 




  Hist("dr_hadr_W_hadr_b")->Fill(deltaR(hadr_W->v4(), hadr_b->v4()), weight);
  Hist("dr_hadr_W_lept_b")->Fill(deltaR(hadr_W->v4(), lept_b->v4()), weight);
  Hist("dr_hadr_W_lepton")->Fill(deltaR(hadr_W->v4(), lepton->v4()), weight);
  Hist("dr_lept_b_hadr_b")->Fill(deltaR(lept_b->v4(), hadr_b->v4()), weight);
  Hist("dr_lept_b_lepton")->Fill(deltaR(lept_b->v4(), lepton->v4()), weight);


  Hist("dphi_hadr_W_hadr_b")->Fill(deltaPhiAbs(hadr_W->phi(), hadr_b->phi()), weight);
  Hist("dphi_hadr_W_lept_b")->Fill(deltaPhiAbs(hadr_W->phi(), lept_b->phi()), weight);
  Hist("dphi_hadr_W_lepton")->Fill(deltaPhiAbs(hadr_W->phi(), lepton->phi()), weight);
  Hist("dphi_lept_b_hadr_b")->Fill(deltaPhiAbs(lept_b->phi(), hadr_b->phi()), weight);
  Hist("dphi_lept_b_lepton")->Fill(deltaPhiAbs(lept_b->phi(), lepton->phi()), weight);
 
  Hist("pt_hadr_b")->Fill(hadr_b->pt(), weight);
  Hist("pt_leponic_b")->Fill(lept_b->pt(), weight);

 
  

}


 


void TruthHists::Finish()
{
  /*
  // final calculations, like division and addition of certain histograms
  */
}
