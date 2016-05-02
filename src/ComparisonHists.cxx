#include "include/ComparisonHists.h"
#include "NtupleWriter/include/JetProps.h"
#include "SFrameTools/include/EventCalc.h"
#include "NtupleWriter/include/GenTopJet.h"
#include "TH2F.h"
#include "SFrameTools/include/TTbarGen.h"
#include <iostream>


using namespace std;

ComparisonHists::ComparisonHists(const char* name) : BaseHists(name)
{
  // named default constructor
   
}

ComparisonHists::~ComparisonHists(){
}

void ComparisonHists::Init()
{

  //  TFile* f = new TFile("~/Substructure_makros/ScaleFactors2.root","READ"); 
  // m_scale_factors = (TH1F*)f->Get("scale_factors_pt_nosel");

  Book( TH1F("g","",50,0,500));
  //masses
  Book( TH1F("gen","m_{jet}",50,0,500));
  Book( TH1F("2nd_jet_mass","m_{2nd jet}",50,0,500));

  Book( TH1F("gen_WP","m_{jet} (genjetswithparts)",50,0,500));
  Book( TH1F("2nd_jet_mass_WP","m_{2nd jet} (genjetswithparts)",50,0,500));

  Book( TH1F("3rd_jet_mass","m_{3rd jet}",50,0,500));
  Book( TH1F("4th_jet_mass","m_{4th jet}",50,0,500));
  Book( TH1F("5th_jet_mass","m_{5th jet}",50,0,500));
  Book( TH1F("mass_diff","m_{leading jet}-m_{2nd jet}",100,-500,500));
  Book( TH1F("leading_plus_nearest_mass","m_{leading jet + nearest jet pt < 150} (#Delta#phi < 1.5)",50,0,500));
  Book( TH1F("mtt", "m_{t#bar{t}} [GeV]", 100, 0, 2000));

  Book( TH1F("gen_boosted","m_{jet}",50,0,500));

  //momenta
  Book( TH1F("gen_pt1","p_{t, leading jet}",50,0,1000));
  Book( TH1F("gen_pt2","p_{t, 2nd jet}",50,0,500));
  Book( TH1F("gen_pt3","p_{t, 3rd jet}",50,0,500));

  Book( TH1F("gen_pt1_WP","p_{t, leading jet} (genjetswithparts)",50,0,1000));
  Book( TH1F("gen_pt2_WP","p_{t, 2nd jet} (genjetswithparts)",50,0,500));
  Book( TH1F("gen_pt3_WP","p_{t, 3rd jet} (genjetswithparts)",50,0,500));

  Book( TH1F("gen_pt4","p_{t, 4th jet}",50,0,500));
  Book( TH1F("gen_pt5","p_{t, 5th jet}",50,0,500));
 

  Book( TH1F("lepton_pt","p_{t,lepton}",50,0,500));
  Book( TH1F("top_pt","p_{t,top}",80,0,800));
 
  Book( TH1F("leading_plus_nearest_pt","p_{t, leading jet + nearest jet pt < 150} (#Delta#phi < 1.5)",50,0,500));
  Book( TH1F("leading_plus_nearest_pt_diff","p_{t, leading jet + nearest jet pt < 150} - p_{t,2} (#Delta#phi < 1.5)",50,0,500));
  //distances
  Book( TH1F("drlep","#DeltaR_{leading jet,lepton}",20,0,10));
  Book( TH1F("drlep2","#DeltaR_{2nd jet,lepton}",20,0,10));
  Book( TH1F("dr_diff","#DeltaR_{leading jet,lepton}-#DeltaR_{2nd jet,lepton}",20,-10,10));

  Book( TH1F("dr_jets","#DeltaR_{leading jet, 2nd jet}",20,0,10));

  Book( TH1F("dr","#DeltaR_{jet, top}",20,0,10));
  Book( TH1F("leading_plus_nearest_dr","#DeltaR(leading jet , nearest jet pt < 150}",20,0,10));

  Book( TH1F("dPhi_2ndjet_lepton", "#Delta#phi_{2nd jet, lepton}", 16, 0, 3.2)); 

  //boost cuts
  Book( TH2F("dPhi_2ndjet_lepton_vs_dR_2ndjet_lepton", "#Delta#phi_{2nd jet, lepton} vs #Delta#R_{2nd jet, lepton}",20,0,5, 16, 0, 3.2));
  Book( TH2F("pTopHad_vs_dR_2ndjet_lepton", "#p_{top_had} vs #Delta#R_{2nd jet, lepton}",20,0,5,100,0,1000));
  Book( TH2F("dR_W_b_vs_dR_2ndjet_lepton", "#Delta#phi_{W had, b} vs #Delta#R_{2nd jet, lepton}",20,0,5,16, 0, 3.2));


  //agles
  Book( TH1F("eta1","#eta_{leading jet}",30,-3,3));
  Book( TH1F("eta2","#eta_{2nd jet}",30,-3,3));
  Book( TH1F("eta3","#eta_{3rd jet}",60,-6,6));

  Book( TH1F("eta1_WP","#eta_{leading jet} (genjetswithparts)",30,-3,3));
  Book( TH1F("eta2_WP","#eta_{2nd jet} (genjetswithparts)",30,-3,3));
  Book( TH1F("eta3_WP","#eta_{3rd jet} (genjetswithparts)",60,-6,6));

  Book( TH1F("eta4","#eta_{4th jet}",60,-6,6));
  Book( TH1F("eta5","#eta_{5th jet}",60,-6,6));
  Book( TH1F("eta_lep","#eta_{lepton}",30,-3,3));

  Book( TH1F("phi1","#phi_{leading jet}",40,-3.5,3.5));
  Book( TH1F("phi2","#phi_{2nd jet}",40,-3.5,3.5));
  Book( TH1F("phi3","#phi_{3rd jet}",40,-3.5,3.5));

  Book( TH1F("phi1_WP","#phi_{leading jet} (genjetswithparts)",40,-3.5,3.5));
  Book( TH1F("phi2_WP","#phi_{2nd jet} (genjetswithparts)",40,-3.5,3.5));
  Book( TH1F("phi3_WP","#phi_{3rd jet} (genjetswithparts)",40,-3.5,3.5));

  Book( TH1F("phi4","#phi_{4th jet}",40,-3.5,3.5));
  Book( TH1F("phi5","#phi_{5th jet}",40,-3.5,3.5));
  Book( TH1F("phi_lep","#phi_{lepton}",40,-3.5,3.5));
  Book( TH1F("leading_plus_nearest_dphi","#Delta#phi(leading jet , nearest jet pt < 150}",40,-3.5,3.5));
  //numbers


  Book( TH1F("hadrlept","hadronic or leptonic",6,0,6));
  Book( TH1F("NCAjets","N_{C/A jets} (p_{t} > 150 GeV, #eta < 2.5)",10,0,10));
  Book( TH1F("N_addCAjets","N_{additional C/A jets, p_{t} < 150, #eta < 2.5}",10,0,10));



  //Scatter Plots
  Book( TH2F("mass_vs_LeadingJetPt","mass_vs_LeadingJetPt",50,0,500,70,0,700));
  Book( TH2F("mass_vs_2ndJetPt","mass_vs_2ndJetPt",50,0,500,70,0,700));
  Book( TH2F("mass_vs_3rdJetPt","mass_vs_3rdJetPt",50,0,500,70,0,700));
  Book( TH2F("mass_vs_TopPt","mass_vs_TopPt",50,0,500,70,0,700));
  Book( TH2F("mass_vs_LeptonPt","mass_vs_LeptonPt",50,0,500,70,0,700));
  Book( TH2F("Pt1/Pt2_vs_TopPt","Pt1/Pt2_vs_TopPt",80,-4,4,80,0,800));
  Book( TH2F("mass_vs_Pt1/Pt2","mass_vs_Pt1/Pt2",80,-4,4,80,0,800));

  //Substructure
  Book( TH1F("Nparts","N_{gen particles} (leading jet)",125,0,350));
  Book( TH1F("Ptparts","p_{T, gen particle} (leading jet)",50,0,50));
  Book( TH1F("DRparts","#Delta R(particle, jet axis) (leading jet)",15,0,1.5));
}

void ComparisonHists::Fill()
{
  // fill the histograms


  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  bool IsRealData = calc->IsRealData();
  LuminosityHandler* lumih = calc->GetLumiHandler();

  // important: get the event weight
  double weight = calc->GetGenWeight();

  int run = calc->GetRunNum();
  int lumiblock = calc->GetLumiBlock();
  int Npvs = calc->GetPrimaryVertices()->size();

 
 
  // Bool_t isMuon = false;
  // Bool_t isEle = false;
  // const GenParticle* electron = NULL;
  // const GenParticle* muon = NULL;

  int hadr = 0;
  int lept = 0;

  GenParticle hadr_top;
  GenParticle lept_top;

  GenParticle lepton;


  if ( bcc->genparticles){
 
    TTbarGen *ttbargen = new TTbarGen(bcc);

   std::vector<GenParticle>* genparts = bcc->genparticles;
    //int Ngenparts = genparts->size();
    if(ttbargen->DecayChannel() == TTbarGen::e_ehad || ttbargen->DecayChannel() == TTbarGen::e_muhad  ){
      hadr_top = ttbargen->TopHad();
      lept_top = ttbargen->TopLep();
      lepton = ttbargen->ChargedLepton();
    }

    double dr_jet_lep = -1;
    double dr_2ndjet_lep = -1;
    double dr_diff;
    double dr_lep_hadr1;
    double dr_lep_hadr2;

 
    dr_jet_lep = deltaR(bcc->cagenjets->at(0).v4(), lepton.v4());
    dr_2ndjet_lep = deltaR(bcc->cagenjets->at(1).v4(), lepton.v4());
    dr_diff =  dr_jet_lep-dr_2ndjet_lep;

    /*

      for(int l=0; l<Ngenparts ; l++){
      GenParticle genpart = genparts->at(l);
      if(genpart.status()<3){break;}
      if( fabs(genpart.pdgId()) == 6){

      const  GenParticle* daughter = genpart.daughter(genparts, 1);
      if((fabs(genpart.daughter(genparts, 1)->pdgId())!=24)) daughter = genpart.daughter(genparts, 2);//W-bosons
      if(daughter){
      if( fabs(daughter->daughter(genparts, 1)->pdgId())<6) {hadr_top=genpart; hadr ++;}
      if( fabs(daughter->daughter(genparts, 1)->pdgId())>6) {lept_top=genpart; lept ++;}
      // cout << daughter->pdgId() << " -> " << daughter->daughter(genparts, 1)->pdgId() << " + " << daughter->daughter(genparts, 2)->pdgId() << endl;
      // cout << daughter->daughter(genparts, 1)->mother(genparts, 1)->pdgId() << "  " << daughter->daughter(genparts, 1)->mother(genparts, 1)->pdgId() << endl;
      // cout << ttbargen->WTop().daughter(genparts, 1)->pdgId() << " + " << ttbargen->WTop().daughter(genparts, 2)->pdgId() << endl;
      // cout << ttbargen->WAntitop().daughter(genparts, 1)->pdgId() << " + " << ttbargen->WAntitop().daughter(genparts, 2)->pdgId() << endl;
      if( fabs(daughter->daughter(genparts, 1)->pdgId()) == 13){
      isMuon=true;
      muon = daughter->daughter(genparts, 1);
      }
      else if( fabs(daughter->daughter(genparts, 2)->pdgId()) == 13){
      isMuon=true;
      muon = daughter->daughter(genparts, 2);
      }
      else if( fabs(daughter->daughter(genparts, 1)->pdgId()) == 11){
      isEle=true;
      electron = daughter->daughter(genparts, 1);
      }
      else if( fabs(daughter->daughter(genparts, 2)->pdgId()) == 11){
      isEle=true;
      electron = daughter->daughter(genparts, 2);
      }
      }
      }
      }

      }

      Bool_t isHadr=false;
      if(bcc->cagenjets->size()>1 && deltaR(bcc->cagenjets->at(0).v4(), hadr_top.v4())<1.2){
      isHadr= true;
      }

  
      weight = calc->GetGenWeight();


      if(bcc->cagenjets->size()>1){
      Hist("top_pt")->Fill(hadr_top.pt(),weight);
      // LorentzVector lep_top1;
      // LorentzVector lep_top2;
      LorentzVector lepton;

      double dr_jet_lep = -1;
      double dr_2ndjet_lep = -1;
      double dr_diff;
      double dr_lep_hadr1;
      double dr_lep_hadr2;

      if(isMuon){
      dr_jet_lep = deltaR(bcc->cagenjets->at(0).v4(), muon->v4());
      dr_2ndjet_lep = deltaR(bcc->cagenjets->at(1).v4(), muon->v4());
 
      lepton= muon->v4();

      }else{
      if(isEle){
      dr_jet_lep = deltaR(bcc->cagenjets->at(0).v4(), electron->v4());
      dr_2ndjet_lep = deltaR(bcc->cagenjets->at(1).v4(), electron->v4());
  
      lepton= electron->v4();
      }     
      }
      dr_diff =  dr_jet_lep-dr_2ndjet_lep;
  
    */
    LorentzVector leading_plus = bcc->cagenjets->at(0).v4();
    LorentzVector nearest; 
    double pt_max = 0;
    double dR_leading_nearest = 10;
    int N_add_jets = 0;
    if(bcc->cagenjets->size()>2){
      for(unsigned int i = 2; i< bcc->cagenjets->size() ; i++){
	if(bcc->cagenjets->at(i).pt()<150 
	
	   && bcc->cagenjets->at(i).eta()<2.5){
	  N_add_jets++;

	  if(deltaR(bcc->cagenjets->at(0).v4(),bcc->cagenjets->at(i).v4())< dR_leading_nearest){
	    dR_leading_nearest = deltaR(bcc->cagenjets->at(0).v4(),bcc->cagenjets->at(i).v4());
	 
	    nearest = bcc->cagenjets->at(i).v4();
	    leading_plus = bcc->cagenjets->at(0).v4()+nearest;
	    //	  }
	  }
	  if(bcc->cagenjets->at(i).pt()>pt_max){
	    pt_max = bcc->cagenjets->at(i).pt();

	  }
	}
      }
    }


 

    bool boosted = false;
    bool very_boosted = false;

  
    if( deltaR(bcc->cagenjets->at(1).v4(), ttbargen->bTop().v4()) < 1.2 
	&& deltaR(bcc->cagenjets->at(1).v4(), ttbargen->Wdecay1().v4()) < 1.2 
	&& deltaR(bcc->cagenjets->at(1).v4(), ttbargen->Wdecay2().v4()) < 1.2 
	){
      boosted = true;
    }else{
      if( deltaR(bcc->cagenjets->at(1).v4(), ttbargen->bAntitop().v4()) < 1.2 
	  && deltaR(bcc->cagenjets->at(1).v4(), ttbargen->WMinusdecay1().v4()) < 1.2 
	  && deltaR(bcc->cagenjets->at(1).v4(), ttbargen->WMinusdecay2().v4()) < 1.2 
	  ){
	boosted = true;
      }
    }

    if( deltaR(bcc->cagenjets->at(1).v4(), ttbargen->bTop().v4()) < 1.2 
	&& deltaR(ttbargen->bTop().v4(), ttbargen->Wdecay1().v4()) < 1.2 
	&& deltaR(ttbargen->bTop().v4(), ttbargen->Wdecay2().v4()) < 1.2  
	){
      very_boosted = true;
    }else{
      if( deltaR(bcc->cagenjets->at(1).v4(), ttbargen->bAntitop().v4()) <  1.2 
	  && deltaR(ttbargen->bAntitop().v4(), ttbargen->WMinusdecay1().v4()) <  1.2  
	  && deltaR(ttbargen->bAntitop().v4(), ttbargen->WMinusdecay2().v4()) <  1.2 
	  ){
	very_boosted = true;
      }
    }



    //masses
    Hist("gen")->Fill( bcc->cagenjets->at(0).v4().mass(),weight);

    if(boosted)Hist("gen_boosted")->Fill( bcc->cagenjets->at(0).v4().mass(),weight);
    if(very_boosted)Hist("g")->Fill( bcc->cagenjets->at(0).v4().mass(),weight);

    if(bcc->cagenjets->size()>1) Hist("2nd_jet_mass")->Fill(bcc->cagenjets->at(1).v4().mass(),weight); 
    if(bcc->cagenjets->size()>2) Hist("3rd_jet_mass")->Fill(bcc->cagenjets->at(2).v4().mass(),weight);
    if(bcc->cagenjets->size()>3) Hist("4th_jet_mass")->Fill(bcc->cagenjets->at(3).v4().mass(),weight);
    if(bcc->cagenjets->size()>4) Hist("5th_jet_mass")->Fill(bcc->cagenjets->at(4).v4().mass(),weight);

    Hist("gen_WP")->Fill( bcc->genjetswithparts->at(0).v4().mass(),weight);
    if(bcc->genjetswithparts->size()>1) Hist("2nd_jet_mass_WP")->Fill(bcc->genjetswithparts->at(1).v4().mass(),weight); 

    Hist("mass_diff")->Fill(bcc->cagenjets->at(0).v4().mass()-bcc->cagenjets->at(1).v4().mass(),weight);
    Hist("mtt")->Fill((hadr_top.v4()+lept_top.v4()).M(), weight);

    Hist("dr")->Fill(deltaR(bcc->cagenjets->at(0).v4(), hadr_top.v4()),weight);
    Hist("drlep")->Fill(dr_jet_lep,weight);
    Hist("drlep2")->Fill(dr_2ndjet_lep,weight);
    Hist("dr_diff")->Fill(dr_diff,weight);

    if(deltaPhiAbs(bcc->cagenjets->at(0).phi(), nearest.Phi())<1.5) Hist("leading_plus_nearest_mass")->Fill(leading_plus.M(),weight);
    else Hist("leading_plus_nearest_mass")->Fill( bcc->cagenjets->at(0).v4().mass(),weight);

    //momenta
    Hist("gen_pt1")->Fill(bcc->cagenjets->at(0).pt(),weight);
    if(bcc->cagenjets->size()>1) Hist("gen_pt2")->Fill(bcc->cagenjets->at(1).pt(),weight);
    if(bcc->cagenjets->size()>2)  Hist("gen_pt3")->Fill(bcc->cagenjets->at(2).pt(),weight);
    if(bcc->cagenjets->size()>3)  Hist("gen_pt4")->Fill(bcc->cagenjets->at(3).pt(),weight);
    if(bcc->cagenjets->size()>4)  Hist("gen_pt5")->Fill(bcc->cagenjets->at(4).pt(),weight);

    Hist("gen_pt1_WP")->Fill(bcc->genjetswithparts->at(0).pt(),weight);
    if(bcc->genjetswithparts->size()>1) Hist("gen_pt2_WP")->Fill(bcc->genjetswithparts->at(1).pt(),weight);
    if(bcc->genjetswithparts->size()>2)  Hist("gen_pt3_WP")->Fill(bcc->genjetswithparts->at(2).pt(),weight);


    Hist("lepton_pt")->Fill(lepton.pt(),weight);

    Hist("top_pt")->Fill(hadr_top.pt(),weight);

    if(deltaPhiAbs(bcc->cagenjets->at(0).phi(), nearest.Phi())<1.5) Hist("leading_plus_nearest_pt")->Fill(leading_plus.Pt(),weight);
 
    if(deltaPhiAbs(bcc->cagenjets->at(0).phi(), nearest.Phi())<1.5) Hist("leading_plus_nearest_pt_diff")->Fill(leading_plus.Pt()-bcc->cagenjets->at(1).pt(),weight);
    //distances
    Hist("dr")->Fill(deltaR(bcc->cagenjets->at(0).v4(), hadr_top.v4()),weight);
    Hist("dr_jets")->Fill(deltaR(bcc->cagenjets->at(0).v4(),bcc->cagenjets->at(1).v4()),weight);
    Hist("drlep")->Fill(dr_jet_lep,weight);
    Hist("drlep2")->Fill(dr_2ndjet_lep,weight);
    Hist("dr_diff")->Fill(dr_diff,weight);
    Hist("leading_plus_nearest_dr")->Fill(dR_leading_nearest, weight);



    Hist("dPhi_2ndjet_lepton")->Fill(deltaPhiAbs(bcc->cagenjets->at(1).phi(),lepton.phi()), weight);

  //boost
    ((TH2F*) Hist("dPhi_2ndjet_lepton_vs_dR_2ndjet_lepton"))->Fill(dr_2ndjet_lep,deltaPhiAbs(bcc->cagenjets->at(1).phi(),lepton.phi()), weight);
    ((TH2F*) Hist("pTopHad_vs_dR_2ndjet_lepton"))->Fill(dr_2ndjet_lep,hadr_top.v4().P(), weight);
    ((TH2F*) Hist("dR_W_b_vs_dR_2ndjet_lepton"))->Fill(dr_2ndjet_lep,deltaR(hadr_top.daughter(genparts,1)->v4(), hadr_top.daughter(genparts,2)->v4()), weight);



    //angles
    Hist("eta1")->Fill(bcc->cagenjets->at(0).eta(),weight);
    if(bcc->cagenjets->size()>1) Hist("eta2")->Fill(bcc->cagenjets->at(1).eta(),weight);
    if(bcc->cagenjets->size()>2)  Hist("eta3")->Fill(bcc->cagenjets->at(2).eta(),weight);
    if(bcc->cagenjets->size()>3)  Hist("eta4")->Fill(bcc->cagenjets->at(3).eta(),weight);
    if(bcc->cagenjets->size()>4)  Hist("eta5")->Fill(bcc->cagenjets->at(4).eta(),weight);

    Hist("eta1_WP")->Fill(bcc->genjetswithparts->at(0).eta(),weight);
    if(bcc->genjetswithparts->size()>1) Hist("eta2_WP")->Fill(bcc->genjetswithparts->at(1).eta(),weight);
    if(bcc->genjetswithparts->size()>2)  Hist("eta3_WP")->Fill(bcc->genjetswithparts->at(2).eta(),weight);

    Hist("eta_lep")->Fill(lepton.eta(),weight);

    Hist("phi1")->Fill(bcc->cagenjets->at(0).phi(),weight);
    if(bcc->cagenjets->size()>1) Hist("phi2")->Fill(bcc->cagenjets->at(1).phi(),weight);
    if(bcc->cagenjets->size()>2)  Hist("phi3")->Fill(bcc->cagenjets->at(2).phi(),weight);
    if(bcc->cagenjets->size()>3)  Hist("phi4")->Fill(bcc->cagenjets->at(3).phi(),weight);
    if(bcc->cagenjets->size()>4)  Hist("phi5")->Fill(bcc->cagenjets->at(4).phi(),weight);

    Hist("phi1_WP")->Fill(bcc->genjetswithparts->at(0).phi(),weight);
    if(bcc->genjetswithparts->size()>1) Hist("phi2_WP")->Fill(bcc->genjetswithparts->at(1).phi(),weight);
    if(bcc->genjetswithparts->size()>2)  Hist("phi3_WP")->Fill(bcc->genjetswithparts->at(2).phi(),weight);

    Hist("phi_lep")->Fill(lepton.phi(),weight);

    Hist("leading_plus_nearest_dphi")->Fill(deltaPhiAbs(bcc->cagenjets->at(0).phi(), nearest.Phi()),weight);

  

    //numbers
    //   Hist("hadrlept")->Fill(hadr,weight);
    //   Hist("hadrlept")->Fill(lept+3,weight);
    Hist("NCAjets")->Fill(N_add_jets+2,weight);
    Hist("N_addCAjets")->Fill(N_add_jets, weight);

    //substructure
    // std::vector<GenParticle> genparts = *bcc->genparticles;
    GenJetWithParts genjet = bcc->genjetswithparts->at(0);
    int NParts = 0;
    for(unsigned int i = 0; i < genjet.genparticles_indices().size(); i++){
      NParts++;
      GenParticle genpart = genparts->at(genjet.genparticles_indices().at(i));
      Hist("Ptparts")->Fill(genpart.pt(),weight);
      Hist("DRparts")->Fill(deltaR(genjet.v4(), genpart.v4()),weight);
    }
    Hist("Nparts")->Fill(NParts,weight);


 
    TString mode= "" ;
    ((TH2F*)Hist("mass_vs_LeadingJetPt"+mode))->Fill(bcc->cagenjets->at(0).v4().mass(),bcc->cagenjets->at(0).pt(),weight);
    if(bcc->cagenjets->size()>1) ((TH2F*)Hist("mass_vs_2ndJetPt"+mode))->Fill(bcc->cagenjets->at(0).v4().mass(),bcc->cagenjets->at(1).pt(),weight);
    if(bcc->cagenjets->size()>2) ((TH2F*)Hist("mass_vs_3rdJetPt"+mode))->Fill(bcc->cagenjets->at(0).v4().mass(),bcc->cagenjets->at(2).pt(),weight);  
    ((TH2F*)Hist("mass_vs_TopPt"+mode))->Fill(bcc->cagenjets->at(0).v4().mass(),hadr_top.pt(),weight);    
    ((TH2F*)Hist("mass_vs_LeptonPt"+mode))->Fill(bcc->cagenjets->at(0).v4().mass(),lepton.pt(),weight);
    //  if(isEle) ((TH2F*)Hist("mass_vs_LeptonPt"+mode))->Fill(bcc->cagenjets->at(0).v4().mass(),electron->pt(),weight);

    delete ttbargen;

  }
}
void ComparisonHists::Finish()
{
  // final calculations, like division and addition of certain histograms
  //  EventCalc* calc = EventCalc::Instance();
  // bool IsRealData = calc->IsRealData();
  // if (IsRealData){
 
  // }

}

