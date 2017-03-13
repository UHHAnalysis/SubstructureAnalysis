#include "../include/JetMassTools.h"
#include "../../NtupleWriter/include/JetProps.h"
#include "fastjet/tools/Filter.hh"

#include <algorithm>
#include <memory>


void ScaleJetMass(double scale){

    EventCalc* calc = EventCalc::Instance();
    BaseCycleContainer *bcc = calc->GetBaseCycleContainer();

    for(unsigned int i = 0; i<bcc->topjets->size();  i++){
      double E = bcc->topjets->at(i).energy();
      double mass = bcc->topjets->at(i).v4().mass();
      bcc->topjets->at(i).set_energy(sqrt(pow(E,2)-(1-pow(scale,2))*pow(mass,2)));
      //cout << "mass1:  " << mass << "  mass2:  " << bcc->topjets->at(i).v4().mass() << "  mass2/mass1:  " << bcc->topjets->at(i).v4().mass()/mass << endl ;
    }

   

   //reset everything in EventCalc except the event weight
    EventCalc* evc = EventCalc::Instance();
    double weight = evc->GetWeight();
    double recweight = evc->GetRecWeight();
    double genweight = evc->GetGenWeight();
    evc->Reset();
    evc->ProduceWeight(weight);
    evc->ProduceRecWeight(recweight);
    evc->ProduceGenWeight(genweight);

}


//IsTagged(jet, btagtype)!!!!!!!!!!!
std::vector<Jet> GetBJets( std::vector<Jet> jets , E_BtagType btagtype, double pt_min, double eta_max){
 std::vector<Jet> bjets; 
  for(unsigned int i = 0; i < jets.size(); ++i){
    if(jets.at(i).pt() > pt_min && fabs(jets.at(i).eta()) < eta_max){
      bool btag = false;
      if(btagtype==e_CSVL && jets.at(i).btag_combinedSecondaryVertex()>0.244) btag = true;
      if(btagtype==e_CSVM && jets.at(i).btag_combinedSecondaryVertex()>0.679) btag = true;
      if(btagtype==e_CSVT && jets.at(i).btag_combinedSecondaryVertex()>0.898) btag = true;
      if(btag) bjets.push_back(jets.at(i));
    }
  }
  return bjets;
}

std::vector<Jet> GetNOBJets( std::vector<Jet> jets , E_BtagType btagtype, double pt_min, double eta_max){
  std::vector<Jet> NObjets; 
  for(unsigned int i = 0; i < jets.size(); ++i){
    if(jets.at(i).pt() > pt_min && fabs(jets.at(i).eta()) < eta_max){
      bool btag = false;
      if(btagtype==e_CSVL && jets.at(i).btag_combinedSecondaryVertex()>0.244) btag = true;
      if(btagtype==e_CSVM && jets.at(i).btag_combinedSecondaryVertex()>0.679) btag = true;
      if(btagtype==e_CSVT && jets.at(i).btag_combinedSecondaryVertex()>0.898) btag = true;
      if(!btag) NObjets.push_back(jets.at(i));
    }
  }
  return NObjets;
}



TTbar_gen_properties::TTbar_gen_properties(){
  EventCalc* calc = EventCalc::Instance();
  bcc = calc->GetBaseCycleContainer();
}

void TTbar_gen_properties::evaluate(){

  // const GenParticle* electron;
  // const GenParticle* muon;

  m_isMuon = false;
  m_isEle = false;
  
  if ( bcc->genparticles){
    std::vector<GenParticle>* genparts = bcc->genparticles;
    int Ngenparts = genparts->size();
    	  
    for(int l=0; l<Ngenparts ; l++){
      GenParticle genpart = genparts->at(l);
      if(genpart.status()<3){break;}
      if( fabs(genpart.pdgId()) == 6){
	//     if((fabs(genpart.daughter(genparts, 1)->pdgId())!=24)) daughter1 = genpart.daughter(genparts, 2);//W-bosons
	const  GenParticle* daughter1 = genpart.daughter(genparts, 1)->daughter(genparts, 1);
	const  GenParticle* daughter2 = genpart.daughter(genparts, 1)->daughter(genparts, 2);
	if((fabs(genpart.daughter(genparts, 1)->pdgId())!=24)){
	  daughter1 = genpart.daughter(genparts, 2)->daughter(genparts, 1);
	  daughter2 = genpart.daughter(genparts, 2)->daughter(genparts, 2);
	}
	if(daughter1 && daughter2){
	  if( fabs(daughter1->pdgId())<6) m_hadronic_top=genpart;
	  if( fabs(daughter1->pdgId()) == 13 ){
	    m_isMuon=true;
	    //   muon = daughter1;
	    m_lepton = daughter1; 
	  }
	  else if ( fabs(daughter2->pdgId()) == 13 ){
	    m_isMuon=true;
	    //  muon = daughter2;
	    m_lepton = daughter2; 
	  }
	  else if( fabs(daughter1->pdgId()) == 11 ){
	    m_isEle=true;
	    //  electron = daughter1;
	    m_lepton = daughter1;
	  }
	  else if( fabs(daughter2->pdgId()) == 11 ){
	    m_isEle=true;
	    //  electron = daughter2;
	    m_lepton = daughter2;
	  }
	}
      }
    }

  }

}


GenCleaner::GenCleaner( BaseCycleContainer* input)
{
    bcc = input;
}

GenCleaner::GenCleaner(){
  calc = EventCalc::Instance();
  bcc = calc->GetBaseCycleContainer();
}

void GenCleaner::CAGenJetCleaner(double ptmin, double etamax){

    std::vector<Particle> cleaned_jets;
    //  cout << "-----------------" << endl;
    for(unsigned int i=0; i < bcc->cagenjets->size(); ++i){
      //  cout << i << endl;
      //  cout << bcc->cagenjets->size() << endl;
      Particle jet = bcc->cagenjets->at(i);
      //  cout << "++++++++++++++" << endl;
      if(jet.pt()>ptmin) {
	if(fabs(jet.eta())<etamax) {          
	  cleaned_jets.push_back(jet);               
	}
      }
    }
   
    bcc->cagenjets->clear();
  
    for(unsigned int i=0; i<cleaned_jets.size(); ++i) {
       bcc->cagenjets->push_back(cleaned_jets[i]);
    }

    //cout <<"size: " << cleaned_jets.size() << endl;
    sort(bcc->cagenjets->begin(), bcc->cagenjets->end(), HigherPt());
 
    resetEventCalc();
}




void GenCleaner::GenCaJetLeptonSubstractor( bool sort){


  // removes all status 3 electrons and muons from the genjets (genjetswithparts and caXXgenjets)

  if(bcc->genjetswithparts){
    for(unsigned int i=0; i<bcc->genjetswithparts->size(); ++i) {

      LorentzVector jet_v4 = bcc->genjetswithparts->at(i).v4();
      std::vector<unsigned int> GENindices =  bcc->genjetswithparts->at(i).genparticles_indices();
      std::vector<GenParticle>* genparts = bcc->genparticles;

      bool LeptonInJet = false;

      for(unsigned int j = 0; j < genparts->size(); j++){
	if(genparts->at(j).status() < 3) continue; //loop only over status three particles
	if(fabs(genparts->at(j).pdgId()) == 11 || fabs(genparts->at(j).pdgId()) == 13 ){
	  const GenParticle lepton = genparts->at(j);
	  for(unsigned int k = 0 ; k < GENindices.size(); k++){
	    GenParticle genpart = genparts->at(GENindices.at(k));
	    if(lepton.pdgId() == genpart.pdgId() 
	       && fabs(1-lepton.pt()/genpart.pt()) < 0.3 
	       && fabs(1-lepton.eta()/genpart.eta()) < 0.3 
	       && fabs(1-lepton.phi()/genpart.phi()) < 0.3 
	       ){ LeptonInJet = true; break;}	
	  }
	  if(LeptonInJet) jet_v4 -= lepton.v4();
	}
      }

      bcc->genjetswithparts->at(i).set_v4(jet_v4);
    }
  }
 

  if(bcc->cagenjets){

    LorentzVector lepton1v4;
    LorentzVector lepton2v4;
    
    /*
    TTbarGen *ttbar = new TTbarGen(bcc);
    //const GenParticle lepton;
    if(ttbar->DecayChannel() == TTbarGen::e_ehad || ttbar->DecayChannel() == TTbarGen::e_muhad){
      const GenParticle lepton = ttbar->ChargedLepton();
      for(unsigned int i=0; i<bcc->cagenjets->size(); ++i) {
	LorentzVector jet_v4 = bcc->cagenjets->at(i).v4();
	if(deltaR(lepton.v4(),jet_v4)<1.2) jet_v4 -= lepton.v4();
	lepton1v4 = lepton.v4();
	bcc->cagenjets->at(i).set_v4(jet_v4);
      }
    }
    */
    //tested to give the same result for Pythia and Herwig samples!

    for(unsigned int i=0; i<bcc->cagenjets->size(); ++i) {
      LorentzVector jet_v4 = bcc->cagenjets->at(i).v4();
      std::vector<GenParticle>* genparts = bcc->genparticles;
      for(unsigned int j = 0; j < genparts->size(); j++){
	if(genparts->at(j).status() < 3) continue; //loop only over status three particles
	if(fabs(genparts->at(j).pdgId()) == 11 || fabs(genparts->at(j).pdgId()) == 13 ){
	  const GenParticle lepton = genparts->at(j);
	  if(deltaR(lepton.v4(),jet_v4)<1.2) jet_v4 -= lepton.v4();
	  lepton2v4 = lepton.v4();
	  //  if(!(lepton1v4 == lepton2v4)){
	  //   calc->PrintGenParticles();
	  //}
	}
      }
      bcc->cagenjets->at(i).set_v4(jet_v4);
    }
    


  }
  

  if(sort){
     if(bcc->cagenjets) std::sort(bcc->cagenjets->begin(), bcc->cagenjets->end(), HigherPt());
    if(bcc->genjetswithparts)std::sort(bcc->genjetswithparts->begin(), bcc->genjetswithparts->end(), HigherPt());
  }

  resetEventCalc();
}

void GenCleaner::resetEventCalc(){

    //reset everything in EventCalc except the event weight
    EventCalc* evc = EventCalc::Instance();
    double weight = evc->GetWeight();
    double recweight = evc->GetRecWeight();
    double genweight = evc->GetGenWeight();
    evc->Reset();
    evc->ProduceWeight(weight);
    evc->ProduceRecWeight(recweight);
    evc->ProduceGenWeight(genweight);
}



int Wjet_index(std::vector<Jet> *jets , std::vector<TopJet> *topjets, LorentzVector lepton, bool massdiff ){

  int Wjet_index = -1; 
  int Ntag = 0;

  /*  std::vector<Jet>* bjets; 

      for(int i = 0; i < jets->size(); i++){
      if(jets->at(i).btag_combinedSecondaryVertex()>0.898){
      Ntag ++;
      //bjets->push_back(jets->at(i));
      }
      else cout<< "WARNING: more than two btags" << endl;
      }
  */

  int Nbjets = 0;

  int indexB1 = -1;
  int indexB2 = -1;

  for(unsigned int i = 0;( i<3 && i<topjets->size()); i++){
    for(unsigned int k = 0; k < jets->size(); k++){
      if(jets->at(k).btag_combinedSecondaryVertex()>0.244 && deltaR(topjets->at(i).v4(), jets->at(k).v4())<1.2)
	//tight: 0.898
	//loose:0.244
	{
	  Nbjets++;
	  if(indexB1 < 0)indexB1 = i;
	  if(indexB1 >= 0)indexB2 = i;
	}
    }
  }

  if(Nbjets == 2) { 
    for(unsigned int i = 0; (i < 3 &&  i<topjets->size()); i++){
      if(i != indexB1 && i != indexB2) Wjet_index = i;
    } 
    return Wjet_index;
  }
  
  /*  if(Nbjets > 2) cout<<"WARNING: more than two bjets" << endl;
      if(Nbjets < 2) cout<<"WARNING: less than two bjets" << endl;
      if(Wjet_index >3) cout<<"WARNING: no Wjet identifyed" <<endl;
  */
  else if(Nbjets == 1 && massdiff){
    int i_jet1 = -1;
    int i_jet2 = -1;
    for(unsigned int i = 0; (i < 3 &&  i<topjets->size()); i++){
      if(i != indexB1 && i_jet1 < 0) i_jet1 = i;
      if(i != indexB1 && i_jet1 >= 0) i_jet2 = i;
      if(i_jet1 >= 0 && i_jet2 >= 0){
	double md1 = (topjets->at(indexB1).v4()+lepton).M()-(topjets->at(i_jet1).v4()+topjets->at(i_jet2).v4()).M();
	double md2 = (topjets->at(i_jet1).v4()+lepton).M()-(topjets->at(indexB1).v4()+topjets->at(i_jet2).v4()).M();
	double md3 = (topjets->at(i_jet2).v4()+lepton).M()-(topjets->at(indexB1).v4()+topjets->at(i_jet1).v4()).M();
	if(md3 < md2 && md3 < md1) Wjet_index = i_jet1;
	if(md2 < md3 && md2 < md1) Wjet_index = i_jet2;
	if(md1 < md3 && md1 < md2){
	  //	  cout << "WARNING: chice of W-jet ambigous!" << endl;
	  if(topjets->at(i_jet1).v4().M()>topjets->at(i_jet2).v4().M()) Wjet_index = i_jet1;
	  else  Wjet_index = i_jet2;
	}
      }
    } 
    return Wjet_index;      
  }
  else return -1;
}


fastjet::PseudoJet FilterJet(TopJet t, std::vector<PFParticle>* pfparticles){
 
   JetProps jp(&t, pfparticles);

   std::vector<fastjet::PseudoJet> jets = jp.GetFastJet(2.0);

   fastjet::Filter trimm_2_003(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.003));

   // first: create a fastjet-jet from the TopJet
   fastjet::PseudoJet trimmed03_jet = trimm_2_003(jets[0]);

   return trimmed03_jet;
}


void FilterNTopJets(unsigned int Njets){
  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  for(unsigned int i=0; (i<bcc->topjets->size() && i< Njets); ++i) {
      
	std::vector<PFParticle> *PFparticles = bcc->pfparticles;

 
	JetProps jp(&bcc->topjets->at(i), PFparticles);

	std::vector<fastjet::PseudoJet> jets = jp.GetFastJet(2.0);

	//	fastjet::Filter filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3)); //recommended
	fastjet::Filter filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03));//recommended

	if(jets.size()){
	fastjet::PseudoJet filtered_jet = filter(jets[0]);

	LorentzVector jet_v4;
	jet_v4.SetPt(filtered_jet.pt());
 	jet_v4.SetPhi(filtered_jet.phi());  
	jet_v4.SetEta(filtered_jet.eta());
 	jet_v4.SetE(filtered_jet.E());
  
	//jet energy correction
	jet_v4 *= 1/bcc->topjets->at(i).JEC_factor_raw();

	bcc->topjets->at(i).set_v4(jet_v4);
	}
    }
}

void SoftDropNTopJets(unsigned int Njets){
  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  for(unsigned int i=0; (i<bcc->topjets->size() && i< Njets); ++i) {
      
    std::vector<PFParticle> *PFparticles = bcc->pfparticles;

    SoftDrop* Softdrop;
    Softdrop = new SoftDrop();

    fastjet::PseudoJet groomed_jet;
    Softdrop->GetJet(bcc->topjets->at(i),PFparticles,0.1,0.0,groomed_jet);

    delete Softdrop;

    LorentzVector jet_v4;
    jet_v4.SetPt(groomed_jet.pt());
    jet_v4.SetPhi(groomed_jet.phi());  
    jet_v4.SetEta(groomed_jet.eta());
    jet_v4.SetE(groomed_jet.E());
  
    //jet energy correction
    jet_v4 *= 1/bcc->topjets->at(i).JEC_factor_raw();

    bcc->topjets->at(i).set_v4(jet_v4);
  }
}


void SoftDropNGenJetsWithParts(unsigned int Njets){
  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  for(unsigned int i=0; (i<bcc->genjetswithparts->size() && i< Njets); ++i) {
      
    std::vector<GenParticle> *GENparticles = bcc->genparticles;

    SoftDrop* Softdrop;
    Softdrop = new SoftDrop();

    fastjet::PseudoJet groomed_jet;
    Softdrop->GetJet(bcc->genjetswithparts->at(i),GENparticles,0.1,0.0,groomed_jet);

    delete Softdrop;

    LorentzVector jet_v4;
    jet_v4.SetPt(groomed_jet.pt());
    jet_v4.SetPhi(groomed_jet.phi());  
    jet_v4.SetEta(groomed_jet.eta());
    jet_v4.SetE(groomed_jet.E());

    bcc->genjetswithparts->at(i).set_v4(jet_v4);
  }
}





void TrimmNTopJets(unsigned int Njets){
  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  for(unsigned int i=0; (i<bcc->topjets->size() && i< Njets); ++i) {
      
	std::vector<PFParticle> *PFparticles = bcc->pfparticles;

 
	JetProps jp(&bcc->topjets->at(i), PFparticles);

	std::vector<fastjet::PseudoJet> jets = jp.GetFastJet(2.0);

	fastjet::Filter trimm_2_003(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.003));


	if(jets.size()){
	fastjet::PseudoJet trimmed03_jet = trimm_2_003(jets[0]);

	LorentzVector jet_v4;
	jet_v4.SetPt(trimmed03_jet.pt());
 	jet_v4.SetPhi(trimmed03_jet.phi());  
	jet_v4.SetEta(trimmed03_jet.eta());
 	jet_v4.SetE(trimmed03_jet.E());
  
	//jet energy correction
	jet_v4 *= 1/bcc->topjets->at(i).JEC_factor_raw();

	bcc->topjets->at(i).set_v4(jet_v4);
	}
    }
    


}

