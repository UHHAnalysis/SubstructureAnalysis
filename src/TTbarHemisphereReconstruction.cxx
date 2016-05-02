#include "include/TTbarHemisphereReconstruction.h"

using namespace std;

TTbarHemisphereReconstruction::TTbarHemisphereReconstruction( BaseCycleContainer* bcc )
{
  m_bcc = bcc;


  // set default values
  m_bJetLepCSV = e_CSVM;
  m_bJetLepPt = 30;
  m_bJetLepEta = 2.1;

  m_bJetHadCSV = e_CSVL;
  m_bJetHadPt = 30;
  m_bJetHadEta = 2.1;

  m_CAbJetPt = 50;
  m_CAbJetEta = 2.4;

  m_WJetPt = 90;//90
  m_WJetEta = 2.4;

  m_LeptonicHemisphereFound = false;
  m_HadronicHemisphereFound = false;

}


TTbarHemisphereReconstruction::~TTbarHemisphereReconstruction()
{

}




bool TTbarHemisphereReconstruction::ReconstructLeptonicHemisphere(){
  
  //get the primary charged lepton
  EventCalc* calc = EventCalc::Instance();
  m_Lepton = *calc->GetPrimaryLepton();

  //neutrino reconstruction
  std::vector<LorentzVector> neutrinos = calc->NeutrinoReconstruction(m_Lepton.v4(), m_bcc->met->v4());
  if(neutrinos.size() == 0) return false; 

  if(!m_bcc->jets){
    cout<< "no jets!! -> no defintion of leptonic hemisphere possible" << endl;
    return false;
  }

  //get b-jets 
  std::vector<Jet> bjets = GetBJets( *m_bcc->jets, m_bJetLepCSV, m_bJetLepPt, m_bJetLepEta); 


  Jet bJetCandidate;
  //find the b-jet as the closest to the lepton but still dR(b,l) < pi/2
  double mindR = FLT_MAX;
  for(unsigned int i = 0; i < bjets.size(); ++i){
    Jet& ijet = bjets.at(i);

    double dR = deltaR(ijet.v4(), m_Lepton.v4());

    // if(dR > TMath::Pi()/2.) continue;

    if(dR < mindR){
      mindR = dR;
      bJetCandidate = ijet;
      if(dR < TMath::Pi()/2.) m_bJetLep = ijet;
    } 
  }

  //debugging
  if(!m_bJetLep.pt()){
     /*//debugging
    cout<< "no bjet found close to the lepton  pt: " << m_bJetLep.pt() << "   dR: "<< mindR << endl;
    cout<< "b-jet candiadate  pt: " << bJetCandidate.pt() << " ,  eta: " << bJetCandidate.eta() << endl;

    cout << "--------------------------------------------------------" << endl;
 cout<< "  run "<< m_bcc->run << " event " << m_bcc->event << endl;
    cout << "--------------------------------------------------------" << endl;
    for(unsigned int i = 0; i < m_bcc->jets->size(); ++i){
      cout<< "  jet "<< i << " pt: " << m_bcc->jets->at(i).pt() << "  eta: "<< m_bcc->jets->at(i).eta()
	  << "  CSV: " << m_bcc->jets->at(i).btag_combinedSecondaryVertex()<< endl;
    }
    cout << "--------------------------------------------------------" << endl;
     */
    return false; 
  }

  //leptonic top reconstruction
  double mtop = 173.34; // !!!!!!!!
  LorentzVector CombinedVector = neutrinos.at(0) + m_bJetLep.v4() + m_Lepton.v4();
  m_Neutrino = neutrinos.at(0);
  if(neutrinos.size() > 1){
    LorentzVector CombinedVector2 = neutrinos.at(1) + m_bJetLep.v4() + m_Lepton.v4();
    if(fabs(CombinedVector2.mass()-mtop) < fabs(CombinedVector.mass()-mtop)){
      CombinedVector = CombinedVector2;
      m_Neutrino = neutrinos.at(1);
    }
  }
  if(CombinedVector.mass() <= 0 || CombinedVector.mass() > 230) return false;
  m_TopLep = CombinedVector;
 
  //just one b-jet on the leptoic hemisphere
  int number = 0; 
  for(unsigned int i=0; i<bjets.size(); ++i){
    Jet& ijet = bjets.at(i);
    if(deltaR(ijet.v4(),CombinedVector) < TMath::Pi()/2.) number ++;
  }
  m_NbJetsLeptonicHemisphere = number;

  if (number == 1){
    m_LeptonicHemisphereFound = true;
    return true;
  }

  return false;
}

bool TTbarHemisphereReconstruction::ReconstructHadronicHemisphere(){
 
  if(!m_LeptonicHemisphereFound) return false; 

  //get b-jets 
  std::vector<Jet> bjets = GetBJets( *m_bcc->jets, m_bJetHadCSV, m_bJetHadPt, m_bJetHadEta); 
  //get non b-tagged jets 
  std::vector<Jet> NObjets = GetNOBJets( *m_bcc->jets,  m_bJetHadCSV, 25, m_bJetHadEta);//25 

  sort(bjets.begin(), bjets.end(), HigherPt()); 
  sort(NObjets.begin(), NObjets.end(), HigherPt()); 

  //find the hadronic b jet (leading b-jet in the hadronic hemisphere)
  for(unsigned int i = 0; i < bjets.size(); i++){
    Jet& bjet = bjets.at(i);
    if(deltaR(bjet.v4(), m_TopLep) > TMath::Pi()/2.){
      m_bJetHad = bjet;
      break;
    }
  }
  if(!m_bJetHad.pt()) return false;

 
 //find the first quark jet close to the bjet
  double mindR = FLT_MAX;
  for(unsigned int i = 0; i < NObjets.size(); i++){
    Jet& jet = NObjets.at(i);
    if(deltaR(jet.v4(), m_TopLep) < TMath::Pi()/2.) continue;
    double dR_jet_bhad = deltaR(jet.v4(), m_bJetHad.v4());
    if( dR_jet_bhad < mindR){
      mindR = dR_jet_bhad;
      m_WJet1 = jet; 
    }
  }
  if(!m_WJet1.pt()) return false;

  //find the second quark jet close to the first one 
  double mindR2 = FLT_MAX;
  for(unsigned int i = 0; i < NObjets.size(); i++){
    Jet& jet = NObjets.at(i);
    if(deltaR(jet.v4(), m_TopLep) < TMath::Pi()/2.) continue;
    double dR_jet_WJet1 = deltaR(jet.v4(), m_WJet1.v4());
    if( dR_jet_WJet1 < mindR2 && dR_jet_WJet1 > 0){
      mindR2 = dR_jet_WJet1;
      m_WJet2 = jet; 
    }
  }
  if(!m_WJet2.pt()) return false;

  //find CA8 jets for the W and the b quark 
  std::vector<TopJet> cajets = *m_bcc->topjets;

  sort(cajets.begin(), cajets.end(), HigherPt());

  int NCAJets = 0;

  for(unsigned int i = 0; i < cajets.size(); i++){
    TopJet& cajet = cajets.at(i); 
    if(!m_WJet.pt()){
      if(cajet.pt() > m_WJetPt && fabs(cajet.eta()) < m_WJetEta 
	 && deltaR(cajet.v4(),  m_TopLep) > TMath::Pi()/2.
       	 && deltaR(cajet.v4(), m_WJet1.v4()) < 1.2
       	 && deltaR(cajet.v4(), m_WJet2.v4()) < 1.2
       	 && deltaR(cajet.v4(), m_bJetHad.v4()) > 1.2
	 ){
	m_WJet = cajet;
	NCAJets++;
	continue;
      }
    }
    if(!m_CAbJet.pt()){
      if(cajet.pt() > m_CAbJetPt && fabs(cajet.eta()) < m_CAbJetEta 
	 && deltaR(cajet.v4(),  m_TopLep) > TMath::Pi()/2.
       	 && deltaR(cajet.v4(), m_WJet1.v4()) > 1.2
	 && deltaR(cajet.v4(), m_WJet2.v4()) > 1.2
	 && deltaR(cajet.v4(), m_bJetHad.v4()) < 1.2
	 ){
	m_CAbJet = cajet;
	NCAJets++;
	continue;
      }
    }

    if(cajet.pt() > m_CAbJetPt && fabs(cajet.eta()) < m_CAbJetEta 
       && deltaR(cajet.v4(),  m_TopLep) > TMath::Pi()/2.) NCAJets++; 
  }
 
  if( m_WJet.pt() && NCAJets < 3){
    m_HadronicHemisphereFound = true;
    return true;
  }

  return false;
}


void TTbarHemisphereReconstruction::SetLeptonicBJetPtEtaCSV( double pt, double eta, E_BtagType btagtype){
  m_bJetLepCSV = btagtype;
  m_bJetLepPt = pt;
  m_bJetLepEta = eta;
}

void TTbarHemisphereReconstruction::SetHadronicBJetPtEtaCSV( double pt, double eta, E_BtagType btagtype ){
  m_bJetHadCSV = btagtype;
  m_bJetHadPt = pt;
  m_bJetHadEta = eta;
}

void TTbarHemisphereReconstruction::SetCABJetPtEta( double pt, double eta){
  m_CAbJetPt = pt; 
  m_CAbJetEta = eta;
}

void TTbarHemisphereReconstruction::SetWJetPtEta(double pt, double eta){
  m_WJetPt = pt;
  m_WJetEta = eta;
}


bool TTbarHemisphereReconstruction::LeptonicHemisphereFound(){
  return m_LeptonicHemisphereFound;
}

bool TTbarHemisphereReconstruction::HadronicHemisphereFound(){
 return m_HadronicHemisphereFound;
}

LorentzVector TTbarHemisphereReconstruction::TopLep(){
  LorentzVector Tlep = m_bJetLep.v4() + m_Lepton.v4() + m_Neutrino;
  return Tlep;
}

LorentzVector TTbarHemisphereReconstruction::TopHad(){
  LorentzVector Thad = m_bJetHad.v4()+ m_WJet1.v4()+ m_WJet2.v4();
  return Thad;
}

LorentzVector TTbarHemisphereReconstruction::Neutrino(){return m_Neutrino;}

Particle TTbarHemisphereReconstruction::Lepton(){return m_Lepton;}

Jet TTbarHemisphereReconstruction::bJetLep(){return m_bJetLep;}

Jet TTbarHemisphereReconstruction::bJetHad(){return m_bJetHad;}

LorentzVector TTbarHemisphereReconstruction::WRecoLep(){
  LorentzVector W =  m_Lepton.v4() + m_Neutrino;
  return W;
}

LorentzVector TTbarHemisphereReconstruction::WRecoHad(){
  LorentzVector W = m_WJet1.v4() + m_WJet2.v4();
  return W;
}

TopJet TTbarHemisphereReconstruction::WJet(){return m_WJet;}

TopJet TTbarHemisphereReconstruction::CAbJet(){return m_CAbJet;}

Jet TTbarHemisphereReconstruction::WJet1(){return m_WJet1;}

Jet TTbarHemisphereReconstruction::WJet2(){return m_WJet2;}
