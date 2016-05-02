//---sframe new---- 
#include "include/SubstructureSelectionModules.h"
#include "TLorentzVector.h"

SubstructureGeneratorPreSelection::SubstructureGeneratorPreSelection(double pt_hadonic_top)
{
  m_pt_hadonic_top = pt_hadonic_top;
}

bool SubstructureGeneratorPreSelection::pass(BaseCycleContainer *bcc)
{

  bool sel = false;
  TTbarGen *ttbar = new TTbarGen(bcc);
  if(ttbar->DecayChannel() == TTbarGen::e_muhad || ttbar->DecayChannel() == TTbarGen::e_ehad ) sel = true;
  delete ttbar;

  if(sel) return true; 
  return false;

}

std::string SubstructureGeneratorPreSelection::description()
{
    return "semileptonic ttbar-decays";
}

SubstructureGeneratorPreSelectionMuon::SubstructureGeneratorPreSelectionMuon(double pt_hadonic_top)
{
  m_pt_hadonic_top = pt_hadonic_top;
}

bool SubstructureGeneratorPreSelectionMuon::pass(BaseCycleContainer *bcc)
{

  bool muhad = false;
  TTbarGen *ttbar = new TTbarGen(bcc);
  if(ttbar->DecayChannel() == TTbarGen::e_muhad) muhad = true;
  delete ttbar;

  if(muhad) return true; 
  return false;

 
}

std::string SubstructureGeneratorPreSelectionMuon::description()
{
    return "semileptonic ttbar-decays including a muon";
}


SubstructureGeneratorPreSelectionElectron::SubstructureGeneratorPreSelectionElectron(double pt_hadonic_top)
{
  m_pt_hadonic_top = pt_hadonic_top;
}

bool SubstructureGeneratorPreSelectionElectron::pass(BaseCycleContainer *bcc)
{

  bool ehad = false;
  TTbarGen *ttbar = new TTbarGen(bcc);
  if(ttbar->DecayChannel() == TTbarGen::e_ehad) ehad = true;
  delete ttbar;

  if(ehad) return true; 
  return false;
 
 
}

std::string SubstructureGeneratorPreSelectionElectron::description()
{
    return "semileptonic ttbar-decays including an electron";
}



SubstructureGeneratorPreSelectionHadron::SubstructureGeneratorPreSelectionHadron()
{

}

bool SubstructureGeneratorPreSelectionHadron::pass(BaseCycleContainer *bcc)
{

  bool sel = false;
  TTbarGen *ttbar = new TTbarGen(bcc);
  if(ttbar->DecayChannel() == TTbarGen::e_had) sel = true;
  delete ttbar;

  if(sel) return true; 
  return false;

}

std::string SubstructureGeneratorPreSelectionHadron::description()
{
    return "hadronic ttbar-decays";
}








GenLeptonSelection::GenLeptonSelection(double pt_lep, double eta_max)
{
  m_pt_lep = pt_lep;
  m_eta_max=eta_max;
}

bool GenLeptonSelection::pass(BaseCycleContainer *bcc)
{

 bool sel = false;
  TTbarGen *ttbar = new TTbarGen(bcc);
  if(ttbar->DecayChannel() == TTbarGen::e_ehad || ttbar->DecayChannel() == TTbarGen::e_muhad  ){
    if(ttbar->ChargedLepton().pt() > m_pt_lep && fabs(ttbar->ChargedLepton().eta()) < m_eta_max) sel = true;
  }
  delete ttbar;

  return sel; 

}

std::string GenLeptonSelection::description()
{
    return "semileptonic ttbar-decays";
}





/*NCAGenJetSelection::NCAGenJetSelection( double ptmin, double ptmin2, double etamax)
{
    m_ptmin2=ptmin2;
    m_ptmin=ptmin;
    m_etamax=etamax;
}

bool NCAGenJetSelection::pass(BaseCycleContainer *bcc)
{
  if(!bcc->cagenjets){return false;}
    int Njetacc=0;
    int Njetacc2=0; 
      for(unsigned int i=0; i<bcc->cagenjets->size(); ++i) {
        if(bcc->cagenjets->at(i).pt()>m_ptmin && bcc->cagenjets->at(i).eta()<m_etamax) Njetacc++;
	if(bcc->cagenjets->at(i).pt()>m_ptmin2 && bcc->cagenjets->at(i).eta()<m_etamax) Njetacc2++;
     }

    if (Njetacc >= 1 && Njetacc2 >= 2){return true;}
    if (Njetacc >= 2 && Njetacc2 >= 2){return true;}

 return false;
 }

std::string NCAGenJetSelection::description()
{
    char s[100];
 sprintf(s, " at least 1 GenTopJet with pt>%.1f GeV, 1 GenTopJet with pt>%.1f GeV, |eta|<%.1f",m_ptmin,m_ptmin2,m_etamax);

      return s;
}*/

NCAGenJetSelection::NCAGenJetSelection(int min_nparticle, int max_nparticle, double ptmin, double etamax)
{
    m_min_nparticle=min_nparticle;
    m_max_nparticle=max_nparticle;
    m_ptmin=ptmin;
    m_etamax=etamax;
}

bool NCAGenJetSelection::pass(BaseCycleContainer *bcc)
{
    int nparticle=0;
    for(unsigned int i=0; i<bcc->cagenjets->size(); ++i) {
        if(bcc->cagenjets->at(i).pt()>m_ptmin && fabs(bcc->cagenjets->at(i).eta())<m_etamax) nparticle++;
    }
    return nparticle >= m_min_nparticle && nparticle<=m_max_nparticle;
}

std::string NCAGenJetSelection::description()
{
    char s[100];
    sprintf(s, "%d <= number of cagenjets <= %d, with pt>%.1f GeV, abs(eta)<%.1f",m_min_nparticle,m_max_nparticle,m_ptmin,m_etamax);

    return s;
}


NCAGenJetveto::NCAGenJetveto( double ptmin, double njets, double etamax)
{

    m_ptmin=ptmin;
    m_njets=njets;
    m_etamax=etamax;
}

bool NCAGenJetveto::pass(BaseCycleContainer *bcc)
{
  if(!bcc->cagenjets){return false;}
    int Njetacc=0;
      for(unsigned int i=0; i<bcc->cagenjets->size(); ++i) {
        if(bcc->cagenjets->at(i).pt()>m_ptmin && fabs(bcc->cagenjets->at(i).eta())<m_etamax) Njetacc++;
     }

    if (Njetacc == m_njets){return true;}


 return false;
}

std::string NCAGenJetveto::description()
{
    char s[100];
    sprintf(s, "<= %.1f cagenjets with pt>%.1f GeV, |eta|<%.1f",m_njets,m_ptmin,m_etamax);

    return s;
}



NTopJetveto::NTopJetveto( double ptmin, double njets, double etamax)
{

    m_ptmin=ptmin;
    m_njets=njets;
    m_etamax=etamax;
}

bool NTopJetveto::pass(BaseCycleContainer *bcc)
{
  if(!bcc->topjets){return false;}
    int Njetacc=0;
      for(unsigned int i=0; i<bcc->topjets->size(); ++i) {
        if(bcc->topjets->at(i).pt()>m_ptmin && fabs(bcc->topjets->at(i).eta())<m_etamax) Njetacc++;
     }

    if (Njetacc == m_njets){return true;}


 return false;
}

std::string NTopJetveto::description()
{
    char s[100];
    sprintf(s, "<= %.1f TopJets with pt>%.1f GeV, |eta|<%.1f",m_njets,m_ptmin,m_etamax);

    return s;
}




LeptonicHemisphereSelection::LeptonicHemisphereSelection(double jetpt, double abseta, E_BtagType btagtype, double dR )
{
  m_jetpt = jetpt;
  m_abseta = abseta;
  m_dR = dR;
  m_type = btagtype;
}

bool LeptonicHemisphereSelection::pass(BaseCycleContainer* bcc)
{  
 bool reconstructed = false;  

 TTbarHemisphereReconstruction* tt = new TTbarHemisphereReconstruction(bcc);
 if(tt->ReconstructLeptonicHemisphere()) reconstructed = true;

 delete tt;
 return reconstructed;
}

std::string LeptonicHemisphereSelection::description()
{
  char s[500];
  sprintf(s, "Selection of events which allow a reconstruction of a leptonic hemisphere");
  return s;
}



BJetCloseToLeptonSelection::BJetCloseToLeptonSelection(double jetpt, double abseta, E_BtagType btagtype, double dR){
  m_jetpt = jetpt;
  m_abseta = abseta;
  m_dR = dR;
  m_type = btagtype;
}

bool BJetCloseToLeptonSelection::pass(BaseCycleContainer* bcc){

  EventCalc* calc = EventCalc::Instance();
  Particle lepton = *calc->GetPrimaryLepton();

  std::vector<Jet> bjets = GetBJets(*bcc->jets,m_type,m_jetpt,m_abseta);

  for(unsigned int i = 0; i < bjets.size(); i++){
    if(deltaR(bjets.at(i).v4(), lepton.v4()) < m_dR){
      return true;
    }
  }
  return false; 
}

std::string BJetCloseToLeptonSelection::description()
{
  char s[500];
  sprintf(s, "Slection of events with a b-jet (pt > %f, |eta| < %f) closer than deltaR = %f to the lepton.", m_jetpt, m_abseta, m_dR);
  return s;
}


HadronicHemisphereSelection::HadronicHemisphereSelection(double jetpt, double abseta, E_BtagType btagtype, double dR )
{
  m_jetpt = jetpt;
  m_abseta = abseta;
  m_dR = dR;
  m_type = btagtype;
}

bool HadronicHemisphereSelection::pass(BaseCycleContainer* bcc)
{  
 bool reconstructed = false;  

  TTbarHemisphereReconstruction* tt = new TTbarHemisphereReconstruction(bcc);
  tt->ReconstructLeptonicHemisphere();
  if(tt->ReconstructHadronicHemisphere()) reconstructed = true;

  delete tt;
  return reconstructed;
}

std::string HadronicHemisphereSelection::description()
{
  char s[500];
  sprintf(s, "Selection of events which allow a reconstruction of a hadronic hemisphere");
  return s;
}


LeadingTopJetHigherMassSelection::LeadingTopJetHigherMassSelection(TString mode)
{
  if(mode != "Default" && mode != "Lepton"){
    cout << "Unknown mode for the LeadingTopJetHigherMassSelection -- shold be Default or Lepton -- Default will be applied!" << endl;
    m_mode = "Default";
  }
  else m_mode = mode;
}

bool LeadingTopJetHigherMassSelection::pass(BaseCycleContainer* bcc)
{
 if(bcc->topjets->size()>1){

   double mass1 = bcc->topjets->at(0).v4().mass();
   double mass2 = int_infinity(); 

   if(m_mode == "Default")
     mass2 = bcc->topjets->at(1).v4().mass();
   if(m_mode == "Lepton")
     {
       EventCalc* calc = EventCalc::Instance();
       Particle Lepton = *calc->GetPrimaryLepton();
       mass2 = (bcc->topjets->at(1).v4()+Lepton.v4()).mass();
     }
   if(mass1 > mass2) return true;
 }
 return false;   
}

std::string LeadingTopJetHigherMassSelection::description()
{
  char s[500];
  if(m_mode == "Default")
    sprintf(s, "Selection of events with m(leading CAJet) > m(2nd CAJet)");
  if(m_mode == "Lepton")
    sprintf(s, "Selection of events with m(leading CAJet) > m(2nd CAJet + lepton)");
  return s;
}



LeadingGenCAJetHigherMassSelection::LeadingGenCAJetHigherMassSelection(TString mode)
{
  if(mode != "Default" && mode != "Lepton" && mode != "LeptonAndNeutrino"){
    cout << "Unknown mode for the LeadingTopJetHigherMassSelection -- shold be Default, Lepton or LeptonAndNeutrino -- Default will be applied!" << endl;
    m_mode = "Default";
  }
  else m_mode = mode;
}

bool LeadingGenCAJetHigherMassSelection::pass(BaseCycleContainer* bcc)
{
  if(bcc->cagenjets->size()>1){

    double mass1 = bcc->cagenjets->at(0).v4().mass();
    double mass2 = int_infinity();

    if(m_mode == "Default")
      mass2 = bcc->cagenjets->at(1).v4().mass(); 
    else
      {
	TTbarGen *ttbar = new TTbarGen(bcc); 
	if(!(ttbar->DecayChannel() == TTbarGen::e_ehad) && !(ttbar->DecayChannel() == TTbarGen::e_muhad)) cout << "wrong channel!!!!!!!!!!"<<endl;

	if(m_mode == "Lepton") mass2 = (bcc->cagenjets->at(1).v4()+ttbar->ChargedLepton().v4()).mass() ;
	if(m_mode == "LeptonAndNeutrino") mass2 = (bcc->cagenjets->at(1).v4()+ttbar->ChargedLepton().v4()+ttbar->Neutrino().v4()).mass();
	delete ttbar;
      }  
    if(mass1 > mass2) return true;
  }
  return false;
}

std::string LeadingGenCAJetHigherMassSelection::description()
{
  char s[500];
  if(m_mode == "Default")
    sprintf(s, "Selection of events with m(leading CAGenJet) > m(2nd CAGenJet)");
  if(m_mode == "Lepton")
    sprintf(s, "Selection of events with m(leading CAGenJet) > m(2nd CAGenJet + lepton)");
  if(m_mode == "LeptonAndNeutrino")
    sprintf(s, "Selection of events with m(leading CAGenJet) > m(2nd CAGenJet + lepton + neutrino)");
  return s;
} 





DR2ndGenCAJetLeptonSelection::DR2ndGenCAJetLeptonSelection( double dR_max)
{
    m_dR_max = dR_max;
}

bool DR2ndGenCAJetLeptonSelection::pass(BaseCycleContainer *bcc)
{
  if(!bcc->cagenjets) return false;

  TTbarGen *ttbar = new TTbarGen(bcc);
  double dR = 100.;
  if(ttbar->DecayChannel() == TTbarGen::e_ehad || ttbar->DecayChannel() == TTbarGen::e_muhad)
    {
      dR = deltaR(bcc->cagenjets->at(1).v4(), ttbar->ChargedLepton().v4());
    }
  delete ttbar;

  if(dR < m_dR_max) return true;
  return false;
}

std::string DR2ndGenCAJetLeptonSelection::description()
{
    char s[100];
    sprintf(s, "deltaR(2nd cagenjet, lepton) < %.1f",m_dR_max);
    return s;
}


DR2ndTopJetLeptonSelection::DR2ndTopJetLeptonSelection( double dR_max)
{
    m_dR_max = dR_max;
}

bool DR2ndTopJetLeptonSelection::pass(BaseCycleContainer *bcc)
{
  if(!bcc->topjets) return false;
  EventCalc* calc = EventCalc::Instance();
  Particle Lepton = *calc->GetPrimaryLepton();

  double dR = deltaR(bcc->topjets->at(1).v4(), Lepton.v4());

  if(dR < m_dR_max)return true;
  return false;
}

std::string DR2ndTopJetLeptonSelection::description()
{
    char s[100];
    sprintf(s, "deltaR(2nd TopJet, lepton) < %.1f",m_dR_max);
    return s;
}


GenCAJetMatchedSelection::GenCAJetMatchedSelection( double dR_max, TString channel )
{
    m_dR_max = dR_max;
    m_channel = channel;
}

bool GenCAJetMatchedSelection::pass(BaseCycleContainer *bcc)
{
  if(!bcc->cagenjets->size()>0) return false;

  TTbarGen *ttbar = new TTbarGen(bcc); 
  double dR = 100;
  if(m_channel == "LeptonPlusJets"){
    if(!ttbar->DecayChannel() == TTbarGen::e_ehad && !ttbar->DecayChannel() == TTbarGen::e_muhad) 
      cout << "wrong channel in matched selection"<<endl;
    dR = deltaR(bcc->cagenjets->at(0).v4(), ttbar->TopHad().v4());
  }
  else if(m_channel == "AllHadronic"){
    dR = deltaR(bcc->cagenjets->at(0).v4(), ttbar->Top().v4());
    double dR2 = deltaR(bcc->cagenjets->at(0).v4(), ttbar->Antitop().v4());
    if(dR2 < dR) dR = dR2;
  }
  else{
    cout << "no channel selected for the matching. All events will be treated as mismatched" << endl;
    return false;
  }
  delete ttbar;

  if(dR < m_dR_max)return true;
  return false;
}

std::string GenCAJetMatchedSelection::description()
{
    char s[100];
    sprintf(s, "deltaR(2nd cagenjet, top) < %.1f",m_dR_max);
    return s;
}




TopJetMatchedSelection::TopJetMatchedSelection( double dR_max)
{
    m_dR_max = dR_max;
}

bool TopJetMatchedSelection::pass(BaseCycleContainer *bcc)
{
  if(!bcc->topjets->size()) return false;

  TTbarGen *ttbar = new TTbarGen(bcc); 
  double dR = 100;
  if(ttbar->DecayChannel() == TTbarGen::e_ehad || ttbar->DecayChannel() == TTbarGen::e_muhad  ){
    dR = deltaR(bcc->topjets->at(0).v4(), ttbar->TopHad().v4());
  }
  delete ttbar;

  if(dR < m_dR_max)return true;
  return false;
}

std::string TopJetMatchedSelection::description()
{
    char s[100];
    sprintf(s, "deltaR(2nd topjet, top) < %.1f",m_dR_max);
    return s;
}



