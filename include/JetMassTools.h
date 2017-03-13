#ifndef JetMassTools_H
#define JetMassTools_H

#include <iostream>
#include <vector>
#include "NtupleWriter/include/JetProps.h"
#include "../../SFrameTools/include/Objects.h"
#include "../../SFrameTools/include/BaseCycleContainer.h"
#include "SFrameTools/include/EventCalc.h"
#include "SoftDrop.h"


using namespace std;

struct HigherMass {
    bool operator() (const Particle& j1, const Particle& j2) const {
      return j1.v4().mass() > j2.v4().mass();
    };
};

class TTbar_gen_properties{

 public:

  TTbar_gen_properties();
  ~TTbar_gen_properties(){};

  void evaluate();
  bool isMuon(){return m_isMuon;};
  bool isEle(){return m_isEle;};
  const GenParticle* GetLepton(){return m_lepton;};
  GenParticle GetHadronicTop(){return m_hadronic_top;};


 private:
  BaseCycleContainer* bcc;
  bool m_isMuon;
  bool m_isEle;
  const GenParticle* m_lepton;
  GenParticle m_hadronic_top;
}; 


class GenCleaner{

 public:
  GenCleaner();
  GenCleaner(BaseCycleContainer*);
  ~GenCleaner(){}; 
  void GenCaJetLeptonSubstractor(bool sort);
  void CAGenJetCleaner(double ptmin, double etamax);

 private:
  EventCalc* calc;
  BaseCycleContainer* bcc;
  void resetEventCalc();
 
};

int Wjet_index(std::vector<Jet> *jets , std::vector<TopJet> *topjets, LorentzVector lepton, bool massdiff = false);

fastjet::PseudoJet FilterJet(fastjet::PseudoJet injet, std::vector<PFParticle>* pfparticles);

  void FilterNTopJets(unsigned int Njets);

void SoftDropNTopJets(unsigned int Njets);
void SoftDropNGenJetsWithParts(unsigned int Njets);
void TrimmNTopJets(unsigned int Njets);
void ScaleJetMass(double scale);
std::vector<Jet> GetBJets( std::vector<Jet> jets, E_BtagType btagtype, double pt_min, double eta_max);
std::vector<Jet> GetNOBJets( std::vector<Jet> jets, E_BtagType btagtype, double pt_min, double eta_max);




#endif
