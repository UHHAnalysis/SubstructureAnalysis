//--sframe new--
#ifndef SubstructureSelectionModules_H
#define SubstructureSelectionModules_H

#include "SFrameTools/include/Objects.h"
#include "SFrameTools/include/Utils.h"
#include "SFrameTools/include/BaseCycleContainer.h"
#include "SFrameTools/include/Selection.h"
#include "SFrameTools/include/EventCalc.h"
#include "SFrameTools/include/HypothesisDiscriminator.h"
#include "SFrameTools/include/TTbarGen.h"
#include "SFrameAnalysis/include/EventFilterFromListStandAlone.h"
#include "SubstructureAnalysis/include/JetMassTools.h"
#include "SubstructureAnalysis/include/TTbarHemisphereReconstruction.h"

/*
#include "SFrameAnalysis/include/CMSTopTagSelectionMods.h"
#include "SFrameAnalysis/include/JetSelectionMods.h"
#include "SFrameAnalysis/include/LeptonSelectionMods.h"
#include "SFrameAnalysis/include/HepSelectionMods.h"
#include "SFrameAnalysis/include/BTagSelectionMods.h"
#include "SFrameAnalysis/include/HTSelectionMods.h"
*/
#include <algorithm>
#include <memory>


//------------------------------------
// The first part of the name of the selection decides in which file it is stored eg:
// MuonBTagSelection --> LeptonSelectionMods
// HEPMuonSelection  --> HepSelectionMods
//-----------------------------------------
class SubstructureGeneratorPreSelection: public SelectionModule {
 public:
  SubstructureGeneratorPreSelection(double pt_hadronic_top);
  ~SubstructureGeneratorPreSelection() {};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_pt_hadonic_top;
};


class SubstructureGeneratorPreSelectionMuon: public SelectionModule {
 public:
  SubstructureGeneratorPreSelectionMuon(double pt_hadronic_top);
  ~SubstructureGeneratorPreSelectionMuon() {};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_pt_hadonic_top;
};


class SubstructureGeneratorPreSelectionElectron: public SelectionModule {
 public:
  SubstructureGeneratorPreSelectionElectron(double pt_hadronic_top);
  ~SubstructureGeneratorPreSelectionElectron() {};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_pt_hadonic_top;
};

class SubstructureGeneratorPreSelectionHadron: public SelectionModule {
 public:
  SubstructureGeneratorPreSelectionHadron();
  ~SubstructureGeneratorPreSelectionHadron() {};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_pt_hadonic_top;
};

class GenLeptonSelection: public SelectionModule {
 public:
  GenLeptonSelection(double pt_lep, double eta_max);
  ~GenLeptonSelection() {};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_pt_lep;
  double m_eta_max;
};


class NCAGenJetSelection: public SelectionModule {
 public:
  NCAGenJetSelection(int min_nparticle, int max_nparticle=int_infinity(), double ptmin=0., double etamax=double_infinity());
  ~NCAGenJetSelection() {};
  
  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();
  
 private:
  int m_min_nparticle;
  int m_max_nparticle;
  double m_ptmin;
  double m_etamax;
};


class NCAGenJetveto: public SelectionModule {
 public:
  NCAGenJetveto( double ptmin, double njets, double etamax);
  ~NCAGenJetveto() {};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_ptmin;
  double m_njets;
  double m_etamax;
};


class NTopJetveto: public SelectionModule {
 public:
  NTopJetveto( double ptmin, double njets, double etamax);
  ~NTopJetveto() {};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_ptmin;
  double m_njets;
  double m_etamax;
};


class TTbarDecayChannelSelectionGen: public SelectionModule {
 public:
  TTbarDecayChannelSelectionGen(std::string DecayChannel);
  ~TTbarDecayChannelSelectionGen() {};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  std::string m_DecayChannel;
};


class LeptonicHemisphereSelection: public SelectionModule {
 public:
  LeptonicHemisphereSelection(double jetpt = 30, double abseta = 2.1, E_BtagType btagtype = e_CSVM, double dR = 1.2);
  ~LeptonicHemisphereSelection() {};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_jetpt;
  double m_abseta;
  double m_dR;
  E_BtagType m_type;
};


class BJetCloseToLeptonSelection: public SelectionModule {
 public:
  BJetCloseToLeptonSelection(double jetpt, double abseta, E_BtagType btagtype, double dR);
  ~BJetCloseToLeptonSelection() {};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_jetpt;
  double m_abseta;
  double m_dR;
  E_BtagType m_type;
};


class HadronicHemisphereSelection: public SelectionModule {
 public:
  HadronicHemisphereSelection(double jetpt = 30, double abseta = 2.1, E_BtagType btagtype = e_CSVL, double dR = 1.2);
  ~HadronicHemisphereSelection() {};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_jetpt;
  double m_abseta;
  double m_dR;
  E_BtagType m_type;
};


class LeadingTopJetHigherMassSelection: public SelectionModule {
 public:
  LeadingTopJetHigherMassSelection( TString mode = "Default"); 
  ~LeadingTopJetHigherMassSelection(){};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  TString m_mode;
};

class LeadingGenCAJetHigherMassSelection: public SelectionModule {
 public:
  LeadingGenCAJetHigherMassSelection( TString mode = "Default"); 
  ~LeadingGenCAJetHigherMassSelection(){};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  TString m_mode;
};

class DR2ndGenCAJetLeptonSelection: public SelectionModule {
 public:
  DR2ndGenCAJetLeptonSelection( double dR_max); 
  ~DR2ndGenCAJetLeptonSelection(){};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_dR_max;
};


class DR2ndTopJetLeptonSelection: public SelectionModule {
 public:
  DR2ndTopJetLeptonSelection( double dR_max); 
  ~DR2ndTopJetLeptonSelection(){};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_dR_max;
};



class GenCAJetMatchedSelection: public SelectionModule {
 public:
  GenCAJetMatchedSelection( double dR_max, TString channel); 
  ~GenCAJetMatchedSelection(){};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_dR_max;
  TString m_channel;
};

class TopJetMatchedSelection: public SelectionModule {
 public:
  TopJetMatchedSelection( double dR_max); 
  ~TopJetMatchedSelection(){};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();

 private:
  double m_dR_max;
};


#endif
