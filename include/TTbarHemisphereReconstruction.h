#ifndef TTbarHemisphereReconstruction_H
#define TTbarHemisphereReconstruction_H

#include "SFrameTools/include/BaseCycleContainer.h"
#include "SubstructureAnalysis/include/JetMassTools.h"
#include "SFrameTools/include/Utils.h"

class TTbarHemisphereReconstruction
{

 public:
  //default constructor recostructs leptonic and hadronic hemispheres if possible
  TTbarHemisphereReconstruction(BaseCycleContainer* bcc);

  //default destructor
  ~TTbarHemisphereReconstruction();


  bool ReconstructLeptonicHemisphere();
  bool ReconstructHadronicHemisphere();

  //check if a leptoic or hadronic hemisphere was reconstructed
  bool LeptonicHemisphereFound();
  bool HadronicHemisphereFound();


  //leptonic parameters: only accessable if a leptonic hemisphere was found
  LorentzVector TopLep();

  LorentzVector Neutrino();
  Particle Lepton();

  Jet bJetLep();

  LorentzVector WRecoLep();


  //hadronic parameters: only accessable if a hadronic hemisphere was found
  LorentzVector TopHad();

  Jet WJet1();
  Jet WJet2();

  Jet bJetHad();

  TopJet WJet();

  TopJet CAbJet();

  LorentzVector WRecoHad();

 
  //
  void SetLeptonicBJetPtEtaCSV( double pt, double eta, E_BtagType btagtype);
  void SetHadronicBJetPtEtaCSV( double pt, double eta, E_BtagType btagtype );
  void SetCABJetPtEta( double pt, double eta);
  void SetWJetPtEta(double pt, double eta);

 private:

  // std::vector<Jet> GetBJets( E_BtagType btagtype, double pt_min, double eta_max);
  //std::vector<Jet> GetNOBJets( E_BtagType btagtype, double pt_min, double eta_max)
  //  m_CSV_WP;
  LorentzVector m_TopLep;

  double m_WJetPt;
  double m_WJetEta;

  double m_CAbJetPt;
  double m_CAbJetEta;

  E_BtagType m_bJetLepCSV;
  double m_bJetLepPt;
  double m_bJetLepEta;

  E_BtagType m_bJetHadCSV;
  double m_bJetHadPt;
  double m_bJetHadEta;


  bool m_LeptonicHemisphereFound;
  bool m_HadronicHemisphereFound;

  LorentzVector m_Neutrino;
  Particle m_Lepton;

  Jet m_bJetLep;
  Jet m_bJetHad;

  TopJet m_WJet;
  TopJet m_CAbJet;

  Jet m_WJet1;
  Jet m_WJet2;

  int m_NbJetsLeptonicHemisphere;
  int m_NbJetsHadronicHemisphere;

  BaseCycleContainer* m_bcc;

}; //class TTbarHemisphereReconstruction

#endif// TTbarHemisphereReconstruction_H

