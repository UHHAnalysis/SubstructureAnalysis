#ifndef SoftDrop_H
#define SoftDrop_H

#include "SFrameTools/include/Objects.h"
#include <fastjet/PseudoJet.hh>
#include "SFrameTools/include/fwd.h"
#include "SFrameTools/include/boost_includes.h" // for shared_array
#include <TMVA/Reader.h>
#include "TVector3.h"
#include <limits>
#include <algorithm>
#include <memory>
#include <TF1.h>
#include "SFrameTools/include/Utils.h"
#include "SFrameTools/include/EventCalc.h"
#include "SFrameAnalysis/include/Cleaner.h"
//#include "FactorizedJetCorrector.h"
//#include "JetCorrectorParameters.h"


//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVectorXYZE;
class SoftDrop{
 private:
  static SoftDrop* m_instance;
  mutable SLogger m_logger;


 public:
SoftDrop();
  ~SoftDrop();
 static SoftDrop* Instance();


 void GetJet(TopJet topjet,std::vector<PFParticle>* allparts, double zcut, double beta,fastjet::PseudoJet & sd_jet);
 void GetJet(GenJetWithParts topjet,std::vector<GenParticle>* allparts, double zcut, double beta,fastjet::PseudoJet & sd_jet);
 void GetJet(fastjet::PseudoJet topjet,double zcut, double beta,fastjet::PseudoJet & sd_jet);

};

#endif
