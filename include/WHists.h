#ifndef WHists_H
#define WHists_H

#include "SFrameTools/include/BaseHists.h"
#include "SubstructureAnalysis/include/TTbarHemisphereReconstruction.h"

#include "SFrameTools/include/EventCalc.h"
#include "TH2F.h"

#include "../include/JetMassTools.h"


class WHists : public BaseHists {

public:
   /// Named constructor
   WHists(const char* name);

   /// Default destructor
   ~WHists();

   void Init();

   void Fill();

   void Finish();

private:

   TTbarHemisphereReconstruction *m_tt;

   void FillPtEtaPhiMass(Jet jet, TString postfix, double weight);
   void FillPtEtaPhiMass(TopJet jet, TString postfix, double weight);
   void FillPtEtaPhiMass(LorentzVector vec, TString postfix, double weight);
   void FillPtEtaPhi(Particle particle, TString postfix, double weight);


}; 


#endif // WHists_H
