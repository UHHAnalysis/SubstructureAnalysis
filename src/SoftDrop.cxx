#include "include/SoftDrop.h"
#include "NtupleWriter/include/JetProps.h"
#include "include/GenJetProps.h"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include "fastjet/contrib/SoftDrop.hh"
#include "SFrameTools/include/EventCalc.h"

//#include "include/TopTagfunctions.h"
//#include "include/topquark.h"

using namespace std;

SoftDrop* SoftDrop::m_instance = NULL;

SoftDrop* SoftDrop::Instance()
{
  // Get a pointer to the object handler.
  // This is the only way to access this class, 
  // since it's a singleton. This method is accessible
  // from everywhere.

  if (m_instance == NULL){
    m_instance = new SoftDrop();
  }
  return m_instance; 
   
}


SoftDrop::SoftDrop() : m_logger("SoftDrop")
{
  // constructor: initialise all variables

  
  m_logger << DEBUG << "Constructor called." << SLogger::endmsg;
  
     //Reset();

  
}



SoftDrop::~SoftDrop()
{
  // default destructor
  
}


void SoftDrop::GetJet(TopJet topjet,std::vector<PFParticle>* allparts , double zcut, double beta, fastjet::PseudoJet &sd_jet){
  fastjet::ClusterSequence* JetFinder;
  fastjet::JetDefinition* JetDef ;

  JetProps jp(&topjet, allparts);
  std::vector<fastjet::PseudoJet> jetpart = jp.GetJetConstituents();
 
  //Clustering definition
  double conesize=3;
  JetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm,conesize); 

  JetFinder = new fastjet::ClusterSequence(jetpart, *JetDef);

  std::vector<fastjet::PseudoJet> jets = JetFinder->inclusive_jets(10.);
  std::vector<fastjet::PseudoJet> SortedJets = sorted_by_pt(jets);
  
  fastjet::contrib::SoftDrop sd(beta, zcut);

  // cout << "SoftDrop groomer is: " << sd.description() << endl;
  // cout<<jets.size()<<endl;
  if(jets.size()>0) {sd_jet = sd(SortedJets[0]);
    //std::cout<<"SOFT DROP "<<sd_jet.perp()<<std::endl;
  }
}


void SoftDrop::GetJet(GenJetWithParts topjet,std::vector<GenParticle>* allparts , double zcut, double beta, fastjet::PseudoJet &sd_jet){
  fastjet::ClusterSequence* JetFinder;
  fastjet::JetDefinition* JetDef ;

  GenJetProps jp(&topjet, allparts);
  std::vector<fastjet::PseudoJet> jetpart = jp.GetJetConstituents();
 
  //Clustering definition
  double conesize=3;
  JetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm,conesize); 

  JetFinder = new fastjet::ClusterSequence(jetpart, *JetDef);

  std::vector<fastjet::PseudoJet> jets = JetFinder->inclusive_jets(10.);
  std::vector<fastjet::PseudoJet> SortedJets = sorted_by_pt(jets);
  
  fastjet::contrib::SoftDrop sd(beta, zcut);

  // cout << "SoftDrop groomer is: " << sd.description() << endl;
  //cout<<jets.size()<<endl;
  if(jets.size()>0) {sd_jet = sd(SortedJets[0]);
    //std::cout<<"SOFT DROP "<<sd_jet.perp()<<std::endl;
  }
}

void SoftDrop::GetJet(fastjet::PseudoJet topjet, double zcut, double beta, fastjet::PseudoJet &sd_jet){
  fastjet::contrib::SoftDrop sd(beta, zcut);
  sd_jet = sd(topjet);
  

}

