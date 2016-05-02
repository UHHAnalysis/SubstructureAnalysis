// $Id: SubstructureCycle.cxx,v 1.10 2012/12/07 14:21:51 peiffer Exp $

#include <iostream>

using namespace std;
#include "include/SubstructureCycle.h"

#include "SFrameAnalysis/include/SelectionModules.h"
#include "SFrameAnalysis/include/HypothesisHists.h"
#include "SFrameAnalysis/include/TopJetHists.h"
#include "SFrameAnalysis/include/JetHists.h"
#include "SFrameTools/include/TopFitCalc.h"
#include "include/SubstructureSelectionModules.h"

ClassImp( SubstructureCycle );

SubstructureCycle::SubstructureCycle()
  : AnalysisCycle() {

  // constructor, declare additional variables that should be 
  // obtained from the steering-xml file
  
  // set the integrated luminosity per bin for the lumi-yield control plots
  SetIntLumiPerBin(500.);
  m_sys_var = e_Default;
  m_sys_unc = e_None;

  m_mttgencut = false;
  DeclareProperty( "ApplyMttbarGenCut", m_mttgencut );
  DeclareProperty( "Electron_Or_Muon_Selection", m_Electron_Or_Muon_Selection );

  //default: no btagging cuts applied, other cuts can be defined in config file
  m_Nbtags_min=0;
  m_Nbtags_max=int_infinity();
  DeclareProperty( "Nbtags_min", m_Nbtags_min);
  DeclareProperty( "Nbtags_max", m_Nbtags_max);

  // steering property for data-driven qcd in electron channel
  m_reversed_electron_selection = false;
  DeclareProperty( "ReversedElectronSelection", m_reversed_electron_selection);

  // veto electron trigger when running on the JetHT data stream
  m_veto_electron_trigger = false;
  DeclareProperty( "vetoElectronTrigger", m_veto_electron_trigger);

  // put the selected trigger in OR with HLT_PFJet320_v* (electron channel)
  m_useORTriggerWithPFJet320 = false;
  DeclareProperty( "useORTriggerWithPFJet320", m_useORTriggerWithPFJet320);
}

SubstructureCycle::~SubstructureCycle() 
{
  // destructor
}

void SubstructureCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void SubstructureCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void SubstructureCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{
  // declaration of histograms and selections.
  // AnalysisCyle expects Selections and HistCollections to be registered here.
  // Their memory will be released in AnalysisCycle::EndInputData.

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );



}

void SubstructureCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );
}

void SubstructureCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer

  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );
}

void SubstructureCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError ) 
{
  // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results

  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset
  AnalysisCycle::ExecuteEvent( id, weight );

  // get the histogram collections. NOTE: this could be done more performant by making
  // all thse BaseHists* vairables private member variables of SubstructureCycle and
  // setting them in BeginInputData. Then, there is no need here to call GetHistColletion ...
  // BaseHists* HistsNoCutsNoChi2 = GetHistCollection("NoCutsNoChi2");
 
  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();


 
  return;
  
}


