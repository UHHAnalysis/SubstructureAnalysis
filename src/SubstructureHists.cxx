#include "include/SubstructureHists.h"
#include "NtupleWriter/include/JetProps.h"
#include "SFrameTools/include/EventCalc.h"
#include "NtupleWriter/include/GenTopJet.h"
#include "TH2F.h"
#include <iostream>


using namespace std;

 SubstructureHists::SubstructureHists(const char* name) : BaseHists(name)
{
  // named default constructor
   
}

SubstructureHists::~SubstructureHists(){
}

void SubstructureHists::Init()
{
  // book all histograms here
  Book( TH1F("jetmass","m_{jet}",20,0,300));
  Book( TH1F("ptratio","p_{t,2}/p_{t,1}",100,0,1)); 
  Book( TH1F("jetmass_mpu","m_{jet}",20,0,300));
  Book( TH1F("jetmass_lpu","m_{jet}",20,0,300));
  Book( TH2F("median","#rho_{m} : NPV",50,0,50,50,0,50));
  Book( TH2F("rho_npv","#rho_{} : NPV",50,0,50,50,0,50));
  Book( TH1F("pt_leading", "p_{t} leading jet", 100,0,600));
  Book( TH1F("pt_2nd", "p_{t} leading jet", 100,0,600));

  /*
  Book( TH1F("tau1","#tau_{1}",50,0,0.5));
  Book( TH1F("tau1_pruned","#tau_{1} (pruned)",50,0,0.5));
  Book( TH1F("tau2","#tau_{2}",50,0,0.5));
  Book( TH1F("tau2_pruned","#tau_{2} (pruned)",50,0,0.5));
  Book( TH1F("tau3","#tau_{3}",50,0,0.5));
  Book( TH1F("tau3_pruned","#tau_{3} (pruned)",50,0,0.5));
  Book( TH1F("tau32","#tau_{3}/#tau_{2}",50,0,1));
  Book( TH1F("tau32_pruned","#tau_{3}/#tau_{2} (pruned)",50,0,1));
  Book( TH1F("tau21","#tau_{2}/#tau_{1}",50,0,1));
  Book( TH1F("tau21_pruned","#tau_{2}/#tau_{1} (pruned)",50,0,1));
  Book( TH1F("qjets","qjets volatility",50,0,1));

 
  Book( TH1F("t1","#tau_{3} (gen)",50,0,0.5));
  Book( TH1F("t2","#tau_{3} (gen)",50,0,0.5));
  Book( TH1F("t3","#tau_{3} (gen)",50,0,0.5));
  Book( TH1F("t32","#tau_{3}/#tau_{2} (gen)",50,0,1));
  Book( TH1F("t21","#tau_{2}/#tau_{1} (gen)",50,0,1));
  Book( TH1F("qjetsg","qjets volatility (gen)",50,0,1));

  Book( TH2F("migt3t2","#tau_{3}(gen)/#tau_{2}",20,0,1,40,0,1));
  Book( TH2F("migt3","#tau_{3}(gen)",20,0,0.3,40,0,0.3));
  */
}

void SubstructureHists::Fill()
{
  // fill the histograms


  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  bool IsRealData = calc->IsRealData();
  LuminosityHandler* lumih = calc->GetLumiHandler();

  // important: get the event weight
  double weight = calc->GetWeight();

  int run = calc->GetRunNum();
  int lumiblock = calc->GetLumiBlock();
  int Npvs = calc->GetPrimaryVertices()->size();

  ////////////////////////////////////////////////////////////////////
  //calculate rho_m
  ////////////////////////////////////////////////////////////////////
  std::vector<Jet>* Jets = calc->GetJets();
  int Njets =Jets->size();

  std::vector<double> m_a;

  for(int i=0; i<Njets; i++){
    if(fabs(Jets->at(i).eta())<2.5){
      m_a.push_back(Jets->at(i).v4().M()/Jets->at(i).jetArea());
    }
  }

  std::sort (m_a.begin(), m_a.end());
  std::reverse(m_a.begin(),m_a.end()); 

  double rhom=0;

  int Nj = m_a.size();

  if( Nj%2 == 0) rhom=(m_a.at(Nj/2-1)+ m_a.at(Nj/2))/2;
  else rhom=m_a.at((Nj-1)/2) ;

  ////////////////////////////////////////////////////////////////////

  ((TH2F*)Hist("median"))->Fill(Npvs,rhom,weight);
  ((TH2F*)Hist("rho_npv"))->Fill(Npvs,bcc->rho,weight);
  /*
  //Migration matrix
 
  if(calc->GetCAJets() &&calc->GetGenTopJets() ){

     std::vector<GenTopJet>* GenTop_Jets = calc->GetGenTopJets();
    std::vector<TopJet>* Top_Jets = calc->GetCAJets();
     int Ntopjets =Top_Jets->size();
      int Ngentopjets =GenTop_Jets->size();

    for(int i=0; i<Ntopjets && i<Ngentopjets ;i++){
      if(i==0){ 
	TopJet t=Top_Jets->at(i);
	//	std::cout<< "tpt  "<< calc->GetCAJets()->size()<<std::endl;
      	GenTopJet gentop=GenTop_Jets->at(i);

       	JetProps jp(&t, calc->GetPFParticles() );
       	double tau3=jp.GetNsubjettiness(3, Njettiness::onepass_kt_axes, 1., 0.8);
       	double tau2=jp.GetNsubjettiness(2, Njettiness::onepass_kt_axes, 1., 0.8);

       	JetProps jpg(&gentop, calc->GetGenParticles() );
       	double tau3g=jpg.GetNsubjettiness(3, Njettiness::onepass_kt_axes, 1., 0.8);
      	double tau2g=jpg.GetNsubjettiness(2, Njettiness::onepass_kt_axes, 1., 0.8);

       	((TH2F*)Hist("migt3"))->Fill(tau3,tau3g,weight);
     	((TH2F*)Hist("migt3t2"))->Fill(tau3/tau2,tau3g/tau2g,weight);
      }
    }
  

  }
  */
 

  if(calc->GetCAJets()->size()>1){

    std::vector<TopJet>* Top_Jets = calc->GetCAJets();
    int Ntopjets =Top_Jets->size();

    // std::cout<<"size of topjets:  "<<Top_Jets->size()<<std::endl;
    for(int i=0; i<Ntopjets;i++){

      TopJet t=Top_Jets->at(i);
      if(i==0){   
	std::vector<unsigned int> indices = t.pfconstituents_indices();
	  int Nind = indices.size();
	LorentzVector tot;
	std::vector<PFParticle>* pfps = calc->GetPFParticles();
	for( int l=0;l<Nind;l++){
	  tot+= pfps->at(t.pfconstituents_indices().at(l)).v4();
	} 
	//	double jetmass = tot.M();
       	double jetmass = t.v4().M();

	std::vector<unsigned int> indices2 = Top_Jets->at(1).pfconstituents_indices();
	int Nind2 = indices2.size();
	LorentzVector tot2;
	for( int l=0;l<Nind2;l++){
	  tot2+= pfps->at(Top_Jets->at(1).pfconstituents_indices().at(l)).v4();
	} 

	Hist("jetmass")->Fill(jetmass,weight);
	Hist("pt_leading")->Fill(Top_Jets->at(0).pt(),weight);
	Hist("pt_2nd")->Fill(Top_Jets->at(1).pt(),weight);
	if(bcc->rho > 15)	Hist("jetmass_mpu")->Fill(jetmass,weight);
	if(bcc->rho < 15)	Hist("jetmass_lpu")->Fill(jetmass,weight);

	//	if(Ntopjets>1) Hist("ptratio")->Fill(Top_Jets->at(1).pt()/Top_Jets->at(0).pt(),weight);
	if(Ntopjets>1) Hist("ptratio")->Fill(Top_Jets->at(1).pt()/t.pt(),weight);
	//	std::cout<<"pt of RecTopJet:  "<<Top_Jets->at(i).pt()<<"   eta of RecTopJet:  "<<Top_Jets->at(i).eta()<<std::endl;
	//    JetProps jp(&t, calc->GetPFParticles() );
	//  double tau1=jp.GetNsubjettiness(1, Njettiness::onepass_kt_axes, 1., 0.8);
	//	double tau1_pruned=jp.GetPrunedNsubjettiness(1, Njettiness::onepass_kt_axes, 1., 0.8);
	//	double tau2=jp.GetNsubjettiness(2, Njettiness::onepass_kt_axes, 1., 0.8);
	//	double tau2_pruned=jp.GetPrunedNsubjettiness(2, Njettiness::onepass_kt_axes, 1., 0.8);
	//	double tau3=jp.GetNsubjettiness(3, Njettiness::onepass_kt_axes, 1., 0.8);
	//	double tau3_pruned=jp.GetPrunedNsubjettiness(3, Njettiness::onepass_kt_axes, 1., 0.8);
	//	double qjv=jp.GetQjetVolatility(12, 0.8);

	//  if(i==0){
	//	Hist("tau1")->Fill(tau1,weight);
	//	Hist("tau1_pruned")->Fill(tau1_pruned,weight);
	//	Hist("tau2")->Fill(tau2,weight);
	// 	Hist("tau2_pruned")->Fill(tau2_pruned,weight);
	//	Hist("tau3")->Fill(tau3,weight);
	//	Hist("tau3_pruned")->Fill(tau3_pruned,weight);
	//	Hist("tau32")->Fill(tau3/tau2,weight); 
	//	Hist("tau32_pruned")->Fill(tau3_pruned/tau2_pruned,weight);
	//	Hist("tau21")->Fill(tau2/tau1,weight);
	//	Hist("tau21_pruned")->Fill(tau2_pruned/tau1_pruned,weight);
	//	Hist("qjets")->Fill(qjv,weight);

      }
    }
  }
  /*
  if(calc->GetGenTopJets()){
    std::vector<GenTopJet>* GenTop_Jets = calc->GetGenTopJets();
    int Ngentopjets =GenTop_Jets->size();
  //  std::cout<<"size of gentojets:  "<<GenTop_Jets->size()<<std::endl;
    
    for(int i=0;i<Ngentopjets && i<100;i++){
      GenTopJet gentop=GenTop_Jets->at(i);

    //  std::cout<<"pt:  "<<GenTop_Jets->at(i).pt()
	// <<""<<GenTop_Jets->at(i).genparticles_indices().size()
      //     <<std::endl;
    
      JetProps jp(&gentop, calc->GetGenParticles() );
        double tau1=jp.GetNsubjettiness(1, Njettiness::onepass_kt_axes, 1., 0.8);
	double tau2=jp.GetNsubjettiness(2, Njettiness::onepass_kt_axes, 1., 0.8);
	double tau3=jp.GetNsubjettiness(3, Njettiness::onepass_kt_axes, 1., 0.8);
	double qjv=jp.GetQjetVolatility(12, 0.8);

      if(i==0){
	//	std::cout<<"pt of GenTopJet:  "<<GenTop_Jets->at(i).pt()<<"   eta of GenTopJet:  "<<GenTop_Jets->at(i).eta()<<std::endl;
	Hist("t1")->Fill(tau3,weight);

	Hist("t2")->Fill(tau2,weight);
	Hist("t3")->Fill(tau3,weight);
	Hist("t32")->Fill(tau3/tau2,weight); 
	Hist("t21")->Fill(tau2/tau1,weight);
	Hist("qjetsg")->Fill(qjv,weight);
      }
    }
 
  }
  */
 
}

void SubstructureHists::Finish()
{
  // final calculations, like division and addition of certain histograms
  EventCalc* calc = EventCalc::Instance();
  bool IsRealData = calc->IsRealData();
  if (IsRealData){
    //  Hist("N_pv_perLumiBin")->Divide( Hist("N_pv_perLumiBin"), Hist("N_events_perLumiBin"));
    // Hist( "N_pv_perLumiBin")->GetYaxis()->SetTitle("Events/Lumi");
  }

}

