#include "EventBase.h"
#include "Analysis.h"
#include "HistSvc.h"
#include "Logger.h"
#include "ConfigSvc.h"
#include "SkimSvc.h"
#include "SparseSvc.h"
#include "CutsBase.h"
#include "CutsCountMSSS.h"
#include "CountMSSS.h"

// Constructor
CountMSSS::CountMSSS()
  : Analysis()
{
  // Set up the expected event structure, include branches required for analysis.
  // List the branches required below, see full list here https://luxzeplin.gitlab.io/docs/softwaredocs/analysis/analysiswithrqs/rqlist.html
  // Load in the Single Scatters Branch
  m_event->IncludeBranch("ss");
  m_event->IncludeBranch("ms");
  m_event->IncludeBranch("kr83m");
  m_event->IncludeBranch("pileUp");
  m_event->IncludeBranch("other");
  m_event->IncludeBranch("mcTruthVertices");
  m_event->IncludeBranch("mcTruthEvent");
  m_event->Initialize();
  ////////

  // Setup logging
  logging::set_program_name("CountMSSS Analysis");
  
  // Setup the analysis specific cuts.
  m_cutsCountMSSS = new CutsCountMSSS(m_event);
  
  // Setup config service so it can be used to get config variables
  m_conf = ConfigSvc::Instance();
  //outfile.open("ParentParticle.txt");
  outfile.open("ParentParticle_EvtView.txt");
}

// Destructor
CountMSSS::~CountMSSS() {
  delete m_cutsCountMSSS;
}

/////////////////////////////////////////////////////////////////
// Start of main analysis code block
// 

// Initialize() -  Called once before the event loop.
void CountMSSS::Initialize() {
  INFO("Initializing CountMSSS Analysis");

  EventNumber = 0;
  ParentParticle = "void";


 //SS
  m_hists->BookScatterPlot("SS_LogS1S2_PTFEAlphaN_neutron_p");
  m_hists->BookScatterPlot("SS_LogS1S2_neutron_p");
  m_hists->BookScatterPlot("SS_LogS1S2_USF_p");

  m_hists->BookScatterPlot("SS_FV_LogS1S2_PTFEAlphaN_neutron_p");
  m_hists->BookScatterPlot("SS_FV_LogS1S2_neutron_p");
  m_hists->BookScatterPlot("SS_FV_LogS1S2_USF_p");

 // MS
  m_hists->BookScatterPlot("MS_LogS1S2_PTFEAlphaN_neutron_p");
  m_hists->BookScatterPlot("MS_LogS1S2_neutron_p");
  m_hists->BookScatterPlot("MS_LogS1S2_USF_p");

  m_hists->BookScatterPlot("MS_FV_LogS1S2_PTFEAlphaN_neutron_p");
  m_hists->BookScatterPlot("MS_FV_LogS1S2_neutron_p");
  m_hists->BookScatterPlot("MS_FV_LogS1S2_USF_p");
}


// Execute() - Called once per event.
void CountMSSS::Execute() {

  int nBinsS1 = 500;
  int nBinsR = 5400;
  //int nBinsR = 350;
  int nBinsXY = 800;
  float R2Max = 5400;
  int nBinsZ = 500;
  int s1_min = 0;
  int s1_max = 500;

  int eventID = (*m_event->m_eventHeader)->eventID;
  int runID =(*m_event->m_eventHeader)->runID;
  string rawFileName = (*m_event->m_eventHeader)->rawFileName;

  bool ETrain = false;
  if (EventNumber == ((*m_event->m_TruthEvent)->eventID - 1) && ParentParticle == (*m_event->m_TruthEvent)->parentParticle){
    ETrain = true;
  }
  EventNumber = (*m_event->m_TruthEvent)->eventID;
  ParentParticle = (*m_event->m_TruthEvent)->parentParticle;


  // Print all the Parent Particle
  //if ((*m_event->m_TruthEvent)->parentParticle.find("neutron") != std::string::npos){
  //outfile << "Parent Particle:  " << (*m_event->m_TruthEvent)->parentParticle << endl <<endl;}
  ///
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Ulate_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Ulate_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Ulate_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Uearly_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Uearly_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Uearly_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Po210_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Po210_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Po210_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Th_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Th_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
    if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Th_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
    if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Th_n3", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Th_n4", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_Th_n5", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }

  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" or 
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5" && ETrain == false) 
    { m_hists->BookFillHist("PEn_PTFEAlphaN_neutron_All", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  
  if ((*m_event->m_TruthEvent)->parentParticle == "neutron" && ETrain == false) 
    { m_hists->BookFillHist("PEn_neutron", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  
  if ((*m_event->m_TruthEvent)->parentParticle == "USF_neutron" && ETrain == false) 
    { m_hists->BookFillHist("PEn_USF", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);  }
/*
 for (int p = 0; p < (*m_event->m_TruthPulses)->nRQMCTruthPulses; p++){
    if ((*m_event->m_TruthPulses)->pmtIndices[p][1] < 600){
       if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0") 
          { m_hists->BookFillHist("TPC_PEn_PTFE", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
       if ((*m_event->m_TruthEvent)->parentParticle == "neutron") 
          { m_hists->BookFillHist("TPC_PEn_neutron", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
       if ((*m_event->m_TruthEvent)->parentParticle == "USF_neutron") 
          { m_hists->BookFillHist("TPC_PEn_USF", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);  }

    }
  }
*/
  // Make a variable from the contents of the data
  //  prevents having to write full ' (*event->Branch)->variableName ' structure many times in the analysis
  int numSS = (*m_event->m_singleScatter)->nSingleScatters;
  int numMS = (*m_event->m_multipleScatter)->nMultipleScatters;
  int numKr83m = (*m_event->m_kr83mScatter)->nKr83mScatters;
  int numpileUp = (*m_event->m_pileupScatter)->nPileUpScatters;
  int numother = (*m_event->m_otherScatter)->nOtherScatters;
/*
  if (numKr83m > 0){
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0") 
        { m_hists->BookFillHist("Kr83m_PEn_PTFE", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      
      if ((*m_event->m_TruthEvent)->parentParticle == "neutron") 
        { m_hists->BookFillHist("Kr83m_PEn_neutron", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      
      if ((*m_event->m_TruthEvent)->parentParticle == "USF_neutron") 
        { m_hists->BookFillHist("Kr83m_PEn_USF", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  }
*/
  if (numpileUp > 0){
      
  if ((*m_event->m_TruthEvent)->parentParticle == "neutron" && ETrain == false) 
      { m_hists->BookFillHist("pileUp_PEn_neutron", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }

  if ((*m_event->m_TruthEvent)->parentParticle == "USF_neutron" && ETrain == false) 
      { m_hists->BookFillHist("pileUp_PEn_USF", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }


  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" or 
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_All", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }

  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Ulate_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Ulate_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Ulate_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Uearly_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Uearly_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Uearly_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Po210_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Po210_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Po210_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Th_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Th_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Th_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Th_n3", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Th_n4", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5" && ETrain == false) 
    { m_hists->BookFillHist("pileUp_PEn_PTFEAlphaN_neutron_Th_n5", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }

  }

  if (numother > 0){
      
  if ((*m_event->m_TruthEvent)->parentParticle == "neutron" && ETrain == false) 
        { m_hists->BookFillHist("other_PEn_neutron", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      
  if ((*m_event->m_TruthEvent)->parentParticle == "USF_neutron" && ETrain == false) 
        { m_hists->BookFillHist("other_PEn_USF", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }

  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" or 
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_All", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }    

  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Ulate_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Ulate_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Ulate_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Uearly_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Uearly_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Uearly_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Po210_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Po210_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Po210_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Th_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Th_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Th_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Th_n3", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Th_n4", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5" && ETrain == false) 
    { m_hists->BookFillHist("other_PEn_PTFEAlphaN_neutron_Th_n5", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }

  }

  // if there is a single scatter in the event then plot the S1 pulse area
  // SS 
  Double_t S2c_ss =log10((*m_event->m_singleScatter)->correctedS2Area_phd);
  Double_t S1c_ss = (*m_event->m_singleScatter)->correctedS1Area_phd;
  Double_t Xc_ss = (*m_event->m_singleScatter)->correctedX_cm;
  Double_t Yc_ss = (*m_event->m_singleScatter)->correctedY_cm;
  Double_t Zc_ss = (*m_event->m_singleScatter)->correctedZ_cm;
  double Rc_ss = sqrt( Xc_ss*Xc_ss + Yc_ss*Yc_ss );
  double R2c_ss = ( Xc_ss*Xc_ss + Yc_ss*Yc_ss );
       
  if (numSS > 0){
      
  if ((*m_event->m_TruthEvent)->parentParticle == "neutron" && ETrain == false) 
        { m_hists->BookFillHist("SS_PEn_neutron", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      
  if ((*m_event->m_TruthEvent)->parentParticle == "USF_neutron" && ETrain == false) 
        { m_hists->BookFillHist("SS_PEn_USF", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }

  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" or 
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" or
      (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_All", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }     

  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Ulate_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Ulate_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Ulate_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Uearly_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Uearly_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Uearly_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Po210_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Po210_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Po210_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Th_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Th_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Th_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Th_n3", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Th_n4", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
  if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5" && ETrain == false) 
    { m_hists->BookFillHist("SS_PEn_PTFEAlphaN_neutron_Th_n5", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }

  }

// SS S1 area < 500 
  if (numSS > 0 && S1c_ss < 500.){

      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Ulate_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Ulate_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Ulate_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Uearly_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Uearly_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Uearly_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Po210_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Po210_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Po210_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Th_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Th_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Th_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Th_n3", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Th_n4", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5") 
        { m_hists->BookFillHist("SS500_PEn_PTFEAlphaN_neutron_Th_n5", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }


      if ( (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" or 
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5") 
          {m_hists->BookFillHist("SS_PEn_PTFE_S1_500", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); 
           m_hists->BookFillHist("SS_LogS1S2_PTFEAlphaNneutron",nBinsS1,s1_min,s1_max, nBinsS1, 2, 7, S1c_ss, S2c_ss);
           m_hists->GetScatterPlot("SS_LogS1S2_PTFEAlphaN_neutron_p")->Fill( S1c_ss, S2c_ss);
           m_hists->BookFillHist("SS_R2vsDrift_PTFEAlphaNneutron", nBinsR,0.,R2Max,nBinsZ,0.,1000.,R2c_ss,Zc_ss);
           m_hists->BookFillHist("SS_XY_PTFEAlphaNneutron",nBinsXY,-80,80,nBinsXY,-80, 80, Xc_ss, Yc_ss);
           
           outfile << "SS(S1500): PTFEAlphaN_neutron_All rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ss << "," << S2c_ss << endl << endl;
           std::cout << "SS(S1500): PTFEAlphaN_neutron_All rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ss << "," << S2c_ss << endl << endl; }
      
      if ((*m_event->m_TruthEvent)->parentParticle == "neutron") 
         {m_hists->BookFillHist("SS_PEn_neutron_S1_500", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);
          m_hists->BookFillHist("SS_LogS1S2_neutron",nBinsS1,s1_min,s1_max, nBinsS1, 2, 7, S1c_ss, S2c_ss);
          m_hists->GetScatterPlot("SS_LogS1S2_neutron_p")->Fill( S1c_ss, S2c_ss); 
          m_hists->BookFillHist("SS_R2vsDrift_neutron", nBinsR,0.,R2Max,nBinsZ,0.,1000.,R2c_ss,Zc_ss);
          m_hists->BookFillHist("SS_XY_neutron",nBinsXY,-80,80,nBinsXY,-80,80,Xc_ss,Yc_ss);
          
          outfile << "SS(S1500): neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ss << "," << S2c_ss << endl << endl;
          std::cout << "SS(S1500): neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ss << "," << S2c_ss << endl << endl;}
      
      if ((*m_event->m_TruthEvent)->parentParticle == "USF_neutron") 
         {m_hists->BookFillHist("SS_PEn_USF_S1_500", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);
          m_hists->BookFillHist("SS_LogS1S2_USF",nBinsS1,s1_min,s1_max, nBinsS1, 2, 7, S1c_ss, S2c_ss);
          m_hists->GetScatterPlot("SS_LogS1S2_USF_p")->Fill( S1c_ss, S2c_ss); 
          m_hists->BookFillHist("SS_R2vsDrift_USF", nBinsR,0.,R2Max,nBinsZ,0.,1000.,R2c_ss,Zc_ss);
          m_hists->BookFillHist("SS_XY_USF",nBinsXY,-80,80,nBinsXY,-80,80, Xc_ss, Yc_ss);

          outfile << "SS(S1500): USF_neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ss << "," << S2c_ss << endl << endl;
          std::cout << "SS(S1500): USF_neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ss << "," << S2c_ss << endl << endl;}
  }

// SS, S1 area < 500, FV (Rc_ss < 68.8 && Zc_ss > 1.5 && Zc_ss < 132.1)
  if (numSS > 0 && S1c_ss < 500. && Rc_ss < 68.8 && Zc_ss > 1.5 && Zc_ss < 132.1){

      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Ulate_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Ulate_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Ulate_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Uearly_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Uearly_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Uearly_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Po210_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Po210_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Po210_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Th_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Th_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Th_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Th_n3", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Th_n4", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5") 
        { m_hists->BookFillHist("SS500_FV_PEn_PTFEAlphaN_neutron_Th_n5", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }


      if ( (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" or 
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5") 
         {m_hists->BookFillHist("SS_FV_PEn_PTFE_S1_500", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);
          m_hists->BookFillHist("SS_FV_LogS1S2_PTFEAlphaNneutron",nBinsS1,s1_min,s1_max, nBinsS1, 2, 7, S1c_ss, S2c_ss); 
          m_hists->GetScatterPlot("SS_FV_LogS1S2_PTFEAlphaN_neutron_p")->Fill( S1c_ss, S2c_ss);
          m_hists->BookFillHist("SS_FV_R2vsDrift_PTFEAlphaNneutron", nBinsR,0.,R2Max,nBinsZ,0.,1000.,R2c_ss,Zc_ss);
          m_hists->BookFillHist("SS_FV_XY_PTFEAlphaNneutron",nBinsXY,-80,80,nBinsXY,-80, 80, Xc_ss, Yc_ss);

         outfile << "SS_FV(S1500): PTFEAlphaN_neutron_All rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ss << "," << S2c_ss << endl << endl;
         std::cout << "SS_FV(S1500): PTFEAlphaN_neutron_All rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ss << "," << S2c_ss << endl << endl;}
      
      if ((*m_event->m_TruthEvent)->parentParticle == "neutron") 
         {m_hists->BookFillHist("SS_FV_PEn_neutron_S1_500", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);
          m_hists->BookFillHist("SS_FV_LogS1S2_neutron",nBinsS1,s1_min,s1_max, nBinsS1, 2, 7, S1c_ss, S2c_ss);
          m_hists->GetScatterPlot("SS_FV_LogS1S2_neutron_p")->Fill( S1c_ss, S2c_ss);
          m_hists->BookFillHist("SS_FV_R2vsDrift_neutron", nBinsR,0.,R2Max,nBinsZ,0.,1000.,R2c_ss,Zc_ss);
          m_hists->BookFillHist("SS_FV_XY_neutron",nBinsXY,-80,80,nBinsXY,-80,80,Xc_ss,Yc_ss); 

         outfile << "SS_FV(S1500): neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ss << "," << S2c_ss << endl << endl;
         std::cout << "SS_FV(S1500): neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ss << "," << S2c_ss << endl << endl;}
      
      if ((*m_event->m_TruthEvent)->parentParticle == "USF_neutron") 
         {m_hists->BookFillHist("SS_FV_PEn_USF_S1_500", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);
          m_hists->BookFillHist("SS_FV_LogS1S2_USF",nBinsS1,s1_min,s1_max, nBinsS1, 2, 7, S1c_ss, S2c_ss);
          m_hists->GetScatterPlot("SS_FV_LogS1S2_USF_p")->Fill( S1c_ss, S2c_ss);
          m_hists->BookFillHist("SS_FV_R2vsDrift_USF", nBinsR,0.,R2Max,nBinsZ,0.,1000.,R2c_ss,Zc_ss);
          m_hists->BookFillHist("SS_FV_XY_USF",nBinsXY,-80,80,nBinsXY,-80,80, Xc_ss, Yc_ss); 

          outfile << "SS_FV(S1500): USF_neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ss << "," << S2c_ss << endl << endl;
          std::cout << "SS_FV(S1500): USF_neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ss << "," << S2c_ss << endl << endl;}
  }

//////////////////////////////
//
// Multiple Scatters ////////
//
/////////////////////////////  

// For MS
      double S2c_ms = 0.; double S1c_ms = 0.; 
      double Xc_ms = 0.; double Yc_ms = 0.; double Zc_ms = 0.;
      double largestS2 = 0.;
      int S2index;
      for (int i = 0; i < ((*m_event->m_multipleScatter)->nS2s); i++ ){
            ////cout << "  in nS2s MS loop... ";
            S2c_ms += ((*m_event->m_multipleScatter)->correctedS2Area_phd)[i];  //add up the S2 areas
            S1c_ms += ((*m_event->m_multipleScatter)->correctedS1Areas_phd)[i] * ((*m_event->m_multipleScatter)->correctedS2Area_phd)[i] ;   //add each correctedS1, weighted by its S2
            Xc_ms += ((*m_event->m_multipleScatter)->correctedX_cm)[i] * ((*m_event->m_multipleScatter)->correctedS2Area_phd)[i];  //add each X-coord., weighted by its S2
            Yc_ms += ((*m_event->m_multipleScatter)->correctedY_cm)[i] * ((*m_event->m_multipleScatter)->correctedS2Area_phd)[i];
            Zc_ms += ((*m_event->m_multipleScatter)->correctedZ_cm)[i] * ((*m_event->m_multipleScatter)->correctedS2Area_phd)[i];
            //double z = ((*m_event->m_multipleScatter)->Zc_ss);  // same for Z
            if ( ((*m_event->m_multipleScatter)->correctedS2Area_phd)[i] > largestS2 ) {
               largestS2 = ((*m_event->m_multipleScatter)->correctedS2Area_phd)[i];    //save the largest S2, just in case
               S2index = i; 
            }
      }
          
          
      S1c_ms /= S2c_ms; //area-weighted average corrected S1
      Xc_ms /= S2c_ms; Yc_ms /= S2c_ms; Zc_ms /= S2c_ms;  //weighted X and Y
      double logS2c_ms = log10( S2c_ms );
      double Rc_ms = sqrt( Xc_ms*Xc_ms + Yc_ms*Yc_ms );
      double R2c_ms = ( Xc_ms*Xc_ms + Yc_ms*Yc_ms );


    if (numMS > 0){
        if ( (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" or 
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5") 
         {m_hists->BookFillHist("MS_PEn_PTFE_All", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }

      if ((*m_event->m_TruthEvent)->parentParticle == "neutron") 
         {m_hists->BookFillHist("MS_PEn_neutron", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }

      if ((*m_event->m_TruthEvent)->parentParticle == "USF_neutron") 
         {m_hists->BookFillHist("MS_PEn_USF", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);}
    //m_hists->BookFillHist("SingleScatters", 300, 0., 3000., (*m_event->m_singleScatter)->s1Area_phd);

      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Ulate_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Ulate_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Ulate_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Uearly_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Uearly_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Uearly_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Po210_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Po210_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Po210_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Th_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Th_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Th_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Th_n3", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Th_n4", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5") 
        { m_hists->BookFillHist("MS_PEn_PTFEAlphaN_neutron_Th_n5", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }

  }

  // MS S1 area < 500 
  if (numMS > 0 && S1c_ms < 500.){


      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Ulate_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Ulate_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Ulate_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Uearly_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Uearly_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Uearly_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Po210_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Po210_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Po210_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Th_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Th_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Th_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Th_n3", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Th_n4", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5") 
        { m_hists->BookFillHist("MS500_PEn_PTFEAlphaN_neutron_Th_n5", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }

      if ( (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" or 
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5") 
        { m_hists->BookFillHist("MS_PEn_PTFE_S1_500", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);
          m_hists->BookFillHist("MS_LogS1S2_PTFEAlphaNneutron",nBinsS1,s1_min,s1_max, nBinsS1, 2, 7, S1c_ms, logS2c_ms); 
          m_hists->GetScatterPlot("MS_LogS1S2_PTFEAlphaN_neutron_p")->Fill( S1c_ms, logS2c_ms);
          m_hists->BookFillHist("MS_R2vsDrift_PTFEAlphaNneutron", nBinsR,0.,R2Max,nBinsZ,0.,1000.,R2c_ms,Zc_ms);
          m_hists->BookFillHist("MS_XY_PTFEAlphaNneutron",nBinsXY,-80,80,nBinsXY,-80, 80, Xc_ms, Yc_ms);

          outfile << "MS(S1500): PTFEAlphaN_neutron_All rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ms << "," << logS2c_ms << endl << endl;
          std::cout << "MS(S1500): PTFEAlphaN_neutron_All rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ms << "," << logS2c_ms << endl << endl;}
      
      if ((*m_event->m_TruthEvent)->parentParticle == "neutron") 
        { m_hists->BookFillHist("MS_PEn_neutron_S1_500", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);
          m_hists->BookFillHist("MS_LogS1S2_neutron",nBinsS1,s1_min,s1_max, nBinsS1, 2, 7, S1c_ms, logS2c_ms);
          m_hists->GetScatterPlot("MS_LogS1S2_neutron_p")->Fill( S1c_ms, logS2c_ms); 
          m_hists->BookFillHist("MS_R2vsDrift_neutron", nBinsR,0.,R2Max,nBinsZ,0.,1000.,R2c_ms,Zc_ms);
          m_hists->BookFillHist("MS_XY_neutron",nBinsXY,-80,80,nBinsXY,-80,80,Xc_ms,Yc_ms);

         outfile << "MS(S1500): neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ms << "," << logS2c_ms << endl << endl;
         std::cout << "MS(S1500): neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ms << "," << logS2c_ms << endl << endl;
       }
      
      if ((*m_event->m_TruthEvent)->parentParticle == "USF_neutron") 
        { m_hists->BookFillHist("MS_PEn_USF_S1_500", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);
          m_hists->BookFillHist("MS_LogS1S2_USF",nBinsS1,s1_min,s1_max, nBinsS1, 2, 7, S1c_ms, logS2c_ms);
          m_hists->GetScatterPlot("MS_LogS1S2_USF_p")->Fill( S1c_ms, logS2c_ms); 
          m_hists->BookFillHist("MS_R2vsDrift_USF", nBinsR,0.,R2Max,nBinsZ,0.,1000.,R2c_ms,Zc_ms);
          m_hists->BookFillHist("MS_XY_USF",nBinsXY,-80,80,nBinsXY,-80,80, Xc_ms, Yc_ms);

          outfile << "MS(S1500): USF_neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ms << "," << logS2c_ms << endl << endl;
          std::cout << "MS(S1500): USF_neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ms << "," << logS2c_ms << endl << endl;}
        }

  // MS, S1 area < 500, FV (Rc_ms < 68.8 && Zc_ms > 1.5 && Zc_ms < 132.1)
  if (numMS > 0 && S1c_ms < 500. && Rc_ms < 68.8 && Zc_ms > 1.5 && Zc_ms < 132.1){


      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Ulate_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Ulate_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Ulate_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Uearly_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Uearly_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Uearly_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Po210_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Po210_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Po210_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Th_n0", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Th_n1", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Th_n2", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Th_n3", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Th_n4", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }
      if ((*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5") 
        { m_hists->BookFillHist("MS500_FV_PEn_PTFEAlphaN_neutron_Th_n5", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV); }


      if ( (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n1" or 
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Ulate_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Uearly_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Po210_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n0" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n1" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n2" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n3" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n4" or
            (*m_event->m_TruthEvent)->parentParticle == "PTFEAlphaN_neutron_Th_n5")  
        { m_hists->BookFillHist("MS_FV_PEn_PTFE_S1_500", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);
          m_hists->BookFillHist("MS_FV_LogS1S2_PTFEAlphaNneutron",nBinsS1,s1_min,s1_max, nBinsS1, 2, 7, S1c_ms, logS2c_ms); 
          m_hists->GetScatterPlot("MS_FV_LogS1S2_PTFEAlphaN_neutron_p")->Fill( S1c_ms, logS2c_ms);
          m_hists->BookFillHist("MS_FV_R2vsDrift_PTFEAlphaNneutron", nBinsR,0.,R2Max,nBinsZ,0.,1000.,R2c_ms,Zc_ms);
          m_hists->BookFillHist("MS_FV_XY_PTFEAlphaNneutron",nBinsXY,-80,80,nBinsXY,-80, 80, Xc_ms, Yc_ms);

          outfile << "MS_FV(S1500): PTFEAlphaN_neutron_All rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ms << "," << logS2c_ms << endl << endl;
          std::cout << "MS_FV(S1500): PTFEAlphaN_neutron_All rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ms << "," << logS2c_ms << endl << endl;}
      
      if ((*m_event->m_TruthEvent)->parentParticle == "neutron") 
        { m_hists->BookFillHist("MS_FV_PEn_neutron_S1_500", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);
          m_hists->BookFillHist("MS_FV_LogS1S2_neutron",nBinsS1,s1_min,s1_max, nBinsS1, 2, 7, S1c_ms, logS2c_ms);
          m_hists->GetScatterPlot("MS_FV_LogS1S2_neutron_p")->Fill( S1c_ms, logS2c_ms); 
          m_hists->BookFillHist("MS_FV_R2vsDrift_neutron", nBinsR,0.,R2Max,nBinsZ,0.,1000.,R2c_ms,Zc_ms);
          m_hists->BookFillHist("MS_FV_XY_neutron",nBinsXY,-80,80,nBinsXY,-80,80,Xc_ms,Yc_ms);

          outfile << "MS_FV(S1500): neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ms << "," << logS2c_ms << endl << endl;
          std::cout << "MS_FV(S1500): neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ms << "," << logS2c_ms << endl << endl;
        }
      
      if ((*m_event->m_TruthEvent)->parentParticle == "USF_neutron") 
        { m_hists->BookFillHist("MS_FV_PEn_USF_S1_500", 100, 0., 10000., (*m_event->m_TruthEvent)->parentEnergy_keV);
          m_hists->BookFillHist("MS_FV_LogS1S2_USF",nBinsS1,s1_min,s1_max, nBinsS1, 2, 7, S1c_ms, logS2c_ms);
          m_hists->GetScatterPlot("MS_FV_LogS1S2_USF_p")->Fill( S1c_ms, logS2c_ms); 
          m_hists->BookFillHist("MS_FV_R2vsDrift_USF", nBinsR,0.,R2Max,nBinsZ,0.,1000.,R2c_ms,Zc_ms);
          m_hists->BookFillHist("MS_FV_XY_USF",nBinsXY,-80,80,nBinsXY,-80,80, Xc_ms, Yc_ms);

          outfile << "MS_FV(S1500): USF_neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ms << "," << logS2c_ms << endl << endl;
          std::cout << "MS_FV(S1500): USF_neutron rawFile:  "<< rawFileName<< "," << runID << "," << eventID<< "," << S1c_ms << "," << logS2c_ms << endl << endl;}
        }

}

// Finalize() - Called once after event loop.
void CountMSSS::Finalize() {
  INFO("Finalizing CountMSSS Analysis");

 //SS
  m_hists->DrawERNRBands("SS_LogS1S2_PTFEAlphaNneutron", "log(s2)");
  m_hists->DrawERNRBands("SS_LogS1S2_neutron", "log(s2)");
  m_hists->DrawERNRBands("SS_LogS1S2_USF", "log(s2)");
  m_hists->DrawERNRBands("SS_FV_LogS1S2_PTFEAlphaNneutron", "log(s2)");
  m_hists->DrawERNRBands("SS_FV_LogS1S2_neutron", "log(s2)");
  m_hists->DrawERNRBands("SS_FV_LogS1S2_USF", "log(s2)");

  m_hists->DrawERNRBandsScatterPlot("SS_LogS1S2_PTFEAlphaN_neutron_p","log(s2)");
  m_hists->DrawERNRBandsScatterPlot("SS_LogS1S2_neutron_p","log(s2)");
  m_hists->DrawERNRBandsScatterPlot("SS_LogS1S2_USF_p","log(s2)");

  m_hists->DrawERNRBandsScatterPlot("SS_FV_LogS1S2_PTFEAlphaN_neutron_p","log(s2)");
  m_hists->DrawERNRBandsScatterPlot("SS_FV_LogS1S2_neutron_p","log(s2)");
  m_hists->DrawERNRBandsScatterPlot("SS_FV_LogS1S2_USF_p","log(s2)");


  // MS

  m_hists->DrawERNRBands("MS_LogS1S2_PTFEAlphaNneutron", "log(s2)");
  m_hists->DrawERNRBands("MS_LogS1S2_neutron", "log(s2)");
  m_hists->DrawERNRBands("MS_LogS1S2_USF", "log(s2)");
  m_hists->DrawERNRBands("MS_FV_LogS1S2_PTFEAlphaNneutron", "log(s2)");
  m_hists->DrawERNRBands("MS_FV_LogS1S2_neutron", "log(s2)");
  m_hists->DrawERNRBands("MS_FV_LogS1S2_USF", "log(s2)");

  m_hists->DrawERNRBandsScatterPlot("MS_LogS1S2_PTFEAlphaN_neutron_p","log(s2)");
  m_hists->DrawERNRBandsScatterPlot("MS_LogS1S2_neutron_p","log(s2)");
  m_hists->DrawERNRBandsScatterPlot("MS_LogS1S2_USF_p","log(s2)");

  m_hists->DrawERNRBandsScatterPlot("MS_FV_LogS1S2_PTFEAlphaN_neutron_p","log(s2)");
  m_hists->DrawERNRBandsScatterPlot("MS_FV_LogS1S2_neutron_p","log(s2)");
  m_hists->DrawERNRBandsScatterPlot("MS_FV_LogS1S2_USF_p","log(s2)");

  outfile.close();
}

