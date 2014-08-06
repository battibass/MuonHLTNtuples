#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TTree.h"
#include "TBranch.h"

#include "MuTree.h"
#include "tdrstyle.C"

#include <iostream>
#include <algorithm>
#include <sstream>
#include <vector>
#include <map>

using namespace ciemat;

// BASE cuts and definitions
// to configute the analysis
const int   PT_THR       = 24;
const int   MAX_ETA      = 2.4;
const float ETA_BINS[13] = {-2.4, -2.1, -1.6, -1.2, -0.8, -0.4, 0., 0.4, 0.8, 1.2, 1.6, 2.1, 2.4};

class RatePlotter 
{

public:

  RatePlotter(float ptThr, float maxEta, std::string filterTag);
  ~RatePlotter();

  void book(TFile* fileOut);
  void fill(ciemat::Event *ev);

  void plotAndSave(float nEntries);

private:

  float m_ptThr;
  float m_maxEta;
  std::string m_filterTag;

  std::map<std::string, TH1F*> m_histoMap;

};

RatePlotter::RatePlotter(float ptThr, float maxEta, std::string filterTag) :
                         m_ptThr(ptThr), m_maxEta(maxEta), m_filterTag(filterTag)
{ 

}

RatePlotter::~RatePlotter()
{

  m_histoMap.clear();

}

void RatePlotter::book(TFile* fileOut)
{

  std::stringstream buildPlotterTag;
  buildPlotterTag << m_filterTag << "_Pt" << m_ptThr << "_Eta" << int(m_maxEta*10); // CB do better here 
  
  std::string plotterTag = buildPlotterTag.str();

  fileOut->cd("/");
  fileOut->mkdir(plotterTag.c_str());
  fileOut->cd(plotterTag.c_str());

  float minPtBin = m_ptThr - 1.;
  int   nPtBins  = (100. - m_ptThr) / 2;
  float maxPtBin = m_ptThr + nPtBins*2. - 1.;

  Float_t pi = TMath::Pi();

  // Declare histograms
  m_histoMap["hHLTvsEta"] = new TH1F(("hHLTvsEta"+plotterTag).c_str(), ("HLT occupancy vs #eta " + plotterTag +"; #eta; A.U.").c_str(), 12, ETA_BINS);
  m_histoMap["hGENvsEta"] = new TH1F(("hGENvsEta"+plotterTag).c_str(), ("GEN occupancy vs #eta " + plotterTag +"; #eta; A.U.").c_str(), 12, ETA_BINS);

  m_histoMap["hHLTvsPhi"] = new TH1F(("hHLTvsPhi"+plotterTag).c_str(), ("HLT occupancy vs #phi " + plotterTag +"; #phi; A.U.").c_str(), 20, -pi, pi);
  m_histoMap["hGENvsPhi"] = new TH1F(("hGENvsPhi"+plotterTag).c_str(), ("GEN occupancy vs #phi " + plotterTag +"; #phi; A.U.").c_str(), 20, -pi, pi);

  m_histoMap["hHLTvsPt"] = new TH1F(("hHLTvsPt"+plotterTag).c_str(), ("HLT occupancy vs p_{T} " + plotterTag +"; p_{T} GeV; A.U.").c_str(), nPtBins, minPtBin, maxPtBin);
  m_histoMap["hGENvsPt"] = new TH1F(("hGENvsPt"+plotterTag).c_str(), ("GEN occupancy vs p_{T} " + plotterTag +"; p_{T} GeV; A.U.").c_str(), nPtBins, minPtBin, maxPtBin);

  m_histoMap["hHLTvsPtThr"] = new TH1F(("hHLTvsPtThr"+plotterTag).c_str(), ("HLT occupancy vs p_{T} cut " + plotterTag +"; p_{T} cut GeV; A.U.").c_str(), nPtBins, minPtBin, maxPtBin );
  m_histoMap["hGENvsPtThr"] = new TH1F(("hGENvsPtThr"+plotterTag).c_str(), ("GEN occupancy vs p_{T} cut " + plotterTag +"; p_{T} cut GeV; A.U.").c_str(), nPtBins, minPtBin, maxPtBin );

  m_histoMap["hHLTvsInter"] = new TH1F(("hHLTvsInter"+plotterTag).c_str(), ("HLT occupancy vs PU" + plotterTag +"; PU; A.U.").c_str(), 35, 0., 70.);
  m_histoMap["hGENvsInter"] = new TH1F(("hGENvsInter"+plotterTag).c_str(), ("GEN occupancy vs PU" + plotterTag +"; PU; A.U.").c_str(), 35, 0., 70.);


}

void RatePlotter::fill(ciemat::Event *ev)
{

  const ciemat::GenParticle *bestGenPart = 0; // bestGenPart is highest pT status 1 mu in eta range

  std::vector<ciemat::GenParticle>::const_iterator genPartIt  = ev->genParticles.begin();
  std::vector<ciemat::GenParticle>::const_iterator genPartEnd = ev->genParticles.end();
      
  for (; genPartIt != genPartEnd; ++genPartIt) 
    {
      if (genPartIt->status == 1           &&
	  abs(genPartIt->pdgId)==13        &&
	  fabs(genPartIt->eta) <= m_maxEta &&
	  genPartIt->pt >= m_ptThr)
	{
	  if (bestGenPart)
	    {
	      if (bestGenPart->pt < genPartIt->pt)
		bestGenPart = &(*genPartIt);
		}
	  else
	    {
	      bestGenPart = &(*genPartIt);
	    }

	}
    }

  if (bestGenPart) 
    {
      m_histoMap["hGENvsEta"]->Fill(bestGenPart->eta);
      m_histoMap["hGENvsPhi"]->Fill(bestGenPart->phi);
      m_histoMap["hGENvsPt"]->Fill(bestGenPart->pt);

      for (int iBin = 0; iBin <= m_histoMap["hGENvsPtThr"]->GetNbinsX(); ++ iBin)
	{
	  float binPt = m_histoMap["hGENvsPtThr"]->GetBinCenter(iBin);	
	  if (bestGenPart->pt >= binPt) m_histoMap["hGENvsPtThr"]-> Fill(binPt);
	}

      Int_t nInt = ev.genInfos.size() > 0 : ev.genInfos.trueNumberOfInteractions : -1;
      m_histoMap["hGENvsInter"]->Fill(nInt);

    }
  

  std::vector<ciemat::HLTObject>::const_iterator hltObjIt  = ev->hlt.objects.begin();
  std::vector<ciemat::HLTObject>::const_iterator hltObjEnd = ev->hlt.objects.end();

  const HLTObject * bestHLT = 0; // bestHLT is highest pT HLT from m_filterTag HLT filter in eta range

  for (; hltObjIt != hltObjEnd; ++hltObjIt)
    {
      if (hltObjIt->filterTag.find(m_filterTag) != std::string::npos &&
	  hltObjIt->pt >= m_ptThr && fabs(hltObjIt->eta) <= m_maxEta )
	{
	  if (bestHLT)
	    {
	      if (bestHLT->pt < hltObjIt->pt)
		bestHLT = &(*hltObjIt);
	    }
	  else
	    {
	      bestHLT = &(*hltObjIt);
	    }
	}
    }
  
  if (bestHLT)
    {
      m_histoMap["hHLTvsEta"]->Fill(bestHLT->eta);
      m_histoMap["hHLTvsPhi"]->Fill(bestHLT->phi);
      m_histoMap["hHLTvsPt"]->Fill(bestHLT->pt);

      for (int iBin = 0; iBin <= m_histoMap["hHLTvsPtThr"]->GetNbinsX(); ++iBin)
	{
	  float binPt = m_histoMap["hHLTvsPtThr"]->GetBinCenter(iBin);
	  if (bestHLT->pt >= binPt) m_histoMap["hHLTvsPtThr"]-> Fill(binPt);
	}

      Int_t nInt = ev.genInfos.size() > 0 : ev.genInfos.trueNumberOfInteractions : -1;
      m_histoMap["hHLTvsInter"]->Fill(nInt);


    }
  
}

void RatePlotter::plotAndSave(float nEntries)
{

  std::map<std::string, TH1F*>::iterator histoMapIt  =  m_histoMap.begin();
  std::map<std::string, TH1F*>::iterator histoMapEnd =  m_histoMap.end();

  for(; histoMapIt!=histoMapEnd; ++histoMapIt)
    {
      histoMapIt->second->Sumw2();
      histoMapIt->second->Scale(1./nEntries);
    }

}


int main(int argc, char* argv[]){

  if (argc < 2) 
    {
      std::cout << "Usage : " << argv[0] << " PATH_TO_FILE \n";
      exit(100);
    }

  // Input root file
  TString fileName = argv[1]; // "/data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_v5_53X_HLT701_25PU_25ns_QCDMu3050.root";

  int tagBegin = fileName.Index("MuonHltTree")+12;
  int tagEnd   = fileName.Index(".root");
  TString fileTag  = fileName(tagBegin,std::max(0,tagEnd-tagBegin));

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("CMS Root Application", &argc, argv);

  setTDRStyle();

  std::vector<RatePlotter> plotters;

  plotters.push_back(RatePlotter(16,2.4,"hltL1sMu16"));
  plotters.push_back(RatePlotter(16,2.4,"hltL2fL1sMu16L1f0L2Filtered16Q"));
  plotters.push_back(RatePlotter(24,2.4,"hltL3fL1sMu16L1f0L2f16QL3Filtered24Q"));
  plotters.push_back(RatePlotter(40,2.4,"hltL3fL1sMu16L1f0L2f16QL3Filtered40Q"));
  plotters.push_back(RatePlotter(24,2.4,"hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15"));
	   
  plotters.push_back(RatePlotter(16,2.1,"hltL1sMu16"));
  plotters.push_back(RatePlotter(16,2.1,"hltL2fL1sMu16L1f0L2Filtered16Q"));
  plotters.push_back(RatePlotter(24,2.1,"hltL3fL1sMu16L1f0L2f16QL3Filtered24Q"));
  plotters.push_back(RatePlotter(40,2.1,"hltL3fL1sMu16L1f0L2f16QL3Filtered40Q"));
  plotters.push_back(RatePlotter(24,2.1,"hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15"));

  // Initialize pointers to summary and full event structure
 
  ciemat::Event* ev = new ciemat::Event();
  TTree* tree;
  TBranch* evBranch;

  // Open file, get tree, set branches

  std::cout << "[" << argv[0] << "] Processing file " << fileName.Data() << std::endl;

  TFile* inputFile = TFile::Open(fileName,"READONLY");
  tree = (TTree*)inputFile->Get("MUHLTTREE");
  if (!tree) inputFile->GetObject("MuonHltTree/MUHLTTREE",tree);

  evBranch = tree->GetBranch("event");
  evBranch->SetAddress(&ev);

  system("mkdir -p results");
  
  TString outputFileName("results/rateResults_"+fileTag+".root");  
  TFile* outputFile = TFile::Open(outputFileName,"RECREATE");  

  std::vector<RatePlotter>::iterator plotterIt  = plotters.begin();
  std::vector<RatePlotter>::iterator plotterEnd = plotters.end();
  
  for(; plotterIt!=plotterEnd; ++ plotterIt)
    {
      plotterIt->book(outputFile);
    }

  // Watch number of entries
  int nEntries = tree->GetEntriesFast();
  std::cout << "[" << argv[0] << "] Number of entries = " << nEntries << std::endl;

  for (Long64_t iEvent=0; iEvent<nEntries; ++iEvent) 
    {
      if (tree->LoadTree(iEvent)<0) break;

      evBranch->GetEntry(iEvent);

      plotterIt  = plotters.begin();
      plotterEnd = plotters.end();

      for(; plotterIt!=plotterEnd; ++ plotterIt)
      	{
      	  plotterIt->fill(ev);
      	}
      
    }

  plotterIt  = plotters.begin();
  plotterEnd = plotters.end();

  for(; plotterIt!=plotterEnd; ++ plotterIt)
    {
      plotterIt->plotAndSave(nEntries);
    }
  
  outputFile->Write();  

  // TCanvas * cGen = new TCanvas("GenVsEta","GenVsEta",500,500);

  // hGENvsEta->Scale(1./nEntries);
  // hGENvsEta->Draw();

  // cGen->SaveAs("genVsEta_"+fileTag+".C");

  // TCanvas * cL3 = new TCanvas("L3VsPt","L3VsPt",500,500);

  // hL3vsEta->Scale(1./nEntries);
  // hL3vsEta->Draw();

  // cL3->SaveAs("hltVsEta_"+fileTag+".C");

  if (!gROOT->IsBatch()) app->Run();

  return 0;
}
