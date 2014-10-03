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

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <vector>
#include <map>

using namespace ciemat;

// BASE cuts and definitions
// to configute the analysis
const float ETA_BINS[17] = {-5., -3., -2.4, -2.1, -1.6, -1.2, -0.8, -0.4, 0., 0.4, 0.8, 1.2, 1.6, 2.1, 2.4, 3., 5.};

class RatePlotter 
{

public:

  RatePlotter(float ptThr, float maxEta, std::string filterTag, bool controlPlots);
  ~RatePlotter();

  void book(TFile* fileOut);
  void fill(ciemat::Event *ev);

  void plotAndSave(float scaleFactor);

private:

  float m_ptThr;
  float m_maxEta;
  std::string m_filterTag;
  bool m_controlPlots;

  std::map<std::string, TH1F*> m_histoMap;

};

RatePlotter::RatePlotter(float ptThr, float maxEta, std::string filterTag, bool controlPlots) :
  m_ptThr(ptThr), m_maxEta(maxEta), m_filterTag(filterTag), m_controlPlots(controlPlots)
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
  m_histoMap["hHLTvsEta"] = new TH1F(("hHLTvsEta"+plotterTag).c_str(), ("HLT occupancy vs #eta " + plotterTag +";HLT mu #eta;A.U.").c_str(), 16, ETA_BINS);
  m_histoMap["hGENvsEta"] = new TH1F(("hGENvsEta"+plotterTag).c_str(), ("GEN occupancy vs #eta " + plotterTag +";GEN mu #eta; A.U.").c_str(), 16, ETA_BINS);

  m_histoMap["hHLTvsPhi"] = new TH1F(("hHLTvsPhi"+plotterTag).c_str(), ("HLT occupancy vs #phi " + plotterTag +";HLT mu #phi [rad];A.U.").c_str(), 10, -pi, pi);
  m_histoMap["hGENvsPhi"] = new TH1F(("hGENvsPhi"+plotterTag).c_str(), ("GEN occupancy vs #phi " + plotterTag +";GEN mu #phi [rad]; A.U.").c_str(), 10, -pi, pi);

  m_histoMap["hHLTvsPt"] = new TH1F(("hHLTvsPt"+plotterTag).c_str(), ("HLT occupancy vs p_{T} " + plotterTag +";HLT mu p_{T} [GeV];A.U.").c_str(), nPtBins, minPtBin, maxPtBin);
  m_histoMap["hGENvsPt"] = new TH1F(("hGENvsPt"+plotterTag).c_str(), ("GEN occupancy vs p_{T} " + plotterTag +";GEN mu p_{T} [GeV]; A.U.").c_str(), nPtBins, minPtBin, maxPtBin);

  m_histoMap["hHLTvsPtThr"] = new TH1F(("hHLTvsPtThr"+plotterTag).c_str(), ("HLT occupancy vs p_{T} cut " + plotterTag +";HLT mu p_{T} cut [GeV]; A.U.").c_str(), nPtBins, minPtBin, maxPtBin );
  m_histoMap["hGENvsPtThr"] = new TH1F(("hGENvsPtThr"+plotterTag).c_str(), ("GEN occupancy vs p_{T} cut " + plotterTag +";GEN mu p_{T} cut [GeV]; A.U.").c_str(), nPtBins, minPtBin, maxPtBin );

  m_histoMap["hHLTvsInter"] = new TH1F(("hHLTvsInter"+plotterTag).c_str(), ("HLT occupancy vs PU" + plotterTag +"; N GEN vtx; A.U.").c_str(), 35, 0., 70.);
  m_histoMap["hGENvsInter"] = new TH1F(("hGENvsInter"+plotterTag).c_str(), ("GEN occupancy vs PU" + plotterTag +"; N GEN vtx; A.U.").c_str(), 35, 0., 70.);

  if (m_controlPlots)
    {
      m_histoMap["hGENMultip"] = new TH1F(("hGENMultip"+plotterTag).c_str(), ("GEN multiplicity of mu+ + mu- " + plotterTag +";# muons;A.U.").c_str(), 30, -0.5, 29.5);
      m_histoMap["hGENMultip1"] = new TH1F(("hGENMultip1"+plotterTag).c_str(), ("GEN multiplicity of mu+ + mu- of status 1" + plotterTag +";# muons;A.U.").c_str(), 30, -0.5, 29.5);

      m_histoMap["hGENvsEtaLowPt"] = new TH1F(("hGENvsEtaLowPt"+plotterTag).c_str(), ("GEN occupancy vs #eta p_{T} < 5 GeV " + plotterTag +";GEN mu #eta;A.U.").c_str(), 16, ETA_BINS);
      m_histoMap["hGENvsPhiLowPt"] = new TH1F(("hGENvsPhiLowPt"+plotterTag).c_str(), ("GEN occupancy vs #phi p_{T} < 5 GeV " + plotterTag +";GEN mu #phi [rad]; A.U.").c_str(), 10, -pi, pi);

      m_histoMap["hGENvsEtaAll"] = new TH1F(("hGENvsEtaAll"+plotterTag).c_str(), ("GEN occupancy vs #eta " + plotterTag +" all muons;GEN mu #eta; A.U.").c_str(), 16, ETA_BINS);
      m_histoMap["hGENvsPhiAll"] = new TH1F(("hGENvsPhiAll"+plotterTag).c_str(), ("GEN occupancy vs #phi " + plotterTag +" all muons;GEN mu #phi [rad]; A.U.").c_str(), 10, -pi, pi);
      m_histoMap["hGENvsPtAll"] = new TH1F(("hGENvsPtAll"+plotterTag).c_str(), ("GEN occupancy vs p_{T} " + plotterTag +" all muons;GEN mu p_{T} [GeV]; A.U.").c_str(), nPtBins, minPtBin, maxPtBin);
    }

}

void RatePlotter::fill(ciemat::Event *ev)
{

  const ciemat::GenParticle *bestGenPart = 0; // bestGenPart is highest pT status 1 mu in eta range

  std::vector<ciemat::GenParticle>::const_iterator genPartIt  = ev->genParticles.begin();
  std::vector<ciemat::GenParticle>::const_iterator genPartEnd = ev->genParticles.end();
      
  // std::cout << "Here is one event:" << std::endl;

  int nMu  = 0;
  int nMu1 = 0;

  for (; genPartIt != genPartEnd; ++genPartIt) 
    {

      if (m_controlPlots)
	{

	  if (abs(genPartIt->pdgId)==13)
	    {
	      nMu++;
	      if (genPartIt->status==1)
		{
		  nMu1++;
		  m_histoMap["hGENvsEtaAll"]->Fill(genPartIt->eta);
		  m_histoMap["hGENvsPhiAll"]->Fill(genPartIt->phi);
		  m_histoMap["hGENvsPtAll"]->Fill(genPartIt->pt);
		  if (genPartIt->pt < 5.)
		    {
		      m_histoMap["hGENvsEtaLowPt"]->Fill(genPartIt->eta);
		      m_histoMap["hGENvsPhiLowPt"]->Fill(genPartIt->phi);
		    }
		}
	    }
	}

      // if (abs(genPartIt->pdgId)==13) 
      // 	{
      // 	  std::cout << genPartIt->pdgId  << " "
      // 		    << genPartIt->status << " "
      // 		    << genPartIt->pt     << " "
      // 		    << genPartIt->eta    << std::endl;
      // 	}

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

  if (m_controlPlots)
    {
      m_histoMap["hGENMultip"]->Fill(nMu);
      m_histoMap["hGENMultip1"]->Fill(nMu1++);
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

      Int_t nInt = ev->genInfos.size() > 0 ? ev->genInfos.at(0).actualNumberOfInteractions : -1;
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

      Int_t nInt = ev->genInfos.size() > 0 ? ev->genInfos.at(0).actualNumberOfInteractions : -1;
      m_histoMap["hHLTvsInter"]->Fill(nInt);


    }
  
}

void RatePlotter::plotAndSave(float scaleFactor)
{

  std::map<std::string, TH1F*>::iterator histoMapIt  =  m_histoMap.begin();
  std::map<std::string, TH1F*>::iterator histoMapEnd =  m_histoMap.end();

  for(; histoMapIt!=histoMapEnd; ++histoMapIt)
    {
      histoMapIt->second->Sumw2();

      if (scaleFactor < 0. )
	{
	  // CB hack 4 simone if SF < 0
	  // normalise by Integral o plot
	  float area = histoMapIt->second->Integral();
	  scaleFactor = area>1 ? 1./area : 1;
	}
      histoMapIt->second->Scale(scaleFactor);
    }

}


int main(int argc, char* argv[]){

  float xSec      = 1;
  float filterEff = 1;

  if (argc < 2) 
    {
      std::cout << "Usage : " << argv[0] << " PATH_TO_FILE <X_SEC> <FILTER_EFF>\n";
      exit(100);
    }

  if (argc >= 3) 
    {
      xSec = std::atof(argv[2]);
      if (fabs(xSec) < 1E-20)
	{
	  std::cout << "atof(<X_SEC>) : " << xSec << " is too small\n";
	  exit(100);
	}
    }

  if (argc >= 4) 
    {
      filterEff = std::atof(argv[3]);
      if (fabs(filterEff) < 1E-20)
	{
	  std::cout << "atof(<FILTER_EFF>) : " << filterEff << " is too small\n";
	  exit(100);
	}
    }
  
  // Input root file
  TString fileName = argv[1];

  int tagBegin = fileName.Index("MuonHltTree")+12;
  int tagEnd   = fileName.Index(".root");
  TString fileTag  = fileName(tagBegin,std::max(0,tagEnd-tagBegin));

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("CMS Root Application", &argc, argv);

  setTDRStyle();
  
  std::vector<RatePlotter> plotters;

  plotters.push_back(RatePlotter(0.,5.,"genReferencePt0",true));
  plotters.push_back(RatePlotter(1.,3.,"genReferencePt11",false));

  plotters.push_back(RatePlotter(16,2.4,"hltL1sMu16",false));
  plotters.push_back(RatePlotter(16,2.4,"hltL2fL1sMu16L1f0L2Filtered16Q",false));
  plotters.push_back(RatePlotter(24,2.4,"hltL3fL1sMu16L1f0L2f16QL3Filtered24Q",false));
  plotters.push_back(RatePlotter(40,2.4,"hltL3fL1sMu16L1f0L2f16QL3Filtered40Q",false));
  plotters.push_back(RatePlotter(24,2.4,"hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15",false));
	   
  plotters.push_back(RatePlotter(16,2.1,"hltL1sMu16",false));
  plotters.push_back(RatePlotter(16,2.1,"hltL2fL1sMu16L1f0L2Filtered16Q",false));
  plotters.push_back(RatePlotter(24,2.1,"hltL3fL1sMu16L1f0L2f16QL3Filtered24Q",false));
  plotters.push_back(RatePlotter(40,2.1,"hltL3fL1sMu16L1f0L2f16QL3Filtered40Q",false));
  plotters.push_back(RatePlotter(24,2.1,"hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15",false));

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

  float scaleFactor = 1./nEntries * xSec * filterEff;

  plotterIt  = plotters.begin();
  plotterEnd = plotters.end();
  
  for(; plotterIt!=plotterEnd; ++ plotterIt)
    {
      plotterIt->plotAndSave(scaleFactor);
    }
  
  outputFile->Write();

  if (!gROOT->IsBatch()) app->Run();

  return 0;
}
