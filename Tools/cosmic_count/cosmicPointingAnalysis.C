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

using namespace muon_hlt;

int main(int argc, char* argv[]){


  if (argc < 2) 
    {
      std::cout << "Usage : " << argv[0] << " PATH_TO_FILE \n";
      exit(100);
    }

  // Input root file
  TString fileName = argv[1];

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 

  // Initialize Root application
  TRint* app = new TRint("CMS Root Application", &argc, argv);

  //setTDRStyle();
  
  // Initialize pointers to summary and full event structure
 
  muon_hlt::Event* ev = new muon_hlt::Event();
  TTree* tree;
  TBranch* evBranch;

  // Open file, get tree, set branches

  std::cout << "[" << argv[0] << "] Processing file " << fileName.Data() << std::endl;

  TFile* inputFile = TFile::Open(fileName,"READONLY");
  tree = (TTree*)inputFile->Get("MUHLTTREE");
  if (!tree) inputFile->GetObject("MuonHltTree/MUHLTTREE",tree);

  evBranch = tree->GetBranch("event");
  evBranch->SetAddress(&ev);

  TFile* outputFile = TFile::Open("results.root","RECREATE");  

  TH1F* hNPixelHits     = new TH1F("hNPixelHits","Number of pixel hits in muon; # pixel hits; # muons per hit bin", 5,-0.5,4.5);
  TH1F* hNTrackerLayers = new TH1F("hNTrackerLayers","Number of tracker layers in muon; # tracker layers; # muons per layer bin", 11,-0.5,10.5);

  TH1F* hNMuonEta = new TH1F("hNMuonEta","Muon #eta; #eta; # muons per #eta bin", 100,-2.5,2.5);
  TH1F* hNMuonPhi = new TH1F("hNMuonPhi","Muon #phi; #phi; # muons per #phi bin", 100,-TMath::Pi(),TMath::Pi());
  TH1F* hNMuonPt  = new TH1F("hNMuonPt","Muon p_{T}; p_{T}; # muons per p_{T} bin", 100,0,700.);

  TH1F* hNMuonDxy = new TH1F("hNMuonDxy","Muon d_{xy}; #d_{xy}; # muons per #d_{xy} bin", 25,-50,50);
  TH1F* hNMuonDz  = new TH1F("hNMuonDz","Muon d_{z}; #d_{z}; # muons per #d_{z} bin", 25,-250,250);

  TH1F* hNMuonDxyCraft = new TH1F("hNMuonDxyCraft","Muon #d_{xy}; #d_{xy}; # muons per #d_{xy} bin", 25,-50,50);
  TH1F* hNMuonDzCraft  = new TH1F("hNMuonDzCraft","Muon #d_{z}; #d_{z}; # muons per #d_{z} bin", 25,-250,250);
  
  // Watch number of entries
  int nEntries = tree->GetEntriesFast();
  std::cout << "[" << argv[0] << "] Number of entries = " << nEntries << std::endl;

  int nFilteredEvents = 0;
  
  for (Long64_t iEvent=0; iEvent<nEntries; ++iEvent) 
    {
      if (tree->LoadTree(iEvent)<0) break;

      evBranch->GetEntry(iEvent);

      std::vector<muon_hlt::Muon>::const_iterator muonIt  = ev->muons.begin();
      std::vector<muon_hlt::Muon>::const_iterator muonEnd = ev->muons.end();

      for (; muonIt != muonEnd; ++muonIt)
	{
	  
	  muon_hlt::Muon const mu = (*muonIt);

	  if (!mu.isGlobal) continue;
	  
	  hNPixelHits->Fill(mu.nHitsPixel);
	  hNTrackerLayers->Fill(mu.nHitsTracker);

	  hNMuonEta->Fill(mu.eta);
	  hNMuonPhi->Fill(mu.phi);
	  hNMuonPt->Fill(mu.pt);

	  hNMuonDxy->Fill(mu.dxy);
	  hNMuonDz->Fill(mu.dz);
	  
	  if (mu.nHitsPixel   > 0 &&
	      mu.nHitsTracker > 7)
	    {	    
	      nFilteredEvents++;
	      hNMuonDxyCraft->Fill(mu.dxy);
	      hNMuonDzCraft->Fill(mu.dz);
	    }
	}
      
    }

  std::cout << "# of fileterd events in sample : " << nFilteredEvents << std::endl;

  outputFile->Write();
  
  if (!gROOT->IsBatch()) app->Run();

  return 0;
}
