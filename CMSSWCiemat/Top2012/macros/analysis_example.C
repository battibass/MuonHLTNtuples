#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TTree.h"
#include "TBranch.h"

#include "MuTree.h"

using namespace ciemat;

int main(int argc, char** argv){

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kFALSE); 

  // Initialize Root application
  TRint* app = new TRint("CMS Root Application", &argc, argv);

  // Canvas
  TCanvas* c1 = new TCanvas("c1","Top2012 analysis");
  gStyle->SetOptStat(1111111);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);

  // Input root file
  TString sample = "Data.root";
  
  // Declare histograms
  TH1D* hbtag = new TH1D("hbtag", "CSV discriminant, most significant jet", 25, 0., 1.);
  hbtag->Sumw2();
  
  // Initialize pointers to summary and full event structure
  Summary summary;
  Summary* pointerToSummary = &summary;
  Event ev;
  Event* pointerToEvent = &ev;
  TTree* tree;
  TBranch* bSummary;
  TBranch* bEvent;

  // Open file, get tree, set branches
  printf("Processing sample '%s' ...\n", sample.Data());
  TFile* pinput_sample = TFile::Open(sample,"READONLY");
  TFile& input_sample = *pinput_sample;
  tree = (TTree*)input_sample.Get("MUTREE");
  if (!tree) input_sample.GetObject("MuTree/MUTREE",tree);
  bSummary = tree->GetBranch("summary");
  bEvent = tree->GetBranch("event");
  bSummary->SetAddress(&pointerToSummary);
  bEvent->SetAddress(&pointerToEvent);

  // Watch number of entries
  int nentries = tree->GetEntriesFast();
  printf("Number of entries = %d\n", nentries);

  // Typical way to normalize the MC to data
  double lumi = 20000.; // 1/pb
  double xsection = 200.; // pb
  double skimFilterEfficiency = 1.; // different from 1 if MC was already skimmed
  double weight = lumi*xsection*skimFilterEfficiency / nentries;

  for (Long64_t iEvent=0; iEvent<nentries; iEvent++) {
      if (tree->LoadTree(iEvent)<0) break;

      // First, access the summary
      bSummary->GetEntry(iEvent);

      // Cut on summary information to save processing time
      if (summary.nMuons==0) continue; // look for events with at least 1 muon
      if (summary.maxMT<30.) continue; // look for events with maxMT>30
      if (summary.nJets<2) continue; // look for events with at least 2 jets

      // Now get the full event, once we know that summary conditions are satisfied
      bEvent->GetEntry(iEvent);

      // Get the most significant b jet
      unsigned int njets = ev.jets.size();
      double btag = 0.;
      for (unsigned int ij=0; ij<njets; ++ij) {
            Jet jet = ev.jets[ij];
            if (jet.btagCSV>btag) {
                  btag = jet.btagCSV;
            }
      }
      
      // Fill histogram
      hbtag->Fill(btag, weight);
  }

  hbtag->Draw("e");

  c1->SaveAs("btag.jpg");

  if (!gROOT->IsBatch()) app->Run();

  return 0;
}
