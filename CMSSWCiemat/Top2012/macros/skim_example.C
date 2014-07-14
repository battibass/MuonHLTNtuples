#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "MuTree.h"

using namespace ciemat;

int main(int argc, char** argv){

  // Input root file name
  TString chinput = "Signal.root";

  // Output root file name
  TString choutput = "Signal_skimmed.root";
  
  // Initialize pointers to summary and full event structure
  Summary summary;
  Summary* pointerToSummary = &summary;
  Event ev;
  Event* pointerToEvent = &ev;
  TTree* tree;
  TBranch* bSummary;
  TBranch* bEvent;

  // Open file, get tree, set branches
  printf("Processing sample '%s' ...\n", chinput.Data());
  TFile* pinput_sample = TFile::Open(chinput,"READONLY");
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

  // Output file and new tree to save selected events 
  TFile* output_file = new TFile(choutput,"RECREATE");
  TTree* output_tree = tree->CloneTree(0);

  for (Long64_t iEvent=0; iEvent<nentries; iEvent++) {
      if (tree->LoadTree(iEvent)<0) break;

      // First, access the summary
      bSummary->GetEntry(iEvent);

      // Cut on summary information to save processing time if event gets rejected
      if (summary.nElectrons==0) continue; // look for events with only 1 electron
      if (summary.maxMT<30.) continue; // look for events with maxMT>30
      if (summary.nJets<2) continue; // look for events with at least 2 jets

      // Now get the full event (mandatory to get all branches to write the full event)
      bEvent->GetEntry(iEvent);

      // Get the most significant b jet (additional cut just as an example)
      unsigned int njets = ev.jets.size();
      double btag = 0.;
      for (unsigned int ij=0; ij<njets; ++ij) {
            Jet jet = ev.jets[ij];
            if (jet.btagCSV>btag) {
                  btag = jet.btagCSV;
            }
      }

      // Select events with a jet with CSV > 0.5
      if (btag<0.5) continue;
      
      // Add event to new tree
      output_file->cd();
      output_tree->Fill();
  }

  // Write selected events and histograms into the output file
  output_file = output_tree->GetCurrentFile();
  output_file->Write();

  // final printout
  int selectedEvents = output_tree->GetEntriesFast();
  printf("Selected events: %d, original: %d, selected fraction: %.2f %%\n", selectedEvents, nentries, selectedEvents*100./nentries);

  // Close output file
  output_file->Close();

  return 0;
}
