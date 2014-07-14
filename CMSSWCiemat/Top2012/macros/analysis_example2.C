#include "TROOT.h"
#include "TSystem.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TTree.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom2.h"

#include "MuTree.h"
#include "MyAnalysis.h"

using namespace ciemat;

class TopAnalysis : public MyAnalysis {
public:
      TopAnalysis(){};
      virtual ~TopAnalysis(){};

      // Preselection
      bool Preselect();

      // Selection
      bool Select();
private:
};

int main(int argc, char** argv){

  // Directory where files are sitting
  TString dir = "./";

  // Initialize analysis structure
  TopAnalysis ana;

  // Add data sample. Input is: (titleId, file, luminosity)
  // One should add this sample first, to define the luminosity for normalizations
  double lumi_invpb = 20000.;
  TString chlumi_invfb = TString::Itoa(int(lumi_invpb/1000+0.5),10);
  ana.AddDataSample("Data, L=" + chlumi_invfb + " fb^{-1}", dir+"Data.root", lumi_invpb);

  // Set maximum number of events in MC
  int maxevents = -1; // use all MC events, final plots
  //int maxevents = 3000000; // for fast tests but reasonabe MC statistics
  //int maxevents = 500000; // for very fast tests with MC
  //int maxevents = 50000; // for very fast tests with MC
 
  // Add signal sample. Input is: (titleId, file, maxevents, xsection in pb)
  double xsec_signal = 286 / lumi_invpb * 0.25*0.75; // not real, just for test
  ana.AddMCSignalSample("Signal", dir+"Signal.root",maxevents, xsec_signal);

  // Add MC samples. Input is: (titleId, file, maxevents, xsection in pb)
  double xsec_bkgd = 286 / lumi_invpb * 0.25*0.25; // not real, just for test
  ana.AddMCSample("Background", dir+"Bkgd.root",maxevents, xsec_bkgd);

  // Initialize 1D histograms
  ana.AddPlot1D("hmet", "Missing E_{T} [GeV]", 50, 0., 200.);

  // Loop on samples
  int nsamples = ana.GetNumberOfSamples();
  for (unsigned int iSample=0; iSample<nsamples; ++iSample) {

      // Set tree
      ana.SetTree(iSample);

      // Loop on events for the current sample
      int nevents = ana.GetNumberOfEvents(iSample);
      for (int iEvent=0; iEvent<nevents; iEvent++) {

            if (iEvent%1000000==0) printf("... event index %d\n", iEvent);

            // Read summary header 
            Summary* psummary = ana.ReadSummary(iEvent);
            if (!psummary) break;
            const Summary& summary = *psummary; // if you prefer not to use the pointer

            // Preselect (summary branch must be read before)
            if (!ana.Preselect()) continue;

            // Now get the full event, once we know that summary conditions are satisfied
            Event* pev = ana.ReadEvent(iEvent);
            if (!pev) break;
            const Event& ev = *pev; // if you prefer not to use the pointer

            // Select (event branch must be read before)
            if (!ana.Select()) continue;

            // Fill histogram
            ana.FillPlot1D("hmet", iSample, ev.met.met);

      }

  }

  // To see things interactively (commnent otherwise) 
  TRint* app = new TRint("Top Analysis", &argc, argv); 

  // Draw histograms
  ana.DrawPlot1D("hmet");

  // To see things interactively (commnent otherwise) 
  if (!gROOT->IsBatch()) app->Run();

  return 0;
}

bool TopAnalysis::Preselect() {
  // Cut on summary information to save processing time
  const Summary& summary = GetSummary();
  if (summary.nMuons+summary.nElectrons<2) return false; // look for events with >= 2 leptons
  if (summary.nJetsPt30<2) return false; // look for events with at least 2 jets with pt>30 GeV
  if (summary.bestCSV<0.2) return false; // at least 1 jet with minimal b-tagging

  return true;    
}

bool TopAnalysis::Select() {
  // Select if there are at least two leptons wit some phase space cuts
  const double ptmuCut = 25.;
  const double etamuCut = 25.;
  const double ptelCut = 2.1;
  const double etaelCut = 2.5;

  const Event& ev = GetEvent();

  int nmusel = 0;
  unsigned int nmuons = ev.muons.size();
  for (unsigned int im=0; im<nmuons; ++im) {
      const Muon& mu = ev.muons[im];
      if (mu.pt>ptmuCut && fabs(mu.eta)<etamuCut) nmusel++;
  }

  int nelsel = 0;
  unsigned int nelectrons = ev.electrons.size();
  for (unsigned int ie=0; ie<nelectrons; ++ie) {
      const Electron& el = ev.electrons[ie];
      if (el.pt>ptelCut && fabs(el.eta)<etaelCut) nelsel++;
  }

  if (nmusel+nelsel>=2) return true; else return false;
}
