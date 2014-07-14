#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TString.h"

class MyAnalysis {
public:
      MyAnalysis();
      virtual ~MyAnalysis();

      // Add new samples to the game
      void AddDataSample(const TString& id, const TString& file, double luminosity);
      void AddMCSignalSample(const TString& id, const TString& file, int maxevents, double xsection);
      void AddMCSample(const TString& id, const TString& file, int maxevents, double xsection);

      // Set tree to the relevant sample
      void SetTree(int i);

      // Read summary, get status
      ciemat::Summary* ReadSummary(int entry);

      // Read event, get status
      ciemat::Event* ReadEvent(int entry);

      // Get summary from the current tree
      const ciemat::Summary& GetSummary(){return *_summary;};
      // Get pointer to summary from the current tree
      ciemat::Summary* GetSummaryPointer(){return _summary;};

      // Get event from the current tree
      const ciemat::Event& GetEvent(){return *_event;};
      // Get pointer to event from the current tree
      ciemat::Event* GetEventPointer(){return _event;};

      // Is data?
      bool IsData(){return _currentIndex==_dataIndex;};
      // Get luminosity
      double GetLumi(){return _lumi;};
      // Get id that represents data
      int GetDataIndex() {return _dataIndex;};
      // Get id that represents MC signal
      int GetSignalIndex() {return _signalIndex;};
      // Get current tree
      TTree* GetTree() {return _tree;};
      // Get current file
      TFile* GetFile() {return _file;};
      // Get current index
      int GetIndex() {return _currentIndex;};

      // Get number of samples
      int GetNumberOfSamples() {return _sampleId.size();};
      // Get number of events to be read in sample
      int GetNumberOfEvents(int i) {return _sampleNevents[i];};
      // Get sampleId
      TString GetSampleId(int i) {return _sampleId[i];};
      // Get sampleFile
      TString GetSampleFile(int i) {return _sampleFile[i];};
      // Get cross section
      double GetSampleXsection(int i) {return _sampleXsection[i];};
      // Get sample weight
      double GetSampleWeight(int i) {return _sampleWeight[i];};

      // Book 1D Plot for data and all bkgd components
      void AddPlot1D(const TString& name, const TString& title, int nbins, double xmin, double xmax);
      // Fill histogram for the 1D Plot
      void FillPlot1D(const TString& name, int isample, double value, double weight=1.);
      // Draw 1D plot with data and all bckg components
      void DrawPlot1D(const TString& name, const TString& suffix="");
      void DrawPlot1D_MC(const TString& name, const TString& suffix="");
      void DrawPlot1DFitData(const TString& name, const TString& suffix="", double width=20., double xmin=0, double xmax=0);
      void DrawPlot1DFitDataDstar(const TString& name, const TString& suffix="", double width=20.);
      void FitDMass(const TString& name, const TString& suffix="", double width=20., double xmin=0, double xmax=0);
      void FitDMassGaussian(const TString& name, const TString& suffix="", double width=20., double xmin=0, double xmax=0);
      void FitDstarMass(const TString& name, const TString& suffix="", double width=20.);

      // Book 2D Plot for data and all bkgd components
      void AddPlot2D(const TString& name, const TString& title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax);
      // Fill histogram for the 2D Plot
      void FillPlot2D(const TString& name, int isample, double xvalue, double yvalue, double weight=1.);
      // Draw 2D plot with data and all bckg components
      void DrawPlot2D(const TString& name, const TString& suffix="");

private:
      double _lumi;
      std::vector<TString> _sampleFile;
      std::vector<TString> _sampleId;
      std::vector<int> _sampleNevents;
      std::vector<double> _sampleXsection;
      std::vector<double> _sampleWeight;

      int _dataIndex;
      int _signalIndex;

      int _currentIndex;
      TFile* _file;
      TTree* _tree;
      ciemat::Summary* _summary;
      ciemat::Event* _event;
      TBranch* _bSummary;
      TBranch* _bEvent;

      std::vector<TH1D*> hists_1D;
};
