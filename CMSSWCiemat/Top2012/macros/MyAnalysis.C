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
#include "TLatex.h"
#include "TColor.h"
#include "TMath.h"
#include "TF1.h"

#include "tdrstyle.C"

#include "MuTree.h"
#include "MyAnalysis.h"

MyAnalysis::MyAnalysis() {
      _file = 0;
      _tree = 0; 
      _currentIndex = -1; 
      _event = new ciemat::Event;
      _summary = new ciemat::Summary;
      _bSummary = 0;
      _bEvent = 0;
      _dataIndex = -1;
      _signalIndex = -1;
};
MyAnalysis::~MyAnalysis() { 
      delete _summary;
      delete _event;
};

void MyAnalysis::AddDataSample(const TString& id, const TString& file, double luminosity) {
  _lumi = luminosity;
  _sampleId.push_back(id);
  _sampleFile.push_back(file);

  TFile file_tmp(file,"READONLY");
  TTree* tree_tmp = 0;
  file_tmp.GetObject("MUTREE",tree_tmp);
  if (!tree_tmp) file_tmp.GetObject("MuTree/MUTREE",tree_tmp);
  if (!tree_tmp) printf("Error reading tree for file %s; crash expected!!\n", file.Data());
  _sampleNevents.push_back(tree_tmp->GetEntriesFast());

  _sampleXsection.push_back(_sampleNevents.back()/_lumi);
  _sampleWeight.push_back(1.);
  _dataIndex = _sampleId.size()-1;
}

void MyAnalysis::AddMCSample(const TString& id, const TString& file, int maxevents, double xsec) {
  if (_dataIndex<0) {
      printf(">>> Warning: you should call AddDataSample first, to define the luminosity!!!\n");
      printf(">>> SAMPLE NOT ADDED!\n");
      return;
  }
  _sampleId.push_back(id);
  _sampleFile.push_back(file);
  _sampleXsection.push_back(xsec);

  TFile file_tmp(file,"READONLY");
  TTree* tree_tmp = 0;
  file_tmp.GetObject("MUTREE",tree_tmp);
  if (!tree_tmp) file_tmp.GetObject("MuTree/MUTREE",tree_tmp);
  if (!tree_tmp) printf("Error reading tree for file %s; crash expected!!\n", file.Data());
  int maxeventsInTree = tree_tmp->GetEntriesFast();
  if (maxevents<0 || maxevents>maxeventsInTree) {
      _sampleNevents.push_back(maxeventsInTree);
  } else {
      _sampleNevents.push_back(maxevents);
  }

  _sampleWeight.push_back(_lumi*_sampleXsection.back()/_sampleNevents.back());
}

void MyAnalysis::AddMCSignalSample(const TString& id, const TString& file, int maxevents, double xsec) {
  if (_dataIndex<0) {
      printf(">>> Warning: you should call AddDataSample first, to define the luminosity!!!\n");
      printf(">>> SAMPLE NOT ADDED!\n");
      return;
  }
  AddMCSample(id, file, maxevents, xsec);
  _signalIndex = _sampleId.size()-1;
}

void MyAnalysis::AddPlot1D(const TString& name, const TString& title, int nbins, double xmin, double xmax) {
      for (unsigned int i=0; i<_sampleId.size(); ++i) {
            bool existing = false;
            for (unsigned int j=0; j<hists_1D.size(); ++j) {
                  TString thisname = hists_1D[j]->GetName();
                  if (thisname == _sampleId[i]+"_"+name) {
                        existing = true; 
                        break;
                  }
            }
            if (existing) continue;

            hists_1D.push_back(new TH1D(_sampleId[i]+"_"+name, title, nbins, xmin, xmax));
            hists_1D[hists_1D.size()-1]->Sumw2();
      }
}
  
void MyAnalysis::FillPlot1D(const TString& name, int isample, double value, double weight) {
      for (unsigned int j=0; j<hists_1D.size(); ++j) {
            if (hists_1D[j]->GetName()==_sampleId[isample]+"_"+name) {
                  hists_1D[j]->Fill(value,_sampleWeight[isample]*weight);
                  return;
            }
      }
}
  
void MyAnalysis::DrawPlot1D(const TString& name, const TString& suffix) {
      setTDRStyle();

      //gROOT->SetStyle("Pub"); 
      gROOT->SetStyle("Plain"); 

      //gStyle->SetPadGridX(true); 
      //gStyle->SetPadGridY(true);
      gStyle->SetOptStat(0);

      TLegend* leg = new TLegend(0.52,0.6,0.85,0.85);
      //TLegend* leg = new TLegend(0.80,0.75,1.00,1.00); // to see the whole histogram
      leg->SetFillColor(0);

      THStack* hMCStack = new THStack(name,name+" histograms");
      TH1D* hData;

      // Get data; stack all MC except signal MC
      unsigned int nhists = hists_1D.size();
      //int colors[10] = {46, 47, 48, 49, 50, 51, 52, 53, 54, 55};
      //int colors[10] = {6, 7, 8, 9, 11, 44, 46, 2, 3, 4}; 
      //int colors[10] = {TColor::GetColor("#ff55ff"), TColor::GetColor("#00ffff"), TColor::GetColor("#ff5500"), 9, 11, 44, 46, 2, 3, 4}; 
      int colors[10] = {TColor::GetColor("#00ddff"), TColor::GetColor("#0055ff"), kMagenta, kRed, 11, 44, 46, 2, 3, 4}; 
      int mcindex = -1;
      for (unsigned int j=0; j<nhists; ++j) {
            TString histname = hists_1D[j]->GetName();
            TString suffix = "_" + name;
            if (!histname.EndsWith(suffix)) continue;

            for (unsigned int i=0; i<_sampleFile.size(); ++i) {
                  TString prefix = _sampleId[i] + "_";
                  if (!histname.BeginsWith(prefix)) continue;
                  if (histname==_sampleId[_dataIndex]+"_"+name) {
                        hData = hists_1D[j];
                        hData->SetMarkerStyle(20);
                        hData->SetMarkerSize(1.0);
                        leg->AddEntry(hData,_sampleId[i].Data(),"P");
                        break;
                  } else {
                        mcindex++;
                        int color = colors[mcindex%10];
                        hists_1D[j]->SetLineWidth(3);
                        hists_1D[j]->SetFillColor(color);
                        hists_1D[j]->SetLineColor(TColor::GetColorDark(color));
                        leg->AddEntry(hists_1D[j],_sampleId[i].Data(),"F");
                        // Do not add the signal component to the stack yet
                        if (histname!=_sampleId[_signalIndex]+"_"+name) {
                              hMCStack->Add(hists_1D[j]);
                        }
                        break;
                  }
            }
      }
      
      // Add the signal component to the stack now
      for (unsigned int j=0; j<nhists; ++j) {
            if (hists_1D[j]->GetName()==_sampleId[_signalIndex]+"_"+name) {
                  hMCStack->Add(hists_1D[j]);
                  break;
            }
      }
      
      TString c1name = "c1_" + name;
      TCanvas* c1 = new TCanvas(c1name.Data(),c1name.Data(),10,10,600,600);
      hData->SetXTitle(hData->GetTitle());
      hData->SetTitle("");
      hData->SetTitleOffset(1.2);
      if (hData->GetMinimum()>0.) hData->SetMinimum(0.);
      hData->Draw("e");
      //TLatex* preliminary = new TLatex(0.56,0.92,"CMS preliminary");
      TLatex* preliminary = new TLatex(0.05,0.92,"CMS preliminary");
      preliminary->SetNDC();
      preliminary->SetTextFont(42);
      preliminary->Draw();
      hMCStack->Draw("samehist");
      hData->Draw("esame");
      leg->Draw();

      gPad->SetTicks(1,1);
      gPad->RedrawAxis();

      if (suffix!="") {
            c1->SaveAs(name+"_"+suffix+".root");
            c1->SaveAs(name+"_"+suffix+".jpg");
            c1->SaveAs(name+"_"+suffix+".pdf");
      } else {
            c1->SaveAs(name+".root");
            c1->SaveAs(name+".jpg");
            c1->SaveAs(name+".pdf");
      }
}
  
void MyAnalysis::SetTree(int i) {
      if (_file) {
            _file->Close();
            _file = NULL;
      }

      printf("Processing sample '%s'...\n", _sampleFile[i].Data());
      _currentIndex = i;
      _file = new TFile(_sampleFile[i],"READONLY");
      _tree = (TTree*)_file->Get("MUTREE");
      _file->GetObject("MUTREE",_tree);
      if (!_tree) _file->GetObject("MuTree/MUTREE",_tree);

      _bSummary = _tree->GetBranch("summary");
      _bEvent = _tree->GetBranch("event");
      _bSummary->SetAddress(&_summary);
      _bEvent->SetAddress(&_event);

      int nentriesInTree = _tree->GetEntriesFast();
      printf("\tReading %d entries from a total of %d\n", _sampleNevents[i], nentriesInTree);

}

ciemat::Summary* MyAnalysis::ReadSummary(int iEvent) {
      if (_tree->LoadTree(iEvent)<0) return 0;
      _bSummary->GetEntry(iEvent);
      return _summary;
}

ciemat::Event* MyAnalysis::ReadEvent(int iEvent) {
      if (_tree->LoadTree(iEvent)<0) return 0;
      _bEvent->GetEntry(iEvent);
      return _event;
}
