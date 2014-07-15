/////////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/EventSetup.h"

#include "MuonHLTNtuples/Tools/src/MuTree.h"
#include "TTree.h"

#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include <iostream>
#include <algorithm>

class MuonHltTreeProducer : public edm::EDAnalyzer {
public:
  MuonHltTreeProducer (const edm::ParameterSet &);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
private:

  edm::InputTag trigResultsTag_;
  edm::InputTag trigSummaryTag_;

  edm::InputTag muonTag_;
  edm::InputTag primaryVertexTag_;
  edm::InputTag beamSpotTag_;

  edm::InputTag genTag_;
  edm::InputTag pileUpInfoTag_;

  ciemat::Event event_;
  std::map<std::string,TTree*> htree_;

};

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"


void fillMuons( const edm::Handle<reco::MuonCollection> & muonCollection,
		const edm::Handle<std::vector<reco::Vertex> > & vertexCollection,
		const edm::Handle<reco::BeamSpot> & beamspot,
		std::vector<ciemat::Muon> & muons);

using namespace trigger;

MuonHltTreeProducer::MuonHltTreeProducer( const edm::ParameterSet & cfg ) :
  // Input collections
  trigResultsTag_(cfg.getUntrackedParameter<edm::InputTag>("TrigResultsTag", edm::InputTag("TriggerResults::HLT"))),
  trigSummaryTag_(cfg.getUntrackedParameter<edm::InputTag>("TrigSummaryTag", edm::InputTag("hltTriggerSummaryAOD::HLT"))),

  muonTag_(cfg.getUntrackedParameter<edm::InputTag>("MuonTag", edm::InputTag("muons"))),
  primaryVertexTag_(cfg.getUntrackedParameter<edm::InputTag>("PrimaryVertexTag", edm::InputTag("offlinePrimaryVertices"))),
  beamSpotTag_(cfg.getUntrackedParameter<edm::InputTag>("BeamSpotTag", edm::InputTag("offlineBeamSpot"))),

  genTag_(cfg.getUntrackedParameter<edm::InputTag>("GenTag", edm::InputTag("prunedGenParticles"))),
  pileUpInfoTag_(cfg.getUntrackedParameter<edm::InputTag>("PileUpInfoTag", edm::InputTag("pileupInfo")))
{

}

void MuonHltTreeProducer::beginJob() {

  edm::Service<TFileService> fs;
  htree_["hMuHltTree"] = fs->make<TTree>("MUHLTTREE","HLT Muon Tree");

  int splitBranches = 2;
  htree_["hMuHltTree"]->Branch("event",&event_,64000,splitBranches);

}


void MuonHltTreeProducer::beginRun(const edm::Run & run, const edm::EventSetup & config )  {
  
}


void MuonHltTreeProducer::endJob() {

}

void MuonHltTreeProducer::analyze (const edm::Event & ev, const edm::EventSetup &)
{

  /// clearing branch variables
  event_.hlt.triggers.clear();
  event_.hlt.objects.clear();

  event_.genParticles.clear();
  event_.genInfos.clear();
  event_.muons.clear();


  // Run, luminosity block, event
  event_.runNumber = ev.id().run();
  event_.luminosityBlockNumber = ev.id().luminosityBlock();
  event_.eventNumber = ev.id().event();


  // Fill GenInfo contents, but only in MC
  ciemat::GenInfo genInfo;
  if (!ev.isRealData()) {
    // Pileup information for MC
    genInfo.trueNumberOfInteractions = -1.;
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    if (pileUpInfoTag_.label() != "none") {

      if (ev.getByLabel(pileUpInfoTag_, PupInfo)) {

	std::vector<PileupSummaryInfo>::const_iterator PVI;
	for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      
	  int BX = PVI->getBunchCrossing();
      
	  if(BX == 0) { 
	    genInfo.trueNumberOfInteractions = PVI->getTrueNumInteractions();
	    continue;
	  }
	}
      } 
      else {
	edm::LogError("") << "[MuonHltTreeProducer]: Pile-Up Info collection does not exist !!!";
      }
    }
    
  }


  // Trigger
  edm::Handle<edm::TriggerResults> triggerResults;
  edm::Handle<TriggerEvent> triggerEvent;
  if (trigResultsTag_.label() != "none" &&
      trigSummaryTag_.label() != "none" ) {
      
    if (ev.getByLabel(trigResultsTag_, triggerResults) &&
	ev.getByLabel(trigSummaryTag_, triggerEvent)) {

      const edm::TriggerNames & triggerNames = ev.triggerNames(*triggerResults);
    
      for (unsigned int itrig=0; itrig<triggerNames.size(); ++itrig) {
	if (triggerResults->accept(itrig)) {
	  std::string thisPathName = triggerNames.triggerName(itrig);
	  event_.hlt.triggers.push_back(thisPathName);
	  
	  // std::cout << "Trigger Results accepts : " 
	  // 	  << thisPathName << std::endl;
	  
	}
      }
      
      const size_type nFilters(triggerEvent->sizeFilters());
      for (size_type iFilter=0; iFilter!=nFilters; ++iFilter) {
	
	std::string filterTag = triggerEvent->filterTag(iFilter).encode();
	
	// std::cout << "TriggerEvent filter tag : "
	// 		<< filterTag << std::endl;
	
	Keys objectKeys = triggerEvent->filterKeys(iFilter);
	const TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
	
	for (size_type iKey=0; iKey<objectKeys.size(); ++iKey) {
	  
	  size_type objKey = objectKeys.at(iKey);
	  const TriggerObject& triggerObj(triggerObjects[objKey]);
	  
	  ciemat::HLTObject hltObj;
	  
	  float trigObjPt = triggerObj.pt();
	  float trigObjEta = triggerObj.eta();
	  float trigObjPhi = triggerObj.phi();
	  
	  // std::cout << "\t pT : "  << trigObjPt
	  // 	  << "   eta : " << trigObjEta
	  // 	  << "   phi : " << trigObjPhi
	  // 	  <<std::endl;
	  
	  hltObj.filterTag = filterTag;
	  hltObj.pt = trigObjPt;
	  hltObj.eta = trigObjEta;
	  hltObj.phi = trigObjPhi;
	  
	  event_.hlt.objects.push_back(hltObj);
	  
	}       
      }
    }
    else {
      edm::LogError("") << "[MuonHltTreeProducer]: Trigger collections do not exist !!!";
    }
  }
  

  // 'Hardest' primary vertex (PV) and #PVs reconstructed in the event
  for (unsigned int ix=0; ix<3; ++ix) {
    event_.primaryVertex[ix] = 0.;
    for (unsigned int iy=0; iy<3; ++iy) {
      event_.cov_primaryVertex[ix][iy] = 0.;
    }
  }
  event_.nvvertex = -1;
  
  edm::Handle<std::vector<reco::Vertex> > vertexCollection;

  if(primaryVertexTag_.label() != "none") {
    if (ev.getByLabel(primaryVertexTag_, vertexCollection)) {
      
      int nvvertex = 0;
      unsigned int vertexCollectionSize = vertexCollection->size();
      
      for (unsigned int i=0; i<vertexCollectionSize; ++i) {
	const reco::Vertex& vertex = vertexCollection->at(i);
	if (!vertex.isValid()) continue;
	++nvvertex;
	if (i==0) {
	  event_.primaryVertex[0] = vertex.x();
	  event_.primaryVertex[1] = vertex.y();
	  event_.primaryVertex[2] = vertex.z();
	  for (unsigned int ix=0; ix<3; ++ix) {
	    for (unsigned int iy=0; iy<3; ++iy) {
	      event_.cov_primaryVertex[ix][iy] = vertex.covariance(ix,iy);
	    }
	  }
	}
      }
      event_.nvvertex = nvvertex;
    }
    else {
      edm::LogError("") << "[MuonHltTreeProducer]: Vertex collection does not exist !!!";
    }
  }


  // Beam spot (used only if there is no primary vertex found)
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  if (beamSpotTag_.label() != "none" ) { 
    if (!ev.getByLabel(beamSpotTag_, beamSpotHandle)) {
      edm::LogError("") << "[MuonHltTreeProducer]: Beam spot collection not found !!!";
    }
  }


  ///////////////////////////////////////////////////////
  // Muon collection
  edm::Handle<reco::MuonCollection> muonCollection;
  if (muonTag_.label() != "none" ) { 
    if (!ev.getByLabel(muonTag_, muonCollection)) {
      edm::LogError("") << "[MuonHltTreeProducer] Muon collection does not exist !!!";
    }
  }
  
  if ( muonCollection.isValid()   && 
       vertexCollection.isValid() &&
       beamSpotHandle.isValid() ) {
    fillMuons( muonCollection, vertexCollection, beamSpotHandle, event_.muons);
  }
  
  // GenParticles
  if (!ev.isRealData()) {
    edm::Handle<reco::GenParticleCollection> genParticles;
    if (genTag_.label() != "none" ) { 
      if (ev.getByLabel(genTag_, genParticles)) {

	unsigned int gensize = genParticles->size();
    
	// Do not record the initial protons
	for (unsigned int i=0; i<gensize; ++i) {
	  const reco::GenParticle& part = genParticles->at(i);
	  
	  ciemat::GenParticle gensel;
	  gensel.pdgId = part.pdgId();
	  gensel.status = part.status();
	  gensel.energy = part.energy();
	  gensel.pt = part.pt();
	  gensel.eta = part.eta();
	  gensel.phi = part.phi();
	  gensel.vx = part.vx();
	  gensel.vy = part.vy();
	  gensel.vz = part.vz();
	  gensel.mothers.clear();
	  for (unsigned int im=0; im<part.numberOfMothers(); ++im) {
	    gensel.mothers.push_back(part.motherRef(im)->pdgId());
	  }
	  // Protect agains bug in genParticles (missing mother => first proton)
	  if (i>=2 && part.numberOfMothers()==0) gensel.mothers.push_back(0);
	  
	  event_.genParticles.push_back(gensel);
	}
      }
      else {
	edm::LogError("") << ">>> GEN collection does not exist !!!";
      }
    }
  }
  
  htree_["hMuHltTree"]->Fill();
  
}


void fillMuons(const edm::Handle<reco::MuonCollection> & muonCollection,
	       const edm::Handle<std::vector<reco::Vertex> > & vertexCollection,
	       const edm::Handle<reco::BeamSpot> & beamspot,
	       std::vector<ciemat::Muon> & muons)
{
  
  size_t muonCollectionSize = muonCollection->size();
  
  for (unsigned int i=0; i<muonCollectionSize; ++i) {
    
    const reco::Muon& mu = muonCollection->at(i);

    bool isGlobal = mu.isGlobalMuon();
    bool hasInnerTrack = !mu.innerTrack().isNull();

    double dxy = isGlobal ? mu.globalTrack()->dxy(vertexCollection->at(0).position()) :
      hasInnerTrack ? mu.innerTrack()->dxy(vertexCollection->at(0).position()) : -1000;
    double dz  = isGlobal ? mu.globalTrack()->dz(vertexCollection->at(0).position()) :
      hasInnerTrack ? mu.innerTrack()->dz(vertexCollection->at(0).position()) : -1000;
    
    double dxybs = isGlobal ? mu.globalTrack()->dxy(beamspot->position()) :
      hasInnerTrack ? mu.innerTrack()->dxy(beamspot->position()) : -1000;
    double dzbs  = isGlobal ? mu.globalTrack()->dz(beamspot->position()) :
      hasInnerTrack ? mu.innerTrack()->dz(beamspot->position()) : -1000;

    ciemat::Muon musel;

    musel.pt  = mu.pt();
    musel.eta = mu.eta();
    musel.phi = mu.phi();

    reco::MuonPFIsolation iso04 = mu.pfIsolationR04();
    reco::MuonPFIsolation iso03 = mu.pfIsolationR03();

    musel.chargedHadronIso = iso04.sumChargedHadronPt;
    musel.neutralHadronIso = iso04.sumNeutralHadronEt;
    musel.photonIso        = iso04.sumPhotonEt;

    
    musel.isGlobal = isGlobal ? 1: 0;
	
    musel.isLoose  = muon::isLooseMuon(mu) ? 1 : 0;	  
    musel.isSoft   = muon::isSoftMuon(mu,vertexCollection->at(0)) ? 1 : 0;	  
    musel.isTight  = muon::isTightMuon(mu,vertexCollection->at(0)) ? 1 : 0;	  
    musel.isHighPt = muon::isHighPtMuon(mu,vertexCollection->at(0)) ? 1 : 0;	  
    
    musel.charge = mu.charge();

    musel.dxy    = dxy;
    musel.dz     = dz;
    musel.edxy   = isGlobal ? mu.globalTrack()->dxyError() : hasInnerTrack ? mu.innerTrack()->dxyError() : -1000;
    musel.edz    = isGlobal ? mu.globalTrack()->dzError()  : hasInnerTrack ? mu.innerTrack()->dzError() : -1000;
    musel.dxybs  = dxybs;
    musel.dzbs   = dzbs;
    
    musel.iso_pflow = (iso04.sumChargedHadronPt+ std::max(0.,iso04.sumPhotonEt+iso04.sumNeutralHadronEt - 0.5*iso04.sumPUPt)) / mu.pt();
    
    musel.iso03_pflow = (iso03.sumChargedHadronPt+ std::max(0.,iso03.sumPhotonEt+iso03.sumNeutralHadronEt - 0.5*iso03.sumPUPt)) / mu.pt();
    
    muons.push_back(musel);

  }

}


// void fillGenPart( const reco::GenParticle * gpart, ciemat::GenParticle gensel )
// {
  
//   gensel.pdgId = gpart->pdgId();
//   gensel.status = gpart->status();
//   gensel.energy = gpart->energy();
//   gensel.pt = gpart->pt();
//   gensel.eta = gpart->eta();
//   gensel.phi = gpart->phi();
//   gensel.vx = gpart->vx();
//   gensel.vy = gpart->vy();
//   gensel.vz = gpart->vz();
//   gensel.mothers.clear();
//   for (unsigned int im=0; im<gpart->numberOfMothers(); ++im) {
//     gensel.mothers.push_back(gpart->motherRef(im)->pdgId());
//   }
//   // Protect agains bug in genParticles (missing mother => first proton)
//   if ( gpart->numberOfMothers() == 0 ) gensel.mothers.push_back(0);

// }

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MuonHltTreeProducer);
