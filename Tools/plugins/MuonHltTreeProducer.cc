//////////////////////////////////////
// Ntuplizer that fills muon_hlt trees
//////////////////////////////////////

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "MuonHLTNtuples/Tools/src/MuTree.h"
#include "TTree.h"

#include <algorithm>
#include <iostream>

class MuonHltTreeProducer : public edm::EDAnalyzer 
{
public:

  MuonHltTreeProducer(const edm::ParameterSet &);
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
  
private:
  
  void fillGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > &);

  void fillGenParticles(const edm::Handle<reco::GenParticleCollection> &);

  void fillHlt(const edm::Handle<edm::TriggerResults> &, 
	       const edm::Handle<trigger::TriggerEvent> &,
	       const edm::TriggerNames &);
  
  void fillPV(const edm::Handle<std::vector<reco::Vertex> > &);
  
  void fillMuons(const edm::Handle<reco::MuonCollection> &,
		 const edm::Handle<std::vector<reco::Vertex> > &,
		 const edm::Handle<reco::BeamSpot> &);
  
  edm::InputTag trigResultsTag_;
  edm::InputTag trigSummaryTag_;

  edm::InputTag muonTag_;
  edm::InputTag primaryVertexTag_;
  edm::InputTag beamSpotTag_;

  edm::InputTag genTag_;
  edm::InputTag pileUpInfoTag_;

  muon_hlt::Event event_;
  std::map<std::string,TTree*> tree_;
  
};


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


void MuonHltTreeProducer::beginJob() 
{
  
  edm::Service<TFileService> fs;
  tree_["muHltTree"] = fs->make<TTree>("MUHLTTREE","HLT Muon Tree");

  int splitBranches = 2;
  tree_["muHltTree"]->Branch("event",&event_,64000,splitBranches);

}


void MuonHltTreeProducer::beginRun(const edm::Run & run, const edm::EventSetup & config )
{
  
}


void MuonHltTreeProducer::endJob() 
{

}


void MuonHltTreeProducer::analyze (const edm::Event & ev, const edm::EventSetup &)
{

  // Clearing branch variables
  // and setting default values
  event_.hlt.triggers.clear();
  event_.hlt.objects.clear();

  event_.genParticles.clear();
  event_.genInfos.clear();
  event_.muons.clear();
  
  for (unsigned int ix=0; ix<3; ++ix) {
    event_.primaryVertex[ix] = 0.;
    for (unsigned int iy=0; iy<3; ++iy) {
      event_.cov_primaryVertex[ix][iy] = 0.;
    }
  }
  event_.nVtx = -1;


  // Fill general information
  // run, luminosity block, event
  event_.runNumber = ev.id().run();
  event_.luminosityBlockNumber = ev.id().luminosityBlock();
  event_.eventNumber = ev.id().event();


  // Fill GEN pile up information
  if (!ev.isRealData()) 
    {
      if (pileUpInfoTag_.label() != "none") 
	{
	  edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
	  if (ev.getByLabel(pileUpInfoTag_, puInfo)) 
	    fillGenInfo(puInfo);
	  else 
	    edm::LogError("") << "[MuonHltTreeProducer]: Pile-Up Info collection does not exist !!!";
	}      
    }
  

  // Fill GEN particles information
  if (!ev.isRealData()) 
    {
      if (genTag_.label() != "none" ) 
	{ 
	  edm::Handle<reco::GenParticleCollection> genParticles;
	  if (ev.getByLabel(genTag_, genParticles)) 
	    fillGenParticles(genParticles);
	  else 
	    edm::LogError("") << ">>> GEN collection does not exist !!!";
	}
    }
  

  // Fill trigger information
  if (trigResultsTag_.label() != "none" &&
      trigSummaryTag_.label() != "none") 
    {
      
      edm::Handle<edm::TriggerResults> triggerResults;
      edm::Handle<trigger::TriggerEvent> triggerEvent;
      
      if (ev.getByLabel(trigResultsTag_, triggerResults) &&
	  ev.getByLabel(trigSummaryTag_, triggerEvent)) 
	fillHlt(triggerResults, triggerEvent,ev.triggerNames(*triggerResults));
      else 
	edm::LogError("") << "[MuonHltTreeProducer]: Trigger collections do not exist !!!";
    }
  
  
  // Fill vertex information
  edm::Handle<std::vector<reco::Vertex> > vertexes;

  if(primaryVertexTag_.label() != "none") 
    {
      if (ev.getByLabel(primaryVertexTag_, vertexes))
	fillPV(vertexes);
      else 
	edm::LogError("") << "[MuonHltTreeProducer]: Vertex collection does not exist !!!";
    }
  

  // Get beam spot for muons
  edm::Handle<reco::BeamSpot> beamSpot;
  if (beamSpotTag_.label() != "none" ) 
    { 
      if (!ev.getByLabel(beamSpotTag_, beamSpot)) 
	edm::LogError("") << "[MuonHltTreeProducer]: Beam spot collection not found !!!";
    }


  // Get muons  
  edm::Handle<reco::MuonCollection> muons;
  if (muonTag_.label() != "none" ) 
    { 
      if (!ev.getByLabel(muonTag_, muons)) 
	edm::LogError("") << "[MuonHltTreeProducer] Muon collection does not exist !!!";
    }
  

  // Fill muon information
  if (muons.isValid() && vertexes.isValid() && beamSpot.isValid()) 
    {
      fillMuons(muons,vertexes,beamSpot);
    }
  
  tree_["muHltTree"]->Fill();
  
}


void MuonHltTreeProducer::fillGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > & puInfo)
{

  muon_hlt::GenInfo genInfo;
  
  genInfo.trueNumberOfInteractions   = -1.;
  genInfo.actualNumberOfInteractions = -1 ;    
  
  std::vector<PileupSummaryInfo>::const_iterator puInfoIt  = puInfo->begin();
  std::vector<PileupSummaryInfo>::const_iterator puInfoEnd = puInfo->end();

  for(; puInfoIt != puInfoEnd; ++puInfoIt) 
    {
    
      int bx = puInfoIt->getBunchCrossing();
	  
      if(bx == 0) 
	{ 
	  genInfo.trueNumberOfInteractions   = puInfoIt->getTrueNumInteractions();
	  genInfo.actualNumberOfInteractions = puInfoIt->getPU_NumInteractions();
	  continue;
	}
    }
  
  event_.genInfos.push_back(genInfo);
  
}


void MuonHltTreeProducer::fillGenParticles(const edm::Handle<reco::GenParticleCollection> & genParticles)
{
  
  unsigned int gensize = genParticles->size();
  
  // Do not record the initial protons
  for (unsigned int i=0; i<gensize; ++i) 
    {

      const reco::GenParticle& part = genParticles->at(i);
    
      muon_hlt::GenParticle gensel;
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
      unsigned int nMothers = part.numberOfMothers();

      for (unsigned int iMother=0; iMother<nMothers; ++iMother) 
	{
	  gensel.mothers.push_back(part.motherRef(iMother)->pdgId());
	}

      // Protect agains bug in genParticles (missing mother => first proton)
      if (i>=2 && nMothers==0) gensel.mothers.push_back(0);
      
      event_.genParticles.push_back(gensel);
    }
  
}


void MuonHltTreeProducer::fillHlt(const edm::Handle<edm::TriggerResults> & triggerResults, 
				  const edm::Handle<trigger::TriggerEvent> & triggerEvent,
				  const edm::TriggerNames & triggerNames)
{    

  for (unsigned int iTrig=0; iTrig<triggerNames.size(); ++iTrig) 
    {
      
      if (triggerResults->accept(iTrig)) 
	{
	  std::string pathName = triggerNames.triggerName(iTrig);
	  event_.hlt.triggers.push_back(pathName);
	}
    }
      
  const trigger::size_type nFilters(triggerEvent->sizeFilters());

  for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter) 
    {
	
      std::string filterTag = triggerEvent->filterTag(iFilter).encode();

      trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
      const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
	
      for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey) 
	{  
	  trigger::size_type objKey = objectKeys.at(iKey);
	  const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);
	  
	  muon_hlt::HLTObject hltObj;
	  
	  float trigObjPt = triggerObj.pt();
	  float trigObjEta = triggerObj.eta();
	  float trigObjPhi = triggerObj.phi();
	  
	  hltObj.filterTag = filterTag;

	  hltObj.pt  = trigObjPt;
	  hltObj.eta = trigObjEta;
	  hltObj.phi = trigObjPhi;
	  
	  event_.hlt.objects.push_back(hltObj);
	  
	}       
    }

}


void MuonHltTreeProducer::fillPV(const edm::Handle<std::vector<reco::Vertex> > & vertexes)
{
      
  int nVtx = 0;

  std::vector<reco::Vertex>::const_iterator vertexIt  = vertexes->begin();
  std::vector<reco::Vertex>::const_iterator vertexEnd = vertexes->end();

  for (; vertexIt != vertexEnd; ++vertexIt) 
    {

      const reco::Vertex& vertex = *vertexIt;

      if (!vertex.isValid()) continue;
      ++nVtx;

      if (vertexIt == vertexes->begin()) 
	{
	  event_.primaryVertex[0] = vertex.x();
	  event_.primaryVertex[1] = vertex.y();
	  event_.primaryVertex[2] = vertex.z();

	  for (unsigned int ix=0; ix<3; ++ix) 
	    {
	      for (unsigned int iy=0; iy<3; ++iy) 
		{
		  event_.cov_primaryVertex[ix][iy] = vertex.covariance(ix,iy);
		}
	    }
	}
    }
  
  event_.nVtx = nVtx;
  
}


void MuonHltTreeProducer::fillMuons(const edm::Handle<reco::MuonCollection> & muons,
				    const edm::Handle<std::vector<reco::Vertex> > & vertexes,
				    const edm::Handle<reco::BeamSpot> & beamSpot)
{

  reco::MuonCollection::const_iterator muonIt  = muons->begin();
  reco::MuonCollection::const_iterator muonEnd = muons->end();

  for (; muonIt != muonEnd; ++muonIt) 
    {
      
      const reco::Muon& mu = (*muonIt);
      const reco::Vertex & vertex = vertexes->at(0); // CB for now vertex is always valid, but add a protection	    

      bool isGlobal      = mu.isGlobalMuon();
      bool isTracker     = mu.isTrackerMuon();
      bool isStandAlone  = mu.isStandAloneMuon();

      bool hasInnerTrack = !mu.innerTrack().isNull();

      double dxy = isGlobal ? mu.globalTrack()->dxy(vertex.position()) :
	hasInnerTrack ? mu.innerTrack()->dxy(vertex.position()) : -1000;
      double dz  = isGlobal ? mu.globalTrack()->dz(vertex.position()) :
	hasInnerTrack ? mu.innerTrack()->dz(vertex.position()) : -1000;
      
      double dxybs = isGlobal ? mu.globalTrack()->dxy(beamSpot->position()) :
	hasInnerTrack ? mu.innerTrack()->dxy(beamSpot->position()) : -1000;
      double dzbs  = isGlobal ? mu.globalTrack()->dz(beamSpot->position()) :
	hasInnerTrack ? mu.innerTrack()->dz(beamSpot->position()) : -1000;
      
      muon_hlt::Muon ntupleMu;
      
      ntupleMu.pt  = mu.pt();
      ntupleMu.eta = mu.eta();
      ntupleMu.phi = mu.phi();
      
      reco::MuonPFIsolation iso04 = mu.pfIsolationR04();
      reco::MuonPFIsolation iso03 = mu.pfIsolationR03();

      ntupleMu.chargedHadronIso = iso04.sumChargedHadronPt;
      ntupleMu.neutralHadronIso = iso04.sumNeutralHadronEt;
      ntupleMu.photonIso        = iso04.sumPhotonEt;

      ntupleMu.isGlobal     = isGlobal ? 1: 0;	
      ntupleMu.isTracker    = isTracker ? 1: 0;	
      ntupleMu.isStandAlone = isStandAlone ? 1: 0;

      ntupleMu.nHitsGlobal     = isGlobal     ? mu.globalTrack()->numberOfValidHits() : -999;	
      ntupleMu.nHitsTracker    = isTracker    ? mu.innerTrack()->numberOfValidHits()  : -999;	
      ntupleMu.nHitsStandAlone = isStandAlone ? mu.outerTrack()->numberOfValidHits()  : -999;
	
      ntupleMu.isLoose  = muon::isLooseMuon(mu)         ? 1 : 0;	  
      ntupleMu.isSoft   = muon::isSoftMuon(mu,vertex)   ? 1 : 0;	  
      ntupleMu.isTight  = muon::isTightMuon(mu,vertex)  ? 1 : 0;	  
      ntupleMu.isHighPt = muon::isHighPtMuon(mu,vertex) ? 1 : 0;	  
    
      ntupleMu.charge = mu.charge();

      ntupleMu.dxy    = dxy;
      ntupleMu.dz     = dz;
      ntupleMu.edxy   = isGlobal ? mu.globalTrack()->dxyError() : hasInnerTrack ? mu.innerTrack()->dxyError() : -1000;
      ntupleMu.edz    = isGlobal ? mu.globalTrack()->dzError()  : hasInnerTrack ? mu.innerTrack()->dzError() : -1000;
      ntupleMu.dxybs  = dxybs;
      ntupleMu.dzbs   = dzbs;
    
      ntupleMu.isoPflow04 = (iso04.sumChargedHadronPt+ std::max(0.,iso04.sumPhotonEt+iso04.sumNeutralHadronEt - 0.5*iso04.sumPUPt)) / mu.pt();
    
      ntupleMu.isoPflow03 = (iso03.sumChargedHadronPt+ std::max(0.,iso03.sumPhotonEt+iso03.sumNeutralHadronEt - 0.5*iso03.sumPUPt)) / mu.pt();
      
      event_.muons.push_back(ntupleMu);

    }

}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MuonHltTreeProducer);
