import FWCore.ParameterSet.Config as cms

MuonHltTree = cms.EDAnalyzer("MuonHltTreeProducer",
                             TrigResultsTag = cms.untracked.InputTag("TriggerResults::HLT"),
                             TrigSummaryTag = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),

                             MuonTag          = cms.untracked.InputTag("muons"),
                             PrimaryVertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
                             BeamSpotTag      = cms.untracked.InputTag("offlineBeamSpot"),
                             
                             GenTag = cms.untracked.InputTag("genParticles"), #CB do we need pruned?
                             PileUpInfoTag = cms.untracked.InputTag("pileupInfo")
                             )



