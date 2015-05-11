import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLES")

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(
             'file:/afs/cern.ch/work/c/calabria/public/step3_1000_1_9NR.root'
             ),
        secondaryFileNames = cms.untracked.vstring()
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "PH1_1K_FB_V3::All"

process.load('Configuration.Geometry.GeometryExtended2019Reco_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2019WithGem 
process = cust_2019WithGem(process)


process.muonHltTree = cms.EDAnalyzer("MuonHltTreeProducer",
                             TrigResultsTag = cms.untracked.InputTag(""),
                             TrigSummaryTag = cms.untracked.InputTag(""),

                             MuonTag          = cms.untracked.InputTag("muons"),
                             PrimaryVertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
                             BeamSpotTag      = cms.untracked.InputTag("offlineBeamSpot"),
                             
                             GenTag = cms.untracked.InputTag("genParticles"), 
                             PileUpInfoTag = cms.untracked.InputTag("addPileupInfo")                                     
                             )
process.TFileService = cms.Service('TFileService',
        fileName = cms.string("hltNtuple.root")
    )

process.AOutput = cms.EndPath(process.muonHltTree)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
