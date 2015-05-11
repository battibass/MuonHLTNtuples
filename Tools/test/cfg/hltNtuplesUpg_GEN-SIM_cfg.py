import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLES")

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(
             'file:/afs/cern.ch/user/b/battilan/work/private/GenQCDUpg/TAU-2023MuonUpg14-00003.root'
             ),
        secondaryFileNames = cms.untracked.vstring()
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "DES23_62_V1::All"

process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaper_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023SHCal 
process = cust_2023SHCal(process)

process.muonHltTree = cms.EDAnalyzer("MuonHltTreeProducer",
                             TrigResultsTag = cms.untracked.InputTag(""),
                             TrigSummaryTag = cms.untracked.InputTag(""),

                             MuonTag          = cms.untracked.InputTag(""),
                             PrimaryVertexTag = cms.untracked.InputTag(""),
                             BeamSpotTag      = cms.untracked.InputTag(""),
                             
                             GenTag = cms.untracked.InputTag("genParticles"), # pruned
                             PileUpInfoTag = cms.untracked.InputTag("")
                             )
process.TFileService = cms.Service('TFileService',
        fileName = cms.string("hltNtuple.root")
    )

process.AOutput = cms.EndPath(process.muonHltTree)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
