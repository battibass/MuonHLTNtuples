import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLES")

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(
             # Z 710pre9 RelVal
             '/store/user/battilan/data/62X_RAW_RECO.root'
             ),
        secondaryFileNames = cms.untracked.vstring()
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "POSTLS170_V3::All"

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

from MuonHLTNtuples.Tools.MuonHltNtuples_cff import appendMuonHltNtuple

appendMuonHltNtuple(process,True,"test.root")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
