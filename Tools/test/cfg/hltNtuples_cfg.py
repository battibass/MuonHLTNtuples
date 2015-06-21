import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLES")

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(
             ' /store/relval/CMSSW_7_5_0_pre5/RelValSingleMuPt100_UP15/GEN-SIM-RECO/MCRUN2_75_V5-v1/00000/9C3F56B9-B80B-E511-A5FC-0025905A605E.root'
             ),
        secondaryFileNames = cms.untracked.vstring()
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "MCRUN2_75_V5"

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

from MuonHLTNtuples.Tools.MuonHltNtuples_cff import appendMuonHltNtuple

appendMuonHltNtuple(process,True,"test.root")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
