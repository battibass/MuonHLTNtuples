import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLES")

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(
             'file:/afs/cern.ch/user/v/vieri/public/cosmics/SPCosmics_GEN-SIM-RECO_step3_20to100.root'
             ),
        secondaryFileNames = cms.untracked.vstring()
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "START50_V13::All"

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

#from SLHCUpgradeSimulations.Configuration.postLS1Customs import *
#process = customise_HLT( process )

from MuonHLTNtuples.Tools.MuonHltNtuples_cff import appendMuonHltNtuple
appendMuonHltNtuple(process,False,"HLT","ntuple_20to100.root")

process.MuonHltTree.TrigResultsTag = cms.untracked.InputTag("TriggerResults::HLT")
process.MuonHltTree.TrigSummaryTag = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT")

process.MuonHltTree.PrimaryVertexTag = cms.untracked.InputTag("none")
process.MuonHltTree.BeamSpotTag      = cms.untracked.InputTag("none")
process.MuonHltTree.PileUpInfoTag    = cms.untracked.InputTag("none")


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
