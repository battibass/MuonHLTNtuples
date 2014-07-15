import FWCore.ParameterSet.Config as cms


def appendMuonHltNtuple(process, runOnMC, processTag="HLT", ntupleFileName="MuonHltTree.root") :

    process.load("CMSSWCiemat.Top2012.MuonHltTreeProducer_cfi")

    process.MuonHltTree.PileUpInfoTag = "none" # CB hack for now

    if processTag != "HLT" :
        print "[MuonHltNtuples]: Customising process tag for TriggerResults / Summary to :", processTag
        process.MuonHltTree.TrigResultsTag = "TriggerResults::"+processTag
        process.MuonHltTree.TrigSummaryTag = "hltTriggerSummaryAOD::"+processTag

    if runOnMC :
        process.load("CMSSWCiemat.Top2012.PrunedGenParticles_cfi")
        process.muonHltNtuple = cms.Sequence(process.prunedGenParticles + process.MuonHltTree)
    else :
        process.muonHltNtuple = cms.Sequence(process.MuonHltTree)

    process.TFileService = cms.Service('TFileService',
        fileName = cms.string(ntupleFileName)
    )

    if hasattr(process,"AOutput") :
        print "[MuonHltNtuples]: EndPath AOutput found, appending ntuples"
        process.AOutput.replace(process.hltOutputA, process.hltOutputA + process.muonHltNtuple)
    else :
        print "[MuonHltNtuples]: EndPath AOuptput not found, creating it for ntuple sequence"
        process.AOutput = cms.EndPath(process.muonHltNtuple)
