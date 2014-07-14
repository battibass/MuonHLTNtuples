import FWCore.ParameterSet.Config as cms


def appendMuonHltNtuple(process, runOnMC, ntupleFileName="MuonHltTree.root") :

    process.load("CMSSWCiemat.Top2012.MuonHltTreeProducer_cfi")

    process.MuonHltTree.PileUpInfoTag = "none" # CB hack for now 

    if runOnMC :
        process.load("CMSSWCiemat.Top2012.PrunedGenParticles_cfi")
        process.muonHltNtuple = cms.Sequence(process.prunedGenParticles + process.MuonHltTree)
    else :
        process.muonHltNtuple = cms.Sequence(process.MuonHltTree)

    process.TFileService = cms.Service('TFileService',
        fileName = cms.string(ntupleFileName)
    )

    if hasattr(process,"out") :
        print "[MuonHltNtuples]: EndPath out found, appending ntuples"
        process.out.append(muonHltNtuple)
    else :
        print "[MuonHltNtuples]: EndPath out not found, creating it for ntuple sequence"
        process.out = cms.EndPath(process.muonHltNtuple)
