import FWCore.ParameterSet.Config as cms

prunedGenParticles = cms.EDProducer("GenParticlePruner",
                                    src = cms.InputTag("genParticles"),
                                    select = cms.vstring("drop *"
                                                         , "keep status = 3"
                                                         , "++keep pdgId=11 & pt>10 & abs(eta)<3"
                                                         , "++keep pdgId=13 & pt>10 & abs(eta)<3"
                                                         )
                                    )
