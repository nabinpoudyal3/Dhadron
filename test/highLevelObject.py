import FWCore.ParameterSet.Config as cms

process = cms.Process("hlo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames  = cms.untracked.vstring(
    "file:/uscms_data/d3/npoudyal/HCC/Dataset/h2cc_miniAOD.root",           
    )
)

process.hlo  = cms.EDAnalyzer("Dhadron",

    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons    = cms.InputTag("slimmedMuons"),
    electrons= cms.InputTag("slimmedElectrons"),
    taus     = cms.InputTag("slimmedTaus"),
    photons  = cms.InputTag("slimmedPhotons"),
    jets     = cms.InputTag("slimmedJets"),
    fatjets  = cms.InputTag("slimmedJetsAK8"),
    mets     = cms.InputTag("slimmedMETs"),
    pfCands  = cms.InputTag("packedPFCandidates"),
)

process.p = cms.Path(process.hlo)