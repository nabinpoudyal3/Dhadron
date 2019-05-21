import FWCore.ParameterSet.Config as cms

process = cms.Process("reClusterONE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    "file:/uscms_data/d3/npoudyal/HCC/Dataset/h2cc_miniAOD.root",
  	# "root://cmsxrootd.fnal.gov//store/mc/RunIISummer17MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/060B1DAB-6FAA-E711-8F2B-0CC47A4D762A.root",
    )
)
from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
jetToolbox( process, 'ak1', 'ak1JetSubs', 'noOutput',
  PUMethod ='CHS',  
  runOnMC=True,miniAOD=True,
  # Cut ='pt > 10. && abs(eta) < 2.1',              
  JETCorrPayload = 'AK4PFchs', JETCorrLevels = ['L1FastJet','L2Relative', 'L3Absolute']
)

# jetToolbox( process, 'ak2', 'ak2JetSubs', 'noOutput',
#   PUMethod ='CHS',       
#   Cut ='pt > 10. && abs(eta) < 2.1',              
#   JETCorrPayload = 'AK4PFchs', JETCorrLevels = ['L1FastJet','L2Relative', 'L3Absolute']
# )

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
  triggerConditions = cms.vstring(
    "HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v5",  
    "HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg_v5",  
    "HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg_v5",  
    "HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg_v8",  
    "HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg_v8",  
    "HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v8",  
    "HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg_v8",  
    "HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg_v8",  
    "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v8",  
    "HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v8",  
    "HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v8",  
    "HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v8",  
    "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v8",  
    "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v8",  
    "HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v8",                                   
  ),
  hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
  l1tResults = cms.InputTag( "" ),
  throw = cms.bool(False)
)

process.reClusterONE = cms.EDAnalyzer('Dhadron',
    vertices  = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons     = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus      = cms.InputTag("slimmedTaus"),
    photons   = cms.InputTag("slimmedPhotons"),
    AK4CHS    = cms.InputTag("slimmedJets"),
    fatjets   = cms.InputTag("slimmedJetsAK8"),
    mets      = cms.InputTag("slimmedMETs"),
    pfCands   = cms.InputTag("packedPFCandidates"),   
)
process.TFileService = cms.Service("TFileService",
      fileName = cms.string("ReClusteredHist_May21_test.root"),
      closeFileFast = cms.untracked.bool(True)
)
# Uncomment the following line if you would like to output the jet collections in a root file
# process.endpath = cms.EndPath(process.out)
# process.p = cms.Path(process.triggerSelection*process.reClusterONE)
process.p = cms.Path(process.reClusterONE)

