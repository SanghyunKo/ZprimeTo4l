import FWCore.ParameterSet.Config as cms

mergedLeptonIDProducer = cms.EDAnalyzer("MergedLeptonIDProducer",
  srcEle=cms.InputTag("slimmedElectrons"),
  srcPv=cms.InputTag("offlineSlimmedPrimaryVertices"),
  trkIsoMap=cms.InputTag("ModifiedHEEPIDVarValueMaps","eleTrkPtIso"),
  ecalIsoMap=cms.InputTag("ModifiedEcalRecHitIsolationScone","EcalRecHitIso"),
  nrSatCrysMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleNrSaturateIn5x5"),
  addGsfTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrk"),
  rho=cms.InputTag("fixedGridRhoFastjetAll"),
  xgbPathDR1EB=cms.untracked.FileInPath("ZprimeTo4l/MergedLepton/data/dr1_EB.xml"),
  meanstdPathDR1EB=cms.untracked.FileInPath("ZprimeTo4l/MergedLepton/data/dr1_EB.csv"),
  cutDR1EB=cms.double(0.5),
  xgbPathDR2EB=cms.untracked.FileInPath("ZprimeTo4l/MergedLepton/data/dr2_EB.xml"),
  meanstdPathDR2EB=cms.untracked.FileInPath("ZprimeTo4l/MergedLepton/data/dr2_EB.csv"),
  cutDR2EB=cms.double(0.5),
  xgbPathDR3EB=cms.untracked.FileInPath("ZprimeTo4l/MergedLepton/data/dr3_EB.xml"),
  meanstdPathDR3EB=cms.untracked.FileInPath("ZprimeTo4l/MergedLepton/data/dr3_EB.csv"),
  cutDR3EB=cms.double(0.5),
  xgbPathDR1EE=cms.untracked.FileInPath("ZprimeTo4l/MergedLepton/data/dr1_EE.xml"),
  meanstdPathDR1EE=cms.untracked.FileInPath("ZprimeTo4l/MergedLepton/data/dr1_EE.csv"),
  cutDR1EE=cms.double(0.5),
  xgbPathDR2EE=cms.untracked.FileInPath("ZprimeTo4l/MergedLepton/data/dr2_EE.xml"),
  meanstdPathDR2EE=cms.untracked.FileInPath("ZprimeTo4l/MergedLepton/data/dr2_EE.csv"),
  cutDR2EE=cms.double(0.5),
  xgbPathDR3EE=cms.untracked.FileInPath("ZprimeTo4l/MergedLepton/data/dr3_EE.xml"),
  meanstdPathDR3EE=cms.untracked.FileInPath("ZprimeTo4l/MergedLepton/data/dr3_EE.csv"),
  cutDR3EE=cms.double(0.5)
)
