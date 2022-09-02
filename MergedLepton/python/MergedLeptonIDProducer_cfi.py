import FWCore.ParameterSet.Config as cms

mergedLeptonIDProducer = cms.EDProducer("MergedLeptonIDProducer",
  srcEle=cms.InputTag("slimmedElectrons"),
  srcPv=cms.InputTag("offlineSlimmedPrimaryVertices"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  trkIsoMap=cms.InputTag("ModifiedHEEPIDVarValueMaps","eleTrkPtIso"),
  ecalIsoMap=cms.InputTag("ModifiedEcalRecHitIsolationScone","EcalRecHitIso"),
  nrSatCrysMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleNrSaturateIn5x5"),
  addGsfTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrk"),
  rho=cms.InputTag("fixedGridRhoFastjetAll"),
  xgbPathDR1EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/dr1_EB.xml"),
  meanstdPathDR1EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/dr1_EB.csv"),
  cutDR1EB=cms.double(0.9),
  xgbPathDR2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/dr2_EB.xml"),
  meanstdPathDR2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/dr2_EB.csv"),
  cutDR2EB=cms.double(0.9),
  xgbPathDR3EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/dr3_EB.xml"),
  meanstdPathDR3EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/dr3_EB.csv"),
  cutDR3EB=cms.double(0.9),
  xgbPathDR1EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/dr1_EE.xml"),
  meanstdPathDR1EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/dr1_EE.csv"),
  cutDR1EE=cms.double(0.9),
  xgbPathDR2EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/dr2_EE.xml"),
  meanstdPathDR2EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/dr2_EE.csv"),
  cutDR2EE=cms.double(0.9),
  xgbPathDR3EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/dr3_EE.xml"),
  meanstdPathDR3EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/dr3_EE.csv"),
  cutDR3EE=cms.double(0.9)
)
