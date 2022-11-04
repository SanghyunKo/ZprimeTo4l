import FWCore.ParameterSet.Config as cms

mergedEleBkgMvaInput = cms.EDAnalyzer("MergedEleBkgMvaInput",
  srcEle=cms.InputTag("slimmedElectrons"),
  srcPv=cms.InputTag("offlineSlimmedPrimaryVertices"),
  trkIsoMap=cms.InputTag("ModifiedHEEPIDVarValueMaps","eleTrkPtIso"),
  ecalIsoMap=cms.InputTag("ModifiedEcalRecHitIsolationScone","EcalRecHitIso"),
  nrSatCrysMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleNrSaturateIn5x5"),
  addGsfTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrk"),
  rho=cms.InputTag("fixedGridRhoFastjetAll"),
  conversions = cms.InputTag("reducedEgamma:reducedConversions"),
  generator = cms.InputTag("generator"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  ptThres=cms.double(20.),
  drThres=cms.double(0.3)
)
