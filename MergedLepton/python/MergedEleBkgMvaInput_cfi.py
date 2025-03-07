import FWCore.ParameterSet.Config as cms

mergedEleBkgMvaInput = cms.EDAnalyzer("MergedEleBkgMvaInput",
  srcEle = cms.InputTag("slimmedElectrons"),
  srcGenPtc = cms.InputTag("prunedGenParticles"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  srcGsfTrack = cms.InputTag("reducedEgamma:reducedGsfTracks"),
  srcPackedCand = cms.InputTag("packedPFCandidates"),
  srcLostTracks = cms.InputTag("lostTracks"),
  trkIsoMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleTrkPtIso"),
  ecalIsoMap = cms.InputTag("modifiedEcalRecHitIsolationScone2nd","EcalRecHitIso"),
  dPerpIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","dPerpIn"),
  dEtaInSeed2nd = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","dEtaInSeed2nd"),
  dPhiInSC2nd = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","dPhiInSC2nd"),
  alphaTrack = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","alphaTrack"),
  alphaCalo = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","alphaCalo"),
  normalizedDParaIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","normalizedDParaIn"),
  union5x5covIeIe = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5covIeIe"),
  union5x5covIeIp = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5covIeIp"),
  union5x5covIpIp = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5covIpIp"),
  union5x5dEtaIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5dEtaIn"),
  union5x5dPhiIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5dPhiIn"),
  union5x5Energy = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5Energy"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  addPackedCandMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddPackedCand"),
  conversions = cms.InputTag("reducedEgamma:reducedConversions"),
  generator = cms.InputTag("generator"),
  lheEvent = cms.InputTag("externalLHEProducer"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  EBrecHits = cms.InputTag("reducedEgamma","reducedEBRecHits"),
  EErecHits = cms.InputTag("reducedEgamma","reducedEERecHits"),
  ptThres = cms.double(20.),
  ptThres2nd = cms.double(10.),
  drThres = cms.double(0.1),
  select0J = cms.bool(False),
  selectHT = cms.bool(False),
  maxHT = cms.double(70.),
  posCalcLog = cms.PSet( T0_barl      = cms.double(7.4),
                         T0_endc      = cms.double(3.1),
                         T0_endcPresh = cms.double(1.2),
                         LogWeighted  = cms.bool(True),
                         W0           = cms.double(4.2),
                         X0           = cms.double(0.89) )
)
