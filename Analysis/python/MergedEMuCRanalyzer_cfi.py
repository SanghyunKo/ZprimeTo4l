import FWCore.ParameterSet.Config as cms

mergedEMuCRanalyzer = cms.EDAnalyzer("MergedEMuCRanalyzer",
  isMC = cms.bool(True),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  generator = cms.InputTag("generator"),
  genptc = cms.InputTag("prunedGenParticles"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  METfilters = cms.InputTag("TriggerResults","","PAT"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  srcMET = cms.InputTag("slimmedMETs"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  addPackedCandMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddPackedCand"),
  union5x5dEtaIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5dEtaIn"),
  union5x5dPhiIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5dPhiIn"),
  union5x5Energy = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5Energy"),
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter"
  ),
  trigList = cms.vstring(
  # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonHLT
    "HLT_Mu50_v*",
    "HLT_TkMu50_v*"
  ),
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016postVFP.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL16.root"),
  FFpathNonPrompt = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_run2.root"),
  MMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MMFF_20UL16.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL16.json"),
  PUname = cms.string("Collisions16_UltraLegacy_goldenJSON"),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016bUL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2016_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2016_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2016_schemaV2.json"),
  muonBoostIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/muIso20UL16.root"),
  ptThres = cms.double(50.),
  ptMuThres = cms.double(20.),
  drThres = cms.double(0.3),
  drThresCR = cms.double(0.6),
  ratioThresLo = cms.double(0.5),
  ratioThresHi = cms.double(1.5),
  mergedEleSFmuHasTrk = cms.double(1.028),
  mergedEleSFmuNoTrk = cms.double(1.004),
  mergedEleSFcl95HasTrk = cms.double(0.06242),
  mergedEleSFcl95NoTrk = cms.double(0.2448),
  mergedEleSFupperHasTrk = cms.double(1.085),
  mergedEleSFupperNoTrk = cms.double(1.785),
  mergedElePolHasTrkStr = cms.string("0.001304*x+0.959"),
  mergedElePolNoTrkStr = cms.string("0.004586*x+0.919"),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2016postVFP.root"),
  modHeepSFmuEB1 = cms.double(0.986),
  modHeepSFmuEB2 = cms.double(1.002),
  modHeepSFmuEE = cms.double(1.022),
  modHeepSFcl95EB1 = cms.double(0.02226),
  modHeepSFcl95EB2 = cms.double(0.03557),
  modHeepSFcl95EE = cms.double(0.05408),
  modHeepSFupperEB1 = cms.double(1.098),
  modHeepSFupperEB2 = cms.double(1.139),
  modHeepSFupperEE = cms.double(1.344),
  modHeepPolEB1str = cms.string("2.805*0.000001*x+0.986"),
  modHeepPolEB2str = cms.string("3.102*0.000001*x+1.001"),
  modHeepPolEEstr = cms.string("0.0001103*x+1.010"),
  muScaleBias = cms.vdouble(0.002, 0.049, -0.045, -0.061, 0.016, -0.108,
                            -0.005, -0.064, 0.023, -0.004, 0.041, -0.023,
                            0.026, -0.142, 0.039, 0.000, 0.003, -0.091 ),
  muSmearFactors = cms.vdouble(0.,0.),
  muSmearParams = cms.vdouble(0.0102, 6.77e-05, -3.72e-08, 8.53e-12, 0.0129, 6.48e-05, -3.04e-08, 6.63e-12),
  year = cms.string("2016nonAPV"),
  abcdScaleAbove = cms.double(0.9964),
  abcdScaleBelow = cms.double(1.032),
  abcdSmear = cms.double(0.0158),
)

mergedEMuCRanalyzer20UL16APV = mergedEMuCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016preVFP.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL16APV.root"),
  MMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MMFF_20UL16APV.root"),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016aUL.txt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL16APV.json"),
  PUname = cms.string("Collisions16_UltraLegacy_goldenJSON"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2016_preVFP_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2016_preVFP_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2016_preVFP_schemaV2.json"),
  muonBoostIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/muIso20UL16APV.root"),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2016preVFP.root"),
  modHeepSFmuEB1 = cms.double(0.996),
  modHeepSFmuEB2 = cms.double(1.000),
  modHeepSFmuEE = cms.double(1.006),
  modHeepSFcl95EB1 = cms.double(0.03649),
  modHeepSFcl95EB2 = cms.double(0.03954),
  modHeepSFcl95EE = cms.double(0.05728),
  modHeepSFupperEB1 = cms.double(1.117),
  modHeepSFupperEB2 = cms.double(1.113),
  modHeepSFupperEE = cms.double(1.213),
  modHeepPolEB1str = cms.string("-7.931*0.000001*x+0.998"),
  modHeepPolEB2str = cms.string("-3.003*0.00001*x+1.010"),
  modHeepPolEEstr = cms.string("1.211*0.000001*x+1.005"),
  muScaleBias = cms.vdouble(-0.075, -0.053, -0.020, -0.009, -0.010, 0.078,
                            -0.049, -0.100, 0.072, 0.057, -0.036, -0.073,
                            -0.007, 0.003, -0.015, -0.003, 0.060, -0.196 ),
  muSmearFactors = cms.vdouble(0.,0.),
  muSmearParams = cms.vdouble(0.011, 6.87e-05, -3.88e-08, 9.03e-12, 0.013, 6.93e-05, -3.46e-08, 7.72e-12),
  year = cms.string("2016APV"),
  abcdScaleAbove = cms.double(0.9962),
  abcdScaleBelow = cms.double(1.024),
  abcdSmear = cms.double(0.0175),
)

mergedEMuCRanalyzer20UL17 = mergedEMuCRanalyzer.clone(
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_ecalBadCalibFilter"
  ),
  trigList = cms.vstring(
    "HLT_Mu50_v*",
    "HLT_OldMu100_v*",
    "HLT_TkMu100_v*"
  ),
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2017.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL17.root"),
  MMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MMFF_20UL17.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL17.json"),
  PUname = cms.string("Collisions17_UltraLegacy_goldenJSON"),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2017UL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2017_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2017_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2017_schemaV2.json"),
  muonBoostIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/muIso20UL17.root"),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2017.root"),
  modHeepSFmuEB1 = cms.double(0.986),
  modHeepSFmuEB2 = cms.double(0.993),
  modHeepSFmuEE = cms.double(0.996),
  modHeepSFcl95EB1 = cms.double(0.01585),
  modHeepSFcl95EB2 = cms.double(0.0199),
  modHeepSFcl95EE = cms.double(0.02101),
  modHeepSFupperEB1 = cms.double(1.117),
  modHeepSFupperEB2 = cms.double(1.092),
  modHeepSFupperEE = cms.double(1.142),
  modHeepPolEB1str = cms.string("2.752*0.00001*x+0.982"),
  modHeepPolEB2str = cms.string("2.542*0.00001*x+0.988"),
  modHeepPolEEstr = cms.string("9.058*0.00001*x+0.986"),
  muScaleBias = cms.vdouble(-0.049, -0.032, 0.027, -0.026, 0.012, -0.151,
                            -0.048, 0.001, 0.039, 0.027, -0.007, 0.058,
                            -0.025, -0.019, -0.049, 0.018, 0.000, 0.004 ),
  muSmearFactors = cms.vdouble(0.,0.3202),
  muSmearParams = cms.vdouble(0.0104, 6.11e-05, -3.31e-08, 6.73e-12, 0.0121, 5.92e-05, -2.61e-08, 5.11e-12),
  year = cms.string("2017"),
  abcdScaleAbove = cms.double(0.9952),
  abcdScaleBelow = cms.double(1.03),
  abcdSmear = cms.double(0.0165),
)

mergedEMuCRanalyzer20UL18 = mergedEMuCRanalyzer.clone(
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_ecalBadCalibFilter"
  ),
  trigList = cms.vstring(
    "HLT_Mu50_v*",
    "HLT_OldMu100_v*",
    "HLT_TkMu100_v*"
  ),
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2018.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL18.root"),
  MMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MMFF_20UL18.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL18.json"),
  PUname = cms.string("Collisions18_UltraLegacy_goldenJSON"),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2018UL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2018_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2018_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2018_schemaV2.json"),
  muonBoostIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/muIso20UL18.root"),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2018.root"),
  modHeepSFmuEB1 = cms.double(0.980),
  modHeepSFmuEB2 = cms.double(0.983),
  modHeepSFmuEE = cms.double(0.998),
  modHeepSFcl95EB1 = cms.double(0.01631),
  modHeepSFcl95EB2 = cms.double(0.02315),
  modHeepSFcl95EE = cms.double(0.0391),
  modHeepSFupperEB1 = cms.double(1.111),
  modHeepSFupperEB2 = cms.double(1.126),
  modHeepSFupperEE = cms.double(1.126),
  modHeepPolEB1str = cms.string("4.257*0.00001*x+0.977"),
  modHeepPolEB2str = cms.string("1.915*0.00001*x+0.981"),
  modHeepPolEEstr = cms.string("2.257*0.00001*x+0.992"),
  muScaleBias = cms.vdouble(-0.089, -0.005, -0.014, 0.020, 0.042, -0.012,
                            0.057, -0.028, -0.009, 0.007, -0.068, -0.235,
                            0.149, 0.021, -0.028, 0.000, -0.059, 0.025),
  muSmearFactors = cms.vdouble(0.,0.46),
  muSmearParams = cms.vdouble(0.0108, 5.93e-05, -3.08e-08, 6.04e-12, 0.0136, 5.47e-05, -2.3e-08, 4.66e-12),
  year = cms.string("2018"),
  abcdScaleAbove = cms.double(0.9965),
  abcdScaleBelow = cms.double(1.03),
  abcdSmear = cms.double(0.0165),
)
