import FWCore.ParameterSet.Config as cms

resolvedEMuCRanalyzer = cms.EDAnalyzer("ResolvedEMuCRanalyzer",
  isMC = cms.bool(True),
  generator = cms.InputTag("generator"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  genptc = cms.InputTag("prunedGenParticles"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  METfilters = cms.InputTag("TriggerResults","","PAT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
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
    "HLT_Mu50_v*",
    "HLT_TkMu50_v*"
  ),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016bUL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2016_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2016_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2016_schemaV2.json"),
  muonBoostIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/muIso20UL16.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  recoEleSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016postVFP.root"),
  RMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL16.root"),
  REFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/REFF_20UL16.root"),
  ptThresTrig = cms.double(52.),
  ptThres = cms.double(20.),
  ffSystCL95 = cms.double(0.5),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2016postVFP.root"),
  modHeepSFmuEB1 = cms.double(0.992),
  modHeepSFmuEB2 = cms.double(1.010),
  modHeepSFmuEE = cms.double(1.035),
  modHeepSFcl95EB1 = cms.double(0.03516),
  modHeepSFcl95EB2 = cms.double(0.03548),
  modHeepSFcl95EE = cms.double(0.06211),
  modHeepSFupperEB1 = cms.double(1.098),
  modHeepSFupperEB2 = cms.double(1.142),
  modHeepSFupperEE = cms.double(1.344),
  modHeepPolEB1str = cms.string("-7.294*0.000001*x+0.995"),
  modHeepPolEB2str = cms.string("-6.237*0.000001*x+1.012"),
  modHeepPolEEstr = cms.string("7.755*0.00001*x+1.026"),
  muScaleBias = cms.vdouble(0.002, 0.049, -0.045, -0.061, 0.016, -0.108,
                            -0.005, -0.064, 0.023, -0.004, 0.041, -0.023,
                            0.026, -0.142, 0.039, 0.000, 0.003, -0.091 ),
  muSmearFactors = cms.vdouble(0.,0.),
  muSmearParams = cms.vdouble(0.0102, 6.77e-05, -3.72e-08, 8.53e-12, 0.0129, 6.48e-05, -3.04e-08, 6.63e-12)
)

resolvedEMuCRanalyzer20UL16APV = resolvedEMuCRanalyzer.clone(
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
    "HLT_Mu50_v*",
    "HLT_TkMu50_v*"
  ),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016aUL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2016_preVFP_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2016_preVFP_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2016_preVFP_schemaV2.json"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  recoEleSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016preVFP.root"),
  RMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL16APV.root"),
  REFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/REFF_20UL16APV.root"),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2016preVFP.root"),
  modHeepSFmuEB1 = cms.double(0.995),
  modHeepSFmuEB2 = cms.double(1.003),
  modHeepSFmuEE = cms.double(1.007),
  modHeepSFcl95EB1 = cms.double(0.03915),
  modHeepSFcl95EB2 = cms.double(0.03809),
  modHeepSFcl95EE = cms.double(0.05239),
  modHeepSFupperEB1 = cms.double(1.116),
  modHeepSFupperEB2 = cms.double(1.112),
  modHeepSFupperEE = cms.double(1.227),
  modHeepPolEB1str = cms.string("-1.039*0.00001*x+0.998"),
  modHeepPolEB2str = cms.string("-3.185*0.00001*x+1.012"),
  modHeepPolEEstr = cms.string("7.773*0.000001*x+1.005"),
  muScaleBias = cms.vdouble(-0.075, -0.053, -0.020, -0.009, -0.010, 0.078,
                            -0.049, -0.100, 0.072, 0.057, -0.036, -0.073,
                            -0.007, 0.003, -0.015, -0.003, 0.060, -0.196 ),
  muSmearFactors = cms.vdouble(0.,0.),
  muSmearParams = cms.vdouble(0.011, 6.87e-05, -3.88e-08, 9.03e-12, 0.013, 6.93e-05, -3.46e-08, 7.72e-12)
)

resolvedEMuCRanalyzer20UL17 = resolvedEMuCRanalyzer.clone(
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
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2017UL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2017_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2017_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2017_schemaV2.json"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL17.root"),
  recoEleSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2017.root"),
  RMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL17.root"),
  REFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/REFF_20UL17.root"),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2017.root"),
  modHeepSFmuEB1 = cms.double(1.015),
  modHeepSFmuEB2 = cms.double(1.002),
  modHeepSFmuEE = cms.double(1.003),
  modHeepSFcl95EB1 = cms.double(0.03528),
  modHeepSFcl95EB2 = cms.double(0.02765),
  modHeepSFcl95EE = cms.double(0.02814),
  modHeepSFupperEB1 = cms.double(1.118),
  modHeepSFupperEB2 = cms.double(1.093),
  modHeepSFupperEE = cms.double(1.144),
  modHeepPolEB1str = cms.string("1.119*0.00001*x+1.012"),
  modHeepPolEB2str = cms.string("2.415*0.00001*x+0.995"),
  modHeepPolEEstr = cms.string("2.153*0.00001*x+1.001"),
  muScaleBias = cms.vdouble(-0.049, -0.032, 0.027, -0.026, 0.012, -0.151,
                            -0.048, 0.001, 0.039, 0.027, -0.007, 0.058,
                            -0.025, -0.019, -0.049, 0.018, 0.000, 0.004 ),
  muSmearFactors = cms.vdouble(0.,0.3202),
  muSmearParams = cms.vdouble(0.0104, 6.11e-05, -3.31e-08, 6.73e-12, 0.0121, 5.92e-05, -2.61e-08, 5.11e-12)
)

resolvedEMuCRanalyzer20UL18 = resolvedEMuCRanalyzer.clone(
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
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2018UL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2018_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2018_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2018_schemaV2.json"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL18.root"),
  recoEleSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2018.root"),
  RMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL18.root"),
  REFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/REFF_20UL18.root"),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2018.root"),
  modHeepSFmuEB1 = cms.double(0.998),
  modHeepSFmuEB2 = cms.double(0.991),
  modHeepSFmuEE = cms.double(1.005),
  modHeepSFcl95EB1 = cms.double(0.02489),
  modHeepSFcl95EB2 = cms.double(0.03772),
  modHeepSFcl95EE = cms.double(0.04466),
  modHeepSFupperEB1 = cms.double(1.110),
  modHeepSFupperEB2 = cms.double(1.128),
  modHeepSFupperEE = cms.double(1.125),
  modHeepPolEB1str = cms.string("-8.873*0.000001*x+1.000"),
  modHeepPolEB2str = cms.string("2.09*0.000001*x+0.990"),
  modHeepPolEEstr = cms.string("2.354*0.00001*x+0.999"),
  muScaleBias = cms.vdouble(-0.089, -0.005, -0.014, 0.020, 0.042, -0.012,
                            0.057, -0.028, -0.009, 0.007, -0.068, -0.235,
                            0.149, 0.021, -0.028, 0.000, -0.059, 0.025),
  muSmearFactors = cms.vdouble(0.,0.46),
  muSmearParams = cms.vdouble(0.0108, 5.93e-05, -3.08e-08, 6.04e-12, 0.0136, 5.47e-05, -2.3e-08, 4.66e-12)
)
