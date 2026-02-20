#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "ZprimeTo4l/Analysis/interface/MuonCorrectionHelper.h"

#include "TF1.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFitResult.h"

#include "correction.h"

class MergedMuNonIsoCRanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedMuNonIsoCRanalyzer(const edm::ParameterSet&);
  ~MergedMuNonIsoCRanalyzer() override {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  bool isMC_;

  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> genptcToken_;
  const edm::EDGetTokenT<double> prefweight_token;
  const edm::EDGetTokenT<double> prefweightUp_token_;
  const edm::EDGetTokenT<double> prefweightDn_token_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;
  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;

  const edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<pat::MET>> metToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;

  const std::vector<std::string> METfilterList_;
  const std::vector<std::string> trigList_;

  const edm::FileInPath MMFFpath_;
  const edm::FileInPath rochesterPath_;
  const edm::FileInPath triggerSFpath_;
  const edm::FileInPath muonIdIsoSFpath_;
  const edm::FileInPath muonBoostIsoSFpath_;
  const edm::FileInPath muonRecoSFpath_;

  std::unique_ptr<correction::CorrectionSet> purwgt_;
  const std::string puname_;

  const std::string year_;

  const double ptThres_;
  const double ptMuThres_;
  const double drThres_;
  const double drThresCR_;
  const double ratioThresLo_;
  const double ratioThresHi_;
  const double mumass_ = 0.1056583745;

  MuonCorrectionHelper mucorrHelper_;

  std::unique_ptr<TFile> MMFFfile_;
  std::unique_ptr<TF1> ffFunc_;
  TFitResultPtr fitResult_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;

  TTree* interestEvtTree_ = nullptr;
  unsigned int runNo_ = 0;
  unsigned int lumiNo_ = 0;
  unsigned long long evtNo_ = 0;

  TTree* nonIsoSRTree_ = nullptr;
  float wgt_ = 1.;
  float m1pt_ = 0;
  float m1eta_ = std::numeric_limits<float>::max();
  float m1phi_ = std::numeric_limits<float>::max();
  float m1iso_ = -1.;
  int m1global_ = -1;
  float m2pt_ = 0;
  float m2eta_ = std::numeric_limits<float>::max();
  float m2phi_ = std::numeric_limits<float>::max();
  float m2iso_ = -1.;
  int m2global_ = -1;
  float mmpt_ = 0;
  float mmeta_ = std::numeric_limits<float>::max();
  float mmphi_ = std::numeric_limits<float>::max();
  float mmiso_ = -1.;
  float metPt_ = 0.;
  float metPhi_ = 0.;
  float mt_ = 0.;
  float mm1m2_ = 0.;
  float ratioPt_ = 0.;
  float dphiMet_ = 0.;
};

MergedMuNonIsoCRanalyzer::MergedMuNonIsoCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getParameter<bool>("isMC")),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
genptcToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genptc"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
prefweightUp_token_(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbUp"))),
prefweightDn_token_(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbDown"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
metToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
MMFFpath_(iConfig.getParameter<edm::FileInPath>("MMFFpath")),
rochesterPath_(iConfig.getParameter<edm::FileInPath>("rochesterPath")),
triggerSFpath_(iConfig.getParameter<edm::FileInPath>("triggerSF")),
muonIdIsoSFpath_(iConfig.getParameter<edm::FileInPath>("muonIdIsoSFpath")),
muonBoostIsoSFpath_(iConfig.getParameter<edm::FileInPath>("muonBoostIsoSFpath")),
muonRecoSFpath_(iConfig.getParameter<edm::FileInPath>("muonRecoSFpath")),
purwgt_(std::move(correction::CorrectionSet::from_file((iConfig.getParameter<edm::FileInPath>("PUrwgt")).fullPath()))),
puname_(iConfig.getParameter<std::string>("PUname")),
year_(iConfig.getParameter<std::string>("year")),
ptThres_(iConfig.getParameter<double>("ptThres")),
ptMuThres_(iConfig.getParameter<double>("ptMuThres")),
drThres_(iConfig.getParameter<double>("drThres")),
drThresCR_(iConfig.getParameter<double>("drThresCR")),
ratioThresLo_(iConfig.getParameter<double>("ratioThresLo")),
ratioThresHi_(iConfig.getParameter<double>("ratioThresHi")),
mucorrHelper_(rochesterPath_,triggerSFpath_,muonIdIsoSFpath_,muonBoostIsoSFpath_,muonRecoSFpath_) {
  usesResource("TFileService");
}

void MergedMuNonIsoCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  MMFFfile_ = std::make_unique<TFile>(MMFFpath_.fullPath().c_str(),"READ");
  ffFunc_ = std::make_unique<TF1>("ffFunc","[0]",0.,2500.);
  // ffFunc_ = static_cast<TF1*>(MMFFfile_->FindObjectAny("MMFF"));
  fitResult_ = (static_cast<TH1D*>(MMFFfile_->Get("1M_MMFF_numer_rebin")))->Fit(ffFunc_.get(),"RS");

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["totWeightedSum_4M"] = fs->make<TH1D>("totWeightedSum_4M","totWeightedSum_4M",1,0.,1.);
  histo1d_["nPV"] = fs->make<TH1D>("nPV","nPV",99,0.,99.);

  histo1d_["GEN_pt_1st"] = fs->make<TH1D>("GEN_pt_1st",";p_{T};",1000,0.,2000.);
  histo1d_["GEN_pt_2nd"] = fs->make<TH1D>("GEN_pt_2nd",";p_{T};",1000,0.,2000.);
  histo1d_["GEN_pt_3rd"] = fs->make<TH1D>("GEN_pt_3rd",";p_{T};",1000,0.,2000.);
  histo1d_["GEN_pt_4th"] = fs->make<TH1D>("GEN_pt_4th",";p_{T};",1000,0.,2000.);

  histo1d_["cutflow"] = fs->make<TH1D>("cutflow","cutflow",20,0.,20.);

  histo1d_["3M_MM_pt"] = fs->make<TH1D>("3M_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_MM_eta"] = fs->make<TH1D>("3M_MM_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["3M_MM_phi"] = fs->make<TH1D>("3M_MM_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["3M_M1_pt"] = fs->make<TH1D>("3M_M1_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_M1_eta"] = fs->make<TH1D>("3M_M1_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["3M_M1_phi"] = fs->make<TH1D>("3M_M1_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["3M_M2_pt"] = fs->make<TH1D>("3M_M2_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_M2_eta"] = fs->make<TH1D>("3M_M2_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["3M_M2_phi"] = fs->make<TH1D>("3M_M2_phi","Phi;#phi;",128,-3.2,3.2);

  // SR mt histograms
  histo1d_["3M_MET_pt"] = fs->make<TH1D>("3M_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_MET_phi"] = fs->make<TH1D>("3M_MET_phi","Phi;#phi;",128,-3.2,3.2);
  histo1d_["3M_MET_dphi"] = fs->make<TH1D>("3M_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["3M_mt"] = fs->make<TH1D>("3M_mt","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_idUp"] = fs->make<TH1D>("3M_mt_idUp","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_idDn"] = fs->make<TH1D>("3M_mt_idDn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_trigUp"] = fs->make<TH1D>("3M_mt_trigUp","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_trigDn"] = fs->make<TH1D>("3M_mt_trigDn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_recoUp"] = fs->make<TH1D>("3M_mt_recoUp","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_recoDn"] = fs->make<TH1D>("3M_mt_recoDn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_prefireUp"] = fs->make<TH1D>("3M_mt_prefireUp","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_prefireDn"] = fs->make<TH1D>("3M_mt_prefireDn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_PUrwgtUp"] = fs->make<TH1D>("3M_mt_PUrwgtUp","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_PUrwgtDn"] = fs->make<TH1D>("3M_mt_PUrwgtDn","m_{T};m_{T};",500,0.,2500.);

  // denominator region mt histograms
  histo1d_["3M_antiRpt_MM_pt"] = fs->make<TH1D>("3M_antiRpt_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiRpt_MM_pt_xFF"] = fs->make<TH1D>("3M_antiRpt_MM_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiRpt_M1_pt"] = fs->make<TH1D>("3M_antiRpt_M1_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiRpt_M2_pt"] = fs->make<TH1D>("3M_antiRpt_M2_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiRpt_MET_pt"] = fs->make<TH1D>("3M_antiRpt_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiRpt_MET_pt_xFF"] = fs->make<TH1D>("3M_antiRpt_MET_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiRpt_MET_dphi"] = fs->make<TH1D>("3M_antiRpt_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["3M_antiRpt_mt"] = fs->make<TH1D>("3M_antiRpt_mt","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiRpt_mt_xFF"] = fs->make<TH1D>("3M_antiRpt_mt_xFF","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiRpt_mt_xFF_up"] = fs->make<TH1D>("3M_antiRpt_mt_xFF_up","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiRpt_mt_xFF_dn"] = fs->make<TH1D>("3M_antiRpt_mt_xFF_dn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiRpt_MET_ratioPt"] = fs->make<TH1D>("3M_antiRpt_MET_ratioPt","p_{T} ratio;R(p_{T});",200,0.,5.);

  // SR m4l histograms with eta hypothesis
  histo1d_["3M_m4l"] = fs->make<TH1D>("3M_m4l","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_m4l_idUp"] = fs->make<TH1D>("3M_m4l_idUp","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_m4l_idDn"] = fs->make<TH1D>("3M_m4l_idDn","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_m4l_trigUp"] = fs->make<TH1D>("3M_m4l_trigUp","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_m4l_trigDn"] = fs->make<TH1D>("3M_m4l_trigDn","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_m4l_recoUp"] = fs->make<TH1D>("3M_m4l_recoUp","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_m4l_recoDn"] = fs->make<TH1D>("3M_m4l_recoDn","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_m4l_prefireUp"] = fs->make<TH1D>("3M_m4l_prefireUp","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_m4l_prefireDn"] = fs->make<TH1D>("3M_m4l_prefireDn","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_m4l_PUrwgtUp"] = fs->make<TH1D>("3M_m4l_PUrwgtUp","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_m4l_PUrwgtDn"] = fs->make<TH1D>("3M_m4l_PUrwgtDn","m_{4l};m_{4l};",500,0.,2500.);

  // denominator region m4l histograms with eta hypothesis
  histo1d_["3M_antiRpt_m4l"] = fs->make<TH1D>("3M_antiRpt_m4l","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_antiRpt_m4l_xFF"] = fs->make<TH1D>("3M_antiRpt_m4l_xFF","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_antiRpt_m4l_xFF_up"] = fs->make<TH1D>("3M_antiRpt_m4l_xFF_up","m_{4l};m_{4l};",500,0.,2500.);
  histo1d_["3M_antiRpt_m4l_xFF_dn"] = fs->make<TH1D>("3M_antiRpt_m4l_xFF_dn","m_{4l};m_{4l};",500,0.,2500.);

  histo1d_["3M_ABCD_MET_dphi"] = fs->make<TH1D>("3M_ABCD_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["3M_ABCD_MET_ratioPt"] = fs->make<TH1D>("3M_ABCD_MET_ratioPt","p_{T} ratio;R(p_{T});",200,0.,5.);

  histo1d_["3M_check_resolved_M1M2_dR"] = fs->make<TH1D>("3M_check_resolved_M1M2_dR","#Delta R(M1,M2);#Delta R;",128,0.,6.4);
  histo1d_["3M_check_resolved_M1MM_dR"] = fs->make<TH1D>("3M_check_resolved_M1MM_dR","#Delta R(M1,MM);#Delta R;",128,0.,6.4);
  histo1d_["3M_check_resolved_M2MM_dR"] = fs->make<TH1D>("3M_check_resolved_M2MM_dR","#Delta R(M2,MM);#Delta R;",128,0.,6.4);
  histo1d_["3M_check_resolved_M1M2_invM"] = fs->make<TH1D>("3M_check_resolved_M1M2_invM","M(M1,M2);M(M1,M2);",400,0.,400.);
  histo1d_["3M_check_resolved_M1M2_invM_zoomed"] = fs->make<TH1D>("3M_check_resolved_M1M2_invM_zoomed","M(M1,M2);M(M1,M2);",400,0.,20.);
  histo1d_["3M_check_resolved_M1M2MM_invM"] = fs->make<TH1D>("3M_check_resolved_M1M2MM_invM","M;GeV;",400,0.,400);

  histo1d_["3M_check_M1M2_dR"] = fs->make<TH1D>("3M_check_M1M2_dR","#Delta R(M1,M2);#Delta R;",100,0.,0.5);
  histo1d_["3M_check_M1MM_dR"] = fs->make<TH1D>("3M_check_M1MM_dR","#Delta R(M1,MM);#Delta R;",128,0.,6.4);
  histo1d_["3M_check_M2MM_dR"] = fs->make<TH1D>("3M_check_M2MM_dR","#Delta R(M2,MM);#Delta R;",128,0.,6.4);
  histo1d_["3M_check_M1M2_invM"] = fs->make<TH1D>("3M_check_M1M2_invM","M(M1,M2);M(M1,M2);",400,0.,200.);
  histo1d_["3M_check_M1M2_invM_zoomed"] = fs->make<TH1D>("3M_check_M1M2_invM_zoomed","M(M1,M2);M(M1,M2);",400,0.,20.);
  histo1d_["3M_check_M1M2MM_invM"] = fs->make<TH1D>("3M_check_M1M2MM_invM","M;GeV;",400,0.,400);
  histo1d_["3M_check_M1_iso"] = fs->make<TH1D>("3M_check_M1_iso","M1 Track iso;#Sigma p_{T iso};",400,0.,100.);
  histo1d_["3M_check_M2_iso"] = fs->make<TH1D>("3M_check_M2_iso","M2 Track iso;#Sigma p_{T iso};",400,0.,100.);
  histo1d_["3M_check_MM_iso"] = fs->make<TH1D>("3M_check_MM_iso","MM Track iso;#Sigma p_{T iso};",400,0.,100.);
  histo1d_["3M_check_M1_iso_musubtract"] = fs->make<TH1D>("3M_check_M1_iso_musubtract","M1 Track iso - M2 p_{T};#Sigma p_{T iso} - p_{T #mu};",400,0.,100.);
  histo1d_["3M_check_M2_iso_musubtract"] = fs->make<TH1D>("3M_check_M2_iso_musubtract","M2 Track iso - M1 p_{T};#Sigma p_{T iso} - p_{T #mu};",400,0.,100.);

  histo2d_["3M_mt_dphi"] = fs->make<TH2D>("3M_mt_dphi","m_{T} vs #Delta#phi;m_{T};#Delta#phi",100,0.,500.,128,-3.2,3.2);
  histo2d_["3M_mt_ratioPt"] = fs->make<TH2D>("3M_mt_ratioPt","m_{T} vs p_{T} ratio;m_{T};p_{T} ratio",100,0.,500.,100,0.,5.);
  histo2d_["3M_dphi_ratioPt"] = fs->make<TH2D>("3M_dphi_ratioPt","#Delta#phi vs p_{T} ratio;#Delta#phi;p_{T} ratio",128,-3.2,3.2,100,0.,5.);

  interestEvtTree_ = fs->make<TTree>("evtTree","evtTree");
  interestEvtTree_->Branch("runNo",&runNo_,"runNo/i");
  interestEvtTree_->Branch("lumiNo",&lumiNo_,"lumiNo/i");
  interestEvtTree_->Branch("evtNo",&evtNo_,"evtNo/l");

  nonIsoSRTree_ = fs->make<TTree>("nonIsoSRTree","nonIsoSRTree");
  nonIsoSRTree_->Branch("wgt",&wgt_,"wgt/F");
  nonIsoSRTree_->Branch("m1pt",&m1pt_,"m1pt/F");
  nonIsoSRTree_->Branch("m1eta",&m1eta_,"m1eta/F");
  nonIsoSRTree_->Branch("m1phi",&m1phi_,"m1phi/F");
  nonIsoSRTree_->Branch("m1iso",&m1iso_,"m1iso/F");
  nonIsoSRTree_->Branch("m1global",&m1global_,"m1global/I");
  nonIsoSRTree_->Branch("m2pt",&m2pt_,"m2pt/F");
  nonIsoSRTree_->Branch("m2eta",&m2eta_,"m2eta/F");
  nonIsoSRTree_->Branch("m2phi",&m2phi_,"m2phi/F");
  nonIsoSRTree_->Branch("m2iso",&m2iso_,"m2iso/F");
  nonIsoSRTree_->Branch("m2global",&m2global_,"m2global/I");
  nonIsoSRTree_->Branch("mmpt",&mmpt_,"mmpt/F");
  nonIsoSRTree_->Branch("mmeta",&mmeta_,"mmeta/F");
  nonIsoSRTree_->Branch("mmphi",&mmphi_,"mmphi/F");
  nonIsoSRTree_->Branch("mmiso",&mmiso_,"mmiso/F");
  nonIsoSRTree_->Branch("metPt",&metPt_,"metPt/F");
  nonIsoSRTree_->Branch("metPhi",&metPhi_,"metPhi/F");
  nonIsoSRTree_->Branch("mt",&mt_,"mt/F");
  nonIsoSRTree_->Branch("mm1m2",&mm1m2_,"mm1m2/F");
  nonIsoSRTree_->Branch("ratioPt",&ratioPt_,"ratioPt/F");
  nonIsoSRTree_->Branch("dphiMet",&dphiMet_,"dphiMet/F");
}

void MergedMuNonIsoCRanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<edm::View<pat::MET>> metHandle;
  iEvent.getByToken(metToken_, metHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  std::vector<reco::GenParticleRef> promptMuons;

  double aWeight = 1.;
  double purwgtNo = 1.;
  double purwgtUp = 1.;
  double purwgtDn = 1.;
  double prefireNo = 1.;
  double prefireUp = 1.;
  double prefireDn = 1.;

  if (isMC_) {
    edm::Handle<double> theprefweight;
    iEvent.getByToken(prefweight_token, theprefweight);
    prefireNo = *theprefweight;

    edm::Handle<double> theprefweightUp;
    iEvent.getByToken(prefweightUp_token_, theprefweightUp);
    prefireUp = *theprefweightUp;

    edm::Handle<double> theprefweightDn;
    iEvent.getByToken(prefweightDn_token_, theprefweightDn);
    prefireDn = *theprefweightDn;

    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(generatorToken_, genInfo);
    double mcweight = genInfo->weight();

    aWeight = prefireNo*mcweight/std::abs(mcweight);

    edm::Handle<edm::View<PileupSummaryInfo>> pusummary;
    iEvent.getByToken(pileupToken_, pusummary);

    for (unsigned int idx = 0; idx < pusummary->size(); ++idx) {
      const auto& apu = pusummary->refAt(idx);

      int bx = apu->getBunchCrossing();

      if (bx==0) { // in-time PU only
        purwgtUp = purwgt_->at(puname_)->evaluate({apu->getTrueNumInteractions(),"up"});
        purwgtDn = purwgt_->at(puname_)->evaluate({apu->getTrueNumInteractions(),"down"});
        purwgtNo = purwgt_->at(puname_)->evaluate({apu->getTrueNumInteractions(),"nominal"});

        aWeight *= purwgtNo;

        break;
      }
    }

    edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
    iEvent.getByToken(genptcToken_, genptcHandle);

    for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
      const auto& genptc = genptcHandle->refAt(idx);

      if ( ( std::abs(genptc->pdgId())==13 ) && genptc->fromHardProcessFinalState() )
        promptMuons.push_back(genptc.castTo<reco::GenParticleRef>());
    }

    auto sortByPt = [](const reco::GenParticleRef& a, const reco::GenParticleRef& b) {
      return a->pt() > b->pt();
    };

    if (promptMuons.size()==4) {
      histo1d_["totWeightedSum_4M"]->Fill(0.5,aWeight);

      std::sort(promptMuons.begin(),promptMuons.end(),sortByPt);

      histo1d_["GEN_pt_1st"]->Fill(promptMuons.at(0)->pt(),aWeight);
      histo1d_["GEN_pt_2nd"]->Fill(promptMuons.at(1)->pt(),aWeight);
      histo1d_["GEN_pt_3rd"]->Fill(promptMuons.at(2)->pt(),aWeight);
      histo1d_["GEN_pt_4th"]->Fill(promptMuons.at(3)->pt(),aWeight);
    }
  }

  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);
  histo1d_["cutflow"]->Fill( 0.5, aWeight );

  edm::Handle<edm::TriggerResults> trigResultHandle;
  iEvent.getByToken(triggerToken_,trigResultHandle);

  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjHandle;
  iEvent.getByToken(triggerobjectsToken_, trigObjHandle);

  const unsigned int nTrig = trigResultHandle.product()->size();
  const edm::TriggerNames trigList = iEvent.triggerNames(*trigResultHandle);

  bool isFired = false;

  for (unsigned int iTrig = 0; iTrig != nTrig; iTrig++) {
    const std::string& trigName = trigList.triggerName(iTrig);
    for (unsigned int jTrig = 0; jTrig != trigList_.size(); jTrig++) {
      if (trigName.find(trigList_.at(jTrig).substr(0, trigList_.at(jTrig).find("*"))) != std::string::npos) {
        if (trigResultHandle.product()->accept(iTrig))
          isFired = true;
      }
    } // wanted triggers
  } // fired triggers

  std::vector<edm::RefToBase<pat::TriggerObjectStandAlone>> trigObjs;

  for (unsigned iTrig=0; iTrig < trigObjHandle->size(); iTrig++) {
    const auto& trigObj = trigObjHandle->refAt(iTrig);
    auto trigObjInst = trigObjHandle->at(iTrig); // workaround for copy
    trigObjInst.unpackPathNames(trigList);
    const auto& pathNames = trigObjInst.pathNames();

    for (const auto name : pathNames) {
      for (unsigned int jTrig = 0; jTrig < trigList_.size(); jTrig++) {
        if ( name.find(trigList_.at(jTrig).substr(0, trigList_.at(jTrig).find("*"))) != std::string::npos &&
             trigObjInst.hasPathName(name,true,true) ) {
          trigObjs.push_back(trigObj);
        }
      } // wanted triggers
    } // fired triggers
  } // trigger objs

  if (!isFired)
    return;

  histo1d_["cutflow"]->Fill( 1.5, aWeight ); // fired trigger

  edm::Handle<edm::TriggerResults> METfilterHandle;
  iEvent.getByToken(METfilterToken_,METfilterHandle);
  edm::TriggerNames METfilters = iEvent.triggerNames(*METfilterHandle);

  unsigned int nPassedFilters = 0;

  for (unsigned int iTrig = 0; iTrig < METfilterHandle.product()->size(); iTrig++) {
    const std::string trigname = METfilters.triggerName(iTrig);

    if (METfilterHandle.product()->accept(iTrig)) {
      for (const auto& filterName : METfilterList_) {
        if (trigname.find(filterName) != std::string::npos)
          nPassedFilters++;
      }
    }
  }

  if (nPassedFilters!=METfilterList_.size())
    return;

  histo1d_["cutflow"]->Fill( 2.5, aWeight ); // MET filter

  reco::Vertex primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->at(0);
  else
    return;

  std::vector<pat::MuonRef> highPtMuons;
  std::vector<pat::MuonRef> highPtTrackerMuons; // but not highPtMuon

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    const auto& aMuon = muonHandle->refAt(idx);

    if ( aMuon->tunePMuonBestTrack()->pt() < ptMuThres_||
         std::abs(aMuon->tunePMuonBestTrack()->eta()) > 2.4 )
      continue;

    const auto casted = aMuon.castTo<pat::MuonRef>();

    if ( muon::isHighPtMuon(*aMuon,primaryVertex) )
      highPtMuons.push_back(casted);
    else if ( muon::isTrackerHighPtMuon(*aMuon,primaryVertex) )
      highPtTrackerMuons.push_back(casted);
    else {}
  }

  std::map<pat::MuonRef,float> highPtMuonMap;
  std::map<pat::MuonRef,float> highPtTrackerMuonMap;

  mucorrHelper_.modifiedIsoWithVal(highPtMuonMap,
                                   highPtTrackerMuonMap,
                                   highPtMuons,
                                   highPtTrackerMuons,
                                   *beamSpotHandle);

  if ( highPtMuonMap.size() < 2 )
    return;

  histo1d_["cutflow"]->Fill( 3.5, aWeight ); // at least two high-pt muons

  unsigned int nHighPtMuons = highPtMuonMap.size() + highPtTrackerMuonMap.size();

  pat::MuonRef leadingMuon;

  for (const auto& mu : highPtMuonMap) {
    if (!leadingMuon.isNonnull() ||
        mu.first->tunePMuonBestTrack()->pt() > leadingMuon->tunePMuonBestTrack()->pt()) {
      leadingMuon = mu.first;
    }
  }

  if ( leadingMuon->tunePMuonBestTrack()->pt() < 52. )
    return;

  histo1d_["cutflow"]->Fill( 4.5, aWeight ); // the leading muon passes trigger pt threshold

  bool trigMatched = false;
  std::vector<double> trigSyst;

  for (unsigned idx = 0; idx < trigObjs.size(); idx++) {
    const auto& trigObj = trigObjs.at(idx);

    const auto& leadMu = leadingMuon;
    const bool isGlobal = true; // true by construction
    // const bool isGlobal = std::find(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),leadMu)!=isolatedHighPtMuons.end();

    if ( reco::deltaR2(trigObj->eta(),
                       trigObj->phi(),
                       leadMu->tunePMuonBestTrack()->eta(),
                       leadMu->tunePMuonBestTrack()->phi()) < 0.01 ) {
      trigMatched = true;

      if (isMC_) {
        if (isGlobal) {
          aWeight *= mucorrHelper_.trigSFglobal(leadMu);
          trigSyst.push_back( mucorrHelper_.trigSFglobalSyst(leadMu) /
                              mucorrHelper_.trigSFglobal(leadMu) );
        } else {
          aWeight *= mucorrHelper_.trigSFtracker(leadMu);
          trigSyst.push_back( mucorrHelper_.trigSFtrackerSyst(leadMu) /
                              mucorrHelper_.trigSFtracker(leadMu) );
        }
      }

      break;
    }
  }

  if ( !trigMatched )
    return;

  histo1d_["cutflow"]->Fill( 5.5, aWeight ); // trigger matching

  std::vector<pat::MuonRef> allHighPtMuons(highPtMuons);
  allHighPtMuons.insert( allHighPtMuons.end(), highPtTrackerMuons.begin(), highPtTrackerMuons.end() );

  std::vector<pat::MuonRef> nonHighPtMuons;
  std::vector<pat::MuonRef> nonHighPtMuonsVLiso;
  std::map<pat::MuonRef,float> nonHighPtIsos;

  mucorrHelper_.nonHighPtMuonIso(nonHighPtMuons,
                                 nonHighPtMuonsVLiso,
                                 nonHighPtIsos,
                                 muonHandle,
                                 allHighPtMuons,
                                 highPtMuons,
                                 highPtTrackerMuons,
                                 primaryVertex,
                                 *beamSpotHandle,
                                 ptMuThres_);

  if (!nonHighPtMuonsVLiso.empty())
    return;

  histo1d_["cutflow"]->Fill( 6.5, aWeight ); // reject non-high-pt muons (to avoid overlap with CR)

  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  bool hasEle = false;

  for (unsigned int idx = 0; idx < eleHandle->size(); idx++) {
    const auto& aEle = eleHandle->refAt(idx);

    if ( std::abs(aEle->superCluster()->eta()) > 2.5 )
      continue;

    // veto EBEE gap
    if ( std::abs(aEle->superCluster()->eta()) > 1.4442 && std::abs(aEle->superCluster()->eta()) < 1.566 )
      continue;

    if ( aEle->electronID("modifiedHeepElectronID") )
      hasEle = true;
  }

  if (hasEle)
    return;

  histo1d_["cutflow"]->Fill( 7.5, aWeight ); // no HEEP electrons

  std::vector<double> idSyst;
  std::vector<double> recoSyst;

  auto systSFratio = [] (const std::vector<double>& vec) -> std::pair<double,double> {
    double up = 1., dn = 1.;

    for (const auto& systOverSF : vec) {
      up *= (1.+systOverSF);
      dn *= std::max(1.-systOverSF,0.);
    }

    return std::make_pair(up,dn); // ratio to the nominal
  };

  if (isMC_) {
    for (const auto& aMu : highPtMuons) {
      aWeight *= mucorrHelper_.highptIdSF(aMu);
      aWeight *= mucorrHelper_.recoSF(aMu);
      idSyst.push_back( mucorrHelper_.highptIdSFsyst(aMu)/mucorrHelper_.highptIdSF(aMu) );
      recoSyst.push_back( mucorrHelper_.recoSFsyst(aMu)/mucorrHelper_.recoSF(aMu) );
    }

    for (const auto& aMu : highPtTrackerMuons) {
      aWeight *= mucorrHelper_.trkHighptIdSF(aMu);
      aWeight *= mucorrHelper_.recoSF(aMu);
      idSyst.push_back( mucorrHelper_.trkHighptIdSFsyst(aMu)/mucorrHelper_.trkHighptIdSF(aMu) );
    }
  }

  if ( nHighPtMuons==3 ) {
    histo1d_["cutflow"]->Fill( 8.5, aWeight ); // 3 high-pt muons

    pat::MuonRef aMu = allHighPtMuons.at(0);
    pat::MuonRef bMu = allHighPtMuons.at(1);
    pat::MuonRef cMu = allHighPtMuons.at(2);
    std::vector<std::pair<pat::MuonRef,pat::MuonRef>> muPairs;
    muPairs.push_back(std::make_pair(aMu,bMu));
    muPairs.push_back(std::make_pair(bMu,cMu));
    muPairs.push_back(std::make_pair(cMu,aMu));
    const auto lvecMa = math::PtEtaPhiMLorentzVector(aMu->tunePMuonBestTrack()->pt(),
                                                     aMu->tunePMuonBestTrack()->eta(),
                                                     aMu->tunePMuonBestTrack()->phi(),
                                                     mumass_);
    const auto lvecMb = math::PtEtaPhiMLorentzVector(bMu->tunePMuonBestTrack()->pt(),
                                                     bMu->tunePMuonBestTrack()->eta(),
                                                     bMu->tunePMuonBestTrack()->phi(),
                                                     mumass_);
    const auto lvecMc = math::PtEtaPhiMLorentzVector(cMu->tunePMuonBestTrack()->pt(),
                                                     cMu->tunePMuonBestTrack()->eta(),
                                                     cMu->tunePMuonBestTrack()->phi(),
                                                     mumass_);
    const double mass3M = (lvecMa+lvecMb+lvecMc).M();

    if ( mass3M < 50. )
      return;

    histo1d_["cutflow"]->Fill( 9.5, aWeight ); // inv mass cut

    auto sortByDR = [] (const std::pair<pat::MuonRef,pat::MuonRef>& apair, const std::pair<pat::MuonRef,pat::MuonRef>& bpair) {
      const auto& a1 = apair.first->tunePMuonBestTrack();
      const auto& a2 = apair.second->tunePMuonBestTrack();
      const auto& b1 = bpair.first->tunePMuonBestTrack();
      const auto& b2 = bpair.second->tunePMuonBestTrack();
      return reco::deltaR2(a1->eta(),a1->phi(),a2->eta(),a2->phi()) < reco::deltaR2(b1->eta(),b1->phi(),b2->eta(),b2->phi());
    };

    std::sort(muPairs.begin(),muPairs.end(),sortByDR);
    const auto& M1M2 = muPairs.at(0);
    const auto& M2MM = muPairs.at(1);
    const auto& M1MM = muPairs.at(2);
    const auto lvecM1 = math::PtEtaPhiMLorentzVector(M1M2.first->tunePMuonBestTrack()->pt(),
                                                     M1M2.first->tunePMuonBestTrack()->eta(),
                                                     M1M2.first->tunePMuonBestTrack()->phi(),
                                                     mumass_);
    const auto lvecM2 = math::PtEtaPhiMLorentzVector(M1M2.second->tunePMuonBestTrack()->pt(),
                                                     M1M2.second->tunePMuonBestTrack()->eta(),
                                                     M1M2.second->tunePMuonBestTrack()->phi(),
                                                     mumass_);
    const auto lvecM1M2 = lvecM1 + lvecM2;

    histo1d_["3M_check_resolved_M1M2_dR"]->Fill( reco::deltaR(M1M2.first->tunePMuonBestTrack()->eta(),
                                                              M1M2.first->tunePMuonBestTrack()->phi(),
                                                              M1M2.second->tunePMuonBestTrack()->eta(),
                                                              M1M2.second->tunePMuonBestTrack()->phi()) , aWeight );
    histo1d_["3M_check_resolved_M2MM_dR"]->Fill( reco::deltaR(M2MM.first->tunePMuonBestTrack()->eta(),
                                                              M2MM.first->tunePMuonBestTrack()->phi(),
                                                              M2MM.second->tunePMuonBestTrack()->eta(),
                                                              M2MM.second->tunePMuonBestTrack()->phi()) , aWeight );
    histo1d_["3M_check_resolved_M1MM_dR"]->Fill( reco::deltaR(M1MM.first->tunePMuonBestTrack()->eta(),
                                                              M1MM.first->tunePMuonBestTrack()->phi(),
                                                              M1MM.second->tunePMuonBestTrack()->eta(),
                                                              M1MM.second->tunePMuonBestTrack()->phi()) , aWeight );
    histo1d_["3M_check_resolved_M1M2_invM"]->Fill( lvecM1M2.M() , aWeight );
    histo1d_["3M_check_resolved_M1M2_invM_zoomed"]->Fill( lvecM1M2.M() , aWeight );
    histo1d_["3M_check_resolved_M1M2MM_invM"]->Fill( mass3M , aWeight );

    double drThres2 = drThres_*drThres_;
    std::vector<std::pair<pat::MuonRef,pat::MuonRef>> cands;

    // find a collimated pair of muons
    for (unsigned idx = 0; idx < highPtMuons.size(); idx++) { // one of them should be global
      const auto& cand1 = highPtMuons.at(idx);

      for (unsigned jdx = idx+1; jdx < allHighPtMuons.size(); jdx++) {
        const auto& cand2 = allHighPtMuons.at(jdx);

        // they shouldn't be the same, throw exception otherwise
        if ( cand1==cand2 )
          throw cms::Exception("LogicError") << "Error: MergedMuNonIsoCRanalyzer::analyze - attempting to compare the same muons " << std::endl;

        // check CWR
        double dr2 = reco::deltaR2(cand1->tunePMuonBestTrack()->eta(),cand1->tunePMuonBestTrack()->phi(),
                                   cand2->tunePMuonBestTrack()->eta(),cand2->tunePMuonBestTrack()->phi());

        // const bool isGlobalHighPt = std::find(highPtMuons.begin(),highPtMuons.end(),cand2)!=highPtMuons.end();

        if ( dr2 < drThres2 ) // || isGlobalHighPt )
          cands.push_back(std::make_pair(cand1,cand2));
      }
    }

    if ( cands.empty() )
      return;

    histo1d_["cutflow"]->Fill( 10.5, aWeight ); // has two collimated muons

    std::sort(cands.begin(),cands.end(),sortByDR);

    // third muon is the merged muon - should be highPtMuon
    const pat::MuonRef& firstMuon = cands.front().first;
    const pat::MuonRef& secondMuon = cands.front().second;
    pat::MuonRef mergedMuon;

    for (unsigned idx = 0; idx < highPtMuons.size(); idx++) {
      const auto& cand = highPtMuons.at(idx);

      if ( cand!=firstMuon && cand!=secondMuon ) {
        mergedMuon = cand;
        break;
      }
    }

    if ( !firstMuon.isNull() && !secondMuon.isNull() && !mergedMuon.isNull() ) {
      histo1d_["cutflow"]->Fill( 11.5, aWeight ); // and the other muons is a high-pt muon

      double momcorr1 = 1.;
      double momcorr2 = 1.;
      double momcorrMM = 1.;

      const auto lvecFirst = math::PtEtaPhiMLorentzVector(firstMuon->tunePMuonBestTrack()->pt(),
                                                          firstMuon->tunePMuonBestTrack()->eta(),
                                                          firstMuon->tunePMuonBestTrack()->phi(),
                                                          mumass_);
      const auto lvecSecond = math::PtEtaPhiMLorentzVector(secondMuon->tunePMuonBestTrack()->pt(),
                                                           secondMuon->tunePMuonBestTrack()->eta(),
                                                           secondMuon->tunePMuonBestTrack()->phi(),
                                                           mumass_);
      const auto lvecMerged = math::PtEtaPhiMLorentzVector(mergedMuon->tunePMuonBestTrack()->pt(),
                                                           mergedMuon->tunePMuonBestTrack()->eta(),
                                                           mergedMuon->tunePMuonBestTrack()->phi(),
                                                           mumass_);

      auto firstP4 = lvecFirst*momcorr1;
      auto secondP4 = lvecSecond*momcorr2;
      auto mergedP4 = lvecMerged*momcorrMM;

      histo1d_["3M_check_M1M2_invM"]->Fill( (firstP4+secondP4).M() , aWeight );
      histo1d_["3M_check_M1M2_invM_zoomed"]->Fill( (firstP4+secondP4).M() , aWeight );
      histo1d_["3M_check_M1M2MM_invM"]->Fill( (firstP4+secondP4+mergedP4).M() , aWeight );

      histo1d_["3M_check_M1M2_dR"]->Fill( reco::deltaR(firstP4.eta(),firstP4.phi(),secondP4.eta(),secondP4.phi()) , aWeight );
      histo1d_["3M_check_M1MM_dR"]->Fill( reco::deltaR(firstP4.eta(),firstP4.phi(),mergedP4.eta(),mergedP4.phi()) , aWeight );
      histo1d_["3M_check_M2MM_dR"]->Fill( reco::deltaR(secondP4.eta(),secondP4.phi(),mergedP4.eta(),mergedP4.phi()) , aWeight );
      histo1d_["3M_check_M1_iso"]->Fill( firstMuon->trackIso() , aWeight );
      histo1d_["3M_check_M2_iso"]->Fill( secondMuon->trackIso() , aWeight );
      histo1d_["3M_check_MM_iso"]->Fill( mergedMuon->trackIso() , aWeight );
      bool insideIsoVeto = reco::deltaR2(firstMuon->innerTrack()->eta(),firstMuon->innerTrack()->phi(),
                                         secondMuon->innerTrack()->eta(),secondMuon->innerTrack()->phi()) < 0.0001;
      histo1d_["3M_check_M1_iso_musubtract"]->Fill( insideIsoVeto ? firstMuon->trackIso() :
                                                                    firstMuon->trackIso() - secondMuon->innerTrack()->pt() , aWeight );
      histo1d_["3M_check_M2_iso_musubtract"]->Fill( insideIsoVeto ? secondMuon->trackIso() :
                                                                    secondMuon->trackIso() - firstMuon->innerTrack()->pt() , aWeight );

      const auto& aMET = metHandle->at(0);
      const auto metp4 = aMET.p4();
      const auto tpMET = metp4 - firstP4 - secondP4 - mergedP4
                               + firstMuon->p4() + secondMuon->p4() + mergedMuon->p4();

      if ( tpMET.pt() > ptThres_ ) {
        histo1d_["cutflow"]->Fill( 12.5, aWeight ); // MET cut

        double dphi = reco::deltaPhi(mergedP4.phi(),tpMET.phi());
        auto p4 = firstP4 + secondP4 + mergedP4 + tpMET;
        double mt = p4.mt();

        // substitute eta hypothesis to check Greg's comment
        const auto metp4_etaHypo = math::PtEtaPhiMLorentzVector(aMET.pt(),mergedMuon->tunePMuonBestTrack()->eta(),aMET.phi(),0.);

        // eta hypothesis variations
        const auto tpMET_etaHypo = metp4_etaHypo - firstP4 - secondP4 - mergedP4
                                                 + firstMuon->p4() + secondMuon->p4() + mergedMuon->p4();
        const auto lvec4lMET_etaHypo = firstP4 + secondP4 + mergedP4 + tpMET_etaHypo;
        const double m4lMET = lvec4lMET_etaHypo.M();

        bool passDPhi = std::abs(dphi) < drThres_;
        bool passDPhiCR = (std::abs(dphi) < drThresCR_) && !passDPhi;
        bool passRatioPt = false;

        double ratioPt = (firstP4+secondP4+mergedP4).pt()/tpMET.pt();

        if ( ratioThresLo_ < ratioPt && ratioPt < ratioThresHi_ )
          passRatioPt = true;

        const double ffMM = std::max(ffFunc_->Eval(mt),0.);
        const double xvalMET[1] = {mt};
        double ciMM[1];
        fitResult_->GetConfidenceIntervals(1,1,0,xvalMET,ciMM,0.95,false);

        histo1d_["nPV"]->Fill( static_cast<float>(pvHandle->size())+0.5, aWeight );
        histo1d_["3M_ABCD_MET_dphi"]->Fill(dphi, aWeight);
        histo1d_["3M_ABCD_MET_ratioPt"]->Fill(ratioPt, aWeight);

        histo2d_["3M_mt_dphi"]->Fill(mt,dphi,aWeight);
        histo2d_["3M_mt_ratioPt"]->Fill(mt,ratioPt,aWeight);
        histo2d_["3M_dphi_ratioPt"]->Fill(dphi,ratioPt,aWeight);

        auto allMuonMaps = highPtMuonMap;
        allMuonMaps.insert(highPtTrackerMuonMap.begin(), highPtTrackerMuonMap.end());

        wgt_ = aWeight;
        m1pt_ = firstP4.pt();
        m1eta_ = firstP4.eta();
        m1phi_ = firstP4.phi();
        m1iso_ = allMuonMaps.at(firstMuon);
        m1global_ = firstMuon->isGlobalMuon();
        m2pt_ = secondP4.pt();
        m2eta_ = secondP4.eta();
        m2phi_ = secondP4.phi();
        m2iso_ = allMuonMaps.at(secondMuon);
        m2global_ = secondMuon->isGlobalMuon();
        mmpt_ = mergedP4.pt();
        mmeta_ = mergedP4.eta();
        mmphi_ = mergedP4.phi();
        mmiso_ = allMuonMaps.at(mergedMuon);
        metPt_ = tpMET.pt();
        metPhi_ = tpMET.phi();
        mt_ = mt;
        mm1m2_ = (firstP4 + secondP4).M();
        ratioPt_ = ratioPt;
        dphiMet_ = dphi;
        nonIsoSRTree_->Fill();

        if (passDPhi)
          histo1d_["cutflow"]->Fill( 13.5, aWeight ); // delta phi cut

        if ( passDPhi && passRatioPt ) {
          histo1d_["cutflow"]->Fill( 14.5, aWeight ); // SR

          if ( true /*mt < 500. || isMC_ (unblinded)*/ ) { // blinded
            histo1d_["3M_MM_pt"]->Fill( mergedP4.pt(), aWeight );
            histo1d_["3M_MM_eta"]->Fill( mergedP4.eta(), aWeight );
            histo1d_["3M_MM_phi"]->Fill( mergedP4.phi(), aWeight );

            histo1d_["3M_MET_pt"]->Fill( tpMET.pt(), aWeight );
            histo1d_["3M_MET_phi"]->Fill( tpMET.phi(), aWeight );

            histo1d_["3M_M1_pt"]->Fill( firstP4.pt(), aWeight );
            histo1d_["3M_M1_eta"]->Fill( firstP4.eta(), aWeight );
            histo1d_["3M_M1_phi"]->Fill( firstP4.phi(), aWeight );

            histo1d_["3M_M2_pt"]->Fill( secondP4.pt(), aWeight );
            histo1d_["3M_M2_eta"]->Fill( secondP4.eta(), aWeight );
            histo1d_["3M_M2_phi"]->Fill( secondP4.phi(), aWeight );

            histo1d_["3M_MET_dphi"]->Fill( dphi, aWeight );
            histo1d_["3M_mt"]->Fill( mt, aWeight );
            histo1d_["3M_mt_idUp"]->Fill( mt, aWeight*systSFratio(idSyst).first );
            histo1d_["3M_mt_idDn"]->Fill( mt, aWeight*systSFratio(idSyst).second );
            histo1d_["3M_mt_trigUp"]->Fill( mt, aWeight*systSFratio(trigSyst).first );
            histo1d_["3M_mt_trigDn"]->Fill( mt, aWeight*systSFratio(trigSyst).second );
            histo1d_["3M_mt_recoUp"]->Fill( mt, aWeight*systSFratio(recoSyst).first );
            histo1d_["3M_mt_recoDn"]->Fill( mt, aWeight*systSFratio(recoSyst).second );
            histo1d_["3M_mt_PUrwgtUp"]->Fill( mt, aWeight*purwgtUp/purwgtNo );
            histo1d_["3M_mt_PUrwgtDn"]->Fill( mt, aWeight*purwgtDn/purwgtNo );
            histo1d_["3M_mt_prefireUp"]->Fill( mt, aWeight*prefireUp/prefireNo );
            histo1d_["3M_mt_prefireDn"]->Fill( mt, aWeight*prefireDn/prefireNo );

            // eta hypothesis variations
            histo1d_["3M_m4l"]->Fill( m4lMET, aWeight );
            histo1d_["3M_m4l_idUp"]->Fill( m4lMET, aWeight*systSFratio(idSyst).first );
            histo1d_["3M_m4l_idDn"]->Fill( m4lMET, aWeight*systSFratio(idSyst).second );
            histo1d_["3M_m4l_trigUp"]->Fill( m4lMET, aWeight*systSFratio(trigSyst).first );
            histo1d_["3M_m4l_trigDn"]->Fill( m4lMET, aWeight*systSFratio(trigSyst).second );
            histo1d_["3M_m4l_recoUp"]->Fill( m4lMET, aWeight*systSFratio(recoSyst).first );
            histo1d_["3M_m4l_recoDn"]->Fill( m4lMET, aWeight*systSFratio(recoSyst).second );
            histo1d_["3M_m4l_PUrwgtUp"]->Fill( m4lMET, aWeight*purwgtUp/purwgtNo );
            histo1d_["3M_m4l_PUrwgtDn"]->Fill( m4lMET, aWeight*purwgtDn/purwgtNo );
            histo1d_["3M_m4l_prefireUp"]->Fill( m4lMET, aWeight*prefireUp/prefireNo );
            histo1d_["3M_m4l_prefireDn"]->Fill( m4lMET, aWeight*prefireDn/prefireNo );

            if (mt > 800.) {
              runNo_ = iEvent.id().run();
              lumiNo_ = iEvent.id().luminosityBlock();
              evtNo_ = iEvent.id().event();
              interestEvtTree_->Fill();

              std::cout << "Interesting event: "
                        << " run " << runNo_
                        << " lumi " << lumiNo_
                        << " event " << evtNo_
                        << std::endl;
              std::cout << "Event details: "
                        << " mt " << mt
                        << " dphi " << dphi
                        << " ratioPt " << ratioPt
                        << std::endl;
              std::cout << "Muons: " << std::endl;
              std::cout << "    First muon:" << std::endl;
              std::cout << "        pt " << firstP4.pt() << " eta " << firstP4.eta() << " phi " << firstP4.phi() << std::endl;
              std::cout << "        ptErr " << firstMuon->tunePMuonBestTrack()->ptError()
                        << "        cov(1,1) " << firstMuon->tunePMuonBestTrack()->covariance(1,1)
                        << "        cov(2,2) " << firstMuon->tunePMuonBestTrack()->covariance(2,2)
                        << "        cov(1,2) " << firstMuon->tunePMuonBestTrack()->covariance(1,2)
                        << std::endl;
              std::cout << "    Second muon:" << std::endl;
              std::cout << "        pt " << secondP4.pt() << " eta " << secondP4.eta() << " phi " << secondP4.phi() << std::endl;
              std::cout << "        ptErr " << secondMuon->tunePMuonBestTrack()->ptError()
                        << "        cov(1,1) " << secondMuon->tunePMuonBestTrack()->covariance(1,1)
                        << "        cov(2,2) " << secondMuon->tunePMuonBestTrack()->covariance(2,2)
                        << "        cov(1,2) " << secondMuon->tunePMuonBestTrack()->covariance(1,2)
                        << std::endl;
              std::cout << "    Merged muon:" << std::endl;
              std::cout << "        pt " << mergedP4.pt() << " eta " << mergedP4.eta() << " phi " << mergedP4.phi() << std::endl;
              std::cout << "        ptErr " << mergedMuon->tunePMuonBestTrack()->ptError()
                        << "        cov(1,1) " << mergedMuon->tunePMuonBestTrack()->covariance(1,1)
                        << "        cov(2,2) " << mergedMuon->tunePMuonBestTrack()->covariance(2,2)
                        << "        cov(1,2) " << mergedMuon->tunePMuonBestTrack()->covariance(1,2)
                        << std::endl;
              std::cout << "MET: " << std::endl;
              std::cout << "    pt " << tpMET.pt() << " phi " << tpMET.phi() << " sumEt " << aMET.sumEt() << std::endl;
              std::cout << "    cov(0,0) " << aMET.getSignificanceMatrix().At(0,0)
                        << "    cov(1,1) " << aMET.getSignificanceMatrix().At(1,1)
                        << "    cov(0,1) " << aMET.getSignificanceMatrix().At(0,1)
                        << std::endl;
            } // mt > 800.
          } // blinded
        } else if ( passDPhi && !passRatioPt ) {
          histo1d_["3M_antiRpt_MM_pt"]->Fill( mergedP4.pt(), aWeight );
          histo1d_["3M_antiRpt_MM_pt_xFF"]->Fill( mergedP4.pt(), aWeight*ffMM );
          histo1d_["3M_antiRpt_M1_pt"]->Fill( firstP4.pt(), aWeight );
          histo1d_["3M_antiRpt_M2_pt"]->Fill( secondP4.pt(), aWeight );
          histo1d_["3M_antiRpt_MET_pt"]->Fill( tpMET.pt(), aWeight );
          histo1d_["3M_antiRpt_MET_pt_xFF"]->Fill( tpMET.pt(), aWeight*ffMM );
          histo1d_["3M_antiRpt_MET_dphi"]->Fill( dphi, aWeight );
          histo1d_["3M_antiRpt_mt"]->Fill( mt, aWeight );
          histo1d_["3M_antiRpt_mt_xFF"]->Fill( mt, aWeight*ffMM );
          histo1d_["3M_antiRpt_mt_xFF_up"]->Fill( mt, aWeight*(ffMM+ciMM[0]) );
          histo1d_["3M_antiRpt_mt_xFF_dn"]->Fill( mt, aWeight*(ffMM-std::min(ciMM[0],ffMM)) );
          histo1d_["3M_antiRpt_MET_ratioPt"]->Fill( ratioPt, aWeight );

          // eta hypothesis variations
          histo1d_["3M_antiRpt_m4l"]->Fill( m4lMET, aWeight );
          histo1d_["3M_antiRpt_m4l_xFF"]->Fill( m4lMET, aWeight*ffMM );
          histo1d_["3M_antiRpt_m4l_xFF_up"]->Fill( m4lMET, aWeight*(ffMM+ciMM[0]) );
          histo1d_["3M_antiRpt_m4l_xFF_dn"]->Fill( m4lMET, aWeight*(ffMM-std::min(ciMM[0],ffMM)) );
        }
      } // MET ptThres
    } // find a pair of collimated muons
  } // nHighPtMuons==3

  return;
}

DEFINE_FWK_MODULE(MergedMuNonIsoCRanalyzer);
