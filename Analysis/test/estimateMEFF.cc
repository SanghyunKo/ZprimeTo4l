#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

static double retrieveLumi(const std::string& anlyzrEra) {
  if (anlyzrEra=="20UL16APV")
    return 19.5;
  else if (anlyzrEra=="20UL16" || anlyzrEra=="")
    return 16.8;
  else if (anlyzrEra=="20UL17")
    return 41.48;
  else if (anlyzrEra=="20UL18")
    return 59.83;

  return 0.;
}

void estimateMEFF(TString era) {
  setTDRStyle();
  gStyle->SetOptFit(0);

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"

  static constexpr double WZxsec = 5.213; // 0.65*62.78;
  static constexpr double ZZxsec = 13.81;
  static TString postfix = era;
  TString fname = era;

  if (era=="20UL16APV") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "19.5 fb^{-1}";
  } else if (era=="20UL16") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "16.8 fb^{-1}";
    postfix = "";
  } else if (era=="20UL17") {
    lumi_sqrtS = "2017 (13 TeV)";
    lumi_13TeV = "41.48 fb^{-1}";
  } else if (era=="20UL18") {
    lumi_sqrtS = "2018 (13 TeV)";
    lumi_13TeV = "59.83 fb^{-1}";
  } else if (era=="run2") {
    lumi_sqrtS = "";
    lumi_13TeV = "137.6 fb^{-1}";
    postfix = "";
    fname = "20UL16";
  } else {
    std::cout << "check era..." << std::endl;
  }

  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 11;

  if( iPos==0 )
    relPosX = 0.12;

  int W = 800;
  int H = 800;

  int H_ref = 800;
  int W_ref = 800;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  TFile* datafile = new TFile("EleAnalyzer_"+fname+"_data.root","READ");
  TFile* WZfile = new TFile("EleAnalyzer_"+fname+"_WZFXFX.root","READ");
  TFile* ZZfile = new TFile("EleAnalyzer_"+fname+"_ZZ.root","READ");

  TFile* datafile1 = new TFile("EleAnalyzer_20UL16APV_data.root","READ");
  TFile* WZfile1 = new TFile("EleAnalyzer_20UL16APV_WZFXFX.root","READ");
  TFile* ZZfile1 = new TFile("EleAnalyzer_20UL16APV_ZZ.root","READ");

  TFile* datafile2 = new TFile("EleAnalyzer_20UL17_data.root","READ");
  TFile* WZfile2 = new TFile("EleAnalyzer_20UL17_WZFXFX.root","READ");
  TFile* ZZfile2 = new TFile("EleAnalyzer_20UL17_ZZ.root","READ");

  TFile* datafile3 = new TFile("EleAnalyzer_20UL18_data.root","READ");
  TFile* WZfile3 = new TFile("EleAnalyzer_20UL18_WZFXFX.root","READ");
  TFile* ZZfile3 = new TFile("EleAnalyzer_20UL18_ZZ.root","READ");

  TFile* datafile4 = new TFile("MuAnalyzer_"+fname+"_data.root","READ");
  TFile* WZfile4 = new TFile("MuAnalyzer_"+fname+"_WZFXFX.root","READ");
  TFile* ZZfile4 = new TFile("MuAnalyzer_"+fname+"_ZZ.root","READ");

  TFile* datafile5 = new TFile("MuAnalyzer_20UL16APV_data.root","READ");
  TFile* WZfile5 = new TFile("MuAnalyzer_20UL16APV_WZFXFX.root","READ");
  TFile* ZZfile5 = new TFile("MuAnalyzer_20UL16APV_ZZ.root","READ");

  TFile* datafile6 = new TFile("MuAnalyzer_20UL17_data.root","READ");
  TFile* WZfile6 = new TFile("MuAnalyzer_20UL17_WZFXFX.root","READ");
  TFile* ZZfile6 = new TFile("MuAnalyzer_20UL17_ZZ.root","READ");

  TFile* datafile7 = new TFile("MuAnalyzer_20UL18_data.root","READ");
  TFile* WZfile7 = new TFile("MuAnalyzer_20UL18_WZFXFX.root","READ");
  TFile* ZZfile7 = new TFile("MuAnalyzer_20UL18_ZZ.root","READ");

  auto estimateCenter = [] (const std::vector<double>& vec) -> std::vector<double> {
    std::vector<double> out;

    for (unsigned idx = 1; idx < vec.size()-1; idx++)
      out.push_back( (vec.at(idx) + vec.at(idx+1) ) / 2. );

    out.push_back(vec.back());

    return std::move(out);
  };

  auto estimateWidth = [] (const std::vector<double>& vec) -> std::vector<double> {
    std::vector<double> out;

    for (unsigned idx = 1; idx < vec.size()-1; idx++)
      out.push_back( ( vec.at(idx+1) - vec.at(idx) ) / 2. );

    out.push_back(0.);

    return std::move(out);
  };

  auto rebinnedHisto = [&] (TFile* afile, const TString& aname, std::vector<double>& binning, const TString& aEra, double xsec=0.) -> TH1D* {
    TH1D* ahist = (TH1D*)afile->Get(aname)->Clone();
    const int nbin = binning.size()-1;
    TH1D* rebin = (TH1D*)ahist->Rebin(nbin,TString(ahist->GetName())+"_rebin",&(binning[0]));

    if (xsec > 0.) {
      const double sumwgt = ((TH1D*)afile->Get("evtCounter/h_sumW"))->GetBinContent(1);
      rebin->Scale(xsec*retrieveLumi(aEra.Data())*1000./sumwgt);
    }

    return rebin;
  };

  auto subtractHist = [] (const TH1D* ahist, const TH1D* bhist) -> TH1D* {
    TH1D* result = (TH1D*)ahist->Clone();

    for (unsigned ibin = 0; ibin < ahist->GetNbinsX()+2; ibin++) {
      if (ahist->GetBinContent(ibin)==0.)
        continue;

      result->SetBinContent(ibin, std::max( ahist->GetBinContent(ibin) - bhist->GetBinContent(ibin), 0.) );
      result->SetBinError(ibin, result->GetBinContent(ibin)==0. ? 0. : std::hypot(ahist->GetBinError(ibin), bhist->GetBinError(ibin)) );
    }

    return result;
  };

  auto subtractPrompt = [&] (const TString& aname, std::vector<double>& binning, TFile* adata, TFile* aWZ, TFile* aZZ, const TString& analyzerName, const TString& aEra) -> TH1D* {
    TH1D* dataHist = rebinnedHisto(adata,analyzerName+"Data/"+aname,binning,aEra);

    TH1D* WZhist = rebinnedHisto(aWZ,analyzerName+aEra+"/"+aname,binning,aEra,WZxsec);
    TH1D* ZZhist = rebinnedHisto(aZZ,analyzerName+aEra+"/"+aname,binning,aEra,ZZxsec);

    return subtractHist( subtractHist( dataHist, WZhist ), ZZhist );
  };

/*
  TH1D* SSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_SSCR_EB_mixedME")->Clone();
  TH1D* SSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_SSCR_EB_antiME")->Clone();
  std::vector<double> xbinsSS = {0, 50, 60, 70, 100, 150, 250, 500, 1000};
  const int nbinsSS = xbinsSS.size()-1;
  std::vector<double> xcenSS = estimateCenter(xbinsSS);
  TH1D* SSnum_rebin = (TH1D*)SSnum->Rebin(nbinsSS, "2E_Et_SSCR_EB_mixedME", &(xbinsSS[0]));
  TH1D* SSdenom_rebin = (TH1D*)SSdenom->Rebin(nbinsSS, "2E_Et_SSCR_EB_antiME", &(xbinsSS[0]));

  SSnum_rebin->Divide( SSdenom_rebin );
*/

  std::vector<double> xbinsSS = {0,50,55,60,75,100,150,300,500}; //{0,50,60,70,80,90,100,120,150,200,500};
  std::vector<double> xcenSS = estimateCenter(xbinsSS);
  const int nbinsSS = xbinsSS.size()-1;

  TH1D* subtracted_EB_SSnumer = subtractPrompt("3E_Et_Ztag_EB_CRME",xbinsSS,datafile,WZfile,ZZfile,"mergedEleCRanalyzer",postfix);
  TH1D* subtracted_EB_SSdenom = subtractPrompt("3E_Et_Ztag_EB_antiME",xbinsSS,datafile,WZfile,ZZfile,"mergedEleCRanalyzer",postfix);

  TH1D* subtracted_EB_SSnumer1 = subtractPrompt("3E_Et_Ztag_EB_CRME",xbinsSS,datafile1,WZfile1,ZZfile1,"mergedEleCRanalyzer","20UL16APV");
  TH1D* subtracted_EB_SSdenom1 = subtractPrompt("3E_Et_Ztag_EB_antiME",xbinsSS,datafile1,WZfile1,ZZfile1,"mergedEleCRanalyzer","20UL16APV");

  TH1D* subtracted_EB_SSnumer2 = subtractPrompt("3E_Et_Ztag_EB_CRME",xbinsSS,datafile2,WZfile2,ZZfile2,"mergedEleCRanalyzer","20UL17");
  TH1D* subtracted_EB_SSdenom2 = subtractPrompt("3E_Et_Ztag_EB_antiME",xbinsSS,datafile2,WZfile2,ZZfile2,"mergedEleCRanalyzer","20UL17");

  TH1D* subtracted_EB_SSnumer3 = subtractPrompt("3E_Et_Ztag_EB_CRME",xbinsSS,datafile3,WZfile3,ZZfile3,"mergedEleCRanalyzer","20UL18");
  TH1D* subtracted_EB_SSdenom3 = subtractPrompt("3E_Et_Ztag_EB_antiME",xbinsSS,datafile3,WZfile3,ZZfile3,"mergedEleCRanalyzer","20UL18");

  subtracted_EB_SSnumer->Add(subtracted_EB_SSnumer1);
  subtracted_EB_SSnumer->Add(subtracted_EB_SSnumer2);
  subtracted_EB_SSnumer->Add(subtracted_EB_SSnumer3);
  subtracted_EB_SSdenom->Add(subtracted_EB_SSdenom1);
  subtracted_EB_SSdenom->Add(subtracted_EB_SSdenom2);
  subtracted_EB_SSdenom->Add(subtracted_EB_SSdenom3);

  TH1D* subtracted_Mu_SSnumer = subtractPrompt("2M_Et_Ztag_EB_CRME",xbinsSS,datafile4,WZfile4,ZZfile4,"mergedEMuCRanalyzer",postfix);
  TH1D* subtracted_Mu_SSdenom = subtractPrompt("2M_Et_Ztag_EB_antiME",xbinsSS,datafile4,WZfile4,ZZfile4,"mergedEMuCRanalyzer",postfix);

  TH1D* subtracted_Mu_SSnumer1 = subtractPrompt("2M_Et_Ztag_EB_CRME",xbinsSS,datafile5,WZfile5,ZZfile5,"mergedEMuCRanalyzer","20UL16APV");
  TH1D* subtracted_Mu_SSdenom1 = subtractPrompt("2M_Et_Ztag_EB_antiME",xbinsSS,datafile5,WZfile5,ZZfile5,"mergedEMuCRanalyzer","20UL16APV");

  TH1D* subtracted_Mu_SSnumer2 = subtractPrompt("2M_Et_Ztag_EB_CRME",xbinsSS,datafile6,WZfile6,ZZfile6,"mergedEMuCRanalyzer","20UL17");
  TH1D* subtracted_Mu_SSdenom2 = subtractPrompt("2M_Et_Ztag_EB_antiME",xbinsSS,datafile6,WZfile6,ZZfile6,"mergedEMuCRanalyzer","20UL17");

  TH1D* subtracted_Mu_SSnumer3 = subtractPrompt("2M_Et_Ztag_EB_CRME",xbinsSS,datafile7,WZfile7,ZZfile7,"mergedEMuCRanalyzer","20UL18");
  TH1D* subtracted_Mu_SSdenom3 = subtractPrompt("2M_Et_Ztag_EB_antiME",xbinsSS,datafile7,WZfile7,ZZfile7,"mergedEMuCRanalyzer","20UL18");

  subtracted_Mu_SSnumer->Add(subtracted_Mu_SSnumer1);
  subtracted_Mu_SSnumer->Add(subtracted_Mu_SSnumer2);
  subtracted_Mu_SSnumer->Add(subtracted_Mu_SSnumer3);
  subtracted_Mu_SSdenom->Add(subtracted_Mu_SSdenom1);
  subtracted_Mu_SSdenom->Add(subtracted_Mu_SSdenom2);
  subtracted_Mu_SSdenom->Add(subtracted_Mu_SSdenom3);

  subtracted_EB_SSnumer->Add(subtracted_Mu_SSnumer);
  subtracted_EB_SSdenom->Add(subtracted_Mu_SSdenom);

  subtracted_EB_SSnumer->Divide( subtracted_EB_SSdenom );

/*  TH1D* OSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_OSCR_EB_mixedME")->Clone();
  TH1D* OSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_OSCR_EB_antiME")->Clone();
  std::vector<double> xbinsOS = {0, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
                                 225, 250, 275, 300, 400, 500, 700, 1000};
  const int nbinsOS = xbinsOS.size()-1;
  std::vector<double> xcenOS = estimateCenter(xbinsOS);
  TH1D* OSnum_rebin = (TH1D*)OSnum->Rebin(nbinsOS, "2E_Et_OSCR_EB_mixedME", &(xbinsOS[0]));
  TH1D* OSdenom_rebin = (TH1D*)OSdenom->Rebin(nbinsOS, "2E_Et_OSCR_EB_antiME", &(xbinsOS[0]));

  OSnum_rebin->Divide( OSdenom_rebin );*/

  TF1* ssboth = new TF1("ssboth","[0]",50,500);
  ssboth->SetLineColor(kBlue);
  ssboth->SetLineWidth(2);
  ssboth->SetLineStyle(2);
  TFitResultPtr fitSS = subtracted_EB_SSnumer->Fit(ssboth,"RS");
  fitSS->SetName("fitSS");
  double ciSS[nbinsSS];
  fitSS->GetConfidenceIntervals(nbinsSS,1,0,&(xcenSS[0]),ciSS,0.68,false); // 0.6827
  std::vector<double> xbinwSS = estimateWidth(xbinsSS);
  double ybinSS[nbinsSS];

  for (unsigned idx = 0; idx < nbinsSS; idx++) {
    ybinSS[idx] = ssboth->Eval(xcenSS[idx]);
  }

  auto errSS = new TGraphErrors(nbinsSS,&(xcenSS[0]),ybinSS,&(xbinwSS[0]),ciSS);
  errSS->SetFillColor(kBlue);
  errSS->SetFillStyle(3004);

/*  TH1D* OSnumEta = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_eta_OSCR_EB_mixedME")->Clone();
  TH1D* OSdenomEta = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_eta_OSCR_EB_antiME")->Clone();

  OSnumEta->Divide( OSdenomEta );*/
  std::vector<double> xbinsEta = {-2.5,-1.5,-1.25,-1.,-0.75,-0.5,-0.25,0.,0.25,0.5,0.75,1.,1.25,1.5,2.5};
  std::vector<double> xcenEta = estimateCenter(xbinsEta);
  std::vector<double> xbinwEta = estimateWidth(xbinsEta);
  const int nbinsEta = xbinsEta.size()-1;

  TH1D* subtracted_etaE_SSnumer = subtractPrompt("3E_eta_Ztag_EB_CRME",xbinsEta,datafile,WZfile,ZZfile,"mergedEleCRanalyzer",postfix);
  TH1D* subtracted_etaE_SSdenom = subtractPrompt("3E_eta_Ztag_EB_antiME",xbinsEta,datafile,WZfile,ZZfile,"mergedEleCRanalyzer",postfix);

  TH1D* subtracted_etaE_SSnumer1 = subtractPrompt("3E_eta_Ztag_EB_CRME",xbinsEta,datafile1,WZfile1,ZZfile1,"mergedEleCRanalyzer","20UL16APV");
  TH1D* subtracted_etaE_SSdenom1 = subtractPrompt("3E_eta_Ztag_EB_antiME",xbinsEta,datafile1,WZfile1,ZZfile1,"mergedEleCRanalyzer","20UL16APV");

  TH1D* subtracted_etaE_SSnumer2 = subtractPrompt("3E_eta_Ztag_EB_CRME",xbinsEta,datafile2,WZfile2,ZZfile2,"mergedEleCRanalyzer","20UL17");
  TH1D* subtracted_etaE_SSdenom2 = subtractPrompt("3E_eta_Ztag_EB_antiME",xbinsEta,datafile2,WZfile2,ZZfile2,"mergedEleCRanalyzer","20UL17");

  TH1D* subtracted_etaE_SSnumer3 = subtractPrompt("3E_eta_Ztag_EB_CRME",xbinsEta,datafile3,WZfile3,ZZfile3,"mergedEleCRanalyzer","20UL18");
  TH1D* subtracted_etaE_SSdenom3 = subtractPrompt("3E_eta_Ztag_EB_antiME",xbinsEta,datafile3,WZfile3,ZZfile3,"mergedEleCRanalyzer","20UL18");

  subtracted_etaE_SSnumer->Add(subtracted_etaE_SSnumer1);
  subtracted_etaE_SSnumer->Add(subtracted_etaE_SSnumer2);
  subtracted_etaE_SSnumer->Add(subtracted_etaE_SSnumer3);
  subtracted_etaE_SSdenom->Add(subtracted_etaE_SSdenom1);
  subtracted_etaE_SSdenom->Add(subtracted_etaE_SSdenom2);
  subtracted_etaE_SSdenom->Add(subtracted_etaE_SSdenom3);

  TH1D* subtracted_etaM_SSnumer = subtractPrompt("2M_eta_Ztag_EB_CRME",xbinsEta,datafile4,WZfile4,ZZfile4,"mergedEMuCRanalyzer",postfix);
  TH1D* subtracted_etaM_SSdenom = subtractPrompt("2M_eta_Ztag_EB_antiME",xbinsEta,datafile4,WZfile4,ZZfile4,"mergedEMuCRanalyzer",postfix);

  TH1D* subtracted_etaM_SSnumer1 = subtractPrompt("2M_eta_Ztag_EB_CRME",xbinsEta,datafile5,WZfile5,ZZfile5,"mergedEMuCRanalyzer","20UL16APV");
  TH1D* subtracted_etaM_SSdenom1 = subtractPrompt("2M_eta_Ztag_EB_antiME",xbinsEta,datafile5,WZfile5,ZZfile5,"mergedEMuCRanalyzer","20UL16APV");

  TH1D* subtracted_etaM_SSnumer2 = subtractPrompt("2M_eta_Ztag_EB_CRME",xbinsEta,datafile6,WZfile6,ZZfile6,"mergedEMuCRanalyzer","20UL17");
  TH1D* subtracted_etaM_SSdenom2 = subtractPrompt("2M_eta_Ztag_EB_antiME",xbinsEta,datafile6,WZfile6,ZZfile6,"mergedEMuCRanalyzer","20UL17");

  TH1D* subtracted_etaM_SSnumer3 = subtractPrompt("2M_eta_Ztag_EB_CRME",xbinsEta,datafile7,WZfile7,ZZfile7,"mergedEMuCRanalyzer","20UL18");
  TH1D* subtracted_etaM_SSdenom3 = subtractPrompt("2M_eta_Ztag_EB_antiME",xbinsEta,datafile7,WZfile7,ZZfile7,"mergedEMuCRanalyzer","20UL18");

  subtracted_etaM_SSnumer->Add(subtracted_etaM_SSnumer1);
  subtracted_etaM_SSnumer->Add(subtracted_etaM_SSnumer2);
  subtracted_etaM_SSnumer->Add(subtracted_etaM_SSnumer3);
  subtracted_etaM_SSdenom->Add(subtracted_etaM_SSdenom1);
  subtracted_etaM_SSdenom->Add(subtracted_etaM_SSdenom2);
  subtracted_etaM_SSdenom->Add(subtracted_etaM_SSdenom3);

  subtracted_etaE_SSnumer->Add(subtracted_etaM_SSnumer);
  subtracted_etaE_SSdenom->Add(subtracted_etaM_SSdenom);

  subtracted_etaE_SSnumer->Divide( subtracted_etaE_SSdenom );

/*  TH2D* OSnum2d = (TH2D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_eta_OSCR_EB_mixedME")->Clone();
  TH2D* OSdenom2d = (TH2D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_eta_OSCR_EB_antiME")->Clone();

  OSdenom2d->RebinX(2);
  OSnum2d->RebinX(2);

  OSnum2d->Divide( OSdenom2d );

  TObjArray aSlices;
  TF1* os2d = new TF1("os2d","[0]+[1]*x+[2]/sqrt(x)",50,1000);
  os2d->SetParLimits(0,-0.5,0.5);
  os2d->SetParLimits(1,-0.001,0.001);
  OSnum2d->FitSlicesY(os2d,0,-1,0,"QNR",&aSlices);

  TF2* osboth = new TF2("osboth","[0]+[1]*y+[2]*sqrt(cosh(x)/y)",-1.4442,1.4442,50.,1000.);
  osboth->SetLineColor(kRed);
  osboth->SetLineWidth(2);
  osboth->SetLineStyle(2);
  TFitResultPtr fitOS = OSnum2d->Fit(osboth,"RS+");
  fitOS->SetName("fitOS");
  double ciOS[nbinsOS];
  std::vector<double> xycenOS;
  const double projectX = 0.75;
  TString formula_projected;
  formula_projected.Form("[0]+[1]*x+[2]*sqrt(cosh(%.2g)/x)",projectX);
  TF1* osboth_projected = new TF1("osboth_projected",formula_projected,50.,1000.);
  osboth_projected->SetLineColor(kRed);
  osboth_projected->SetLineWidth(2);
  osboth_projected->SetLineStyle(2);
  osboth_projected->SetParameters(osboth->GetParameters());

  for (const auto& xcen : xcenOS) {
    xycenOS.push_back(projectX);
    xycenOS.push_back(xcen);
  }

  fitOS->GetConfidenceIntervals(nbinsOS,2,1,&(xycenOS[0]),ciOS,0.95,false);
  std::vector<double> xbinwOS = estimateWidth(xbinsOS);
  double ybinOS[nbinsOS];

  for (unsigned idx = 0; idx < nbinsOS; idx++) {
    ybinOS[idx] = osboth->Eval(projectX,xcenOS[idx]);
  }

  auto errOS = new TGraphErrors(nbinsOS,&(xcenOS[0]),ybinOS,&(xbinwOS[0]),ciOS);
  errOS->SetFillColor(kRed);
  errOS->SetFillStyle(3005);*/

  // save file
  TFile* outfile = new TFile("MEFF_"+era+".root","RECREATE");
  subtracted_EB_SSnumer->Write();
//  OSnum2d->Write();
  ssboth->Write();
//  osboth->Write();
  fitSS->Write();
//  fitOS->Write();
  outfile->Close();

  auto* canvas = new TCanvas("canvas","canvas",50,50,W,H);

  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin( L/W );
  canvas->SetRightMargin( R/W );
  canvas->SetTopMargin( T/H );
  canvas->SetBottomMargin( B/H );
  canvas->SetTickx(0);
  canvas->SetTicky(0);

  // EB
//  auto legend = std::make_unique<TLegend>(0.8,0.75,0.95,0.9);
//  legend->SetBorderSize(0);
//  legend->AddEntry(subtracted_EB_SSnumer,"SS");
//  legend->AddEntry(OSnum_rebin,"OS");

  subtracted_EB_SSnumer->GetYaxis()->SetRangeUser(0.,2.0);
  subtracted_EB_SSnumer->SetLineWidth(2);
  subtracted_EB_SSnumer->GetYaxis()->SetTitle("Fake factor");
  subtracted_EB_SSnumer->GetXaxis()->SetTitle("E_{T} [GeV]");
  subtracted_EB_SSnumer->SetStats(0);
  subtracted_EB_SSnumer->SetLineColor(kBlue);
  subtracted_EB_SSnumer->Draw("E1");
  errSS->Draw("3");
//  errOS->Draw("3");
//  legend->Draw();

  TPaveText* textlow = new TPaveText(0.75,0.8,0.95,0.95,"NDC");
  textlow->SetBorderSize(0);
  textlow->SetFillStyle(3025);
  textlow->SetFillColor(0);
  TString textsslow;
  textsslow.Form(" %.3f#pm%.3f", ssboth->GetParameter(0), ssboth->GetParError(0));
  textlow->AddText(textsslow);
  ((TText*)textlow->GetListOfLines()->Last())->SetTextColor(kBlue);
  ((TText*)textlow->GetListOfLines()->Last())->SetTextAlign(12);

  textlow->Draw();

  canvas->Update();

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("MEFF_Et_"+era+".pdf");

  subtracted_etaE_SSnumer->GetYaxis()->SetRangeUser(0.,2.0);
  subtracted_etaE_SSnumer->GetYaxis()->SetTitle("Fake factor");
  subtracted_etaE_SSnumer->GetXaxis()->SetTitle("#eta_{SC}");
  subtracted_etaE_SSnumer->SetLineWidth(2);
  subtracted_etaE_SSnumer->SetLineColor(kBlue);
  subtracted_etaE_SSnumer->Draw("E1");
//  legend->Draw();

  TF1* ssEta = new TF1("ssEta","[0]",-1.4442,1.4442);
  ssEta->SetLineColor(kBlue);
  ssEta->SetLineWidth(2);
  ssEta->SetLineStyle(2);
  TFitResultPtr fitSSEta = subtracted_etaE_SSnumer->Fit(ssEta,"RS");
  ssEta->Draw("same");

  double etaSS[2] = {-1.4442,1.4442};
  double ciSSeta[nbinsEta];
  fitSSEta->GetConfidenceIntervals(nbinsEta,1,0,&(xcenEta[0]),ciSSeta,0.68,false); // 0.6827

  double ybinSSeta[2] = {ssEta->GetParameter(0),ssEta->GetParameter(0)};
  auto errSSeta = new TGraphErrors(2,etaSS,ybinSSeta,&(xbinwEta[0]),ciSSeta);
  errSSeta->SetFillColor(kBlue);
  errSSeta->SetFillStyle(3004);
  errSSeta->Draw("3");

  TPaveText* textlowEta = new TPaveText(0.75,0.8,0.95,0.95,"NDC");
  textlowEta->SetBorderSize(0);
  textlowEta->SetFillStyle(3025);
  textlowEta->SetFillColor(0);
  TString textsslowEta;
  textsslowEta.Form(" %.3f#pm%.3f", ssEta->GetParameter(0), ssEta->GetParError(0));
  textlowEta->AddText(textsslowEta);
  ((TText*)textlowEta->GetListOfLines()->Last())->SetTextColor(kBlue);
  ((TText*)textlowEta->GetListOfLines()->Last())->SetTextAlign(12);

  textlowEta->Draw();

  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("MEFF_eta_"+era+".pdf");

  CMS_lumi( canvas, iPeriod, iPos );

  return;
}
