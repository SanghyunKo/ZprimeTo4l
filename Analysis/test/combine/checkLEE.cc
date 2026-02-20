#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

#include <algorithm>

#include "/u/user/sako/ModHEEP/CMSSW_10_6_29/work/tdrstyle.C"
#include "/u/user/sako/ModHEEP/CMSSW_10_6_29/work/CMS_lumi.C"

void checkLEE() {
  setTDRStyle();
  // gStyle->SetLineWidth(2);
  TString era = "run2";

  writeExtraText = true;       // if extra text
  extraText  = "   Internal";  // default extra text is "Preliminary"

  if (era=="20UL16APV") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "19.5 fb^{-1}";
  } else if (era=="20UL16") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "16.8 fb^{-1}";
  } else if (era=="20UL17") {
    lumi_sqrtS = "2017 (13 TeV)";
    lumi_13TeV = "41.48 fb^{-1}";
  } else if (era=="20UL18") {
    lumi_sqrtS = "2018 (13 TeV)";
    lumi_13TeV = "59.83 fb^{-1}";
  } else if (era=="run2") {
    lumi_sqrtS = "";
    lumi_13TeV = "137.6 fb^{-1}";
  } else {
    std::cout << "check era..." << std::endl;
  }

  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 0;

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

  auto* canvas1 = new TCanvas("canvas","canvas",50,50,W,H);
  canvas1->SetFillColor(0);
  canvas1->SetBorderMode(0);
  canvas1->SetFrameFillStyle(0);
  canvas1->SetFrameBorderMode(0);
  canvas1->SetLeftMargin( L/W );
  canvas1->SetRightMargin( R/W );
  canvas1->SetTopMargin( T/H );
  canvas1->SetBottomMargin( B/H );
  canvas1->SetTickx(0);
  canvas1->SetTicky(0);

  auto SaveAs = [&] (TCanvas* canvas, const std::string& name, TPad* pad = nullptr) {
    canvas->Update();

    // writing the lumi information and the CMS "logo"
    CMS_lumi( canvas, iPeriod, iPos );

    if (pad) {
      pad->RedrawAxis();
      pad->GetFrame()->Draw();
    } else {
      canvas->Update();
      canvas->RedrawAxis();
      canvas->GetFrame()->Draw();
    }

    canvas->SaveAs(name.c_str());
  };

  auto retrieve = [] (int mX, int ntoy) {
    TString strMX = TString(std::to_string(mX));
    TString strNtoy = TString(std::to_string(ntoy));
    TString filename = TString("./filesLEE/")+"higgsCombine.LEEmX"+strMX+"toy"+strNtoy+".HybridNew.mH"+strMX+".123456.root";

    auto afile = std::make_unique<TFile>(filename,"READ");

    if (afile->IsZombie() || afile->TestBit(TFile::kRecovered))
      return 1.;

    TTree* atree = (TTree*)afile->Get("limit");
    double val;
    atree->SetBranchAddress("limit",&val);

    unsigned nentry = atree->GetEntries();

    if (nentry==0)
      return 1.;

    atree->GetEntry(0);
    double valFinal = val;

    return valFinal;
  };

  std::vector<int> massvec = {250,275,300,325,350,375,400,425,450,500,550,650,750,850,1000,1250,1500,1750,2000};
  std::vector<double> massvecDb = {250.,275.,300.,325.,350.,375.,400.,425.,450.,500.,550.,650.,750.,850.,1000.,1250.,1500.,1750.,2000.};

  auto toyPvalues = [&retrieve,&massvec] (int ntoy) -> std::vector<double> {
    std::vector<double> pvalues;

    for (unsigned idx=0; idx<massvec.size(); idx++) {
      double val = retrieve(massvec.at(idx),ntoy);

      pvalues.push_back(val);
    }

    return pvalues;
  };

  int totalNtoys = 2000;
  int lowPvalToy = 0;

  for (int idx=1; idx<=totalNtoys; idx++) {
    auto vals = toyPvalues(idx);

    double minPval = *std::min_element(vals.begin(),vals.end());

    if (minPval < 0.009025) {
      std::cout << "Ntoy = " << idx << ", min p-val = " << minPval << std::endl;
      lowPvalToy++;

      auto gr = std::make_unique<TGraph>(vals.size(),&(massvecDb[0]),&(vals[0]));
      gr->SetMaximum(1.5);
      gr->SetMinimum(10e-5);

      gr->Draw("AL");
      canvas1->SetLogy();

      gr->SetLineColor(kBlack);
      gr->SetMarkerSize(0);
      gr->GetXaxis()->SetTitle("M_{X} [GeV]");
      gr->GetXaxis()->SetTitleSize(0.05);
      gr->GetYaxis()->SetTitle("p-value w.r.t. bkg-only toy");
      gr->GetYaxis()->SetTitleOffset(1.3);
      gr->GetXaxis()->SetTitleOffset(1.0);
      gr->GetYaxis()->SetTitleSize(0.04);
      gr->GetYaxis()->SetLabelSize(0.04);
      gr->GetXaxis()->SetLabelSize(0.04);
      //SaveAs(canvas1,std::string("plotLEE_")+std::to_string(idx)+".png");
    }
  }

  std::cout << "Finished! Total " << lowPvalToy << " / 2000 toys have p-value lower than 0.009025" << std::endl;
}
