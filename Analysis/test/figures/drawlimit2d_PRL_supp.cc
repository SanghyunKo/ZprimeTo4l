#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

// #include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
// #include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

#include "/u/user/sako/ModHEEP/CMSSW_10_6_29/work/tdrstyle.C"
#include "/u/user/sako/ModHEEP/CMSSW_10_6_29/work/CMS_lumi.C"

void drawlimit2d_PRL_supp() {
  TString era="run2";
  setTDRStyle();
  // gStyle->SetLineWidth(2);

  writeExtraText = true;       // if extra text
  extraText  = "";  // default extra text is "Preliminary"
  customLumiOffset = 0.01;
  customCmsTextOffset = 0.01;

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
    lumi_13TeV = "138 fb^{-1}";
  } else {
    std::cout << "check era..." << std::endl;
  }

  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 0;

  if( iPos==0 )
    relPosX = 0.25; // doesn't work

  int W = 1600;
  int H = 1400;

  int H_ref = 1050;
  int W_ref = 1200;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.1*H_ref;
  float L = 0.1*W_ref;
  float R = 0.16*W_ref;

  auto* canvas1 = new TCanvas("canvas","canvas",1600,1600,W,H);
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

  auto retrieve = [] (const TString& hmass, const TString& amass, const TString& isZA, const TString& category, TString quant="") -> double {
    TString filename = TString("./")+isZA+amass+"_"+category+"/"+hmass+"/higgsCombine.X"+hmass+isZA+amass+".HybridNew.mH"+hmass;

    if (quant!="")
      filename += ".quant"+quant;

    filename += ".root";

    auto afile = std::make_unique<TFile>(filename,"READ");

    if (afile->IsZombie() || afile->TestBit(TFile::kRecovered))
      return -1.;

    TTree* atree = (TTree*)afile->Get("limit");
    double val;
    atree->SetBranchAddress("limit",&val);

    unsigned nentry = atree->GetEntries();

    if (nentry==0)
      return 9999.;

    atree->GetEntry(0);
    double valFinal = val;

    return 0.1*valFinal;
  };

  auto retrieveSignificance = [] (const TString& hmass, const TString& amass, const TString& isZA, const TString& category) -> double {
    TString filename = TString("./")+isZA+amass+"_"+category+"/"+hmass+"/higgsCombine.significanceCWR.X"+hmass+isZA+amass+".readSignificance.HybridNew.mH"+hmass+".root";

    auto afile = std::make_unique<TFile>(filename,"READ");

    if (afile->IsZombie() || afile->TestBit(TFile::kRecovered))
      return -9999.;

    TTree* atree = (TTree*)afile->Get("limit");
    double val;
    atree->SetBranchAddress("limit",&val);

    unsigned nentry = atree->GetEntries();

    if (nentry==0)
      return -9999.;

    atree->GetEntry(0);
    double valFinal = val;

    return valFinal;
  };

  auto checkKinematic = [] (double xmass, double ymass, const TString& isZA) -> bool {
    if (isZA=="A") {
      if (ymass < xmass/2.)
        return true;
    } else if (isZA=="Z") {
      if (ymass < xmass - 91.1876)
        return true;
    }

    return false;
  };

  auto fillY = [&checkKinematic] (TH2D* hist, const TString& isZA) {
    for (int idx=2; idx<=hist->GetNbinsX(); idx++) {
      for (int jdx=2; jdx<=hist->GetNbinsY(); jdx++) {
        double lowX = hist->GetXaxis()->GetBinLowEdge(idx);
        double lowY = hist->GetYaxis()->GetBinLowEdge(jdx);

        if ( !checkKinematic(lowX,lowY,isZA) )
          continue;

        double con = hist->GetBinContent(idx,jdx);

        double conp1 = hist->GetBinContent(idx,jdx+1);
        double conp2 = hist->GetBinContent(idx,jdx+2);
        double conm1 = hist->GetBinContent(idx,jdx-1);

        if (con == 0. && (conp1 > 0. || conp2 > 0.) && conm1 > 0.) {
          hist->SetBinContent( hist->GetBin(idx,jdx), conm1 );
        }
      }
    }
  };

  auto fillX = [&checkKinematic] (TH2D* hist, const TString& isZA) {
    for (int idx=2; idx<=hist->GetNbinsX(); idx++) {
      for (int jdx=2; jdx<=hist->GetNbinsY(); jdx++) {
        double lowX = hist->GetXaxis()->GetBinLowEdge(idx);
        double lowY = hist->GetYaxis()->GetBinLowEdge(jdx);

        if ( !checkKinematic(lowX,lowY,isZA) )
          continue;

        double con = hist->GetBinContent(idx,jdx);

        double conp1 = hist->GetBinContent(idx+1,jdx);
        double conp2 = hist->GetBinContent(idx+2,jdx);
        double conm1 = hist->GetBinContent(idx-1,jdx);

        if (con == 0. && (conp1 > 0. || conp2 > 0.) && conm1 > 0.) {
          hist->SetBinContent( hist->GetBin(idx,jdx), conm1 );
        }
      }
    }
  };

  std::vector<int> massvec = {250,275,300,325,350,375,400,425,450,500,550,650,750,850,1000,1250,1500,1750,2000};
  std::vector<TString> avecA = {"0p4","0p6","0p8","1","1p5","2","5","10","50","100","250","500","750"};
  std::vector<TString> avecZ = {"0p4","0p6","0p8","1","1p5","2","5","10","50","100","250","500","750","1000","1500"};

  // std::vector<float> xbin = {250., 275., 300., 325., 350., 375., 400., 425., 450., 500., 550., 650., 750., 850., 1000., 1250., 1500., 1750., 2000., 2500.};
  std::vector<float> xbin = {237.5,262.5,287.5,312.5,337.5,362.5,387.5,412.5,437.5,475., 525., 575., 700., 800., 925.,  1125., 1375., 1625., 1875., 2125.};
  std::vector<float> ybinA = {0.4,0.6,0.8,1.,1.5,2.,5.,10.,50.,100.,250.,500.,750.,1250.};
  std::vector<float> ybinZ = {0.4,0.6,0.8,1.,1.5,2.,5.,10.,50.,100.,250.,500.,750.,1000.,1500.,2000.};
  // std::vector<float> ybinForHist = {0.4, 0.6, 0.8, 1., 1.5, 2.,  5., 10., 50., 100., 250., 500., 750., 1000., 1500., 2000., 5000.};
  std::vector<float> ybinForHist = {0.35,0.5, 0.7, 0.9,1.25,1.75,2.5,7.5, 25., 75.,  125., 350., 600., 875.,  1250., 1750., 5000.};

  class SubFigure {
  public:
    SubFigure()=default;
    ~SubFigure()=default;

    TH2D* limit2d_ = nullptr;
    TH2D* sig2d_ = nullptr;
    TH2D* sig2d2sig_ = nullptr;
    TString isZA_;
    TString isMu_;
    double limit_max_ = -1.;
    double limit_min_ = std::numeric_limits<double>::max();
    double sig_resolved_max_ = -1.;
    double xmax_resolved_;
    double ymax_resolved_;
    double sig_merged_max_ = -1.;
    double xmax_merged_;
    double ymax_merged_;

    // for HepData
    TTree* limitTree_ = nullptr;
    float mX_, mY_, limit_, significance_;
    float nominal_, quant0p025_, quant0p975_, quant0p16_, quant0p84_;

    SubFigure(const std::vector<float>& xbin, const std::vector<float>& ybin, const TString& isZA, const TString& isMu) {
      limit2d_ = new TH2D("limit2d"+isZA+"_"+isMu,"limit2d",xbin.size()-1,&(xbin[0]),ybin.size()-1,&(ybin[0]));
      sig2d_ = new TH2D("sig2d"+isZA+"_"+isMu,"sig2d",xbin.size()-1,&(xbin[0]),ybin.size()-1,&(ybin[0]));
      isZA_ = isZA;
      isMu_ = isMu;

      limitTree_ = new TTree("limitTree"+isZA+"_"+isMu,"limitTree");
    }
    
    void setupBranches() { /// to be called after vector construction to ensure objects are in final memory locations
      limitTree_->Branch("mX",&mX_,"mX/F");
      limitTree_->Branch("mY",&mY_,"mY/F");
      limitTree_->Branch("limit",&limit_,"limit/F");
      limitTree_->Branch("significance",&significance_,"significance/F");
      limitTree_->Branch("nominal",&nominal_,"nominal/F");
      limitTree_->Branch("quant0p025",&quant0p025_,"quant0p025/F");
      limitTree_->Branch("quant0p975",&quant0p975_,"quant0p975/F");
      limitTree_->Branch("quant0p16",&quant0p16_,"quant0p16/F");
      limitTree_->Branch("quant0p84",&quant0p84_,"quant0p84/F");
    }
  };

  std::vector<SubFigure> subFigures = {
    SubFigure(xbin, ybinForHist, "A", "el"),
    SubFigure(xbin, ybinForHist, "A", "mu"),
    SubFigure(xbin, ybinForHist, "Z", "el"),
    SubFigure(xbin, ybinForHist, "Z", "mu")
  };

  // Setup branch addresses after vector construction to ensure objects are in final memory locations
  for (auto& subFigure : subFigures) {
    subFigure.setupBranches();
  }

  // for the sake of sorting HepData tables in mX
  for (auto& subFigure : subFigures) {
    std::vector<TString> avec = (subFigure.isZA_=="A") ? avecA : avecZ;

    for (int idx = 0; idx<avec.size(); idx++) {
      for (const auto hmass : massvec) {
        auto ymassStr = avec.at(idx);
        double ymass = (subFigure.isZA_=="A") ? ybinA.at(idx) : ybinZ.at(idx);

        if (!checkKinematic(static_cast<double>(hmass),ymass,subFigure.isZA_))
          continue;

        double val = retrieve(std::to_string(hmass),ymassStr,subFigure.isZA_, subFigure.isMu_);
        double valSignificance = retrieveSignificance(std::to_string(hmass),ymassStr,subFigure.isZA_, subFigure.isMu_);
        int ibin = subFigure.limit2d_->FindFixBin(static_cast<float>(hmass)+1.,ymass+0.001);

        if (valSignificance > 2.)
          std::cout << subFigure.isZA_ << " significance > 2 sigma found at Mx=" << hmass << " My=" << ymass << " : " << valSignificance << std::endl;

        if (val > 0.) {
          subFigure.limit2d_->SetBinContent(ibin, val );
          subFigure.limit_max_ = val>subFigure.limit_max_ ? val : subFigure.limit_max_;
          subFigure.limit_min_ = val<subFigure.limit_min_ ? val : subFigure.limit_min_;

          subFigure.sig2d_->SetBinContent(ibin, std::max(valSignificance, 0.));

          if (ymass >= 5.) {
            if (valSignificance > subFigure.sig_resolved_max_) {
              subFigure.sig_resolved_max_ = valSignificance;
              subFigure.xmax_resolved_ = hmass;
              subFigure.ymax_resolved_ = ymass;
            }
          }

          if (ymass < 5.) {
            if (valSignificance > subFigure.sig_merged_max_) {
              subFigure.sig_merged_max_ = valSignificance;
              subFigure.xmax_merged_ = hmass;
              subFigure.ymax_merged_ = ymass;
            }
          }

          // fill HepData tree
          subFigure.mX_ = static_cast<float>(hmass);
          subFigure.mY_ = static_cast<float>(ymass);
          subFigure.limit_ = static_cast<float>(val);
          subFigure.significance_ = std::max(valSignificance, 0.);
          subFigure.nominal_ = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,subFigure.isZA_, subFigure.isMu_,"0.500"));
          subFigure.quant0p025_ = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,subFigure.isZA_, subFigure.isMu_,"0.025")) - subFigure.nominal_;
          subFigure.quant0p975_ = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,subFigure.isZA_, subFigure.isMu_,"0.975")) - subFigure.nominal_;
          subFigure.quant0p16_ = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,subFigure.isZA_, subFigure.isMu_,"0.160")) - subFigure.nominal_;
          subFigure.quant0p84_ = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,subFigure.isZA_, subFigure.isMu_,"0.840")) - subFigure.nominal_;
          subFigure.limitTree_->Fill();
        } // if val > 0.
      }
    }
  }

  double zmin = 0.03, zmax = 200.;

  for (auto& subFigure : subFigures) {
    std::cout << subFigure.isZA_ << "_" << subFigure.isMu_ << " limit range: " << subFigure.limit_min_ << " - " << subFigure.limit_max_ << std::endl;
    std::cout << "resolved sig max " << subFigure.sig_resolved_max_ << " at Mx: " << subFigure.xmax_resolved_ << " My: " << subFigure.ymax_resolved_ << std::endl;
    std::cout << "merged sig max " << subFigure.sig_merged_max_ << " at Mx: " << subFigure.xmax_merged_ << " My: " << subFigure.ymax_merged_ << std::endl;

    fillY(subFigure.limit2d_, subFigure.isZA_);
    fillX(subFigure.limit2d_, subFigure.isZA_);

    fillY(subFigure.sig2d_, subFigure.isZA_);
    fillX(subFigure.sig2d_, subFigure.isZA_);

    subFigure.limit2d_->GetXaxis()->SetLabelSize(0.07);
    subFigure.limit2d_->GetYaxis()->SetLabelSize(0.07);
    subFigure.limit2d_->GetZaxis()->SetRangeUser(zmin,zmax); // define z-scale here
  }

  gStyle->SetPalette(kViridis);

  // multipad
  // Define the rows and columns for subpads
  int rows = 2;
  int cols = 2;

  // Create subpads and divide the canvas
  double leftmargin = L/W;
  double rightmargin = R/W;
  double topmargin = T/H;
  double bottommargin = 1.2*B/H;

  canvas1->SetLogy();
  // canvas1->SetLogx();
  canvas1->SetLogz();
  
  // loop over pads
  for (int j = 0; j < rows; ++j) {
    for (int i = 0; i < cols; ++i) {
      unsigned iPad = j*cols + i;
      unsigned iFig = iPad;

      // Calculate the position of each subpad
      double xmin = leftmargin + static_cast<double>(i)*(1.-leftmargin-rightmargin)/cols;
      double xmax = leftmargin + static_cast<double>(i+1)*(1.-leftmargin-rightmargin)/cols;
      double ymin = bottommargin + static_cast<double>(rows-j-1)*(1.-topmargin-bottommargin)/rows;
      double ymax = bottommargin + static_cast<double>(rows-j)*(1.-topmargin-bottommargin)/rows;

      if (i==0)
        xmin -= 0.6*leftmargin;
      if (j==rows-1)
        ymin = 0.01;

      if (i==0)
        xmax -= 0.0;
      if (i==1)
        xmin -= 0.0;

      if (j==0)
        ymin -= 0.02;
      if (j==1)
        ymax -= 0.02;

      // Create subpad
      TPad *padUp = new TPad(Form("padUp%d", iPad), Form("SubpadUp %d", iPad), xmin, ymin, xmax, ymax);
      padUp->SetLeftMargin( 0. );
      padUp->SetRightMargin( 0. );
      padUp->SetTopMargin( 0.0 );
      padUp->SetBottomMargin( j==0 ? 0.0 : 1.2*bottommargin );
      padUp->SetFillColor(0);
      padUp->SetBorderMode(0);
      padUp->SetFrameFillStyle(0);
      padUp->SetFrameBorderMode(0);

      if (i==0)
        padUp->SetLeftMargin(1.6*leftmargin);

      padUp->Draw();
      padUp->cd();
      padUp->SetLogy();
      padUp->SetLogz();

      auto* limitHist = subFigures.at(iFig).limit2d_;
      auto* significanceHist1sig = subFigures.at(iFig).sig2d_;
      const auto& isZA = subFigures.at(iFig).isZA_;
      const auto& isMu = subFigures.at(iFig).isMu_;

      auto* significanceHist2sig = (TH2D*)significanceHist1sig->Clone("significance2sig"+isZA+"_"+isMu);
      significanceHist1sig->SetContour(1,&(std::vector<double>{1.} )[0]);
      significanceHist1sig->SetLineColor(TColor::GetColor("#f89c20"));
      significanceHist1sig->SetLineWidth(2);
      significanceHist2sig->SetContour(1,&(std::vector<double>{2.} )[0]);
      significanceHist2sig->SetLineColor(TColor::GetColor("#e42536"));
      significanceHist2sig->SetLineWidth(2);
      subFigures.at(iFig).sig2d2sig_ = significanceHist2sig;

      // TString limitDrawOption = (iPad==3) ? "colz" : "col";

      if (j==1)
        limitHist->GetXaxis()->ChangeLabel(-1,-1,0.);

      limitHist->GetXaxis()->SetNdivisions(505);

      limitHist->Draw("col");
      significanceHist1sig->Draw("cont3 same");
      significanceHist2sig->Draw("cont3 same");

      double xoffset = (i==0) ? 0.18 : 0.18;

      auto* atext = new TPaveText(0.6+xoffset,0.87,0.8+xoffset,0.99,"NDC");
      atext->SetBorderSize(0);
      atext->SetFillColor(0);
      atext->SetFillStyle(0);
      TString channelName = (isZA=="A") ? "\\mathrm{X}\\rightarrow\\mathrm{YY}" : "\\mathrm{X}\\rightarrow\\mathrm{ZY}";
      channelName += (isMu=="mu") ? "\\:(\\mathrm{Y}\\rightarrow\\mu\\mu)" : "\\:(\\mathrm{Y}\\rightarrow\\mathrm{ee})";
      atext->AddText(channelName);
      ((TText*)atext->GetListOfLines()->Last())->SetTextColor(kBlack);
      ((TText*)atext->GetListOfLines()->Last())->SetTextAlign(31);
      ((TText*)atext->GetListOfLines()->Last())->SetTextSize(0.08);
      ((TText*)atext->GetListOfLines()->Last())->SetTextFont(42);
      atext->Draw();

      canvas1->cd();  // Go back to the main canvas
    }
  }

  canvas1->cd();
  canvas1->Update();  // Refresh canvas

  // create a dedicated small pad on the right that reserves space for labels/title
  canvas1->cd();
  canvas1->Update();
  TPad *padPal = new TPad("padPalette","padPalette",0.88,0.06,0.95,0.94,0);
  padPal->SetFillStyle(0);
  padPal->SetBorderMode(0);
  padPal->SetFrameBorderMode(0);
  padPal->SetLeftMargin(0.0);
  padPal->SetRightMargin(0.5);
  padPal->SetTopMargin(0.0);
  padPal->SetBottomMargin(0.0);
  padPal->Draw();
  padPal->cd();
  gPad->Update();
  // TPaletteAxis *palette = (TPaletteAxis*)subFigures.back().limit2d_->GetListOfFunctions()->FindObject("palette");
  TPaletteAxis *palette = new TPaletteAxis(0.05,0.02,0.3,0.98, zmin, zmax); //subFigures.back().limit2d_);
  canvas1->Modified();
  canvas1->Update();
  // coordinates here are normalized to padPal (0..1)
  palette->SetX1NDC(0.05);
  palette->SetX2NDC(0.3);
  palette->SetY1NDC(0.02);
  palette->SetY2NDC(0.98);
  palette->Draw();

  double px1 = palette->GetX1NDC();
  double px2 = palette->GetX2NDC();
  double py1 = palette->GetY1NDC();
  double py2 = palette->GetY2NDC();
  double ax = px2;
  // create a vertical axis mapping zmin..zmax to visually match the palette
  TGaxis *gaxis = new TGaxis(ax, py1, ax, py2, zmin, zmax, 50510, "+GLS");
  gaxis->SetLabelFont(42);
  gaxis->SetLabelSize(0.4);    // tune to taste
  // gaxis->SetTitle("Z scale");
  // gaxis->SetTitleSize(0.5);
  // gaxis->SetTitleOffset(0.05);
  gaxis->SetTickSize(0.15);
  gaxis->Draw();
  padPal->Update();
  canvas1->cd();
  canvas1->Modified();
  canvas1->Update();

  auto legend = new TLegend(0.088,0.83,0.8,0.93);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetMargin(0.12);
  legend->SetNColumns(2);
  legend->AddEntry(subFigures.front().sig2d_,"1 s.d. obs. local significance         ","l");
  legend->AddEntry(subFigures.front().sig2d2sig_,"2 s.d. obs. local significance","l");
  canvas1->cd();
  legend->Draw();

  auto ytext = std::make_unique<TPaveText>(0.03,0.68,0.06,0.9,"NDC");
  ytext->SetBorderSize(0);
  ytext->SetFillColor(0);
  ytext->SetFillStyle(0);
  TString astr = "M_{Y} [GeV]";
  ytext->AddText(astr);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextFont(42);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextSize(0.04);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextAngle(90);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextAlign(11);
  ytext->Draw();

  auto xtext = std::make_unique<TPaveText>(0.8,0.03,0.95,0.08,"NDC");
  xtext->SetBorderSize(0);
  xtext->SetFillColor(0);
  xtext->SetFillStyle(0);
  TString xstr = "M_{X} [GeV]";
  xtext->AddText(xstr);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextFont(42);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextAlign(33);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextSize(0.04);
  xtext->Draw();

  auto ztext = std::make_unique<TPaveText>(0.975,0.35,0.99,0.85,"NDC");
  ztext->SetBorderSize(0);
  ztext->SetFillColor(0);
  ztext->SetFillStyle(0);
  TString zstr = "\\sigma(\\mathrm{pp}\\rightarrow\\mathrm{X})\\times\\mathrm{B}\\:[\\mathrm{fb}]";
  ztext->AddText(zstr);
  ((TText*)ztext->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)ztext->GetListOfLines()->Last())->SetTextFont(42);
  ((TText*)ztext->GetListOfLines()->Last())->SetTextAlign(11);
  ((TText*)ztext->GetListOfLines()->Last())->SetTextSize(0.04);
  ((TText*)ztext->GetListOfLines()->Last())->SetTextAngle(90);
  ztext->Draw();

  auto ztext2 = std::make_unique<TPaveText>(0.9,0.4,0.952,0.77,"NDC");
  ztext2->SetBorderSize(0);
  ztext2->SetFillColor(0);
  ztext2->SetFillStyle(0);
  TString zstr2 = "Upper limits on";
  ztext2->AddText(zstr2);
  ((TText*)ztext2->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)ztext2->GetListOfLines()->Last())->SetTextFont(42);
  ((TText*)ztext2->GetListOfLines()->Last())->SetTextAlign(33);
  ((TText*)ztext2->GetListOfLines()->Last())->SetTextSize(0.045);
  ((TText*)ztext2->GetListOfLines()->Last())->SetTextAngle(90);
  ztext2->Draw();

  SaveAs(canvas1,std::string(TString("limit2D_PRL_supp.png").Data()));

  // HepData output
  auto outFile = std::make_unique<TFile>("HepData_limit2D_PRL_supp.root","RECREATE");

  for (auto& subFigure : subFigures)
    subFigure.limitTree_->Write();

  outFile->Close();

  // save individual canvases for each channel
  // create new canvas
  customCmsTextOffset = 0.;
  customLumiOffset = 0.;

  iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  iPos = 0;
  relPosX = 0.045;

  if( iPos==0 )
    relPosX = 0.12;

  W = 800;
  H = 800;

  H_ref = 800;
  W_ref = 800;

  // references for T, B, L, R
  T = 0.08*H_ref;
  B = 0.08*H_ref;
  L = 0.12*W_ref; // 0.1
  R = 0.16*W_ref; // 0.02

  auto* canvas2 = new TCanvas("canvas2","canvas2",W,H,W,H);
  canvas2->SetFillColor(0);
  canvas2->SetBorderMode(0);
  canvas2->SetFrameFillStyle(0);
  canvas2->SetFrameBorderMode(0);
  canvas2->SetLeftMargin( L/W );
  canvas2->SetRightMargin( R/W );
  canvas2->SetTopMargin( T/H );
  canvas2->SetBottomMargin( B/H );
  canvas2->SetTickx(0);
  canvas2->SetTicky(0);
  canvas2->SetLogy();
  canvas2->SetLogz();
  canvas2->cd();

  for (unsigned iFig=0; iFig<subFigures.size(); iFig++) {
    auto subFigure = subFigures.at(iFig);
    auto limitHist = subFigure.limit2d_;
    auto significanceHist1sig = subFigure.sig2d_;
    auto significanceHist2sig = subFigure.sig2d2sig_;
    auto isZA = subFigure.isZA_;
    auto isMu = subFigure.isMu_;

    // if (i==0)
    // limitHist->GetXaxis()->ChangeLabel(-1,-1,0.);
    limitHist->GetXaxis()->SetNdivisions(505);

    limitHist->GetXaxis()->SetLabelSize(0.035);
    limitHist->GetXaxis()->SetTitle("M_{X} [GeV]");
    limitHist->GetXaxis()->SetTitleSize(0.04);
    limitHist->GetYaxis()->SetLabelSize(0.04);
    limitHist->GetYaxis()->SetTitle("M_{Y} [GeV]");
    limitHist->GetYaxis()->SetTitleSize(0.04);
    limitHist->GetZaxis()->SetLabelSize(0.035);
    limitHist->GetZaxis()->SetTitleSize(0.04);
    limitHist->GetZaxis()->SetTitle("\\sigma(\\mathrm{pp}\\rightarrow\\mathrm{X}) \\times \\mathrm{B} [\\mathrm{fb}]");
    limitHist->GetZaxis()->SetTitleOffset(1.25);

    limitHist->Draw("colz");
    significanceHist1sig->Draw("cont3 same");
    significanceHist2sig->Draw("cont3 same");

    auto legend2 = new TLegend(0.13,0.8,0.7,0.92);
    legend2->SetBorderSize(0);
    legend2->SetFillColor(0);
    legend2->SetFillStyle(0);
    legend2->SetTextFont(42);
    legend2->SetTextSize(0.035);
    legend2->SetMargin(0.12);
    legend2->AddEntry(significanceHist1sig,"1 s.d. obs. local sig.","l");
    legend2->AddEntry(significanceHist2sig,"2 s.d. obs. local sig.","l");
    legend2->Draw();

    auto* atext = new TPaveText(0.73,0.85,0.83,0.9,"NDC");
    atext->SetBorderSize(0);
    atext->SetFillColor(0);
    atext->SetFillStyle(0);
    TString channelName = (isZA=="A") ? "\\mathrm{X}\\rightarrow\\mathrm{YY}" : "\\mathrm{X}\\rightarrow\\mathrm{ZY}";
    channelName += (isMu=="mu") ? "\\:(\\mathrm{Y}\\rightarrow\\mu\\mu)" : "\\:(\\mathrm{Y}\\rightarrow\\mathrm{ee})";
    atext->AddText(channelName);
    ((TText*)atext->GetListOfLines()->Last())->SetTextColor(kBlack);
    ((TText*)atext->GetListOfLines()->Last())->SetTextAlign(31);
    ((TText*)atext->GetListOfLines()->Last())->SetTextSize(0.045);
    ((TText*)atext->GetListOfLines()->Last())->SetTextFont(42);
    atext->Draw();

    auto btext = std::make_unique<TPaveText>(0.9,0.4,0.956,0.8,"NDC");
    btext->SetBorderSize(0);
    btext->SetFillColor(0);
    btext->SetFillStyle(0);
    TString bstr = "Upper limits on";
    btext->AddText(bstr);
    ((TText*)btext->GetListOfLines()->Last())->SetTextColor(kBlack);
    ((TText*)btext->GetListOfLines()->Last())->SetTextFont(42);
    ((TText*)btext->GetListOfLines()->Last())->SetTextAlign(33);
    ((TText*)btext->GetListOfLines()->Last())->SetTextSize(0.04);
    ((TText*)btext->GetListOfLines()->Last())->SetTextAngle(90);
    btext->Draw();

    canvas2->Update();  // Refresh canvas
    SaveAs(canvas2,std::string((TString("HepData_limit2D_")+isZA+"_"+isMu+"_PRL_supp.png").Data()));
  }

  return;
}
