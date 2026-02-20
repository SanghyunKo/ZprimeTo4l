#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

// #include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
// #include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

#include "/u/user/sako/ModHEEP/CMSSW_10_6_29/work/tdrstyle.C"
#include "/u/user/sako/ModHEEP/CMSSW_10_6_29/work/CMS_lumi.C"

void drawlimit2d_PRL() {
  TString era="run2";
  setTDRStyle();
  // gStyle->SetLineWidth(2);

  writeExtraText = true;       // if extra text
  extraText  = "";  // default extra text is "Preliminary"
  customLumiOffset = 0.06;
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
  int H = 800;

  int H_ref = 800;
  int W_ref = 1600;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.08*H_ref;
  float L = 0.07*W_ref;
  float R = 0.01*W_ref;

  auto* canvas1 = new TCanvas("canvas","canvas",800,800,W,H);
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

  auto retrieve = [] (const TString& hmass, const TString& amass, const TString& isZA, TString quant="") -> double {
    TString filename = TString("./")+isZA+amass+"/"+hmass+"/higgsCombine.X"+hmass+isZA+amass+".HybridNew.mH"+hmass;

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

  auto retrieveSignificance = [] (const TString& hmass, const TString& amass, const TString& isZA) -> double {
    TString filename = TString("./")+isZA+amass+"/"+hmass+"/higgsCombine.significanceCWR.X"+hmass+isZA+amass+".readSignificance.HybridNew.mH"+hmass+".root";

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

  TH2D* customHistA = new TH2D("limitA","limitA",xbin.size()-1,&(xbin[0]),ybinForHist.size()-1,&(ybinForHist[0]));
  TH2D* customHistSigA = new TH2D("significanceA","significanceA",xbin.size()-1,&(xbin[0]),ybinForHist.size()-1,&(ybinForHist[0]));

  TH2D* customHistZ = new TH2D("limitZ","limitZ",xbin.size()-1,&(xbin[0]),ybinForHist.size()-1,&(ybinForHist[0]));
  TH2D* customHistSigZ = new TH2D("significanceZ","significanceZ",xbin.size()-1,&(xbin[0]),ybinForHist.size()-1,&(ybinForHist[0]));

  double maxLimitA = -1., minLimitA = std::numeric_limits<double>::max();
  double maxLimitZ = -1., minLimitZ = std::numeric_limits<double>::max();

  // for HepData Tree
  TTree* limitTreeA = new TTree("limitTreeA","limitTreeA");
  float mX, mY, limit, significance;
  float nominal, quant0p025, quant0p975, quant0p16, quant0p84;
  limitTreeA->Branch("mX",&mX,"mX/F");
  limitTreeA->Branch("mY",&mY,"mY/F");
  limitTreeA->Branch("limit",&limit,"limit/F");
  limitTreeA->Branch("significance",&significance,"significance/F");
  limitTreeA->Branch("nominal",&nominal,"nominal/F");
  limitTreeA->Branch("quant0p025",&quant0p025,"quant0p025/F");
  limitTreeA->Branch("quant0p975",&quant0p975,"quant0p975/F");
  limitTreeA->Branch("quant0p16",&quant0p16,"quant0p16/F");
  limitTreeA->Branch("quant0p84",&quant0p84,"quant0p84/F");
  TTree* limitTreeZ = new TTree("limitTreeZ","limitTreeZ");
  limitTreeZ->Branch("mX",&mX,"mX/F");
  limitTreeZ->Branch("mY",&mY,"mY/F");
  limitTreeZ->Branch("limit",&limit,"limit/F");
  limitTreeZ->Branch("significance",&significance,"significance/F");
  limitTreeZ->Branch("nominal",&nominal,"nominal/F");
  limitTreeZ->Branch("quant0p025",&quant0p025,"quant0p025/F");
  limitTreeZ->Branch("quant0p975",&quant0p975,"quant0p975/F");
  limitTreeZ->Branch("quant0p16",&quant0p16,"quant0p16/F");
  limitTreeZ->Branch("quant0p84",&quant0p84,"quant0p84/F");

  // for the sake of sorting HepData tables in mX
  for (int idx = 0; idx<avecA.size(); idx++) {
    for (const auto hmass : massvec) {
      auto ymassStr = avecA.at(idx);
      double ymass = ybinA.at(idx);

      if (!checkKinematic(static_cast<double>(hmass),ymass,"A"))
        continue;

      double val = retrieve(std::to_string(hmass),ymassStr,"A");
      int ibin = customHistA->FindFixBin(static_cast<float>(hmass)+1.,ymass+0.001);

      double valSignificance = retrieveSignificance(std::to_string(hmass),ymassStr,"A");

      if (valSignificance > 2.)
        std::cout << "A significance > 2 sigma found at Mx=" << hmass << " My=" << ymass << " : " << valSignificance << std::endl;

      if (val > 0.) {
        customHistA->SetBinContent(ibin, val);
        maxLimitA = val>maxLimitA ? val : maxLimitA;
        minLimitA = val<minLimitA ? val : minLimitA;
        mX = static_cast<float>(hmass);
        mY = static_cast<float>(ymass);
        limit = static_cast<float>(val);
        nominal = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,"A","0.500"));
        quant0p025 = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,"A","0.025")) - nominal;
        quant0p975 = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,"A","0.975")) - nominal;
        quant0p16 = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,"A","0.160")) - nominal;
        quant0p84 = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,"A","0.840")) - nominal;

        significance = std::max(valSignificance, 0.);
        customHistSigA->SetBinContent(ibin, significance);
        limitTreeA->Fill();
      }
    }
  }

  for (int idx = 0; idx<avecZ.size(); idx++) {
    for (const auto hmass : massvec) {
      auto ymassStr = avecZ.at(idx);
      double ymass = ybinZ.at(idx);

      if (!checkKinematic(static_cast<double>(hmass),ymass,"Z"))
        continue;

      double val = retrieve(std::to_string(hmass),ymassStr,"Z");
      int ibin = customHistZ->FindFixBin(static_cast<float>(hmass)+1.,ymass+0.001);

      double valSignificance = retrieveSignificance(std::to_string(hmass),ymassStr,"Z");

      if (valSignificance > 2.)
        std::cout << "Z significance > 2 sigma found at Mx=" << hmass << " My=" << ymass << " : " << valSignificance << std::endl;

      if (val > 0.) {
        customHistZ->SetBinContent(ibin, val);
        maxLimitZ = val>maxLimitZ ? val : maxLimitZ;
        minLimitZ = val<minLimitZ ? val : minLimitZ;
        mX = static_cast<float>(hmass);
        mY = static_cast<float>(ymass);
        limit = static_cast<float>(val);
        nominal = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,"Z","0.500"));
        quant0p025 = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,"Z","0.025")) - nominal;
        quant0p975 = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,"Z","0.975")) - nominal;
        quant0p16 = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,"Z","0.160")) - nominal;
        quant0p84 = static_cast<float>(retrieve(std::to_string(hmass),ymassStr,"Z","0.840")) - nominal;

        significance = std::max(valSignificance, 0.);
        customHistSigZ->SetBinContent(ibin, significance);
        limitTreeZ->Fill();
      }
    }
  }

  std::cout << "A limit range: " << minLimitA << " - " << maxLimitA << std::endl;
  std::cout << "Z limit range: " << minLimitZ << " - " << maxLimitZ << std::endl;

  fillY(customHistA,"A");
  fillX(customHistA,"A");

  fillY(customHistSigA,"A");
  fillX(customHistSigA,"A");

  fillY(customHistZ,"Z");
  fillX(customHistZ,"Z");

  fillY(customHistSigZ,"Z");
  fillX(customHistSigZ,"Z");

  std::vector<TH2D*> histList = {customHistA, customHistZ};
  std::vector<TH2D*> histListForHepData = {(TH2D*)customHistA->Clone(), (TH2D*)customHistZ->Clone()};
  std::vector<TH2D*> histSigList = {customHistSigA, customHistSigZ};
  std::vector<TH2D*> histSigList2sig = {};
  std::vector<TString> isZAList = {"A","Z"};

  gStyle->SetPalette(kViridis);
  customHistA->GetXaxis()->SetLabelSize(0.05);
  customHistA->GetYaxis()->SetLabelSize(0.05);
  customHistA->GetZaxis()->SetRangeUser(0.04,12.);
  // customHistA->GetYaxis()->SetTitle("M_{Y} [GeV]");
  // customHistA->GetYaxis()->SetTitleSize(0.06);

  customHistZ->GetXaxis()->SetLabelSize(0.05);
  customHistZ->GetZaxis()->SetLabelSize(0.045);
  customHistZ->GetZaxis()->SetTitleSize(0.05);
  customHistZ->GetZaxis()->SetTitle("\\sigma(\\mathrm{pp}\\rightarrow\\mathrm{X}) \\times \\mathrm{B}(\\mathrm{X}\\rightarrow4\\ell) [\\mathrm{fb}]");
  customHistZ->GetZaxis()->SetTitleOffset(1.25);
  customHistZ->GetZaxis()->SetRangeUser(0.04,12.);

  // multipad
  // Define the rows and columns for subpads
  int rows = 1;
  int cols = 2;

  // Create subpads and divide the canvas
  double leftmargin = L/W;
  double rightmargin = R/W;
  double topmargin = T/H;
  double bottommargin = 1.2*B/H;

  canvas1->SetLogy();
  // canvas1->SetLogx();
  canvas1->SetLogz();
  TPad* padWithPalette = nullptr;
  
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
        xmax -= 0.03;
      if (i==1)
        xmin -= 0.03;

      // Create subpad
      TPad *padUp = new TPad(Form("padUp%d", iPad), Form("SubpadUp %d", iPad), xmin, ymin, xmax, ymax);
      padUp->SetLeftMargin( 0. );
      padUp->SetRightMargin( i==0 ? 0.0 : 0.18 );
      padUp->SetTopMargin( 0.0 );
      padUp->SetBottomMargin( 1.2*bottommargin );
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

      auto limitHist = histList.at(iFig);
      auto significanceHist1sig = histSigList.at(iFig);
      auto isZA = isZAList.at(iFig);

      auto significanceHist2sig = (TH2D*)significanceHist1sig->Clone("significance2sig"+isZA);
      significanceHist1sig->SetContour(1,&(std::vector<double>{1.} )[0]);
      significanceHist1sig->SetLineColor(TColor::GetColor("#f89c20"));
      significanceHist1sig->SetLineWidth(2);
      // significanceHist1sig->SetLineStyle(kDashed);
      significanceHist2sig->SetContour(1,&(std::vector<double>{2.} )[0]);
      significanceHist2sig->SetLineColor(TColor::GetColor("#e42536"));
      significanceHist2sig->SetLineWidth(2);
      histSigList2sig.push_back(significanceHist2sig);

      TString limitDrawOption = (isZA=="A") ? "col" : "colz";
      padWithPalette = padUp; // Store the pad with the palette for later use

      // if (i==0)
      // limitHist->GetXaxis()->ChangeLabel(-1,-1,0.);
      limitHist->GetXaxis()->SetNdivisions(505);

      limitHist->Draw(limitDrawOption);
      significanceHist1sig->Draw("cont3 same");
      significanceHist2sig->Draw("cont3 same");

      double xoffset = (i==0) ? 0.18 : 0.;

      auto* atext = new TPaveText(0.6+xoffset,0.88,0.8+xoffset,0.99,"NDC");
      atext->SetBorderSize(0);
      atext->SetFillColor(0);
      atext->SetFillStyle(0);
      const TString channelName = (isZA=="A") ? "\\mathrm{X}\\rightarrow\\mathrm{YY}" : "\\mathrm{X}\\rightarrow\\mathrm{ZY}";
      atext->AddText(channelName);
      ((TText*)atext->GetListOfLines()->Last())->SetTextColor(kBlack);
      ((TText*)atext->GetListOfLines()->Last())->SetTextAlign(31);
      ((TText*)atext->GetListOfLines()->Last())->SetTextSize(0.06);
      ((TText*)atext->GetListOfLines()->Last())->SetTextFont(42);
      atext->Draw();

      canvas1->cd();  // Go back to the main canvas
    }
  }

  canvas1->Update();  // Refresh canvas

  auto legend = new TLegend(0.088,0.77,0.4,0.92);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->SetMargin(0.12);
  legend->AddEntry(histSigList.front(),"1 s.d. obs. local significance","l");
  legend->AddEntry(histSigList2sig.front(),"2 s.d. obs. local significance","l");
  canvas1->cd();
  legend->Draw();

  auto ytext = std::make_unique<TPaveText>(0.03,0.63,0.06,0.77,"NDC");
  ytext->SetBorderSize(0);
  ytext->SetFillColor(0);
  ytext->SetFillStyle(0);
  TString astr = "M_{Y} [GeV]";
  ytext->AddText(astr);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextFont(42);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextSize(0.055);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextAngle(90);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextAlign(11);
  ytext->Draw();

  auto xtext = std::make_unique<TPaveText>(0.8,0.05,0.92,0.088,"NDC");
  xtext->SetBorderSize(0);
  xtext->SetFillColor(0);
  xtext->SetFillStyle(0);
  TString xstr = "M_{X} [GeV]";
  xtext->AddText(xstr);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextFont(42);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextAlign(33);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextSize(0.055);
  xtext->Draw();

  auto ztext = std::make_unique<TPaveText>(0.9,0.15,0.97,0.6,"NDC");
  ztext->SetBorderSize(0);
  ztext->SetFillColor(0);
  ztext->SetFillStyle(0);
  TString zstr = "Upper limits on";
  ztext->AddText(zstr);
  ((TText*)ztext->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)ztext->GetListOfLines()->Last())->SetTextFont(42);
  ((TText*)ztext->GetListOfLines()->Last())->SetTextAlign(33);
  ((TText*)ztext->GetListOfLines()->Last())->SetTextSize(0.045);
  ((TText*)ztext->GetListOfLines()->Last())->SetTextAngle(90);
  ztext->Draw();

  SaveAs(canvas1,std::string(TString("limit2D_PRL.png").Data()));
  
  // HepData output
  auto outFile = std::make_unique<TFile>("HepData_limit2D_PRL.root","RECREATE");
  limitTreeA->Write();
  limitTreeZ->Write();
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

  for (unsigned iFig=0; iFig<histListForHepData.size(); iFig++) {
    auto limitHist = histListForHepData.at(iFig);
    auto significanceHist1sig = histSigList.at(iFig);
    auto significanceHist2sig = histSigList2sig.at(iFig);
    auto isZA = isZAList.at(iFig);

    // if (i==0)
    // limitHist->GetXaxis()->ChangeLabel(-1,-1,0.);
    limitHist->GetXaxis()->SetNdivisions(505);

    limitHist->GetXaxis()->SetLabelSize(0.035);
    limitHist->GetXaxis()->SetTitle("M_{X} [GeV]");
    limitHist->GetXaxis()->SetTitleSize(0.04);
    limitHist->GetYaxis()->SetLabelSize(0.04);
    limitHist->GetYaxis()->SetTitle("M_{Y} [GeV]");
    limitHist->GetYaxis()->SetTitleSize(0.04);
    limitHist->GetZaxis()->SetLabelSize(0.032);
    limitHist->GetZaxis()->SetTitleSize(0.04);
    limitHist->GetZaxis()->SetTitle("\\sigma(\\mathrm{pp}\\rightarrow\\mathrm{X}) \\times \\mathrm{B}(\\mathrm{X}\\rightarrow4\\ell) [\\mathrm{fb}]");
    limitHist->GetZaxis()->SetTitleOffset(1.25);
    limitHist->GetZaxis()->SetRangeUser(0.04,12.);

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
    legend2->AddEntry(histSigList.front(),"1 s.d. obs. local significance","l");
    legend2->AddEntry(histSigList2sig.front(),"2 s.d. obs. local significance","l");
    legend2->Draw();

    auto* atext = new TPaveText(0.73,0.82,0.8,0.88,"NDC");
    atext->SetBorderSize(0);
    atext->SetFillColor(0);
    atext->SetFillStyle(0);
    const TString channelName = (isZA=="A") ? "\\mathrm{X}\\rightarrow\\mathrm{YY}" : "\\mathrm{X}\\rightarrow\\mathrm{ZY}";
    atext->AddText(channelName);
    ((TText*)atext->GetListOfLines()->Last())->SetTextColor(kBlack);
    ((TText*)atext->GetListOfLines()->Last())->SetTextAlign(31);
    ((TText*)atext->GetListOfLines()->Last())->SetTextSize(0.05);
    ((TText*)atext->GetListOfLines()->Last())->SetTextFont(42);
    atext->Draw();

    auto btext = std::make_unique<TPaveText>(0.9,0.28,0.956,0.6,"NDC");
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
    SaveAs(canvas2,std::string((TString("HepData_limit2D_")+isZA+TString("_PRL.png")).Data()));
  }

  return;
}
