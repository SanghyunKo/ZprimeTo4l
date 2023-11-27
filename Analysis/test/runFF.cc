#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void runFF(TString era) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "   Work in progress";  // default extra text is "Preliminary"

  static double valLumi = 0.;
  static constexpr double WZxsec_ = 0.65*62.78;
  static constexpr double ZZxsec_ = 13.81;
  static TString postfix = era;

  if (era=="20UL16APV") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "19.5 fb^{-1}";
    valLumi = 19.5;
  } else if (era=="20UL16") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "16.8 fb^{-1}";
    valLumi = 16.8;
    postfix = "";
  } else if (era=="20UL17") {
    lumi_sqrtS = "2017 (13 TeV)";
    lumi_13TeV = "41.48 fb^{-1}";
    valLumi = 41.48;
  } else if (era=="20UL18") {
    lumi_sqrtS = "2018 (13 TeV)";
    lumi_13TeV = "59.83 fb^{-1}";
    valLumi = 59.83;
  } else {
    std::cout << "check era..." << std::endl;
  }

  static TString anlyzrMC = "mergedEleCRanalyzer"+postfix;
  static TString anlyzrData = "mergedEleCRanalyzerData";

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

  TFile* datafile = new TFile("MergedEleCR_"+era+"_data.root","READ");
  TFile* WZfile = new TFile("MergedEleCR_"+era+"_WZ.root","READ");
  TFile* ZZfile = new TFile("MergedEleCR_"+era+"_ZZ.root","READ");

  //TFile* H250A1file = new TFile("MergedEleCR_"+era+"_H250A1.root","READ");
  //TFile* H750A1file = new TFile("MergedEleCR_"+era+"_H750A1.root","READ");
  //TFile* H2000A1file = new TFile("MergedEleCR_"+era+"_H2000A1.root","READ");
  //TFile* H250A10file = new TFile("MergedEleCR_"+era+"_H250A10.root","READ");
  //TFile* H750A10file = new TFile("MergedEleCR_"+era+"_H750A10.root","READ");
  //TFile* H2000A10file = new TFile("MergedEleCR_"+era+"_H2000A10.root","READ");

  class SigSample {
  public:
    SigSample(TFile* afile, TString name)
    : file_(afile), name_(name) {}

    ~SigSample()=default;

  public:
    TFile* file_;
    TString name_;
  };

  auto H750A1sample = SigSample(new TFile("MergedEleCR_"+era+"_H750A1.root","READ"),"H750A1");

  std::vector<SigSample> sigsamples = {H750A1sample};
  std::vector<SigSample> sigsamples3E = {H750A1sample};

  auto* canvas_2 = new TCanvas("canvas_2","canvas_2",50,50,W,H);
  canvas_2->SetFillColor(0);
  canvas_2->SetBorderMode(0);
  canvas_2->SetFrameFillStyle(0);
  canvas_2->SetFrameBorderMode(0);
  canvas_2->SetLeftMargin( L/W );
  canvas_2->SetRightMargin( R/W );
  canvas_2->SetTopMargin( T/H );
  canvas_2->SetBottomMargin( B/H );
  canvas_2->SetTickx(0);
  canvas_2->SetTicky(0);

  TPad* p1 = new TPad("p1","",0,0.3,1,0.95);
  p1->SetFillColor(0);
  p1->SetFrameBorderSize(0);
  p1->SetBorderMode(0);
  p1->SetFrameFillStyle(0);
  p1->SetFrameBorderMode(0);
  p1->SetTickx(0);
  p1->SetTicky(0);
  p1->SetBottomMargin(0.00);
  p1->SetLeftMargin( L/W );
  p1->SetRightMargin( R/W );
  p1->Draw();

  TPad* p2 = new TPad("p2","",0,0,1,0.3);
  p2->SetTickx(0);
  p2->SetTicky(0);
  p2->SetTopMargin(0.0);
  p2->SetBottomMargin(0.3);
  p2->SetLeftMargin( L/W );
  p2->SetRightMargin( R/W );
  p2->Draw();

  auto* canvas_1 = new TCanvas("canvas_1","canvas_1",50,50,W,H);
  canvas_1->SetFillColor(0);
  canvas_1->SetBorderMode(0);
  canvas_1->SetFrameFillStyle(0);
  canvas_1->SetFrameBorderMode(0);
  canvas_1->SetLeftMargin( L/W );
  canvas_1->SetRightMargin( R/W );
  canvas_1->SetTopMargin( T/H );
  canvas_1->SetBottomMargin( B/H );
  canvas_1->SetTickx(0);
  canvas_1->SetTicky(0);

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

  class HistLoaderBase {
  protected:
    TFile* datafile_ = nullptr;
    TFile* WZfile_ = nullptr;
    TFile* ZZfile_ = nullptr;

    TFile* datacard_ = nullptr;
    TDirectory* dir_ = nullptr;

  public:
    HistLoaderBase(TFile* adatafile, TFile* aWZfile, TFile* aZZfile) {
      datafile_ = adatafile;
      WZfile_ = aWZfile;
      ZZfile_ = aZZfile;
    }

    HistLoaderBase(TFile* adatafile) {
      datafile_ = adatafile;
    }

    ~HistLoaderBase()=default;

    void preparecard(TString name, TString dirname) {
      datacard_ = new TFile(name,"RECREATE");
      dir_ = datacard_->mkdir(dirname);
      dir_->cd();
    }

    void close() {
      datacard_->Close();
      datacard_ = nullptr;
      dir_ = nullptr;
    }
  };

  class HistLoader2E : public HistLoaderBase {
  public:
    HistLoader2E(TFile* adatafile, std::vector<SigSample> sigFiles)
    : HistLoaderBase(adatafile),
      sigFiles_(sigFiles) {}

    ~HistLoader2E() {
      delete dataHist_, FFHist_;

      if (FFHistUp_) {
        delete FFHistUp_, FFHistDn_;
        FFHistUp_ = nullptr;
      }

      if (dataHist2_) {
        delete dataHist2_, FFHist2_;
        dataHist2_ = nullptr;

        if (FFHist2Up_) {
          delete FFHist2Up_, FFHist2Dn_;
          FFHist2Up_ = nullptr;
        }
      }

    }

  private:
    TH1D* dataHist_ = nullptr;
    TH1D* FFHist_ = nullptr;

    TH1D* dataHist2_ = nullptr;
    TH1D* FFHist2_ = nullptr;

    TH1D* FFHistUp_ = nullptr;
    TH1D* FFHistDn_ = nullptr;
    TH1D* FFHist2Up_ = nullptr;
    TH1D* FFHist2Dn_ = nullptr;

    std::vector<SigSample> sigFiles_;
    std::vector<TH1D*> sigHist_;
    std::vector<TH1D*> sigHist_modHeepIdUp_;
    std::vector<TH1D*> sigHist_modHeepIdDn_;
    std::vector<TH1D*> sigHist_mergedEleIdUp_;
    std::vector<TH1D*> sigHist_mergedEleIdDn_;

  public:
    void load(const std::string& nameNum, const std::string& name, int color, const std::string& nameNum2="", const std::string& name2="", int color2=0) {
      if (dataHist_)
        delete dataHist_, FFHist_;

      if (FFHistUp_) {
        delete FFHistUp_, FFHistDn_;
        FFHistUp_ = nullptr;
      }

      if (dataHist2_) {
        delete dataHist2_, FFHist2_;
        dataHist2_ = nullptr;

        if (FFHist2Up_) {
          delete FFHist2Up_, FFHist2Dn_;
          FFHist2Up_ = nullptr;
        }
      }

      dataHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+nameNum).c_str() )->Clone();
      FFHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name).c_str() )->Clone();

      FFHist_->SetLineWidth(0);
      FFHist_->SetFillColor(color);

      if ( TString(name).Contains("invM") || TString(name).Contains("eta") ) {
        FFHistUp_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_up").c_str() )->Clone();
        FFHistDn_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_dn").c_str() )->Clone();
      }

      const double lumi = valLumi;
      sigHist_.clear();
      sigHist_modHeepIdUp_.clear();
      sigHist_modHeepIdDn_.clear();
      sigHist_mergedEleIdUp_.clear();
      sigHist_mergedEleIdDn_.clear();


      if (nameNum2!="") {
        dataHist2_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+nameNum2).c_str() )->Clone();
        FFHist2_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name2).c_str() )->Clone();
        FFHist2_->SetLineWidth(0);
        FFHist2_->SetFillColor(color2);

        if ( TString(name2).Contains("invM") ) {
          FFHist2Up_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name2+"_up").c_str() )->Clone();
          FFHist2Dn_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name2+"_dn").c_str() )->Clone();
        }

        auto retrieveSigHist = [this,&nameNum,&nameNum2,&lumi] (TFile* afile, const std::string& systName) -> TH1D* {
          TH1D* ahist = (TH1D*)afile->Get( (std::string(anlyzrMC+"/")+nameNum+systName).c_str() )->Clone();
          ahist->Add( (TH1D*)afile->Get( (std::string(anlyzrMC+"/")+nameNum2+systName).c_str() ) );
          ahist->Scale( lumi*1000.*0.001 / ( (TH1D*)afile->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );

          return ahist;
        };

        for (unsigned idx=0; idx<sigFiles_.size(); idx++) {
          TH1D* asigHist = retrieveSigHist(sigFiles_.at(idx).file_,"");
          sigHist_.push_back( asigHist );
          sigHist_.back()->SetLineWidth(2);
          int color = (idx > 2) ? kOrange : kRed;
          sigHist_.back()->SetLineColor(color);

          if ( TString(name2).Contains("invM") ) {
            TH1D* sigModHeepUp = retrieveSigHist(sigFiles_.at(idx).file_,"_heepIdUp");
            sigHist_modHeepIdUp_.push_back( sigModHeepUp );
            TH1D* sigModHeepDn = retrieveSigHist(sigFiles_.at(idx).file_,"_heepIdDn");
            sigHist_modHeepIdDn_.push_back( sigModHeepDn );
            TH1D* sigMergedEleUp = retrieveSigHist(sigFiles_.at(idx).file_,"_mergedEleIdUp");
            sigHist_mergedEleIdUp_.push_back( sigMergedEleUp );
            TH1D* sigMergedEleDn = retrieveSigHist(sigFiles_.at(idx).file_,"_mergedEleIdDn");
            sigHist_mergedEleIdDn_.push_back( sigMergedEleDn );
          }
        }
      }
    }

    void compare(TPad* padUp, bool removeZero=false, TPad* padDn=nullptr, int rebin=1) {
      if ( FFHist_->GetNbinsX() % dataHist_->GetNbinsX()!=0 ) {
        // hardcode (rebin to GCD)
        FFHist_->Rebin( 5 );
        dataHist_->Rebin( 2 );

        if (FFHistUp_) {
          FFHistUp_->Rebin( 5 );
          FFHistDn_->Rebin( 5 );
        }

        if (dataHist2_) {
          FFHist2_->Rebin( 5 );
          dataHist2_->Rebin( 2 );

          if (FFHist2Up_) {
            FFHist2Up_->Rebin( 5 );
            FFHist2Dn_->Rebin( 5 );
          }
        }
      } else {
        FFHist_->Rebin( FFHist_->GetNbinsX()/dataHist_->GetNbinsX() );

        if (FFHistUp_) {
          FFHistUp_->Rebin( FFHistUp_->GetNbinsX()/dataHist_->GetNbinsX() );
          FFHistDn_->Rebin( FFHistDn_->GetNbinsX()/dataHist_->GetNbinsX() );
        }

        if (dataHist2_) {
          FFHist2_->Rebin( FFHist2_->GetNbinsX()/dataHist2_->GetNbinsX() );

          if (FFHist2Up_) {
            FFHist2Up_->Rebin( FFHist2Up_->GetNbinsX()/dataHist2_->GetNbinsX() );
            FFHist2Dn_->Rebin( FFHist2Dn_->GetNbinsX()/dataHist2_->GetNbinsX() );
          }
        }
      }

      dataHist_->Rebin(rebin);
      FFHist_->Rebin(rebin);

      if (FFHistUp_) {
        FFHistUp_->Rebin(rebin);
        FFHistDn_->Rebin(rebin);
      }

      if (dataHist2_) {
        dataHist2_->Rebin(rebin);
        FFHist2_->Rebin(rebin);

        if (FFHist2Up_) {
          FFHist2Up_->Rebin(rebin);
          FFHist2Dn_->Rebin(rebin);
        }
      }

      TH1D* dataHist = dataHist_;
      TH1D* FFadded = (TH1D*)FFHist_->Clone();
      TNamed* FFobj = FFHist_;

      if (dataHist2_) {
        dataHist->Add(dataHist2_);
        FFadded->Add(FFHist2_);
        auto* stack = new THStack("stack",";GeV;");
        stack->Add(FFHist_);
        stack->Add(FFHist2_);
        FFobj = stack;
      }

      if ( TString(dataHist_->GetName()).Contains("invM") && TString(dataHist_->GetName()).Contains("mixedME") )
        dataHist->GetXaxis()->SetRangeUser(0.,500.);
      else if ( TString(dataHist_->GetName()).Contains("invM") && TString(dataHist_->GetName()).Contains("CRME") ) {
        if (dataHist2_)
          dataHist->GetXaxis()->SetRangeUser(0.,1000.);
        else
          dataHist->GetXaxis()->SetRangeUser(0.,250.);
      }

      padUp->cd();

      if ( padUp->GetLogy() )
        dataHist->SetMaximum( 10.*std::max(dataHist_->GetMaximum(),FFadded->GetMaximum()) );
      else
        dataHist->SetMaximum( 1.2*std::max(dataHist_->GetMaximum(),FFadded->GetMaximum()) );

      if (removeZero)
        dataHist->SetMinimum( 0.001 );

      dataHist->SetLineWidth(2);
      dataHist->SetLineColor(kBlack);
      dataHist->Draw("E1");
      FFobj->Draw("hist&same");
      dataHist->Draw("E1&same");

      if (dataHist2_) {
        TLegend* legend = new TLegend(0.45,0.72,0.95,0.93);
        legend->SetBorderSize(0);
        legend->SetNColumns(2);
        legend->AddEntry(dataHist,"Data");

        legend->AddEntry(FFHist_,"Data-driven (X #rightarrow ME)");
        legend->AddEntry(FFHist2_,"Data-driven (e #rightarrow ME)");

        if ( TString(dataHist_->GetName()).Contains("CRME") && TString(dataHist_->GetName()).Contains("invM") ) {
          for (unsigned idx=0; idx<sigHist_.size(); idx++) {
            sigHist_.at(idx)->Rebin(rebin);
            sigHist_modHeepIdUp_.at(idx)->Rebin(rebin);
            sigHist_modHeepIdDn_.at(idx)->Rebin(rebin);
            sigHist_mergedEleIdUp_.at(idx)->Rebin(rebin);
            sigHist_mergedEleIdDn_.at(idx)->Rebin(rebin);
            sigHist_.at(idx)->Draw("hist&same");
          }

          //legend->AddEntry(sigHist_.at(3),"H250/750/2000A1");
          legend->AddEntry(sigHist_.at(0),"M_{A} = 1 GeV (#sigma = 1 fb)");
        }

        legend->Draw();
      }

      TH1D* ratio = nullptr;

      if (padDn) {
        ratio = (TH1D*)dataHist->Clone();
        ratio->SetStats(0);
        ratio->SetTitle("");
        ratio->Divide(FFadded);
        ratio->GetYaxis()->SetTitle("Obs/Exp");
        ratio->GetYaxis()->SetTitleSize(0.1);
        ratio->GetYaxis()->SetTitleOffset(0.4);
        ratio->GetXaxis()->SetLabelSize(0.1);
        ratio->GetYaxis()->SetLabelSize(0.1);
        ratio->GetXaxis()->SetLabelOffset(0.01);
        ratio->GetYaxis()->SetLabelOffset(0.01);
        ratio->GetYaxis()->SetRangeUser(0.2,1.8);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetTitleOffset(0.75);
        ratio->SetLineColor(kBlack);

        padDn->cd();
        ratio->Draw("E1");
      }

      if (FFHistUp_) {
        std::vector<double> x0, y0, errx, erryDn, erryUp;
        std::vector<double> r0, errRdn, errRup;

        TH1D* FFHistUpAdded = (TH1D*)FFHistUp_->Clone();
        TH1D* FFHistDnAdded = (TH1D*)FFHistDn_->Clone();
        TH1D* FFHist2UpAdded = nullptr;
        TH1D* FFHist2DnAdded = nullptr;

        if (FFHist2Up_) {
          FFHistUpAdded->Add(FFHist2_);
          FFHistDnAdded->Add(FFHist2_);
          FFHist2UpAdded = (TH1D*)FFHist2Up_->Clone();
          FFHist2DnAdded = (TH1D*)FFHist2Dn_->Clone();
          FFHist2UpAdded->Add(FFHist_);
          FFHist2DnAdded->Add(FFHist_);
        }

        for (unsigned idx = 1; idx <= FFadded->GetNbinsX(); idx++) {
          x0.push_back(FFadded->GetBinCenter(idx));
          y0.push_back(FFadded->GetBinContent(idx));
          errx.push_back(FFadded->GetBinWidth(idx)/2.);

          double valFFUp = FFHistUp_->GetBinContent(idx) - FFHist_->GetBinContent(idx);
          double valFFDn = FFHist_->GetBinContent(idx) - FFHistDn_->GetBinContent(idx);
          double valFF2Up = 0.;
          double valFF2Dn = 0.;

          if (FFHist2Up_) {
            valFF2Up = FFHist2Up_->GetBinContent(idx) - FFHist2_->GetBinContent(idx);
            valFF2Dn = FFHist2_->GetBinContent(idx) - FFHist2Dn_->GetBinContent(idx);
          }

          erryUp.push_back( std::sqrt( valFFUp*valFFUp + valFF2Up*valFF2Up ) );
          erryDn.push_back( std::sqrt( valFFDn*valFFDn + valFF2Dn*valFF2Dn ) );

          if (ratio) {
            r0.push_back(1.);

            double rvalFFUp = FFadded->GetBinContent(idx)/FFHistDnAdded->GetBinContent(idx) - 1.;
            double rvalFFDn = 1. - FFadded->GetBinContent(idx)/FFHistUpAdded->GetBinContent(idx);
            double rvalFF2Up = 0.;
            double rvalFF2Dn = 0.;

            if (FFHist2Up_) {
              rvalFF2Up = FFadded->GetBinContent(idx)/FFHist2DnAdded->GetBinContent(idx) - 1.;
              rvalFF2Dn = 1. - FFadded->GetBinContent(idx)/FFHist2UpAdded->GetBinContent(idx);
            }

            errRup.push_back( std::sqrt( rvalFFUp*rvalFFUp + rvalFF2Up*rvalFF2Up ) );
            errRdn.push_back( std::sqrt( rvalFFDn*rvalFFDn + rvalFF2Dn*rvalFF2Dn ) );
          }
        }

        auto gr = new TGraphAsymmErrors(dataHist->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
        gr->SetFillColor(kGray+2);
        gr->SetLineColor(kGray+2);
        gr->SetFillStyle(3004);
        padUp->cd();
        gr->Draw("2");

        if (padDn) {
          auto rgr = new TGraphAsymmErrors(dataHist->GetNbinsX(),&(x0[0]),&(r0[0]),&(errx[0]),&(errx[0]),&(errRdn[0]),&(errRup[0]));
          rgr->SetFillColor(kGray+2);
          rgr->SetLineColor(kGray+2);
          rgr->SetFillStyle(3004);
          padDn->cd();
          rgr->Draw("2");
        }

        if (dir_) {
          dir_->WriteTObject(dataHist_,"data_obs");
          dir_->WriteTObject(FFHist_,"SS");
          dir_->WriteTObject(FFHist2_,"OS");

          dir_->WriteTObject(FFHistUp_,"SS_mergedEleFakeFactorSSUp");
          dir_->WriteTObject(FFHistDn_,"SS_mergedEleFakeFactorSSDown");
          dir_->WriteTObject(FFHist2Up_,"OS_mergedEleFakeFactorOSUp");
          dir_->WriteTObject(FFHist2Dn_,"OS_mergedEleFakeFactorOSDown");

          for (unsigned idx=0; idx<sigHist_.size(); idx++) {
            dir_->WriteTObject(sigHist_.at(idx),sigFiles_.at(idx).name_);
            dir_->WriteTObject(sigHist_modHeepIdUp_.at(idx),sigFiles_.at(idx).name_+"_modHeepIdUp");
            dir_->WriteTObject(sigHist_modHeepIdDn_.at(idx),sigFiles_.at(idx).name_+"_modHeepIdDown");
            dir_->WriteTObject(sigHist_mergedEleIdUp_.at(idx),sigFiles_.at(idx).name_+"_mergedEleIdUp");
            dir_->WriteTObject(sigHist_mergedEleIdDn_.at(idx),sigFiles_.at(idx).name_+"_mergedEleIdDown");
          }
        }
      }
    }
  };

  // invM
  auto aloader2E = HistLoader2E(datafile,sigsamples);

  aloader2E.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1);
  aloader2E.compare(p1,true,p2,2);
  SaveAs(canvas_2,"FF_2E_invM_denom_antiSS.png",p1);

  aloader2E.load("2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange);
  aloader2E.compare(p1,true,p2);
  SaveAs(canvas_2,"FF_2E_invM_denom_mixedOS.png",p1);

  p1->SetLogy();
  aloader2E.load("2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange);
  aloader2E.compare(p1,false,p2);
  SaveAs(canvas_2,"FF_2E_invM_denom_antiOS.png",p1);

  //p1->SetLogx();
  //p2->SetLogx();
  p1->SetLogy();
  aloader2E.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,"2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange);
  aloader2E.preparecard("ME2E_"+era+"_datacard.root","mergedEle2E");
  aloader2E.compare(p1,false,p2,2);
  SaveAs(canvas_2,"FF_2E_invM_denom_CR.png",p1);
  aloader2E.close();

  p1->SetLogy();
  aloader2E.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,"2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange);
  aloader2E.compare(p1,false,p2,2);
  SaveAs(canvas_2,"FF_2E_invM_mixed_CR.png",p1);

  //p1->SetLogx(0);
  //p2->SetLogx(0);

  aloader2E.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1);
  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_invM_denom_mixedSS.png");

  // Et
  aloader2E.load("2E_Et_SSCR_EB_CRME","2E_Et_SSCR_EB_mixedAntiME_xFF",kCyan+1);
  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_Et_denom_mixedSS.png");

  aloader2E.load("2E_Et_SSCR_EB_mixedME","2E_Et_SSCR_EB_antiME_xFF",kCyan+1);
  aloader2E.compare(p1,false,p2);
  SaveAs(canvas_2,"FF_2E_Et_denom_antiSS.png",p1);

  aloader2E.load("2E_Et_OSCR_EB_CRME","2E_Et_OSCR_EB_mixedAntiME_xFF",kOrange);
  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_Et_denom_mixedOS.png");

  aloader2E.load("2E_Et_OSCR_EB_mixedME","2E_Et_OSCR_EB_antiME_xFF",kOrange);
  aloader2E.compare(p1,false,p2);
  SaveAs(canvas_2,"FF_2E_Et_denom_antiOS.png",p1);

  // eta
  p1->SetLogy(0);
  aloader2E.load("2E_CRME_SSll_eta","2E_mixedAntiME_SSll_eta_xFF",kCyan+1);
  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_eta_denom_mixedSS.png");

  aloader2E.load("2E_mixedME_SSll_eta","2E_antiME_SSll_eta_xFF",kCyan+1);
  aloader2E.compare(p1,true,p2);
  SaveAs(canvas_2,"FF_2E_eta_denom_antiSS.png",p1);

  aloader2E.load("2E_CRME_OSll_eta","2E_mixedAntiME_OSll_eta_xFF",kOrange);
  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_eta_denom_mixedOS.png");

  aloader2E.load("2E_mixedME_OSll_eta","2E_antiME_OSll_eta_xFF",kOrange);
  aloader2E.compare(p1,true,p2);
  SaveAs(canvas_2,"FF_2E_eta_denom_antiOS.png",p1);

  class HistLoader3E : public HistLoaderBase {
  public:
    HistLoader3E(TFile* adatafile, TFile* aWZfile, TFile* aZZfile, std::vector<SigSample> sigFiles)
    : HistLoaderBase(adatafile,aWZfile,aZZfile),
      sigFiles_(sigFiles) {}

    ~HistLoader3E() {
      delete dataHist_;
    }

    class denomHists {
    public:
      denomHists() {}
      ~denomHists()=default;

    private:
      TH1D* SSFFHist_ = nullptr;
      TH1D* SSWZHist_ = nullptr;
      TH1D* SSZZHist_ = nullptr;
      TH1D* OSWZHist_ = nullptr;
      TH1D* OSZZHist_ = nullptr;

    public:
      void load(TString denomNameSS, TString denomNameOS, TFile* datafile, TFile* WZfile, TFile* ZZfile) {
        SSFFHist_ = (TH1D*)datafile->Get(anlyzrData+"/"+denomNameSS)->Clone();

        const double WZsumwgt = ((TH1D*)WZfile->Get("evtCounter/h_sumW"))->GetBinContent(1);
        const double ZZsumwgt = ((TH1D*)ZZfile->Get("evtCounter/h_sumW"))->GetBinContent(1);

        SSWZHist_ = (TH1D*)WZfile->Get(anlyzrMC+"/"+denomNameSS)->Clone();
        SSWZHist_->Scale( valLumi*1000.*WZxsec_/WZsumwgt );
        SSZZHist_ = (TH1D*)ZZfile->Get(anlyzrMC+"/"+denomNameSS)->Clone();
        SSZZHist_->Scale( valLumi*1000.*ZZxsec_/ZZsumwgt );
        OSWZHist_ = (TH1D*)WZfile->Get(anlyzrMC+"/"+denomNameOS)->Clone();
        OSWZHist_->Scale( valLumi*1000.*WZxsec_/WZsumwgt );
        OSZZHist_ = (TH1D*)ZZfile->Get(anlyzrMC+"/"+denomNameOS)->Clone();
        OSZZHist_->Scale( valLumi*1000.*ZZxsec_/ZZsumwgt );
      }

      void rebin(int r) {
        SSFFHist_->Rebin(r);
        SSWZHist_->Rebin(r);
        SSZZHist_->Rebin(r);
        OSWZHist_->Rebin(r);
        OSZZHist_->Rebin(r);
      }

      TH1D* SSFFHist() { return SSFFHist_; }

      std::unique_ptr<TH1D> subtractHist(const TH1D* denom, const TH1D* denom_prompt) {
        auto cloned = std::unique_ptr<TH1D>((TH1D*)denom->Clone());

        for (unsigned idx = 0; idx < cloned->GetNbinsX()+2; idx++) {
          cloned->SetBinContent(idx,std::max(cloned->GetBinContent(idx)-denom_prompt->GetBinContent(idx),0.));
          cloned->SetBinError(idx,std::hypot(cloned->GetBinError(idx),denom_prompt->GetBinError(idx)));
        }

        return std::move(cloned);
      };

      std::unique_ptr<TH1D> returnSS() {
        auto tempSubtract = subtractHist( SSFFHist_, SSWZHist_ );
        auto denomSSfinal = subtractHist( tempSubtract.get(), SSZZHist_ );
        denomSSfinal->SetFillColor(kCyan+1);

        return std::move(denomSSfinal);
      }

      std::unique_ptr<TH1D> returnOS() {
        auto denomOSfinal = std::unique_ptr<TH1D>((TH1D*)OSWZHist_->Clone());
        denomOSfinal->Add(OSZZHist_);
        denomOSfinal->SetFillColor(kOrange);

        return std::move(denomOSfinal);
      }

      std::unique_ptr<TH1D> returnAdded() {
        auto ss = returnSS();
        auto os = returnOS();
        ss->Add( os.get() );

        return std::move(ss);
      }
    }; // class denomHists

  private:
    TH1D* dataHist_ = nullptr;
    denomHists nominal_;
    denomHists SSFFup_;
    denomHists SSFFdn_;
    denomHists OSFFup_;
    denomHists OSFFdn_;
    denomHists heepIdUp_;
    denomHists heepIdDn_;

    std::vector<SigSample> sigFiles_;
    std::vector<TH1D*> sigHist_;
    std::vector<TH1D*> sigHist_modHeepIdUp_;
    std::vector<TH1D*> sigHist_modHeepIdDn_;
    std::vector<TH1D*> sigHist_mergedEleIdUp_;
    std::vector<TH1D*> sigHist_mergedEleIdDn_;

  public:
    void load(TString numName, TString denomName) {
      if (dataHist_) {
        delete dataHist_;
      }

      dataHist_ = (TH1D*)datafile_->Get(anlyzrData+"/"+numName)->Clone();
      nominal_.load(denomName+"_xSSFF",denomName+"_xOSFF",datafile_,WZfile_,ZZfile_);
      SSFFup_.load(denomName+"_xSSFF_up",denomName+"_xOSFF",datafile_,WZfile_,ZZfile_);
      SSFFdn_.load(denomName+"_xSSFF_dn",denomName+"_xOSFF",datafile_,WZfile_,ZZfile_);
      OSFFup_.load(denomName+"_xSSFF",denomName+"_xOSFF_up",datafile_,WZfile_,ZZfile_);
      OSFFdn_.load(denomName+"_xSSFF",denomName+"_xOSFF_dn",datafile_,WZfile_,ZZfile_);

      sigHist_.clear();
      sigHist_modHeepIdUp_.clear();
      sigHist_modHeepIdDn_.clear();
      sigHist_mergedEleIdUp_.clear();
      sigHist_mergedEleIdDn_.clear();

      if (numName.Contains("invM")) {
        const unsigned isigDiv = 3;

        auto retrieveSigHist = [this,&numName] (TFile* afile, const TString& systName) -> TH1D* {
          TH1D* ahist = (TH1D*)afile->Get( anlyzrMC+"/"+numName+systName )->Clone();
          ahist->Scale( valLumi*1000.*0.001 / ( (TH1D*)afile->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );

          return ahist;
        };

        for (unsigned isig = 0; isig < sigFiles_.size(); isig++) {
          sigHist_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"") );
          sigHist_.back()->SetLineWidth(2);
          sigHist_.back()->SetLineColor(kRed-3);

          sigHist_modHeepIdUp_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"_heepIdUp") );
          sigHist_modHeepIdDn_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"_heepIdDn") );
          sigHist_mergedEleIdUp_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"_mergedEleIdUp") );
          sigHist_mergedEleIdDn_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"_mergedEleIdDn") );

          if (isig > isigDiv-1)
            sigHist_.back()->SetLineColor(kBlue-3);
        }

        heepIdUp_.load(denomName+"_xSSFF_heepIdUp",denomName+"_xOSFF_heepIdUp",datafile_,WZfile_,ZZfile_);
        heepIdDn_.load(denomName+"_xSSFF_heepIdDn",denomName+"_xOSFF_heepIdDn",datafile_,WZfile_,ZZfile_);
      }
    }

    void compare(TPad* pad, int rebin=1) {
      if (rebin!=1) {
        dataHist_->Rebin(rebin);

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx)->Rebin(rebin);
          sigHist_modHeepIdUp_.at(idx)->Rebin(rebin);
          sigHist_modHeepIdDn_.at(idx)->Rebin(rebin);
          sigHist_mergedEleIdUp_.at(idx)->Rebin(rebin);
          sigHist_mergedEleIdDn_.at(idx)->Rebin(rebin);
        }
      }

      nominal_.rebin(nominal_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      SSFFup_.rebin(SSFFup_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      SSFFdn_.rebin(SSFFdn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      OSFFup_.rebin(OSFFup_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      OSFFdn_.rebin(OSFFdn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());

      auto denomSSfinal = nominal_.returnSS();
      auto denomOSfinal = nominal_.returnOS();

      auto denomfinal_nominal = nominal_.returnAdded();
      auto denomfinal_SSFFup = SSFFup_.returnAdded();
      auto denomfinal_SSFFdn = SSFFdn_.returnAdded();
      auto denomfinal_OSFFup = OSFFup_.returnAdded();
      auto denomfinal_OSFFdn = OSFFdn_.returnAdded();
      std::unique_ptr<TH1D> denomfinal_heepIdUp, denomfinal_heepIdDn;

      if ( TString(dataHist_->GetName()).Contains("invM") ) {
        heepIdUp_.rebin(heepIdUp_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        heepIdDn_.rebin(heepIdDn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        denomfinal_heepIdUp = heepIdUp_.returnAdded();
        denomfinal_heepIdDn = heepIdDn_.returnAdded();
      }

      THStack* denomFinal = new THStack("final",";GeV;");
      denomOSfinal->SetLineWidth(0);
      denomSSfinal->SetLineWidth(0);
      denomFinal->Add(denomOSfinal.get());
      denomFinal->Add(denomSSfinal.get());

      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->SetMaximum(1.5*dataHist_->GetMaximum());

      if (TString(dataHist_->GetName()).Contains("invM")) {
        dataHist_->GetXaxis()->SetRangeUser(0.,1000.);
        dataHist_->SetMaximum(5.*dataHist_->GetMaximum());
        //dataHist_->SetMinimum(0.2);
      }

      dataHist_->Draw("E1");
      denomFinal->Draw("hist&same");
      dataHist_->Draw("E1&same");

      if (TString(dataHist_->GetName()).Contains("invM")) {
        for (auto* sigNum : sigHist_)
          sigNum->Draw("hist&same");
      }

      if (!TString(dataHist_->GetName()).Contains("Et")) {
        TLegend* legend = new TLegend(0.55,0.7,0.95,0.9);
        legend->SetBorderSize(0);
        legend->AddEntry(dataHist_,"Data");
        legend->AddEntry(denomSSfinal.get(),"Data-driven (X #rightarrow ME)");
        legend->AddEntry(denomOSfinal.get(),"Data-driven (e #rightarrow ME)");

        if (TString(dataHist_->GetName()).Contains("invM")) {
          legend->AddEntry(sigHist_.at(0),"M_{A} = 1 GeV (#sigma = 1 fb)");
          //legend_left->AddEntry(sigNums.at(isigDiv),"M_{A} = 10 GeV (#sigma = 10 fb)");
        }

        legend->Draw();
      }

      std::vector<double> x0, y0, errx, erryDn, erryUp;

      for (unsigned idx = 1; idx <= denomfinal_nominal->GetNbinsX(); idx++) {
        x0.push_back(denomfinal_nominal->GetBinCenter(idx));
        y0.push_back(denomfinal_nominal->GetBinContent(idx));
        errx.push_back(denomfinal_nominal->GetBinWidth(idx)/2.);

        double valSSFFup = denomfinal_SSFFup->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
        double valSSFFdn = denomfinal_nominal->GetBinContent(idx) - denomfinal_SSFFdn->GetBinContent(idx);
        double valOSFFup = denomfinal_OSFFup->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
        double valOSFFdn = denomfinal_nominal->GetBinContent(idx) - denomfinal_OSFFdn->GetBinContent(idx);

        double valHeepIdUp = 0., valHeepIdDn = 0.;

        if ( TString(dataHist_->GetName()).Contains("invM") ) {
          valHeepIdUp = denomfinal_heepIdUp->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
          valHeepIdDn = denomfinal_nominal->GetBinContent(idx) - denomfinal_heepIdDn->GetBinContent(idx);
        }

        erryUp.push_back( std::sqrt( valSSFFup*valSSFFup + valOSFFup*valOSFFup + valHeepIdUp*valHeepIdUp ) );
        erryDn.push_back( std::sqrt( valSSFFdn*valSSFFdn + valOSFFdn*valOSFFdn + valHeepIdDn*valHeepIdDn ) );
      }

      auto gr = new TGraphAsymmErrors(denomfinal_nominal->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
      gr->SetFillColor(kGray+2);
      gr->SetLineColor(kGray+2);
      gr->SetFillStyle(3004);
      gr->Draw("2");

      if (dir_) {
        dir_->WriteTObject(dataHist_,"data_obs");
        dir_->WriteTObject(nominal_.returnSS().release(),"SS");
        dir_->WriteTObject(nominal_.returnOS().release(),"OS");

        dir_->WriteTObject(SSFFup_.returnSS().release(),"SS_mergedEleFakeFactorSSUp");
        dir_->WriteTObject(SSFFdn_.returnSS().release(),"SS_mergedEleFakeFactorSSDown");
        dir_->WriteTObject(OSFFup_.returnOS().release(),"OS_mergedEleFakeFactorOSUp");
        dir_->WriteTObject(OSFFdn_.returnOS().release(),"OS_mergedEleFakeFactorOSDown");
        dir_->WriteTObject(heepIdUp_.returnSS().release(),"SS_modHeepIdUp");
        dir_->WriteTObject(heepIdUp_.returnOS().release(),"OS_modHeepIdUp");
        dir_->WriteTObject(heepIdDn_.returnSS().release(),"SS_modHeepIdDown");
        dir_->WriteTObject(heepIdDn_.returnOS().release(),"OS_modHeepIdDown");

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          dir_->WriteTObject(sigHist_.at(idx),sigFiles_.at(idx).name_);
          dir_->WriteTObject(sigHist_modHeepIdUp_.at(idx),sigFiles_.at(idx).name_+"_modHeepIdUp");
          dir_->WriteTObject(sigHist_modHeepIdDn_.at(idx),sigFiles_.at(idx).name_+"_modHeepIdDown");
          dir_->WriteTObject(sigHist_mergedEleIdUp_.at(idx),sigFiles_.at(idx).name_+"_mergedEleIdUp");
          dir_->WriteTObject(sigHist_mergedEleIdDn_.at(idx),sigFiles_.at(idx).name_+"_mergedEleIdDown");
        }
      }

      denomSSfinal.release();
      denomOSfinal.release();
    }
  };

  canvas_1->cd();

  auto aloader3E = HistLoader3E(datafile,WZfile,ZZfile,sigsamples3E);

  canvas_1->SetLogy();
  aloader3E.load("3E_CRME_lll_invM","3E_antiME_lll_invM_CR");
  aloader3E.preparecard("ME3E_"+era+"_datacard.root","mergedEle3E");
  aloader3E.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_3E_invM.png");
  aloader3E.close();
  canvas_1->SetLogy(0);

  aloader3E.load("3E_Et_CR_EB_CRME","3E_Et_CR_EB_antiME");
  aloader3E.compare(canvas_1);
  SaveAs(canvas_1,"FF_3E_Et.png");

  aloader3E.load("3E_CRME_all_eta","3E_antiME_eta");
  aloader3E.compare(canvas_1);
  SaveAs(canvas_1,"FF_3E_eta.png");

  return;
}
