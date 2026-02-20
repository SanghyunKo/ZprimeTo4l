#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"

#include "../verARC_v5/tdrstyle.C"
#include "../verARC_v5/CMS_lumi.C"

void combineFlavor_PRL_postFit() {
  TString era = "run2";
  setTDRStyle();

  customLeft = "";
  customRight = ""; // to adjust right margin
  customCmsTextOffset = 0.;
  customLumiOffset = 0.04;

  writeExtraText = true;       // if extra text
  extraText  = "";  // default extra text is "Preliminary"

  lumi_sqrtS = "";
  lumi_13TeV = "138 fb^{-1}";

  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 0;

  if( iPos==0 )
    relPosX = 0.12;

  int W = 3000;
  int H = 2100;

  int H_ref = 2100;
  int W_ref = 2000;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.08*H_ref;
  float L = 0.08*W_ref; // 0.1
  float R = 0.03*W_ref; // 0.02

  // gROOT->SetWebDisplay("firefox");

  auto* canvas1 = new TCanvas("canvas","canvas",W,H,W,H);
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
      pad->cd();
      pad->RedrawAxis();
      pad->GetFrame()->Draw();
    } else {
      canvas->Update();
      // canvas->RedrawAxis();
      // canvas->GetFrame()->Draw();
    }

    canvas->SaveAs(name.c_str());
  };

  class SignalRegion {
  public:
    class process {
    public:
      process(TH1F* ashape, TString aname, int acolor)
      : shape_(ashape), name_(aname), color_(acolor) {}
      ~process()=default;

      TH1F* shape_;
      TString name_;
      int color_;
    };

    TH1F* data_;
    TH1F* totBkg_;
    TH1F* sig_;
    std::vector<process> shapes_;
    TString SRname_;

  protected:
    TFile* file_;
    TFile* sigFile_;

    int sigColor_ = kRed;

  public:
    void divideByBinWidth(TH1F* ahist) {
      for (int i=1; i<ahist->GetNbinsX()+1; i++) {
        ahist->SetBinContent( i, ahist->GetBinContent(i)/ahist->GetBinWidth(i) );
        ahist->SetBinError( i , ahist->GetBinError(i)/ahist->GetBinWidth(i) );
      }
    }

    void init(TString aSRname) {
      SRname_ = aSRname;
      TGraphAsymmErrors* in_data = (TGraphAsymmErrors*)file_->Get("shapes_fit_b/"+aSRname+"/data");
      TH1F* in_totBkg = (TH1F*)file_->Get("shapes_fit_b/"+aSRname+"/total_background");

      std::vector<process> in_bkg;
      std::vector<double> binEdges;

      if (aSRname.Contains("mergedEMu2M") || aSRname.Contains("mergedEle3E")) {
        in_bkg.push_back(process((TH1F*)file_->Get("shapes_fit_b/"+aSRname+"/OS"),"\\mathrm{e_{ME}}\\hbox{ misID (Prompt)}",TColor::GetColor("#f89c20")));
        in_bkg.push_back(process((TH1F*)file_->Get("shapes_fit_b/"+aSRname+"/SS"),"\\mathrm{e_{ME}}\\hbox{ misID (Nonprompt)}",TColor::GetColor("#9c9ca1")));
        binEdges = {200.,220.,240.,260.,280.,300.,325.,350.,375.,400.,425.,450.,500.,550.,600.,700.,1000.,2500.};
        sig_ = (TH1F*)sigFile_->Get(aSRname+"/H550Z1");
      } else if ( aSRname.Contains("mergedEle2E") ) {
        in_bkg.push_back(process((TH1F*)file_->Get("shapes_fit_b/"+aSRname+"/OS"),"\\mathrm{e_{ME}}\\hbox{ misID (Prompt)}",TColor::GetColor("#f89c20")));
        binEdges ={200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,
                        400.,425.,450.,475.,500.,550.,600.,700.,800.,
                        900.,1000.,1200.,1500.,2500.};
        sig_ = (TH1F*)sigFile_->Get(aSRname+"/H550A1");
      } else if ( aSRname.Contains("resolved") ) {
        in_bkg.push_back(process((TH1F*)file_->Get("shapes_fit_b/"+aSRname+"/Nonprompt"),"\\hbox{Resolved }\\ell\\hbox{ misID }(\\Delta\\mathrm{R} > 0.3)",TColor::GetColor("#832db6")));
        in_bkg.push_back(process((TH1F*)file_->Get("shapes_fit_b/"+aSRname+"/NonpromptDR03"),"\\hbox{Resolved }\\ell\\hbox{ misID }(\\Delta\\mathrm{R} < 0.3)",TColor::GetColor("#92dadd")));
        in_bkg.push_back(process((TH1F*)file_->Get("shapes_fit_b/"+aSRname+"/ZZ"),"\\mathrm{ZZ}\\rightarrow4\\ell",TColor::GetColor("#5790fc")));
        binEdges = {200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,
                        400.,425.,450.,475.,500.,550.,600.,700.,800.,
                        900.,1000.,1200.,1500.,2500.};
        sig_ = (TH1F*)sigFile_->Get(aSRname+"/H550A100");
      } else if ( aSRname.Contains("mergedEMu1M") ) {
        in_bkg.push_back(process((TH1F*)file_->Get("shapes_fit_b/"+aSRname+"/MM"),"\\mu_{\\mathrm{MM}}\\hbox{ misID}",TColor::GetColor("#b9ac70")));
        binEdges = {250.,275.,300.,325.,350.,375.,400.,425.,450.,475.,500.,550.,600.,700.,800.,1000.,2500.};
        sig_ = (TH1F*)sigFile_->Get(aSRname+"/H2000A1");
      } else if ( aSRname.Contains("mergedMu3M") ) {
        in_bkg.push_back(process((TH1F*)file_->Get("shapes_fit_b/"+aSRname+"/MM"),"\\mu_{\\mathrm{MM}}\\hbox{ misID}",TColor::GetColor("#b9ac70")));
        binEdges = {250.,300.,350.,400.,450.,500.,600.,800.,2500.};
        sig_ = (TH1F*)sigFile_->Get(aSRname+"/H2000Z1");
      } else if ( aSRname.Contains("mergedMu2E") ) {
        in_bkg.push_back(process((TH1F*)file_->Get("shapes_fit_b/"+aSRname+"/MM"),"\\mu_{\\mathrm{MM}}\\hbox{ misID}",TColor::GetColor("#b9ac70")));
        binEdges = {250.,300.,350.,400.,450.,500.,600.,800.,1000.,2500.};
        sig_ = (TH1F*)sigFile_->Get(aSRname+"/H2000Z1");
      }

      data_ = new TH1F("out_data_"+aSRname,";M [GeV];",binEdges.size()-1,&(binEdges[0]));
      totBkg_ = new TH1F("out_totBkg_"+aSRname,";M [GeV];",binEdges.size()-1,&(binEdges[0]));
      sig_->SetLineColor(kRed);
      sig_->SetLineWidth(1);
      sig_->Scale(100.); // 100 = 10 fb
      sig_->SetMarkerColor(kRed);
      sig_->SetMarkerSize(0);

      for (const auto& aproc : in_bkg) {
        shapes_.push_back(process(new TH1F("out_bkg_"+aSRname+"_"+aproc.shape_->GetName(),";M [GeV];",binEdges.size()-1,&(binEdges[0])),aproc.name_,aproc.color_));
        shapes_.back().shape_->SetFillColor(aproc.color_);
        shapes_.back().shape_->SetLineColor(aproc.color_);
        shapes_.back().shape_->SetMarkerSize(0);
        shapes_.back().shape_->SetTitle(aproc.name_);
      }

      for (int ibin=1; ibin<=data_->GetNbinsX(); ibin++) {
        double x=0, y=0;
        in_data->GetPoint(ibin-1,x,y);
        data_->SetBinContent(ibin,y);
        data_->SetBinError(ibin,(in_data->GetErrorYhigh(ibin-1)+in_data->GetErrorYlow(ibin-1))/2.);

        totBkg_->SetBinContent(ibin,in_totBkg->GetBinContent(ibin));
        totBkg_->SetBinError(ibin,in_totBkg->GetBinError(ibin));

        for (unsigned idx=0; idx<in_bkg.size(); idx++) {
          shapes_.at(idx).shape_->SetBinContent(ibin,in_bkg.at(idx).shape_->GetBinContent(ibin));
          shapes_.at(idx).shape_->SetBinError(ibin,in_bkg.at(idx).shape_->GetBinError(ibin));
        }
      }
    } // init

  public:
    SignalRegion(TFile* afile, TFile* aSigFile)
    : file_(afile), sigFile_(aSigFile) {}

    ~SignalRegion()=default;

    void add(SignalRegion& other, SignalRegion* other2=nullptr) {
      data_->Add(other.data_);
      sig_->Add(other.sig_);

      for (unsigned idx=0; idx<shapes_.size(); idx++)
        shapes_.at(idx).shape_->Add(other.shapes_.at(idx).shape_);

      if (other2) {
        data_->Add(other2->data_);
        sig_->Add(other2->sig_);

        for (unsigned idx=0; idx<shapes_.size(); idx++)
          shapes_.at(idx).shape_->Add(other2->shapes_.at(idx).shape_);
      }

      // retrieve covariance matrix
      TH2D* covMat = (TH2D*)file_->Get("shapes_fit_b/overall_total_covar");

      // recompute total bkg
      for (unsigned ibin=1; ibin<totBkg_->GetNbinsX()+1; ibin++) {
        double binContent = totBkg_->GetBinContent(ibin) + other.totBkg_->GetBinContent(ibin);
        double err1 = totBkg_->GetBinError(ibin);
        double err2 = other.totBkg_->GetBinError(ibin);

        TString binName1 = SRname_ + "_" + TString( std::to_string(ibin-1) );
        TString binName2 = other.SRname_ + "_" + TString( std::to_string(ibin-1) );
        int binIdx1 = covMat->GetXaxis()->FindFixBin(binName1);
        int binIdx2 = covMat->GetYaxis()->FindFixBin(binName2);
        int binGlobal = covMat->GetBin(binIdx1,binIdx2);

        double cov = covMat->GetBinContent(binGlobal);
        double uncSq = err1*err1 + err2*err2 + 2.*cov;

        if (other2) {
          binContent += other2->totBkg_->GetBinContent(ibin);
          double err3 = other2->totBkg_->GetBinError(ibin);

          TString binName3 = other2->SRname_ + "_" + TString( std::to_string(ibin-1) );
          int binIdx3 = covMat->GetXaxis()->FindFixBin(binName3);
          int binGlobal13 = covMat->GetBin(binIdx1,binIdx3);
          int binGlobal23 = covMat->GetBin(binIdx2,binIdx3);

          double cov13 = covMat->GetBinContent(binGlobal13);
          double cov23 = covMat->GetBinContent(binGlobal23);

          uncSq += err3*err3 + 2.*(cov13 + cov23);
        }

        totBkg_->SetBinContent( ibin, binContent );
        totBkg_->SetBinError( ibin, uncSq > 0. ? sqrt(uncSq) : 0. );
      }
    } // add

    void addMergedMu(SignalRegion& other) {
      // assume this is mergedMu3M and the other is mergedMu2E
      std::vector<double> binEdges = {250.,300.,350.,400.,450.,500.,600.,800.,2500.};

      std::unique_ptr<TH1F> other_data = std::make_unique<TH1F>("mergedMu2E_dataRebinned","",binEdges.size()-1,&(binEdges[0]));
      std::unique_ptr<TH1F> other_totBkg = std::make_unique<TH1F>("mergedMu2E_totBkgRebinned","",binEdges.size()-1,&(binEdges[0]));
      std::unique_ptr<TH1F> other_sig = std::make_unique<TH1F>("mergedMu2E_sigRebinned","",binEdges.size()-1,&(binEdges[0]));
      std::vector<process> other_shapes;

      // rebin background processes (ignore error)
      for (unsigned ishape = 0; ishape<other.shapes_.size(); ishape++) {
        const auto& aproc = other.shapes_.at(ishape);
        other_shapes.push_back(process(new TH1F("mergedMu2E_shapeRebinned_"+aproc.name_,"",binEdges.size()-1,&(binEdges[0])),aproc.name_,aproc.color_));

        for (unsigned ibin = 0; ibin<other_data->GetNbinsX()+1; ibin++) {
          double binContent = aproc.shape_->GetBinContent(ibin);
          
          if (ibin==other_data->GetNbinsX()) // merge the last bin
            binContent += aproc.shape_->GetBinContent(ibin+1);

          other_shapes.back().shape_->SetBinContent(ibin, binContent);
          other_shapes.back().shape_->SetBinError(ibin,0.);
        }
      } // loop shapes

      // rebin data and signal (only data error matters)
      for (unsigned ibin = 0; ibin<other_data->GetNbinsX()+1; ibin++) {
        double binContent = other.data_->GetBinContent(ibin);
        double binError2 = other.data_->GetBinError(ibin)*other.data_->GetBinError(ibin);
        
        if (ibin==other_data->GetNbinsX()) { // merge the last bin
          binContent += other.data_->GetBinContent(ibin+1);
          binError2 += other.data_->GetBinError(ibin+1)*other.data_->GetBinError(ibin+1);
        }

        other_data->SetBinContent(ibin, binContent);
        other_data->SetBinError(ibin, sqrt(binError2));

        binContent = other.sig_->GetBinContent(ibin);
        
        if (ibin==other_sig->GetNbinsX()) // merge the last bin
          binContent += other.sig_->GetBinContent(ibin+1);

        other_sig->SetBinContent(ibin, binContent);
        other_sig->SetBinError(ibin,0.);
      } // loop bins of data and signal

      data_->Add(other_data.get());
      sig_->Add(other_sig.get());

      for (unsigned idx=0; idx<shapes_.size(); idx++)
        shapes_.at(idx).shape_->Add(other_shapes.at(idx).shape_);

      // retrieve covariance matrix
      TH2D* covMat = (TH2D*)file_->Get("shapes_fit_b/overall_total_covar");

      // recompute total bkg
      for (unsigned ibin=1; ibin<totBkg_->GetNbinsX()+1; ibin++) {
        double binContent = totBkg_->GetBinContent(ibin) + other.totBkg_->GetBinContent(ibin);
        double err1 = totBkg_->GetBinError(ibin);
        double err2 = other.totBkg_->GetBinError(ibin);

        TString binName1 = SRname_ + "_" + TString( std::to_string(ibin-1) );
        TString binName2 = other.SRname_ + "_" + TString( std::to_string(ibin-1) );
        int binIdx1 = covMat->GetXaxis()->FindFixBin(binName1);
        int binIdx2 = covMat->GetYaxis()->FindFixBin(binName2);
        int binGlobal = covMat->GetBin(binIdx1,binIdx2);

        double cov = covMat->GetBinContent(binGlobal);
        double uncSq = err1*err1 + err2*err2 + 2.*cov;

        if (ibin==totBkg_->GetNbinsX()) { // merge the last bin
          double err3 = other.totBkg_->GetBinError(ibin+1);

          TString binName3 = other.SRname_ + "_" + TString( std::to_string(ibin) );
          int binIdx3 = covMat->GetXaxis()->FindFixBin(binName3);
          int binGlobal13 = covMat->GetBin(binIdx1,binIdx3);
          int binGlobal23 = covMat->GetBin(binIdx2,binIdx3);

          double cov13 = covMat->GetBinContent(binGlobal13);
          double cov23 = covMat->GetBinContent(binGlobal23);

          binContent += other.totBkg_->GetBinContent(ibin+1);
          uncSq += err3*err3 + 2.*(cov13 + cov23);
        }

        totBkg_->SetBinContent( ibin, binContent );
        totBkg_->SetBinError( ibin, uncSq > 0. ? sqrt(uncSq) : 0. );
      }
    } // addMergedMu

    void divideByBinWidth() {
      divideByBinWidth(data_);
      divideByBinWidth(totBkg_);
      divideByBinWidth(sig_);

      for (unsigned idx=0; idx<shapes_.size(); idx++)
        divideByBinWidth(shapes_.at(idx).shape_);
    }

    std::unique_ptr<THStack> returnStack() {
      auto astack = std::make_unique<THStack>("stack_"+TString(data_->GetName()),data_->GetTitle());

      for (auto& x : shapes_)
        astack->Add(x.shape_);

      return std::move(astack);
    }

    std::unique_ptr<TH1F> returnRatio() {
      auto ratio = std::unique_ptr<TH1F>((TH1F*)data_->Clone());
      auto norm = std::unique_ptr<TH1F>((TH1F*)totBkg_->Clone());

      for(unsigned ibin=0; ibin < norm->GetNbinsX()+2; ibin++)
        norm->SetBinError(ibin, 0.); // stat uncertainty is part of nuisances

      ratio->Divide(norm.get());

      ratio->SetStats(0);
      ratio->SetTitle("");
      ratio->SetLineColor(kBlack);

      return std::move(ratio);
    }

    std::unique_ptr<TGraphAsymmErrors> returnError() {
      std::vector<double> x0, y0, errx, erryDn, erryUp;

      for (int idx=1; idx<=totBkg_->GetNbinsX(); idx++) {
        x0.push_back(totBkg_->GetBinCenter(idx));
        y0.push_back(totBkg_->GetBinContent(idx));
        errx.push_back(totBkg_->GetBinWidth(idx)/2.);

        erryUp.push_back(totBkg_->GetBinError(idx));
        erryDn.push_back(totBkg_->GetBinError(idx));
      }

      auto gr = std::make_unique<TGraphAsymmErrors>(totBkg_->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
      gr->SetFillColor(kGray+3);
      gr->SetLineColor(kGray+3);
      gr->SetFillStyle(3002);

      return std::move(gr);
    }

    std::unique_ptr<TGraphAsymmErrors> returnRatioError() {
      auto nominal = std::unique_ptr<TH1F>((TH1F*)totBkg_->Clone());
      std::vector<double> x0, y0, errx, erryDn, erryUp;

      for (int idx=1; idx<=nominal->GetNbinsX(); idx++) {
        x0.push_back(nominal->GetBinCenter(idx));
        y0.push_back(1.);
        errx.push_back(nominal->GetBinWidth(idx)/2.);

        double errUpSqr = 0.;
        double errDnSqr = 0.;
        double valNom = nominal->GetBinContent(idx);

        erryUp.push_back(valNom > 0. ? nominal->GetBinError(idx)/valNom : 0.);
        erryDn.push_back(valNom > 0. ? nominal->GetBinError(idx)/valNom : 0.);        
      }

      auto gr = std::make_unique<TGraphAsymmErrors>(nominal->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
      gr->SetFillColor(kGray+3);
      gr->SetLineColor(kGray+3);
      gr->SetFillStyle(3002);

      return std::move(gr);
    }
  };

  class ResolvedEleSR : public SignalRegion {
  public:
    ResolvedEleSR(TFile* afile, TFile* aSigFile)
    : SignalRegion(afile,aSigFile) {

      init("resolvedEle");
    }
  };

  class ResolvedEMuSR : public SignalRegion {
  public:
    ResolvedEMuSR(TFile* afile, TFile* aSigFile)
    : SignalRegion(afile,aSigFile) {

      init("resolvedEMu");
    }
  };

  class ResolvedMuSR : public SignalRegion {
  public:
    ResolvedMuSR(TFile* afile, TFile* aSigFile)
    : SignalRegion(afile,aSigFile) {

      init("resolvedMu");
    }
  };

  class MergedEle3ESR : public SignalRegion {
  public:
    MergedEle3ESR(TFile* afile, TFile* aSigFile)
    : SignalRegion(afile,aSigFile) {

      init("mergedEle3E");
    }
  };

  class MergedEMu2MSR : public SignalRegion {
  public:
    MergedEMu2MSR(TFile* afile, TFile* aSigFile)
    : SignalRegion(afile,aSigFile) {

      init("mergedEMu2M");
    }
  };

  class ME2ESR : public SignalRegion {
  public:
    ME2ESR(TFile* afile, TFile* aSigFile)
    : SignalRegion(afile,aSigFile) {

      init("mergedEle2E");
    }
  };

  class MergedMu3MSR : public SignalRegion {
  public:
    MergedMu3MSR(TFile* afile, TFile* aSigFile)
    : SignalRegion(afile,aSigFile) {

      init("mergedMu3M");
    }
  };

  class MergedMu2ESR : public SignalRegion {
  public:
    MergedMu2ESR(TFile* afile, TFile* aSigFile)
    : SignalRegion(afile,aSigFile) {

      init("mergedMu2E");
    }
  };

  class MEMu1MSR : public SignalRegion {
  public:
    MEMu1MSR(TFile* afile, TFile* aSigFile)
    : SignalRegion(afile,aSigFile) {

      init("mergedEMu1M");
    }
  };

  struct SubFigure {
  public:
    SubFigure()=default;
    SubFigure(TH1F* data, THStack* expected, TGraphAsymmErrors* err, TH1F* sig,
              TH1F* ratio, TGraphAsymmErrors* ratioErr)
        : data_(data), expected_(expected), err_(err), sig_(sig),
          ratio_(ratio), ratioErr_(ratioErr) {}

    TH1F* data_;
    THStack* expected_;
    TGraphAsymmErrors* err_;
    TH1F* sig_;
    TH1F* ratio_;
    TGraphAsymmErrors* ratioErr_;

    void writeHepDataFile(TString aname) {
      TFile* afile = new TFile(aname,"RECREATE");
      auto* data = (TH1D*)data_->Clone();
      auto* err = (TH1D*)err_->Clone();
      auto* sig = (TH1D*)sig_->Clone();
      afile->WriteTObject(data, "Data");
      afile->WriteTObject(err, "Total background");

      for (unsigned idx=0; idx<expected_->GetHists()->GetEntries(); idx++) {
        auto* hist = (TH1F*)expected_->GetHists()->At(idx)->Clone();
        afile->WriteTObject(hist, hist->GetName());
      }

      afile->WriteTObject(sig,"signal");
      afile->Close();
    }
  };

  // create hists for resolved SRs
  TFile* resolvedFile = new TFile("fitDiagnostics.fitDiagCWRX500Y100.root","READ");
  TFile* resolvedSigFile1 = new TFile("REFF_run2_datacard_rebin.root","READ");
  TFile* resolvedSigFile2 = new TFile("REMuFF_run2_datacard_rebin.root","READ");
  TFile* resolvedSigFile3 = new TFile("RMFF_run2_datacard_rebin.root","READ");
  auto resolvedEleSR = ResolvedEleSR(resolvedFile,resolvedSigFile1);
  auto resolvedEMuSR = ResolvedEMuSR(resolvedFile,resolvedSigFile2);
  auto resolvedMuSR = ResolvedMuSR(resolvedFile,resolvedSigFile3);
  resolvedEleSR.add(resolvedEMuSR,&resolvedMuSR);
  resolvedEleSR.divideByBinWidth();

  auto resolvedStack = resolvedEleSR.returnStack();
  auto resolvedRatio = resolvedEleSR.returnRatio();
  auto resolvedError = resolvedEleSR.returnError();
  auto resolvedRatioError = resolvedEleSR.returnRatioError();

  auto resolvedSubFigure = SubFigure(resolvedEleSR.data_,
                                     resolvedStack.get(),
                                     resolvedError.get(),
                                     resolvedEleSR.sig_,
                                     resolvedRatio.get(),
                                     resolvedRatioError.get());

  // create hists for merged SRs
  TFile* mergedFile = new TFile("fitDiagnostics.fitDiagCWRX500Z0p8.root","READ");
  TFile* mergedSigFile1 = new TFile("ME3E_run2_datacard_rebin.root","READ");
  TFile* mergedSigFile2 = new TFile("MEMu2M_run2_datacard_rebin.root","READ");
  auto mergedEle3ESR = MergedEle3ESR(mergedFile,mergedSigFile1);
  auto mergedEMu2MSR = MergedEMu2MSR(mergedFile,mergedSigFile2);
  mergedEle3ESR.add(mergedEMu2MSR);
  mergedEle3ESR.divideByBinWidth();

  auto mergedStack = mergedEle3ESR.returnStack();
  auto mergedRatio = mergedEle3ESR.returnRatio();
  auto mergedError = mergedEle3ESR.returnError();
  auto mergedRatioError = mergedEle3ESR.returnRatioError();

  auto mergedEleSubFigure = SubFigure(mergedEle3ESR.data_,
                                      mergedStack.get(),
                                      mergedError.get(),
                                      mergedEle3ESR.sig_,
                                      mergedRatio.get(),
                                      mergedRatioError.get());

  // create hists for ME2E SRs
  TFile* meeFile = mergedFile; // new TFile("fitDiagnostics.fitDiagCWRX500Y0p8.root","READ");
  TFile* meeSigFile = new TFile("ME2E_run2_datacard_rebin.root","READ");
  auto me2eSR = ME2ESR(meeFile,meeSigFile);
  me2eSR.divideByBinWidth();

  auto me2eStack = me2eSR.returnStack();
  auto me2eRatio = me2eSR.returnRatio();
  auto me2eError = me2eSR.returnError();
  auto me2eRatioError = me2eSR.returnRatioError();

  auto me2eSubFigure = SubFigure(me2eSR.data_,
                                 me2eStack.get(),
                                 me2eError.get(),
                                 me2eSR.sig_,
                                 me2eRatio.get(),
                                 me2eRatioError.get());

  // create hists for cleanedMu SRs
  TFile* mmffSigFile = new TFile("MMFF_run2_datacard_rebin.root","READ");
  TFile* mm2eSigFile = new TFile("MM2E_run2_datacard_rebin.root","READ");
  auto mergedMu3MSR = MergedMu3MSR(mergedFile,mmffSigFile);
  auto mergedMu2ESR = MergedMu2ESR(mergedFile,mm2eSigFile);
  mergedMu3MSR.addMergedMu(mergedMu2ESR);
  mergedMu3MSR.divideByBinWidth();

  auto mergedMuStack = mergedMu3MSR.returnStack();
  auto mergedMuRatio = mergedMu3MSR.returnRatio();
  auto mergedMuError = mergedMu3MSR.returnError();
  auto mergedMuRatioError = mergedMu3MSR.returnRatioError();

  auto mergedMuSubFigure = SubFigure(mergedMu3MSR.data_,
                                     mergedMuStack.get(),
                                     mergedMuError.get(),
                                     mergedMu3MSR.sig_,
                                     mergedMuRatio.get(),
                                     mergedMuRatioError.get());

  // create hists for MEMu1M SRs
  TFile* cmeSigFile = new TFile("MEMu1M_run2_datacard_rebin.root","READ");
  auto cleanedMu1MSR = MEMu1MSR(meeFile,cmeSigFile);
  cleanedMu1MSR.divideByBinWidth();
  auto cmeStack = cleanedMu1MSR.returnStack();
  auto cmeRatio = cleanedMu1MSR.returnRatio();
  auto cmeError = cleanedMu1MSR.returnError();
  auto cmeRatioError = cleanedMu1MSR.returnRatioError();

  auto cmeSubFigure = SubFigure(cleanedMu1MSR.data_,
                                cmeStack.get(),
                                cmeError.get(),
                                cleanedMu1MSR.sig_,
                                cmeRatio.get(),
                                cmeRatioError.get());

  // Define the rows and columns for subpads
  int rows = 2;
  int cols = 3;

  // Create subpads and divide the canvas
  double leftmargin = 1.5*L/W;
  double rightmargin = R/W;
  double topmargin = T/H;
  double bottommargin = 1.2*B/H;

  std::vector<SubFigure> subfigureList = {mergedEleSubFigure,mergedMuSubFigure,resolvedSubFigure,me2eSubFigure,cmeSubFigure};
  std::vector<TString> figureNames = {"2\\ell\\mathrm{e_{ME}}",
                                      "2\\ell\\mu\\mu_{\\mathrm{MM}}",
                                      "4\\ell",
                                      "2\\mathrm{e_{ME}}",
                                      "\\mathrm{e_{ME}}\\mu\\mu_{\\mathrm{MM}}"};
  std::vector<TString> figureFileNames = {"mergedEle","mergedMu","resolved","mergedEle2E","mergedEMu1M"};

  auto finalLegend = new TLegend(0.02,0.52,0.25,0.92);
  finalLegend->AddEntry(resolvedEleSR.data_,"Data","LEP");
  finalLegend->AddEntry(resolvedEleSR.shapes_.at(2).shape_,"\\mathrm{ZZ}\\rightarrow4\\ell");
  finalLegend->AddEntry(resolvedEleSR.shapes_.at(1).shape_,"\\hbox{Resolved }\\ell\\hbox{ misID }(\\Delta\\mathrm{R} < 0.3)");
  finalLegend->AddEntry(resolvedEleSR.shapes_.at(0).shape_,"\\hbox{Resolved }\\ell\\hbox{ misID }(\\Delta\\mathrm{R} > 0.3)");
  finalLegend->AddEntry(mergedEle3ESR.shapes_.at(0).shape_,"\\mathrm{e_{ME}}\\hbox{ misID (Prompt)}");
  finalLegend->AddEntry(mergedEle3ESR.shapes_.at(1).shape_,"\\mathrm{e_{ME}}\\hbox{ misID (Nonprompt)}");
  finalLegend->AddEntry(mergedMu3MSR.shapes_.at(0).shape_,"\\mu_{\\mathrm{MM}}\\hbox{ misID}");
  finalLegend->AddEntry(resolvedEleSR.sig_,"Signal");
  finalLegend->SetTextSize(0.03);
  finalLegend->SetBorderSize(0);
  finalLegend->SetFillStyle(0);
  finalLegend->SetMargin(0.2);
  finalLegend->SetTextFont(42);

  for (int j = 0; j < rows; ++j) {
    for (int i = 0; i < cols; ++i) {
      unsigned iPad = j*cols + i;
      unsigned iFig = iPad > 0 ? iPad-1 : iPad; // skip the (0,1) pad

      if (iFig >= subfigureList.size())
        continue;

      if (iPad == 0) {
        finalLegend->Draw();
        continue;
      }

      // Calculate the position of each subpad
      double xmin = leftmargin + static_cast<double>(i)*(1.-leftmargin-rightmargin)/cols;
      double xmax = leftmargin + static_cast<double>(i+1)*(1.-leftmargin-rightmargin)/cols;
      double ymin = bottommargin + static_cast<double>(rows-j-1)*(1.-topmargin-bottommargin)/rows;
      double ymax = bottommargin + static_cast<double>(rows-j)*(1.-topmargin-bottommargin)/rows;
      double ydiff = ymax - ymin;

      if (i==0)
        xmin -= 0.6*leftmargin;
      if (j==rows-1)
        ymin -= 0.3*bottommargin;

      if (iFig==0)
        xmin -= 0.61*leftmargin;

      if (i==1)
        xmax -= 0.15*leftmargin;
      if (i==2)
        xmin -= 0.15*leftmargin;

      // Create subpad
      TPad *padUp = new TPad(Form("padUp%d", iPad), Form("SubpadUp %d", iPad), xmin, ymin + ydiff*0.3, xmax, ymax);
      padUp->SetLeftMargin( 0. );
      padUp->SetRightMargin( 0.0 );
      padUp->SetTopMargin( 0.0 );
      padUp->SetBottomMargin( 0.0 );
      padUp->SetFillColor(0);
      padUp->SetBorderMode(0);
      padUp->SetFrameFillStyle(0);
      padUp->SetFrameBorderMode(0);
      TPad* padDn = new TPad(Form("padDn%d", iPad), Form("SubpadDn %d", iPad), xmin, ymin, xmax, ymin + ydiff*0.3);
      padDn->SetLeftMargin( 0. );
      padDn->SetRightMargin( 0.0 );
      padDn->SetTopMargin( 0.0 );
      padDn->SetBottomMargin( 0.0 );
      padDn->SetFillColor(0);
      padDn->SetBorderMode(0);
      padDn->SetFrameFillStyle(0);
      padDn->SetFrameBorderMode(0);

      auto aSubFigure = subfigureList.at(iFig);

      if (aSubFigure.data_ == nullptr)
        continue;

      // draw options
      TString axisOpt = "E1";

      if ( iFig==2 || iFig==0 ) {
        padUp->SetLeftMargin(1.8*leftmargin);
        padDn->SetLeftMargin(1.8*leftmargin);
        aSubFigure.data_->GetYaxis()->SetNdivisions(505);
        axisOpt += "Y-";
      } else if (iFig!=3) {
        padUp->SetRightMargin(1.7*leftmargin); // for labels
        padDn->SetRightMargin(1.7*leftmargin);
        axisOpt += "Y-"; // Y+
        // aSubFigure.data_->GetYaxis()->SetLabelOffset(-0.1);
      }
      //if (i==cols-1)
      //  pad->SetRightMargin(rightmargin);

      if (j==rows-1) {
        padDn->SetBottomMargin(2.2*bottommargin);
      }
      //if (j==0)
      //  pad->SetTopMargin(topmargin);

      padUp->Draw();
      padDn->Draw();

      padUp->cd();

      if (iFig == 2 || iFig == 3 || iFig==0 || iFig==1 || iFig==4) {
        padUp->SetLogy();
        padUp->SetLogx();
        padDn->SetLogx();
      } else {
        padUp->SetLogy(0);
      }

      // some cosmetics      
      aSubFigure.ratio_->GetYaxis()->SetTitle("");
      aSubFigure.data_->GetYaxis()->SetTitle("");

      if (iFig==2 || iFig==3 || iFig==0 || iFig==1 || iFig==4) {
        aSubFigure.data_->SetMaximum(500.);
        aSubFigure.data_->SetMinimum(50e-5);
      }

      double xmaxUser = 2500.; // (iFig==1 || iFig==4) ? 2500. : 1500.;
      double xminUser = (iFig==1 || iFig==4) ? 250. : 200.;

      aSubFigure.data_->GetXaxis()->SetRangeUser(xminUser,xmaxUser);
      aSubFigure.ratio_->GetXaxis()->SetRangeUser(xminUser,xmaxUser);
      aSubFigure.ratio_->GetXaxis()->SetTickLength(0.08);
      // aSubFigure.ratio_->GetYaxis()->SetTickLength(0.05);
      aSubFigure.ratio_->GetXaxis()->SetTitle("");
      aSubFigure.ratio_->GetXaxis()->SetLabelSize(0.25);
      aSubFigure.ratio_->GetYaxis()->SetLabelSize(0.25);
      aSubFigure.ratio_->GetYaxis()->SetLabelOffset(0.01);
      aSubFigure.ratio_->GetYaxis()->SetNdivisions(505); // ratio Y
      // aSubFigure.ratio_->GetXaxis()->SetNdivisions(707);
      // aSubFigure.data_->GetXaxis()->SetNdivisions(707);
      aSubFigure.data_->GetYaxis()->SetLabelSize(0.1);
      aSubFigure.data_->GetYaxis()->SetLabelOffset(0.01);

      // if (iFig==1 || iFig==4) {
      //   aSubFigure.ratio_->GetXaxis()->SetNdivisions(504);
      //   aSubFigure.data_->GetXaxis()->SetNdivisions(504);
      //   aSubFigure.data_->GetYaxis()->SetNdivisions(iFig==4 ? 508 : 505);
      // }

      // if (iFig==3 || iFig==4)
        // aSubFigure.ratio_->GetXaxis()->ChangeLabel(1,-1,0.);

      // if (iFig==2 || iFig==3)
        // aSubFigure.ratio_->GetXaxis()->ChangeLabel(-1,-1,0.);

      aSubFigure.ratio_->GetYaxis()->SetRangeUser(0.2,1.8);
      aSubFigure.ratio_->GetXaxis()->SetMoreLogLabels();
      aSubFigure.ratio_->GetXaxis()->SetNoExponent();

      aSubFigure.data_->Draw(axisOpt);
      aSubFigure.expected_->Draw("hist&same");
      aSubFigure.data_->Draw(axisOpt+"&same");
      aSubFigure.sig_->Draw("hist&same");
      aSubFigure.err_->Draw("2");

      padDn->cd();
      aSubFigure.ratio_->Draw(axisOpt);
      aSubFigure.ratioErr_->Draw("2");

      padUp->RedrawAxis();
      padDn->RedrawAxis();

      double xoffset = (i==2) ? -0.1 : 0.;

      if (iFig==0)
        xoffset += 0.1;

      if (iFig==2)
        xoffset += 0.05;

      padUp->cd();

      double offsetCustom = 0.;

      if (iFig==2)
        offsetCustom = -0.05;
      else if (iFig==0)
        offsetCustom = -0.1;

      if (i==2)
        offsetCustom -= 0.01;

      auto atext = std::make_unique<TPaveText>(0.5+xoffset,0.8,0.95+xoffset+offsetCustom,0.95,"NDC");
      atext->SetBorderSize(0);
      atext->SetFillColor(0);
      atext->SetFillStyle(0);
      atext->AddText(figureNames.at(iFig));
      ((TText*)atext->GetListOfLines()->Last())->SetTextColor(kBlack);
      ((TText*)atext->GetListOfLines()->Last())->SetTextAlign(31);
      ((TText*)atext->GetListOfLines()->Last())->SetTextSize(0.11);
      atext->Draw();
      atext.release();

      offsetCustom = 0.;

      TString model = (j==1) ? "\\mathrm{X}\\rightarrow\\mathrm{YY}" : "\\mathrm{X}\\rightarrow\\mathrm{ZY}";

      if (iFig==2)
        offsetCustom = -0.05;
      else if (iFig==0)
        offsetCustom = -0.1;

      if (i==2)
        offsetCustom -= 0.01;

      auto dtext = std::make_unique<TPaveText>(0.5+xoffset,0.5,0.95+xoffset+offsetCustom,0.68,"NDC");
      dtext->SetBorderSize(0);
      dtext->SetFillColor(0);
      dtext->SetFillStyle(0);

      dtext->AddText(model);
      ((TText*)dtext->GetListOfLines()->Last())->SetTextColor(kBlack);
      ((TText*)dtext->GetListOfLines()->Last())->SetTextAlign(31);
      ((TText*)dtext->GetListOfLines()->Last())->SetTextSize(0.1);
      ((TText*)dtext->GetListOfLines()->Last())->SetTextFont(42);
      // ((TText*)atext->GetListOfLines()->Last())->SetTextAlign(12);
      dtext->Draw();
      dtext.release();

      TString massPoint = "(M_{X}, M_{Y}) = (550, 1) GeV";

      if (iFig==2)
        massPoint = "(M_{X}, M_{Y}) = (550, 100) GeV";
      else if (i==2)
        massPoint = "(M_{X}, M_{Y}) = (2000, 1) GeV";

      offsetCustom = 0.;

      if (iFig==0)
        offsetCustom = -0.08;
      else if (iFig==2)
        offsetCustom = -0.035;
      else
        offsetCustom = 0.02;

      if (i==2)
      offsetCustom -= 0.01;

      auto etext = std::make_unique<TPaveText>(0.5+xoffset-0.35,0.68,0.95+xoffset+offsetCustom,0.8,"NDC");
      etext->SetBorderSize(0);
      etext->SetFillColor(0);
      etext->SetFillStyle(0);
      etext->AddText(massPoint);
      ((TText*)etext->GetListOfLines()->Last())->SetTextColor(kBlack);
      ((TText*)etext->GetListOfLines()->Last())->SetTextAlign(31);
      ((TText*)etext->GetListOfLines()->Last())->SetTextSize(0.09);
      ((TText*)etext->GetListOfLines()->Last())->SetTextFont(42);
      etext->Draw();
      etext.release();

      canvas1->cd();  // Go back to the main canvas
    }
  }

  canvas1->cd();

  auto ytext = std::make_unique<TPaveText>(0.0,0.13,0.2,0.48,"NDC");
  ytext->SetBorderSize(0);
  ytext->SetFillColor(0);
  ytext->SetFillStyle(0);
  TString astr = "N_{events}/GeV";
  ytext->AddText(astr);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextFont(42);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextSize(0.045);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextAngle(90);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextAlign(13);
  ytext->Draw();

  auto ytextRatio = std::make_unique<TPaveText>(0.01,0.02,0.2,0.12,"NDC");
  ytextRatio->SetBorderSize(0);
  ytextRatio->SetFillColor(0);
  ytextRatio->SetFillStyle(0);
  astr = "Obs/Exp";
  ytextRatio->AddText(astr);
  ((TText*)ytextRatio->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)ytextRatio->GetListOfLines()->Last())->SetTextFont(42);
  ((TText*)ytextRatio->GetListOfLines()->Last())->SetTextSize(0.035);
  ((TText*)ytextRatio->GetListOfLines()->Last())->SetTextAngle(90);
  ((TText*)ytextRatio->GetListOfLines()->Last())->SetTextAlign(13);
  ytextRatio->Draw();

  auto xtext = std::make_unique<TPaveText>(0.82,0.01,0.99,0.06,"NDC");
  xtext->SetBorderSize(0);
  xtext->SetFillColor(0);
  xtext->SetFillStyle(0);
  TString xstr = "\\mathrm{M}_{\\mathrm{T}} [\\mathrm{GeV}]";
  xtext->AddText(xstr);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextFont(42);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextAlign(11);
  xtext->Draw();

  auto xtextInv = std::make_unique<TPaveText>(0.53,0.01,0.85,0.06,"NDC");
  xtextInv->SetBorderSize(0);
  xtextInv->SetFillColor(0);
  xtextInv->SetFillStyle(0);
  TString xstrInv = "\\mathrm{M}_{4\\ell} [\\mathrm{GeV}]";
  xtextInv->AddText(xstrInv);
  ((TText*)xtextInv->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)xtextInv->GetListOfLines()->Last())->SetTextFont(42);
  ((TText*)xtextInv->GetListOfLines()->Last())->SetTextAlign(11);
  xtextInv->Draw();

  canvas1->Update();  // Refresh canvas

  SaveAs(canvas1,std::string((TString("signalRegions")+"_PRL_"+era+"_postFit.png").Data()));

  // save hepdata files
  std::cout << "Saving HepData files..." << std::endl;
  resolvedSubFigure.writeHepDataFile("HepData_figure1_resolvedSR.root");
  mergedEleSubFigure.writeHepDataFile("HepData_figure1_mergedEleSR.root");
  me2eSubFigure.writeHepDataFile("HepData_figure1_me2eSR.root");
  mergedMuSubFigure.writeHepDataFile("HepData_figure1_mergedMuSR.root");
  cmeSubFigure.writeHepDataFile("HepData_figure1_cmeSR.root");

  // save individual canvases for each SR
  // create new canvas
  customCmsTextOffset = 0.;
  customLumiOffset = 0.;

  iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  iPos = 11;
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
  R = 0.04*W_ref; // 0.02

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

  // create subpads for each SR
  TPad* padUp = new TPad("padUp","",0,0.3,1,0.95);
  padUp->SetFillColor(0);
  padUp->SetFrameBorderSize(0);
  padUp->SetBorderMode(0);
  padUp->SetFrameFillStyle(0);
  padUp->SetFrameBorderMode(0);
  padUp->SetTickx(0);
  padUp->SetTicky(0);
  padUp->SetBottomMargin(0.00);
  padUp->SetLeftMargin( L/W );
  padUp->SetRightMargin( R/W );
  padUp->Draw();

  TPad* padDn = new TPad("padDn","",0,0,1,0.3);
  padDn->SetTickx(0);
  padDn->SetTicky(0);
  padDn->SetTopMargin(0.0);
  padDn->SetBottomMargin(0.3);
  padDn->SetLeftMargin( L/W );
  padDn->SetRightMargin( R/W );
  padDn->Draw();

  // now draw each SR in a separate canvas and save
  for (unsigned idx=0; idx<subfigureList.size(); idx++) {
    auto aSubFigure = subfigureList.at(idx);

    if (aSubFigure.data_ == nullptr)
      continue;

    padUp->cd();
    padUp->SetLogy();
    padUp->SetLogx();
    padDn->SetLogx();

    aSubFigure.data_->Draw("E1");
    aSubFigure.data_->GetYaxis()->SetTitleSize(0.08);
    aSubFigure.data_->GetYaxis()->SetTitle("N_{events}/GeV");
    aSubFigure.data_->GetYaxis()->SetTitleOffset(0.7);
    aSubFigure.data_->GetYaxis()->SetLabelSize(0.07);
    aSubFigure.expected_->Draw("hist&same");
    aSubFigure.data_->Draw("E1&same");
    aSubFigure.sig_->SetLineWidth(2);
    aSubFigure.sig_->Draw("hist&same");
    aSubFigure.err_->Draw("2");

    TLegend* leg = nullptr;

    if (figureFileNames.at(idx).Contains("resolved")) {
      leg = new TLegend(0.4,0.6,0.9,0.92);
      aSubFigure.data_->SetMaximum(aSubFigure.data_->GetMaximum()*2.);
    } else if (figureFileNames.at(idx).Contains("mergedEle2E")) {
      leg = new TLegend(0.55,0.65,0.9,0.92);
    } else if (figureFileNames.at(idx).Contains("Mu")) {
      leg = new TLegend(0.7,0.65,0.9,0.92);
    } else
      leg = new TLegend(0.45,0.6,0.9,0.92);

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.055);
    leg->SetMargin(0.2);
    leg->AddEntry(aSubFigure.data_,"Data","LEP");

    for (unsigned i=0; i<aSubFigure.expected_->GetHists()->GetEntries(); i++) {
      auto* hist = (TH1F*)aSubFigure.expected_->GetHists()->At(i);
      leg->AddEntry(hist, hist->GetTitle(), "f");
    }

    leg->AddEntry(aSubFigure.sig_,"Signal","L");
    leg->Draw();

    padDn->cd();
    aSubFigure.ratio_->Draw("E1");
    aSubFigure.ratio_->GetXaxis()->SetTitleSize(0.16);
    aSubFigure.ratio_->GetXaxis()->SetLabelSize(0.14);
    TString xTitle = figureFileNames.at(idx).Contains("Mu") ? "M_{T} [GeV]" : "\\mathrm{M}_{4\\ell}~[\\mathrm{GeV}]";
    aSubFigure.ratio_->GetXaxis()->SetTitle(xTitle);
    aSubFigure.ratio_->GetXaxis()->SetTitleOffset(0.8);
    aSubFigure.ratio_->GetYaxis()->SetTitleSize(0.16);
    aSubFigure.ratio_->GetYaxis()->SetLabelSize(0.14);
    aSubFigure.ratio_->GetYaxis()->SetTitle("Obs/Exp");
    aSubFigure.ratio_->GetYaxis()->SetTitleOffset(0.35);
    aSubFigure.ratioErr_->Draw("2");

    padUp->RedrawAxis();
    padDn->RedrawAxis();

    SaveAs(canvas2,std::string((TString("HepData_figure1_") + figureFileNames.at(idx) + "_" + era + "_postFit.png").Data()));
  }

  return;
}
