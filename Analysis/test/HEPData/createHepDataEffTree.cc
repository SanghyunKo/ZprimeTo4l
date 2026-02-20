#include "TFile.h"
#include "TH1D.h"

#include <regex>

void createHepDataEffTree(TString filename) {
  // HepData output
  TTree* effTree = new TTree("effTree","effTree");
  float mX, logmY, mY, eff;
  effTree->Branch("mX",&mX,"mX/F");
  effTree->Branch("logmY",&logmY,"logmY/F");
  effTree->Branch("mY",&mY,"mY/F");
  effTree->Branch("eff",&eff,"eff/F");

  TFile* afile = TFile::Open(filename,"READ");
  auto keyList = afile->GetListOfKeys();

  // mass vectors
  std::vector<int> massvec = {250,275,300,325,350,375,400,425,450,500,550,650,750,850,1000,1250,1500,1750,2000};
  std::vector<double> avecDouble = {0.4,0.6,0.8,1.,1.5,2.,5.,10.,50.,100.};
  std::vector<std::pair<int, double>> filled;

  for (auto keyObj : *keyList) {
    auto key = (TKey*)keyObj;

    if (TString(key->GetClassName())==TString("TDirectoryFile")) {
      TDirectory* adir = (TDirectory*)key->ReadObj();
      TString dirname = adir->GetName();
      auto keyList2 = adir->GetListOfKeys();

      for (auto keyObj2 : *keyList2) {
        auto key2 = (TKey*)keyObj2;

        if (TString(key2->GetClassName()).Contains("TH1")) {
          std::string keyname = std::string(key2->GetName());
          std::regex sigPattern("H([0-9]+)([AZ])([0-9]*p*[0-9]+).*");
          std::regex sigSystPattern("H[0-9]+[AZ][0-9]*p*[0-9]+_(.*)");
          std::smatch XYmass, systName;

          bool isSig = std::regex_match(keyname, XYmass, sigPattern);
          bool isSyst = std::regex_match(keyname, systName, sigSystPattern);
          std::string ymassStr = std::regex_replace( XYmass[3].str(), std::regex("p"), "." );

          double mx = 0, my = 0;
          bool isZA = false;

          if (isSig) {
            mx = std::stod(XYmass[1].str());
            my = std::stod(ymassStr);
            isZA = XYmass[2].str()=="Z";
          }

          if (isSig && !isSyst && mx <= 2500. && isZA) {
            if (std::find(filled.begin(), filled.end(), std::make_pair(static_cast<int>(mx), my)) != filled.end()) {
              std::cout << "Already filled: " << keyname << std::endl;
              continue;
            }

            TH1* ahist = (TH1*)key2->ReadObj();
            eff = ahist->Integral()/(137.6*0.1);
            mX = mx;
            mY = my;
            logmY = std::log10(my);
            effTree->Fill();

            filled.push_back(std::make_pair(static_cast<int>(mx), my));
          }
        }
      } // loop hists


    } // if TDirectory
  } // key loop

  for (auto my : avecDouble) {
    double myDouble = my;

    for (auto mx : massvec) {
      if (std::find(filled.begin(), filled.end(), std::make_pair(mx, myDouble)) != filled.end())
        continue;

      // fill with zero for missing points
      mX = mx;
      mY = myDouble;
      logmY = std::log10(myDouble);
      eff = 0.;
      effTree->Fill();
    }
  }

  TFile* outFile = TFile::Open("HepDataEff_Z_"+filename,"RECREATE");
  effTree->Write();
  outFile->Close();
}
