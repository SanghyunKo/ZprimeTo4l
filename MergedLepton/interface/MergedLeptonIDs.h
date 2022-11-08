#ifndef MergedLeptonIDs_H
#define MergedLeptonIDs_H 1

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

namespace MergedLeptonIDs {
  enum cutflowElectron {
    baseline = -1,
    minEnergy,
    ecalDriven,
    dPhiIn,
    saturatedXtal,
    missingInnerHits,
    trackIso,
    HoE,
    caloIso,
    dxy, // passed modified HEEP
    dEtaIn,
    showerShape, // passed HEEP
    failedHEEP = 20,
    no2ndGsf,
    has2ndGsf,
    failConvVeto,
    hasMatchedElectron, // has 2nd GSF track but it's from another electron
    passedMVA1, // merged electron MVA selection when there exists 2nd GSF track
    outsideDef, // outside Et or dR of interest
    passedMVA2, // merged electron MVA selection when there is no 2nd GSF track
    has2ndGsfCRpass, // CR with minEt < Et < etThresEB && pass MVA
    has2ndGsfCRfail, // CR with minEt < Et < etThresEB && fail MVA
    no2ndGsfCRpass, // CR with minEt < Et < etThresEB && pass MVA
    no2ndGsfCRfail  // CR with minEt < Et < etThresEB && fail MVA
  };

  enum GSFtype {
    nulltype = -1,
    DR1Et2EB,
    DR2Et1EB,
    DR2Et2EB,
    DR1Et2EE,
    DR2Et1EE,
    DR2Et2EE,
    extendedCR
  };

  inline bool isSameGsfTrack(const reco::GsfTrackRef& aTrk, const reco::GsfTrackRef& bTrk) {
    return ( aTrk->pt()==bTrk->pt() &&
             aTrk->eta()==bTrk->eta() &&
             aTrk->phi()==bTrk->phi() );
  }

  GSFtype checkElectronGSFtype( const pat::ElectronRef& aEle,
                                const reco::GsfTrackRef& orgGsfTrk,
                                const reco::GsfTrackRef& addGsfTrk,
                                const double etThresEB,
                                const double etThresEE,
                                const double minEt );

  void fillCutflow(TH1* ahist, const int acutflow, const double aWeight);

  bool isModifiedHEEP(const reco::GsfElectron& el,
                      const reco::Vertex& primaryVertex,
                      const float& trkIso,
                      const float& ecalIso,
                      const int& nrSatCrys,
                      const double& rho,
                      cutflowElectron& cutflow);

  bool hasPassedHEEP(const reco::GsfElectron& el,
                     const reco::Vertex& primaryVertex,
                     const int& nrSatCrys,
                     const double& rho,
                     cutflowElectron& cutflow);

  bool isHighPtTrackerMuon(const reco::Muon& muon,
                           const reco::Vertex& vtx);
} /* MergedLeptonIDs */

#endif
