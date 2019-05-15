#ifndef _Dhadron_h
#define _Dhadron_h

// system include files
#include <memory>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <string>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// To access the high level physics objects:
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
// Packed Candidates, pfs
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
// MC Truth
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// TFile:
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
// root
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
// Reclustering
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
// hlt triggers
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
// l1 trigger
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/L1Trigger/interface/Tau.h"

// class declaration

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
class Dhadron : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit Dhadron(const edm::ParameterSet&);
    ~Dhadron();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT <std::vector<reco::Vertex>>      vtxToken_;
    edm::EDGetTokenT <pat::ElectronCollection>        electronToken_;
    edm::EDGetTokenT <pat::TauCollection>             tauToken_;
    edm::EDGetTokenT <pat::MuonCollection>            muonToken_;
    edm::EDGetTokenT <pat::PhotonCollection>          photonToken_;
    edm::EDGetTokenT <pat::JetCollection>             jetToken_;
    edm::EDGetTokenT <pat::METCollection>             metToken_;
    edm::EDGetTokenT <pat::JetCollection>             fatjetToken_;
    edm::EDGetTokenT <pat::PackedCandidateCollection> pfToken_;

    edm::EDGetTokenT <pat::JetCollection>             ak1PFCHSjetToken_;
    
    edm::EDGetTokenT <edm::TriggerResults>                        triggerBits_;
    edm::EDGetTokenT <std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;
    edm::EDGetTokenT <pat::PackedTriggerPrescales>                triggerPrescales_;
    
    edm::EDGetTokenT <l1t::EGammaBxCollection> l1EGammasToken_;
    edm::EDGetTokenT <l1t::MuonBxCollection  > l1MuonsToken_;
    edm::EDGetTokenT <l1t::JetBxCollection   > l1JetsToken_;
    edm::EDGetTokenT <l1t::EtSumBxCollection > l1EtSumsToken_ ;
    edm::EDGetTokenT <l1t::TauBxCollection   > l1TauToken_;

    edm::EDGetTokenT <reco::GenJetCollection      > ak4genjetToken_;     //slimmedGenJets
    edm::EDGetTokenT <reco::GenParticleCollection > prunedGenParticles_; //prunedGenParticles --> all pythia level gen particles with c an c bars, higgs
    edm::EDGetTokenT <pat::PackedGenParticleCollection > packedGenParticles_; //packedGenParticles --> all pythia level final state particles
  
	TH1F *h_l1Tau_n                             ;
	TH1F *h_l1Tau_pt                            ;
	TH1F *h_l1Tau_eta                           ;
	TH1F *h_l1Tau_phi                           ;
	TH1F *h_l1Tau_mt                            ;
	TH1F *h_l1Tau_energy                        ;
	TH1F *h_l1Tau_isoEt                         ;
	TH1F *h_l1Tau_hwIso                         ;

	TH1F *h_matched_ak4chs_n                    ;
	TH1F *h_matched_ak4chs_pt                   ;
	TH1F *h_matched_ak4chs_mass                 ;
	TH1F *h_matched_ak4chs_neutralMultiplicity  ;
	TH1F *h_matched_ak4chs_chargedMultiplicity  ;
	TH1F *h_matched_ak4chs_numberOfDaughters    ;
	TH1F *h_matched_ak4chs_partonFlavour        ;
	TH1F *h_matched_ak4chs_jetArea              ;
	TH1F *h_matched_ak4chs_charge               ;
	TH1F *h_matched_ak4chs_isolation            ;

	TH1F *h_selected_ak4chs_n                   ;
	TH1F *h_selected_ak4chs_pt                  ;
	TH1F *h_selected_ak4chs_mass                ;
	TH1F *h_selected_ak4chs_neutralMultiplicity ;
	TH1F *h_selected_ak4chs_chargedMultiplicity ;
	TH1F *h_selected_ak4chs_numberOfDaughters   ;
	TH1F *h_selected_ak4chs_partonFlavour       ;
	TH1F *h_selected_ak4chs_jetArea             ;
	TH1F *h_selected_ak4chs_charge              ;
	TH1F *h_selected_ak4chs_isolation           ;

	TH1F *h_matched_ak1chs_n                    ;
	TH1F *h_matched_ak1chs_pt                   ;
	TH1F *h_matched_ak1chs_mass                 ;
	TH1F *h_matched_ak1chs_neutralMultiplicity  ;
	TH1F *h_matched_ak1chs_chargedMultiplicity  ;
	TH1F *h_matched_ak1chs_numberOfDaughters    ;
	TH1F *h_matched_ak1chs_partonFlavour        ;
	TH1F *h_matched_ak1chs_jetArea              ;
	TH1F *h_matched_ak1chs_charge               ;
	TH1F *h_matched_ak1chs_isolation            ;

	TH1F *h_isolated_ak1chs_n                   ;
	// TH1F *h_isolated_ak1chs_pt                  ;
	// TH1F *h_isolated_ak1chs_mass                ;
	// TH1F *h_isolated_ak1chs_neutralMultiplicity ;
	// TH1F *h_isolated_ak1chs_chargedMultiplicity ;
	// TH1F *h_isolated_ak1chs_numberOfDaughters   ;
	// TH1F *h_isolated_ak1chs_partonFlavour       ;
	// TH1F *h_isolated_ak1chs_jetArea             ;
	// TH1F *h_isolated_ak1chs_charge              ;
	// TH1F *h_isolated_ak1chs_isolation           ;
	// TH1F *h_isolated_ak1chs_higgsMass           ;

	TH1F *h_selected_ak1chs_n                   ;
	TH1F *h_selected_ak1chs_pt                  ;
	TH1F *h_selected_ak1chs_mass                ;
	TH1F *h_selected_ak1chs_neutralMultiplicity ;
	TH1F *h_selected_ak1chs_chargedMultiplicity ;
	TH1F *h_selected_ak1chs_numberOfDaughters   ;
	TH1F *h_selected_ak1chs_partonFlavour       ;
	TH1F *h_selected_ak1chs_jetArea             ;
	TH1F *h_selected_ak1chs_charge              ;
	TH1F *h_selected_ak1chs_isolation           ;
	TH1F *h_selected_ak1chs_higgsMass           ;
};
#endif
