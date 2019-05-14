// -*- C++ -*-
//
// Package:    hTocc/Dhadron
// Class:      Dhadron
// 
/**\class Dhadron Dhadron.cc hTocc/Dhadron/plugins/Dhadron.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nabin Poudyal
//         Created:  Fri, 03 May 2019 17:30:43 GMT
//
//

// system include files
#include <memory>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include "hTocc/Dhadron/plugins/Dhadron.h"
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

// constants, enums and typedefs

// static data member definitions

// constructors and destructor
Dhadron::Dhadron(const edm::ParameterSet& iConfig)
  :
  vtxToken_     (consumes<std::vector<reco::Vertex>>      (iConfig.getParameter<edm::InputTag>("vertices"))),
  electronToken_(consumes<pat::ElectronCollection>        (iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_     (consumes<pat::TauCollection>             (iConfig.getParameter<edm::InputTag>("taus"))),
  muonToken_    (consumes<pat::MuonCollection>            (iConfig.getParameter<edm::InputTag>("muons"))),
  photonToken_  (consumes<pat::PhotonCollection>          (iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_     (consumes<pat::JetCollection>             (iConfig.getParameter<edm::InputTag>("AK4CHS"))),
  metToken_     (consumes<pat::METCollection>             (iConfig.getParameter<edm::InputTag>("mets"))),
  fatjetToken_  (consumes<pat::JetCollection>             (iConfig.getParameter<edm::InputTag>("fatjets"))),
  pfToken_      (consumes<pat::PackedCandidateCollection> (iConfig.getParameter<edm::InputTag>("pfCands"))),

  ak1PFCHSjetToken_(consumes<pat::JetCollection>(edm::InputTag("selectedPatJetsAK1PFCHS"))),

  triggerBits_     (consumes<edm::TriggerResults>                       (edm::InputTag("TriggerResults","","HLT"))),
  triggerObjects_  (consumes<std::vector<pat::TriggerObjectStandAlone>> (edm::InputTag("selectedPatTrigger"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>               (edm::InputTag("patTrigger"))),

  l1EGammasToken_ (consumes<l1t::EGammaBxCollection>(edm::InputTag("caloStage2Digis:EGamma"))),
  l1MuonsToken_   (consumes<l1t::MuonBxCollection  >(edm::InputTag("gmtStage2Digis:Muon"))),
  l1JetsToken_    (consumes<l1t::JetBxCollection   >(edm::InputTag("caloStage2Digis:Jet"))),
  l1EtSumsToken_  (consumes<l1t::EtSumBxCollection >(edm::InputTag("caloStage2Digis:EtSum"))),
  l1TauToken_     (consumes<l1t::TauBxCollection   >(edm::InputTag("caloStage2Digis:Tau")))

  // ak4genjetToken_(consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJets")))
{
    //now do what ever initialization is needed
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    // all info about L1 tau object per event
    h_l1Tau_n                             = fs->make<TH1F>("h_l1Tau_n",     "h_l1Tau_n",     20,0.,20.)  ;
    h_l1Tau_pt                            = fs->make<TH1F>("h_l1Tau_pt",    "h_l1Tau_pt" ,   200,0.,200.);
    h_l1Tau_eta                           = fs->make<TH1F>("h_l1Tau_eta",   "h_l1Tau_eta",   32,-4.,4.)  ;
    h_l1Tau_phi                           = fs->make<TH1F>("h_l1Tau_phi",   "h_l1Tau_phi",   32,-4.,4.)  ;
    h_l1Tau_mt                            = fs->make<TH1F>("h_l1Tau_mt",    "h_l1Tau_mt",    200,0.,200) ;
    h_l1Tau_energy                        = fs->make<TH1F>("h_l1Tau_energy","h_l1Tau_energy",200,0.,200) ;
    h_l1Tau_isoEt                         = fs->make<TH1F>("h_l1Tau_isoEt", "h_l1Tau_isoEt", 200,0.,200) ;
    h_l1Tau_hwIso                         = fs->make<TH1F>("h_l1Tau_hwIso", "h_l1Tau_hwIso", 200,0.,2)   ;
    // if an event has 2 or 3 L1 tau, then below is the info about matched ak4chs jets with L1 tau
    h_matched_ak4chs_n                    = fs->make<TH1F>("h_matched_ak4chs_n",                    "h_matched_ak4chs_n",                  20,0.,20.)  ;
    h_matched_ak4chs_pt                   = fs->make<TH1F>("h_matched_ak4chs_pt",                   "h_matched_ak4chs_pt" ,                200,0.,200.); 
    h_matched_ak4chs_mass                 = fs->make<TH1F>("h_matched_ak4chs_mass",                 "h_matched_ak4chs_mass" ,              200,0.,200.); 
    h_matched_ak4chs_neutralMultiplicity  = fs->make<TH1F>("h_matched_ak4chs_neutralMultiplicity",  "h_matched_ak4chs_neutralMultiplicity",50,0.,50.)  ;
    h_matched_ak4chs_chargedMultiplicity  = fs->make<TH1F>("h_matched_ak4chs_chargedMultiplicity",  "h_matched_ak4chs_chargedMultiplicity",50,0.,50.)  ;
    h_matched_ak4chs_numberOfDaughters    = fs->make<TH1F>("h_matched_ak4chs_numberOfDaughters",    "h_matched_ak4chs_numberOfDaughters",  50,0.,50.)  ;
    h_matched_ak4chs_partonFlavour        = fs->make<TH1F>("h_matched_ak4chs_partonFlavour",        "h_matched_ak4chs_partonFlavour",      25,0.,25.)  ;
    h_matched_ak4chs_jetArea              = fs->make<TH1F>("h_matched_ak4chs_jetArea",              "h_matched_ak4chs_jetArea",            100,0.,1.)  ;
    h_matched_ak4chs_charge               = fs->make<TH1F>("h_matched_ak4chs_charge",               "h_matched_ak4chs_charge",             40,-20.,20.);
    h_matched_ak4chs_isolation            = fs->make<TH1F>("h_matched_ak4chs_isolation",            "h_matched_ak4chs_isolation",          200,0.,2.)  ;
    h_matched_ak4chs_higgsMass            = fs->make<TH1F>("h_matched_ak4chs_higgsMass",            "h_matched_ak4chs_higgsMass",          200,0.,200.);
    // if an event has 2 or 3 matched ak4csh jet then below is the info about matched ak1chs jets with matched ak4chs with matched L1 tau
    h_selected_ak4chs_n                   = fs->make<TH1F>("h_selected_ak4chs_n",                    "h_selected_ak4chs_n",                  20,0.,20.)  ;
    h_selected_ak4chs_pt                  = fs->make<TH1F>("h_selected_ak4chs_pt",                   "h_selected_ak4chs_pt" ,                200,0.,200.); 
    h_selected_ak4chs_mass                = fs->make<TH1F>("h_selected_ak4chs_mass",                 "h_selected_ak4chs_mass" ,              200,0.,200.); 
    h_selected_ak4chs_neutralMultiplicity = fs->make<TH1F>("h_selected_ak4chs_neutralMultiplicity",  "h_selected_ak4chs_neutralMultiplicity",50,0.,50.)  ;
    h_selected_ak4chs_chargedMultiplicity = fs->make<TH1F>("h_selected_ak4chs_chargedMultiplicity",  "h_selected_ak4chs_chargedMultiplicity",50,0.,50.)  ;
    h_selected_ak4chs_numberOfDaughters   = fs->make<TH1F>("h_selected_ak4chs_numberOfDaughters",    "h_selected_ak4chs_numberOfDaughters",  50,0.,50.)  ;
    h_selected_ak4chs_partonFlavour       = fs->make<TH1F>("h_selected_ak4chs_partonFlavour",        "h_selected_ak4chs_partonFlavour",      25,0.,25.)  ;
    h_selected_ak4chs_jetArea             = fs->make<TH1F>("h_selected_ak4chs_jetArea",              "h_selected_ak4chs_jetArea",            100,0.,1.)  ;
    h_selected_ak4chs_charge              = fs->make<TH1F>("h_selected_ak4chs_charge",               "h_selected_ak4chs_charge",             40,-20.,20.);
    h_selected_ak4chs_isolation           = fs->make<TH1F>("h_selected_ak4chs_isolation",            "h_selected_ak4chs_isolation",          200,0.,2.)  ;
    h_selected_ak4chs_higgsMass           = fs->make<TH1F>("h_selected_ak4chs_higgsMass",            "h_selected_ak4chs_higgsMass",          200,0.,200.);
    // matched with selected ak4chs
    h_matched_ak1chs_n                    = fs->make<TH1F>("h_matched_ak1chs_n",                    "h_matched_ak1chs_n",                  20,0.,20.)  ;
    h_matched_ak1chs_pt                   = fs->make<TH1F>("h_matched_ak1chs_pt",                   "h_matched_ak1chs_pt" ,                200,0.,200.); 
    h_matched_ak1chs_mass                 = fs->make<TH1F>("h_matched_ak1chs_mass",                 "h_matched_ak1chs_mass" ,              200,0.,200.); 
    h_matched_ak1chs_neutralMultiplicity  = fs->make<TH1F>("h_matched_ak1chs_neutralMultiplicity",  "h_matched_ak1chs_neutralMultiplicity",50,0.,50.)  ;
    h_matched_ak1chs_chargedMultiplicity  = fs->make<TH1F>("h_matched_ak1chs_chargedMultiplicity",  "h_matched_ak1chs_chargedMultiplicity",50,0.,50.)  ;
    h_matched_ak1chs_numberOfDaughters    = fs->make<TH1F>("h_matched_ak1chs_numberOfDaughters",    "h_matched_ak1chs_numberOfDaughters",  50,0.,50.)  ;
    h_matched_ak1chs_partonFlavour        = fs->make<TH1F>("h_matched_ak1chs_partonFlavour",        "h_matched_ak1chs_partonFlavour",      25,0.,25.)  ;
    h_matched_ak1chs_jetArea              = fs->make<TH1F>("h_matched_ak1chs_jetArea",              "h_matched_ak1chs_jetArea",            100,0.,1.)  ;
    h_matched_ak1chs_charge               = fs->make<TH1F>("h_matched_ak1chs_charge",               "h_matched_ak1chs_charge",             40,-20.,20.);
    h_matched_ak1chs_isolation            = fs->make<TH1F>("h_matched_ak1chs_isolation",            "h_matched_ak1chs_isolation",          200,0.,2.)  ;
    h_matched_ak1chs_higgsMass            = fs->make<TH1F>("h_matched_ak1chs_higgsMass",            "h_matched_ak1chs_higgsMass",          200,0.,200.);
    // if the matched ak1 is isolated with 10%
    h_isolated_ak1chs_n                   = fs->make<TH1F>("h_isolated_ak1chs_n",                    "h_isolated_ak1chs_n",                  20,0.,20.)  ;
    h_isolated_ak1chs_pt                  = fs->make<TH1F>("h_isolated_ak1chs_pt",                   "h_isolated_ak1chs_pt" ,                200,0.,200.); 
    h_isolated_ak1chs_mass                = fs->make<TH1F>("h_isolated_ak1chs_mass",                 "h_isolated_ak1chs_mass" ,              200,0.,200.); 
    h_isolated_ak1chs_neutralMultiplicity = fs->make<TH1F>("h_isolated_ak1chs_neutralMultiplicity",  "h_isolated_ak1chs_neutralMultiplicity",50,0.,50.)  ;
    h_isolated_ak1chs_chargedMultiplicity = fs->make<TH1F>("h_isolated_ak1chs_chargedMultiplicity",  "h_isolated_ak1chs_chargedMultiplicity",50,0.,50.)  ;
    h_isolated_ak1chs_numberOfDaughters   = fs->make<TH1F>("h_isolated_ak1chs_numberOfDaughters",    "h_isolated_ak1chs_numberOfDaughters",  50,0.,50.)  ;
    h_isolated_ak1chs_partonFlavour       = fs->make<TH1F>("h_isolated_ak1chs_partonFlavour",        "h_isolated_ak1chs_partonFlavour",      25,0.,25.)  ;
    h_isolated_ak1chs_jetArea             = fs->make<TH1F>("h_isolated_ak1chs_jetArea",              "h_isolated_ak1chs_jetArea",            100,0.,1.)  ;
    h_isolated_ak1chs_charge              = fs->make<TH1F>("h_isolated_ak1chs_charge",               "h_isolated_ak1chs_charge",             40,-20.,20.);
    h_isolated_ak1chs_isolation           = fs->make<TH1F>("h_isolated_ak1chs_isolation",            "h_isolated_ak1chs_isolation",          200,0.,2.)  ;
    h_isolated_ak1chs_higgsMass           = fs->make<TH1F>("h_isolated_ak1chs_higgsMass",            "h_isolated_ak1chs_higgsMass",          200,0.,200.);
    // if isolated ak1 has at least two jets
    h_selected_ak1chs_n                   = fs->make<TH1F>("h_selected_ak1chs_n",                    "h_selected_ak1chs_n",                  20,0.,20.)  ;
    h_selected_ak1chs_pt                  = fs->make<TH1F>("h_selected_ak1chs_pt",                   "h_selected_ak1chs_pt" ,                200,0.,200.); 
    h_selected_ak1chs_mass                = fs->make<TH1F>("h_selected_ak1chs_mass",                 "h_selected_ak1chs_mass" ,              200,0.,200.); 
    h_selected_ak1chs_neutralMultiplicity = fs->make<TH1F>("h_selected_ak1chs_neutralMultiplicity",  "h_selected_ak1chs_neutralMultiplicity",50,0.,50.)  ;
    h_selected_ak1chs_chargedMultiplicity = fs->make<TH1F>("h_selected_ak1chs_chargedMultiplicity",  "h_selected_ak1chs_chargedMultiplicity",50,0.,50.)  ;
    h_selected_ak1chs_numberOfDaughters   = fs->make<TH1F>("h_selected_ak1chs_numberOfDaughters",    "h_selected_ak1chs_numberOfDaughters",  50,0.,50.)  ;
    h_selected_ak1chs_partonFlavour       = fs->make<TH1F>("h_selected_ak1chs_partonFlavour",        "h_selected_ak1chs_partonFlavour",      25,0.,25.)  ;
    h_selected_ak1chs_jetArea             = fs->make<TH1F>("h_selected_ak1chs_jetArea",              "h_selected_ak1chs_jetArea",            100,0.,1.)  ;
    h_selected_ak1chs_charge              = fs->make<TH1F>("h_selected_ak1chs_charge",               "h_selected_ak1chs_charge",             40,-20.,20.);
    h_selected_ak1chs_isolation           = fs->make<TH1F>("h_selected_ak1chs_isolation",            "h_selected_ak1chs_isolation",          200,0.,2.)  ;
    h_selected_ak1chs_higgsMass           = fs->make<TH1F>("h_selected_ak1chs_higgsMass",            "h_selected_ak1chs_higgsMass",          200,0.,200.);
}

Dhadron::~Dhadron()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void Dhadron::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace pat;
  //
  edm::Handle<pat::JetCollection> AK1CHS; 
  iEvent.getByToken(ak1PFCHSjetToken_, AK1CHS);

  edm::Handle<pat::JetCollection> AK4CHS;
  iEvent.getByToken(jetToken_, AK4CHS); // jets are stored in decending order of pt

  edm::Handle<pat::TauCollection>  taus;
  iEvent.getByToken(tauToken_, taus);
  
  edm::Handle<pat::PackedCandidateCollection> pfCands;
  iEvent.getByToken(pfToken_, pfCands);
  
  edm::Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  //
  edm::Handle<l1t::EGammaBxCollection> L1EGammas; 
  iEvent.getByToken(l1EGammasToken_, L1EGammas);
  
  edm::Handle<l1t::MuonBxCollection> L1Muons; 
  iEvent.getByToken(l1MuonsToken_, L1Muons);
  
  edm::Handle<l1t::JetBxCollection > L1Jets; 
  iEvent.getByToken(l1JetsToken_, L1Jets);
  
  edm::Handle<l1t::EtSumBxCollection> L1EtSums; 
  iEvent.getByToken(l1EtSumsToken_, L1EtSums);
  
  edm::Handle<l1t::TauBxCollection> L1Taus; 
  iEvent.getByToken(l1TauToken_, L1Taus);

  //
  if (!taus.isValid()) {
    edm::LogWarning("Dhadron") << "no pat::Tau in event";
    return;
  }

  if (!AK1CHS.isValid()) {
    edm::LogWarning("Dhadron") << "no pat::AK1CHS in event";
    return;
  }

  if (!AK4CHS.isValid()) {
    edm::LogWarning("Dhadron") << "no pat::AK4CHS in event";
    return;
  }

  // working with L1 tau objects
  h_l1Tau_n->Fill((*L1Taus.product()).size());     
  // vector<l1t::Tau> myTauObjVector; myTauObjVector.clear();
  for (unsigned int bx=0; bx<1; bx++){
    for(unsigned int i=0; i < (*L1Taus.product()).size(bx); i++ ){
        const l1t::Tau &myTauObj = (*L1Taus.product()).at(bx,i);
        // if (myTauObj.pt() >= 28.0 && abs(myTauObj.eta()) < 2.1){
          h_l1Tau_pt     ->Fill(myTauObj.pt());
          h_l1Tau_eta    ->Fill(myTauObj.eta());
          h_l1Tau_phi    ->Fill(myTauObj.phi());
          h_l1Tau_mt     ->Fill(myTauObj.mt());
          h_l1Tau_energy ->Fill(myTauObj.energy());
          h_l1Tau_isoEt  ->Fill(myTauObj.isoEt());
          h_l1Tau_hwIso  ->Fill(myTauObj.hwIso());
          // myTauObjVector.push_back(myTauObj);
        //}
    }
  }
  //
  if((*AK1CHS.product()).size() <=1 || (*AK4CHS.product()).size() <=1 || (*L1Taus.product()).size() <=1 || (*L1Taus.product()).size() >=4) return;
  // finding the matched ak4 and filling them
  vector<pat::Jet> matched_ak4chsJets; //matched_ak4chsJets.clear();
  for (const pat::Jet &ijet : *AK4CHS) {  
    // if(ijet.pt() <= 28. || fabs(ijet.eta()) >= 2.1 || ijet.numberOfDaughters() <= 1 || ijet.numberOfDaughters() >=8 || ijet.mass() <= 1.8 || ijet.mass() >= 4) continue;
    if(ijet.pt() <= 28. || fabs(ijet.eta()) >= 2.1 || ijet.mass() <= 1.8 ) continue;
    for (const l1t::Tau &itau : *L1Taus) {
      if (deltaR(ijet.eta(), ijet.phi(), itau.eta(), itau.phi()) <= 0.15) {
        matched_ak4chsJets.push_back(ijet); // matched
        double charged = 0, neutral = 0, pileup  = 0;
        std::vector< reco::CandidatePtr> constituents(ijet.daughterPtrVector());  // it contains all the daughters from the single jet
        if (constituents.size() == 0) continue;
        // std::sort(constituents.begin(), constituents.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // sorting the daughters by decending pt order
        for (unsigned int i = 0, n = pfCands->size(); i < n; ++i) {
          const pat::PackedCandidate &pf = (*pfCands)[i];
          // cout << "test: " << pf.pdgId() << endl;
          if (deltaR(pf.eta(), pf.phi(), ijet.eta(), ijet.phi()) < 1) {
            // pfcandidate-based constituents removal
            if (std::find(constituents.begin(), constituents.end(), reco::CandidatePtr(pfCands,i)) != constituents.end()) continue;
            if (pf.charge() == 0) {
              if (pf.pt() > 0.5) neutral += pf.pt();
            } else if (pf.fromPV() >= 2) {
              charged += pf.pt();
            } else {
              if (pf.pt() > 0.5) pileup += pf.pt();
            }
          }
        }
        double iso = (charged + std::max(0.0, neutral-0.5*pileup))/ijet.pt();
        h_matched_ak4chs_pt                  ->Fill(ijet.pt());
        h_matched_ak4chs_mass                ->Fill(ijet.mass()); 
        h_matched_ak4chs_neutralMultiplicity ->Fill(ijet.neutralMultiplicity()); 
        h_matched_ak4chs_chargedMultiplicity ->Fill(ijet.chargedMultiplicity ()); 
        h_matched_ak4chs_numberOfDaughters   ->Fill(ijet.numberOfDaughters()); 
        h_matched_ak4chs_partonFlavour       ->Fill(ijet.partonFlavour()); 
        h_matched_ak4chs_jetArea             ->Fill(ijet.jetArea()); 
        h_matched_ak4chs_charge              ->Fill(ijet.charge()); 
        h_matched_ak4chs_isolation           ->Fill(iso);  
        break;
      }
    }
  }
  // filling the matched jets and selecting depending on isolation
  h_matched_ak4chs_n->Fill(matched_ak4chsJets.size());                    
  // for (const pat::Jet &ijet : matched_ak4chsJets) {
        
  // }
  //
    // h_matched_ak4chs_higgsMass               

  // ak1jets
  // vector<TLorentzVector> higgsMass1Total; higgsMass1Total.clear();
  // vector<TLorentzVector> higgsMass1Total1; higgsMass1Total1.clear();
  
  // vector<pat::Jet> myselectedAK1jets; myselectedAK1jets.clear();

  // for (const pat::Jet &ijet : *AK1CHS) {
  //   if(ijet.pt() <= 28. || fabs(ijet.eta()) >= 2.1 || ijet.numberOfDaughters() <= 1 || ijet.numberOfDaughters() >=8 || ijet.mass() <= 1.5 || ijet.mass() >= 4) continue;
  //   for (const l1t::Tau &itau : *L1Taus) {
  //     if (deltaR(ijet.eta(), ijet.phi(), itau.eta(), itau.phi()) >= 0.15) continue; // did not match
    
  //     double charged = 0, neutral = 0, pileup  = 0;
  //     std::vector< reco::CandidatePtr> constituents(ijet.daughterPtrVector());  // it contains all the daughters from the single jet
  //     if (constituents.size() == 0) continue;
  //     std::sort(constituents.begin(), constituents.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // sorting the daughters by decending pt order
  //     for (unsigned int i = 0; i < constituents.size(); ++i) {
  //       const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*constituents[i]); // object daughter object
  //       // printf(" constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", i,cand.pt(),cand.dz(PV.position()),cand.pdgId());
  //     }

  //     for (unsigned int i = 0, n = pfCands->size(); i < n; ++i) {
  //       const pat::PackedCandidate &pf = (*pfCands)[i];
  //       // cout << "test: " << pf.pdgId() << endl;
  //       if (deltaR(pf.eta(), pf.phi(), ijet.eta(), ijet.phi()) < 0.5) {
  //         // pfcandidate-based constituents removal
  //         if (std::find(constituents.begin(), constituents.end(), reco::CandidatePtr(pfCands,i)) != constituents.end()) continue;
  //         if (pf.charge() == 0) {
  //           if (pf.pt() > 0.5) neutral += pf.pt();
  //         } else if (pf.fromPV() >= 2) {
  //           charged += pf.pt();
  //         } else {
  //           if (pf.pt() > 0.5) pileup += pf.pt();
  //         }
  //       }
  //     }
  //     double iso = (charged + std::max(0.0, neutral-0.5*pileup))/ijet.pt();
  //     h_ak1chs_Iso_dR5->Fill(iso);
  //     if (iso >= 0.15) continue; // do things below for isolated jets only
  //     // working with numbers of daughters in the jet
  //     // for (unsigned int i = 0; i < constituents.size(); ++i) {
  //     //   const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*constituents[i]); // object daughter object
  //     //   // printf(" constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", i,cand.pt(),cand.dz(PV.position()),cand.pdgId());
  //     //   cout << cand.pdgId() << endl;

  //     // }
  //     myselectedAK1jets.push_back(ijet);
  //     constituents.clear(); 
  //     TLorentzVector higgsMass1(0,0,0,0);
  //     higgsMass1.SetPtEtaPhiM(ijet.pt(),ijet.eta(),ijet.phi(),ijet.mass());
  //     higgsMass1Total.push_back(higgsMass1);
  //     h_ak1chs_pt               ->Fill(ijet.pt());
  //     h_ak1chs_mass             ->Fill(ijet.mass());
  //     h_ak1chs_numberOfDaughters->Fill(ijet.numberOfDaughters());
  //     h_ak1chs_partonFlavour    ->Fill(ijet.partonFlavour());

  //     if (iso <= 0.01) {
  //       TLorentzVector higgsMass1(0,0,0,0);
  //       higgsMass1.SetPtEtaPhiM(ijet.pt(),ijet.eta(),ijet.phi(),ijet.mass());
  //       higgsMass1Total1.push_back(higgsMass1);
  //       h_ak4_pt_01               ->Fill(ijet.pt());
  //       h_ak4_mass_01             ->Fill(ijet.mass());
  //       h_ak4_numberOfDaughters_01->Fill(ijet.numberOfDaughters());
  //       h_ak4chs_partonFlavour_01 ->Fill(ijet.partonFlavour());
  //     }
  //     break;
  //   }
  // }
  // cout << "test, size of the myselectedAK1jets " << myselectedAK1jets.size() << endl;
  // cout << "first pt of the myselectedAK1jets: " << myselectedAK1jets[0].pt() << endl;

  // if (higgsMass1Total.size() >= 2) h_ak1chs_higgsMass->Fill((higgsMass1Total[0]+higgsMass1Total[1]).M());
  // if (higgsMass1Total1.size() >= 2) h_ak1chs_higgsMass1->Fill((higgsMass1Total1[0]+higgsMass1Total1[1]).M());





  // edm::Handle<pat::MuonCollection>  muons;
  // iEvent.getByToken(muonToken_, muons);
  // // for (const pat::Muon &mu : *muons) {
  // //   if (mu.pt() < 5 || !mu.isLooseMuon()) continue;
  // //   printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
  // //           mu.pt(), mu.muonBestTrack()->dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV));
  // // }

  // edm::Handle<pat::ElectronCollection>  electrons;
  // iEvent.getByToken(electronToken_, electrons);
  // // for (const pat::Electron &el : *electrons) {
  // //   if (el.pt() < 5) continue;
  // //   printf("elec with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), lost hits %d, pass conv veto %d\n",
  // //         el.pt(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(), el.passConversionVeto());
  // // }

  // edm::Handle<pat::PhotonCollection>  photons;
  // iEvent.getByToken(photonToken_, photons);
  // // for (const pat::Photon &pho : *photons) {
  // //   if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
  // //   printf("phot with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes)\n",
  // //           pho.pt(), pho.superCluster()->eta(), pho.sigmaIetaIeta(), pho.full5x5_sigmaIetaIeta());
  // // }

 

  // // int ijet = 0;
  // // for (const pat::Jet &j : *jets) {
  // //   if (j.pt() < 20) continue;
  // //   printf("jet  with pt %5.1f (raw pt %5.1f), eta %+4.2f, btag CSV %.3f, CISV %.3f, pileup mva disc %+.2f\n",
  // //         j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")), j.userFloat("pileupJetId:fullDiscriminant"));
  // //   // if ((++ijet) == 1) { // for the first jet, let's print the leading constituents
  // //   //   std::vector daus(j.daughterPtrVector());
  // //   //   std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // the joys of C++11
  // //   //   for (unsigned int i2 = 0, n = daus.size(); i2 < n && i2 <= 3; ++i2) {
  // //   //     const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
  // //   //     printf("         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", i2,cand.pt(),cand.dz(PV.position()),cand.pdgId());
  // //   //   }
  // //   // }
  // // }

  // edm::Handle<pat::JetCollection>  fatjets;
  // iEvent.getByToken(fatjetToken_, fatjets);
  // // for (const pat::Jet &j : *fatjets) {
  // //   printf("AK8j with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed, %5.1f softdrop, %5.1f pruned CHS\n",
  // //         j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), j.mass(), j.userFloat("ak8PFJetsPuppiSoftDropMass"), j.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass") );

  // //   // To get the constituents of the AK8 jets, you have to loop over all of the
  // //   // daughters recursively. To save space, the first two constituents are actually
  // //   // the Soft Drop SUBJETS, which will then point to their daughters.
  // //   // The remaining constituents are those constituents removed by soft drop but
  // //   // still in the AK8 jet.
  // //   std::vector constituents;
  // //   for ( unsigned ida = 0; ida < j.numberOfDaughters(); ++ida ) {
  // //     reco::Candidate const * cand = j.daughter(ida);
  // //     if ( cand->numberOfDaughters() == 0 )
  // //         constituents.push_back( cand ) ;
  // //     else {
  // //       for ( unsigned jda = 0; jda < cand->numberOfDaughters(); ++jda ) {
  // //         reco::Candidate const * cand2 = cand->daughter(jda);
  // //         constituents.push_back( cand2 );
  // //       }
  // //     }
  // //   }
  // //   std::sort( constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda){return ida->pt() > jda->pt();} );

  // //   for ( unsigned int ida = 0; ida < constituents.size(); ++ida ) {
  // //     const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ida]);
  // //     printf("         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", ida,cand.pt(),cand.dz(PV.position()),cand.pdgId());
  // //   }

  // //   auto sdSubjets = j.subjets("SoftDropPuppi");
  // //   for ( auto const & isd : sdSubjets ) {
  // //     printf("  sd subjet with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed\n",
  // //     isd->pt(), isd->pt()*isd->jecFactor("Uncorrected"), isd->eta(), isd->mass() );

  // //   }
  // // }

  // edm::Handle<pat::METCollection>  mets;
  // iEvent.getByToken(metToken_, mets);
  // // const pat::MET &met = mets->front();
  // // printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
  // //   met.pt(), met.phi(), met.sumEt(),
  // //   met.genMET()->pt(),
  // //   met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));

  // // printf("\n");
}
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif

// ------------ method called once each job just before starting event loop  ------------
void Dhadron::beginJob(){
}

// ------------ method called once each job just after ending the event loop  ------------
void Dhadron::endJob(){
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Dhadron::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Dhadron);
