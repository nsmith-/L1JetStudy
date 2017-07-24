// -*- C++ -*-
//
// Package:    Analysis/L1JetTupler
// Class:      L1JetTupler
// 
/**\class L1JetTupler L1JetTupler.cc Analysis/L1JetTupler/plugins/L1JetTupler.cc

 Description: 

 Implementation:
     
*/
//
// Original Author:  Nicholas Charles Smith
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include "DataFormats/Math/interface/deltaR.h"

namespace {
  struct L1JetStruct {
    Long64_t run;
    Long64_t lumi;
    Long64_t event;
    int iJet;
    float recoJet_pt;
    float recoJet_eta;
    float recoJet_phi;
    float recoJet_energy;
    float recoJet_chargedHadronEnergy;
    float recoJet_neutralHadronEnergy;
    float recoJet_electronEnergy;
    float recoJet_photonEnergy;
    float recoJet_muonEnergy;
    float recoJet_hfHadEnergy;
    float recoJet_hfEmEnergy;
    float l1Jet_matchDeltaR;
    float l1Jet_et;
    float l1Jet_eta;
    float l1Jet_phi;
    float l1Jet_hwPt;
    float l1Jet_hwEta;
    float l1Jet_hwPhi;
    float l1Jet_towerIEta;
    float l1Jet_towerIPhi;
    float l1Jet_rawEt;
    float l1Jet_seedEt;
    float l1Jet_puEt;
    std::vector<float> caloTower_hwEt;
    std::vector<float> caloTower_ieta;
    std::vector<float> caloTower_iphi;
    std::vector<float> hcalTower_et;
    std::vector<float> hcalTower_eta;
    std::vector<float> hcalTower_phi;
    std::vector<float> hcalTower_hwEt;
    std::vector<float> hcalTower_ieta;
    std::vector<float> hcalTower_iphi;
    std::vector<float> ecalTower_et;
    std::vector<float> ecalTower_eta;
    std::vector<float> ecalTower_phi;
    std::vector<float> ecalTower_hwEt;
    std::vector<float> ecalTower_ieta;
    std::vector<float> ecalTower_iphi;
  };
}

class L1JetTupler : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>  {
  public:
    explicit L1JetTupler(const edm::ParameterSet&);
    ~L1JetTupler();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
    virtual void endJob() override;

    edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalTPs_;
    edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalTPs_;
    edm::EDGetTokenT<l1t::CaloTowerBxCollection> caloTowers_;
    edm::EDGetTokenT<l1t::JetBxCollection> l1Jets_;
    edm::EDGetTokenT<pat::JetCollection> recoJets_;

    TTree * l1JetTree_;
    L1JetStruct l1Jet_;
    edm::ESHandle<CaloGeometry> caloGeometry_;
};

L1JetTupler::L1JetTupler(const edm::ParameterSet& iConfig):
  ecalTPs_(consumes<EcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("ecalTPs"))),
  hcalTPs_(consumes<HcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("hcalTPs"))),
  caloTowers_(consumes<l1t::CaloTowerBxCollection>(iConfig.getParameter<edm::InputTag>("caloTowers"))),
  l1Jets_(consumes<l1t::JetBxCollection>(iConfig.getParameter<edm::InputTag>("l1Jets"))),
  recoJets_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("patJets")))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  l1JetTree_ = fs->make<TTree>("l1Jet", "All gen electrons, nearest reco");
  l1JetTree_->Branch("run", &l1Jet_.run);
  l1JetTree_->Branch("lumi", &l1Jet_.lumi);
  l1JetTree_->Branch("event", &l1Jet_.event);
  l1JetTree_->Branch("iJet", &l1Jet_.iJet);
  l1JetTree_->Branch("recoJet_pt", &l1Jet_.recoJet_pt);
  l1JetTree_->Branch("recoJet_eta", &l1Jet_.recoJet_eta);
  l1JetTree_->Branch("recoJet_phi", &l1Jet_.recoJet_phi);
  l1JetTree_->Branch("recoJet_energy", &l1Jet_.recoJet_energy);
  l1JetTree_->Branch("recoJet_chargedHadronEnergy", &l1Jet_.recoJet_chargedHadronEnergy);
  l1JetTree_->Branch("recoJet_neutralHadronEnergy", &l1Jet_.recoJet_neutralHadronEnergy);
  l1JetTree_->Branch("recoJet_electronEnergy", &l1Jet_.recoJet_electronEnergy);
  l1JetTree_->Branch("recoJet_photonEnergy", &l1Jet_.recoJet_photonEnergy);
  l1JetTree_->Branch("recoJet_muonEnergy", &l1Jet_.recoJet_muonEnergy);
  l1JetTree_->Branch("recoJet_hfHadEnergy", &l1Jet_.recoJet_hfHadEnergy);
  l1JetTree_->Branch("recoJet_hfEmEnergy", &l1Jet_.recoJet_hfEmEnergy);
  l1JetTree_->Branch("l1Jet_matchDeltaR", &l1Jet_.l1Jet_matchDeltaR);
  l1JetTree_->Branch("l1Jet_et", &l1Jet_.l1Jet_et);
  l1JetTree_->Branch("l1Jet_eta", &l1Jet_.l1Jet_eta);
  l1JetTree_->Branch("l1Jet_phi", &l1Jet_.l1Jet_phi);
  l1JetTree_->Branch("l1Jet_hwPt", &l1Jet_.l1Jet_hwPt);
  l1JetTree_->Branch("l1Jet_hwEta", &l1Jet_.l1Jet_hwEta);
  l1JetTree_->Branch("l1Jet_hwPhi", &l1Jet_.l1Jet_hwPhi);
  l1JetTree_->Branch("l1Jet_towerIEta", &l1Jet_.l1Jet_towerIEta);
  l1JetTree_->Branch("l1Jet_towerIPhi", &l1Jet_.l1Jet_towerIPhi);
  l1JetTree_->Branch("l1Jet_rawEt", &l1Jet_.l1Jet_rawEt);
  l1JetTree_->Branch("l1Jet_seedEt", &l1Jet_.l1Jet_seedEt);
  l1JetTree_->Branch("l1Jet_puEt", &l1Jet_.l1Jet_puEt);
  l1JetTree_->Branch("caloTower_hwEt", &l1Jet_.caloTower_hwEt);
  l1JetTree_->Branch("caloTower_ieta", &l1Jet_.caloTower_ieta);
  l1JetTree_->Branch("caloTower_iphi", &l1Jet_.caloTower_iphi);
  l1JetTree_->Branch("hcalTower_et", &l1Jet_.hcalTower_et);
  l1JetTree_->Branch("hcalTower_eta", &l1Jet_.hcalTower_eta);
  l1JetTree_->Branch("hcalTower_phi", &l1Jet_.hcalTower_phi);
  l1JetTree_->Branch("hcalTower_hwEt", &l1Jet_.hcalTower_hwEt);
  l1JetTree_->Branch("hcalTower_ieta", &l1Jet_.hcalTower_ieta);
  l1JetTree_->Branch("hcalTower_iphi", &l1Jet_.hcalTower_iphi);
  l1JetTree_->Branch("ecalTower_et", &l1Jet_.ecalTower_et);
  l1JetTree_->Branch("ecalTower_eta", &l1Jet_.ecalTower_eta);
  l1JetTree_->Branch("ecalTower_phi", &l1Jet_.ecalTower_phi);
  l1JetTree_->Branch("ecalTower_hwEt", &l1Jet_.ecalTower_hwEt);
  l1JetTree_->Branch("ecalTower_ieta", &l1Jet_.ecalTower_ieta);
  l1JetTree_->Branch("ecalTower_iphi", &l1Jet_.ecalTower_iphi);
}


L1JetTupler::~L1JetTupler()
{
}

void
L1JetTupler::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  iSetup.get<CaloGeometryRecord>().get(caloGeometry_);
}

void
L1JetTupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::ESHandle<CaloTPGTranscoder> decoder;
  iSetup.get<CaloTPGRecord>().get(decoder);

  l1Jet_.run = iEvent.run();
  l1Jet_.lumi = iEvent.luminosityBlock();
  l1Jet_.event = iEvent.id().event();

  Handle<EcalTrigPrimDigiCollection> ecalTPs;
  iEvent.getByToken(ecalTPs_, ecalTPs);

  Handle<HcalTrigPrimDigiCollection> hcalTPs;
  iEvent.getByToken(hcalTPs_, hcalTPs);

  Handle<l1t::CaloTowerBxCollection> caloTowers;
  iEvent.getByToken(caloTowers_, caloTowers);

  Handle<l1t::JetBxCollection> l1Jets;
  iEvent.getByToken(l1Jets_, l1Jets);

  Handle<pat::JetCollection> recoJets;
  iEvent.getByToken(recoJets_, recoJets);

  int iJet = 0;
  for (const auto& recoJet : *recoJets) {
    if ( recoJet.pt() < 250 or std::abs(recoJet.eta()) < 2.5 ) continue;
    // std::cout << "New jet" << std::endl;

    l1Jet_.iJet = iJet;
    l1Jet_.recoJet_pt = recoJet.pt();
    l1Jet_.recoJet_eta = recoJet.eta();
    l1Jet_.recoJet_phi = recoJet.phi();
    l1Jet_.recoJet_energy = recoJet.energy();
    l1Jet_.recoJet_chargedHadronEnergy = recoJet.chargedHadronEnergy();
    l1Jet_.recoJet_neutralHadronEnergy = recoJet.neutralHadronEnergy();
    l1Jet_.recoJet_electronEnergy = recoJet.electronEnergy();
    l1Jet_.recoJet_photonEnergy = recoJet.photonEnergy();
    l1Jet_.recoJet_muonEnergy = recoJet.muonEnergy();
    l1Jet_.recoJet_hfHadEnergy = recoJet.HFHadronEnergy();
    l1Jet_.recoJet_hfEmEnergy = recoJet.HFEMEnergy();

    l1t::Jet matchedJet;
    float matchDr = 999.;
    for (auto l1jet=l1Jets->begin(0); l1jet!=l1Jets->end(0); ++l1jet) {
      if ( reco::deltaR(*l1jet, recoJet) < matchDr ) {
        matchDr = reco::deltaR(*l1jet, recoJet);
        matchedJet = *l1jet;
      }
    }

    // std::cout << "hw at L1: " << matchedJet.hwEta() << "," << matchedJet.hwPhi() << std::endl;
    // hwEta from unpacked jets in GT coordinates
    // Cheap conversion trick: find matching tower
    int ieta{0}, iphi{0};
    for (const auto& hcalTp : *hcalTPs) {
      int ietafix = l1t::CaloTools::gtEta(hcalTp.id().ieta());
      if ( ietafix == matchedJet.hwEta() and l1t::CaloTools::gtPhi(0, hcalTp.id().iphi()) == matchedJet.hwPhi() ) {
        // std::cout << "gt matched tp" << hcalTp.id() << std::endl;
        ieta = hcalTp.id().ieta();
        iphi = hcalTp.id().iphi();
      }
    }
    // Not sure why it doesn't always work
    if ( ieta == 0 ) {
      float matchTowDr = 999.;
      for (const auto& hcalTp : *hcalTPs) {
        float deta = std::fabs(l1t::CaloTools::towerEta(hcalTp.id().ieta()) - matchedJet.eta());
        float dphi = reco::deltaPhi(l1t::CaloTools::towerPhi(0, hcalTp.id().iphi()), matchedJet.phi());
        float dr = std::hypot(deta, dphi);
        if ( dr < matchTowDr ) {
          matchTowDr = dr;
          ieta = hcalTp.id().ieta();
          iphi = hcalTp.id().iphi();
        }
      }
      // std::cout << "gt matched tp ieta " << ieta << " iphi " << iphi << ", dR=" << matchTowDr << std::endl;
      // std::cout << "gt coord " << l1t::CaloTools::gtEta(ieta) << ", " << l1t::CaloTools::gtPhi(0, iphi) << std::endl;
    }

    l1Jet_.l1Jet_matchDeltaR = matchDr;
    l1Jet_.l1Jet_et = matchedJet.et();
    l1Jet_.l1Jet_eta = matchedJet.eta();
    l1Jet_.l1Jet_phi = matchedJet.phi();
    l1Jet_.l1Jet_hwPt = matchedJet.hwPt();
    l1Jet_.l1Jet_hwEta = matchedJet.hwEta();
    l1Jet_.l1Jet_hwPhi = matchedJet.hwPhi();
    l1Jet_.l1Jet_towerIEta = matchedJet.towerIEta();
    l1Jet_.l1Jet_towerIPhi = matchedJet.towerIPhi();
    l1Jet_.l1Jet_rawEt = matchedJet.rawEt();
    l1Jet_.l1Jet_seedEt = matchedJet.seedEt();
    l1Jet_.l1Jet_puEt = matchedJet.puEt();

    l1Jet_.caloTower_hwEt.clear();
    l1Jet_.caloTower_ieta.clear();
    l1Jet_.caloTower_iphi.clear();
    for (const auto& caloTp : *caloTowers) {
      int deta = std::abs(caloTp.hwEta() - ieta);
      if ( deta > 4 ) continue;
      int dphi = std::abs((caloTp.hwPhi() - iphi + 36) % 72 - 36);
      if ( dphi > 4 ) continue;
      l1Jet_.caloTower_hwEt.push_back(caloTp.hwPt());
      l1Jet_.caloTower_ieta.push_back(caloTp.hwEta());
      l1Jet_.caloTower_iphi.push_back(caloTp.hwPhi());
    }

    l1Jet_.hcalTower_et.clear();
    l1Jet_.hcalTower_eta.clear();
    l1Jet_.hcalTower_phi.clear();
    l1Jet_.hcalTower_hwEt.clear();
    l1Jet_.hcalTower_ieta.clear();
    l1Jet_.hcalTower_iphi.clear();
    for (const auto& hcalTp : *hcalTPs) {
      int deta = std::abs(hcalTp.id().ieta() - ieta);
      if ( deta > 4 ) continue;
      int dphi = std::abs((hcalTp.id().iphi() - iphi + 36) % 72 - 36);
      if ( dphi > 4 ) continue;
      //auto cell = caloGeometry_->getGeometry(hcalTp.id());
      l1Jet_.hcalTower_et.push_back(decoder->hcaletValue(hcalTp.id(), hcalTp.t0()));
      //l1Jet_.hcalTower_eta.push_back(cell->etaPos());
      //l1Jet_.hcalTower_phi.push_back(cell->phiPos());
      l1Jet_.hcalTower_hwEt.push_back(hcalTp.SOI_compressedEt());
      l1Jet_.hcalTower_ieta.push_back(hcalTp.id().ieta());
      l1Jet_.hcalTower_iphi.push_back(hcalTp.id().iphi());
    }

    l1Jet_.ecalTower_et.clear();
    l1Jet_.ecalTower_eta.clear();
    l1Jet_.ecalTower_phi.clear();
    l1Jet_.ecalTower_hwEt.clear();
    l1Jet_.ecalTower_ieta.clear();
    l1Jet_.ecalTower_iphi.clear();
    for (const auto& ecalTp : *ecalTPs) {
      int deta = std::abs(ecalTp.id().ieta() - ieta);
      if ( deta > 4 ) continue;
      int dphi = std::abs((ecalTp.id().iphi() - iphi + 36) % 72 - 36);
      if ( dphi > 4 ) continue;
      //auto cell = caloGeometry_->getGeometry(ecalTp.id());
      l1Jet_.ecalTower_et.push_back(ecalTp.compressedEt()*0.5);
      //l1Jet_.ecalTower_eta.push_back(cell->etaPos());
      //l1Jet_.ecalTower_phi.push_back(cell->phiPos());
      l1Jet_.ecalTower_hwEt.push_back(ecalTp.compressedEt());
      l1Jet_.ecalTower_ieta.push_back(ecalTp.id().ieta());
      l1Jet_.ecalTower_iphi.push_back(ecalTp.id().iphi());
    }

    // std::cout << "hcal tower size: " << l1Jet_.hcalTower_et.size() << std::endl;
    // std::cout << "ecal tower size: " << l1Jet_.ecalTower_et.size() << std::endl;
    float jetet = 0.;
    for(auto et : l1Jet_.hcalTower_et) jetet += et;
    for(auto et : l1Jet_.ecalTower_et) jetet += et;
    // std::cout << "jet sum et: " << jetet << std::endl;
    // std::cout << "l1jet sum et: " << matchedJet.et() << std::endl;


    l1JetTree_->Fill();
    iJet++;
  }
}


void 
L1JetTupler::beginJob()
{
}


void 
L1JetTupler::endJob() 
{
}


void
L1JetTupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(L1JetTupler);
