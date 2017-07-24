// -*- C++ -*-
//
// Package:    Analysis/L1Saturation
// Class:      L1Saturation
// 
/**\class L1Saturation L1Saturation.cc Analysis/L1Saturation/plugins/L1Saturation.cc

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

class L1Saturation : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>  {
  public:
    explicit L1Saturation(const edm::ParameterSet&);
    ~L1Saturation();

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

    TTree * l1JetTree_;
    L1JetStruct l1Jet_;
    edm::ESHandle<CaloGeometry> caloGeometry_;
};

L1Saturation::L1Saturation(const edm::ParameterSet& iConfig):
  ecalTPs_(consumes<EcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("ecalTPs"))),
  hcalTPs_(consumes<HcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("hcalTPs"))),
  caloTowers_(consumes<l1t::CaloTowerBxCollection>(iConfig.getParameter<edm::InputTag>("caloTowers")))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  l1JetTree_ = fs->make<TTree>("l1Jet", "All gen electrons, nearest reco");
  l1JetTree_->Branch("run", &l1Jet_.run);
  l1JetTree_->Branch("lumi", &l1Jet_.lumi);
  l1JetTree_->Branch("event", &l1Jet_.event);
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


L1Saturation::~L1Saturation()
{
}

void
L1Saturation::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  iSetup.get<CaloGeometryRecord>().get(caloGeometry_);
}

void
L1Saturation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  l1Jet_.caloTower_hwEt.clear();
  l1Jet_.caloTower_ieta.clear();
  l1Jet_.caloTower_iphi.clear();
  for (const auto& caloTp : *caloTowers) {
    if ( caloTp.hwPt() <= 508 ) continue;
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
    if ( hcalTp.SOI_compressedEt() < 255 ) continue;
    double physEt = decoder->hcaletValue(hcalTp.id(), hcalTp.t0());
    if ( physEt < 127.5 ) continue;
    l1Jet_.hcalTower_et.push_back(physEt);
    //auto cell = caloGeometry_->getGeometry(hcalTp.id());
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
    if ( ecalTp.compressedEt() < 255 ) continue;
    //auto cell = caloGeometry_->getGeometry(ecalTp.id());
    l1Jet_.ecalTower_et.push_back(ecalTp.compressedEt()*0.5);
    //l1Jet_.ecalTower_eta.push_back(cell->etaPos());
    //l1Jet_.ecalTower_phi.push_back(cell->phiPos());
    l1Jet_.ecalTower_hwEt.push_back(ecalTp.compressedEt());
    l1Jet_.ecalTower_ieta.push_back(ecalTp.id().ieta());
    l1Jet_.ecalTower_iphi.push_back(ecalTp.id().iphi());
  }

  l1JetTree_->Fill();
}


void 
L1Saturation::beginJob()
{
}


void 
L1Saturation::endJob() 
{
}


void
L1Saturation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(L1Saturation);
