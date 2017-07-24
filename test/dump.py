#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
from DataFormats.FWLite import Handle, Runs, Lumis, Events
import sys
import math

ecalDigiH = Handle('edm::SortedCollection<EcalTriggerPrimitiveDigi,edm::StrictWeakOrdering<EcalTriggerPrimitiveDigi> >')
hcalDigiH = Handle('edm::SortedCollection<HcalTriggerPrimitiveDigi,edm::StrictWeakOrdering<HcalTriggerPrimitiveDigi> >')
caloTowerH = Handle('BXVector<l1t::CaloTower>')

events = Events(["file:L1TEST_RAW2DIGI.root"])

for iev, event in enumerate(events):
    print "\nEvent %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
    event.getByLabel('ecalDigis:EcalTriggerPrimitives', ecalDigiH)
    event.getByLabel('hcalDigis', hcalDigiH)
    event.getByLabel('caloStage2Digis:CaloTower', caloTowerH)

    ecalTPs = ecalDigiH.product()
    hcalTPs = hcalDigiH.product()
    caloTPs = caloTowerH.product()

    print "ecal size", ecalTPs.size()
    print "hcal size", hcalTPs.size()
    print "calo size", caloTPs.size()

    for sat in filter(lambda tp: tp.compressedEt()==255, ecalTPs):
        print "Saturated ecal tp @ (% 3d, %3d)" % (sat.id().ieta(), sat.id().iphi())

    for sat in filter(lambda tp: tp.SOI_compressedEt()==255, hcalTPs):
        print "Saturated hcal tp @ (% 3d, %3d)" % (sat.id().ieta(), sat.id().iphi())

    for sat in filter(lambda tp: tp.hwPt()>=509, caloTPs):
        print "Saturated calo tp @ (% 3d, %3d), code = %3d" % (sat.hwEta(), sat.hwPhi(), sat.hwPt())

    for tp in hcalTPs:
        if tp.id().ieta() == 18 and tp.id().iphi() == 8:
            print tp.SOI_compressedEt()
