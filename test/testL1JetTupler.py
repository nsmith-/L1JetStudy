import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Ntupler", eras.Run2_2016)
options = VarParsing('analysis')
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring(options.inputFiles) )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v16', '')

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.HcalTrigTowerGeometryESProducer = cms.ESProducer("HcalTrigTowerGeometryESProducer")

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)

process.load("Configuration.StandardSequences.RawToDigi_Data_cff")    
process.load("EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi")

process.ntupler = cms.EDAnalyzer("L1JetTupler",
    ecalTPs = cms.InputTag("l1tCaloLayer1Digis"),
    hcalTPs = cms.InputTag("l1tCaloLayer1Digis"),
    # l1Jets = cms.InputTag("caloStage2Digis:Jet"),
    l1Jets = cms.InputTag("simCaloStage2Digis"),
    patJets = cms.InputTag("slimmedJets"),
)

process.p = cms.Path(
    #process.l1tCaloLayer1Digis +
    process.ntupler
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

