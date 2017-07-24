import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.parseArguments()

from Configuration.StandardSequences.Eras import eras
process = cms.Process('SKIM', eras.Run2_2016)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(False),
    wantSummary = cms.untracked.bool(True),
)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        options.inputFiles
    ),
    secondaryFileNames = cms.untracked.vstring()
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v16', '')


process.load("Configuration.StandardSequences.RawToDigi_Data_cff")    
process.load("EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi")

process.skimPath = cms.Path(
    process.RawToDigi +
    process.l1tCaloLayer1Digis
)

process.skimOutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dropMetaData = cms.untracked.string('NONE'),
    fastCloning = cms.untracked.bool(True),
    fileName = cms.untracked.string(options.outputFile),
    overrideInputFileSplitLevels = cms.untracked.bool(False),
    outputCommands = cms.untracked.vstring(
        'keep *_slimmedJets_*_*',
        'keep *_l1tCaloLayer1Digis_*_*',
        'keep *_caloStage2Digis_Jet_*',
        'keep  l1tJetBXVector_simCaloStage2Digis_*_*',
    ),
)
process.output_step = cms.EndPath(process.skimOutput)

process.schedule = cms.Schedule(
    process.skimPath,
    process.output_step,
)

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAW 

#call to customisation function L1TReEmulFromRAW imported from L1Trigger.Configuration.customiseReEmul
process = L1TReEmulFromRAW(process)

#process.remove(process.TotemDAQMappingESSourceXML)
#ctppsDiamondRawToDigi
#totemRPRawToDigi
#totemTriggerRawToDigi
