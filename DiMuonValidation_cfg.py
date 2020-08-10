import FWCore.ParameterSet.Config as cms
process = cms.Process("TkAlignmentDiMuonValidation")
process.load("FWCore.MessageService.MessageLogger_cfi")

import CalibTracker.Configuration.Common.PoolDBESSource_cfi  


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
               'file:FC020B1E-E8B4-E811-9A24-FA163E5A0019.root'
                )
                            )
process.TFileService = cms.Service("TFileService",
      fileName = cms.string("DiMuonValidation_output.root"),
      #closeFileFast = cms.untracked.bool(True)
  )

process.myanalysis = cms.EDAnalyzer('DiMuonValidation',
                                     TkTag = cms.string ('ALCARECOTkAlZMuMu'),#TrackRefitter1
                              )

process.p = cms.Path(process.myanalysis)

