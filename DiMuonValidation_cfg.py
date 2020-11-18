import FWCore.ParameterSet.Config as cms
#process = cms.Process("TkAlignmentDiMuonValidation")
process = cms.Process("AlcarecoAnalysis")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
import CalibTracker.Configuration.Common.PoolDBESSource_cfi  
#process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v10'#92X_dataRun2_Prompt_v11 110X_dataRun2_v13
process.GlobalTag.globaltag = '92X_dataRun2_Prompt_v11'
"""
process.GlobalTag.toGet = cms.VPSet(
   cms.PSet(record = cms.string('TrackerAlignmentRcd'),
            tag = cms.string('TrackerAlignment_MC2018_v1'),
            connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
            ),
)
"""
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
               'file:FC020B1E-E8B4-E811-9A24-FA163E5A0019.root'
                )
                            )
###############################################################################
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
import RecoTracker.TrackProducer.TrackRefitters_cff
process.TrackRefitter1 = process.TrackRefitterP5.clone(
    src =  'ALCARECOTkAlZMuMu', #'AliMomConstraint1',
    TrajectoryInEvent = True,
    TTRHBuilder = "WithAngleAndTemplate",
    NavigationSchool = "",
    #constraint = 'momentum', ### SPECIFIC FOR CRUZET
    #srcConstr='AliMomConstraint1' ### SPECIFIC FOR CRUZET$works only with tag V02-10-02 TrackingTools/PatternTools / or CMSSW >=31X
    )
###############################################################################

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.myanalysis = cms.EDAnalyzer('DiMuonValidation',
                                     TkTag = cms.string ('TrackRefitter1'),
                                     #TkTag = cms.string ('ALCARECOTkAlZMuMu'),#TrackRefitter1
                                     Pair_mass_min = cms.double(80),
                                     Pair_mass_max = cms.double(120),
                                     Pair_eta_min = cms.double(-2.5),
                                     Pair_eta_max = cms.double(+2.5),
                                     Pair_pt_min = cms.double(0),
                                     Pair_pt_max = cms.double(1000),
                                     #Variable_CosThetaCS_nbins = cms.int(20)
                                     Variable_CosThetaCS_xmin = cms.double(-1),
                                     Variable_CosThetaCS_xmax = cms.double(+1),
                                     #Variable_PairPt_nbins= cms.int(100),
                                     Variable_PairPt_xmin = cms.double(0),
                                     Variable_PairPt_xmax = cms.double(100),
                                     Variable_DeltaEta_xmin = cms.double(-4.8),
                                     Variable_DeltaEta_xmax = cms.double(+4.8),
                                     Variable_EtaPlus_xmin = cms.double(-2.4),
                                     Variable_EtaPlus_xmax = cms.double(+2.4),
                                     Variable_EtaMinus_xmin = cms.double(-2.4),
                                     Variable_EtaMinus_xmax = cms.double(+2.4)
                              )
process.TFileService = cms.Service("TFileService",
      fileName = cms.string("DiMuonValidation_output.root"),
      #closeFileFast = cms.untracked.bool(True)
  )
process.p = cms.Path(process.offlineBeamSpot*process.TrackRefitter1*process.myanalysis)

