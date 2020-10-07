from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
config = config()
config.JobType.allowUndistributedCMSSW = True
config.General.requestName = 'DiMuonValidation_jobs_v6_TrackRefitter_92X_dataRun2_Prompt_v11'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../test/DiMuonValidation_cfg.py'
config.Data.inputDataset = '/DoubleMuon/Run2018A-TkAlZMuMu-PromptReco-v3/ALCARECO' #'/DoubleMuon/Run2018A-TkAlZMuMu-PromptReco-v1/ALCARECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20 
config.Data.lumiMask = 'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
config.Data.publication = False
#config.Data.useParent = True

config.Site.storageSite = 'T2_US_Purdue'
config.Data.outLFNDirBase = '/store/user/agu/EPR_Tracker_Alignment/v6_TrackRefitter_92X_dataRun2_Prompt_v11'
