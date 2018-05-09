#!/usr/bin/env python




import os,glob,time

versionDir = '/pnfs/iihe/cms/store/user/tomc/tnp/muons/Moriond18_v2/'

samples2016 = {
  #         'DY_LO'   : 'mc/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYToLL_madgraph/*/*/*.root',
  #         'runB'    : 'data/SingleMuon/crab_Run2016B/*/*/*.root',
  #         'runC'    : 'data/SingleMuon/crab_Run2016C/*/*/*.root',
  #         'runD'    : 'data/SingleMuon/crab_Run2016D/*/*/*.root',
  #         'runE'    : 'data/SingleMuon/crab_Run2016E/*/*/*.root',
  #         'runF'    : 'data/SingleMuon/crab_Run2016F/*/*/*.root',
  #         'runG'    : 'data/SingleMuon/crab_Run2016G/*/*/*.root',
  #         'runH-v2' : 'data/SingleMuon/crab_Run2016H-v2/*/*/*.root',
  #         'runH-v3' : 'data/SingleMuon/crab_Run2016H-v3/*/*/*.root'

}

samples2017 = {'DY'      : 'mc/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_amcatnlo/*/*/*.root',
  }


cut  = "tag_IsoMu24==1 && tag_combRelIsoPF04dBeta<0.2 && tag_pt>25 && pair_probeMultiplicity>0.5 && pair_probeMultiplicity<1.5 && pt>15 && abseta<2.4"
sel  = "mass pt eta abseta SIP logdxy logdz sip3d combRelIsoPF03 segmentCompatibility JetNDauCharged miniIsoCharged miniIsoNeutral relIso Loose Medium JetPtRel JetPtRatio JetBTagCSV JetDeepBTagCSV tkSigmaPtOverPt tag_nVertices fixedGridRhoFastjetCentralNeutral pair_nJets30"


def launch(command, logfile):
    os.system("qsub -v dir=" + os.getcwd() + ",command=\"" + command + "\" -q localgrid@cream02 -o " + logfile + " -e " + logfile + " -l walltime=00:30:00 $CMSSW_BASE/src/MuonAnalysis/TagAndProbe/test/runOnCream02.sh &> .qsub.log")
    with open('.qsub.log','r') as qsublog:
      for l in qsublog:
        if 'Invalid credential' in l:
          time.sleep(10)
          launch(command, logfile)


for samples, year in [(samples2016, '2016'), (samples2017, '2017')]:
  for sample, path in samples.iteritems():
    print versionDir + '/' + path
    inputFiles = glob.glob(versionDir + '/' + path)
    for inputFile in inputFiles:
      outFile = '/user/tomc/public/tagAndProbe/' + year + '/' + sample + '/' + inputFile.split('/')[-1]
      outTag  = outFile.replace('.root','')
      logFile = outFile.replace('.root','.log')
      try:    os.makedirs(os.path.dirname(logFile))
      except: pass
      command  = './skimTree ' + inputFile + ' ' + outFile + ' -d tpTree -t fitter_tree -c \'' + cut + '\' -r \'all\' -k \'' + sel + '\';'
      command  = '\\\"' + command + '\\\"'
      print inputFile
      launch(command, logFile)
    time.sleep(300)
