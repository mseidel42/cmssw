import FWCore.ParameterSet.Config as cms

#
# produce hitFit hypothesis with all necessary
# ingredients
#

## std sequence to perform kinematic fit
import TopQuarkAnalysis.TopHitFit.hitFitTtSemiLepEventMuons_cfi
hitFitTtSemiLepEventHypothesis = TopQuarkAnalysis.TopHitFit.hitFitTtSemiLepEventMuons_cfi.hitFitTtSemiLepEventMuons.clone()

## configure hitFit hypothesis
from TopQuarkAnalysis.TopJetCombination.TtSemiLepHypHitFit_cfi import *

## make hypothesis
makeHypothesis_hitFitTask = cms.Task(
  hitFitTtSemiLepEventHypothesis,
  ttSemiLepHypHitFit
)
makeHypothesis_hitFit = cms.Sequence(makeHypothesis_hitFitTask)
