import FWCore.ParameterSet.Config as cms

hitFitTtSemiLepEvent = cms.EDProducer("TtSemiLepHitFitProducerElectron",
    jets = cms.InputTag("selectedPatJets"),
    leps = cms.InputTag("selectedPatElectrons"),
    mets = cms.InputTag("patMETs"),

    # ------------------------------------------------
    # maximum number of jets to be considered in the
    # jet combinatorics (has to be >= 4, can be set to
    # -1 if you want to take all)
    # ------------------------------------------------
    maxNJets = cms.int32(4),

    #-------------------------------------------------
    # maximum number of jet combinations finally
    # written into the event, starting from the "best"
    # (has to be >= 1, can be set to -1 if you want to 
    # take all)
    #-------------------------------------------------
    maxNComb = cms.int32(1),
    
    # ------------------------------------------------
    # option to use b-tagging
    # ------------------------------------------------
    bTagAlgo          = cms.string("trackCountingHighEffBJetTags"),
    minBDiscBJets     = cms.double(0.0),
    maxBDiscLightJets = cms.double(1.0),
    useBTagging       = cms.bool(False),
    useFiveJets       = cms.bool(False),
    useSixJets        = cms.bool(False),
    useBTagEmulation  = cms.bool(False),
    
    # ------------------------------------------------
    # set mass values used in the constraints
    # set mass to 0 for no constraint
    # ------------------------------------------------    
    mW   = cms.double(80.4),
    mTop = cms.double(0.),
    
    # ------------------------------------------------
    # specify jet correction level as, Uncorrected, L1Offset, L2Relative, L3Absolute, L4Emf,
    # L5Hadron, L6UE, L7Parton. This is currently only used for debugging. The jet collection
    # should be "ready" when it arrives here.
    # ------------------------------------------------
    jetCorrectionLevel = cms.string("L3Absolute"),

    # ------------------------------------------------
    # check that the jet collection behaves as expected, including e.g. the JEC level.
    # ------------------------------------------------
    runJetTests = cms.bool(True),
)
