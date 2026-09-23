# Pythia8DeadConeHooks_cfi.py
#
# Example configuration showing how to switch off or soften the b-quark dead
# cone when generating CMS samples with the Pythia8Interface, in particular
# when showering Powheg (or other) LHE files that contain MASSIVE b quarks.
#
# Three user hooks are provided (all registered with the CustomHookFactory of
# the GeneratorInterface/Pythia8Interface package and attached via the
# "UserCustomization" parameter of the Pythia8HadronizerFilter /
# Pythia8HepMC3HadronizerFilter EDFilter):
#
#   * MasslessLHEInputHook : runs at the process level (after the LHE is read,
#                            before the shower) and makes every final-state b
#                            quark from the LHE massless. For b quarks that come
#                            from a resonance decay (e.g. t -> b W in a Powheg
#                            LHE file with full spin-correlated top decays), the
#                            hook rescales ALL the resonance daughters by a
#                            common factor alpha in the resonance rest frame,
#                            chosen so that (i) the b is massless and (ii) the
#                            total energy of the daughters equals the resonance
#                            mass. This preserves the resonance (top) mass
#                            EXACTLY. For two-body t -> b W the formula reduces
#                            to the standard two-body kinematics with m_b = 0.
#                            For b quarks not from a resonance (e.g. b-initiated
#                            hard processes, b from g -> bb in the ME), the
#                            b's 3-momentum is preserved and its energy is set
#                            to |p|.
#
#                            This is needed because the Pythia setting
#                                5:m0 = 0.0
#                            only changes the *default* b mass used for
#                            NEWLY CREATED b quarks (e.g. g -> bb in the
#                            shower); it does NOT override the mass of b quarks
#                            that are already in the event record from the
#                            LHE.
#
#   * MassRestorationHook   : runs at the end of the parton level (after the
#                            shower, before hadronization) and restores the
#                            b mass to a configurable value (default 4.8 GeV,
#                            the typical B-hadron mass) by finding the nearest
#                            gluon color partner and putting the (b, g) pair on
#                            their new mass shells while preserving the pair
#                            3-momentum. If no kinematically allowed partner can
#                            be found the event is vetoed and Pythia re-tries
#                            the parton level.
#
#                            The Pythia setting "5:m0 = 0.0" must also be set
#                            so the shower actually runs the b quarks as
#                            massless.
#
#   * MasslessBWeightHook    : softens (rather than removes) the dead cone by
#                            reweighting the trial emission probability so that
#                            emissions inside the dead cone are enhanced. The b
#                            is kept massive in the shower. The Pythia settings
#                                Enhancements:doEnhanceTrial = on
#                                Enhancements:overSampleFSR  = 10
#                            must also be set so the over-sampling actually
#                            generates the enhanced emissions.
#
# The b mass used by MassRestorationHook and MasslessBWeightHook is
# configurable via the "bMass" entry of the UserCustomization PSet and defaults
# to 4.8 GeV.
#
# -----------------------------------------------------------------------------
# Recommended combinations for showering MASSIVE-LHE b quarks (e.g. Powheg LHE
# files generated with bmass = 4.8):
#
#   * Option 6 -- massless b in shower, hadronization produces the B mass:
#
#       Pythia settings: 5:m0 = 0.0
#       UserCustomization:
#           pluginName = "MasslessLHEInputHook"
#
#   * Option 5 -- massless b in shower, mass explicitly restored before
#                hadronization:
#
#       Pythia settings: 5:m0 = 0.0
#       UserCustomization (both hooks):
#           pluginName = "MasslessLHEInputHook"
#           pluginName = "MassRestorationHook"
#               bMass = 4.8
#
#   * Option 7 -- soften (do not remove) the dead cone by trial-probability
#                reweighting, keeping the b massive:
#
#       Pythia settings: Enhancements:doEnhanceTrial = on
#                        Enhancements:overSampleFSR  = 10
#       UserCustomization:
#           pluginName = "MasslessBWeightHook"
#               bMass = 4.8
#
# -----------------------------------------------------------------------------
# Note for showering MASSLESS-LHE b quarks (e.g. a 5FS Powheg LHE generated
# with bmass = 0, or a Pythia-internal hard process):
# In that case the LHE b quarks are already massless, so MasslessLHEInputHook
# is a no-op and can be omitted. Just use "5:m0 = 0.0" for Options 5 and 6,
# or the Enhancements:* settings for Option 7.

import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------------------------
# Pythia8 settings blocks for the three dead-cone methods.
# ---------------------------------------------------------------------------

# Settings common to all methods (CP5 tune as an example; replace with the
# tune appropriate for your sample).
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

# Method 1 (Option 6): massless b in the shower. The B-hadron mass is then
# produced by the string fragmentation. To be used together with the
# MasslessLHEInputHook (for massive-LHE input).
pythia8DeadConeMassRestorationSettings = cms.vstring(
    '5:m0 = 0.0',
)

# Method 2 (Option 5): same Pythia settings as Method 1, plus the
# MassRestorationHook to explicitly restore the b mass before hadronization.
# (Same settings block is reused.)

# Method 3 (Option 7): soften the dead cone by reweighting the trial emission
# probability. The b is kept massive in the shower. The enhancement
# machinery must be switched on so the over-sampling actually produces the
# enhanced emissions.
pythia8DeadConeMasslessBWeightSettings = cms.vstring(
    'Enhancements:doEnhanceTrial = on',
    'Enhancements:overSampleFSR  = 10',
)

# ---------------------------------------------------------------------------
# Example generator block (Pythia8HadronizerFilter).
#
# This is an example for an internal-Pythia8 ttbar sample. Adapt the process,
# beam, energy, and tune to your needs. Pick the "UserCustomization" entry
# that corresponds to the method you want to use.
# ---------------------------------------------------------------------------
exampleDeadConeGenerator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring(
            'Top:all = on',
        ),
        parameterSets = cms.vstring(
            'pythia8CommonSettings',
            'pythia8CP5Settings',
            'processParameters',
        )
    ),

    # -------------------------------------------------------------------
    # Pick ONE of the UserCustomization blocks below.
    # Remember to also add the corresponding Pythia settings block
    # (pythia8DeadConeMassRestorationSettings for Options 5/6,
    #  pythia8DeadConeMasslessBWeightSettings for Option 7) to the
    # parameterSets above.
    # -------------------------------------------------------------------

    # Option 6: massless b in shower, B mass produced by hadronization.
    UserCustomization = cms.VPSet(
        cms.PSet(
            pluginName = cms.string("MasslessLHEInputHook")
        )
    ),

    # Option 5: massless b in shower + explicit mass restoration before
    # hadronization. Both hooks must be attached.
    # UserCustomization = cms.VPSet(
    #     cms.PSet(
    #         pluginName = cms.string("MasslessLHEInputHook")
    #     ),
    #     cms.PSet(
    #         pluginName = cms.string("MassRestorationHook"),
    #         bMass = cms.double(4.8)
    #     )
    # ),

    # Option 7: soften the dead cone by trial-probability reweighting. The b
    # is kept massive; MasslessLHEInputHook is NOT used here.
    # UserCustomization = cms.VPSet(
    #     cms.PSet(
    #         pluginName = cms.string("MasslessBWeightHook"),
    #         bMass = cms.double(4.8)
    #     )
    # ),
)
