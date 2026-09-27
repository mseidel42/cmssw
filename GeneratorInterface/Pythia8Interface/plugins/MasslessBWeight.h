// MasslessBWeight.h
//
// CMSSW integration of the MasslessBWeightHook developed for the standalone
// Pythia8 main115.cc example that softens the b-quark dead cone.
//
// The b is kept massive in the shower, but the emission probability is
// reweighted so that the dead cone is effectively shifted from the nominal
// b mass (mB) to a smaller target mass (mTgt). The reweighting factor is
//
//     W = ((theta^2 + theta0^2) / (theta^2 + thetaTgt^2))^2
//
// where theta is the emission angle, theta0 = mB / E (nominal dead-cone
// angle) and thetaTgt = mTgt / E (target dead-cone angle). This shifts
// the dead cone from the nominal mass to the target mass, partially
// softening the dead cone for 0 < mTgt < mB (e.g. mTgt = mB/2 gives
// roughly half the dead-cone angle).
//
// To be used together with the Pythia settings
//     UncertaintyBands:doVariations  = on
//     UncertaintyBands:List           = { dummy fsr:muRfac=1.0 }
//     UncertaintyBands:overSampleFSR  = 10.0
// so that the "delayed acceptance" path in SimpleTimeShower::pTnext() is
// active, the trial is over-sampled, and the physical accept probability
// (including the dead cone) is stored in dip.pAccept. The enhanceEmission
// hook then multiplies this by W to give the reweighted accept probability.
//
// The nominal b mass and the target mass are configurable via the
// ParameterSet entries "bMass" (default 4.8 GeV) and "bTarget" (default 1.5
// GeV, the charm quark mass, shifting the dead cone from the b to the c
// scale). With overSampleFSR = 100, the worst-case weight is W ~ 55,
// safely below 100.
//
// The hook warns if the weight exceeds the oversampling factor, since in
// that case the trial over-sampling is insufficient and the shower
// efficiency drops (consider increasing overSampleFSR or raising mTgt).

#ifndef GeneratorInterface_Pythia8Interface_MasslessBWeight_h
#define GeneratorInterface_Pythia8Interface_MasslessBWeight_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "Pythia8/Pythia.h"
#include <memory>
#include "GeneratorInterface/Pythia8Interface/interface/CustomHook.h"

class MasslessBWeightHook : public Pythia8::UserHooks {
public:
  MasslessBWeightHook(const edm::ParameterSet& iConfig);
  ~MasslessBWeightHook() override {}

  // Tell Pythia we want to enhance trial emissions.
  bool canEnhanceEmission() override { return true; }
  // The core reweighting logic. iRad, iRec, iEmt are indices in workEvent.
  double enhanceEmission(int beamKind, int iRad, int iRec, int iEmt) override;

private:
  double mB;    // nominal b mass (e.g. 4.8 GeV)
  double mTgt; // target b mass (e.g. 2.7 GeV, half nominal)
};

REGISTER_USERHOOK(MasslessBWeightHook);
#endif
