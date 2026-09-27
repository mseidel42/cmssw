// MasslessBWeight.cc
//
// Implementation of the MasslessBWeightHook that reweights trial emissions
// in the FSR shower to soften the b-quark dead cone. The logic mirrors the
// standalone Pythia8 main115.cc example.

#include "MasslessBWeight.h"

#include <algorithm>
#include <cmath>
#include <string>
#include "Pythia8/Logger.h"

MasslessBWeightHook::MasslessBWeightHook(const edm::ParameterSet& iConfig)
    // "bMass" is the nominal b mass (default 4.8 GeV).
    // "bTarget" is the target b mass used in the reweighting weight (default
    // 2.7 GeV, the lowest we can go safely with overSampleFSR = 10, the
    // Pythia maximum and CMS standard for FSR variations).
    : mB(iConfig.exists("bMass") ? iConfig.getParameter<double>("bMass") : 4.8),
      mTgt(iConfig.exists("bTarget") ? iConfig.getParameter<double>("bTarget") : 1.5) {}

double MasslessBWeightHook::enhanceEmission(int /*beamKind*/,
                                            int iRad,
                                            int /*iRec*/,
                                            int iEmt) {
  // Only target b quarks (PDG id = 5).
  if (std::abs(workEvent[iRad].id()) != 5)
    return 1.0;

  // Radiator energy after the splitting.
  double energy = workEvent[iRad].e();
  if (energy <= mB)
    return 1.0;  // avoid unphysical kinematics

  // Splitting angle between the radiator and the emitted parton.
  Pythia8::Vec4 pRad = workEvent[iRad].p();
  Pythia8::Vec4 pEmt = workEvent[iEmt].p();
  double thetaEmt = Pythia8::theta(pRad, pEmt);

  // Protect against perfectly collinear emissions.
  if (thetaEmt < 1e-9)
    return 1.0;

  // Regulated weight to shift the dead cone from mB to mTgt:
  //   W = ((theta^2 + theta0^2) / (theta^2 + thetaTgt^2))^2
  // where theta0 = mB / E (nominal dead-cone angle) and
  // thetaTgt = mTgt / E (target dead-cone angle).
  // For mTgt = 0 this reduces to the full dead-cone removal
  //   W = (1 + (theta0/theta)^2)^2,
  // and for mTgt = mB the weight is 1 (no change).
  double thetaEmt2 = Pythia8::pow2(thetaEmt);
  double thetaZero2 = Pythia8::pow2(mB / energy);
  double thetaTgt2 = Pythia8::pow2(mTgt / energy);
  double weight = Pythia8::pow2((thetaEmt2 + thetaZero2) /
                                (thetaEmt2 + thetaTgt2));

  // Warn if the weight exceeds the over-sampling factor, since then the
  // trial over-sampling is insufficient and the shower efficiency drops.
  double overSampleFSR = settingsPtr->parm("UncertaintyBands:overSampleFSR");
  if (weight > overSampleFSR) {
    loggerPtr->warningMsg("MasslessBWeightHook::enhanceEmission",
      "Calculated enhancement weight " + std::to_string(weight) +
      " exceeds the oversampling factor of " + std::to_string(overSampleFSR) +
      " for this emission. Consider increasing either "
      "UncertaintyBands:overSampleFSR or the target mass.");
  }

  // Cap at 1000 to avoid numerical spikes in the deep collinear limit.
  return std::min(weight, 1000.0);
}
