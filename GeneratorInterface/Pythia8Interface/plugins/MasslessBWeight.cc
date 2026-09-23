// MasslessBWeight.cc
//
// Implementation of the MasslessBWeightHook that reweights trial emissions
// in the FSR shower to soften the b-quark dead cone. The logic mirrors the
// standalone Pythia8 main115.cc example.

#include "MasslessBWeight.h"

#include <algorithm>
#include <cmath>

MasslessBWeightHook::MasslessBWeightHook(const edm::ParameterSet& iConfig)
    // "bMass" is optional in the UserCustomization PSet; default to the
    // typical B-hadron mass of 4.8 GeV.
    : mB(iConfig.exists("bMass") ? iConfig.getParameter<double>("bMass") : 4.8) {}

double MasslessBWeightHook::enhanceEmission(int /*beamKind*/,
                                             int iRad,
                                             int /*iRec*/,
                                             int iEmt) {
  // Only target b quarks (PDG id = 5).
  if (std::abs(workEvent[iRad].id()) != 5)
    return 1.0;

  // Dead-cone angle theta0 = m_b / E, using the radiator energy after the
  // splitting.
  double energy = workEvent[iRad].e();
  if (energy <= mB)
    return 1.0;  // avoid unphysical kinematics
  double theta0 = mB / energy;

  // Splitting angle between the radiator and the emitted parton.
  Pythia8::Vec4 pRad = workEvent[iRad].p();
  Pythia8::Vec4 pEmt = workEvent[iEmt].p();
  double thetaEmt = Pythia8::theta(pRad, pEmt);

  // Protect against perfectly collinear emissions.
  if (thetaEmt < 1e-9)
    return 1.0;

  // W = (1 + (theta0/theta)^2)^2, capped at 1000 to avoid numerical spikes
  // in the deep collinear / infrared limit.
  double ratioSq = Pythia8::pow2(theta0 / thetaEmt);
  double weight = Pythia8::pow2(1.0 + ratioSq);

  return std::min(weight, 1000.0);
}
