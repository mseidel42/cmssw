// MassRestoration.cc
//
// Implementation of the MassRestorationHook that restores the b-quark mass
// after a massless-b shower, just before hadronization. The logic mirrors
// the standalone Pythia8 main115.cc example.

#include "MassRestoration.h"

#include <cmath>

MassRestorationHook::MassRestorationHook(const edm::ParameterSet& iConfig)
    // "bMass" is optional in the UserCustomization PSet; default to the
    // typical B-hadron mass of 4.8 GeV.
    : bMass_(iConfig.exists("bMass") ? iConfig.getParameter<double>("bMass") : 4.8) {}

bool MassRestorationHook::doVetoPartonLevel(const Pythia8::Event& constEvent) {
  // Cast away constness to modify the event record before hadronization.
  Pythia8::Event& event = const_cast<Pythia8::Event&>(constEvent);

  for (int i = 0; i < event.size(); ++i) {
    if (std::abs(event[i].id()) == 5 && event[i].isFinal()) {
      int colPartner = -1;
      double minAngle = 1e9;

      // Find the nearest gluon (smallest opening angle) that is
      // kinematically allowed, i.e. s > (m_b + m_g)^2.
      for (int j = 0; j < event.size(); ++j) {
        if (i == j || !event[j].isFinal() || event[j].id() != 21)
          continue;

        double s = (event[i].p() + event[j].p()).m2Calc();
        if (s > Pythia8::pow2(bMass_ + event[j].m())) {
          double angle = Pythia8::theta(event[i].p(), event[j].p());
          if (angle < minAngle) {
            minAngle = angle;
            colPartner = j;
          }
        }
      }

      if (colPartner != -1) {
        Pythia8::Vec4 pB = event[i].p();
        Pythia8::Vec4 pP = event[colPartner].p();
        double mB = bMass_;
        double mP = event[colPartner].m();

        Pythia8::Vec4 pSum = pB + pP;
        double s = pSum.m2Calc();

        if (s > Pythia8::pow2(mB + mP)) {
          // 3-momentum magnitude of either parton in the pair CM frame
          // for the new (on-shell) masses.
          double pCM =
              std::sqrt((s - Pythia8::pow2(mB + mP)) * (s - Pythia8::pow2(mB - mP)) / (2.0 * std::sqrt(s)));

          Pythia8::RotBstMatrix rotBst;
          rotBst.toCMframe(pB, pP);
          Pythia8::Vec4 pB_CM = pB;
          pB_CM.rotbst(rotBst);

          double pOld = pB_CM.pAbs();
          if (pOld > 0) {
            // Rescale to the new on-shell momenta and energies.
            pB_CM *= (pCM / pOld);
            pB_CM.e(std::sqrt(Pythia8::pow2(pCM) + Pythia8::pow2(mB)));

            Pythia8::Vec4 pP_CM = pP;
            pP_CM.rotbst(rotBst);
            pP_CM *= (pCM / pOld);
            pP_CM.e(std::sqrt(Pythia8::pow2(pCM) + Pythia8::pow2(mP)));

            // Boost back to the lab frame.
            rotBst.invert();
            pB_CM.rotbst(rotBst);
            pP_CM.rotbst(rotBst);

            event[i].p(pB_CM);
            event[i].m(mB);
            event[colPartner].p(pP_CM);
            event[colPartner].m(mP);
          }
        } else {
          // The massless shower populated a kinematically forbidden region
          // (s < (m_b + m_g)^2). Veto so Pythia retries the parton level.
          return true;
        }
      }
    }
  }
  // Let Pythia proceed to string fragmentation with the modified event.
  return false;
}
