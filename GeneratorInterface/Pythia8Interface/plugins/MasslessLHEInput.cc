// MasslessLHEInput.cc
//
// Implementation of the MasslessLHEInputHook that makes the final-state
// b quarks coming from the LHE massless before the shower runs, while
// preserving the mass of any resonance (top, Z, W, Higgs, ...) that the b
// came from.
//
// Two cases are handled:
//  (a) b from a resonance decay (e.g. t -> b W): All of the resonance's
//      direct daughters are rescaled by a common factor alpha in the
//      resonance rest frame, so that (i) every final-state b/bbar daughter
//      becomes massless and (ii) the total energy of the daughters equals
//      the resonance mass. This preserves the resonance mass exactly.
//      Additionally, for every rescaled daughter that itself has descendants
//      (e.g. the W from t -> b W), the daughter's 4-momentum change is
//      propagated to ALL of the daughter's descendants (recursively)
//      using the Lorentz boost that maps the daughter's original
//      4-momentum to its new 4-momentum. This preserves momentum
//      conservation at every level of the decay chain.
//
//      Without this descendant propagation, the resonance's resonance
//      daughter (e.g. the W) would have a new 4-momentum that no longer
//      equals the sum of its (unchanged) descendants' 4-momenta, breaking
//      momentum conservation in the W decay. The broken kinematics can
//      later produce NaN values in PDF/pT evaluations, which was the
//      cause of the "Unphysical x given: -nan" error in the Powheg-matched
//      shower.
//
//  (b) b that is not from a resonance decay (e.g. a b-initiated hard process
//      outgoing parton, or b from g -> bb in the ME): The b's 3-momentum is
//      preserved, its energy is set to |p|, and its stored mass is set to 0.
//      (No resonance mass to preserve in this case; the small energy
//      loss ~ m_b^2 / (2 |p|) is absorbed into the rest of the event.)
//
// Incoming b quarks (e.g. from a b-initiated hard process in the 5FS) are NOT
// modified -- changing their 4-momentum would distort the hard ME PDF
// kinematics.

#include "MasslessLHEInput.h"

#include <cmath>
#include <vector>

namespace {

// Standard model "resonance" PDG IDs whose mass we want to preserve when
// their decay produces a b quark. Absolute values are matched.
bool isResonanceId(int idAbs) {
  return (idAbs == 6 ||   // top
          idAbs == 23 ||  // Z0
          idAbs == 24 ||  // W+/W-
          idAbs == 25 ||  // Higgs
          idAbs == 35 ||  // H0
          idAbs == 36 ||  // A0
          idAbs == 37 ||  // H+
          idAbs == 32 ||  // Z'
          idAbs == 33);    // W'
}

// Walk up the mother chain from index iB, return the index of the first
// resonance ancestor, or -1 if none found.
int findResonanceAncestor(const Pythia8::Event& event, int iB) {
  int i = iB;
  for (int step = 0; step < event.size() && i > 0; ++step) {
    int m1 = event[i].mother1();
    int m2 = event[i].mother2();
    int iMother = (m1 > 0) ? m1 : m2;
    if (iMother <= 0 || iMother >= event.size())
      break;
    if (isResonanceId(std::abs(event[iMother].id()))) {
      return iMother;
    }
    i = iMother;
  }
  return -1;
}

// Recursively apply a Lorentz transformation T to all descendants of
// particle i (modifying their 4-momenta). T must be a Lorentz transformation,
// so the "resonance = sum of descendants" relations are preserved at every
// level of the decay tree.
void transformDescendants(Pythia8::Event& event, int i,
                           const Pythia8::RotBstMatrix& T) {
  std::vector<int> daughters = event[i].daughterList();
  for (int iD : daughters) {
    if (iD <= 0 || iD >= event.size())
      continue;
    Pythia8::Vec4 p = event[iD].p();
    p.rotbst(T);
    event[iD].p(p);
    transformDescendants(event, iD, T);
  }
}

}  // namespace

MasslessLHEInputHook::MasslessLHEInputHook(const edm::ParameterSet& iConfig) {}

bool MasslessLHEInputHook::doVetoProcessLevel(Pythia8::Event& event) {
  // We may visit the same resonance more than once if it has two b daughters
  // (e.g. H -> b bbar). Track which resonance indices have already been
  // rescaled to avoid double-processing.
  std::vector<int> alreadyDone;

  for (int i = 0; i < event.size(); ++i) {
    if (std::abs(event[i].id()) != 5 || !event[i].isFinal())
      continue;

    int iRes = findResonanceAncestor(event, i);
    if (iRes < 0) {
      // Case (b): b not from a resonance. Just make it massless, preserving
      // the 3-momentum.
      double pAbs = event[i].pAbs();
      event[i].e(pAbs);
      event[i].m(0.0);
      continue;
    }

    // Check we haven't already rescaled this resonance.
    bool skip = false;
    for (int idx : alreadyDone) {
      if (idx == iRes) {
        skip = true;
        break;
      }
    }
    if (skip)
      continue;
    alreadyDone.push_back(iRes);

    // Case (a): b from a resonance. Rescale all the resonance daughters by a
    // common factor alpha in the resonance rest frame.
    Pythia8::Vec4 pRes = event[iRes].p();
    double mRes = pRes.mCalc();
    if (mRes <= 0.)
      continue;

    std::vector<int> daughters = event[iRes].daughterList();
    if (daughters.empty())
      continue;

    Pythia8::RotBstMatrix toRes;
    toRes.bstback(pRes);

    // Pre-compute, per daughter, the 3-momentum magnitude and the mass, in
    // the resonance rest frame.
    struct DaughterKin {
      int idx;
      double pMag;  // 3-momentum magnitude in the resonance rest frame
      double mass;  // 0 for b's being massless-ized, original mass otherwise
      bool isB;    // is this a b/bbar we want to make massless?
    };
    std::vector<DaughterKin> dks;
    for (int iD : daughters) {
      if (iD <= 0 || iD >= event.size())
        continue;
      Pythia8::Vec4 pD = event[iD].p();
      pD.rotbst(toRes);
      DaughterKin dk;
      dk.idx = iD;
      dk.pMag = pD.pAbs();
      dk.isB = (std::abs(event[iD].id()) == 5 && event[iD].isFinal());
      dk.mass = dk.isB ? 0.0 : event[iD].m();
      dks.push_back(dk);
    }
    if (dks.empty())
      continue;

    // f(alpha) = sum_d E_d'(alpha) - mRes, where for b daughters
    // E_d' = alpha * |p_d| (massless) and for other daughters
    // E_d' = sqrt(alpha^2 * p_d^2 + m_d^2) (on-shell with original mass).
    // f is monotonically increasing in alpha, so a unique root exists and
    // can be found by bisection.
    auto totalEnergy = [&](double alpha) -> double {
      double E = 0.0;
      for (const auto& dk : dks) {
        double pAlpha = alpha * dk.pMag;
        E += std::sqrt(pAlpha * pAlpha + dk.mass * dk.mass);
      }
      return E;
    };
    auto f = [&](double alpha) { return totalEnergy(alpha) - mRes; };

    double aLo = 0.0, aHi = 1.0;
    while (f(aHi) < 0.0 && aHi < 1e6) {
      aHi *= 2.0;
    }
    for (int iter = 0; iter < 200; ++iter) {
      double mid = 0.5 * (aLo + aHi);
      if (f(mid) < 0.0)
        aLo = mid;
      else
        aHi = mid;
      if (aHi - aLo < 1e-12 * aHi)
        break;
    }
    double alpha = 0.5 * (aLo + aHi);
    if (alpha <= 0.0 || alpha > 1e6)
      continue;  // numerical failure -- give up on this resonance

    // Apply the rescaling: for each daughter, scale the 3-momentum by alpha
    // (in the resonance rest frame), recompute the energy with the (new) mass,
    // and boost back to the lab frame.
    Pythia8::RotBstMatrix fromRes;
    fromRes.bst(pRes);

    for (const auto& dk : dks) {
      int iD = dk.idx;

      // Save the daughter's original lab-frame 4-momentum before rescaling.
      Pythia8::Vec4 pDOrig = event[iD].p();

      // Apply the rescaling in the resonance rest frame.
      Pythia8::Vec4 pD = pDOrig;
      pD.rotbst(toRes);
      double oldMag = pD.pAbs();
      double newMass = dk.isB ? 0.0 : event[iD].m();
      if (oldMag > 0.0) {
        pD.rescale3(alpha);
        double newPMag = alpha * oldMag;
        pD.e(std::sqrt(newPMag * newPMag + newMass * newMass));
      } else {
        pD.e(newMass);
      }
      pD.rotbst(fromRes);

      // Install the daughter's new 4-momentum.
      event[iD].p(pD);
      if (dk.isB) {
        event[iD].m(0.0);
      }

      // Propagate the daughter's 4-momentum change to ALL of its
      // descendants (recursively). This is essential for kinematic
      // consistency: when a daughter's 4-momentum changes (e.g. the W from
      // t -> b W), the daughter's descendants' 4-momenta must also change
      // so that the daughter still equals the sum of its descendants.
      // Otherwise, the broken momentum conservation in the W decay produces
      // NaN values in later PDF / pT evaluations in the Powheg-matched
      // shower, manifesting as the "Unphysical x given: -nan" error.
      //
      // We use the Lorentz boost that maps the daughter's original
      // 4-momentum to its new 4-momentum: T = bst(pDOrig, pDNew). Because T
      // is a Lorentz transformation, applying T to every descendant (at
      // every level of the decay tree, recursively) preserves every
      // "resonance = sum of its descendants" relation.
      //
      // The boost formula requires pOrig^2 == pNew^2, i.e. the daughter's
      // invariant mass must be preserved. This is true for non-b daughters
      // (whose mass is preserved by the rescaling). For b daughters, the mass
      // is changed to 0 -- but the b is a leaf of the decay tree at this
      // point (it hadronizes later), so it has no descendants and the
      // propagation is not needed for it.
      if (alpha != 1.0) {
        Pythia8::Vec4 pDNew = event[iD].p();
        double m2Orig = pDOrig.m2Calc();
        double m2New  = pDNew.m2Calc();
        if (std::abs(m2Orig - m2New) <=
            1e-9 * (m2Orig + m2New + 1.0)) {
          Pythia8::RotBstMatrix T;
          T.bst(pDOrig, pDNew);
          transformDescendants(event, iD, T);
        }
      }
    }
  }

  // Do not veto; let the shower proceed with the modified b quarks.
  return false;
}
