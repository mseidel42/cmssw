// MasslessLHEInput.cc
//
// Implementation of the MasslessLHEInputHook that makes the final-state
// b quarks coming from the LHE massless before the shower runs, while
// preserving the mass of any resonance (top, Z, W, Higgs, ...) that the b
// came from.
//
// Two cases are handled:
//  (a) b from a resonance decay (e.g. t -> b W): The b's momentum is rescaled
//      in the resonance rest frame, together with the other resonance
//      daughters, by a common factor alpha chosen so that (i) the b is
//      massless and (ii) the total energy equals the resonance mass. This
//      preserves the resonance mass exactly. For a two-body t -> b W decay
//      the formula reduces to the standard two-body kinematics with m_b = 0.
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
  // Guard against cycles / running past the start of the event record.
  for (int step = 0; step < event.size() && i > 0; ++step) {
    int m1 = event[i].mother1();
    int m2 = event[i].mother2();
    // Pick the first valid mother.
    int iMother = (m1 > 0) ? m1 : m2;
    if (iMother <= 0 || iMother >= event.size())
      break;
    if (isResonanceId(std::abs(event[iMother].id()))) {
      // The b came from this resonance's decay.
      return iMother;
    }
    i = iMother;
  }
  return -1;
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

    // Skip incoming b quarks (e.g. 5FS b-initiated hard processes).
    // Outgoing final-state b's have status > 0 (isFinal) and the typical
    // LHE status codes for outgoing ME particles; incoming particles have
    // negative status and are NOT isFinal. So the isFinal() check above
    // already excludes incoming b's.

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

    // Collect the (direct) daughters of the resonance.
    std::vector<int> daughters = event[iRes].daughterList();
    if (daughters.empty())
      continue;

    // For each daughter, get the 4-momentum in the resonance rest frame.
    // Identify which daughter is the b we are processing (or a bbar).
    // We rescale ALL daughters by alpha, making any final-state b massless.
    Pythia8::RotBstMatrix toRes;
    toRes.bstback(pRes);

    // Compute the uniform scaling factor alpha by solving
    //   f(alpha) = sum_over_daughters E_i'(alpha) - mRes = 0,
    // where for b daughters E_i' = alpha * |p_i| (massless) and for other
    // daughters E_i' = sqrt(alpha^2 * p_i^2 + m_i^2) (on-shell with original
    // mass). f is monotonically increasing in alpha, so a unique root exists
    // and can be found by bisection.
    //
    // Pre-compute, per daughter, the 3-momentum magnitude and the mass, in
    // the resonance rest frame.
    struct DaughterKin {
      int idx;
      double pMag;  // 3-momentum magnitude in the resonance rest frame
      double mass;  // original mass
      bool isB;      // is this a b/bbar we want to make massless?
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
      dk.mass = event[iD].m();
      dk.isB = (std::abs(event[iD].id()) == 5 && event[iD].isFinal());
      if (dk.isB) {
        // Will be made massless: energy becomes alpha * pMag.
        // Note: mass is what we're going to set to 0; we use mass = 0 in f.
        dk.mass = 0.0;
      }
      dks.push_back(dk);
    }
    if (dks.empty())
      continue;

    // f(alpha) = sum_d E_d'(alpha) - mRes
    auto totalEnergy = [&](double alpha) -> double {
      double E = 0.0;
      for (const auto& dk : dks) {
        double pAlpha = alpha * dk.pMag;
        E += std::sqrt(pAlpha * pAlpha + dk.mass * dk.mass);
      }
      return E;
    };
    auto f = [&](double alpha) { return totalEnergy(alpha) - mRes; };

    // Bisection: find alpha > 0 such that f(alpha) = 0.
    // f is monotonically increasing. At alpha = 0, f = sum_d m_d - mRes (a
    // large negative number). As alpha -> infinity, f -> +infinity. So a
    // unique root exists. Start with a bracket and tighten.
    double aLo = 0.0, aHi = 1.0;
    // Expand aHi until f(aHi) > 0.
    while (f(aHi) < 0.0 && aHi < 1e6) {
      aHi *= 2.0;
    }
    // Bisect.
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
      Pythia8::Vec4 pD = event[iD].p();
      pD.rotbst(toRes);
      double oldMag = pD.pAbs();
      double newMass = dk.isB ? 0.0 : event[iD].m();
      if (oldMag > 0.0) {
        // Scale the 3-momentum by alpha, keep the direction, and set the
        // energy from the (new) on-shell mass.
        double scale = alpha;  // 3-momentum scaling factor
        pD.rescale3(scale);
        double newPMag = alpha * oldMag;
        pD.e(std::sqrt(newPMag * newPMag + newMass * newMass));
      } else {
        // Daughter at rest in the resonance rest frame: only energy changes.
        pD.e(newMass);
      }
      pD.rotbst(fromRes);
      event[iD].p(pD);
      if (dk.isB) {
        event[iD].m(0.0);
      }
    }
  }

  // Do not veto; let the shower proceed with the modified b quarks.
  return false;
}
