// test273.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Includes a user hook that corrects emission in top decay for dipole
// from gluon to W, to instead be from gluon to top.

// Important: the top mass shift analysis encoded here is very primitive,
// does not perform well at all, and should not be taken seriously.
// The important part is that you see how the different scenarios
// should be set up to operate as intended.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// Write own derived UserHooks class for modified emission in top decay.

class TopRecoilHook : public UserHooks {

 public:
  
  // Constructor.
  //  doTopRecoil : eikonal correction in GW dipole on/off when no MEC applied.
  //  useOldDipole  : in GW dipole, use partons before or after branching.
  //  doList        : diagnostic output; set false for production runs.
  TopRecoilHook(bool doTopRecoilIn=true, bool useOldDipoleIn=false,
    bool doListIn = false) {
    doTopRecoil = doTopRecoilIn;
    useOldDipole = useOldDipoleIn;
    doList = doListIn; 
    // Constructor also creates some histograms for analysis inside User Hook.
    wtCorr = new Hist("corrective weight", 100, 0., 2.);
  }

  // Destructor prints histogram.
  ~TopRecoilHook() override {
    if (doTopRecoil) cout << *wtCorr;
    delete wtCorr;
  }

  // Initialise. Only use hook for simple showers with recoilToColoured = off.
  virtual bool initAfterBeams() override {
    int showerModel  = settingsPtr->mode("PartonShowers:Model");
    // Switch off if not using simple showers or if recoilToColoured = on.
    bool recoilToColoured = settingsPtr->flag("TimeShower:recoilToColoured");
    if (showerModel != 1 || recoilToColoured) doTopRecoil=false;
    // Flag if W mass term is already accounted for (true) or not (false).
    recoilDeadCone        = settingsPtr->flag("TimeShower:recoilDeadCone");
    // All ok.
    return true;
  }

  // Allow a veto after an FSR emission
  virtual bool canVetoFSREmission() override {return doTopRecoil;}
  
  // Access the event after an FSR emission, specifically inside top decay.
  virtual bool doVetoFSREmission( int sizeOld, const Event& event, int iSys,
    bool inResonance) override {
    
    // Check that we are inside a resonance decay.
    if (!inResonance) return false;

    // Check that it is a top decay. 
    int iTop = partonSystemsPtr->getInRes(iSys);
    if (iTop == 0 || event[iTop].idAbs() != 6) return false;
    
    // Skip first emission, where ME corrections are already made.
    int sizeOut = partonSystemsPtr->sizeOut(iSys);
    if (sizeOut == 2) return false;
    
    // Location of trial new particles: radiator, emitted, recoiler.
    int iRad = sizeOld;
    int iEmt = sizeOld + 1;
    int iRec = sizeOld + 2;
    
    // The above partons are after emission;
    // alternatively use the ones before.
    if (useOldDipole) {
      iRad = event[iRad].mother1();
      iRec = event[iRec].mother1();
    }
    
    // Check if newly emitted gluon matches (anti)top colour line.
    if (event[iEmt].id() != 21) return false;
    if (event[iTop].id() == 6) {
      if (event[iEmt].col() != event[iTop].col()) return false;
    } else {
      if (event[iEmt].acol() != event[iTop].acol()) return false;
    }
    
    // Recoiler should now be a W, else something is wrong.
    if (event[iRec].idAbs() != 24) {
      cout << " ERROR: recoiler is " << event[iRec].id() << endl;
      return false;
    }
    
    // Denominator: eikonal weight with W as recoiler.
    double pRadRec = event[iRad].p() * event[iRec].p();
    double pRadEmt = event[iRad].p() * event[iEmt].p();
    double pRecEmt = event[iRec].p() * event[iEmt].p();
    double wtW = 2. * pRadRec / (pRadEmt * pRecEmt)
      - pow2(event[iRad].m() / pRadEmt);
    // If recoilDeadCone = on, include W mass term in denominator.
    if (recoilDeadCone) wtW -= pow2(event[iRec].m() / pRecEmt);

    // Numerator: eikonal weight with top as recoiler.
    double pRadTop = event[iRad].p() * event[iTop].p();
    double pTopEmt = event[iTop].p() * event[iEmt].p();
    double wtT = 2. * pRadTop / (pRadEmt * pTopEmt)
      - pow2(event[iRad].m() / pRadEmt) - pow2(event[iTop].m() / pTopEmt);
    
    // Histogram weight ratio.
    wtCorr->fill( wtT / wtW );
    
    // List relevant properties.
    if (doList) {
      cout << "\n now event with sizeOld = " << sizeOld << ", iSys = "
           << iSys << ", sizeOut = " << sizeOut << scientific
           << setprecision(3)
           << ", weight with W = " << wtW << " and with t = " << wtT << endl;
      partonSystemsPtr->list();
      event.list();
    }
    
    // Accept/reject emission. Smooth suppression or step function.
    return (wtT < wtW * rndmPtr->flat());
  }

 private:

  // Options and Histograms.
  bool  doTopRecoil, useOldDipole, doList, recoilDeadCone;
  Hist *wtCorr;
  
};

//==========================================================================

int main() {

  // Number of events to generate.
  // Warning: much statistics is needed for significant results,
  // so this is just an appetizer. Anyway, the reconstruction is
  // pretty lousy, so not useful for real studies.
  int nEvent = 1000000;
  map<int,UserHooksPtr> topRecoilHooks;
  
  // Loop over different scenarios. Use OpenMP parallelisation if enabled.
#ifdef OPENMP
  #pragma omp parallel for
#endif
  for (int mLoop = 0; mLoop <= 3; ++mLoop) {

    // Set up anti-kT jet finder.
    double Rjet = 0.5;
    double pTjetMin = 20.;
    SlowJet sJet( -1, Rjet, pTjetMin);
    
    // Generator at LHC.
    Pythia pythia;
    Event& event = pythia.event;
    pythia.readString("Beams:eCM = 13000.");
    pythia.readString("6:m0 = 173.3");
    pythia.readString("24:m0 = 80.385");
    
    // Reduce printout.
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:numberCount = 100000");
    
    // Switch off MPI, hadronisation, and W decay (for super-simplicity!)
    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("HadronLevel:all = off");
    pythia.readString("24:mayDecay = off");
    
    // q qbar, g g -> t tbar.
    pythia.readString("Top:qqbar2ttbar = on");
    pythia.readString("Top:gg2ttbar = on");
    
    // mLoop  = 0 : baseline with recoilToColoured = on.
    // mLoop >= 1 : baseline with recoilToColoured = off.
    // mLoop  = 2 : GW->GT eikonal correction (with post-branching dipole).
    // mLoop  = 3 : same as =3 but recoilDeadCone = off as starting point.
    bool doTopRecoil  = false;
    bool useOldDipole = false;
    bool doList       = false;
    if (mLoop == 0) pythia.readString("TimeShower:recoilToColoured = on");
    else pythia.readString("TimeShower:recoilToColoured = off");
    if (mLoop >= 2) doTopRecoil = true;
    if (mLoop == 3) pythia.readString("TimeShower:recoilDeadCone = off");
    
    // Set up to do eikonal 
    topRecoilHooks[mLoop] = make_shared<TopRecoilHook>(doTopRecoil, 
      useOldDipole, doList);
    if (mLoop >= 2) pythia.setUserHooksPtr( topRecoilHooks[mLoop] );
    
    // Do the initialisation sequentially (for nicer output).
#ifdef OPENMP
  #pragma omp critical
#endif
    {
      cout<<"\nINITIALISING MODE: "<<mLoop<<"\n";
      pythia.init();
    } // End omp critical section.
    
    // Target top mass.
    double mT = pythia.particleData.m0(6);

    // Histograms for current scenario.
    Hist nH(     "multiplicity",            100,  -0.5, 199.5);
    Hist nJetH(  "jet multiplicity",              20, -0.5, 19.5);
    Hist mTH(    "reconstructed t mass",         100, 130., 210.);
    Hist mTpTH(  "reconstructed delta-m_t(pT_t)", 11,   0., 275.);
    Hist mTerrH( "reconstructed t mass error",   100, -5,  5.0);
    Hist pTTH(   "reconstructed pT_t",            11,   0., 275.);
    Hist ebH(    "energy of b quark",            100,   0., 500.);

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Multiplicity. Also find W bosons.
      int nF = 0;
      vector<int> iW;
      for (int i = 0; i < event.size(); ++i) {
        if (!event[i].isFinal()) continue;
        ++nF;
        if (event[i].idAbs() == 24) {
          // Set status < 0 so not included in jet clustering below.
          event[i].statusNeg();
          iW.push_back(i);
        }
        if (event[i].idAbs() == 5) ebH.fill(event[i].e());
      }
      nH.fill(nF);

      // Find number of jets. At least two to keep going.
      sJet.analyze(event);
      int nJet = sJet.sizeJet();
      nJetH.fill( nJet);
      if (nJet < 2) continue;

      // Find best jet-W mass pair, closest to mT.
      double diff   = 1e10;      
      double mWJrec = 0.;
      Vec4 pRec; 
      for (int i = 0; i <= 1; ++i) {        
        for (int j=0; j<nJet; ++j) {
          double mWJ = (event[iW[i]].p() + sJet.p(j)).mCalc();
          if (abs(mWJ - mT) < diff) {
            mWJrec  = mWJ;
            pRec    = event[iW[i]].p() + sJet.p(j);
            diff    = abs(mWJ - mT);
          }
        }
      }
      mTH.fill( mWJrec);
      mTerrH.fill( mWJrec - mT );

      // Only keep going if within +-20 GeV.
      if (abs(mWJrec - mT) > 20.) continue;

      // Study top pT and dependence of top mass error.
      double pTT = pRec.pT();
      if (pTT > 250.) pTT = 260.;
      pTTH.fill( pTT);
      mTpTH.fill( pTT, mWJrec - mT);

    // End of event loop. Statistics. Histograms.
    }

#ifdef OPENMP
  #pragma omp critical
#endif
    {
      cout<<"\nFINISHED MODE: "<<mLoop<<"\n";
      pythia.stat();
      mTpTH /= pTTH;
      cout <<  nH << nJetH << mTH 
           << mTerrH << pTTH << mTpTH << ebH;
      stringstream ss;
      ss<<"mTerr."<<setw(1)<<mLoop;
      mTerrH.table(ss.str());
      ss.str("");
      ss<<"eb."<<setw(1)<<mLoop;
      ebH.table(ss.str());
    }
    // End loop over models.
  }

  // Done.
  return 0;
}
