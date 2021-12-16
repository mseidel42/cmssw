/**
  @file PatJetHitFitTranslator.cc

  @brief Specialization of template class JetTranslatorBase in the
  package HitFit

  @author Haryo Sumowidagdo <Suharyo.Sumowidagdo@cern.ch>

  @par Created
  Sat Jun 27 17:49:21 2009 UTC

  @version $Id: PatJetHitFitTranslator.cc,v 1.2 2011/05/26 12:57:18 mseidel Exp $
 */


#include <TopQuarkAnalysis/TopHitFit/interface/JetTranslatorBase.h>
#include <DataFormats/PatCandidates/interface/Jet.h>

using std::string;
using std::cout;
using std::endl;

namespace hitfit {

template<>
JetTranslatorBase<pat::Jet>::JetTranslatorBase()
{
  string CMSSW_BASE(getenv("CMSSW_BASE"));
  string resolution_filename = CMSSW_BASE + string("/src/TopQuarkAnalysis/PatHitFit/data/exampleJetResolution.txt");
  udscResolution_ = EtaDepResolution(resolution_filename);
  bResolution_    = EtaDepResolution(resolution_filename);
} // JetTranslatorBase<pat::Jet>::JetTranslatorBase()


template<>
JetTranslatorBase<pat::Jet>::JetTranslatorBase(const string& udscFile, const string& bFile)
{
  string CMSSW_BASE(getenv("CMSSW_BASE"));

  string udscResolution_filename = udscFile;
  if (udscFile.empty()) udscResolution_filename = CMSSW_BASE + string("/src/TopQuarkAnalysis/TopHitFit/data/resolution/tqafUdscJetResolution.txt");

  string bResolution_filename = bFile;
  if (bFile.empty()) bResolution_filename = CMSSW_BASE + string("/src/TopQuarkAnalysis/TopHitFit/data/resolution/tqafBJetResolution.txt");

  udscResolution_ = EtaDepResolution(udscResolution_filename);
  bResolution_    = EtaDepResolution(bResolution_filename);
} // JetTranslatorBase<pat::Jet>::JetTranslatorBase(const string& ifile)


template<>
JetTranslatorBase<pat::Jet>::~JetTranslatorBase()
{
} // JetTranslatorBase<pat::Jet>::~JetTranslatorBase()


template<>
Lepjets_Event_Jet
JetTranslatorBase<pat::Jet>::operator()(const pat::Jet& jet, int type /*= hitfit::unknown_label */)
{
  const bool bCase = type == hitfit::hadb_label || type == hitfit::lepb_label || type == hitfit::higgs_label;
  // If Calo jet usage is resurrected: jet->isCaloJet() ? ((reco::CaloJet*) jet->originalObject())->detectorP4().eta()
  const double jet_eta = jet.eta();

  Vector_Resolution jet_resolution = bCase ? bResolution_.GetResolution(jet_eta) : udscResolution_.GetResolution(jet_eta);
  
  // In the past, flavor-tags have been used in the higher-level JECs.
  // This is not done anymore, and the feature has been disabled, but one could resurrect it e.g. by:
  // pat::Jet partonCorrJet(jet.correctedJet(jetCorrectionLevel_, bCase ? "BOTTOM" : "UDS"));
  const Fourvec p(jet.px(),jet.py(),jet.pz(),jet.energy());

  return Lepjets_Event_Jet(p, type, jet_resolution);
} // Lepjets_Event_Jet JetTranslatorBase<pat::Jet>::operator()(const pat::Jet& j,int type)


template<>
const EtaDepResolution&
JetTranslatorBase<pat::Jet>::udscResolution() const
{
  return udscResolution_;
}


template<>
const EtaDepResolution&
JetTranslatorBase<pat::Jet>::bResolution() const
{
  return bResolution_;
}


template<>
bool
JetTranslatorBase<pat::Jet>::CheckEta(const pat::Jet& jet) const
{
  double jet_eta = jet.eta();

  //if (jet.isPFJet()) do nothing at the moment!
  if (jet.isCaloJet()) jet_eta = static_cast<const reco::CaloJet*>(jet.originalObject())->detectorP4().eta();
  return bResolution_.CheckEta(jet_eta) and udscResolution_.CheckEta(jet_eta);
}

} // namespace hitfit
