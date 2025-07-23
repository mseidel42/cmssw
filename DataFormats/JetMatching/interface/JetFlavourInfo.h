#ifndef DataFormats_JetMatching_JetFlavourInfo_H
#define DataFormats_JetMatching_JetFlavourInfo_H

#include <vector>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

/// For fastjet::contrib::FlavInfo handling
#include "fastjet/contrib/FlavInfo.hh"

namespace reco {
  /**\class JetFlavourInfo JetFlavourInfo.h DataFormats/JetMatching/interface/JetFlavourInfo.h
 * \brief Class storing the jet flavour information
 *
 * JetFlavourInfo class stores the jet flavour information based on hadrons
 * and partons clustered inside the jet. It also provides vectors of
 * EDM references to clustered hadrons and partons. The hadron- and parton-based
 * flavours are defined in the JetFlavourClustering producer.
 *
 * The flavour definition algorithms introduced in fastjet::constrib outputs
 * the flavours in the type of fastjet::contrib::FlavInfo, which is fully
 * characterized by its member _flav_content, an int[7] array. To preserve the
 * info inside a FlavInfo object, we define a type reco::FJContribFlavArray
 * which is a typedef of std::vector<int>.
 *  * (For full info see https://github.com/cms-externals/fastjet-contrib/blob/cms/v1.101/IFNPlugin/include/fastjet/contrib/IFNPlugin.hh)
 * 
 * For the result of each flavour definition algorithm, an std::vector is
 * always assigned and will remain empty if no flavour information given.
 *
 * 
 */
  /// typedef for 
  /// enumeration class for the fastjet::contrib flavour definition algorithms
  enum class FJContribFlavDef {
    kCMP,
    kGHS,
    kIFN,
    kSDF
  };
  // Number of fastjet::contrib flavour definition algorithms
  static constexpr kFJContribFlavAlgoCount = 4;

  class JetFlavourInfo {
  public:
    JetFlavourInfo(void) : m_hadronFlavour(0), m_partonFlavour(0) {}
    JetFlavourInfo(const int hadronFlavour, const int partonFlavour)
        : m_hadronFlavour(hadronFlavour), m_partonFlavour(partonFlavour) {}
    JetFlavourInfo(const GenParticleRefVector& bHadrons,
                   const GenParticleRefVector& cHadrons,
                   const GenParticleRefVector& partons,
                   const GenParticleRefVector& leptons,
                   const int hadronFlavour,
                   const int partonFlavour)
        : m_bHadrons(bHadrons),
          m_cHadrons(cHadrons),
          m_partons(partons),
          m_leptons(leptons),
          m_hadronFlavour(hadronFlavour),
          m_partonFlavour(partonFlavour) {}

    /// Return a vector of GenParticleRef's to b hadrons clustered inside the jet
    const GenParticleRefVector& getbHadrons() const { return m_bHadrons; }
    /// Return a vector of GenParticleRef's to c hadrons clustered inside the jet
    const GenParticleRefVector& getcHadrons() const { return m_cHadrons; }
    /// Return a vector of GenParticleRef's to partons clustered inside the jet
    const GenParticleRefVector& getPartons() const { return m_partons; }
    /// Return a vector of GenParticleRef's to leptons clustered inside the jet
    const GenParticleRefVector& getLeptons() const { return m_leptons; }
    /// Return the hadron-based flavour
    const int getHadronFlavour() const { return m_hadronFlavour; }
    /// Return the parton-based flavour
    const int getPartonFlavour() const { return m_partonFlavour; }

    /// Set the hadron-based flavour
    void setHadronFlavour(const int hadronFlavour) { m_hadronFlavour = hadronFlavour; }
    /// Set the parton-based flavour
    void setPartonFlavour(const int partonFlavour) { m_partonFlavour = partonFlavour; }

    /// view if flavour by some fastjet::contrib algorithm is defined
    bool haveFJContribFlavAlgo(const FJContribFlavDef& algo) const { return !m_fjContribFlav[static_cast<int>(algo)].empty(); }
    /// Set the flavour defined by some fastjet::contrib algorithm, on one digit, only if an array already exists.
    void setFJContribFlavAlgo(const FJContribFlavDef& algo, const fastjet::contrib::FlavInfo& fjFlavInfo) { m_fjContribFlav[static_cast<int>(algo)].assign(fjFlavInfo._flav_content, fjFlavInfo._flav_content + 7); }
    /// Obtain the flavour defined by some fastjet::contrib algorithm, either by digit or by array. 
    int getFJContribFlavAlgo(const FJContribFlavDef& algo, const unsigned int& iflav) const { return m_fjContribFlav[static_cast<int>(algo)][iflav]; }
    std::vector<int>& getFJContribFlavAlgo(const FJContribFlavDef& algo) { return m_fjContribFlav[static_cast<int>(algo)]; }
    /// When needed, clear the flavour defined by some fastjet::contrib algorithm.
    void clearFJContribFlavAlgo(const FJContribFlavDef& algo) { m_fjContribFlav[static_cast<int>(algo)].clear(); }

  private:
    GenParticleRefVector m_bHadrons;
    GenParticleRefVector m_cHadrons;
    GenParticleRefVector m_partons;
    GenParticleRefVector m_leptons;
    /// Original hadronFlavour and partonFlavour definition.
    int m_hadronFlavour;
    int m_partonFlavour;
    /// fastjet::contrib algorithm flavour definition (arXiv:2205.01109)
    std::vector<int> m_fjContribFlav[kFJContribFlavAlgoCount];
  };

}  // namespace reco
#endif
