#ifndef FCCANALYSES_PPP_H
#define FCCANALYSES_PPP_H

#include "FCCAnalyses/JetConstituentsUtils.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/ReconstructedParticle2Track.h"
#include "FCCAnalyses/ReconstructedParticle2MC.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterData.h"
#include "edm4hep/CalorimeterHitData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/EDM4hepVersion.h"
#include "FCCAnalyses/JetClusteringUtils.h"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

#include "ROOT/RVec.hxx"

namespace FCCAnalyses {
namespace JetConstituentsUtils {

class PPP {
public:

    // 1️⃣ cast_constituent template must be defined before get_px/get_py/get_pz
    template <typename F>
    static rv::RVec<FCCAnalysesJetConstituentsData> cast_constituent(
        const rv::RVec<FCCAnalysesJetConstituents>& jcs, F&& meth)
    {
        rv::RVec<FCCAnalysesJetConstituentsData> out;
        for (const auto &jc : jcs)
            out.push_back(meth(jc)); // use push_back, simpler
        return out;
    }

    // 2️⃣ Cartesian momentum components as static functions
    static rv::RVec<FCCAnalysesJetConstituentsData> get_px(
        const rv::RVec<FCCAnalysesJetConstituents>& jcs)
    {
        return cast_constituent(jcs, ReconstructedParticle::get_px);
    }

    static rv::RVec<FCCAnalysesJetConstituentsData> get_py(
        const rv::RVec<FCCAnalysesJetConstituents>& jcs)
    {
        return cast_constituent(jcs, ReconstructedParticle::get_py);
    }

    static rv::RVec<FCCAnalysesJetConstituentsData> get_pz(
        const rv::RVec<FCCAnalysesJetConstituents>& jcs)
    {
        return cast_constituent(jcs, ReconstructedParticle::get_pz);
    }

    static rv::RVec<int> mask_pt_positive(
    const rv::RVec<FCCAnalysesJetConstituents>& jcs)
    {
    rv::RVec<int> mask;
    mask.reserve(jcs.size());
    for (const auto &jc : jcs) {
        auto px = ReconstructedParticle::get_px(jc)[0];
        auto py = ReconstructedParticle::get_py(jc)[0];
        float pt2 = px*px + py*py;
        mask.push_back(pt2 > 0.0 ? 1 : 0);

    }
    return mask;
    }

    
    

};

} // namespace JetConstituentsUtils
} // namespace FCCAnalyses

#endif // FCCANALYSES_PPP_H
