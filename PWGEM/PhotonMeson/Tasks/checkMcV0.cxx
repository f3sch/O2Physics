// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file checkMcV0.cxx
/// \brief Study the reconstruction and reconstructability of photon conversions in MC
/// \author daiki.sekihata@cern.ch felix.schlepper@cern.ch

#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"
#include "PWGMM/Mult/DataModel/Index.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TMCProcess.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::constants::physics;
using namespace o2::pwgem::photonmeson;

struct CheckMcV0 {
  using TracksMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels>;
  using FullMcParticles = soa::Join<aod::McParticles, aod::ParticlesToTracks>;
  using CollisionsMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;

  // Shared truth-photon fiducial selection
  Configurable<bool> requirePrimaryOrGeneratorGamma{"requirePrimaryOrGeneratorGamma", true, "require the photon to be a physical primary or produced by the generator"};
  Configurable<int> maxTPCDriftTime{"maxTPCDriftTime", 4000, "maximum TPC drift time in BC"};
  Configurable<int> minNCont{"minNCont", 10, "Minimum number of contributor to account collision"}; // surpress UPC cont.
  Configurable<float> minPtGamma{"minPtGamma", 0.01f, "minimum truth photon pT in GeV/c"};
  Configurable<float> maxPtGamma{"maxPtGamma", 9999999.f, "maximum truth photon pT in GeV/c"};
  Configurable<float> maxEtaGamma{"maxEtaGamma", 0.9f, "maximum absolute truth photon eta"};
  Configurable<float> minConversionRadius{"minConversionRadius", 2.f, "minimum truth conversion radius in cm"};
  Configurable<float> maxConversionRadius{"maxConversionRadius", 90.f, "maximum truth conversion radius in cm"};
  Configurable<float> maxAbsConversionZ{"maxAbsConversionZ", -1.f, "maximum absolute truth conversion z in cm; negative disables the cut"};
  Configurable<bool> applyRZLineCut{"applyRZLineCut", true, "apply the photon conversion R-Z fiducial line cut"};
  Configurable<float> marginZMC{"marginZMC", 7.f, "margin for the photon conversion R-Z fiducial line cut in cm"};
  // Shared reconstructed-leg selection
  Configurable<bool> requireTPCForLeg{"requireTPCForLeg", true, "require a TPC-containing reconstructed track for each conversion leg"};
  Configurable<bool> requireCorrectMcTrackLabel{"requireCorrectMcTrackLabel", true, "require conversion-leg tracks without MC label mismatch or fake bits"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  // Histograms
  HistogramRegistry registry{"output", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis cfgAxisGammaPt{"cfgAxisGammaPt", {VARIABLE_WIDTH, 0.01f, 0.02f, 0.03f, 0.04f, 0.05f, 0.06f, 0.07f, 0.08f, 0.09f, 0.1f, 0.15f, 0.2f, 0.25f, 0.3f, 0.35f, 0.4f, 0.45f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.f, 2.2f, 2.5f, 3.f, 4.f, 5.f, 6.f, 8.f, 10.f, 15.f, 20.f}, "truth photon pT axis"};
  AxisSpec axisGammaPt{cfgAxisGammaPt, "#it{p}_{T}^{#gamma} (GeV/#it{c})"};
  AxisSpec axisLegPt{cfgAxisGammaPt, "#it{p}_{T}^{leg,true} (GeV/#it{c})"};
  AxisSpec axisElectronPt{cfgAxisGammaPt, "#it{p}_{T}^{e^{-},true} (GeV/#it{c})"};
  AxisSpec axisPositronPt{cfgAxisGammaPt, "#it{p}_{T}^{e^{+},true} (GeV/#it{c})"};
  AxisSpec axisPairPt{cfgAxisGammaPt, "|#vec{p}_{T}^{e^{-}} + #vec{p}_{T}^{e^{+}}| (GeV/#it{c})"};
  AxisSpec axisLegPtSum{cfgAxisGammaPt, "#it{p}_{T}^{e^{-}} + #it{p}_{T}^{e^{+}} (GeV/#it{c})"};
  ConfigurableAxis cfgAxisMinLegPt{"cfgAxisMinLegPt", {VARIABLE_WIDTH, 0.001f, 0.005f, 0.01f, 0.02f, 0.03f, 0.04f, 0.05f, 0.06f, 0.07f, 0.08f, 0.09f, 0.1f, 0.15f, 0.2f, 0.25f, 0.3f, 0.35f, 0.4f, 0.45f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.f, 2.2f, 2.5f, 3.f, 4.f, 5.f, 6.f, 8.f, 10.f, 15.f, 20.f}, "minimum truth daughter pT axis"};
  AxisSpec axisMinLegPt{cfgAxisMinLegPt, "min(#it{p}_{T}^{e^{-}}, #it{p}_{T}^{e^{+}}) (GeV/#it{c})"};
  ConfigurableAxis cfgAxisMomentumFraction{"cfgAxisMomentumFraction", {100, 0.f, 1.f}, "positron momentum fraction axis"};
  AxisSpec axisMomentumFraction{cfgAxisMomentumFraction, "#it{p}_{T}^{e^{+}}/(#it{p}_{T}^{e^{-}} + #it{p}_{T}^{e^{+}})"};
  ConfigurableAxis cfgAxisConversionRadius{"cfgAxisConversionRadius", {400, 0.f, 100.f}, "truth conversion radius axis"};
  AxisSpec axisConversionRadius{cfgAxisConversionRadius, "#it{R}_{conv}^{true} (cm)"};
  ConfigurableAxis cfgAxisConversionPhi{"cfgAxisConversionPhi", {400, 0.f, o2::constants::math::TwoPI}, "truth conversion phi axis"};
  AxisSpec axisConversionPhi{cfgAxisConversionPhi, "#varphi_{conv}^{true} (rad.)"};
  ConfigurableAxis cfgAxisConversionZ{"cfgAxisConversionZ", {400, -100.f, 100.f}, "truth conversion z axis"};
  AxisSpec axisConversionZ{cfgAxisConversionZ, "#it{z}_{conv}^{true} (cm)"};
  ConfigurableAxis cfgAxisOccupancyTracks{"cfgAxisOccupancyTracks", {VARIABLE_WIDTH, -1.f, 0.f, 100.f, 500.f, 1500.f, 3700.f, 8000.f, 15000.f, 50000.f}, "track occupancy axis"};
  AxisSpec axisOccupancyTracks{cfgAxisOccupancyTracks, "track occupancy in time range"};
  ConfigurableAxis cfgAxisCentrality{"cfgAxisCentrality", {110, 0.f, 110.f}, "centrality axis"};
  AxisSpec axisCentrality{cfgAxisCentrality, "centrality (%)"};
  ConfigurableAxis cfgAxisNCol{"cfgAxisNCol", {1001, -1, 1000}, "number of other collisions in drift time"};
  AxisSpec axisNCol{cfgAxisNCol, "number of other collisions in drift time"};
  ConfigurableAxis cfgAxisMultiplicity{"cfgAxisMultiplicity", {101, -0.5f, 100.5f}, "candidate or track multiplicity axis"};
  AxisSpec axisMultiplicity{cfgAxisMultiplicity, "multiplicity"};
  AxisSpec axisV0Type{16, -0.5f, 15.5f, "raw V0 type"};
  AxisSpec axisTrackType{5, -0.5f, 4.5f, "raw V0 daughter track type"};
  AxisSpec axisCollisionAssociation{4, -0.5f, 3.5f, "raw V0 collision association"};

  enum RawV0TrackType : int {
    kITSTPCITSTPC = 0,
    kITSTPCTPCOnly,
    kTPCOnlyTPCOnly,
    kContainsITSOnly,
    kOther
  };

  enum RawV0CollisionAssociation : int {
    kExactRecoCollision = 0,
    kOtherRecoSameMcCollision,
    kOtherMcCollision,
    kNoMcCollisionLabel
  };

  struct RawV0Match {
    int64_t collisionId{-1};
    int64_t mcCollisionId{-1};
    int64_t posTrackId{-1};
    int64_t negTrackId{-1};
    uint8_t v0Type{0};
    int trackType{kOther};
    bool hasCleanLabels{false};
  };

  struct ReconstructedTrackSummary {
    int nSameCollision{0};
    std::vector<int64_t> eligibleTrackIds;
  };

  void init(InitContext const& /*unused*/)
  {
    if (doprocessTrueConversion) {
      registry.add("ConversionProbability/hGeneratedGammaCount", "selected generated photons", HistType::kTH1D, {{1, 0.5f, 1.5f}});
      registry.add("ConversionProbability/hGeneratedGammaPt", "selected generated photons", HistType::kTH1D, {axisGammaPt});
      registry.add("ConversionProbability/hConvertedGammaPt", "selected converted photons in the fiducial volume", HistType::kTH1D, {axisGammaPt});
      registry.add("ConversionProbability/hConversionRadius", "conversion-radius numerator for the cumulative conversion probability", HistType::kTH1D, {axisConversionRadius});
      registry.add("ConversionProbability/hRPhiTrueConversion", "true photon conversion points", HistType::kTH2D, {axisConversionRadius, axisConversionPhi});
      registry.add("ConversionProbability/hRZTrueConversion", "true photon conversion points", HistType::kTH2D, {axisConversionRadius, axisConversionZ});
    }
    if (doprocessMCV0) {
      const std::vector<AxisSpec> legEfficiencyAxes{axisConversionRadius, axisLegPt, axisOccupancyTracks, axisCentrality, axisNCol};
      const std::vector<AxisSpec> partnerEfficiencyAxes{axisConversionRadius, axisElectronPt, axisPositronPt, axisOccupancyTracks, axisCentrality, axisNCol};
      const std::vector<AxisSpec> conversionEfficiencyAxes{axisConversionRadius, axisGammaPt, axisOccupancyTracks, axisCentrality, axisNCol};
      const std::vector<AxisSpec> conversionDaughterKinematicsAxes{axisConversionRadius, axisGammaPt, axisMinLegPt, axisMomentumFraction, axisOccupancyTracks};
      const std::vector<AxisSpec> eligiblePairTypeAxes{axisConversionRadius, axisGammaPt, axisOccupancyTracks, axisTrackType};
      const std::vector<AxisSpec> rawV0MultiplicityAxes{axisConversionRadius, axisGammaPt, axisOccupancyTracks, axisMultiplicity};
      const std::vector<AxisSpec> rawV0TypeAxes{axisConversionRadius, axisGammaPt, axisOccupancyTracks, axisV0Type};
      const std::vector<AxisSpec> rawV0TrackTypeAxes{axisConversionRadius, axisGammaPt, axisOccupancyTracks, axisTrackType};
      const std::vector<AxisSpec> rawV0CollisionAssociationAxes{axisConversionRadius, axisGammaPt, axisOccupancyTracks, axisCollisionAssociation};
      const std::vector<AxisSpec> electronTrackMultiplicityAxes{axisConversionRadius, axisElectronPt, axisOccupancyTracks, axisMultiplicity};
      const std::vector<AxisSpec> positronTrackMultiplicityAxes{axisConversionRadius, axisPositronPt, axisOccupancyTracks, axisMultiplicity};

      registry.add("Event/hCentrality", "collision centrality", HistType::kTH1F, {axisCentrality});
      registry.add("Event/hOccupancy", "track occupancy in the time range around the collision", HistType::kTH1F, {axisOccupancyTracks});
      registry.add("Event/hCentralityOccupancy", "collision centrality versus track occupancy", HistType::kTH2F, {axisCentrality, axisOccupancyTracks});
      registry.add("Event/hNCol", "collision in Drift time", HistType::kTH1F, {axisNCol});

      registry.add("SingleLeg/hElectronEligible", "eligible true conversion electrons", HistType::kTHnSparseF, legEfficiencyAxes);
      registry.add("SingleLeg/hElectronFound", "eligible true conversion electrons with a reconstructed track", HistType::kTHnSparseF, legEfficiencyAxes);
      registry.add("SingleLeg/hPositronEligible", "eligible true conversion positrons", HistType::kTHnSparseF, legEfficiencyAxes);
      registry.add("SingleLeg/hPositronFound", "eligible true conversion positrons with a reconstructed track", HistType::kTHnSparseF, legEfficiencyAxes);

      registry.add("Partner/ElectronTag/hProbeEligible", "positron probes with a reconstructed electron tag", HistType::kTHnSparseF, partnerEfficiencyAxes);
      registry.add("Partner/ElectronTag/hProbeFound", "reconstructed positron probes with a reconstructed electron tag", HistType::kTHnSparseF, partnerEfficiencyAxes);
      registry.add("Partner/PositronTag/hProbeEligible", "electron probes with a reconstructed positron tag", HistType::kTHnSparseF, partnerEfficiencyAxes);
      registry.add("Partner/PositronTag/hProbeFound", "reconstructed electron probes with a reconstructed positron tag", HistType::kTHnSparseF, partnerEfficiencyAxes);

      registry.add("Kinematics/hGammaPtPairPt", "photon and daughter vector-sum transverse momenta", HistType::kTH2F, {axisGammaPt, axisPairPt});
      registry.add("Kinematics/hGammaPtLegPtSum", "photon transverse momentum and scalar daughter transverse-momentum sum", HistType::kTH2F, {axisGammaPt, axisLegPtSum});
      registry.add("Kinematics/hGammaPtPositronFraction", "photon transverse momentum and positron momentum fraction", HistType::kTH2F, {axisGammaPt, axisMomentumFraction});

      registry.add("Conversion/hBothLegsEligible", "fiducial true conversions with both legs eligible", HistType::kTHnSparseF, conversionEfficiencyAxes);
      registry.add("Conversion/hBothLegsFound", "fiducial true conversions with both legs reconstructed", HistType::kTHnSparseF, conversionEfficiencyAxes);
      registry.add("Conversion/hRawV0FoundLabelAnywhere", "fiducial true conversions with a label-matched raw V0 in any collision", HistType::kTHnSparseF, conversionEfficiencyAxes);
      registry.add("Conversion/hRawV0FoundAnywhere", "fiducial true conversions with a clean truth-matched raw V0 in any collision", HistType::kTHnSparseF, conversionEfficiencyAxes);
      registry.add("Conversion/hRawV0FoundSameMcCollision", "fiducial true conversions with a clean raw V0 assigned to the correct MC collision", HistType::kTHnSparseF, conversionEfficiencyAxes);
      registry.add("Conversion/hRawV0FoundSameCollision", "fiducial true conversions with a clean raw V0 assigned to the reference reconstructed collision", HistType::kTHnSparseF, conversionEfficiencyAxes);
      registry.add("Conversion/hRawV0FoundUsingEligiblePair", "fiducial true conversions with a clean raw V0 using a denominator-eligible track pair", HistType::kTHnSparseF, conversionEfficiencyAxes);
      registry.add("Conversion/hRawV0FoundSameCollisionUsingEligiblePair", "fiducial true conversions with a clean raw V0 in the reference collision using a denominator-eligible track pair", HistType::kTHnSparseF, conversionEfficiencyAxes);
      registry.add("Conversion/hV0Found", "fiducial true conversions reconstructed as a raw V0 in the reference collision", HistType::kTHnSparseF, conversionEfficiencyAxes);

      registry.add("ConversionDaughterKinematics/hBothLegsFound", "fiducial true conversions with both legs reconstructed", HistType::kTHnSparseF, conversionDaughterKinematicsAxes);
      registry.add("ConversionDaughterKinematics/hRawV0FoundAnywhere", "fiducial true conversions with a clean truth-matched raw V0 in any collision", HistType::kTHnSparseF, conversionDaughterKinematicsAxes);

      registry.add("EligiblePairType/hBothLegsFound", "denominator-eligible reconstructed track-pair types", HistType::kTHnSparseF, eligiblePairTypeAxes);
      registry.add("EligiblePairType/hRawV0FoundUsingEligiblePair", "denominator-eligible track-pair types used by a clean raw V0", HistType::kTHnSparseF, eligiblePairTypeAxes);

      registry.add("RawV0/hCandidateMultiplicity", "number of clean raw V0 rows per reconstructable truth conversion", HistType::kTHnSparseF, rawV0MultiplicityAxes);
      registry.add("RawV0/hDistinctCollisionMultiplicity", "number of distinct raw V0 collision assignments per reconstructable truth conversion", HistType::kTHnSparseF, rawV0MultiplicityAxes);
      registry.add("RawV0/hDistinctTrackPairMultiplicity", "number of distinct raw V0 daughter track pairs per reconstructable truth conversion", HistType::kTHnSparseF, rawV0MultiplicityAxes);
      registry.add("RawV0/hV0Type", "type of each clean truth-matched raw V0 row", HistType::kTHnSparseF, rawV0TypeAxes);
      registry.add("RawV0/hTrackType", "daughter track type of each clean truth-matched raw V0 row", HistType::kTHnSparseF, rawV0TrackTypeAxes);
      registry.add("RawV0/hCollisionAssociation", "collision association of each clean truth-matched raw V0 row", HistType::kTHnSparseF, rawV0CollisionAssociationAxes);

      registry.add("TrackClones/hElectronSameCollision", "number of reconstructed electron tracks in the reference collision", HistType::kTHnSparseF, electronTrackMultiplicityAxes);
      registry.add("TrackClones/hElectronEligible", "number of denominator-eligible electron tracks in the reference collision", HistType::kTHnSparseF, electronTrackMultiplicityAxes);
      registry.add("TrackClones/hPositronSameCollision", "number of reconstructed positron tracks in the reference collision", HistType::kTHnSparseF, positronTrackMultiplicityAxes);
      registry.add("TrackClones/hPositronEligible", "number of denominator-eligible positron tracks in the reference collision", HistType::kTHnSparseF, positronTrackMultiplicityAxes);
    }
  }

  template <typename TTrack>
  bool passesReconstructedLegSelection(TTrack const& track, int expectedSign)
  {
    if (requireCorrectMcTrackLabel && track.mcMask() != 0) {
      return false;
    }
    if (track.sign() * expectedSign <= 0) {
      return false;
    }
    if (requireTPCForLeg && !track.hasTPC()) {
      return false;
    }
    return track.hasITS() || track.hasTPC();
  }

  template <typename TMCPhoton>
  bool passesGeneratedPhotonSelection(TMCPhoton const& photon)
  {
    if (requirePrimaryOrGeneratorGamma && !photon.isPhysicalPrimary() && !photon.producedByGenerator()) {
      return false;
    }
    if (photon.pt() < minPtGamma || photon.pt() > maxPtGamma) {
      return false;
    }
    return std::abs(photon.eta()) <= maxEtaGamma;
  }

  template <typename TMCPhoton>
  bool isInPhotonFiducialVolume(TMCPhoton const& photon, float conversionRadius, float conversionZ)
  {
    if (!passesGeneratedPhotonSelection(photon)) {
      return false;
    }
    if (conversionRadius < minConversionRadius || conversionRadius > maxConversionRadius) {
      return false;
    }
    if (maxAbsConversionZ >= 0.f && std::abs(conversionZ) > maxAbsConversionZ) {
      return false;
    }
    if (applyRZLineCut && conversionRadius < (std::abs(conversionZ) * std::tan(2.f * std::atan(std::exp(-maxEtaGamma)))) - marginZMC) {
      return false;
    }
    return true;
  }

  template <typename TMCParticle>
  ReconstructedTrackSummary getReconstructedTrackSummary(TMCParticle const& particle, int expectedSign, int64_t collisionId)
  {
    ReconstructedTrackSummary summary;
    for (const auto& track : particle.template tracks_as<TracksMC>()) {
      if (track.collisionId() != collisionId) {
        continue;
      }
      ++summary.nSameCollision;
      if (passesReconstructedLegSelection(track, expectedSign)) {
        summary.eligibleTrackIds.push_back(track.globalIndex());
      }
    }
    return summary;
  }

  template <typename TTrack>
  int getRawV0TrackType(TTrack const& pos, TTrack const& ele)
  {
    const bool posITSTPC = pos.hasITS() && pos.hasTPC();
    const bool eleITSTPC = ele.hasITS() && ele.hasTPC();
    const bool posTPCOnly = !pos.hasITS() && pos.hasTPC();
    const bool eleTPCOnly = !ele.hasITS() && ele.hasTPC();
    const bool posITSOnly = pos.hasITS() && !pos.hasTPC();
    const bool eleITSOnly = ele.hasITS() && !ele.hasTPC();

    if (posITSTPC && eleITSTPC) {
      return kITSTPCITSTPC;
    }
    if ((posITSTPC && eleTPCOnly) || (posTPCOnly && eleITSTPC)) {
      return kITSTPCTPCOnly;
    }
    if (posTPCOnly && eleTPCOnly) {
      return kTPCOnlyTPCOnly;
    }
    if (posITSOnly || eleITSOnly) {
      return kContainsITSOnly;
    }
    return kOther;
  }

  std::vector<int> countNeighbors(CollisionsMC const& collisions)
  {
    struct ValidCollision {
      int64_t bc;
      int64_t originalIndex;
    };
    std::vector<ValidCollision> valid;
    valid.reserve(collisions.size());
    std::vector<int> counts(collisions.size(), -1);
    for (int64_t i = 0; i < collisions.size(); ++i) {
      const auto& collision = collisions.rawIteratorAt(i);
      if (!collision.has_foundBC() || collision.numContrib() < minNCont) {
        continue;
      }
      valid.push_back({static_cast<int64_t>(collision.foundBC_as<aod::BCs>().globalBC()), i});
    }
    size_t left = 0;
    size_t right = 0;
    for (size_t i = 0; i < valid.size(); ++i) {
      const int64_t currentBC = valid[i].bc;
      while (valid[left].bc < currentBC - maxTPCDriftTime) {
        ++left;
      }
      while (right < valid.size() && valid[right].bc <= currentBC + maxTPCDriftTime) {
        ++right;
      }
      counts[valid[i].originalIndex] = static_cast<int>(right - left - 1);
    }
    return counts;
  }

  void processTrueConversion(FullMcParticles const& mcParticles, TracksMC const& /*tracks*/)
  {
    for (const auto& mcPhoton : mcParticles) {
      if (mcPhoton.pdgCode() != kGamma) {
        continue;
      }
      if (!passesGeneratedPhotonSelection(mcPhoton)) {
        continue;
      }

      registry.fill(HIST("ConversionProbability/hGeneratedGammaCount"), 1.f);
      registry.fill(HIST("ConversionProbability/hGeneratedGammaPt"), mcPhoton.pt());
      if (!mcPhoton.has_daughters()) {
        continue;
      }

      bool hasPositron{false};
      bool hasElectron{false};
      float conversionX{0.f};
      float conversionY{0.f};
      float conversionZ{0.f};

      for (const auto& daughter : mcPhoton.template daughters_as<FullMcParticles>()) {
        if (daughter.getProcess() != TMCProcess::kPPair) {
          continue;
        }
        if (daughter.pdgCode() == kPositron) {
          if (!hasPositron && !hasElectron) {
            conversionX = daughter.vx();
            conversionY = daughter.vy();
            conversionZ = daughter.vz();
          }
          hasPositron = true;
        } else if (daughter.pdgCode() == kElectron) {
          if (!hasPositron && !hasElectron) {
            conversionX = daughter.vx();
            conversionY = daughter.vy();
            conversionZ = daughter.vz();
          }
          hasElectron = true;
        }
      }

      if (!hasPositron || !hasElectron) {
        continue;
      }

      const float conversionRadius = std::hypot(conversionX, conversionY);
      if (!isInPhotonFiducialVolume(mcPhoton, conversionRadius, conversionZ)) {
        continue;
      }

      registry.fill(HIST("ConversionProbability/hConvertedGammaPt"), mcPhoton.pt());
      registry.fill(HIST("ConversionProbability/hConversionRadius"), conversionRadius);
      registry.fill(HIST("ConversionProbability/hRPhiTrueConversion"), conversionRadius, mcPhoton.phi());
      registry.fill(HIST("ConversionProbability/hRZTrueConversion"), conversionRadius, conversionZ);
    }
  }
  PROCESS_SWITCH(CheckMcV0, processTrueConversion, "process true photon conversions", true);

  Preslice<FullMcParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  void processMCV0(CollisionsMC const& collisions, aod::BCs const& /*bcs*/, aod::V0s const& rawV0s, FullMcParticles const& mcParticles, aod::McCollisions const& /*mcCollisions*/, TracksMC const& tracks)
  {
    // Index raw production V0 rows by their truth photon. Efficiency numerators
    // below are filled once per truth-photon/reference-collision entry, while
    // this vector retains all cloned rows for multiplicity diagnostics.
    std::unordered_map<int64_t, std::vector<RawV0Match>> rawV0MatchesByPhoton;
    for (const auto& rawV0 : rawV0s) {
      const auto& pos = rawV0.template posTrack_as<TracksMC>();
      const auto& ele = rawV0.template negTrack_as<TracksMC>();
      if (!pos.has_mcParticle() || !ele.has_mcParticle()) {
        continue;
      }

      const auto& posMC = pos.template mcParticle_as<FullMcParticles>();
      const auto& eleMC = ele.template mcParticle_as<FullMcParticles>();
      if (posMC.globalIndex() == eleMC.globalIndex() || posMC.pdgCode() != kPositron || eleMC.pdgCode() != kElectron) {
        continue;
      }
      if (posMC.getProcess() != TMCProcess::kPPair || eleMC.getProcess() != TMCProcess::kPPair || !checkMCParticles<kGamma>(posMC, eleMC)) {
        continue;
      }

      const auto& mother = posMC.template mothers_first_as<FullMcParticles>();
      int64_t assignedMcCollisionId{-1};
      if (rawV0.collisionId() >= 0) {
        const auto& assignedCollision = rawV0.template collision_as<CollisionsMC>();
        if (assignedCollision.has_mcCollision()) {
          assignedMcCollisionId = assignedCollision.mcCollisionId();
        }
      }

      const bool hasCleanLabels = !requireCorrectMcTrackLabel || (pos.mcMask() == 0 && ele.mcMask() == 0);
      rawV0MatchesByPhoton[mother.globalIndex()].push_back({rawV0.collisionId(),
                                                            assignedMcCollisionId,
                                                            pos.globalIndex(),
                                                            ele.globalIndex(),
                                                            rawV0.v0Type(),
                                                            getRawV0TrackType(pos, ele),
                                                            hasCleanLabels});
    }

    const auto nCols = countNeighbors(collisions);
    for (int64_t iCol = 0; iCol < collisions.size(); ++iCol) {
      const auto& collision = collisions.rawIteratorAt(iCol);
      if (!collision.has_mcCollision()) {
        continue;
      }

      const float centrality = [&] {
        switch (cfgCentEstimator) {
          case 2:
            return collision.centFT0C();
          case 1:
            return collision.centFT0A();
          default:
            return collision.centFT0M();
        }
      }();
      const auto nCol = nCols[iCol];
      const auto occupancy = static_cast<float>(collision.trackOccupancyInTimeRange());
      registry.fill(HIST("Event/hCentrality"), centrality);
      registry.fill(HIST("Event/hOccupancy"), occupancy);
      registry.fill(HIST("Event/hCentralityOccupancy"), centrality, occupancy);
      registry.fill(HIST("Event/hNCol"), nCol);

      const auto mcParticlesThisCollision = mcParticles.sliceBy(perMcCollision, collision.mcCollisionId());
      for (const auto& mcPhoton : mcParticlesThisCollision) {
        if (mcPhoton.pdgCode() != kGamma || !mcPhoton.has_daughters()) {
          continue;
        }

        int64_t positronId{-1};
        int64_t electronId{-1};
        float conversionX{0.f};
        float conversionY{0.f};
        float conversionZ{0.f};
        for (const auto& daughter : mcPhoton.template daughters_as<FullMcParticles>()) {
          if (daughter.getProcess() != TMCProcess::kPPair) {
            continue;
          }
          if (daughter.pdgCode() == kPositron && positronId < 0) {
            positronId = daughter.globalIndex();
            conversionX = daughter.vx();
            conversionY = daughter.vy();
            conversionZ = daughter.vz();
          } else if (daughter.pdgCode() == kElectron && electronId < 0) {
            electronId = daughter.globalIndex();
            conversionX = daughter.vx();
            conversionY = daughter.vy();
            conversionZ = daughter.vz();
          }
        }

        if (positronId < 0 || electronId < 0) {
          continue;
        }

        const float conversionRadius = std::hypot(conversionX, conversionY);
        if (!isInPhotonFiducialVolume(mcPhoton, conversionRadius, conversionZ)) {
          continue;
        }

        const auto& positron = mcParticles.iteratorAt(positronId);
        const auto& electron = mcParticles.iteratorAt(electronId);
        const float legPtSum = electron.pt() + positron.pt();
        const float minLegPt = std::min(electron.pt(), positron.pt());
        const float pairPt = std::hypot(electron.px() + positron.px(), electron.py() + positron.py());
        const float positronFraction = legPtSum > 0.f ? positron.pt() / legPtSum : 0.f;
        registry.fill(HIST("Kinematics/hGammaPtPairPt"), mcPhoton.pt(), pairPt);
        registry.fill(HIST("Kinematics/hGammaPtLegPtSum"), mcPhoton.pt(), legPtSum);
        registry.fill(HIST("Kinematics/hGammaPtPositronFraction"), mcPhoton.pt(), positronFraction);

        const auto positronTracks = getReconstructedTrackSummary(positron, +1, collision.globalIndex());
        const auto electronTracks = getReconstructedTrackSummary(electron, -1, collision.globalIndex());
        const bool positronFound = !positronTracks.eligibleTrackIds.empty();
        const bool electronFound = !electronTracks.eligibleTrackIds.empty();

        registry.fill(HIST("TrackClones/hElectronSameCollision"), conversionRadius, electron.pt(), occupancy, electronTracks.nSameCollision);
        registry.fill(HIST("TrackClones/hElectronEligible"), conversionRadius, electron.pt(), occupancy, electronTracks.eligibleTrackIds.size());
        registry.fill(HIST("TrackClones/hPositronSameCollision"), conversionRadius, positron.pt(), occupancy, positronTracks.nSameCollision);
        registry.fill(HIST("TrackClones/hPositronEligible"), conversionRadius, positron.pt(), occupancy, positronTracks.eligibleTrackIds.size());

        registry.fill(HIST("SingleLeg/hElectronEligible"), conversionRadius, electron.pt(), occupancy, centrality, nCol);
        registry.fill(HIST("SingleLeg/hPositronEligible"), conversionRadius, positron.pt(), occupancy, centrality, nCol);
        if (electronFound) {
          registry.fill(HIST("SingleLeg/hElectronFound"), conversionRadius, electron.pt(), occupancy, centrality, nCol);
          registry.fill(HIST("Partner/ElectronTag/hProbeEligible"), conversionRadius, electron.pt(), positron.pt(), occupancy, centrality, nCol);
          if (positronFound) {
            registry.fill(HIST("Partner/ElectronTag/hProbeFound"), conversionRadius, electron.pt(), positron.pt(), occupancy, centrality, nCol);
          }
        }
        if (positronFound) {
          registry.fill(HIST("SingleLeg/hPositronFound"), conversionRadius, positron.pt(), occupancy, centrality, nCol);
          registry.fill(HIST("Partner/PositronTag/hProbeEligible"), conversionRadius, electron.pt(), positron.pt(), occupancy, centrality, nCol);
          if (electronFound) {
            registry.fill(HIST("Partner/PositronTag/hProbeFound"), conversionRadius, electron.pt(), positron.pt(), occupancy, centrality, nCol);
          }
        }

        registry.fill(HIST("Conversion/hBothLegsEligible"), conversionRadius, mcPhoton.pt(), occupancy, centrality, nCol);
        if (!electronFound || !positronFound) {
          continue;
        }

        registry.fill(HIST("Conversion/hBothLegsFound"), conversionRadius, mcPhoton.pt(), occupancy, centrality, nCol);
        registry.fill(HIST("ConversionDaughterKinematics/hBothLegsFound"), conversionRadius, mcPhoton.pt(), minLegPt, positronFraction, occupancy);

        std::array<bool, kOther + 1> hasEligiblePairType{};
        for (const auto positronTrackId : positronTracks.eligibleTrackIds) {
          const auto& positronTrack = tracks.iteratorAt(positronTrackId);
          for (const auto electronTrackId : electronTracks.eligibleTrackIds) {
            const auto& electronTrack = tracks.iteratorAt(electronTrackId);
            hasEligiblePairType[getRawV0TrackType(positronTrack, electronTrack)] = true;
          }
        }
        for (int trackType = 0; trackType <= kOther; ++trackType) {
          if (hasEligiblePairType[trackType]) {
            registry.fill(HIST("EligiblePairType/hBothLegsFound"), conversionRadius, mcPhoton.pt(), occupancy, trackType);
          }
        }

        const auto rawV0MatchesIt = rawV0MatchesByPhoton.find(mcPhoton.globalIndex());
        const bool hasLabelMatchAnywhere = rawV0MatchesIt != rawV0MatchesByPhoton.end() && !rawV0MatchesIt->second.empty();
        if (hasLabelMatchAnywhere) {
          registry.fill(HIST("Conversion/hRawV0FoundLabelAnywhere"), conversionRadius, mcPhoton.pt(), occupancy, centrality, nCol);
        }

        bool hasCleanMatchAnywhere{false};
        bool hasCleanMatchSameMcCollision{false};
        bool hasCleanMatchSameCollision{false};
        bool hasCleanMatchUsingEligiblePair{false};
        bool hasCleanMatchSameCollisionUsingEligiblePair{false};
        int nCleanMatches{0};
        std::array<bool, kOther + 1> hasRawV0UsingEligiblePairType{};
        std::unordered_set<int64_t> distinctCollisionIds;
        std::set<std::pair<int64_t, int64_t>> distinctTrackPairs;
        const std::unordered_set<int64_t> eligiblePositronTrackIds{positronTracks.eligibleTrackIds.begin(), positronTracks.eligibleTrackIds.end()};
        const std::unordered_set<int64_t> eligibleElectronTrackIds{electronTracks.eligibleTrackIds.begin(), electronTracks.eligibleTrackIds.end()};

        if (rawV0MatchesIt != rawV0MatchesByPhoton.end()) {
          for (const auto& rawV0Match : rawV0MatchesIt->second) {
            if (!rawV0Match.hasCleanLabels) {
              continue;
            }
            hasCleanMatchAnywhere = true;
            ++nCleanMatches;
            distinctCollisionIds.insert(rawV0Match.collisionId);
            distinctTrackPairs.emplace(rawV0Match.posTrackId, rawV0Match.negTrackId);
            registry.fill(HIST("RawV0/hV0Type"), conversionRadius, mcPhoton.pt(), occupancy, rawV0Match.v0Type);
            registry.fill(HIST("RawV0/hTrackType"), conversionRadius, mcPhoton.pt(), occupancy, rawV0Match.trackType);

            int collisionAssociation = kOtherMcCollision;
            if (rawV0Match.collisionId == collision.globalIndex()) {
              collisionAssociation = kExactRecoCollision;
            } else if (rawV0Match.mcCollisionId < 0) {
              collisionAssociation = kNoMcCollisionLabel;
            } else if (rawV0Match.mcCollisionId == collision.mcCollisionId()) {
              collisionAssociation = kOtherRecoSameMcCollision;
            }
            registry.fill(HIST("RawV0/hCollisionAssociation"), conversionRadius, mcPhoton.pt(), occupancy, collisionAssociation);

            const bool usesEligiblePair = eligiblePositronTrackIds.contains(rawV0Match.posTrackId) && eligibleElectronTrackIds.contains(rawV0Match.negTrackId);
            if (usesEligiblePair) {
              hasCleanMatchUsingEligiblePair = true;
              hasRawV0UsingEligiblePairType[rawV0Match.trackType] = true;
            }
            if (rawV0Match.mcCollisionId == collision.mcCollisionId()) {
              hasCleanMatchSameMcCollision = true;
            }
            if (rawV0Match.collisionId == collision.globalIndex()) {
              hasCleanMatchSameCollision = true;
              if (usesEligiblePair) {
                hasCleanMatchSameCollisionUsingEligiblePair = true;
              }
            }
          }
        }

        registry.fill(HIST("RawV0/hCandidateMultiplicity"), conversionRadius, mcPhoton.pt(), occupancy, nCleanMatches);
        registry.fill(HIST("RawV0/hDistinctCollisionMultiplicity"), conversionRadius, mcPhoton.pt(), occupancy, distinctCollisionIds.size());
        registry.fill(HIST("RawV0/hDistinctTrackPairMultiplicity"), conversionRadius, mcPhoton.pt(), occupancy, distinctTrackPairs.size());

        if (hasCleanMatchAnywhere) {
          registry.fill(HIST("Conversion/hRawV0FoundAnywhere"), conversionRadius, mcPhoton.pt(), occupancy, centrality, nCol);
          registry.fill(HIST("ConversionDaughterKinematics/hRawV0FoundAnywhere"), conversionRadius, mcPhoton.pt(), minLegPt, positronFraction, occupancy);
        }
        if (hasCleanMatchSameMcCollision) {
          registry.fill(HIST("Conversion/hRawV0FoundSameMcCollision"), conversionRadius, mcPhoton.pt(), occupancy, centrality, nCol);
        }
        if (hasCleanMatchUsingEligiblePair) {
          registry.fill(HIST("Conversion/hRawV0FoundUsingEligiblePair"), conversionRadius, mcPhoton.pt(), occupancy, centrality, nCol);
        }
        if (hasCleanMatchSameCollisionUsingEligiblePair) {
          registry.fill(HIST("Conversion/hRawV0FoundSameCollisionUsingEligiblePair"), conversionRadius, mcPhoton.pt(), occupancy, centrality, nCol);
        }
        if (hasCleanMatchSameCollision) {
          registry.fill(HIST("Conversion/hRawV0FoundSameCollision"), conversionRadius, mcPhoton.pt(), occupancy, centrality, nCol);
          registry.fill(HIST("Conversion/hV0Found"), conversionRadius, mcPhoton.pt(), occupancy, centrality, nCol);
        }
        for (int trackType = 0; trackType <= kOther; ++trackType) {
          if (hasRawV0UsingEligiblePairType[trackType]) {
            registry.fill(HIST("EligiblePairType/hRawV0FoundUsingEligiblePair"), conversionRadius, mcPhoton.pt(), occupancy, trackType);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(CheckMcV0, processMCV0, "process truth-matched reconstructed photon conversions", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<CheckMcV0>(cfgc)};
}
