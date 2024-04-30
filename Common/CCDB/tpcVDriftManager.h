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

#ifndef TPCVDRIFTMANAGER_H_
#define TPCVDRIFTMANAGER_H_

#include <string>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/Logger.h"
#include "DataFormatsTPC/VDriftCorrFact.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

namespace o2::aod::common
{

// Thin wrapper for vdrift ccdb queries
// should eventually mirror VDriftHelper class
class TPCVDriftManager
{
 public:
  void init(o2::ccdb::BasicCCDBManager* ccdb, int run, long timestamp)
  {
    mCCDB = ccdb;

    /// Copied from eventselection
    // access orbit-reset timestamp
    const auto ctpx = mCCDB->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", timestamp);
    int64_t tsOrbitReset = ctpx->at(0); // us
    // access TF duration, start-of-run and end-of-run timestamps from ECS GRP
    std::map<std::string, std::string> metadata;
    metadata["runNumber"] = Form("%d", run);
    const auto grpecs = mCCDB->getSpecific<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", timestamp, metadata);
    uint32_t nOrbitsPerTF = grpecs->getNHBFPerTF(); // assuming 1 orbit = 1 HBF;  nOrbitsPerTF=128 in 2022, 32 in 2023
    int64_t tsSOR = grpecs->getTimeStart();         // ms
    // calculate SOR orbit
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    // adjust to the nearest TF edge
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    // first bc of the first orbit (should coincide with TF start)
    mBCSOR = orbitSOR * o2::constants::lhc::LHCMaxBunches;
    // duration of TF in bcs
    mNBCsPerTF = nOrbitsPerTF * o2::constants::lhc::LHCMaxBunches;
    LOGP(info, "tsOrbitReset={} us, SOR = {} ms, orbitSOR = {}, nBCsPerTF = {}", tsOrbitReset, tsSOR, orbitSOR, mNBCsPerTF);
  }

  void update(long timestamp)
  {
    // Check validity of already present obj, otherwise update
    if (mVD != nullptr && (timestamp > mVD->firstTime || timestamp < mVD->lastTime)) {
      return;
    }

    // Update Obj
    mVD = mCCDB->getForTimeStamp<o2::tpc::VDriftCorrFact>("TPC/Calib/VDriftTgl", timestamp);
    if (mVD == nullptr) {
      LOGP(error, "Got nullptr from ccdb for VDriftCorrFact for {}", timestamp);
      return;
    }

    // Update factors
    mTPCVDrift = mVD->refVDrift * mVD->corrFact;
    mTPCVDriftCorrFact = mVD->corrFact;
    mTPCBin2Z = mTPCVDrift / mMUS2TPCBin;

    LOGP(info, "Updated VDrift for {} with vdrift={} and bin2z={}", mVD->creationTime, mTPCVDrift, mTPCBin2Z);
  }

  template <typename BC, typename Collision, typename TrackExtra, typename Track>
  [[nodiscard]] auto correctTPCTrack(BC bc, const Collision& col, const TrackExtra& trackExtra, Track& track) -> bool
  {
    // Check if there is a good object available otherwise pretend everything is fine
    if (mVD == nullptr) {
      LOGP(debug, "No VDrift object available, pretending track to be correct");
      return true;
    }

    // TPC time is given relative to the TF start
    int64_t bcInTF = (bc.globalBC() - mBCSOR) % mNBCsPerTF;                                       // BC in TF
    float timeCol = bcInTF * o2::constants::lhc::LHCBunchSpacingMUS + col.collisionTime() * 1e-3; // Collision time is given relative to the bc [ns]
    float tTB = timeCol * mMUS2TPCBin;
    float tTBErr = col.collisionTimeRes() * mMUS2TPCBin;

    float dDrift = (tTB - trackExtra.tpcTime0()) * mTPCBin2Z;
    float dDriftErr = tTBErr * mTPCBin2Z;
    if (dDriftErr < 0.f || dDrift > 150.f) { // Generous cut on Z correction
      LOGP(warn, "Skipping faulty correction with dDrift={} +- {}", dDrift, dDriftErr);
      return false;
    }
    // TODO how to check for constrained tracks?
    track.setZ(track.getZ() + ((track.getZ() < 0.) ? -dDrift : dDrift));
    track.setCov(track.getSigmaZ2() + dDriftErr*dDriftErr, o2::track::kSigZ2);

    return true;
  }

 private:
  static constexpr float mMUS2TPCBin{1.f / (8 * o2::constants::lhc::LHCBunchSpacingMUS)};
  float mTPCBin2Z{};
  float mTPCVDrift{};
  float mTPCVDriftCorrFact{};
  int64_t mBCSOR = -1;     // global bc of the start of the first orbit
  int64_t mNBCsPerTF = -1; // duration of TF in bcs, should be 128*3564 or 32*3564

  const o2::tpc::VDriftCorrFact* mVD{nullptr};
  o2::ccdb::BasicCCDBManager* mCCDB;
};

} // namespace o2::aod::common

#endif // TPCVDRIFTMANAGER_H_
