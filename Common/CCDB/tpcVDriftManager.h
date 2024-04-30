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
#include "ReconstructionDataFormats/TrackParametrization.h"

namespace o2::aod::common
{

// Thin wrapper for vdrift ccdb queries
// should eventually mirror VDriftHelper class
class TPCVDriftManager
{
 public:
  void update(long timestamp)
  {
    // Check validity of already present obj, otherwise update
    if (mVD != nullptr && (timestamp > mVD->firstTime || timestamp < mVD->lastTime)) {
      return;
    }

    // Update Obj
    mVD = mCCDB->getForTimeStamp<o2::tpc::VDriftCorrFact>(mVDriftTglPath, timestamp);

    // Update factors
    mTPCVDrift = mVD->refVDrift * mVD->corrFact;
    mTPCVDriftCorrFact = mVD->corrFact;
    mTPCBin2Z = mTPCVDrift / mMUS2TPCBin;

    LOGP(info, "Updated VDrift for {} with vdrift={} and bin2z={}", mVD->creationTime, mTPCVDrift, mTPCBin2Z);
  }

  [[nodiscard]] auto getVDrift() const -> float
  {
    return mTPCVDrift;
  }

  [[nodiscard]] auto getVDriftErr() const -> float
  {
    return mTPCVDriftCorrFact;
  };

  [[nodiscard]] auto getTPCBin2Z() const -> float
  {
    return mTPCBin2Z;
  };

  static constexpr float mMUS2TPCBin{1.f / (8 * o2::constants::lhc::LHCBunchSpacingMUS)};

  template <typename Collision, typename TrackExtra, typename Track>
  [[nodiscard]] auto correctTPCTrack(const Collision& col, const TrackExtra& trackExtra, Track& track) const -> bool
  {
      // TODO Does this Track need a correction?
      float tTB, tTBErr;
      if(col.collisionTimeRes() < 0){
          tTB = trackExtra.tpcTime0();
          tTBErr = trackExtra.trackTimeRes();
      } else{
          tTB = col.collisionTime() * mMUS2TPCBin;
          tTBErr = col.collisionTimeRes() * mMUS2TPCBin;
      }

      float dDrift = (tTB - trackExtra.tpcTime0()) * mTPCBin2Z;
      float dDriftErr = tTBErr * mTPCBin2Z;
      if(dDriftErr < 0.f){
          return false;
      }
      // TODO how to check for constrained tracks?
      track.setZ(track.getZ() + ((track.getZ() < 0.) ? -dDrift : dDrift));
      /* track.setCov(track.getSigmaZ2() + dDriftErr*dDriftErr, o2::track::kSigZ2); */

      return true;
  }

 private:
  float mTPCBin2Z{};
  float mTPCVDrift{};
  float mTPCVDriftCorrFact{};

  const o2::tpc::VDriftCorrFact* mVD{nullptr};
  std::string mVDriftTglPath{"TPC/Calib/VDriftTgl"};
  o2::framework::Service<o2::ccdb::BasicCCDBManager> mCCDB;
};

} // namespace o2::aod::common

#endif // TPCVDRIFTMANAGER_H_
