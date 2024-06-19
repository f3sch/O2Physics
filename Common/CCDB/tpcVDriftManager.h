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
#include "Framework/Logger.h"
#include "Framework/DataTypes.h"
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
  void init(o2::ccdb::BasicCCDBManager* ccdb)
  {
    mCCDB = ccdb;
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

    // TODO Laser calib

    // Update factors
    mTPCVDrift = mVD->refVDrift * mVD->corrFact;
    mTPCVDriftCorrFact = mVD->corrFact;
    mTPCVDriftNS = mTPCVDrift * 1e-3;

    LOGP(info, "Updated VDrift for timestamp {} with vdrift={:.4f} (cm/us)", mVD->creationTime, mTPCVDrift);
  }

  template <typename Collision, typename TrackExtra, typename Track>
  [[nodiscard]] auto moveTPCTrack(const Collision& col, const TrackExtra& trackExtra, Track& track) -> bool
  {
    // track is fine, or cannot be moved has information is not available
    if (!(trackExtra.flags() & o2::aod::track::TrackFlags::TrackTimeAsym)) {
      return true;
    }

    // Check if there is a good object available otherwise pretend everything is fine
    if (mVD == nullptr) {
      LOGP(warn, "No VDrift object available, pretending track to be correct");
      return true;
    }

    // TPC time is given relative to the closest BC in ns
    float tTB, tTBErr;
    if (col.collisionTimeRes() < 0.f) { // use track data
      tTB = trackExtra.trackTime();
      o2::aod::track::extensions::TPCTimeErrEncoding enc;
      enc.encoding.timeErr = trackExtra.trackTimeRes();
      tTBErr = 0.5f * (enc.getDeltaTFwd() + enc.getDeltaTBwd());
    } else {
      tTB = col.collisionTime();
      tTBErr = col.collisionTimeRes();
    }

    float dDrift = (tTB - trackExtra.trackTime()) * mTPCVDriftNS;
    float dDriftErr = tTBErr * mTPCVDriftNS;
    if (dDriftErr < 0.f || dDrift > 250.f) { // we cannot move a track outside the drift volume
      LOGP(warn, "Skipping faulty correction with dDrift={} +- {}", dDrift, dDriftErr);
      return false;
    }
    // TODO how to check for constrained tracks?
    track.setZ(track.getZ() + ((track.getTgl() < 0.) ? -dDrift : dDrift));
    track.setCov(track.getSigmaZ2() + dDriftErr * dDriftErr, o2::track::kSigZ2);

    return true;
  }

 private:
  float mTPCVDrift{};
  float mTPCVDriftNS{};
  float mTPCVDriftCorrFact{};

  const o2::tpc::VDriftCorrFact* mVD{nullptr};
  o2::ccdb::BasicCCDBManager* mCCDB;
};

} // namespace o2::aod::common

#endif // TPCVDRIFTMANAGER_H_
