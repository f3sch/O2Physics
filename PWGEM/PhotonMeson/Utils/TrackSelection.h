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

/// \brief helper functions for pair track selection
/// \author felix.schlepper@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_TRACKSELECTION_H_
#define PWGEM_PHOTONMESON_UTILS_TRACKSELECTION_H_

namespace o2::pwgem::photonmeson
{

/**
 * @brief Track has ITS and TPC
 *
 * @tparam TTrack track
 * @param track track
 * @return true if has both
 */
template <typename TTrack>
bool isITSTPCMatchedTrack(TTrack const& track)
{
  return track.hasITS() && track.hasTPC();
}

/**
 * @brief Track is TPC-only
 *
 * @tparam TTrack track
 * @param track track
 * @return true if tracks is TPC-only
 */
template <typename TTrack>
bool isTPConlyTrack(TTrack const& track)
{
  return !track.hasITS() && track.hasTPC();
}

/**
 * @brief Track is ITS-only
 *
 * @tparam TTrack track
 * @param track track
 * @return true if tracks is ITS-only
 */
template <typename TTrack>
bool isITSonlyTrack(TTrack const& track)
{
  return track.hasITS() && !track.hasTPC();
}

/**
 * @brief If both V0 pairs are ITSTPC-tracks
 *
 * @tparam TTrack track
 * @param track0 track from daughter 0
 * @param track1 track from daughter 1
 * @return true if V0 pairs are ITSTPC-tracks
 */
template <typename TTrack>
bool isITSTPC_ITSTPC(TTrack const& track0, TTrack const& track1)
{
  return isITSTPCMatchedTrack(track0) && isITSTPCMatchedTrack(track1);
}

/**
 * @brief If one track is TPC-only the other ITSTPC
 *
 * @tparam TTrack track
 * @param track0 track from daughter 0
 * @param track1 track from daughter 1
 * @return true if one is TPC-only and the other ITSTPC
 */
template <typename TTrack>
bool isITSTPC_TPConly(TTrack const& track0, TTrack const& track1)
{
  return (isITSTPCMatchedTrack(track0) && isTPConlyTrack(track1)) || (isITSTPCMatchedTrack(track1) && isTPConlyTrack(track0));
}

/**
 * @brief If one track is ITS-only the other ITSTPC
 *
 * @tparam TTrack track
 * @param track0 track from daughter 0
 * @param track1 track from daughter 1
 * @return true if one is ITS-only and the other ITSTPC
 */
template <typename TTrack>
bool isITSTPC_ITSonly(TTrack const& track0, TTrack const& track1)
{
  return (isITSTPCMatchedTrack(track0) && isITSonlyTrack(track1)) || (isITSTPCMatchedTrack(track1) && isITSonlyTrack(track0));
}

/**
 * @brief If V0 pairs are TPC-only tracks
 *
 * @tparam TTrack track
 * @param track0 track from daughter 0
 * @param track1 track from daughter 1
 * @return true if both are TPC-only tracks
 */
template <typename TTrack>
bool isTPConly_TPConly(TTrack const& track0, TTrack const& track1)
{
  return isTPConlyTrack(track0) && isTPConlyTrack(track1);
}

/**
 * @brief If V0 pairs are ITS-only tracks
 *
 * @tparam TTrack track
 * @param track0 track from daughter 0
 * @param track1 track from daughter 1
 * @return true if both are ITS-only tracks
 */
template <typename TTrack>
bool isITSonly_ITSonly(TTrack const& track0, TTrack const& track1)
{
  return isITSonlyTrack(track0) && isITSonlyTrack(track1);
}

/**
 * @brief If one V0 pair is ITS-only and the other TPC-only
 *
 * @tparam TTrack track
 * @param track0 track from daughter 0
 * @param track1 track from daughter 1
 * @return true if either one is ITS-only while the other one is TPC-only
 */
template <typename TTrack>
bool isTPConly_ITSonly(TTrack const& track0, TTrack const& track1)
{
  return (isTPConlyTrack(track0) && isITSonlyTrack(track1)) || (isTPConlyTrack(track1) && isITSonlyTrack(track0));
}
} // namespace o2::pwgem::photonmeson

#endif // PWGEM_PHOTONMESON_UTILS_TRACKSELECTION_H_
