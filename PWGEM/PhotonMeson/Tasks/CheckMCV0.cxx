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

/// \brief check the v0 phase-space
/// dependencies: o2-analysis-lf-lambdakzeromcfinder
/// \author daiki.sekihata@cern.ch felix.schlepper@cern.ch

#include "TTree.h"
#include "TDatabasePDG.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "ReconstructionDataFormats/Track.h"
#include "DetectorsBase/Propagator.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::pwgem::photonmeson;

enum v0TypesEnum : unsigned int {
  ITSTPC_ITSTPC = 0,
  TPConly_TPConly,
  ITSonly_ITSonly,
  ITSTPC_TPConly,
  ITSTPC_ITSonly,
  TPConly_ITSonly,
};

struct CheckMCV0Tree {
  // Track
  float trkPosPt;
  float trkElePt;
  float trkPosEta;
  float trkEleEta;
  float trkPosTgl;
  float trkEleTgl;
  float trkPosX;
  float trkPosY;
  float trkPosZ;
  float trkEleX;
  float trkEleY;
  float trkEleZ;
  // MC
  float trkPosMCPt;
  float trkEleMCPt;
  float trkPosMCEta;
  float trkEleMCEta;
  float trkPosMCX;
  float trkPosMCY;
  float trkPosMCZ;
  float trkEleMCX;
  float trkEleMCY;
  float trkEleMCZ;
  // MC Prop
  float trkPosMCXProp;
  float trkPosMCYProp;
  float trkPosMCZProp;
  float trkEleMCXProp;
  float trkEleMCYProp;
  float trkEleMCZProp;
  // Type
  unsigned int type;
};

struct CheckMCV0 {
  // Output Objects
  HistogramRegistry registry{"output", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TTree> tree{"checkMC"};
  CheckMCV0Tree treeData;

  // Track selection
  Configurable<bool> writeTree{"writeTree", false, "write tree output with MC and track info"};
  Configurable<bool> ptLogAxis{"ptLogAxis", false, "Flag to use a log momentum axis"};
  Configurable<float> minpt{"minpt", 0.001, "min pt for track in GeV/c"};
  Configurable<float> maxpt{"maxpt", 20.0, "max pt for track in GeV/c"};
  Configurable<float> maxeta{"maxeta", 999.0, "eta acceptance"};
  Configurable<float> dcamin{"dcamin", 0.1, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};
  Configurable<float> maxZ{"maxZ", 200.0, "max z for track"};
  Configurable<float> maxX{"maxX", 200.0, "maximum X (starting point of track X)"};
  Configurable<float> maxY{"maxY", 200.0, "maximum Y (starting point of track Y)"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 100.0, "max chi2/NclsTPC"};
  Configurable<int> minCrossedRowsTPC{"minCrossedRowsTPC", 10, "min crossed rows tpc"};

  // Filters
  Filter trackPt = aod::track::pt > minpt&& aod::track::pt < maxpt;
  Filter trackZ = nabs(aod::track::z) < maxZ;
  Filter trackX = aod::track::x < maxX;
  Filter trackY = nabs(aod::track::y) < maxY;
  Filter trackEta = nabs(aod::track::eta) < maxeta;
  Filter trackDCA = nabs(aod::track::dcaXY) > dcamin&& nabs(aod::track::dcaXY) < dcamax;
  using TracksMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TracksCovIU, aod::McTrackLabels>;
  using FilteredTracksMC = soa::Filtered<TracksMC>;

  // Histogram Parameters
  Configurable<int> tglNBins{"tglNBins", 500, "nBins for tgl"};
  Configurable<int> zNBins{"zNBins", 2 * static_cast<int>(maxZ), "nBins for z"};
  Configurable<int> xNBins{"xNBins", static_cast<int>(maxX), "nBins for x"};
  Configurable<int> yNBins{"yNBins", static_cast<int>(maxY), "nBins for y"};
  Configurable<int> ptNBins{"ptNBins", 200, "nBins for pt"};
  Configurable<float> ptLowCut{"ptLowCut", 0.1, "low pt cut"};

  static constexpr std::array<std::string_view, 7> cutsBinLabels{"Pre", "checkV0leg", "sign", "propagationFailed", "lowPt", "highPt", "survived"};
  enum cutsBinEnum {
    PRE = 1,
    CHECKV0LEG,
    SIGN,
    PROPFAIL,
    LOWPT,
    HIGHPT,
    SURVIVED,
  };
  static_assert(cutsBinLabels.size() == cutsBinEnum::SURVIVED);

  static constexpr std::array<std::string_view, 8> checkV0legLabels{"Pt<minPt", "Pt>maxPt", "dca<dcamin", "dca>dcamax", "tpcChi2NCl>maxchi2tpc", "no ITS||TPC", "eta>maxEta", "tpcCrossedRows<minCrossedRowsTPC"};
  enum checkV0legEnum {
    PTMIN = 1,
    PTMAX,
    DCAMIN,
    DCAMAX,
    TPCCHI2NCL,
    NOITSTPC,
    MAXETA,
    MINCROSSEDROWSTPC,
  };
  static_assert(checkV0legLabels.size() == checkV0legEnum::MINCROSSEDROWSTPC);

  // CCDB
  Configurable<std::string> ccdbPath{"ccdb-path", "GLO/GRP/GRP", "path to the ccdb object"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "path to the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int runNumber = -1;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::parameters::GRPMagField* grpmag = nullptr;

  // params
  std::array<float, 3> mcPosXYZProp{};
  std::array<float, 3> mcEleXYZProp{};

  // Track Types
  static constexpr std::array<std::string_view, 6> v0Types = {"ITSTPC_ITSTPC/", "TPConly_TPConly/", "ITSonly_ITSonly/", "ITSTPC_TPConly/", "ITSTPC_ITSonly/", "TPConly_ITSonly/"};
  std::array<bool, v0Types.size()> v0TypesPassed{};
  static constexpr std::array<std::string_view, 3> ptBins = {"lowPt/", "highPt/", "all/"};

  void init(InitContext const& /*unused*/)
  {
    if (writeTree) {
      setTree();
    }

    // CCDB
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));

    const AxisSpec axisTgl{tglNBins, -5.f, +5.f, "tan(#lambda)"};
    const AxisSpec axisTglEle{tglNBins, -5.f, +5.f, "tan(#lambda) of e^{-}"};
    const AxisSpec axisTglPos{tglNBins, -5.f, +5.f, "tan(#lambda) of e^{+}"};
    const AxisSpec axisEta{300, -1.6, +1.6, "#eta"};
    const AxisSpec axisX{xNBins, 0.0, maxX, "reco. x"};
    const AxisSpec axisXMC{xNBins, 0.0, maxX, "MC prop. x"};
    const AxisSpec axisXEle{xNBins, 0.0, maxX, "reco. x of e^{-}"};
    const AxisSpec axisXPos{xNBins, 0.0, maxX, "reco. x of e^{+}"};
    const AxisSpec axisXMCEle{xNBins, 0.0, maxX, "MC prop. x of e^{-}"};
    const AxisSpec axisXMCPos{xNBins, 0.0, maxX, "MC prop. x of e^{+}"};
    const AxisSpec axisY{yNBins, -maxY, maxY, "reco. y"};
    const AxisSpec axisYMC{yNBins, -maxY, maxY, "MC prop. y"};
    const AxisSpec axisYEle{yNBins, -maxY, maxY, "reco. y of e^{-}"};
    const AxisSpec axisYPos{yNBins, -maxY, maxY, "reco. y of e^{+}"};
    const AxisSpec axisYMCEle{yNBins, -maxY, maxY, "MC prop. y of e^{-}"};
    const AxisSpec axisYMCPos{yNBins, -maxY, maxY, "MC prop. y of e^{+}"};
    const AxisSpec axisZ{zNBins, -maxZ, maxZ, "reco. z"};
    const AxisSpec axisZMC{zNBins, -maxZ, maxZ, "MC prop.z"};
    const AxisSpec axisZEle{zNBins, -maxZ, maxZ, "reco. z of e^{-}"};
    const AxisSpec axisZPos{zNBins, -maxZ, maxZ, "reco. z of e^{+}"};
    const AxisSpec axisZMCEle{zNBins, -maxZ, maxZ, "MC prop. z of e^{-}"};
    const AxisSpec axisZMCPos{zNBins, -maxZ, maxZ, "MC prop. z of e^{+}"};
    const AxisSpec axisYDiff{yNBins, -maxY, maxY, "reco. y of e^{+} - reco. y of e^{-}"};
    const AxisSpec axisZDiff{zNBins, -maxZ, maxZ, "reco. z of e^{+} - reco. z of e^{-}"};
    AxisSpec axisPt{ptNBins, minpt, maxpt, "#it{p}_{T} GeV/#it{c}"};
    AxisSpec axisPtPos{ptNBins, minpt, maxpt, "#it{p}_{T} GeV/#it{c} of e^{+}"};
    AxisSpec axisPtEle{ptNBins, minpt, maxpt, "#it{p}_{T} GeV/#it{c} of e^{-}"};
    if (ptLogAxis) {
      axisPt.makeLogarithmic();
      axisPtPos.makeLogarithmic();
      axisPtEle.makeLogarithmic();
    }

    for (const auto& v0type : v0Types) {
      for (const auto& ptbin : ptBins) {
        registry.add(Form("%s%sTglTgl", v0type.data(), ptbin.data()), "tan(#lambda) vs. tan(#lambda)", HistType::kTH2F, {axisTglPos, axisTglEle});
        registry.add(Form("%s%sPtPt", v0type.data(), ptbin.data()), "p_{T} vs. p_{T}", HistType::kTH2F, {axisPtPos, axisPtEle});
        registry.add(Form("%s%sXX", v0type.data(), ptbin.data()), "x vs. x", HistType::kTH2F, {axisXPos, axisXEle});
        registry.add(Form("%s%sYY", v0type.data(), ptbin.data()), "y vs. y", HistType::kTH2F, {axisYPos, axisYEle});
        registry.add(Form("%s%sZZ", v0type.data(), ptbin.data()), "z vs. z", HistType::kTH2F, {axisZPos, axisZEle});
        registry.add(Form("%s%sXZ", v0type.data(), ptbin.data()), "x vs. z", HistType::kTH2F, {axisX, axisZ});
        registry.add(Form("%s%sZPt", v0type.data(), ptbin.data()), "z vs. p_{T}", HistType::kTH2F, {axisZ, axisPt});
        registry.add(Form("%s%sXYDiff", v0type.data(), ptbin.data()), "x vs. y", HistType::kTH2F, {axisXPos, axisYDiff});
        registry.add(Form("%s%sXZDiff", v0type.data(), ptbin.data()), "x vs. Z", HistType::kTH2F, {axisXPos, axisZDiff});
        registry.add(Form("%s%sZTgl", v0type.data(), ptbin.data()), "z vs. tan(#lambda)", HistType::kTH2F, {axisZ, axisTgl});
        registry.add(Form("%s%sPtTgl", v0type.data(), ptbin.data()), "p_{T} vs. tan(#lambda)", HistType::kTH2F, {axisPt, axisTgl});
        registry.add(Form("%s%sEta", v0type.data(), ptbin.data()), "#eta", HistType::kTH1F, {axisEta});
        registry.add(Form("%s%sMC/XX", v0type.data(), ptbin.data()), Form("MC x vs. reconstructed x (ptLowCut=%.2f)", ptLowCut.value), HistType::kTH2F, {axisXMC, axisX});
        registry.add(Form("%s%sMC/YY", v0type.data(), ptbin.data()), Form("MC y vs. rconstructed y (ptLowCut=%.2f)", ptLowCut.value), HistType::kTH2F, {axisYMC, axisY});
        registry.add(Form("%s%sMC/ZZ", v0type.data(), ptbin.data()), Form("MC z vs. reconstructed z (ptLowCut=%.2f)", ptLowCut.value), HistType::kTH2F, {axisZMC, axisZ});
        registry.add(Form("%s%sMC/VertexPropagationX", v0type.data(), ptbin.data()), Form("MC vertex X propagated to track (ptLowCut=%.2f)", ptLowCut.value), HistType::kTH2F, {{axisXMCPos, axisXMCEle}});
        registry.add(Form("%s%sMC/VertexPropagationY", v0type.data(), ptbin.data()), Form("MC vertex Y propagated to track (ptLowCut=%.2f)", ptLowCut.value), HistType::kTH2F, {{axisYMCPos, axisYMCEle}});
        registry.add(Form("%s%sMC/VertexPropagationZ", v0type.data(), ptbin.data()), Form("MC vertex Z propagated to track (ptLowCut=%.2f)", ptLowCut.value), HistType::kTH2F, {{axisZMCPos, axisZMCEle}});
      }
    }

    registry.add("V0Counter", "V0 counter", HistType::kTH1F, {{cutsBinLabels.size(), 0.5, 0.5 + cutsBinLabels.size()}});
    for (int iBin = 0; iBin < cutsBinLabels.size(); ++iBin) {
      registry.get<TH1>(HIST("V0Counter"))->GetXaxis()->SetBinLabel(iBin + 1, cutsBinLabels[iBin].data());
    }

    registry.add("V0TypeCounter", "V0 Type counter", HistType::kTH1F, {{v0Types.size(), 0.5, 0.5 + v0Types.size()}});
    for (int iBin = 0; iBin < v0Types.size(); ++iBin) {
      registry.get<TH1>(HIST("V0TypeCounter"))->GetXaxis()->SetBinLabel(iBin + 1, v0Types[iBin].data());
    }

    registry.add("CheckV0Leg", "CheckV0Leg", HistType::kTH1F, {{checkV0legLabels.size(), 0.5, 0.5 + checkV0legLabels.size()}});
    for (int iBin = 0; iBin < checkV0legLabels.size(); ++iBin) {
      registry.get<TH1>(HIST("CheckV0Leg"))->GetXaxis()->SetBinLabel(iBin + 1, checkV0legLabels[iBin].data());
    }
  }

  void processV0(soa::Join<aod::McCollisionLabels, aod::Collisions> const& collisions, aod::V0s const& v0s, FilteredTracksMC const& /*unused*/, aod::McParticles const& /*unused*/, aod::BCsWithTimestamps const& /*unused*/, aod::McCollisions const& /*unused*/)
  {
    initCCDB(collisions.begin().bc_as<aod::BCsWithTimestamps>());

    for (const auto& v0 : v0s) {
      registry.fill(HIST("V0Counter"), PRE);

      // tracks
      const auto pos = v0.template posTrack_as<TracksMC>(); // positive daughter
      const auto ele = v0.template negTrack_as<TracksMC>(); // negative daughter
      if (!checkV0leg(pos) || !checkV0leg(ele)) {
        registry.fill(HIST("V0Counter"), CHECKV0LEG);
        continue;
      }

      const auto posMC = pos.template mcParticle_as<aod::McParticles>();
      const auto eleMC = ele.template mcParticle_as<aod::McParticles>();
      if (std::signbit(pos.sign()) || !std::signbit(ele.sign())) { // wrong sign and same sign reject
        registry.fill(HIST("V0Counter"), SIGN);
        continue;
      }
      checkPassed(pos, ele);

      // Propagate the vertex position of the MC track to the reconstructed vertex position of the track
      if (!recacluateMCVertex(ele, eleMC, mcEleXYZProp) || !recacluateMCVertex(pos, posMC, mcPosXYZProp)) {
        registry.fill(HIST("V0Counter"), PROPFAIL);
        continue;
      }

      registry.fill(HIST("V0Counter"), SURVIVED);

      if (writeTree) {
        fillTree(pos, posMC, ele, eleMC);
      }

      static_for<0, v0Types.size() - 1>([&](auto i) {
        static_for<0, ptBins.size() - 1>([&](auto j) {
          fillHistograms<i, j>(pos, ele);
        });
      });
    }
  }
  PROCESS_SWITCH(CheckMCV0, processV0, "process reconstructed info", true);

  template <unsigned int v0Type, unsigned int ptBin, typename TTrack>
  void fillHistograms(TTrack const& pos, TTrack const& ele)
  {
    if (!v0TypesPassed[v0Type]) {
      return;
    }
    if constexpr (ptBin == 0) {
      if (pos.pt() > ptLowCut && ele.pt() > ptLowCut) {
        return;
      }
    } else if (ptBin == 1) {
      if (pos.pt() < ptLowCut && ele.pt() < ptLowCut) {
        return;
      }
    }

    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("TglTgl"), pos.tgl(), ele.tgl());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("PtPt"), pos.pt(), ele.pt());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("XX"), pos.x(), ele.x());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("YY"), pos.y(), ele.y());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("ZZ"), pos.z(), ele.z());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("XZ"), pos.x(), pos.z());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("XZ"), ele.x(), ele.z());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("ZPt"), pos.z(), pos.pt());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("ZPt"), ele.z(), ele.pt());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("XYDiff"), pos.x(), pos.y() - ele.y());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("XZDiff"), pos.x(), pos.z() - ele.z());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("ZTgl"), pos.z(), pos.tgl());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("ZTgl"), ele.z(), ele.tgl());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("PtTgl"), pos.pt(), pos.tgl());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("PtTgl"), ele.pt(), ele.tgl());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("Eta"), pos.eta());
    registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("Eta"), ele.eta());
    if constexpr (ptBin == 0) {
      if (pos.pt() < ptLowCut || ele.pt() < ptLowCut) {
        registry.fill(HIST("V0Counter"), cutsBinEnum::LOWPT);
      }
      if (pos.pt() < ptLowCut) {
        registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("XX"), mcPosXYZProp[0], pos.x());
        registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("YY"), mcPosXYZProp[1], pos.y());
        registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("ZZ"), mcPosXYZProp[2], pos.z());
      }
      if (ele.pt() < ptLowCut) {
        registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("XX"), mcEleXYZProp[0], ele.x());
        registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("YY"), mcEleXYZProp[1], ele.y());
        registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("ZZ"), mcEleXYZProp[2], ele.z());
      }
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("VertexPropagationX"), mcPosXYZProp[0], mcEleXYZProp[0]);
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("VertexPropagationY"), mcPosXYZProp[1], mcEleXYZProp[1]);
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("VertexPropagationZ"), mcPosXYZProp[2], mcEleXYZProp[2]);
    } else if (ptBin == 1) {
      if (pos.pt() > ptLowCut || ele.pt() > ptLowCut) {
        registry.fill(HIST("V0Counter"), cutsBinEnum::HIGHPT);
      }
      if (pos.pt() > ptLowCut) {
        registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("XX"), mcPosXYZProp[0], pos.x());
        registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("YY"), mcPosXYZProp[1], pos.y());
        registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("ZZ"), mcPosXYZProp[2], pos.z());
      }
      if (ele.pt() > ptLowCut) {
        registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("XX"), mcEleXYZProp[0], ele.x());
        registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("YY"), mcEleXYZProp[1], ele.y());
        registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("ZZ"), mcEleXYZProp[2], ele.z());
      }
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("VertexPropagationX"), mcPosXYZProp[0], mcEleXYZProp[0]);
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("VertexPropagationY"), mcPosXYZProp[1], mcEleXYZProp[1]);
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("VertexPropagationZ"), mcPosXYZProp[2], mcEleXYZProp[2]);
    } else {
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("XX"), mcPosXYZProp[0], pos.x());
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("YY"), mcPosXYZProp[1], pos.y());
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("ZZ"), mcPosXYZProp[2], pos.z());
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("XX"), mcEleXYZProp[0], ele.x());
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("YY"), mcEleXYZProp[1], ele.y());
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("ZZ"), mcEleXYZProp[2], ele.z());
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("VertexPropagationX"), mcPosXYZProp[0], mcEleXYZProp[0]);
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("VertexPropagationY"), mcPosXYZProp[1], mcEleXYZProp[1]);
      registry.fill(HIST(v0Types[v0Type]) + HIST(ptBins[ptBin]) + HIST("MC/") + HIST("VertexPropagationZ"), mcPosXYZProp[2], mcEleXYZProp[2]);
    }
  }

  void processDummy(aod::V0s const& v0s) {}
  PROCESS_SWITCH(CheckMCV0, processDummy, "process dummy", false);

  // Templates
  template <typename TTrack>
  bool checkV0leg(TTrack const& track)
  {
    if (track.pt() < minpt) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::PTMIN);
      return false;
    }
    if (track.pt() > maxpt) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::PTMAX);
      return false;
    }
    if (!track.hasITS() && !track.hasTPC()) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::NOITSTPC);
      return false;
    }
    if (abs(track.eta()) > maxeta) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::MAXETA);
      return false;
    }
    if (abs(track.dcaXY()) < dcamin) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::DCAMIN);
      return false;
    }
    if (abs(track.dcaXY()) > dcamax) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::DCAMAX);
      return false;
    }
    if (track.tpcChi2NCl() > maxChi2TPC) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::TPCCHI2NCL);
      return false;
    }
    if (track.tpcNClsCrossedRows() < minCrossedRowsTPC) {
      registry.fill(HIST("CheckV0Leg"), checkV0legEnum::MINCROSSEDROWSTPC);
      return false;
    }
    return true;
  }

  template <typename TTrack>
  void checkPassed(TTrack const& track0, TTrack const& track1)
  {
    v0TypesPassed.fill(false); // reset
    v0TypesPassed[0] = isITSTPC_ITSTPC(track0, track1);
    v0TypesPassed[1] = isTPConly_TPConly(track0, track1);
    v0TypesPassed[2] = isITSonly_ITSonly(track0, track1);
    v0TypesPassed[3] = isITSTPC_TPConly(track0, track1);
    v0TypesPassed[4] = isITSTPC_ITSonly(track0, track1);
    v0TypesPassed[5] = isTPConly_ITSonly(track0, track1);
    for (int i = 0; i < v0TypesPassed.size(); ++i) {
      if (v0TypesPassed[i]) {
        registry.fill(HIST("V0TypeCounter"), i + 1);
        treeData.type = v0TypesEnum(i);
      }
    }
  }

  void setTree()
  {

    tree.setObject(new TTree("checkMC", "checkMC"));

    // Reconstructed data
    tree->Branch("trkPosPt", &treeData.trkPosPt, "trkPosPt/F");
    tree->Branch("trkElePt", &treeData.trkElePt, "trkElePt/F");
    tree->Branch("trkPosEta", &treeData.trkPosEta, "trkPosEta/F");
    tree->Branch("trkEleEta", &treeData.trkEleEta, "trkEleEta/F");
    tree->Branch("trkPosTgl", &treeData.trkPosTgl, "trkPosTgl/F");
    tree->Branch("trkEleTgl", &treeData.trkEleTgl, "trkEleTgl/F");
    tree->Branch("trkPosX", &treeData.trkPosX, "trkPosX/F");
    tree->Branch("trkPosY", &treeData.trkPosY, "trkPosY/F");
    tree->Branch("trkPosZ", &treeData.trkPosZ, "trkPosZ/F");
    tree->Branch("trkEleX", &treeData.trkEleX, "trkEleX/F");
    tree->Branch("trkEleY", &treeData.trkEleY, "trkEleY/F");
    tree->Branch("trkEleZ", &treeData.trkEleZ, "trkEleZ/F");
    // MC particle
    tree->Branch("trkPosMCPt", &treeData.trkPosMCPt, "trkPosMCPt/F");
    tree->Branch("trkEleMCPt", &treeData.trkEleMCPt, "trkEleMCPt/F");
    tree->Branch("trkPosMCEta", &treeData.trkPosMCEta, "trkPosMCEta/F");
    tree->Branch("trkEleMCEta", &treeData.trkEleMCEta, "trkEleMCEta/F");
    tree->Branch("trkPosMCX", &treeData.trkPosMCX, "trkPosMCX/F");
    tree->Branch("trkPosMCY", &treeData.trkPosMCY, "trkPosMCY/F");
    tree->Branch("trkPosMCZ", &treeData.trkPosMCZ, "trkPosMCZ/F");
    tree->Branch("trkEleMCX", &treeData.trkEleMCX, "trkEleMCX/F");
    tree->Branch("trkEleMCY", &treeData.trkEleMCY, "trkEleMCY/F");
    tree->Branch("trkEleMCZ", &treeData.trkEleMCZ, "trkEleMCZ/F");
    // Propagated Track position
    tree->Branch("trkPosMCXProp", &treeData.trkPosMCXProp, "trkPosMCXProp/F");
    tree->Branch("trkPosMCYProp", &treeData.trkPosMCYProp, "trkPosMCYProp/F");
    tree->Branch("trkPosMCZProp", &treeData.trkPosMCZProp, "trkPosMCZProp/F");
    tree->Branch("trkEleMCXProp", &treeData.trkEleMCXProp, "trkEleMCXProp/F");
    tree->Branch("trkEleMCYProp", &treeData.trkEleMCYProp, "trkEleMCYProp/F");
    tree->Branch("trkEleMCZProp", &treeData.trkEleMCZProp, "trkEleMCZProp/F");
    // Track Type
    tree->Branch("trkType", &treeData.type, "trkType/i");
  }

  template <typename TTrack, typename MCTrack>
  void fillTree(TTrack const& pos, MCTrack const& posMC, TTrack const& ele, MCTrack const& eleMC)
  {
    // Reconstructed data
    treeData.trkPosPt = pos.pt();
    treeData.trkElePt = ele.pt();
    treeData.trkPosEta = pos.eta();
    treeData.trkEleEta = ele.eta();
    treeData.trkPosTgl = pos.tgl();
    treeData.trkEleTgl = ele.tgl();
    treeData.trkPosX = pos.x();
    treeData.trkPosY = pos.y();
    treeData.trkPosZ = pos.z();
    treeData.trkEleX = ele.x();
    treeData.trkEleY = ele.y();
    treeData.trkEleZ = ele.z();
    // MC particle
    treeData.trkPosMCPt = posMC.pt();
    treeData.trkEleMCPt = eleMC.pt();
    treeData.trkPosMCEta = posMC.eta();
    treeData.trkEleMCEta = eleMC.eta();
    treeData.trkPosMCX = posMC.vx();
    treeData.trkPosMCY = posMC.vy();
    treeData.trkPosMCZ = posMC.vz();
    treeData.trkEleMCX = eleMC.vx();
    treeData.trkEleMCY = eleMC.vy();
    treeData.trkEleMCZ = eleMC.vz();
    // Propagated Track position
    treeData.trkPosMCXProp = mcPosXYZProp[0];
    treeData.trkPosMCYProp = mcPosXYZProp[1];
    treeData.trkPosMCZProp = mcPosXYZProp[2];
    treeData.trkEleMCXProp = mcEleXYZProp[0];
    treeData.trkEleMCYProp = mcEleXYZProp[1];
    treeData.trkEleMCZProp = mcEleXYZProp[2];

    tree->Fill();
  }

  template <typename TTrack, typename MCTrack>
  bool recacluateMCVertex(TTrack const& track, MCTrack const& mcTrack, std::array<float, 3>& xyz)
  {
    std::array<float, 3> xyzMC{mcTrack.vx(), mcTrack.vy(), mcTrack.vz()};
    std::array<float, 3> pxyzMC{mcTrack.px(), mcTrack.py(), mcTrack.pz()};
    auto pPDG = TDatabasePDG::Instance()->GetParticle(mcTrack.pdgCode());
    if (!pPDG) {
      return false;
    }
    o2::track::TrackPar mctrO2(xyzMC, pxyzMC, TMath::Nint(pPDG->Charge() / 3), false);
    if (!mctrO2.rotate(track.alpha()) || !o2::base::Propagator::Instance()->PropagateToXBxByBz(mctrO2, track.x())) {
      return false;
    }
    xyz = {mctrO2.getX(), mctrO2.getY(), mctrO2.getZ()};
    return true;
  }

  template <typename BC>
  void initCCDB(BC const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    o2::base::Propagator::initFieldFromGRP(grpmag);
    o2::base::Propagator::Instance()->setMatLUT(lut);
    runNumber = bc.runNumber();
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<CheckMCV0>(cfgc, TaskName{"check-mc-v0"})};
}
