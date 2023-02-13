#include "MpdTpcHit.h"
#include "MpdTpcInputHit.h"
#include "MpdTpcRunner.h"
#include "PlotUtils.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TMultiGraph.h>

#include <limits>

void buildHistograms(const Mpd::Tpc::Statistics &statistics,
                     Int_t nTracks,
                     Int_t eventCounter) {
  std::string fileName = "stat_" + std::to_string(eventCounter) + ".root";
  TFile rootFile(fileName.c_str(), "RECREATE");
  std::string histoName = "The number of real tracks: " +
      std::to_string(nTracks);
  TH1I h1("statistics", histoName.c_str(), 10, 0.0, 100.0);

  for (const auto &[key, val]: statistics) {
    h1.Fill(val.quality);
  }
  h1.Write();
  rootFile.Close();
}

size_t realTrackLen(const Mpd::Tpc::InputHitContainer &hits,
                    size_t realTrackId) {
  size_t length = 0;
  for (const auto hit : hits) {
    if (hit.trackId == realTrackId) {
      length++;
    }
  }
  return length;
}

void getRecoTrackParams(const Mpd::Tpc::InputHitContainer &hits,
                        const ActsExamples::ProtoTrack &track,
                        Int_t &recoLen,
                        Int_t &realLen,
                        Int_t &realTrackId) {
  std::unordered_map<Int_t, size_t> counts;
  counts.reserve(track.size());

  for (auto iHit : track) {
    counts[hits.at(iHit).trackId]++;
  }
  size_t maxCount = 0;
  realTrackId = -1;
  for (const auto &[id, count] : counts) {
    if (count > maxCount) {
      maxCount = count;
      realTrackId = id;
    }
  }
  recoLen = maxCount;
  realLen = realTrackLen(hits, realTrackId);
}

// Track parameters (needed for .png graphs)
struct TrackParam {
  // The transversal momentum of the first hit in the track
  Double_t p;

  // For reconstructed track only:
  // the length of the maximum part of the reconstructed track,
  // that corresponds to some real track
  Int_t recoLen;

  // For recontructed track:
  // the length of real track that corresponds to reconstructed track
  // for real track: its length
  Int_t realLen;

  // For recontructed track only:
  // trackId of real track that corresponds to reconstructed track
  Int_t realTrackId;

  Int_t eventNumber;
};

using RealTracksMap = std::map<Int_t, std::pair<TrackParam, bool>>;

void drawGraph(
    Int_t nIntervals,
    // Real tracks (input)
    const RealTracksMap &realTracksMap,

    // Reconstructed tracks (output)
    const std::vector<TrackParam> &recoTrackParams,

    std::string fName,
    // If eventNumber == -1 then draw png for all events
    Int_t eventNumber,

    Int_t startRecoIndex,
    Int_t kTrackId,
    Double_t rightBound) {

  Double_t pMin = std::numeric_limits<Double_t>::max();
  Double_t pMax = 0;

  for (Int_t trackI = startRecoIndex;
      trackI < recoTrackParams.size(); trackI++) {

    Double_t p = recoTrackParams.at(trackI).p;
    if (p < pMin) {pMin = p;}
    if (p > pMax) {pMax = p;}
  }
  Int_t nRealTracks = 0;
  Int_t nReconstructedTracks = 0;
  for (auto & [trackIdExtended, pair] : realTracksMap) {
    Int_t eventNumberC = trackIdExtended / kTrackId;

    if ( (eventNumber >= 0) && (eventNumber!= eventNumberC) ) {
      continue;
    }
    Double_t p = pair.first.p;
    if (p < pMin) {pMin = p;}
    if (p > pMax) {pMax = p;}
    nRealTracks++;
    if (pair.second) {
      nReconstructedTracks++;
    }
  }
  if ((pMin == 0) && (pMax == 0)) {
    std::string ev;
    if (eventNumber == -1) {
      ev = std::string("all events");
    } else {
      ev = std::string("event ") + std::to_string(eventNumber);
    }
    std::cout <<  std::endl << "drawGraph(): " <<
      "; can not draw .png graph for " << ev <<
      "; pMin == 0 and pMax == 0" << std::endl;
    return;
  }
  // Patch boundaries
  if (pMax == pMin) {
    pMax += 1;
  }
  Int_t sumRealAmongReconstructed[nIntervals];
  Int_t sumRealAmongReal         [nIntervals];
  Int_t sumsReconstructed        [nIntervals];
  std::fill_n(sumRealAmongReconstructed, nIntervals, 0);
  std::fill_n(sumRealAmongReal,          nIntervals, 0);
  std::fill_n(sumsReconstructed,         nIntervals, 0);

  Double_t appearanceAll  [nIntervals + 1];
  Double_t appearanceReco [nIntervals + 1];
  std::fill_n(appearanceAll,  nIntervals + 1, 0);
  std::fill_n(appearanceReco, nIntervals + 1, 0);

  Double_t intervalLen = (pMax - pMin) / nIntervals;

  // Loop over real tracks
  for (auto & [trackIdExtended, pair] : realTracksMap) {
    Int_t eventNumberC = trackIdExtended / kTrackId;
    if ( (eventNumber >= 0) && (eventNumber!= eventNumberC) ) {
      continue;
    }
    Double_t p = pair.first.p;
    Int_t ind = static_cast<Int_t>((p - pMin) / intervalLen);
    if (ind == nIntervals) {
      ind--;
    }
    appearanceAll[ind]--;

    Int_t realTrackId = trackIdExtended % kTrackId;
    if (pair.second) { // track was reconstructed
      appearanceReco[ind]--;

      Int_t maxRecoLen =  0;
      Int_t trackIndex = -1;
      for (Int_t trackI = startRecoIndex; trackI < recoTrackParams.size();
           trackI++) {

        if ((recoTrackParams.at(trackI).realTrackId != realTrackId) ||
            (recoTrackParams.at(trackI).eventNumber != eventNumberC)) {
          continue;
        }
        Int_t curRecoLen = recoTrackParams.at(trackI).recoLen;
        if (curRecoLen > maxRecoLen) {
          maxRecoLen = curRecoLen;
          trackIndex = trackI;
        }
      }
      if (trackIndex == -1) {
        std::cout << "drawGraph() ERROR: " <<
            "can't find a reconstructed track for the real track; " <<
            "real trackId = " << realTrackId << std::endl;
        continue;
      }
      // Consistency check
      if (pair.first.realLen != recoTrackParams.at(trackIndex).realLen) {
        std::cout << "drawGraph() WARNING: " <<
            "lentgh(real track) != length(real track " <<
            "corresponding to reco track); " <<
            "real trackId = " << realTrackId <<
            "; len(real track) = " << pair.first.realLen <<
            "; length(real track corresponding to reco track) = " <<
                recoTrackParams.at(trackIndex).realLen << std::endl;
      }
      sumRealAmongReal[ind]          += pair.first.realLen;
      sumRealAmongReconstructed[ind] += pair.first.realLen;
      sumsReconstructed[ind]         += recoTrackParams.at(trackIndex).recoLen;

      // Consistency check
      if (pair.first.realLen < recoTrackParams.at(trackIndex).recoLen) {
        std::cout << "drawGraph() WARNING: for reconstructed track " <<
            trackIndex << " real len = " << pair.first.realLen <<
            " < reconstructed len = "  <<
            recoTrackParams.at(trackIndex).recoLen << std::endl;
      }
    } else { // track was not reconstructed
      sumRealAmongReal[ind]          += pair.first.realLen;
    }
  }
  // Normalization of the number of tracks
  for (Int_t iInterval = 0; iInterval < nIntervals; iInterval++) {
    appearanceAll [iInterval] = 100 * appearanceAll [iInterval] / nRealTracks;
    appearanceReco[iInterval] = 100 * appearanceReco[iInterval] / nRealTracks;
  }

  Double_t pI[nIntervals + 1];

  for (Int_t i = 0; i < nIntervals; i++) {
    pI[i] = pMin + i*intervalLen;
  }
  pI[nIntervals] = rightBound;

  Double_t qualitiesAmongReco[nIntervals + 1];
  Double_t qualities         [nIntervals + 1];

  qualitiesAmongReco[nIntervals] = 0;
  qualities[nIntervals] = 0;

  for (Int_t i = 0; i < nIntervals; i++) {
    if (sumRealAmongReconstructed[i] != 0) {
      qualitiesAmongReco[i] = 100.*sumsReconstructed[i] /
                                   sumRealAmongReconstructed[i];
    } else {
      qualitiesAmongReco[i] = 0;
    }
    if (sumRealAmongReal[i] != 0) {
      qualities[i] = 100.*sumsReconstructed[i] / sumRealAmongReal[i];
    } else {
      qualities[i] = 0;
    }
  }
  Int_t canvasY = 3000;
  Int_t canvasX = Int_t(canvasY * 1.161); // golden ratio for beauty
  TCanvas canvas("Quality_on_momentum", "Quality on momentum",
                 canvasX, canvasY);

  TMultiGraph mg;

  TGraph qualitiesAmongRecoGraph(nIntervals + 1, pI, qualitiesAmongReco);
  qualitiesAmongRecoGraph.SetFillColor(38);
  mg.Add(&qualitiesAmongRecoGraph, "B");

  TGraph qualitiesGraph(nIntervals + 1, pI, qualities);
  qualitiesGraph.SetFillColor(kBlue + 2);
  mg.Add(&qualitiesGraph, "B");

  TGraph appearanceAllGraph(nIntervals + 1, pI, appearanceAll);
  appearanceAllGraph.SetFillColor(kRed - 9);
  mg.Add(&appearanceAllGraph, "B");

  TGraph appearanceRecoGraph(nIntervals + 1, pI, appearanceReco);
  appearanceRecoGraph.SetFillColor(kRed);
  mg.Add(&appearanceRecoGraph, "B");

  std::string title = "The number of reconstructed tracks: " +
                      std::to_string(nReconstructedTracks) + " / " +
                      std::to_string(nRealTracks);
  mg.SetTitle(title.c_str());
  mg.GetXaxis()->SetTitle("Pt (GeV)");
  mg.GetYaxis()->SetTitle("Quality");
  mg.Draw("A");
  canvas.Print(fName.c_str());
}

// Draw png graphics quality depending on momentum
//   for every event and
//   for all events
void plotQualityOnP(const Mpd::Tpc::InputHitContainer &hits,
                    const Mpd::Tpc::ProtoTrackContainer &trajectories,
                    Int_t eventCounter) {
  const Double_t rightBound = 2.5; // GeV
  const Int_t kTrackId = 1000000;
  if (eventCounter == 0) {
    return;
  }
  // Map for all real tracks:
  // trackId -> std::pair(trackParam         : TrackParam
  //                      is track was reconstructed: bool)
  static RealTracksMap realTracksMap;

  // Collecting info about real tracks
  Int_t nRealTracks = 0;
  for (auto hit : hits) {
    Int_t trackId = hit.trackId;
    Int_t trackIdExtended = eventCounter * kTrackId + trackId;

    if (realTracksMap.find(trackIdExtended) == realTracksMap.end()) {
      nRealTracks++;
      TrackParam trackParam;
      trackParam.p = std::hypot(hit.momentum[0], hit.momentum[1]);
      trackParam.recoLen = 0;
      trackParam.realLen = 0;
      trackParam.eventNumber = eventCounter;

      bool isReconstructed = false;
      realTracksMap[trackIdExtended] = std::make_pair(trackParam, isReconstructed);
    }
    realTracksMap[trackIdExtended].first.realLen++;
  }
  static std::vector<TrackParam> recoTrackParams;
  Int_t recoLastIndex = recoTrackParams.size();

// Collecting info about reconstructed tracks
  Int_t tI = -1;
  for (auto recoTrack : trajectories) {
    tI++;
    auto hit0Index = recoTrack.at(0);
    auto pX = hits.at(hit0Index).momentum[0];
    auto pY = hits.at(hit0Index).momentum[1];
    Double_t p = std::hypot(pX, pY);

    struct TrackParam trackParam;
    trackParam.p = p;
    Int_t recoLen;
    Int_t realLen;
    Int_t realTrackId;
    getRecoTrackParams(hits, recoTrack, recoLen, realLen, realTrackId);

    bool isReconstructed = true;
    realTracksMap[eventCounter * kTrackId + realTrackId].second = isReconstructed;

    trackParam.recoLen = recoLen;
    trackParam.realLen = realLen;
    trackParam.realTrackId = realTrackId;
    trackParam.eventNumber = eventCounter;
    recoTrackParams.push_back(trackParam);

    // Consistency check
    if (trackParam.realLen !=
        realTracksMap[eventCounter * kTrackId + realTrackId].first.realLen) {
      std::cout << "plotQualityOnP() WARNING: " <<
          "trackParam.realLen != realTracksMap[eventCounter * kTrackId + " <<
          "realTrackId].first.realLen" <<
          "reco track index = " << tI <<
          "; realTrackId = " << realTrackId <<
          "; realLen = "     << trackParam.realLen <<
          "; realTracksMap[eventCounter * kTrackId + " <<
          "realTrackId].first.realLen = " <<
          realTracksMap[eventCounter * kTrackId +
          realTrackId].first.realLen << std::endl;
    }
    if (trackParam.recoLen > trackParam.realLen) {
      std::cout << "plotQualityOnP() WARNING: " <<
          "trackParam.recoLen > trackParam.realLen" <<
          "reco track index = " << tI <<
          "; realTrackId = " << realTrackId <<
          "; recoLen = " << recoLen <<
          "; realLen = " << trackParam.realLen << std::endl;
    }
  }

  const Int_t nIntervals = 50;
  drawGraph(nIntervals,
            realTracksMap,
            recoTrackParams,
            "quality_momentum_ev_" + std::to_string(eventCounter) + ".png",
            eventCounter,
            recoLastIndex,
            kTrackId,
            rightBound);

  drawGraph(nIntervals,
            realTracksMap,
            recoTrackParams,
            std::string("quality_momentum_ev_all.png"),
            -1,
            0,
            kTrackId,
            rightBound);
}

struct PointP {
  Double_t x, y, z;

  PointP(Double_t x, Double_t y, Double_t z):
      x(x), y(y), z(z) {}

  Int_t transform(const CoordSystem coordSystem) {
    Int_t res = 0;
    Double_t xTmp = 0;
    Double_t yTmp = 0;
    Double_t zTmp = 0;
    switch (coordSystem) {
      case XY:
        xTmp = x;
        yTmp = y;
        break;
      case ZY:
        xTmp = z;
        yTmp = y;
        break;
      case ZR:
        xTmp = z;
        yTmp = std::hypot(x, y);
        break;
      default:
        std::cout << "PointP::transform() ERROR: " <<
            "unknown coordinate system type!" << std::endl;
        xTmp = 0;
        yTmp = 0;
        res = -1;
        break;
    }
    x = xTmp;
    y = yTmp;
    z = zTmp;
    return res;
  }
};

void plotOutputTracks(
    Int_t canvasX,
    Int_t canvasY,
    const std::shared_ptr<const Acts::TrackingGeometry> &geometry,
    const Mpd::Tpc::SpacePointContainer &spacePoints,
    const Mpd::Tpc::InputHitContainer &hits,
    const Mpd::Tpc::ProtoTrackContainer &trajectories,
    Int_t eventCounter,
    bool multicoloured,
    Int_t lineWidth,
    CoordSystem coordSystem) {

  TCanvas canvas("outputTrajectories", "Output trajectories", canvasX, canvasY);
  TMultiGraph multiGraph;

  // Detector's sensitive surfaces graph
  std::vector<TGraph*> surfGraphs;
  Acts::GeometryContext actsContext;
  if (coordSystem == CoordSystem::XY) {
    geometry->visitSurfaces([&actsContext, &surfGraphs, &multiGraph](
        const Acts::Surface *surface) {
      std::vector<Double_t> boundsValues = surface->bounds().values();

      Double_t xLeftDown = boundsValues.at(0);
      Double_t yLeftDown = boundsValues.at(1);
      Double_t xRightUp  = boundsValues.at(2);
      Double_t yRightUp  = boundsValues.at(3);

      if (yLeftDown + yRightUp != 0) {
        std::cout << "boundsValues.at(1) + boundsValues.at(3) != 0" <<
            std::endl;
      }
      if (xLeftDown + xRightUp != 0) {
        std::cout << "boundsValues.at(0) + boundsValues.at(2) != 0" <<
            std::endl;
      }
      Acts::Vector2 locBound0(xLeftDown, yLeftDown);
      Acts::Vector2 locBound1(xRightUp,  yRightUp);
      Acts::Vector3 unused;

      Acts::Vector3 globalBound0 = surface->localToGlobal(
          actsContext, locBound0, unused);
      Acts::Vector3 globalBound1 = surface->localToGlobal(
          actsContext, locBound1, unused);

      TGraph *surfGraph = new TGraph;
      surfGraphs.push_back(surfGraph);

      surfGraph->SetPoint(0, globalBound0.x(), globalBound0.y());
      surfGraph->SetPoint(1, globalBound1.x(), globalBound1.y());

      surfGraph->SetMarkerStyle(kFullDotMedium);
      surfGraph->SetLineColor(kYellow);
      surfGraph->SetMarkerColor(kYellow);
      multiGraph.Add(surfGraph, "PL");
    });
  }

  // Space points graph
  TGraph spacePointsGraph;
  spacePointsGraph.SetMarkerStyle(kFullDotMedium);
  spacePointsGraph.SetMarkerColor(kSpring);

  size_t pIndex = 0;
  for (auto spacePoint : spacePoints) {
    PointP point(spacePoint.x(), spacePoint.y(), spacePoint.z());
    point.transform(coordSystem);
    Double_t x = point.x;
    Double_t y = point.y;

    spacePointsGraph.SetPoint(pIndex++, x, y);
  }
  multiGraph.Add(&spacePointsGraph, "P");

  // Input hits graph
  TGraph inputHitsGraph;
  inputHitsGraph.SetMarkerStyle(kFullDotMedium);

  pIndex = 0;
  for (auto hit : hits) {
    PointP point(hit.position[0], hit.position[1], hit.position[2]);
    point.transform(coordSystem);
    Double_t x = point.x;
    Double_t y = point.y;

    inputHitsGraph.SetPoint(pIndex++, x, y);
  }
  multiGraph.Add(&inputHitsGraph, "P");

  // Reconstructed tracks graph
  TGraph outTrajectoryGraphs[trajectories.size()];
  Int_t trackIndex = -1;
  std::vector<Int_t> colors = {kRed, kCyan, kGreen - 3, kBlue - 4,
      kMagenta, kOrange + 7, kOrange - 7};

  for (ActsExamples::ProtoTrack reconstructedTrack : trajectories) {
    trackIndex++;

    TGraph &outTrajectoryGraph = outTrajectoryGraphs[trackIndex];
    outTrajectoryGraph.SetMarkerStyle(kFullDotMedium);
    outTrajectoryGraph.SetLineWidth(lineWidth);
    Int_t color = kRed;
    if (multicoloured) {
      Int_t iColor = trackIndex % colors.size();
      color = colors[iColor];
    }
    outTrajectoryGraph.SetLineColor(color);

    pIndex = 0;
    for (uint32_t hitIndex : reconstructedTrack) {
      PointP point(hits.at(hitIndex).position[0],
                   hits.at(hitIndex).position[1],
                   hits.at(hitIndex).position[2]);
      point.transform(coordSystem);
      Double_t x = point.x;
      Double_t y = point.y;

      outTrajectoryGraph.SetPoint(pIndex++, x, y);
    }
    multiGraph.Add(&outTrajectoryGraph, "PL");
  }
  multiGraph.Draw("A");

  std::string dotPlotFileName = "event_" + std::to_string(eventCounter);
  if (lineWidth == 0) {
    dotPlotFileName += "_w0";
  }
  dotPlotFileName += ".png";
  canvas.Print(dotPlotFileName.c_str());

  for (auto surfGraph : surfGraphs) {
    delete surfGraph;
  }
}
