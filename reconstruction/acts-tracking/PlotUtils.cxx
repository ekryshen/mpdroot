#include "MpdTpcHit.h"
#include "MpdTpcInputHit.h"
#include "MpdTpcRunner.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TMultiGraph.h"

#include <limits>

void buildHistograms(const Mpd::Tpc::Statistics &statistics,
                     const int nTracks,
                     const int eventCounter) {
  std::string fileName = "stat_" + std::to_string(eventCounter) + ".root";
  TFile rootFile(fileName.c_str(), "RECREATE");
  std::string histoName = "The number of real tracks: " + std::to_string(nTracks);
  TH1I h1("statistics", histoName.c_str(), 10, 0.0, 100.0);

  for (auto const &[key, val]: statistics) {
    h1.Fill(val.quality);
  }
  h1.Write();
  rootFile.Close();
}

size_t realTrackLen(const Mpd::Tpc::InputHitContainer &hits,
                    const size_t realTrackId) {
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
                        int *recoLen,
                        int *realLen,
                        int *realTrackId) {
  std::unordered_map<int, size_t> counts;
  counts.reserve(track.size());

  for (auto iHit : track) {
    counts[hits.at(iHit).trackId]++;
  }
  size_t maxCount = 0;
  *realTrackId = -1;
  for (const auto &[id, count] : counts) {
    if (count > maxCount) {
      maxCount = count;
      *realTrackId = id;
    }
  }
  *recoLen = maxCount;
  *realLen = realTrackLen(hits, *realTrackId);
}

// Track parameters (needed for .png graphs)
struct TrackParam {
  double p;    // the transversal momentum of the first hit in the track

               // for reconstructed track only:
  int recoLen; // the length of the maximum part of the reconstructed track,
               // that corresponds to some real track

               // for recontructed track:
  int realLen; // the length of real track that corresponds to reconstructed track
               // for real track: its length

  int realTrackId; // for recontructed track only:
                   // trackId of real track that corresponds to reconstructed track
};

using RealTracksMap = std::map<int, std::pair<TrackParam, bool> >;

void drawGraph(const int nIntervals,
               const RealTracksMap &realTracksMap,              // real tracks (input)
               const std::vector<TrackParam> &recoTrackParams,  // reconstructed tracks (output)
               const std::string fName,

               // if eventNumber == -1 then draw png for all events
               const int eventNumber,

               const int startRecoIndex,
               const int kTrackId,
               const double rightBound) {
  const bool debugPrint = false;
  double pMin = std::numeric_limits<double>::max();
  double pMax = 0;

  for (int trackI = startRecoIndex; trackI < recoTrackParams.size(); trackI++) {
    double p = recoTrackParams.at(trackI).p;
    if (p < pMin) {pMin = p;}
    if (p > pMax) {pMax = p;}
  }
  int nRealTracks = 0;
  int nReconstructedTracks = 0;
  for (auto & [trackIdExtended, pair] : realTracksMap) {
    int eventNumberC = trackIdExtended / kTrackId;

    if ( (eventNumber >= 0) && (eventNumber!= eventNumberC) ) {
      continue;
    }
    double p = pair.first.p;
    if (p < pMin) {pMin = p;}
    if (p > pMax) {pMax = p;}
    nRealTracks++;
    if (pair.second) {  // true
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
  // patch boundaries
  if (pMax == pMin) {
    pMax += 1;
  }
  int sumRealAmongReconstructed[nIntervals];
  int sumRealAmongReal         [nIntervals];
  int sumsReconstructed        [nIntervals];
  std::fill_n(sumRealAmongReconstructed, nIntervals, 0);
  std::fill_n(sumRealAmongReal,          nIntervals, 0);
  std::fill_n(sumsReconstructed,         nIntervals, 0);

  double appearanceAll  [nIntervals + 1];
  double appearanceReco [nIntervals + 1];
  std::fill_n(appearanceAll,  nIntervals + 1, 0);
  std::fill_n(appearanceReco, nIntervals + 1, 0);

  double intervalLen = (pMax - pMin) / nIntervals;

  // loop over real tracks
  for (auto & [trackIdExtended, pair] : realTracksMap) {
    int eventNumberC = trackIdExtended / kTrackId;
    if ( (eventNumber >= 0) && (eventNumber!= eventNumberC) ) {
      continue;
    }
    double p = pair.first.p;
    int ind = int( (p - pMin) / intervalLen );
    if (ind == nIntervals) {
      ind--;
    }
    appearanceAll[ind]--;

    int realTrackId = trackIdExtended % kTrackId;
    if (pair.second) { // track was reconstructed
      appearanceReco[ind]--;

      int maxRecoLen = 0;
      int trackIndex = -1;
      for (int trackI = startRecoIndex; trackI < recoTrackParams.size(); trackI++) {
        if (recoTrackParams.at(trackI).realTrackId != realTrackId) {
          continue;
        }
        int curRecoLen = recoTrackParams.at(trackI).recoLen;
        if (curRecoLen > maxRecoLen) {
          maxRecoLen = curRecoLen;
          trackIndex = trackI;
        }
      }
      if (trackIndex == -1) {
        std::cout << "drawGraph() ERROR: can't find a reconstructed track for the real track; " <<
                     "real trackId = " << realTrackId << std::endl;
        continue;
      }
      if (pair.first.realLen != recoTrackParams.at(trackIndex).realLen) {
         if (debugPrint) {
           std::cout << "drawGraph() WARNING: " <<
               "lentgh(real track) != length(real track corresponding to reco track); " <<
               "real trackId = " << realTrackId <<
               "; len(real track) = " << pair.first.realLen <<
               "; length(real track corresponding to reco track) = recoTrackParams.at(trackIndex).realLen" <<
               "; reco track index = " << trackIndex << std::endl;
        }
      }
      sumRealAmongReal[ind]          += pair.first.realLen;
      sumRealAmongReconstructed[ind] += pair.first.realLen;
      sumsReconstructed[ind]         += recoTrackParams.at(trackIndex).recoLen;
    } else { // track was not reconstructed
      sumRealAmongReal[ind]          += pair.first.realLen;
    }
  }
  double pI[nIntervals + 1];

  for (int i = 0; i < nIntervals; i++) {
    pI[i] = pMin + i*intervalLen;
  }
  pI[nIntervals] = rightBound;

  double qualitiesAmongReco[nIntervals + 1];
  double qualities         [nIntervals + 1];

  qualitiesAmongReco[nIntervals] = 0;
  qualities[nIntervals] = 0;

  for (int i = 0; i < nIntervals; i++) {
    if (sumRealAmongReconstructed[i] != 0) {
      qualitiesAmongReco[i] = 100.*sumsReconstructed[i] / sumRealAmongReconstructed[i];
    } else {
      qualitiesAmongReco[i] = 0;
    }
    if (sumRealAmongReal[i] != 0) {
      qualities[i] = 100.*sumsReconstructed[i] / sumRealAmongReal[i];
    } else {
      qualities[i] = 0;
    }
  }
  int canvasY = 3000;
  int canvasX = int(canvasY * 1.161); // golden ratio for beauty
  TCanvas canvas("Quality_on_momentum", "Quality on momentum", canvasX, canvasY);

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
void drawQualityOnP(const Mpd::Tpc::InputHitContainer &hits,
                    const Mpd::Tpc::ProtoTrackContainer &trajectories,
                    const int eventCounter) {
  double rightBound = 2.5; // GeV
  int kTrackId = 1000000;
  if (eventCounter == 0) {
    return;
  }
  // map for all real tracks:
  // trackId -> std::pair(trackParam         : TrackParam
  //                      is track was reconstructed: bool)
  static RealTracksMap realTracksMap;

// collecting info about real tracks
  int nRealTracks = 0;
  int iHit = -1;
  for (auto hit : hits) {
    iHit++;
    int trackId = hit.trackId;
    int trackIdExtended = eventCounter * kTrackId + trackId;

    if (realTracksMap.find(trackIdExtended) == realTracksMap.end()) {
      nRealTracks++;
      TrackParam trackParam;
      trackParam.p = std::hypot(hit.momentum[0], hit.momentum[1]);
      trackParam.recoLen = 0;
      trackParam.realLen = 0;

      bool isReconstructed = false;
      realTracksMap[trackIdExtended] = std::make_pair(trackParam, isReconstructed);
    }
    realTracksMap[trackIdExtended].first.realLen++;
  }
  static std::vector<TrackParam> recoTrackParams;
  int recoLastIndex = recoTrackParams.size();

// collecting info about reconstructed tracks
  int tI = -1;
  for (auto recoTrack : trajectories) {
    tI++;
    auto hit0Index = recoTrack.at(0);
    auto pX = hits.at(hit0Index).momentum[0];
    auto pY = hits.at(hit0Index).momentum[1];
    double p = std::hypot(pX, pY);

    struct TrackParam trackParam;
    trackParam.p = p;
    int recoLen;
    int realLen;
    int realTrackId;
    getRecoTrackParams(hits, recoTrack, &recoLen, &realLen, &realTrackId);

    bool isReconstructed = true;
    realTracksMap[eventCounter * kTrackId + realTrackId].second = isReconstructed;

    trackParam.recoLen = recoLen;
    trackParam.realLen = realLen;
    trackParam.realTrackId = realTrackId;
    recoTrackParams.push_back(trackParam);
  }

  const int nIntervals = 50;
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

void plotOutputTracks(const int canvasX,
                      const int canvasY,
                      const std::shared_ptr<const Acts::TrackingGeometry> &geometry,
                      const Mpd::Tpc::SpacePointContainer &spacePoints,
                      const Mpd::Tpc::InputHitContainer &hits,
                      const Mpd::Tpc::ProtoTrackContainer &trajectories,
                      const int eventCounter) {
  TCanvas canvas("outputTrajectories", "Output trajectories", canvasX, canvasY);
  TMultiGraph multiGraph;

// detector's sensitive surfaces graph
  std::vector<TGraph*> surfGraphs;
  Acts::GeometryContext actsContext;
  geometry->visitSurfaces([&actsContext, &surfGraphs, &multiGraph](const Acts::Surface *surface) {
    std::vector<double> boundsValues = surface->bounds().values();

    double xLeftDown = boundsValues.at(0);
    double yLeftDown = boundsValues.at(1);
    double xRightUp  = boundsValues.at(2);
    double yRightUp  = boundsValues.at(3);

    if (yLeftDown + yRightUp != 0) {
      std::cout << "boundsValues.at(1) + boundsValues.at(3) != 0" << std::endl;
    }
    if (xLeftDown + xRightUp != 0) {
      std::cout << "boundsValues.at(0) + boundsValues.at(2) != 0" << std::endl;
    }
    Acts::Vector2 locBound0(xLeftDown, yLeftDown);
    Acts::Vector2 locBound1(xRightUp,  yRightUp);
    Acts::Vector3 unused;

    Acts::Vector3 globalBound0 = surface->localToGlobal(actsContext, locBound0, unused);
    Acts::Vector3 globalBound1 = surface->localToGlobal(actsContext, locBound1, unused);

    TGraph *surfGraph = new TGraph;
    surfGraphs.push_back(surfGraph);

    surfGraph->SetPoint(0, globalBound0.x(), globalBound0.y());
    surfGraph->SetPoint(1, globalBound1.x(), globalBound1.y());

    surfGraph->SetMarkerStyle(kFullDotMedium);
    surfGraph->SetLineColor(kYellow);
    surfGraph->SetMarkerColor(kYellow);
    multiGraph.Add(surfGraph, "PL");
  });

// space points graph
  TGraph spacePointsGraph;
  spacePointsGraph.SetMarkerStyle(kFullDotMedium);
  spacePointsGraph.SetMarkerColor(kSpring);

  size_t pIndex = 0;
  for (auto spacePoint : spacePoints) {
    double x = spacePoint.x();
    double y = spacePoint.y();
    spacePointsGraph.SetPoint(pIndex++, x, y);
  }
  multiGraph.Add(&spacePointsGraph, "P");

// input hits graph
  TGraph inputHitsGraph;
  inputHitsGraph.SetMarkerStyle(kFullDotMedium);

  pIndex = 0;
  for (auto hit : hits) {
    double x = hit.position[0];
    double y = hit.position[1];
    inputHitsGraph.SetPoint(pIndex++, x, y);
  }
  multiGraph.Add(&inputHitsGraph, "P");

// reconstructed tracks graph
  TGraph outTrajectoryGraphs[trajectories.size()];
  int trackIndex = -1;
  for (ActsExamples::ProtoTrack reconstructedTrack : trajectories) {
    trackIndex++;

    TGraph &outTrajectoryGraph = outTrajectoryGraphs[trackIndex];
    outTrajectoryGraph.SetMarkerStyle(kFullDotMedium);
    outTrajectoryGraph.SetLineWidth(3);
    outTrajectoryGraph.SetLineColor(kRed);

    pIndex = 0;
    for (uint32_t hitIndex : reconstructedTrack) {
      double x = hits.at(hitIndex).position[0];
      double y = hits.at(hitIndex).position[1];
      outTrajectoryGraph.SetPoint(pIndex++, x, y);
    }
    multiGraph.Add(&outTrajectoryGraph, "PL");
  }
  multiGraph.Draw("A");
  std::string dotPlotFileName = "event_" + std::to_string(eventCounter) + ".png";
  canvas.Print(dotPlotFileName.c_str());

  for (auto surfGraph : surfGraphs) {
    delete surfGraph;
  }
  std::cout << "File created: " << dotPlotFileName << std::endl;
}
