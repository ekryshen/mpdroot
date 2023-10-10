#include "MpdTpcHit.h"
#include "MpdTpcInputHit.h"
#include "MpdTpcRunner.h"
#include "PlotUtils.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TMultiGraph.h>
#include <TText.h>

#include <limits>

void buildHistograms(const Mpd::Tpc::Statistics &statistics,
                     Int_t nTracks,
                     Int_t eventCounter,
                     std::string outPath) {
  auto fname = outPath + "/stat_" + std::to_string(eventCounter) + ".root";
  TFile rootFile(fname.c_str(), "RECREATE");
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

using RealTracksMap = std::map<Int_t, std::pair<TrackParam, Bool_t>>;

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

  for (size_t trackI = startRecoIndex;
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
      for (size_t trackI = startRecoIndex; trackI < recoTrackParams.size();
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
  for (size_t iInterval = 0; iInterval < nIntervals; iInterval++) {
    appearanceAll [iInterval] = 100 * appearanceAll [iInterval] / nRealTracks;
    appearanceReco[iInterval] = 100 * appearanceReco[iInterval] / nRealTracks;
  }

  Double_t pI[nIntervals + 1];

  for (size_t i = 0; i < nIntervals; i++) {
    pI[i] = pMin + i*intervalLen;
  }
  pI[nIntervals] = rightBound;

  Double_t qualitiesAmongReco[nIntervals + 1];
  Double_t qualities         [nIntervals + 1];

  qualitiesAmongReco[nIntervals] = 0;
  qualities[nIntervals] = 0;

  for (size_t i = 0; i < nIntervals; i++) {
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
                    Int_t eventCounter,
                    std::string outPath) {
  const Double_t rightBound = 2.5; // GeV
  const Int_t kTrackId = 1000000;
  if (eventCounter == 0) {
    return;
  }
  // Map for all real tracks:
  // trackId -> std::pair(trackParam         : TrackParam
  //                      is track was reconstructed: Bool_t)
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

      Bool_t isReconstructed = false;
      realTracksMap[trackIdExtended] = std::make_pair(trackParam,
                                                      isReconstructed);
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

    Bool_t isReconstructed = true;
    realTracksMap[eventCounter * kTrackId + realTrackId].second =
        isReconstructed;

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

  Int_t nIntervals = 50;
  auto fnamePrefix = outPath + "/quality_momentum_ev_";
  drawGraph(nIntervals,
            realTracksMap,
            recoTrackParams,
            fnamePrefix + std::to_string(eventCounter) + ".png",
            eventCounter,
            recoLastIndex,
            kTrackId,
            rightBound);

  drawGraph(nIntervals,
            realTracksMap,
            recoTrackParams,
            fnamePrefix + "all.png",
            -1,
            0,
            kTrackId,
            rightBound);
}

struct PointP {
  Double_t x, y, z;

  PointP(Double_t x, Double_t y, Double_t z):
      x(x), y(y), z(z) {}

  Int_t makeProjection(const Projection &projection) {
    Int_t res = 0;
    Double_t xTmp = 0;
    Double_t yTmp = 0;
    Double_t zTmp = 0;
    switch (projection) {
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
        std::cout << "PointP::makeProjection() ERROR: " <<
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

void plotRealTracks(
    Int_t canvasX,
    Int_t canvasY,
    const std::shared_ptr<const Acts::TrackingGeometry> &geometry,
    const Mpd::Tpc::InputHitContainer &hits,
    Int_t eventCounter,
    std::string outPath,
    Bool_t multicoloured,
    Int_t lineWidth,
    Projection projection,
    Double_t zmin,
    Double_t zmax,
    Double_t rmax,
    Bool_t plotLines,
    std::string namePostfix,
    Bool_t grid,
    Bool_t realAspectRatio,
    Int_t txtSize,
    Int_t txtStep) {

  TCanvas canvas("inputTrajectories", "Input trajectories", canvasX, canvasY);
  TMultiGraph multiGraph;

  if (grid) {
    canvas.SetGrid();
  }

  // Detector's sensitive surfaces graph
  std::vector<TGraph*> surfGraphs;
  Acts::GeometryContext actsContext;
  if (projection == Projection::XY) {
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

  std::map<Int_t, std::vector<Int_t>> hitsOnTracks;
  size_t nHits = hits.size();
  size_t nTracks = 0;
  for (size_t iHit = 0; iHit < nHits; iHit++) {
    auto &hit = hits.at(iHit);
    auto trackId = hit.trackId;
    auto &hitsVector = hitsOnTracks[trackId];
    hitsVector.push_back(iHit);
    nTracks++;
  }

  for (auto &[trackId, hitsVector] : hitsOnTracks) {
    std::sort(hitsVector.begin(), hitsVector.end());
  }

  // Lines graphs
  std::vector<TGraph*> linesGraphs;
  if (plotLines) {
    for (Int_t i = 0; i < nHits; i++) {
      auto &hit0 = hits.at(i);
      PointP point0(
          hit0.position[0],
          hit0.position[1],
          hit0.position[2]);
      point0.makeProjection(projection);
      Double_t x0 = point0.x;
      Double_t y0 = point0.y;

      for (Int_t j = 0; j < i; j++) {
        auto &hit1 = hits.at(j);
        PointP point1(
            hit1.position[0],
            hit1.position[1],
            hit1.position[2]);
        point1.makeProjection(projection);
        Double_t x1 = point1.x;
        Double_t y1 = point1.y;

        // pass gorizontal or vertical lines
        if ((x0 == x1) || (y0 == y1)) {
          continue;
        }

        Double_t p0x;
        Double_t p0y;

        if ((fabs(x1) >= fabs(x0)) && (fabs(y1) >= fabs(y0))) {
          p0x = x1;
          p0y = y1;
        } else if ((fabs(x1) <= fabs(x0)) && (fabs(y1) <= fabs(y0))) {
          p0x = x0;
          p0y = y0;
        } else {
          continue;
        }

        Double_t p1x;
        Double_t p1y;

        Double_t k = (y1 - y0) / (x1 - x0);
        Double_t b = (y0 * x1 - y1 * x0) / (x1 - x0);

        if (p0y >= 0) {
          if (b >= 0) {
            p1x = -b/k;
            p1y = 0;
          } else {
            p1x = 0;
            p1y = b;
          }
        } else if (p0y < 0) {
          if (b >= 0) {
            p1x = 0;
            p1y = b;
          } else {
            p1x = -b/k;
            p1y = 0;
          }
        }
        auto *lineGraph = new TGraph;
        lineGraph->SetPoint(0, p0x, p0y);
        lineGraph->SetPoint(1, p1x, p1y);

        linesGraphs.push_back(lineGraph);

        multiGraph.Add(lineGraph, "AC");
      }
    }
  }

  // Text labels
  std::vector<TText*> labels;
  auto txtOffset  = 1. * canvasX / 3000. * 2;
  auto txtColor   = kBlack;
  auto txtFont    = 43;

  size_t pIndex = 0;

  // Real tracks graph and text labels
  std::vector<TGraph*> inputTrajectoriesGraph(nTracks);

  std::vector<Int_t> colors = {kRed, kCyan, kGreen - 3, kBlue - 4,
      kMagenta, kOrange + 7, kOrange - 7};

  for (const auto &[trackId, hitsVector] : hitsOnTracks) {
    auto *trajectoryGraph = new TGraph;
    inputTrajectoriesGraph.push_back(trajectoryGraph);

    trajectoryGraph->SetMarkerStyle(kFullDotMedium);
    trajectoryGraph->SetLineWidth(lineWidth);
    Int_t color = kRed;
    if (multicoloured) {
      Int_t iColor = trackId % colors.size();
      color = colors[iColor];
    }
    trajectoryGraph->SetLineColor(color);

    pIndex = 0;
    for (size_t hitIdx : hitsVector) {
      auto &hit = hits.at(hitIdx);
      PointP point(
          hit.position[0],
          hit.position[1],
          hit.position[2]);
      point.makeProjection(projection);
      Double_t x = point.x;
      Double_t y = point.y;
      trajectoryGraph->SetPoint(
          pIndex++, x, y);

      // Text labels
      auto plotTxt = false;
      if (pIndex % txtStep == 0) plotTxt = true;
      Int_t m = hitsVector.size() / 2;
      if ((hitsVector.size() < txtStep) && (
          (pIndex == 1) || (pIndex == hitsVector.size()) ||
          (pIndex == m + 1))) {
        plotTxt = true;
      }

      if (plotTxt) {
        auto txt = new TText(x, y, std::to_string(trackId).c_str());
        txt->SetTextAlign(12);
        txt->SetTextColor(txtColor);
        txt->SetTextFont(txtFont);
        txt->SetTextSizePixels(txtSize);
        labels.push_back(txt);
      }
    }
    if (pIndex) {
      multiGraph.Add(trajectoryGraph, "PL");
    }
  }

  // Input hits graph
  TGraph hitsGraph;
  hitsGraph.SetMarkerStyle(kFullDotMedium);

  pIndex = 0;
  for (auto hit : hits) {
    PointP point(
        hit.position[0],
        hit.position[1],
        hit.position[2]);
    point.makeProjection(projection);
    hitsGraph.SetPoint(
        pIndex++,
        point.x,
        point.y);
  }
  if (pIndex) {
    multiGraph.Add(&hitsGraph, "P");
  }

  if ((projection == Projection::ZY) ||
      (projection == Projection::ZR)) {
    multiGraph.GetXaxis()->SetLimits(zmin, zmax);
    multiGraph.SetMinimum(-rmax);
    multiGraph.SetMaximum( rmax);
  }
  multiGraph.Draw("A");

  for (const auto txt : labels) {
    txt->Draw();
  }

  auto fname = outPath + "/event_" + std::to_string(eventCounter) +
      namePostfix + ".png";
  if (realAspectRatio) {
    canvas.SetRealAspectRatio();
  }
  canvas.Print(fname.c_str());

  for (auto lineGraph: linesGraphs) {
    delete lineGraph;
  }
  for (auto item : inputTrajectoriesGraph) {
    delete item;
  }
  for (auto item : surfGraphs) {
    delete item;
  }
  for (auto item : labels) {
    delete item;
  }
}

void plotOutputTracks(
    Int_t canvasX,
    Int_t canvasY,
    const std::shared_ptr<const Acts::TrackingGeometry> &geometry,
    const Mpd::Tpc::SpacePointContainer &spacePoints,
    const Mpd::Tpc::InputHitContainer &hits,
    const Mpd::Tpc::ProtoTrackContainer &trajectories,
    Int_t eventCounter,
    std::string outPath,
    Bool_t multicoloured,
    Int_t lineWidth,
    Projection projection,
    std::string namePostfix,
    Bool_t plotLabels,
    Int_t txtSize,
    Int_t txtStep) {

  TCanvas canvas("outputTrajectories", "Output trajectories", canvasX, canvasY);
  TMultiGraph multiGraph;

  // Detector's sensitive surfaces graph
  std::vector<TGraph*> surfGraphs;
  Acts::GeometryContext actsContext;
  if (projection == Projection::XY) {
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
    point.makeProjection(projection);
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
    point.makeProjection(projection);
    Double_t x = point.x;
    Double_t y = point.y;

    inputHitsGraph.SetPoint(pIndex++, x, y);
  }
  multiGraph.Add(&inputHitsGraph, "P");

  // Reconstructed tracks graph with text labels

  // Text labels
  std::vector<TText*> labels;
  auto txtOffset  = 1. * canvasX / 3000. * 2;
  auto txtColor   = kBlack;
  auto txtFont    = 43;

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
    for (auto hitIndex : reconstructedTrack) {
      auto hit = hits.at(hitIndex);
      PointP point(hit.position[0],
                   hit.position[1],
                   hit.position[2]);
      point.makeProjection(projection);
      Double_t x = point.x;
      Double_t y = point.y;

      outTrajectoryGraph.SetPoint(pIndex++, x, y);

      // Text labels
      if (!plotLabels) {
        continue;
      }
      auto plotTxt = false;
      if (pIndex % txtStep == 0) plotTxt = true;
      Int_t m = reconstructedTrack.size() / 2;
      if ((2*m < txtStep) && (
          (pIndex == 1)     ||
          (pIndex == m + 1) ||
          (pIndex == 2*m)
         )) {
        plotTxt = true;
      }
      auto trackId = hit.trackId;
      if (plotTxt) {
        auto txt = new TText(x, y, std::to_string(trackId).c_str());
        txt->SetTextAlign(12);
        txt->SetTextColor(txtColor);
        txt->SetTextFont(txtFont);
        txt->SetTextSizePixels(txtSize);
        labels.push_back(txt);
      }
    }
    multiGraph.Add(&outTrajectoryGraph, "PL");
  }
  multiGraph.Draw("A");

  for (const auto txt : labels) {
    txt->Draw();
  }

  auto fname = outPath + "/event_" + std::to_string(eventCounter) +
      namePostfix + ".png";

  canvas.Print(fname.c_str());

  for (auto surfGraph : surfGraphs) {
    delete surfGraph;
  }
  for (auto item : labels) {
    delete item;
  }
}
