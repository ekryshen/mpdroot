//--------------------------------------------------------------------
//
// Description:
//      MPD TPC-EMC Matching
//
//
// Author List:
//      D.Peresunko RRCKI 2019
//
//--------------------------------------------------------------------

#include "MpdEmcMatchingKI.h"
#include "MpdEmcClusterKI.h"
#include "MpdEmcGeoUtils.h"
#include "MpdEmcTrackExtrap.h"
#include "MpdKalmanFilter.h"
#include "MpdTpcHit.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdVertex.h"

#include "FairLogger.h"
#include "MpdMCTrack.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "MpdEmcSimParams.h"

#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TMath.h>
#include <TVector3.h>

#include <iterator>

// -----   Default constructor   -------------------------------------------

MpdEmcMatchingKI::MpdEmcMatchingKI() : FairTask("TPC-EMC matching") {}

// -----   Destructor   ----------------------------------------------------

MpdEmcMatchingKI::~MpdEmcMatchingKI() {}
// -------------------------------------------------------------------------

// -----   Public method Init   --------------------------------------------

InitStatus MpdEmcMatchingKI::Init()
{
  // Get RootManager
  FairRootManager* ioman = FairRootManager::Instance();
  if (!ioman) {
    LOG(ERROR) << "MpdEmcMatchingKI::Init: "
               << "RootManager not instantiated!" << endl;
    return kFATAL;
  }

  // Get input array
  fClusterArray = (TObjArray*)ioman->GetObject("EmcCluster");
  if (!fClusterArray) {
    LOG(ERROR) << "MpdEmcMatchingKI::Init: "
               << "No EMC clusters array!" << endl;
    return kERROR;
  }

  fTpcTracks = (TClonesArray*)ioman->GetObject("TpcKalmanTrack");
  if (!fTpcTracks) {
    LOG(ERROR) << "MpdEmcMatchingKI::Init: "
               << "No TPC track array!" << endl;
    return kERROR;
  }

  fvtx = (TClonesArray*)ioman->GetObject("Vertex");
  if (!fvtx) {
    LOG(INFO) << "Vertex not found, assume (0,0,0)" << endl;
  }

  double rMax;
  MpdEmcGeoUtils* geom = MpdEmcGeoUtils::GetInstance();
  geom->GetECALTubeSize(frMin, rMax, fzMax);
  return kSUCCESS;
}

//__________________________________________________________________________

void MpdEmcMatchingKI::Finish() {}

//__________________________________________________________________________

void MpdEmcMatchingKI::Exec(Option_t* opt)
{
  // Main processing engine
  // Prepare list of cluster coordinates
  // Prepare list of coordinates of tracks extrapolated to inner ECAL radius
  // For each cluster find N closest tracks and extrapolate them to the 1/2 depth of ECAL - look if distance decrease
  // For each track find M closest clusters and extrapolate to 1/2 depth of ECAL - check if distance decrease
  // Calculate distance along track and perpendicular to track and (dz,dphi)

  Int_t nClu = fClusterArray->GetEntriesFast();
  Int_t ntpc = fTpcTracks->GetEntriesFast();

  LOG(INFO) << "MpdEmcMatchingKI::Exec started: number of clusters " << nClu << ", number of tracks " << ntpc << endl;

  if (ntpc == 0 || nClu == 0) {
    return;
  }

  CorrectClustersVtx();

  ExtrapolateTracks();

  MakeMatchToClusters();

  //  MakeMatchToTracks() ;
}
//__________________________________________________________________________
void MpdEmcMatchingKI::CorrectClustersVtx()
{
  // Get vertex
  if (!fvtx) { // no vertex, no corrections
    return;
  }

  MpdEmcSimParams* simParams = MpdEmcSimParams::GetInstance();

  TVector3 primVert;
  MpdVertex* vtx = (MpdVertex*)fvtx->First();
  vtx->Position(primVert);
  if (primVert.Z() == 0. || !simParams->CorrectZ()) { // no need to correct
    return;
  }

  Int_t nClu = fClusterArray->GetEntriesFast();
  for (int i = 0; i < nClu; i++) {
    MpdEmcClusterKI* clu = static_cast<MpdEmcClusterKI*>(fClusterArray->UncheckedAt(i));
    clu->CorrectVertex(primVert.Z());
  }
}

//__________________________________________________________________________
void MpdEmcMatchingKI::ExtrapolateTracks()
{
  // Clean containers
  // Sort Cluster extrapolations according to grid
  fTrackPoints.clear();

  MpdKalmanFilter* pKF = MpdKalmanFilter::Instance("KF", "KF");
  MpdKalmanHit hEnd;
  hEnd.SetType(MpdKalmanHit::kFixedR);

  int nTracks = fTpcTracks->GetEntriesFast();
  for (int itrack = 0; itrack < nTracks; itrack++) {
    MpdTpcKalmanTrack* tr = (MpdTpcKalmanTrack*)fTpcTracks->UncheckedAt(itrack);
    if (TMath::Abs((*tr->GetParamAtHit())(1, 0)) > fzMax)
      continue;
    MpdTpcKalmanTrack tr1(*tr);
    tr1.SetParam(*tr1.GetParamAtHit());
    tr1.SetParamNew(*tr1.GetParamAtHit());
    tr1.SetWeight(*tr1.GetWeightAtHit());
    tr1.SetPos(tr1.GetPosAtHit());
    tr1.SetPosNew(tr1.GetPos());
    tr1.SetLength(tr1.GetLengAtHit());

    // Propagate to EMC inner radius
    hEnd.SetPos(frMin);
    if (!pKF->PropagateToHit(&tr1, &hEnd, kTRUE)) {
      continue;
    }
    double phi = tr1.GetParamNew(0) / tr1.GetPosNew();
    double z = tr1.GetParamNew(1);
    double r = tr1.GetPosNew();

    // Remember position at ECAL surface
    double surfX = r * TMath::Cos(phi);
    double surfY = r * TMath::Sin(phi);
    double surfZ = z;

    // Extrapolate further to the reconstruction radius of clusters
    double rClu = MpdEmcGeoUtils::GetInstance()->Rperp(z) ;
    hEnd.SetPos(rClu);
    bool iok = pKF->PropagateToHit(&tr1, &hEnd, kTRUE);
    double depthX, depthY, depthZ;
    if (iok) {
      r = tr1.GetPosNew();
      phi = tr1.GetParamNew(0) / r;
      depthX = r * TMath::Cos(phi);
      depthY = r * TMath::Sin(phi);
      depthZ = tr1.GetParamNew(1);
    } else { // old coordinates
      depthX = surfX;
      depthY = surfY;
      depthZ = surfZ;
    }
    fTrackPoints.emplace_back(itrack, surfX, surfY, surfZ, depthX, depthY, depthZ);
  }
}

//__________________________________________________________________________
void MpdEmcMatchingKI::MakeMatchToClusters()
{
  // Look over track extrapolations and find closest

  Int_t nClu = fClusterArray->GetEntriesFast();
  for (int i = 0; i < nClu; i++) {
    MpdEmcClusterKI* clu = static_cast<MpdEmcClusterKI*>(fClusterArray->UncheckedAt(i));

    std::vector<MpdEmcTrackExtrap>::iterator minTr;
    double minD = 9999.;
    std::vector<MpdEmcTrackExtrap>::iterator a = fTrackPoints.begin();
    while (a != fTrackPoints.end()) {
      double distance = a->Distance(clu);
      if (distance < minD) {
        minD = distance;
        minTr = a;
      }
      a++;
    }

    // Evaluate distances in (phi,z) (alongTr, perpTr)
    if (minD < 999) { // found some tracks
      double dphi, dz;
      minTr->DistanceDphiDz(clu, dphi, dz);
      clu->SetTrackIndex(minTr->GetTrackIndex());
      clu->SetTrackDxDz(dphi, dz);
    } else {
      // Distance to track extrapolation
      // Track index already set to -1 at initialization
      clu->SetTrackDxDz(9999., 9999.);
    }
  }
}
//__________________________________________________________________________
void MpdEmcMatchingKI::MakeMatchToTracks()
{
  // Look over track extrapolations and find closest

  Int_t nClu = fClusterArray->GetEntriesFast();
  std::vector<MpdEmcTrackExtrap>::iterator a = fTrackPoints.begin();
  while (a != fTrackPoints.end()) {
    double minD = 9999.;
    //int icluMin;
    for (int i = 0; i < nClu; i++) {
      MpdEmcClusterKI* clu = static_cast<MpdEmcClusterKI*>(fClusterArray->UncheckedAt(i));
      double distance = a->Distance(clu);
      if (distance < minD) {
        minD = distance;
        //icluMin = i;
      }
    }

    // Evaluate distances in (phi,z) (alongTr, perpTr)
    /*    if(minD<999){ //found some cluster
          a->SetDistance(icluMin) ;
        }
    */
  }
}

ClassImp(MpdEmcMatchingKI)
