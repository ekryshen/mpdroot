// -------------------------------------------------------------------------
// -----                  MpdEctTrackFinderTpc source file             -----
// -----                 Created 28/03/08  by A. Zinchenko             -----
// -------------------------------------------------------------------------

/**  MpdEctTrackFinderTpc.h
 *@author A.Zinchenko <Alexander.Zinchenko@jinr.ru>
 **
 ** Track finder in MPD End-Cap Tracker (ECT) using seeds from TPC
 **/

#include "MpdEctTrackFinderTpc.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanHit.h"
#include "MpdEctKalmanTrack.h"
#include "MpdStrawendcapGeoPar.h"
#include "MpdStrawendcapPoint.h"
#include "MpdDchPoint.h"
#include "MpdTpcKalmanTrack.h"

#include "FairGeoNode.h"
#include "MpdMCTrack.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include <TGeoManager.h>
#include <TGeoTube.h>
#include "TMath.h"
#include "TVector2.h"
#include <TRandom.h>
#include <TFile.h>
#include <iostream>
//#include <vector>

using std::cout;
using std::endl;
//using std::set;

const Double_t MpdEctTrackFinderTpc::fgkChi2Cut = 10; //10 //20; //100;

//__________________________________________________________________________
MpdEctTrackFinderTpc::MpdEctTrackFinderTpc(const char *name, Int_t iVerbose )
  :FairTask(name, iVerbose)
{
  fKHits = new TClonesArray("MpdKalmanHit", 100);
  fTracks = new TClonesArray("MpdEctKalmanTrack", 100);
  fTrackCand = new TClonesArray("MpdEctKalmanTrack", 100);
  fHistoDir = 0x0;
  fhLays = new TH1F("hLays0","ECT layers",fgkNlays+1,0,fgkNlays+1);
  fLayPointers = 0x0;
  fDetType = 1;
  fMirror = kTRUE; //kFALSE;
  fExact = 0;
  if (fExact) fMirror = kFALSE;
}


//__________________________________________________________________________
MpdEctTrackFinderTpc::~MpdEctTrackFinderTpc()
{
  //delete fKHits;
  //delete fTracks;
  //delete fTrackCand;
  delete [] fLayPointers;
  delete fhLays;
}

//__________________________________________________________________________
InitStatus MpdEctTrackFinderTpc::Init()
{
  ReInit();
  if (fDetType == 1) InitGeo();

  TGeoVolume *inW = gGeoManager->GetVolume("tpc01InWall");
  TGeoTube *tube = (TGeoTube*) inW->GetShape();
  fZtpc = tube->GetDZ();

  return kSUCCESS;
}

//__________________________________________________________________________
InitStatus MpdEctTrackFinderTpc::ReInit()
{
  fEctPoints = 0x0; //(TClonesArray *) FairRootManager::Instance()->GetObject("STRAWPoint");
  fEctHits =(TClonesArray *) FairRootManager::Instance()->GetObject("STRAWPoint");
  if (fEctHits == 0x0) {
    fEctHits =(TClonesArray *) FairRootManager::Instance()->GetObject("TgemPoint");
    fDetType = 0;
  }
  if (fEctHits == 0x0) {
    fEctHits =(TClonesArray *) FairRootManager::Instance()->GetObject("DchPoint");
    fDetType = 2;
  }
  if (fEctHits == 0x0) Fatal("ReInit", "No ECT points found !!!");
  fTpcTracks =(TClonesArray *) FairRootManager::Instance()->GetObject("TpcKalmanTrack");
  fMCTracks =(TClonesArray *) FairRootManager::Instance()->GetObject("MCTrack");
  //fSTSTrackMatch = (TClonesArray*) FairRootManager::Instance()->GetObject("STSTrackMatch");
  //fPrimVtx =  (FairVertex *) FairRootManager::Instance() ->GetObject("PrimaryVertex");

  FairRootManager::Instance()->Register("EctTrack", "Ect", fTracks, kTRUE);
  FairRootManager::Instance()->Register("EctKalmanHit", "EctHit", fKHits, kFALSE);

  fNPass = 2; //2;
  return kSUCCESS;
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::Reset() 
{
  ///
  cout << " MpdEctTrackFinderTpc::Reset  " << endl;

  //fKHits->Clear("C");
  fKHits->Delete();
  fTracks->Delete();
  fTrackCand->Delete();
  if (fLayPointers) delete [] fLayPointers;
  fLayPointers = 0x0;
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::SetParContainers()
{
  // Get run and runtime database
  FairRunAna* run = FairRunAna::Instance();
  if ( ! run ) Fatal("SetParContainers", "No analysis run");

  FairRuntimeDb* db = run->GetRuntimeDb();
  if ( ! db ) Fatal("SetParContainers", "No runtime database");

  // Get ECT geometry parameter container
  db->getContainer("MpdStrawendcapGeoPar");
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::Finish()
{
  //Write();
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::Exec(Option_t * option)
{

  static int eventCounter = 0;    
  cout << " EctRec event " << ++eventCounter << endl;

  Reset();

  // Create Kalman hits
  MakeKalmanHits();
  if (fKHits->GetEntriesFast() == 0) return;

  for (Int_t i = 0; i < fNPass; ++i) {
    fTrackCand->Delete();
    GetTrackSeeds(i);

    Int_t nHitsOk = 0, nHits = fKHits->GetEntriesFast();
    for (Int_t j = 0; j < nHits; ++j) {
      MpdKalmanHit *hit = (MpdKalmanHit*) fKHits->UncheckedAt(j);
      if (hit->GetFlag() >= 0) ++nHitsOk;
    }
    
    cout << "  Total number of hits for tracking: " << nHits << ", good: " << nHitsOk << endl;
    cout << "  Total number of track seeds: " << fTrackCand->GetEntriesFast() << endl;

    DoTracking(i);
    RemoveDoubles(); // remove double (clone) tracks
    StoreTracks(i);
    cout << "  Total number of found tracks: " << fTracks->GetEntriesFast() << endl;
    //if (i != fNPass - 1) ExcludeHits(); // exclude used hits
    ExcludeHits(); // exclude used hits
  }
  SelectTracks(0); // do track selection and compute shared hit multiplicities
  AddHits(); // add hit objects to tracks
  cout << "  Total number of found tracks: " << fTracks->GetEntriesFast() << endl;
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::InitGeo()
{
  // Initialize detector geometry

  // Create configuration
  cout << " Creating configuration ... " << endl;

  map<TString, FairGeoNode*> layers;
  map<TString, FairGeoNode*> tubes;
  map<TString, FairGeoNode*> modules;

  FairRuntimeDb* rtdb = FairRun::Instance()->GetRuntimeDb();
  // Create mother volume
  FairGeoNode *cave = new FairGeoNode();
  cave->setVolumeType(kFairGeoTopNode);

  MpdStrawendcapGeoPar *geoPar = (MpdStrawendcapGeoPar*) rtdb->getContainer("MpdStrawendcapGeoPar");
  TObjArray* passNodes = geoPar->GetGeoPassiveNodes();
  Int_t nVol = passNodes->GetEntriesFast();

  // Get volumes
  Int_t nTripl = 0, nMod = 0, nTube = 0;
  frMinMax[0] = frMinMax[2] = -1.;
  for (Int_t i = 0; i < nVol; ++i) {
    FairGeoNode* passVol = (FairGeoNode*) (passNodes->UncheckedAt(i));
    TString name = passVol->GetName();
    TArrayD* params = passVol->getParameters();
    //cout << i << " " << name << " " << (*params)[0] << " " << (*params)[1] << " " 
    // << (*params)[2] << " " << passVol->GetListOfDaughters()->GetEntriesFast() << endl;
    if (name.Contains("layer")) {
      layers.insert(pair<TString, FairGeoNode*>(name,passVol));
      if (frMinMax[0] < 0) {
	frMinMax[0] = (*params)[0];
	frMinMax[1] = (*params)[1];
      }
      Ssiz_t pos = name.First("#") + 1;
      nTripl = TMath::Max (nTripl, (TString(name(pos,name.Length()))).Atoi());
    } else if (name.Contains("tube")) {
      tubes.insert(pair<TString, FairGeoNode*>(name,passVol));
      if (frMinMax[2] < 0) frMinMax[2] = (*params)[1] * 10.; // mm
      Ssiz_t pos = name.First("#") + 1;
      nTube = TMath::Max (nTube, (TString(name(pos,name.Length()))).Atoi());
    } else if (name.Contains("mod")) {
      modules.insert(pair<TString, FairGeoNode*>(name,passVol));
      Ssiz_t pos = name.First("#") + 1;
      nMod = TMath::Max (nMod, (TString(name(pos,name.Length()))).Atoi());
    }
    /*
    CbmGeoTransform *transf = passVol->getLabTransform();
    CbmGeoVector trans = transf->getTranslation();
    //cout << trans << endl;
    CbmGeoRotation rot = transf->getRotMatrix();
    //rot.print();
    passVol->setMother(cave);
    passVol->calcLabTransform();
    */
  }
  //fNofLays = nTripl * 3;
  //fNofLays *= 2;
  fNofTubes = nTube;

  // Set mother volumes and find out views ordering
  map<TString,FairGeoNode*>::iterator it;
  Int_t ishft[3] = {0};
  Double_t zzz[3] = {0.};
  for (it = layers.begin(); it != layers.end(); ++it) {
    FairGeoNode *lay = it->second;
    lay->setMother(cave);
    TString name = lay->GetName();
    if (name == "stt01layerradial#1") zzz[0] = lay->getLabTransform()->getTranslation().getZ();
    else if (name == "stt01layerleft#1") zzz[1] = lay->getLabTransform()->getTranslation().getZ();
    else if (name == "stt01layerright#1") zzz[2] = lay->getLabTransform()->getTranslation().getZ();
    //cout << lay->getLabTransform()->getTranslation() << endl;
    //lay->getLabTransform()->getRotMatrix().print();
  }
  TMath::Sort(3,zzz,ishft,kFALSE);

  // Get transforms
  for (it = tubes.begin(); it != tubes.end(); ++it) {
    FairGeoNode *tube = it->second;
    TString layName = tube->GetName();

    //if (layName.Contains("stt03") || layName.Contains("stt04")) continue; // skip z<0 for the moment !!!

    layName.ReplaceAll("tube","layer");
    Ssiz_t pos = layName.First("#") + 1;
    Int_t group, itube, ind = 0;
    sscanf(&layName[4],"%d",&group);
    --group;
    itube = (TString(layName(pos,layName.Length()))).Atoi() - 1;
    //cout << lay << " " << module << " " << itube << endl;
    if (itube >= fgkNtube) { 
      cout << " Too many tubes !!! " << itube << endl;
      exit(0);
    }
    if (layName.Contains("left")) ind = 1;
    else if (layName.Contains("right")) ind = 2;
    // Loop over triplets
    for (Int_t i = 0; i < nTripl; ++i) {
      layName.Remove(pos);
      layName += (i+1);
      tube->setMother(layers[layName]);
      Int_t ilay = i * 3 + group * nTripl * 3 + ishft[ind];
      if (ilay >= fgkNlays) {
	cout << " Too many layers !!! " << ilay << endl;
	exit(0);
      }
      //cout << tube->GetName() << " " << lay << endl;
      FairGeoTransform *tran = tube->calcLabTransform();
      //cout << tran->getTranslation() << endl;
      //tran->getRotMatrix().print();
      fTransf[ilay][itube] = tran;
      //transf[ilay][itube] = layers[lay]->getLabTransform();
      if (itube == 0) fZplanes[ilay] = layers[layName]->getLabTransform()->getTranslation().getZ();
    }
  } // for (it = tubes.begin();

}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::MakeKalmanHits()
{
  /// Create Kalman hits from ECT points.

  if (fDetType == 0) {
    // TGEM
    // JANO TGEM removed 2021.08.16 - MakeKalmanHitsTgem(); // TGEM
    return;
  } else if (fDetType == 2) {
    // DCH
    MakeKalmanHitsDch(); // DCH
    return;
  }

  Int_t nHits = fEctHits->GetEntriesFast(), lay = 0, nKH = 0, itube = 0;
  Double_t phi, coord, errR = 0.2, errRphi = 0.02; //r, 2mm in R, 200um in R-Phi
  /*Double_t firstHit[100][1000] = {{0},{0}};
  Int_t indxTube[100][1000];
  for (Int_t i = 0; i < 100; ++i) {
    for (Int_t j = 0; j < 1000; ++j) indxTube[i][j] = -1;
    }*/
  Double_t firstHit[fgkNlays+1][fgkNtube] = {{0},{0}};
  Int_t indxTube[fgkNlays+1][fgkNtube];
  Int_t iend = fgkNlays + 1;
  for (Int_t i = 0; i < iend; ++i) {
    for (Int_t j = 0; j < fgkNtube; ++j) indxTube[i][j] = -1;
  }

  for (Int_t ih = 0; ih < nHits; ++ih) {
    MpdStrawendcapPoint *h = (MpdStrawendcapPoint*) fEctHits->UncheckedAt(ih);

    //if (h->GetZ() < 0) continue; // !!! exclude one side for now

    //if (h->GetTrackID() != 1063) continue; // !!!
    phi = h->GetPhi(); // tube Phi
    lay = h->GetDetectorID() / 1000; // + 1; 
    if (h->GetZ() < 0) lay += fgkNlays2;;
    itube = h->GetDetectorID() % 1000; // tube number
    // Extrapolate track to Z = Ztube
    Double_t dZ = h->GetZ() - h->GetTrackZ();
    Double_t dt = 0.; // dZ / h->GetPz();
    if (TMath::Abs(h->GetPz()) > 1.e-6 && h->GetPz() * h->GetZ() > 0) dt = dZ / h->GetPz();
    Double_t xNew = h->GetTrackX() + dt * h->GetPx();
    Double_t yNew = h->GetTrackY() + dt * h->GetPy();
    //cout << dZ << " " << h->GetTrackX() << " " << xNew << " " << h->GetTrackY() << " " << yNew << " " << lay << endl;
    //Double_t zNew = h->GetTrackZ() + dt * h->GetPz(); // just for cross-check
    // Transform to the rotated coordinate system
    Double_t cosSin[2] = {TMath::Cos(phi), TMath::Sin(phi)}; // phi - tube direction
    //Double_t xLoc = h->GetX() * cosPhi + h->GetY() * sinPhi; // cross-check
    //Double_t yLoc = -h->GetX() * sinPhi + h->GetY() * cosPhi;
    //Double_t xRot = xNew * cosSin[0] + yNew * cosSin[1];
    Double_t yRot = -xNew * cosSin[1] + yNew * cosSin[0];
    //Double_t xLoc = (xNew - h->GetX()) * cosPhi + (yNew - h->GetY()) * sinPhi;
    //Double_t yLoc = -(xNew - h->GetX()) * sinPhi + (yNew - h->GetY())  * cosPhi;
    //r = xNew * xNew + yNew * yNew;
    //r = TMath::Sqrt (r);
    //r = TMath::Abs(xLoc);
    //r = xRot;
    //cout << xRot << " " << yRot << " " << r << " " << h->GetPz() << endl;
    coord = yRot;
    Double_t yWire = -h->GetX() * cosSin[1] + h->GetY() * cosSin[0];
    Double_t dY = TMath::Abs (yWire - coord);

    // Add error                                            
    Double_t dRphi = 0, dR = 0;
    gRandom->Rannor(dRphi,dR);
    Double_t meas[2] = {coord+dRphi*errRphi, 0};
    Double_t err[2] = {errRphi, errR};
    MpdKalmanHit *hit = new ((*fKHits)[nKH]) MpdKalmanHit(lay*1000000+itube, 1, MpdKalmanHit::kFixedZ, 
							  meas, err, cosSin, 9., h->GetZ(), ih);
    //hit->SetXY(h->GetX(), h->GetY());
    // Add second measurement - just for test at the moment
    //!!!
    //hit->SetNofDim(2);
    //!!!

    // Check if multiple hits per tube
    if (1 && indxTube[lay][itube] >= 0) {
      // Multiple hits
      MpdKalmanHit *hitOld = (MpdKalmanHit*) fKHits->UncheckedAt(indxTube[lay][itube]); // previous hit
      if (dY < firstHit[lay][itube]) {
	// Closer to anode wire - remove previous hit
	firstHit[lay][itube] = dY;
	const TArrayI *inds = hitOld->Index();
	Int_t nOld = inds->GetSize();
	for (Int_t io = 0; io < nOld; ++io) hit->SetIndex((*inds)[io]);
	fKHits->RemoveAt(indxTube[lay][itube]);
	if (fMirror) fKHits->RemoveAt(indxTube[lay][itube]+1);
	indxTube[lay][itube] = nKH;
      } else {
	hitOld->SetIndex(hit->GetIndex());
	fKHits->RemoveAt(nKH); // remove last hit
      }
    } else {
      firstHit[lay][itube] = dY;
      indxTube[lay][itube] = nKH;
    }
    if (fMirror && fKHits->UncheckedAt(nKH)) {
      // Add mirror hit
      // wire point transverse coordinate in rotated system:
      //yRot = -h->GetX() * cosSin[1] + h->GetY() * cosSin[0];
      //Double_t dy = hit->GetRphi() - yRot;
      Double_t dy = hit->GetMeas(0) - yWire;
      MpdKalmanHit *hitM = new ((*fKHits)[++nKH]) MpdKalmanHit(*hit);
      //hitM->SetRphi(yRot-dy);
      hitM->SetMeas(0, yWire-dy);
      hitM->SetMirror();
    }
    ++nKH;
  }
  if (nKH == 0) return;
  fKHits->Compress();
  fKHits->Sort(); // in descending order in abs(Z)
  Int_t layMax = ((MpdKalmanHit*)fKHits->First())->GetLayer();
  cout << " Max layer = " << layMax << " " << fKHits->GetEntriesFast() << endl;

  fLayPointers = new Int_t [layMax+1];
  fhLays->Reset();
  /*
  Int_t ipos = 0;
  for (Int_t i = 0; i <= layMax; ++i) {
    //cout << i << " " << fhLays->GetCellContent(i+1,0) << endl;
    //if (ipos) cout << ((TpcLheHit*)fHits->UncheckedAt(ipos))->GetLayer() << " "
    //     << ((TpcLheHit*)fHits->UncheckedAt(ipos-1))->GetLayer() << endl;
    Int_t cont = TMath::Nint (fhLays->GetCellContent(i+1,0));
    if (cont == 0) {
      fLayPointers[i] = -1;
      continue;
    }
    fLayPointers[i] = ipos;
    ipos += cont;
  }
  */
  for (Int_t i = 0; i <= layMax; ++i) fLayPointers[i] = -1;
  nKH = fKHits->GetEntriesFast();
  for (Int_t i = 0; i < nKH; ++i) {
    lay = ((MpdKalmanHit*) fKHits->UncheckedAt(i))->GetLayer();
    fhLays->Fill(lay+0.1);
    if (fLayPointers[lay] < 0) fLayPointers[lay] = i;
  }
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::MakeKalmanHitsDch()
{
  /// Create Kalman hits from DCH points.

  fhLays->Reset();
  Int_t nHits = fEctHits->GetEntriesFast(), layMax = 0, lay = 0, nKH = 0, itube = 0;
  if (nHits == 0) return;
  static const Double_t errR = 0.01, errRphi = 0.01; // 100um in R, 100um in R-Phi
  //static const Double_t errR = 0.01, errRphi = 0.02; // 100um in R, 200um in R-Phi
  Double_t firstHit[100][1000] = {{0},{0}};
  Int_t indxTube[100][1000];
  for (Int_t i = 0; i < 100; ++i) {
    for (Int_t j = 0; j < 1000; ++j) indxTube[i][j] = -1;
  }

  for (Int_t ih = 0; ih < nHits; ++ih) {
    MpdDchPoint *h = (MpdDchPoint*) fEctHits->UncheckedAt(ih);
    if (h->GetZ() < 0) continue; // !!! exclude one side for now
    lay = h->GetDetectorID();

    Double_t cosSin[2] = {TMath::Cos(h->GetPhi()), TMath::Sin(h->GetPhi())}; // phi - wire direction
    Double_t u = -h->GetX() * cosSin[1] + h->GetY() * cosSin[0];

    // Add error                                            
    Double_t dRphi = 0, dR = 0;
    gRandom->Rannor(dRphi,dR);
    Double_t meas[2] = {u+dRphi*errRphi, 0};
    Double_t err[2] = {errRphi, errR};
    //MpdKalmanHit *hit = 
    new ((*fKHits)[nKH]) MpdKalmanHit(h->GetDetectorID()*1000000, 1, 
      MpdKalmanHit::kFixedZ, meas, err, cosSin, 9., h->GetZ(), ih);
    // Add second measurement - just for test at the moment
    //!!!
    //hit->SetNofDim(2);
    //!!!
    //lay = hit.GetLayer();
    layMax = TMath::Max (lay, layMax);
    //fhLays->Fill(lay+0.1);

    // Check if multiple hits per tube
    if (0 && indxTube[lay][itube] >= 0) {
      // Multiple hits
      if (TMath::Abs(u) < firstHit[lay][itube]) {
	// Closer to anode wire - remove previous hit
	firstHit[lay][itube] = TMath::Abs(u);
	fKHits->RemoveAt(indxTube[lay][itube]);
	indxTube[lay][itube] = nKH;
      } else fKHits->RemoveAt(nKH); // remove last hit
    } else {
      firstHit[lay][itube] = TMath::Abs(u);
      indxTube[lay][itube] = nKH;
      fhLays->Fill(lay+0.1);
    }
    ++nKH;
  }
  fKHits->Compress();
  fKHits->Sort(); // in descending order in abs(Z)
  cout << " Max layer = " << layMax << " " << fKHits->GetEntriesFast() << " " 
       << ((MpdKalmanHit*)fKHits->First())->GetDist() << endl;

  fLayPointers = new Int_t [layMax+1];
  for (Int_t i = 0; i <= layMax; ++i) fLayPointers[i] = -1;
  nKH = fKHits->GetEntriesFast();
  for (Int_t i = 0; i < nKH; ++i) {
    lay = ((MpdKalmanHit*) fKHits->UncheckedAt(i))->GetLayer();
    if (fLayPointers[lay] < 0) fLayPointers[lay] = i;
  }
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::GetTrackSeeds(Int_t iPass)
{
  /// Build ECT track seeds from TPC tracks (apply some angular selection cuts?)

  Int_t nTpcTracks = fTpcTracks->GetEntriesFast();
  cout << " Seed tracks: " << nTpcTracks << endl;
 
  Int_t nCand = 0;
  for (Int_t itr = 0; itr < nTpcTracks; ++itr) {
    MpdTpcKalmanTrack *tpc = (MpdTpcKalmanTrack*) fTpcTracks->UncheckedAt(itr);
    if (tpc == 0x0 || tpc->GetUniqueID() == 1) continue;
    if (TMath::Abs((*tpc->GetParamAtHit())(3,0)) < 0.3) continue; // track in central region
    MpdKalmanTrack tmp(*tpc);
    tmp.SetParam(*tpc->GetParamAtHit());
    if (tmp.Momentum() < 0.07) continue; // !!! track with low momentum

    //if (tpc->GetParam(3) < 0) continue; // !!! going in backward hemisphere

    MpdEctKalmanTrack *track = new ((*fTrackCand)[nCand++]) MpdEctKalmanTrack(itr, *tpc);
    //cout << nCand-1 << " " << track->GetTrackID() << endl;
    //cout << track->GetHits()->GetEntriesFast() << endl;
    //track->GetParam()->Print();
    //track->GetCovariance()->Print();
    //track->GetWeight()->Print();
    /*AZ 27-12-2011
    Double_t distFirst = ((MpdKalmanHit*) track->GetHits()->First())->GetDist();
    Double_t distLast = ((MpdKalmanHit*) track->GetHits()->Last())->GetDist();
    track->SetPos(TMath::Min(distFirst,distLast));
    track->SetPosNew(track->GetPos());
    Bool_t ok = MpdKalmanFilter::Instance()->Refit(track, -1); // from last to first hit

    //if (!ok || track->GetParamNew(3) < 0) { 
    if (!ok) { 
      // propagation failure or going in backward hemisphere
      fTrackCand->RemoveLast();
      --nCand;
      continue;
    }
    */
    //AZ 27-12-2011
    /*track->SetPos(((MpdKalmanHit*) track->GetHits()->First())->GetDist());
    track->SetPosNew(track->GetPos());
    track->SetLength(track->GetLengAtHit());*/
    //
    track->GetHits()->Clear();
    track->SetChi2(0.);
  }
  cout << " Number of ECT track candidates: " << nCand << endl;
}  
  
//__________________________________________________________________________
void MpdEctTrackFinderTpc::DoTracking(Int_t iPass)
{
  /// Run Kalman tracking
  
  Int_t nCand = fTrackCand->GetEntriesFast(), iok = 0;
 
  for (Int_t i = 0; i < nCand; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTrackCand->UncheckedAt(i);
    //Int_t id = track->GetTrackID();
    //cout << " Track seed No. " << i << ", ID: " << track->GetTrackID() << endl;
    //if (track->GetTrackID() != 77) continue;
    //(*(track->GetParamNew()))(4,0) = -0.5; // just for check
    /*  "Standalone" tracking
    for (Int_t k = 0; k < 5; ++k) {
      for (Int_t j = k; j < 5; ++j) {
	if (j == k) (*track->GetWeight())(k,j) /= 100.;
	else (*track->GetWeight())(k,j) = (*track->GetWeight())(j,k) = 0.;
      }
    }
    */
    
    // Propagate to TPC end-plate
    MpdKalmanHit hit;
    hit.SetType(MpdKalmanHit::kFixedZ);
    //Double_t zEnd = 150.; // get it from geometry !!!
    //hit.SetPos(TMath::Sign(zEnd,track->GetParam(3)));
    hit.SetPos(TMath::Sign(fZtpc,track->GetParam(3)));
    //Double_t thick = 0.06; // material thickness in rad. lengths
    Double_t thick = 0.20; // material thickness in rad. lengths - v7
    MpdKalmanFilter::Instance()->PropagateToHit(track, &hit, kTRUE);
    PassWall(track,thick);

    // Start ECT tracking from different layers to account for missing hits
    MpdEctKalmanTrack tracks[10];
    Int_t layMax = 1; // TGEM or DCH
    //Int_t layMax = 2; // TGEM or DCH
    if (fDetType == 1) layMax = 10; // ECT

    for (Int_t j = 0; j < layMax; ++j) {
      tracks[j] = *track;
      //cout << track->GetParamNew(0) << endl;
      //cout << i << " " << lay << " " << tracks[lay].GetNofHits() << " " << tracks[lay].GetChi2() << " " << tracks[lay].GetParam(0) << endl;
      //Int_t lay = 1;
      Int_t lay = j + 1;
      if (hit.GetPos() < 0) lay += fgkNlays2;
      if (j > 0 && tracks[j-1].GetNofHits() > 0) 
	lay = ((MpdKalmanHit*)tracks[j-1].GetHits()->First())->GetLayer() + 1;
      iok = RunKalmanFilter(&tracks[j], lay);
      //iok = RunKalmanFilter(&tracks[lay], 0);
      //cout << i << " " << lay << " " << tracks[j].GetNofHits() << " " << tracks[j].GetChi2() << endl;
    }

    // Select the best track (with max number of hits)
    Int_t nHitsMax = tracks[0].GetNofHits(), iMax = 0;
    for (Int_t j = 1; j < layMax; ++j) {
      Int_t nhits = tracks[j].GetNofHits();
      if (nhits > nHitsMax) {
	nHitsMax = nhits;
	iMax = j;
      } else if (nhits == nHitsMax) {
	if (tracks[j].GetChi2() < tracks[iMax].GetChi2()) {
	  iMax = j;
	  nHitsMax = nhits;
	}
      }
    }
    //*track = tracks[iMax];
    fTrackCand->RemoveAt(i);
    new ((*fTrackCand)[i]) MpdEctKalmanTrack(tracks[iMax]);

    if (0 && iok == 0 && track->GetNofHits() < 0) {
      //*
      MpdKalmanHit *pHit = (MpdKalmanHit*) track->GetHits()->UncheckedAt(0);
      Double_t dir = TMath::Sign(1.,track->GetPosNew()); 
      track->SetPos(pHit->GetPos()-dir);
      track->SetPosNew(track->GetPos());
      //for (Int_t j = 0; j < 5; ++j) track->SetParam(j,track->GetParam1()[j]);
      //track->SetParamNew(*track->GetParam());
      //MpdKalmanHitZ hitTmp(*pHit);
      //hitTmp.SetZ(track->GetPosNew()-dir);
      //MpdKalmanFilter::Instance()->PropagateToHit(track,&hitTmp);
      MpdKalmanFilter::Instance()->Refit(track,1);
      //*/
      /*
      track->Print("");
      Double_t dir = TMath::Sign(1.,track->GetPosNew()); 
      MpdKalmanHitZ *hitZ = (MpdKalmanHitZ*) track->GetHits()->UncheckedAt(0);
      MpdKalmanHitZ hitTmp(*hitZ);
      hitTmp.SetZ(track->GetPosNew()+dir);
      MpdKalmanFilter::Instance()->PropagateToHit(track,&hitTmp);
      track->Print("");
      MpdKalmanFilter::Instance()->Refit(track,-1);
      track->Print("");
      hitTmp.SetZ(track->GetPosNew()-dir);
      MpdKalmanFilter::Instance()->PropagateToHit(track,&hitTmp);
      track->Print("");
      MpdKalmanFilter::Instance()->Refit(track,1);
      track->Print("");
      */
    }
  } // for (Int_t i = 0; i < nCand;
}
    
//__________________________________________________________________________
void MpdEctTrackFinderTpc::PassWall(MpdEctKalmanTrack *track, Double_t thick)
{
  /// Propagate track thru TPC end plate (include MS)

  // Add multiple scattering
  TMatrixDSym *cov = track->Weight2Cov();
  Double_t th = track->GetParamNew(3);
  Double_t cosTh = TMath::Cos(th);
  Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, thick);
  (*cov)(2,2) += (angle2 / cosTh / cosTh );
  (*cov)(3,3) += angle2;
  Int_t iok = 0;
  MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
  track->SetWeight(*cov);
}

//__________________________________________________________________________
Int_t MpdEctTrackFinderTpc::RunKalmanFilter(MpdEctKalmanTrack *track, Int_t layBeg)
{
  /// Run Kalman filter

  Double_t rMin = 49.2, rMax = 110.5; // min and max radial size of TGEM - to be taken elsewhere
  Double_t dX = 0.009; // 0.9% of X0 (TGEM layer thickness)
  if (fDetType == 1) {
    // ECT
    rMin = 49.2;
    rMax = 110.5;
    dX = 0.001;
  } else if (fDetType == 2) {
    // DCH
    rMin = 9.2;
    rMax = 135.4;
    dX = 0.001;
  }
  if (track->GetRadNew() > rMax) return -1;

  //cout << fHits->GetEntriesFast() << endl;
  //Int_t layMax = ((MpdKalmanHit*)fKHits->First())->GetLayer();
  Int_t layMax = fgkNlays2;
  MpdKalmanHit *hitOK = 0x0;
  MpdKalmanTrack::TrackDir trackDir = track->GetDirection();
  //Int_t layBeg = 0, layEnd = -1, dLay = -1, layOK = -1;
  Int_t layEnd = -1, dLay = -1;//, layOK = -1;
  //if (track->GetPosNew() < 0) layEnd = fgkNlays2;
  if (track->GetPosNew() < 0) layMax = 2*fgkNlays2;
  if (trackDir == MpdKalmanTrack::kOutward) {
    layEnd = layMax + 1;
    dLay = 1;
    //if (track->GetPosNew() > 0 && layMax > fgkNlays2) layEnd -= fgkNlays2;
  }
  
  //Int_t indxOK = hits->IndexOf(hitOK);
  //Int_t nHits = hits->GetEntriesFast();
  Int_t miss = 0;
  TMatrixDSym pointWeight(5), pointWeightTmp(5), saveWeight(5);
  TMatrixD param(5,1), paramTmp(5,1);

  Double_t saveZ = 0., saveLeng = 0.;
  //cout << " Starting hit: " << hitOK->GetLayer() << " " << hitOK->GetTrackID() << " " << hitOK->GetUsage() << endl;
  
  for (Int_t lay = layBeg; lay != layEnd; lay+=dLay) {
    //if (lay > 30) rMin = 38.;
    //else rMin = 33.;
    Int_t nLay = GetNofHitsInLayer(lay);
    Int_t indx0 = GetHitsInLayer(lay);
    Double_t dChi2Min = 1.e+6;
    MpdKalmanHit *hitMin = 0x0;
    //cout << " lay, nLay: " << lay << " " << nLay << " " << indx0 << endl;
    Int_t indxBeg = 0, indxEnd = nLay, dIndx = 1;
    /*
    if (trackDir == MpdKalmanTrack::kOutward) {
      // !!! This part is due to the simplified hit merging (digitization) -
      // different hit position inside one layer - should be removed later 
      indxBeg = nLay - 1;
      indxEnd = -1;
      dIndx = -1;
    }
    */
    //for (Int_t indx = 0; indx < nLay; ++indx) {
    Double_t radNew = -1;
    for (Int_t indx = indxBeg; indx != indxEnd; indx+=dIndx) {
      MpdKalmanHit *hit = (MpdKalmanHit*) fKHits->UncheckedAt(indx0+indx);
      // Exclude used hits
      if (hit->GetFlag() < 0) continue;

      // !!! Exact ID match
      Int_t itmp = 0;
      //if (fExact && ((FairMCPoint*)fEctHits->UncheckedAt(hit->GetIndex()))->GetTrackID() != track->GetTrackID()) continue;
      if (fExact && HitMotherId(hit,track->GetTrackID(),itmp) != track->GetTrackID()) continue;

      if (TMath::Abs(hit->GetPos()-track->GetPosNew()) > 1.e-3) {
	//cout << " window:" << /*hit->GetTrackID() <<*/ " " << hit->GetRphi() << " " << track->GetParamNew(0) << " " << hit->GetZ() << " " << track->GetParamNew(1) << endl;
	// Check if the hit within some window (15x25cm for the moment - check!!!)
	//if (TMath::Abs(hit->GetRphi()-track->GetParamNew(0)) > 9) continue;
	//if (TMath::Abs(Proxim(track,hit)-track->GetParamNew(0)) > 15) continue;
	TVector2 dist = GetDistance(track, hit);
	//if (track->GetTrackID() == 1335) cout << ((TpcLheKalmanTrack*)fTpcTracks->UncheckedAt(track->GetTpcIndex()))->GetNofTrHits() << " " << dist.X() << " " << dist.Y() << " " << hit->GetLayer() << endl;
	if (dist.X() > 15.) continue; // distance in transverse to the tube direction
	if (hit->GetNofDim() > 1 && dist.Y() > 25.) continue; // distance in R 
	//*if (hit->GetTrackID() == 186)*/ cout << " Hit: " << ((FairMCPoint*)fEctHits->UncheckedAt(hit->GetIndex()))->GetTrackID() << " " << hit->GetLayer() << " " << ((FairMCPoint*)fEctHits->UncheckedAt(hit->GetIndex()))->GetDetectorID() << " " << hit->GetRphi() << " " << hit->GetR() << " " << hit->GetZ() << " " << dist.X() << " " << dist.Y() << " " << track->GetParamNew(0) << " " << track->GetParamNew(1) << endl;
	//track->Print("");
	//hit->Print("");
	if (!MpdKalmanFilter::Instance()->PropagateToHit(track,hit)) return -1;
	dist = GetDistance(track, hit);
	//if (track->GetTrackID() == 1335) cout << track->GetTrackID() << " " << dist.X() << " " << dist.Y() << " " << endl;
	
	//*if (hit->GetTrackID() == -607)*/ cout << /*hit->GetTrackID() <<*/ " " << hit->GetRphi() << " " << track->GetParamNew(0) << " " << track->GetParamNew(1) << " " << hit->GetZ() << " " << track->GetPosNew() << endl;
	//
	// For debugging
	/*
	  TVector2 dist0 = GetDistance(track, hit);
	  cout << dist0.X() << " ";
	  MpdKalmanHitZ hitDbg = *hit;
	  Double_t xDbg = hit->GetXY(0) * TMath::Cos(hit->GetPhi()) + hit->GetXY(1) * TMath::Sin(hit->GetPhi());
	  Double_t yDbg = -hit->GetXY(0) * TMath::Sin(hit->GetPhi()) + hit->GetXY(1) * TMath::Cos(hit->GetPhi());
	  hitDbg.SetRphi(yDbg);
	  hitDbg.SetR(xDbg);
	  dist = GetDistance(track, &hitDbg);
	  cout << dist.X() << endl;
	*/
	//
	//Double_t radNew = track->GetRadNew();
	if (radNew < 0) radNew = track->GetRadNew();
	//if (radNew < rMin || radNew > rMax) return -1;
	if (radNew < rMin) break; // next layer
	else if (radNew > rMax) return -1;

	//Double_t err = hit->GetRphiErr();
	//if (track->GetNofHits() == 0) hit->SetRphiErr(0.04);

	// Add multiple scattering
	if (TMath::Abs(hit->GetPos()-track->GetPos()) > 1.e-1) PassWall(track, dX);
      }

      if (fDetType == 1) {
	// ECT
	Double_t phi = TMath::ATan2 (track->GetParamNew(1),track->GetParamNew(0));
	Double_t dPhi = MpdKalmanFilter::Instance()->Proxim(phi,hit->GetPhi()) - phi;
	if (TMath::Abs(dPhi) > 1) continue;
      }

      Double_t dChi2 = MpdKalmanFilter::Instance()->FilterHitZ(track,hit,pointWeight,param);
      //if (track->GetTrackID() == 1335) cout << dChi2 << endl;
      //if (track->GetNofHits() == 0) hit->SetRphiErr(err);
      //if (param(3,0) < 0) { cout << " Check: " << param(3,0) << " " << dChi2 << " " << (*fParamNew)(3,0) << " " << hit->GetRphi() << " " << hit->GetZ() << endl; fParamNew->Print();}
      if (dChi2 < dChi2Min) {
      //if (dChi2 < dChi2Min && dist0.X() <= dist.X()) {
      //if (dChi2 < dChi2Min && dist.X() <= 0.2) {
        dChi2Min = dChi2;
        hitMin = hit;
        saveWeight = *track->GetWeight();
        saveZ = track->GetPosNew();
	saveLeng = track->GetLength();
        // temporary storage for the current track
        paramTmp = param;
        pointWeightTmp = pointWeight;
        //cout << " New min dChi2 = " << dChi2 << " " << hitMin->GetRphi() << " " << hitMin->GetR() << endl;
	//cout << track->GetParamNew(0) << " " << track->GetParamNew(1) << endl;
	//cout << hit->GetRphi() << " " << hit->GetZ() << endl;
	//cout << param(0,0) << " " << param(1,0) << endl;
        //paramTmp.Print();
      }
    } //for (Int_t indx = indxBeg; indx != indxEnd;
    Double_t cut = fgkChi2Cut;
    //if (track->GetNofHits() == 0) cut /= 2.;
    //if (dChi2Min < fgkChi2Cut) {
    if (dChi2Min < cut) {
      //layOK = lay;
      hitOK = hitMin;
      track->GetHits()->Add(hitOK);
      //miss = 0;
      // Restore track params at the best hit
      track->SetChi2(track->GetChi2()+dChi2Min);
      saveWeight += pointWeightTmp;
      track->SetWeight(saveWeight);
      track->SetPosNew(saveZ);
      track->SetLength(saveLeng);
      track->SetParamNew(paramTmp);
      //if (track->GetTrackID() == 1335) cout << " *** Adding hit: " << hitOK->GetLayer() << " " << ((FairMCPoint*)fEctHits->UncheckedAt(hitOK->GetIndex()))->GetTrackID() << " " << ((FairMCPoint*)fEctHits->UncheckedAt(hitOK->GetIndex()))->GetDetectorID() << " " << dChi2Min << " " << track->GetChi2() << endl;
      //paramTmp.Print();
      // Check if the accepted hit is the same as the seeded hit
      //if (hitOK->GetLayer() == f2ndHit->GetLayer() && hitOK != f2ndHit) return -1; // abandon track
      if (track->GetNofHits() == 1) track->SetParam1();
      //if (((FairMCPoint*)fEctHits->UncheckedAt(hitOK->GetIndex()))->GetTrackID() == track->GetTrackID()) track->fLengAtHit = saveLeng;
    } else if (fDetType == 1 && radNew > rMin) {
      // ECT
      //++miss;
      //if (miss > 1) return -1;
      // Check whether tracks goes through a tube
      Int_t layM1 = lay - 1;
      if (nLay == 0) {
	// No hits in the layer - extrapolate track
	MpdKalmanHit hitTmp;
	hitTmp.SetType(MpdKalmanHit::kFixedZ);
	hitTmp.SetDist(fZplanes[layM1]);
	if (!MpdKalmanFilter::Instance()->PropagateToHit(track,&hitTmp)) return -1;
      }
      FairGeoVector vec(track->GetParamNew(0)*10,track->GetParamNew(1)*10,track->GetPosNew()*10); // mm
      //cout << vec.getX() << " " << vec.getY() << " " << vec.getZ() << endl;
      Int_t iTube1 = 0, iter = 0, idTubes = 0;//, iTubeOk = 0;
      Int_t iTube = Int_t (TMath::ATan2 (vec.getY(),vec.getX()) / TMath::TwoPi() * fNofTubes);
      if (iTube < 0) iTube += fNofTubes;
      iTube = iTube % fNofTubes;
      FairGeoVector vecLoc = fTransf[layM1][iTube]->transTo(vec);
      //cout << vecLoc.getX() << " " << vecLoc.getY() << " " << vecLoc.getZ() << " " << iTube << endl;
      //cout << transf[lay][0]->getTranslation() << endl;
      //transf[lay][0]->getRotMatrix().print();
      /*if (TMath::Abs(vecLoc.getZ()) > z2) cout << vecLoc << " " << ph << endl;
	if (TMath::Abs(vecLoc.getZ()) > z2) {
	iTube = nTubes / 2;
	vecLoc = transf[lay][iTube]->transTo(vec);
	cout << vecLoc;
	}*/
      Double_t dyMin = 0;
      while (1) {
	FairGeoVector vecLoc1 = fTransf[layM1][(iTube+1)%fNofTubes]->transTo(vec);
	Double_t dy = vecLoc.getY() - vecLoc1.getY();
	idTubes = TMath::Nint (vecLoc.getY() / dy);
	iTube1 = iTube + idTubes;
	if (iTube1 < 0) iTube1 += fNofTubes;
	iTube1 = iTube1 % fNofTubes;
	vecLoc1 = fTransf[layM1][iTube1]->transTo(vec);
	dyMin = TMath::Abs(vecLoc1.getY());
	//iTubeOk = iTube1;
	//cout << vecLoc1.getX() << " " << vecLoc1.getY() << " " << vecLoc1.getZ() << " " << iTubeOk << " " << dyMin << " " << iter << endl;
	if (dyMin <= frMinMax[2] || (iter && TMath::Abs(idTubes) < 2 && TMath::Abs(dyMin/dy) < 0.51)) break;
	iTube = iTube1;
	vecLoc = vecLoc1;
	++iter;
      }
      //if (dyMin > 4.) ++miss; // too far from the tube
      //if (dyMin < frMinMax[2]+0.5) ++miss; // too far from the tube
      if (dyMin < frMinMax[2]-0.2) ++miss; // too far from the tube
      //cout << dyMin << " " << dChi2Min << endl;
    }
    //cout << " lay, miss: " << lay << " " << miss << " " << dChi2Min << " " << fChi2 << endl;
  } // for (Int_t lay = layOK-1; lay >= 0;
  track->SetMisses(miss);
  return 0;
}

//__________________________________________________________________________
TVector2 MpdEctTrackFinderTpc::GetDistance(MpdEctKalmanTrack *track, MpdKalmanHit *hit)
{
  /// Compute distance between track and hit

  Double_t xTr, yTr;
  if (track->GetType() == MpdKalmanTrack::kEndcap) {
    xTr = track->GetParamNew(0);
    yTr = track->GetParamNew(1);
  } else {
    Double_t rPhi = track->GetParamNew(0);
    Double_t r = track->GetPosNew();
    Double_t ph = rPhi / r;
    xTr = r * TMath::Cos(ph);
    yTr = r * TMath::Sin(ph);
  }
  Double_t phi = hit->GetPhi();
  Double_t cosPhi = TMath::Cos(phi);
  Double_t sinPhi = TMath::Sin(phi);
  // Transform track coordinates to local tube coordinates
  Double_t xLoc = xTr * cosPhi + yTr * sinPhi;
  Double_t yLoc = -xTr * sinPhi + yTr * cosPhi;
  //Double_t xLoc = (xTr - hit->GetXY(0)) * cosPhi + (yTr - hit->GetXY(1)) * sinPhi;
  //Double_t yLoc = -(xTr - hit->GetXY(0)) * sinPhi + (yTr - hit->GetXY(1)) * cosPhi;
  TVector2 dist(TMath::Abs(yLoc-hit->GetMeas(0)), TMath::Abs(xLoc-hit->GetMeas(1)));
  //cout << dist.X() << " " << dist.Y() << endl;
  return dist;
}

//__________________________________________________________________________
Double_t MpdEctTrackFinderTpc::Proxim(Double_t phi0, Double_t phi)
{
  /// Adjust angle phi to be "around" phi0 - to avoid discontinuity around +- Pi

  Double_t dPhi = phi0 - phi;
  if (TMath::Abs(dPhi) > TMath::Pi()) phi += TMath::Pi() * 2 * TMath::Sign(1.,dPhi);
  return phi;
}

//__________________________________________________________________________
Double_t MpdEctTrackFinderTpc::Proxim(MpdKalmanTrack *track, const MpdKalmanHit *hit)
{
  /// Adjust hit coord. R-Phi to be "around" track R-Phi - to avoid 
  /// discontinuity around +- Pi

  Fatal("Proxim", " !!! Not implemented !!!");
  /*
  Double_t hitPhi = hit->GetRphi() / hit->GetR();
  Double_t phi0 = track->GetParamNew(0) / track->GetPosNew();
  return hit->GetR() * Proxim(phi0, hitPhi);
  */
  return 0.;
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::Write()
{
  /// Write

  TFile histoFile("EctRec.root","RECREATE");
  Writedir2current(fHistoDir);
  histoFile.Close();
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::Writedir2current( TObject *obj )
{
  /// Write

  if( !obj->IsFolder() ) obj->Write();
  else{
    TDirectory *cur = gDirectory;
    TDirectory *sub = cur->mkdir(obj->GetName());
    sub->cd();
    TList *listSub = ((TDirectory*)obj)->GetList();
    TIter it(listSub);
    while( TObject *obj1=it() ) Writedir2current(obj1);
    cur->cd();
  }
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::SelectTracks(Int_t iPass)
{
  /// Do track selection and compute shared hit multiplicities

  //if (iPass) return;
  Int_t nFound = fTracks->GetEntriesFast(), nHitsTot = fKHits->GetEntriesFast();
  Int_t *index = new Int_t [nHitsTot];
  for (Int_t i = 0; i < nHitsTot; ++i) index[i] = 0;

  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTracks->UncheckedAt(i);
    Int_t nHits = track->GetNofHits();
    //if (nHits == 0) { fTracks->RemoveAt(i); continue; }
    if (nHits == 0) continue;
    TObjArray *hits = track->GetHits();
    for (Int_t j = 0; j < nHits; ++j) {
      MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(j);
      Int_t ind = fKHits->IndexOf(hit); 
      if (index[ind]) hit->SetFlag(TMath::Abs(hit->GetFlag())+10);
      ++index[ind];
    }
  }
  //fTracks->Compress();

  // Remove ghost tracks (with too many shared hits)
  nFound = fTracks->GetEntriesFast();
  //Int_t *nh = new Int_t [nFound];
  Double_t *etas = new Double_t [nFound];
  Int_t *order = new Int_t [nFound];
  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTracks->UncheckedAt(i);
    //nh[i] = track->GetNofHits();
    etas[i] = TMath::Abs (track->Momentum3().Eta());
  }
  //TMath::Sort(nFound,nh,order,kFALSE); // in ascending order
  TMath::Sort(nFound,etas,order,kTRUE); // in descending order

  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTracks->UncheckedAt(order[i]);
    Int_t nHits = track->GetNofHits();
    TObjArray *hits = track->GetHits();
    Int_t nShared = 0;
    for (Int_t j = 0; j < nHits; ++j) {
      MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(j);
      if (hit->GetFlag()%1000 >= 10) ++nShared;
    }
    if (1.*nShared / nHits >= 0.5) {
      // Remove track
      // Update hit sharing
      for (Int_t j = 0; j < nHits; ++j) {
	MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(j);
	if (hit->GetFlag()%1000 >= 10)  hit->SetFlag(hit->GetFlag()-10);
	else {
	  if (!fMirror) {
	    hit->SetFlag(TMath::Abs(hit->GetFlag())); // use hit again
	  } else {
	    // Check mirror hit 
	    Int_t nLay = GetNofHitsInLayer(hit->GetLayer());
	    Int_t indx0 = GetHitsInLayer(hit->GetLayer());
	    for (Int_t indx = 0; indx < nLay; ++indx) {
	      MpdKalmanHit *h = (MpdKalmanHit*) fKHits->UncheckedAt(indx0+indx);
	      if (h == hit) continue;
	      if (h->GetDetectorID() == hit->GetDetectorID())  {
		if (index[indx0+indx]) break; // do not use hit again
		hit->SetFlag(TMath::Abs(hit->GetFlag())); // use hit again
		h->SetFlag(TMath::Abs(h->GetFlag())); // use hit again
		break;
	      }
	    }
	  }
	}
      }
      fTracks->RemoveAt(order[i]);
      cout << " Removing track: " << nHits << " " << nShared << endl;
    }
  }
  fTracks->Compress();
  //delete [] nh;
  delete [] etas;
  delete [] index;
  delete [] order;
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::RemoveDoubles()
{
  /// Remove double tracks (if number of common hits greater than 50% of hits on track)

  Int_t nFound = fTrackCand->GetEntriesFast(), nOK = 0;
  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTrackCand->UncheckedAt(i);
    if (track == 0x0) continue;
    Int_t nh = track->GetNofHits();
    if (nh == 0) continue;
    //TObjArray *hits = track->GetHits();
    //cout << i << " " << nh << " " << ++nOK << endl;
    for (Int_t j = i + 1; j < nFound; ++j) {
      MpdEctKalmanTrack *track1 = (MpdEctKalmanTrack*) fTrackCand->UncheckedAt(j);
      if (track1 == 0x0) continue;
      Int_t nh1 = track1->GetNofHits();
      if (nh1 == 0) continue;
      //TObjArray *hits1 = track1->GetHits();
      Int_t nc = NofCommonHits(track, track1);
      if (1.*nc/TMath::Min(nh,nh1) < 0.5) continue;
      if (nh > nh1+5) fTrackCand->RemoveAt(j);
      else if (nh < nh1-5) {
        fTrackCand->RemoveAt(i);
        --nOK;
        break;
      } else {
        if (track->GetChi2()/nh > track1->GetChi2()/nh1) {
          fTrackCand->RemoveAt(i);
          --nOK;
          break;
        }
        fTrackCand->RemoveAt(j);
      }
    }
  }
  fTrackCand->Compress();
}

//__________________________________________________________________________
Int_t MpdEctTrackFinderTpc::NofCommonHits(MpdEctKalmanTrack* track, MpdEctKalmanTrack* track1, Int_t dir)
{
  /// Compute number of common hits on 2 tracks

  TObjArray *hits = track->GetHits(), *hits1 = track1->GetHits();

  MpdKalmanHit *hit = (MpdKalmanHit*) hits->First();
  MpdKalmanHit *hit1 = (MpdKalmanHit*) hits1->First();

  Int_t nCom = 0;
  while (hit && hit1) {
    if (hit->GetLayer() > hit1->GetLayer()) {
      hit1 = (MpdKalmanHit*) hits1->After(hit1);
      continue;
    }
    if (hit == hit1) ++nCom;
    hit = (MpdKalmanHit*) hits->After(hit);
  }
  return nCom;
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::AddHits()
{
  /// Add hit objects to tracks and compute number of wrongly assigned hits
  /// (hits with ID different from ID of starting TPC track)

  Int_t nFound = fTracks->GetEntriesFast();
  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTracks->UncheckedAt(i);
    Int_t nHits = track->GetNofHits();
    if (nHits == 0) { fTracks->RemoveAt(i); continue; }
    cout << i << " " << nHits << " " << track->GetTrackID() << " " << track->GetChi2() << " " 
	 << track->Phi()*TMath::RadToDeg() << " " << track->Momentum3().Eta() << " " 
	 << 1./track->GetParam(4) << " " << track->GetMisses() << endl;
    track->SetNofHits(nHits);
    TClonesArray &trHits = *track->GetTrHits();
    TObjArray *hits = track->GetHits();
    Int_t nWrong = 0, nMirr = 0, motherID = track->GetTrackID(), motherID1 = 0;
    // Get track mother ID 
    MpdMCTrack *mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(motherID);
    while (mctrack->GetMotherId() >= 0) {
      motherID = mctrack->GetMotherId();
      mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(mctrack->GetMotherId());
    }

    for (Int_t j = 0; j < nHits; ++j) {
      MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(j);
      new (trHits[j]) MpdKalmanHit(*hit);
      Int_t iproj = (hit->GetLayer() - 1) % 3;
      if (iproj == 0) cout << " R";
      else if (iproj == 1) cout << " U";
      else cout << " V";
      cout << (hit->GetLayer()-1) % fgkNlays2 + 1;
      Int_t hitId = HitMotherId(hit,motherID,motherID1);
      //cout << "-" << motherID1;
      cout << "-" << hitId;
      if (hit->Index()->GetSize() > 1) cout << "-" << hit->Index()->GetSize();
      if (hit->IsMirror()) cout << "(M)";
      if (motherID1 != motherID) ++nWrong;
      else if (hit->IsMirror()) ++nMirr;
      if (hit->GetFlag() > -1) hit->SetFlag(-hit->GetFlag());
    }
    cout << "\n" << nWrong << " " << nMirr << endl;
    track->SetNofWrong(nWrong);
    track->SetMirrors(nMirr);
    track->SetLastLay();
    //track->GetParam()->Print();
  }
  fTracks->Compress();
}

//__________________________________________________________________________
Int_t MpdEctTrackFinderTpc::HitMotherId(MpdKalmanHit* hit, Int_t idM, Int_t &id1)
{
  /// Check if hit has the same mother ID as idM

  Int_t nOver = hit->Index()->GetSize(), idHit = 0;
  for (Int_t i = 0; i < nOver; ++i) {
    FairMCPoint *h = (FairMCPoint*) fEctHits->UncheckedAt(hit->GetIndex(i));
    id1 = idHit = h->GetTrackID();
    // Get point mother ID 
    MpdMCTrack* mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(id1);
    while (mctrack->GetMotherId() >= 0) {
      id1 = mctrack->GetMotherId();
      mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(mctrack->GetMotherId());
    }
    if (id1 == idM) return idHit;
  }
  return idHit;
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::StoreTracks(Int_t iPass)
{
  /// Transfer tracks from fTrackCand to fTracks

  static const Int_t nHitMin = 10, nHitMinAbs = 5;
  static const Double_t etaMin = 1.3, etaMax = 1.5;

  Int_t nFound = fTracks->GetEntriesFast(), nCand = fTrackCand->GetEntriesFast();
  for (Int_t i = 0; i < nCand; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTrackCand->UncheckedAt(i);
    if (fDetType == 1) {
      // Selection criteria for ECT
      if (track->GetNofHits() < nHitMinAbs) continue;
      Double_t eta = TMath::Abs(track->Momentum3().Eta());
      //if (iPass == 0 && track->GetChi2() / (track->GetNofHits()-5) > 3.0) continue;
      //if (iPass == 0 && eta < etaMin) continue;
      if (eta < etaMin) continue;
      if (eta > etaMax && track->GetNofHits() < nHitMin) continue;
      // !!!!
      //if (track->GetMisses() >= track->GetNofHits()) continue;
      // !!!
    } else if (fDetType == 2 && iPass == 0) {
      // Selection criteria for DCH
      if (track->GetNofHits() < 11) continue;
    }
    track->Weight2Cov();
    //SetTrackID(track);
    new ((*fTracks)[nFound++]) MpdEctKalmanTrack(*track);
    //fTpcTracks->RemoveAt(track->GetTpcIndex());
    fTpcTracks->UncheckedAt(track->GetTpcIndex())->SetUniqueID(1); // for one-job running of TPC and ECT tracking 
  }
}

//__________________________________________________________________________
void MpdEctTrackFinderTpc::ExcludeHits()
{
  /// Exclude hits, already used for tracking, from consideration during the next passes

  Int_t nReco = fTracks->GetEntriesFast();
  cout << " nReco: " << nReco << endl;
  for (Int_t i = 0; i < nReco; ++i) {
    MpdKalmanTrack *track = (MpdKalmanTrack*) fTracks->UncheckedAt(i);
    Int_t nhitsKF = track->GetNofHits();
    TObjArray *hits = track->GetHits();
    for (Int_t j = 0; j < nhitsKF; ++j) {
      MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(j);
      hit->SetFlag(-TMath::Abs(hit->GetFlag()));
      if (fMirror) {
	// Exclude mirror hits as well
	Int_t nLay = GetNofHitsInLayer(hit->GetLayer());
	Int_t indx0 = GetHitsInLayer(hit->GetLayer());
	for (Int_t indx = 0; indx < nLay; ++indx) {
	  MpdKalmanHit *h = (MpdKalmanHit*) fKHits->UncheckedAt(indx0+indx);
	  if (h == hit) continue;
	  if (h->GetDetectorID() == hit->GetDetectorID()) { h->SetFlag(-TMath::Abs(h->GetFlag())); break; }
	}
      }
    }
  }
}

ClassImp(MpdEctTrackFinderTpc);
