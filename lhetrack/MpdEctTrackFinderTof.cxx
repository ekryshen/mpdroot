// -------------------------------------------------------------------------
// -----                  MpdEctTrackFinderTof source file             -----
// -----                 Created 28/07/08  by A. Zinchenko             -----
// -------------------------------------------------------------------------

/**  MpdEctTrackFinderTof.h
 *@author A.Zinchenko <Alexander.Zinchenko@jinr.ru>
 **
 ** Track finder in MPD End-Cap Tracker (ECT) using seeds built from ECT, 
 ** ETOF and vertex points
 **/

#include "MpdEctTrackFinderTof.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanHit.h"
#include "MpdEctKalmanTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdStrawendcapGeoPar.h"
#include "MpdStrawendcapPoint.h"
#include "MpdTofHit.h"
#include "MpdTofPoint.h"
#include "MpdTofUtils.h"
#include "MpdVertex.h"
#include "MpdKfPrimaryVertexFinder.h"

#include "FairMCEventHeader.h"
#include "FairGeoNode.h"
#include "MpdMCTrack.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "TMath.h"
#include "TFile.h"
//#include "TLorentzVector.h"
#include "TVector2.h"
//#include "TClonesArray.h"
#include <TRandom.h>

#include <iostream>
#include <set>

using std::cout;
using std::endl;
//using std::vector;

const Double_t MpdEctTrackFinderTof::fgkChi2Cut = 10; //10; //20; //100;
FILE *lun1 = 0x0; //fopen("error1.dat","w");
FILE *lun2 = 0x0; //fopen("match.dat","w");

//__________________________________________________________________________
MpdEctTrackFinderTof::MpdEctTrackFinderTof(const char *name, Int_t iVerbose )
  :FairTask(name, iVerbose),
   fNPass(1), 
   fTpc(kFALSE),
   fHistoDir(0x0),
   fExact(0)

{
  //fKHits = new TClonesArray("MpdKalmanHit", 100);
  fTrackCand = new TClonesArray("MpdEctKalmanTrack", 100);
  fhLays = new TH1F("hLays1","ECT layers",fgkNlays+1,0,fgkNlays+1);
  fLayPointers = 0x0;
}


//__________________________________________________________________________
MpdEctTrackFinderTof::~MpdEctTrackFinderTof()
{
  //delete fKHits;
  //delete fTracks;
  //delete fTrackCand;
  delete [] fLayPointers;
  delete fhLays;
}

//__________________________________________________________________________
InitStatus MpdEctTrackFinderTof::Init()
{
  ReInit();
  InitGeo();
  return kSUCCESS;  
}

//__________________________________________________________________________
InitStatus MpdEctTrackFinderTof::ReInit()
{
  if (fTpc) {
    fKHits =(TClonesArray *) FairRootManager::Instance()->GetObject("EctKalmanHit");
    fTracks =(TClonesArray *) FairRootManager::Instance()->GetObject("EctTrack");
    //fTpcTracks =(TClonesArray *) FairRootManager::Instance()->GetObject("TpcKalmanTrack");
  } else {
    fKHits = new TClonesArray("MpdKalmanHit", 100);
    fTracks = new TClonesArray("MpdEctKalmanTrack", 100);
    FairRootManager::Instance()->Register("EctTrack", "Ect", fTracks, kTRUE);
  }
  fTpcTracks =(TClonesArray *) FairRootManager::Instance()->GetObject("TpcKalmanTrack");
  fEctHits =(TClonesArray *) FairRootManager::Instance()->GetObject("STRAWPoint");
  fTofHits =(TClonesArray *) FairRootManager::Instance()->GetObject("ETOFHit");
  fTofPoints =(TClonesArray *) FairRootManager::Instance()->GetObject("ETOFPoint");
  fMCTracks =(TClonesArray *) FairRootManager::Instance()->GetObject("MCTrack");
  //fSTSTrackMatch = (TClonesArray*) FairRootManager::Instance()->GetObject("STSTrackMatch");
  fPrimVtx =  (TClonesArray *) FairRootManager::Instance() ->GetObject("Vertex");

  fNPass = 1;
  return kSUCCESS;
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::Reset() 
{
  ///
  cout << " MpdEctTrackFinderTof::Reset  " << endl;

  //fKHits->Clear();
  //fTracks->Clear("C");
  if (!fTpc) {
    fKHits->Clear();
    fTracks->Delete();
  }
  fTrackCand->Delete();
  delete [] fLayPointers;
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::SetParContainers()
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
void MpdEctTrackFinderTof::Finish()
{
  //Write();
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::Exec(Option_t * option)
{

  static int eventCounter = 0;    
  cout << " EctRec event " << ++eventCounter << endl;

  Reset();

  // Create Kalman hits
  MakeKalmanHits();
  if (fKHits->GetEntriesFast() == 0) return;

  // Match with ETOF hits for tracks extrapolated from TPC 
  if (fTpc) MatchEtof();
    
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
    RemoveDoubles(); // remove double tracks
    StoreTracks(i);
    cout << "  Total number of found tracks: " << fTracks->GetEntriesFast() << endl;
    if (i != fNPass - 1) ExcludeHits(); // exclude used hits
  }
  SelectTracks(0); // do track selection and compute shared hit multiplicities
  AddHits(); // add hit objects to tracks
  cout << "  Total number of found tracks: " << fTracks->GetEntriesFast() << endl;
  if (fTracks->GetEntriesFast() < 1) return;

  if (fTpcTracks) MatchTpc();
  GoToBeamLine();
  Smooth(); // apply constraints due to primary vertex
  GoOutward(); // for matching with outer detectors
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::InitGeo()
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
    FairGeoNode *lay = (*it).second;
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
    FairGeoNode *tube = (*it).second;
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
void MpdEctTrackFinderTof::MakeKalmanHits()
{
  /// Create Kalman hits from ECT points.

  if (!fTpc) {
    // Without TPC tracking
    Bool_t mirror = kTRUE; //kFALSE; //kTRUE; // add mirror hits if TRUE
    if (fExact) mirror = kFALSE;

    Int_t nHits = fEctHits->GetEntriesFast(), lay = 0, nKH = 0, itube = 0;
    //Double_t phi, r, coord, errR = 0.2, errRphi = 0.02; // 2mm in R, 200um in R-Phi
    Double_t phi, coord, errR = 0.02, errRphi = 0.02; //r, 200um in R, 200um in R-Phi
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

      phi = h->GetPhi(); // tube Phi
      lay = h->GetDetectorID() / 1000;
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
      // Transform to the tube local coordinate system
      Double_t cosSin[2] = {TMath::Cos(phi), TMath::Sin(phi)}; // phi - tube direction
      //Double_t xLoc = h->GetX() * cosPhi + h->GetY() * sinPhi; // cross-check
      //Double_t yLoc = -h->GetX() * sinPhi + h->GetY() * cosPhi;
      //Double_t xLoc = xNew * cosSin[0] + yNew * cosSin[1];
      Double_t yLoc = -xNew * cosSin[1] + yNew * cosSin[0];
      //Double_t xLoc = (xNew - h->GetX()) * cosPhi + (yNew - h->GetY()) * sinPhi;
      //Double_t yLoc = -(xNew - h->GetX()) * sinPhi + (yNew - h->GetY())  * cosPhi;
      //r = xNew * xNew + yNew * yNew;
      //r = TMath::Sqrt (r);
      //r = TMath::Abs(xLoc);
      //r = xLoc;
      //cout << xLoc << " " << yLoc << " " << r << " " << h->GetPz() << endl;
      coord = yLoc;
      Double_t yWire = -h->GetX() * cosSin[1] + h->GetY() * cosSin[0];
      Double_t dY = TMath::Abs (yWire - coord);
      
      // Add error                                            
      Double_t dRphi = 0, dR = 0;
      gRandom->Rannor(dRphi,dR); // add errors
      Double_t meas[2] = {coord+dRphi*errRphi, 0};
      Double_t err[2] = {errRphi, errR};
      MpdKalmanHit *hit = new ((*fKHits)[nKH]) MpdKalmanHit(lay*1000000, 1, MpdKalmanHit::kFixedZ, 
			      //MpdKalmanHit *hit = new ((*fKHits)[nKH]) MpdKalmanHit(h->GetDetectorID()*1000, 1, MpdKalmanHit::kFixedZ, 
							    meas, err, cosSin, 9., h->GetZ(), ih);
      //hit->SetXY(h->GetX(), h->GetY()); // wire coordinate
      // Add second measurement - just for test at the moment
      //!!!
      //hit->SetNofDim(2);
      //!!!
      //lay = hit.GetLayer();
      //layMax = TMath::Max (lay, layMax);
      //fhLays->Fill(lay+0.1);
      
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
	  if (mirror) fKHits->RemoveAt(indxTube[lay][itube]+1);
	  indxTube[lay][itube] = nKH;
	} else {
	  hitOld->SetIndex(hit->GetIndex());
	  fKHits->RemoveAt(nKH); // remove last hit
	}
      } else {
	firstHit[lay][itube] = dY;
	indxTube[lay][itube] = nKH;
	//fhLays->Fill(lay+0.1);
      }
      if (mirror && fKHits->UncheckedAt(nKH)) {
	// Add mirror hit
	// wire point transverse coordinate in rotated system:
	//yLoc = -h->GetX() * cosSin[1] + h->GetY() * cosSin[0];
	//Double_t dy = hit->GetRphi() - yRot;
	Double_t dy = hit->GetMeas(0) - yWire;
	MpdKalmanHit *hitM = new ((*fKHits)[++nKH]) MpdKalmanHit(*hit);
	//hitM->SetRphi(yRot-dy);
	hitM->SetMeas(0, yWire-dy);
	hitM->SetMirror();
	//fhLays->Fill(lay+0.1);
      }
      ++nKH;
    }
    fKHits->Compress();
    fKHits->Sort(); // in descending order in abs(Z)
  } // if (!fTpc)

  if (fKHits->GetEntriesFast() == 0) return;
  Int_t layMax = ((MpdKalmanHit*)fKHits->First())->GetLayer();
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
  Int_t nKHits = fKHits->GetEntriesFast(), nOK = 0;
  for (Int_t i = 0; i < nKHits; ++i) {
    MpdKalmanHit *hit = (MpdKalmanHit*) fKHits->UncheckedAt(i);
    Int_t lay = hit->GetLayer();
    if (hit->GetFlag() > -1) { hit->SetFlag(1+hit->IsMirror()*1000); ++nOK; }
    fhLays->Fill(lay+0.1);
    if (fLayPointers[lay] < 0) fLayPointers[lay] = i;
  }
  cout << " Max layer = " << layMax << " " << fKHits->GetEntriesFast() << " " << nOK << endl;
  //exit(0);
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::MatchEtof()
{
  /// Match tracks, extended from TPC, with ETOF hits.

  static const Double_t rMax = 2.;

  MpdKalmanHit hTmp;
  hTmp.SetType(MpdKalmanHit::kFixedZ);
  hTmp.SetDist(((MpdTofHit*) fTofHits->First())->GetZ());
  Int_t nTof = fTofHits->GetEntriesFast();

  Int_t nTpc = fTracks->GetEntriesFast(), nMatch = 0;
  for (Int_t i = 0; i < nTpc; ++i) {
    MpdEctKalmanTrack *tr = (MpdEctKalmanTrack*) fTracks->UncheckedAt(i);
    MpdEctKalmanTrack track = *tr;
    //MpdKalmanFilter::Instance()->PropagateToHit(&track,&hTmp);
    track.SetParam(*track.GetParamNew());
    track.SetPos(track.GetPosNew());
    hTmp.SetDist(TMath::Abs(hTmp.GetDist())*TMath::Sign(1.,track.GetPos()));
    //MpdKalmanFilter::Instance()->PropagateParamZ(&track,&hTmp,kFALSE);
    MpdKalmanFilter::Instance()->PropagateParamZ(&track,&hTmp,kTRUE);
    Int_t iMin = 0; //id = track.GetTrackID(),
    Double_t rMin = 999999., xTr = track.GetParamNew(0), yTr = track.GetParamNew(1); //rID = 999,

    for (Int_t itof = 0; itof < nTof; ++itof) {
      MpdTofHit *tof = (MpdTofHit*) fTofHits->UncheckedAt(itof);
      //MpdTofPoint *tofP = (MpdTofPoint*) fTofPoints->UncheckedAt(tof->GetRefIndex());

      //if (tof->GetZ() < 0) continue; // !!! exclude one side for now

      Double_t dx = tof->GetX() - xTr, dy = tof->GetY() - yTr;
      Double_t r = dx * dx + dy * dy;
      //if (id == tofP->GetTrackID()) rID = r;
      if (r < rMin) { rMin = r; iMin = itof; }
    }
    //if (lun2) fprintf(lun2,"%d %10.3e %10.3e\n", id, TMath::Sqrt(rID), TMath::Sqrt(rMin));
    if (TMath::Sqrt(rMin) < rMax) {
      // Exclude ETOF hit
      //if (((MpdTofHit*) fTofHits->UncheckedAt(iMin))->GetFlag() != -1) {
	++nMatch;
	MpdTofHit *tof = (MpdTofHit*) fTofHits->UncheckedAt(iMin);
	tof->SetFlag(tof->GetFlag() | MpdTofUtils::IsSelected);
	//3-may-12 tr->SetTofIndex(iMin);
	tr->SetMatchEtof();
	//}
    }
    // Modify track length
    tr->SetLength(track.GetLength());
  }
  cout << " *** Found " << nMatch << " matches with ETOF hits *** " << endl;
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::GetTrackSeeds(Int_t iPass)
{
  /// Build track seeds from TOF, ECT hits and vertex

  Int_t tmp[99999] = {0}, idMax = 0, nID = 0;
  Int_t nTof = fTofHits->GetEntriesFast();
 
  Int_t nCand = 0, layMax = ((MpdKalmanHit*)fKHits->First())->GetLayer();
  const Double_t dphiMax = 0.1; //0.2; // cut on Phi between ECT and ETOF
  TVector3 vert(0.,0.,0.), pmom;

  cout << " Number of ETOF hits: " << nTof << endl;
  TVector3 posTof;
  for (Int_t itof = 0; itof < nTof; ++itof) {
    MpdTofHit *tof = (MpdTofHit*) fTofHits->UncheckedAt(itof);
    if (tof->GetFlag() & MpdTofUtils::IsSelected) continue; 
    MpdTofPoint *tofP = (MpdTofPoint*) fTofPoints->UncheckedAt(tof->GetRefIndex());

    //if (tof->GetZ() < 0) continue; // !!! exclude one side for now

    tof->Position(posTof);
    Double_t rTof = posTof.Pt();
    Double_t drdz = rTof / tof->GetZ();
    Double_t phTof = posTof.Phi(), dph = 0, phEct = 0;
    tofP->Momentum(pmom);

    Int_t layBeg = layMax;
    if (tof->GetZ() > 0 && layMax > fgkNlays2) layBeg = layMax - fgkNlays2;
    Int_t layEnd = layBeg - 10;
    //for (Int_t lay = layMax; lay > 50; --lay) {
    for (Int_t lay = layBeg; lay > layEnd; --lay) {
      Int_t nLay = GetNofHitsInLayer(lay);
      Int_t indx0 = GetHitsInLayer(lay);
      Int_t view = (lay - 1) % 3; // view

      for (Int_t indx = 0; indx < nLay; ++indx) {
	MpdKalmanHit *hit = (MpdKalmanHit*) fKHits->UncheckedAt(indx0+indx);
	// Exclude used hits
	if (hit->GetFlag() < 0) continue;
	FairMCPoint *ectP =  (FairMCPoint*) fEctHits->UncheckedAt(hit->GetIndex());
        Int_t idEct = ectP->GetTrackID();

	// !!! Exact match of TOF and ECT
	//if (fExact) {
	  //Int_t itmp = 0;
	  //if (fExact && idEct != tofP->GetTrackID()) continue;
	  //if (HitMotherId(hit,tofP->GetTrackID(),itmp) != tofP->GetTrackID()) continue;
	//}

	Int_t detID = ectP->GetDetectorID();
	if (fWireMap.find(detID) == fWireMap.end()) {
	  // Add position of wire point
	  fWireMap.insert(pair<Int_t,pair<Double_t,Double_t> >(detID,
							       pair<Double_t,Double_t>(ectP->GetX(),ectP->GetY())));
	}

	// Check Phi between ECT and ETOF
	Double_t rEct = drdz * hit->GetPos();
	if (view == 0) {
	  // radial view
	  phEct = hit->GetPhi();
	} else {
	  // stereo view - find PHI of wire point at radius, extrapolated from ETOF 
	  //Double_t b = hit->GetXY(0) * hit->GetCos() + hit->GetXY(1) * hit->GetSin();
	  //Double_t d = hit->GetXY(0) * hit->GetXY(0) + hit->GetXY(1) * hit->GetXY(1) 
	  Double_t b = fWireMap[detID].first * hit->GetCosSin(0) + fWireMap[detID].second * hit->GetCosSin(1);
	  Double_t d = fWireMap[detID].first * fWireMap[detID].first + fWireMap[detID].second * fWireMap[detID].second 
	             - rEct * rEct;
	  Double_t l1 = -b + TMath::Sqrt (b * b - d);
	  Double_t l2 = -b - TMath::Sqrt (b * b - d);
	  if (TMath::Abs(l1) > TMath::Abs(l2)) l1 = l2;
	  //phEct = TMath::ATan2 (hit->GetXY(1) + l1 * hit->GetSin(),
	  //		hit->GetXY(0) + l1 * hit->GetCos());
	  phEct = TMath::ATan2 (fWireMap[detID].second + l1 * hit->GetCosSin(1),
				fWireMap[detID].first + l1 * hit->GetCosSin(0));
	}
	dph = MpdKalmanFilter::Instance()->Proxim(phTof, phEct) - phTof;
	if (TMath::Abs(dph) > dphiMax) continue;
	// Add track seed
	//if (lun1 && idTpc == tof->GetTrackID()) fprintf(lun1,"%10.3e %10.3e %10.3e %10.3e\n",rTpc,rTof,dphi,slope-slope0);
	if (tofP->GetTrackID() < 99999) { tmp[tofP->GetTrackID()]++; idMax = TMath::Max(idMax,tofP->GetTrackID()); }
	MpdEctKalmanTrack *track = 
	  //new ((*fTracks)[nCand++]) MpdEctKalmanTrack(itof, indx0+indx, tof, hit, vert);
	  new ((*fTrackCand)[nCand++]) MpdEctKalmanTrack(itof, idEct, tof, hit, vert);
	track->SetTrackID(idEct);
	EvalParams(tof, hit, track, rEct, phEct);

	// Eval. track Phi - linear extrapolation of 2 estimates
	//track->SetParam(2, phTof); // phi
	Double_t dx = tof->GetX() - rEct * TMath::Cos(phEct);
	Double_t dy = tof->GetY() - rEct * TMath::Sin(phEct);
	Double_t dr = TMath::Sqrt (dx * dx + dy * dy);
	Double_t phTofEct = TMath::ATan2 (dy, dx);
	Double_t ph = phTofEct;
	//Double_t ph = phEct - 1.0 * (MpdKalmanFilter::Instance()->Proxim(phTofEct,phTof)-phTofEct)/(rTof-dr) * rEct;
	track->SetParam(2, ph); // phi
	Double_t th = TMath::ATan2 (tof->GetZ()-hit->GetPos(), dr);
	track->SetParam(3, th); // dip angle
	Double_t pt = 1. / TMath::Abs(track->GetParam(4));
	/*
	track->SetParam(2, pmom.Phi()); // exact value - for test
	track->SetParam(3, TMath::PiOver2()-pmom.Theta()); // exact value - for test
	track->SetParam(4, TMath::Sign(1.,track->GetParam(4))/pmom.Pt()); // exact value - for test
	*/
	/*
	if (track->GetParam(3) < 0) {
	  // Going in backward direction
	  --nCand;
	  fTrackCand->RemoveAt(nCand);
	  continue;
	}
	*/
	track->SetParamNew(*track->GetParam());
	//track->GetParam()->Print();
	EvalCovar(tof, hit, track, rEct, phEct);
	//if (lun1 && idEct == tofP->GetTrackID()) fprintf(lun1,"%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",track->GetParam(2),MpdKalmanFilter::Instance()->Proxim(track->GetParam(2),ph),track->GetParam(3),MpdKalmanFilter::Instance()->Proxim(track->GetParam(3),th),1./TMath::Sqrt((*track->GetWeight())(2,2)),1./TMath::Sqrt((*track->GetWeight())(3,3)),pmom.Mag(),rTof,pmom.Pt(),pt);
	if (lun1 && idEct == tofP->GetTrackID()) fprintf(lun1,"%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",track->GetParam(2),MpdKalmanFilter::Instance()->Proxim(track->GetParam(2),pmom.Phi()),track->GetParam(3),MpdKalmanFilter::Instance()->Proxim(track->GetParam(3),TMath::PiOver2()-pmom.Theta()),1./TMath::Sqrt((*track->GetWeight())(2,2)),1./TMath::Sqrt((*track->GetWeight())(3,3)),pmom.Mag(),rTof,pmom.Pt(),pt);

	//cout << nCand-1 << " " << lay << " " << track->GetTrackID() << " " << tofP->GetTrackID() << endl;
      } // for (Int_t indx = 0; indx < nLay;
    } // for (Int_t lay = layMax; lay > 50;
  } // for (Int_t itof = 0; itof < nTof;

  for (Int_t j = 0; j <= idMax; ++j) if (tmp[j] > 0) ++nID;
  cout << " Number of ECT track candidates: " << nCand << " " << nID << endl;
  //if (lun1) fclose(lun1);
  //exit(0);
}  
  
//__________________________________________________________________________
void MpdEctTrackFinderTof::EvalParams(const MpdTofHit *tof, const MpdKalmanHit *ect, 
				      MpdEctKalmanTrack *track, Double_t rEct, Double_t phEct) 
{
  /// Evaluate track parameters

  // Evaluate signed track Pt (curvature) assuming the track coming from the 
  // primary vertex
  //Double_t rTof = TMath::Sqrt (tof->GetX() * tof->GetX() + tof->GetY() * tof->GetY());
  Double_t phTof = TMath::ATan2 (tof->GetY(), tof->GetX());
  //TVector2 vecTof(rTof*TMath::Cos(phTof)-0.,rTof*TMath::Sin(phTof)-0.);
  TVector2 vecTof(tof->GetX()-0.,tof->GetY()-0.);
  TVector2 vecEct(rEct*TMath::Cos(phEct)-0.,rEct*TMath::Sin(phEct)-0.);
  TVector2 dvec = vecTof - vecEct;
  Double_t cosAlpha = vecEct * dvec / vecEct.Mod() / dvec.Mod();
  Double_t rad = vecTof.Mod() / 2. / TMath::Sin(TMath::ACos(cosAlpha));
  Double_t factor = 0.0015; // 0.3 * 0.5T * 0.01
  //Double_t charge = phTof - MpdKalmanFilter::Instance()->Proxim(phTof,phTpc);
  Double_t charge = MpdKalmanFilter::Instance()->Proxim(phTof,phEct) - phTof;
  Double_t pt = factor * TMath::Abs(rad);
  pt = TMath::Min (pt, 2.5) * TMath::Sign(1., -charge);

  track->SetParam(0, vecTof.X()); // X
  track->SetParam(1, vecTof.Y()); // Y
  track->SetParam(4, 1./pt); // 1/pt
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::EvalCovar(const MpdTofHit *tof, const MpdKalmanHit *ect, 
				     MpdEctKalmanTrack *track, Double_t rEct, Double_t phEct)
{
  /// Evaluate covariance matrix for the track seed

  static const Double_t xErrTof = 2. / TMath::Sqrt(12.),  yErrTof = xErrTof;
  //static const Double_t xErrTof = 0.1, yErrTof = xErrTof;
  TMatrixDSym weight(5);
  weight(0,0) = xErrTof * xErrTof; // <xx>
  //weight(0,0) *= 4.; // extra factor of 2

  weight(1,1) = yErrTof * yErrTof; // <yy>
  //weight(1,1) *= 4.; // extra factor of 2

  Double_t xEct = rEct * TMath::Cos(phEct);
  Double_t yEct = rEct * TMath::Sin(phEct);
  Double_t xTof = tof->GetX(), yTof = tof->GetY();
  Double_t dx = xTof - xEct, dy = yTof - yEct;
  Double_t dist2 = dx * dx + dy * dy;
  weight(2,2) = xErrTof * xErrTof / dist2; // <PhiPhi>
  //weight(2,2) *= 2.; // extra factor of sqrt(2)

  //Double_t tanThe = TMath::Tan (track->GetParam(3));
  Double_t rTof = TMath::Sqrt (tof->GetX() * tof->GetX() + tof->GetY() * tof->GetY());
  //Double_t dR = rTof - rEct;
  //Double_t dz = tof->GetZ() - ect->GetPos();
  Double_t dR = rTof;
  Double_t dz = tof->GetZ();
  Double_t dThdR = dz / (dR * dR + dz * dz);
  weight(3,3) = dThdR * dThdR * xErrTof * xErrTof; // <TheThe>
  //weight(3,3) /= 4.; // extra factor of 0.5
  weight(3,3) *= 4.; // extra factor of 2

  //weight(4,4) = (track->GetParam(4)*0.5) * (track->GetParam(4)*0.5); // error 50%
  weight(4,4) = (track->GetParam(4)*1.) * (track->GetParam(4)*1.); // error 100%
  //weight(4,4) = TMath::Max (1.5 * 1.5, weight(4,4));
  //weight(4,4) = (track->GetParam(4)*2.) * (track->GetParam(4)*2.); // error 200%
  //weight(4,4) = (track->GetParam(4)*4.) * (track->GetParam(4)*4.); // error 400%

  //weight.Print();
  //fWeight->Invert(); // weight matrix
  Int_t iok = 0;
  MpdKalmanFilter::Instance()->MnvertLocal(weight.GetMatrixArray(), 5, 5, 5, iok);
  //fWeight->Print();
  track->SetWeight(weight);
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::DoTracking(Int_t iPass)
{
  /// Run Kalman tracking
  
  Int_t nCand = fTrackCand->GetEntriesFast(), iok = -1;
  TMatrixD param(5,1);
  TMatrixDSym weight(5), pointWeight(5);
 
  for (Int_t i = 0; i < nCand; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTrackCand->UncheckedAt(i);
    //cout << " Track seed No. " << i << " " << track->GetTrackID() << " " ;//<< endl;
    //if (track->GetTrackID() != 117) continue;
    //if (track->GetTrackID() != 10) continue;
    //if (track->GetTrackID() != 1105) continue;
    //if (track->GetTrackID() != 77) continue;
    //(*(track->GetParamNew()))(4,0) = -0.5; // just for check
    /*
    for (Int_t k = 0; k < 5; ++k) {
      for (Int_t j = i; j < 5; ++j) {
	if (j == i) (*track->GetWeight())(i,j) /= 100.;
	else (*track->GetWeight())(i,j) = (*track->GetWeight())(j,i) = 0.;
      }
    }
    */

    // Add first ECT hit
    //track->Weight2Cov()->Print();
    //track->GetParamNew()->Print();
    MpdKalmanHit *hit = (MpdKalmanHit*) track->GetHits()->UncheckedAt(0);
    MpdKalmanFilter::Instance()->PropagateToHit(track, hit);
    //cout << track->GetLength() << endl;
    Double_t dChi2 = MpdKalmanFilter::Instance()->FilterHitZ(track,hit,pointWeight,param);
    track->SetChi2(track->GetChi2()+dChi2);
    weight = *track->GetWeight();
    weight += pointWeight;
    track->SetWeight(weight);
    track->SetParamNew(param);
    //if (i == 0) track->SetLength(0.); // for correct track length
    //track->Weight2Cov()->Print();
    //track->GetParamNew()->Print();

    // Start ECT tracking from different layers to account for missing hits
    Int_t laySkip = 2, lay0 = hit->GetLayer() - 1;
    //Int_t layMax = ((MpdKalmanHit*)fKHits->Last())->GetLayer();
    MpdEctKalmanTrack tracks[10];

    for (Int_t j = 0; j < laySkip; ++j) {
      tracks[j] = *track;
      //cout << track->GetParamNew(0) << endl;
      //cout << i << " " << lay << " " << tracks[lay].GetNofHits() << " " << tracks[lay].GetChi2() << " " << tracks[lay].GetParam(0) << endl;
      if (j > 0) {
	lay0 = -1;
	// Find first layer after the missing one
	Int_t nh = tracks[j-1].GetNofHits(), lastLay = hit->GetLayer();
	TObjArray *hits = tracks[j-1].GetHits();
	for (Int_t jh = 1; jh < nh; ++jh) {
	  hit = (MpdKalmanHit*) hits->UncheckedAt(jh);
	  if (TMath::Abs(hit->GetLayer()-lastLay) == 1) { 
	    lastLay = hit->GetLayer();
	    tracks[j].GetHits()->Add(hit);
	    continue; 
	  }
	  lay0 = hit->GetLayer() - 1;
	  break;
	}
	if (lay0 >= 0) {
	  // Refit track to the last kept hit
	  nh = tracks[j].GetNofHits();
	  for (Int_t jh = 1; jh < nh; ++jh) {
	    hit = (MpdKalmanHit*) tracks[j].GetHits()->UncheckedAt(jh);
	    MpdKalmanFilter::Instance()->PropagateToHit(&tracks[j], hit);
	    dChi2 = MpdKalmanFilter::Instance()->FilterHitZ(&tracks[j],hit,pointWeight,param);
	    tracks[j].SetChi2(tracks[j].GetChi2()+dChi2);
	    weight = *tracks[j].GetWeight();
	    weight += pointWeight;
	    tracks[j].SetWeight(weight);
	    tracks[j].SetParamNew(param);
	    if (jh == nh - 1) {
	      tracks[j].SetWeightAtHit(weight);
	      tracks[j].SetParamAtHit(param);
	      tracks[j].SetLengAtHit(tracks[j].GetLength());
	    }
	  }
	}
      } // if (j > 0)
      //if (lay0 >= 0) iok = RunKalmanFilter(&tracks[j], lay0);
      if (lay0 >= 0) iok &= RunKalmanFilter(&tracks[j], lay0);
      else --laySkip;
      //cout << iok << " " << tracks[j].GetMisses() << " ";//<< endl;
    }
    //cout << endl;

    // Select the best track (with max number of hits)
    Int_t nHitsMax = tracks[0].GetNofHits(), iMax = 0;
    //if (iok > -9) {
    if (iok >= 0) {
      for (Int_t j = 1; j < laySkip; ++j) {
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
    }

    fTrackCand->RemoveAt(i);
    if (iok < 0) continue;
    track = new ((*fTrackCand)[i]) MpdEctKalmanTrack(tracks[iMax]);

    // Refit
    if (1) {
      //track->Weight2Cov()->Print();
      //track->GetParamNew()->Print();
      Double_t pos = track->GetPosNew();
      track->SetDirection(MpdKalmanTrack::kOutward);
      MpdKalmanFilter::Instance()->Refit(track, -1);
      //track->GetParamNew()->Print();
      track->SetDirection(MpdKalmanTrack::kInward);
      MpdKalmanFilter::Instance()->Refit(track, 1);
      track->SetParamAtHit(*track->GetParamNew());
      track->SetWeightAtHit(*track->GetWeight());
      //track->SetLengAtHit(track->GetLength());
      MpdKalmanHit hitTmp;
      hitTmp.SetType(MpdKalmanHit::kFixedZ);
      hitTmp.SetPos(pos);
      MpdKalmanFilter::Instance()->PropagateToHit(track,&hitTmp,kFALSE);
    }
    //iok = RunKalmanFilter(track, hit->GetLayer()-1);

    //cout << i << " " << track->GetNofHits() << endl;
    if (0 && iok == 0 && track->GetNofHits() < 0) {
      //*
      hit = (MpdKalmanHit*) track->GetHits()->UncheckedAt(0);
      Double_t dir = TMath::Sign(1.,track->GetPosNew()); 
      track->SetPos(hit->GetPos()-dir);
      track->SetPosNew(track->GetPos());
      //for (Int_t j = 0; j < 5; ++j) track->SetParam(j,track->GetParam1()[j]);
      //track->SetParamNew(*track->GetParam());
      //MpdKalmanHitZ hitTmp(*hit);
      //hitTmp.SetZ(track->GetPosNew()-dir);
      //MpdKalmanFilter::Instance()->PropagateToHit(track,&hitTmp);
      MpdKalmanFilter::Instance()->Refit(track,1);
      //*/
      /*
      track->Print("");
      Double_t dir = TMath::Sign(1.,track->GetPosNew()); 
      MpdKalmanHitZ *hit = (MpdKalmanHitZ*) track->GetHits()->UncheckedAt(0);
      MpdKalmanHitZ hitTmp(*hit);
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
  }
}
    
//__________________________________________________________________________
Int_t MpdEctTrackFinderTof::RunKalmanFilter(MpdEctKalmanTrack *track, Int_t layBeg)
{
  /// Run Kalman filter

  const Double_t rMin = 49.2, rMax = 110.5; // min and max radial size of ECT - to be taken elsewhere
  //cout << fHits->GetEntriesFast() << endl;
  Int_t layMax = ((MpdKalmanHit*)fKHits->Last())->GetLayer();
  MpdKalmanHit *hitOK = 0x0;
  MpdKalmanTrack::TrackDir trackDir = track->GetDirection();
  //Int_t layBeg = layMax, layEnd = -1, dLay = -1, layOK = -1;
  Int_t layEnd = -1, dLay = -1;//, layOK = -1;
  if (track->GetPosNew() < 0) layEnd = fgkNlays2;
  if (trackDir == MpdKalmanTrack::kOutward) {
    layEnd = layMax + 1;
    dLay = 1;
    if (track->GetPosNew() > 0 && layMax > fgkNlays2) layEnd -= fgkNlays2;
  }
  
  //Int_t indxOK = hits->IndexOf(hitOK);
  //Int_t nHits = hits->GetEntriesFast();
  Int_t miss = 0, missSec = 0;
  TMatrixDSym pointWeight(5), pointWeightTmp(5), saveWeight(5);
  TMatrixD param(5,1), paramTmp(5,1);

  Double_t saveZ = 0;
  //cout << " Starting hit: " << hitOK->GetLayer() << " " << hitOK->GetTrackID() << " " << hitOK->GetUsage() << endl;
  
  for (Int_t lay = layBeg; lay != layEnd; lay+=dLay) {
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
    Double_t radNew = -1;
    //for (Int_t indx = 0; indx < nLay; ++indx) {
    for (Int_t indx = indxBeg; indx != indxEnd; indx+=dIndx) {
      MpdKalmanHit *hit = (MpdKalmanHit*) fKHits->UncheckedAt(indx0+indx);

      // !!! Exact ID match
      Int_t itmp = 0;
      //if (fExact && ((FairMCPoint*)fEctHits->UncheckedAt(hit->GetIndex()))->GetTrackID() != track->GetTrackID()) continue;
      if (fExact && HitMotherId(hit,track->GetTrackID(),itmp) != track->GetTrackID()) continue;

      // Exclude used hits
      if (hit->GetFlag() < 0) continue;
      //cout << " window:" << /*hit->GetTrackID() <<*/ " " << hit->GetRphi() << " " << track->GetParamNew(0) << " " << hit->GetZ() << " " << track->GetParamNew(1) << endl;
      // Check if the hit within some window (15x15cm for the moment - check!!!)
      //if (TMath::Abs(hit->GetRphi()-track->GetParamNew(0)) > 9) continue;
      //if (TMath::Abs(Proxim(track,hit)-track->GetParamNew(0)) > 15) continue;
      TVector2 dist = GetDistance(track, hit);
      if (dist.X() > 15.) continue; // distance in transverse to the tube direction
      if (hit->GetNofDim() > 1 && dist.Y() > 25.) continue; // distance in R 
      //*if (hit->GetTrackID() == 186)*/ cout << " Hit: " << ((FairMCPoint*)fEctHits->UncheckedAt(hit->GetIndex()))->GetTrackID() << " " << hit->GetLayer() << " " << ((FairMCPoint*)fEctHits->UncheckedAt(hit->GetIndex()))->GetDetectorID() << " " << hit->GetRphi() << " " << hit->GetR() << " " << hit->GetZ() << " " << dist.X() << " " << dist.Y() << " " << track->GetParamNew(0) << " " << track->GetParamNew(1) << endl;
      //track->Print("");
      //hit->Print("");
      if (!MpdKalmanFilter::Instance()->PropagateToHit(track,hit)) { track->SetMisses(miss); return -1; }

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
      //Double_t radNew = track->GetRadNew();
      if (radNew < 0) radNew = track->GetRadNew();
      //if (radNew < rMin || radNew > rMax) return -1;
      // Going to interaction point
      if (radNew > rMax) break; // next layer
      //else if (radNew < rMin) return -1;
      else if (radNew < rMin) { track->SetMisses(miss); return 0; }

      //Double_t err = hit->GetRphiErr();
      //if (track->GetNofHits() == 0) hit->SetRphiErr(0.04);

      Double_t dChi2 = MpdKalmanFilter::Instance()->FilterHitZ(track,hit,pointWeight,param);
      //if (track->GetNofHits() == 0) hit->SetRphiErr(err);
      //if (param(3,0) < 0) { cout << " Check: " << param(3,0) << " " << dChi2 << " " << (*fParamNew)(3,0) << " " << hit->GetRphi() << " " << hit->GetZ() << endl; fParamNew->Print();}
      if (dChi2 < dChi2Min) {
      //if (dChi2 < dChi2Min && dist0.X() <= dist.X()) {
      //if (dChi2 < dChi2Min && dist.X() <= 0.2) {
        dChi2Min = dChi2;
        hitMin = hit;
        saveWeight = *track->GetWeight();
        saveZ = track->GetPosNew();
        // temporary storage for the current track
        paramTmp = param;
        pointWeightTmp = pointWeight;
        //cout << " New min dChi2 = " << dChi2 << " " << hitMin->GetRphi() << " " << hitMin->GetR() << " " << ((FairMCPoint*)fEctHits->UncheckedAt(hit->GetIndex()))->GetTrackID() << endl;
	//cout << track->GetParamNew(0) << " " << track->GetParamNew(1) << endl;
	//cout << hit->GetRphi() << " " << hit->GetZ() << endl;
	//cout << param(0,0) << " " << param(1,0) << endl;
        //paramTmp.Print();
      }
    } // for (Int_t indx = indxBeg; indx != indxEnd;
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
      track->SetParamNew(paramTmp);
      //cout << " *** Adding hit: " << hitOK->GetLayer() << " " << ((FairMCPoint*)fEctHits->UncheckedAt(hitOK->GetIndex()))->GetTrackID() << " " << ((FairMCPoint*)fEctHits->UncheckedAt(hitOK->GetIndex()))->GetDetectorID() << " " << dChi2Min << " " << track->GetChi2() << endl;
      //paramTmp.Print();
      // Check if the accepted hit is the same as the seeded hit
      //if (hitOK->GetLayer() == f2ndHit->GetLayer() && hitOK != f2ndHit) return -1; // abandon track
      if (track->GetNofHits() == 1) track->SetParam1();
      if (track->GetDirection() == MpdKalmanTrack::kInward) {
	// Save track params at last hit
	track->SetLengAtHit(track->GetLength()); 
	track->SetParamAtHit(paramTmp);
	track->SetWeightAtHit(saveWeight);
	track->SetPosAtHit(saveZ);
      }
      missSec = 0;
    } else if (lay > 0 && radNew < rMax && radNew > rMin) {
      // Check whether tracks goes through a tube
      Int_t layM1 = lay - 1;
      if (nLay == 0) {
	// No hits in the layer - extrapolate track
	MpdKalmanHit hitTmp;
	hitTmp.SetType(MpdKalmanHit::kFixedZ);
	hitTmp.SetDist(fZplanes[layM1]);
	if (!MpdKalmanFilter::Instance()->PropagateToHit(track,&hitTmp)) { track->SetMisses(miss); return -1; }
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
      if (dyMin < frMinMax[2]-0.2) { ++miss; ++missSec; } // going thru the tube
      //cout << dyMin << " " << dChi2Min << endl;
    }
    //cout << " lay, miss: " << lay << " " << miss << " " << dChi2Min << " " << fChi2 << endl;
    if (missSec > 5) { track->SetMisses(miss); return -9; }
  } // for (Int_t lay = layBeg; lay != layEnd;
  track->SetMisses(miss);
  return 0;
}

//__________________________________________________________________________
TVector2 MpdEctTrackFinderTof::GetDistance(MpdEctKalmanTrack *track, MpdKalmanHit *hit)
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
  return dist;
}

//__________________________________________________________________________
Double_t MpdEctTrackFinderTof::Proxim(MpdKalmanTrack *track, const MpdKalmanHit *hit)
{
  /// Adjust hit coord. R-Phi to be "around" track R-Phi - to avoid 
  /// discontinuity around +- Pi

  Fatal("Proxim", " !!! Not implemented !!!");
  /*
  Double_t hitPhi = hit->GetRphi() / hit->GetR();
  Double_t phi0 = track->GetParamNew(0) / track->GetPosNew();
  return hit->GetR() * MpdKalmanFilter::Instance()->Proxim(phi0, hitPhi);
  */
  return 0;
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::Write()
{
  /// Write

  TFile histoFile("EctRec.root","RECREATE");
  Writedir2current(fHistoDir);
  histoFile.Close();
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::Writedir2current( TObject *obj )
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
void MpdEctTrackFinderTof::RemoveDoubles()
{
  /// Remove double tracks (if number of common hits greater than 50% of hits on track)

  Int_t nFound = fTrackCand->GetEntriesFast(), nOK = 0;
  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTrackCand->UncheckedAt(i);
    if (track == 0x0) continue;
    Int_t nh = track->GetNofHits();
    //cout << i << " " << nh << " " << ++nOK << endl;
    for (Int_t j = i + 1; j < nFound; ++j) {
      MpdEctKalmanTrack *track1 = (MpdEctKalmanTrack*) fTrackCand->UncheckedAt(j);
      if (track1 == 0x0) continue;
      Int_t nh1 = track1->GetNofHits();
      Int_t nc = NofCommonHits(track, track1);
      if (1.*nc/TMath::Min(nh,nh1) < 0.5) continue;
      if (nh > nh1) fTrackCand->RemoveAt(j);
      //if (nh > nh1+5) fTrackCand->RemoveAt(j);
      else if (nh < nh1) {
      //else if (nh < nh1-5) {
	fTrackCand->RemoveAt(i);
	--nOK;
	break;
      } else {
	if (track->GetChi2() > track1->GetChi2()) {
        //if (track->GetChi2()/nh > track1->GetChi2()/nh1) {
	  fTrackCand->RemoveAt(i);
	  --nOK;
	  break;
	}
	fTrackCand->RemoveAt(j);
      }
    }
  }
  fTrackCand->Compress();

  // Remove double tracks (originated from the same ETOF hit)
  //*
  nFound = fTrackCand->GetEntriesFast();
  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTrackCand->UncheckedAt(i);
    if (track == 0x0) continue;
    Int_t iTof = track->GetTofIndex();
    Int_t nh = track->GetNofHits();
    for (Int_t j = i + 1; j < nFound; ++j) {
      MpdEctKalmanTrack *track1 = (MpdEctKalmanTrack*) fTrackCand->UncheckedAt(j);
      if (track1 == 0x0) continue;
      Int_t iTof1 = track1->GetTofIndex();
      if (iTof1 != iTof) continue;
      Int_t nh1 = track1->GetNofHits();
      if (nh > nh1) fTrackCand->RemoveAt(j);
      else if (nh < nh1) {
	fTrackCand->RemoveAt(i);
	break;
      } else {
	if (track->GetChi2() > track1->GetChi2()) {
	  fTrackCand->RemoveAt(i);
	  break;
	}
	fTrackCand->RemoveAt(j);
      }
    }
  }
  fTrackCand->Compress();
  //*/
}
	  
//__________________________________________________________________________
Int_t MpdEctTrackFinderTof::NofCommonHits(MpdEctKalmanTrack* track, MpdEctKalmanTrack* track1)
{
  /// Compute number of common hits on 2 tracks

  TObjArray *hits = track->GetHits(), *hits1 = track1->GetHits();
  
  MpdKalmanHit *hit = (MpdKalmanHit*) hits->First();
  MpdKalmanHit *hit1 = (MpdKalmanHit*) hits1->First();
  Int_t dir = 1;
  if (hit->GetLayer() < ((MpdKalmanHit*) hits->Last())->GetLayer()) dir = -1;

  Int_t nCom = 0;
  while (hit && hit1) {
    if (dir * hit->GetLayer() < dir * hit1->GetLayer()) {
      hit1 = (MpdKalmanHit*) hits1->After(hit1);
      continue;
    }
    if (hit == hit1) ++nCom;
    hit = (MpdKalmanHit*) hits->After(hit);
  }
  return nCom;
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::StoreTracks(Int_t iPass)
{
  /// Transfer tracks from fTrackCand to fTracks

  static const Int_t nHitMin = 10;
  //static const Double_t etaMin = 1.3, etaMax = 1.5;
  static const Double_t etaMin = 1.6, etaMax = 1.5;

  Int_t nFound = fTracks->GetEntriesFast(), nCand = fTrackCand->GetEntriesFast();
  for (Int_t i = 0; i < nCand; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTrackCand->UncheckedAt(i);
    Double_t eta = TMath::Abs(track->Momentum3().Eta());
    //if (iPass == 0 && track->GetChi2() / (track->GetNofHits()-5) > 3.0) continue;
    //if (iPass == 0 && eta < etaMin) continue;
    if (eta < etaMin) continue;
    if (eta > etaMax && track->GetNofHits() < nHitMin) continue;
    track->Weight2Cov();
    //SetTrackID(track);
    new ((*fTracks)[nFound++]) MpdEctKalmanTrack(*track);
    //fTpcTracks->RemoveAt(track->GetTpcIndex());
    //fTrackCand->RemoveAt(i);
  }
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::AddHits()
{
  /// Add hit objects to tracks and compute number of wrongly assigned hits
  /// (hits with ID different from ID of starting track)

  Int_t nFound = fTracks->GetEntriesFast();
  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTracks->UncheckedAt(i);
    if (track->GetTofIndex() < 0) continue; // track from TPC
    Int_t nHits = track->GetNofHits();
    if (nHits == 0) { fTracks->RemoveAt(i); continue; }
    track->SetNofHits(nHits);
    TClonesArray &trHits = *track->GetTrHits();
    SetTrackID(track); // set track ID as ID of majority of its hits
    TObjArray *hits = track->GetHits();
    Int_t nWrong = 0, motherID = track->GetTrackID(), motherID1 = 0, nMirror = 0;
    cout << i << " " << nHits << " " << track->GetTrackID() << " " << track->GetChi2() << " " << track->Phi()*TMath::RadToDeg() 
	 << " " << track->Momentum3().Eta() << " " << 1./track->GetParam(4) << " " 
	 << track->Momentum() << " " << track->GetMisses() << endl;
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
      cout << hit->GetLayer();
      Int_t hitId = HitMotherId(hit,motherID,motherID1);
      cout << "-" << hitId;
      if (hit->Index()->GetSize() > 1) cout << "-" << hit->Index()->GetSize();
      if (hit->IsMirror()) cout << "(M)";
      if (motherID1 != motherID) nWrong++;
      else if (hit->IsMirror()) ++nMirror;
    }
    cout << "\n" << nWrong << " " << track->GetTrackID() << " " << motherID << " " << nMirror << endl;
    track->SetNofWrong(nWrong);
    track->SetLastLay(((MpdKalmanHit*)track->GetHits()->First())->GetLayer());
  }
  fTracks->Compress();
}

//__________________________________________________________________________
Int_t MpdEctTrackFinderTof::HitMotherId(MpdKalmanHit* hit, Int_t idM, Int_t &id1)
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
void MpdEctTrackFinderTof::SetTrackID(MpdEctKalmanTrack* track)
{
  /// Set track ID as ID of majority of its hits

  multiset<Int_t> ids;
  multiset<Int_t>::iterator it;
  Int_t nHits = track->GetNofHits(), nMax = 0, idMax = -1;
  TObjArray *hits = track->GetHits();

  for (Int_t i = 0; i < nHits; ++i) {
    MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(i);
    MpdStrawendcapPoint *p = (MpdStrawendcapPoint*) fEctHits->UncheckedAt(hit->GetIndex());
    Int_t id = p->GetTrackID();
    ids.insert(id);
  }

  for (it = ids.begin(); it != ids.end(); ++it) {
    Int_t nC = ids.count(*it);
    if (nC > nMax) {
      nMax = nC;
      idMax = *it;
    }
  }

  if (nMax > 1) track->SetTrackID(idMax);
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::SelectTracks(Int_t iPass)
{
  /// Do track selection and compute shared hit multiplicities

  //if (iPass) return;
  Int_t nFound = fTracks->GetEntriesFast(), nHitsTot = fKHits->GetEntriesFast();
  Int_t *index = new Int_t [nHitsTot];
  for (Int_t i = 0; i < nHitsTot; ++i) index[i] = 0;

  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTracks->UncheckedAt(i);
    if (track->GetTofIndex() < 0) continue; // track from TPC
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
  delete [] index;

  // Remove ghost tracks (with too many shared hits)
  nFound = fTracks->GetEntriesFast();
  Int_t *nh = new Int_t [nFound];
  Double_t *etas = new Double_t [nFound];
  index = new Int_t [nFound];
  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTracks->UncheckedAt(i);
    if (track->GetTofIndex() < 0) continue; // track from TPC
    nh[i] = track->GetNofHits();
    etas[i] = TMath::Abs (track->Momentum3().Eta());
  }
  TMath::Sort(nFound,nh,index,kFALSE); // in ascending order
  //TMath::Sort(nFound,etas,index,kTRUE); // in descending order

  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTracks->UncheckedAt(index[i]);
    if (track->GetTofIndex() < 0) continue; // track from TPC
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
      }
      fTracks->RemoveAt(index[i]);
      cout << " *** Removing track (too many shared hits): " << i << " " 
	   << track->GetTrackID() << " " << nHits << " " << nShared << endl;
    }
  }
  fTracks->Compress();

  // Apply rejection based on the number of missing hits
  nFound = fTracks->GetEntriesFast();
  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTracks->UncheckedAt(i);
    if (track->GetTofIndex() < 0) continue; // track from TPC
    Int_t nHits = track->GetNofHits();
    if (track->GetMisses() >= nHits) {
      fTracks->RemoveAt(i);
      cout << " *** Removing track (too many misses): " << i << " "
	   << track->GetTrackID() << " " << nHits << " " << track->GetMisses() << endl;
    }
  }
  fTracks->Compress();

  delete [] nh;
  delete [] etas;
  delete [] index;
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::MatchTpc()
{
  /// Match ETOF-ECT track with tracks from TPC.

  Int_t nTpc = fTpcTracks->GetEntriesFast();//, nMatch = 0;
  map<Int_t,pair<Double_t,Double_t> > tpcM, ectM;
  set<Int_t> tpcEctS; 

  cout << " TPC tracks: " << nTpc << endl;
  for (Int_t i = 0; i < nTpc; ++i) {
    MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) fTpcTracks->UncheckedAt(i);

    // Exclude TPC tracks used for ECT tracking
    //if (tr == 0x0 || tr->GetUniqueID() == 1) continue;
    if (tr == 0x0) continue;
    if (tr->GetUniqueID() == 1) tpcEctS.insert(i);

    //if (tr->GetParam(3) < 0) continue; // !!! going in backward hemisphere
    if (TMath::Abs(tr->GetParam(3)) < 0.3 ||
	tr->Momentum() < 0.07) continue; // !!! track in central region or with low momentum

    MpdEctKalmanTrack tpc = MpdEctKalmanTrack(i, *tr);
    /*2-may-12
    Double_t distFirst = ((MpdKalmanHit*) tpc.GetHits()->First())->GetDist();
    Double_t distLast = ((MpdKalmanHit*) tpc.GetHits()->Last())->GetDist();
    tpc.SetPos(TMath::Min(distFirst,distLast));
    tpc.SetPosNew(tpc.GetPos());
    Bool_t ok = MpdKalmanFilter::Instance()->Refit(&tpc, -1); // from last to first hit
    //cout << i << " " << tpc.GetHits()->GetEntriesFast() << " " << tpc.GetParamNew(3) << " " << tpc.GetParamNew(1) << " " << tpc.GetPosNew() << endl;
    //if (!ok || tpc.GetParamNew(3) < 0) { 
    if (!ok) { 
      // propagation failure or going in backward hemisphere
      //fTrackCand->RemoveLast();
      //--nCand;
      continue;
    }
    */
    //cout << i << " " << tpc.GetHits()->GetEntriesFast() << " " << tpc.GetType() << " " << tpc.GetParamNew(1) 
    // << " " << tpc.GetParamNew(3) << " " << tpc.GetPosNew() << endl;

    // Propagate to TPC end-plate
    MpdKalmanHit hit;
    hit.SetType(MpdKalmanHit::kFixedZ);
    Double_t zEnd = 150.; // get it from geometry !!!
    hit.SetPos(TMath::Sign(zEnd,tpc.GetParam(3)));
    MpdKalmanFilter::Instance()->PropagateToHit(&tpc, &hit, kFALSE);
    Double_t rad = tpc.GetParamNew(0) * tpc.GetParamNew(0) + tpc.GetParamNew(1) * tpc.GetParamNew(1);
    rad = TMath::Sqrt(rad);
    //cout << tpc.GetTrackID() << " " << tpc.GetType() << " " << tpc.GetPosNew() << " " << tpc.GetParamNew(3) << " " << rad << endl;
    if (rad > 100.) continue; // outer radius of TPC
    tpcM.insert(pair<Int_t,pair<Double_t,Double_t> > (i,pair<Double_t,Double_t>(tpc.GetParamNew(0),tpc.GetParamNew(1))));
  }

  Int_t nFound = fTracks->GetEntriesFast();
  cout << " ECT tracks: " << nFound << endl;
  for (Int_t i = 0; i < nFound; ++i) {
    MpdEctKalmanTrack *tr = (MpdEctKalmanTrack*) fTracks->UncheckedAt(i);

    if (tr->IsFromTpc()) continue; // track from TPC

    // Propagate to TPC end-plate
    MpdKalmanHit hit;
    hit.SetType(MpdKalmanHit::kFixedZ);
    Double_t zEnd = 150.; // get it from geometry !!!
    //hit.SetPos(TMath::Sign(zEnd,tr->GetPosNew()));
    hit.SetPos(TMath::Sign(zEnd,tr->GetParam(3)));
    MpdKalmanFilter::Instance()->PropagateToHit(tr, &hit, kTRUE);
    Double_t rad = tr->GetParamNew(0) * tr->GetParamNew(0) + tr->GetParamNew(1) * tr->GetParamNew(1);
    rad = TMath::Sqrt(rad);
    //cout << tr->GetTrackID() << " " << tr->GetPosNew() << " " << tr->GetParamNew(3) << " " << rad << endl;
    //if (rad < 34.) continue; // inner radius of TPC
    Double_t thick = 0.06; // TPC end-plate material thickness in rad. lengths
    PassWall(tr, thick);
    ectM.insert(pair<Int_t,pair<Double_t,Double_t> > (i,pair<Double_t,Double_t>(tr->GetParamNew(0),tr->GetParamNew(1))));
  }  

  map<Int_t,pair<Double_t,Double_t> >::iterator it, it1;
  multimap<Int_t,Int_t> matchM;
  Double_t rMax = 0.5; // 5 mm matching radius

  for (it = ectM.begin(); it != ectM.end(); ++it) {
    MpdEctKalmanTrack *tr = (MpdEctKalmanTrack*) fTracks->UncheckedAt(it->first);
    Int_t id = tr->GetTrackID(), dir = Int_t(tr->GetPosNew());
    for (it1 = tpcM.begin(); it1 != tpcM.end(); ++it1) {
      MpdTpcKalmanTrack *tr1 = (MpdTpcKalmanTrack*) fTpcTracks->UncheckedAt(it1->first);
      Int_t id1 = tr1->GetTrackID(), dir1 = Int_t(tr1->GetParam(3));
      if (dir * dir1 < 0) continue; // different arms
      Double_t dx = it->second.first - it1->second.first;
      Double_t dy = it->second.second - it1->second.second;
      Double_t dr = dx * dx + dy * dy;
      dr = TMath::Sqrt (dr);
      if (id == id1) {
	// Exact ID matching
	Double_t r = it->second.first * it->second.first + it->second.second * it->second.second;
	if (lun2) fprintf(lun2,"%6d %10.3e %10.3e %1d %3d %3d\n", id, TMath::Sqrt(r), dr, tr1->GetUniqueID(),
			  tr->GetNofHits(), tr1->GetNofHits());
      }

      if (dr < rMax) matchM.insert(pair<Int_t,Int_t>(it->first,it1->first));
    }
  }

  cout << " Number of ETOF tracks to TPC: " << ectM.size() << ", potential matches: " << matchM.size() << endl;

  // Find matches
  multimap<Int_t,Int_t>::iterator itm, itm1;
  pair<multimap<Int_t,Int_t>::iterator,multimap<Int_t,Int_t>::iterator> ret;
  map<Int_t,MpdEctKalmanTrack> trackM;
  TMatrixDSym pointWeight(5);
  TMatrixD param(5,1);
  Int_t indx0 = -1;

  for (itm = matchM.begin(); itm != matchM.end(); ++itm) {
    Int_t indx = itm->first;
    if (indx == indx0) continue;
    indx0 = indx;
    MpdEctKalmanTrack *tr = (MpdEctKalmanTrack*) fTracks->UncheckedAt(indx);
    MpdEctKalmanTrack ect(*tr);

    ret = matchM.equal_range(indx);
    trackM.clear();
    for (itm1 = ret.first; itm1 != ret.second; ++itm1) {
      MpdTpcKalmanTrack *tr1 = (MpdTpcKalmanTrack*) fTpcTracks->UncheckedAt(itm1->second);
      TClonesArray *hits = tr1->GetTrHits();
      Int_t nHits = hits->GetEntriesFast(), lay = -1, lay0 = -1;

      MpdEctKalmanTrack trTMP = ect;
      trTMP.GetHits()->Clear();

      // Follow track thru TPC
      map<Int_t,Double_t>& stepMap = tr1->GetSteps();
      for (Int_t j = 0; j < nHits; ++j) {
	MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(j);
	lay0 = lay;
	lay = hit->GetLayer();
	Double_t stepBack = -1.0;
	if (stepMap.find(lay) != stepMap.end() && stepMap.find(lay0) != stepMap.end()) 
	  stepBack = TMath::Abs (stepMap[lay] - stepMap[lay0]);
	if (!MpdKalmanFilter::Instance()->PropagateToHit(&trTMP,hit,kTRUE,kTRUE,stepBack)) {
	  trTMP.SetNodeNew(trTMP.GetNode()); // restore node
	  continue;
	}
	Double_t dChi2 = MpdKalmanFilter::Instance()->FilterHit(&trTMP,hit,pointWeight,param);
	//cout << j << " " << dChi2 << " " << hit->GetPos() << endl;
	if (dChi2 > fgkChi2Cut * 4) continue;
	trTMP.SetParam(param);
	trTMP.SetParamNew(param);
	(*trTMP.GetWeight()) += pointWeight;
	trTMP.GetHits()->Add(hit);
	trTMP.SetChi2(trTMP.GetChi2()+dChi2);
	// Set params at hit
	trTMP.SetParamAtHit(*trTMP.GetParamNew());
	trTMP.SetWeightAtHit(*trTMP.GetWeight());
	trTMP.SetLengAtHit(trTMP.GetLength());
	trTMP.SetPosAtHit(trTMP.GetPosNew());
      }
      trTMP.SetPartID(tr1->GetPartID());
      trackM.insert(pair<Int_t,MpdEctKalmanTrack>(itm1->second,trTMP));
    }
    
    // Select the best track (with max number of hits)
    map<Int_t,MpdEctKalmanTrack>::iterator iter;
    Int_t nHitsMax = 0, iMax = -1;
    for (iter = trackM.begin(); iter != trackM.end(); ++iter) {
      Int_t nhits = iter->second.GetNofHits();
      if (nhits == 0) continue;
      if (nhits > nHitsMax) {
	nHitsMax = nhits;
	iMax = iter->first;
      } else if (nhits == nHitsMax) {
	if (iter->second.GetChi2() < trackM[iMax].GetChi2()) {
	  iMax = iter->first;
	  nHitsMax = nhits;
	}
      }
    }
    if (iMax >= 0) {
      cout << " ECT-TPC match found: " << tr->GetTrackID() << " " << trackM[iMax].GetTrackID() << " " << iMax << endl;
      trackM[iMax].SetMatchTpc();
      trackM[iMax].SetTpcIndex(iMax);
      // Check if there is a TPC-ECT track (clone)
      Bool_t removeThis = kFALSE;
      if (tpcEctS.find(iMax) != tpcEctS.end()) {
	for (Int_t j = 0; j < nFound; ++j) {
	  MpdEctKalmanTrack *tr1 = (MpdEctKalmanTrack*) fTracks->UncheckedAt(j);
	  if (tr1 == 0x0) continue;
	  if (!tr1->IsFromTpc()) continue; // track not from TPC
	  if (tr1->GetTpcIndex() != iMax) continue;
	  cout << tr->GetNofHits() << " " << tr1->GetNofHits() << " " << tr1->GetChi2() << endl;
	  if (tr->GetNofHits() > tr1->GetNofHits()) fTracks->RemoveAt(j);
	  else if (tr->GetNofHits() < tr1->GetNofHits()) removeThis = kTRUE;
	  else if (tr->GetChi2() > tr1->GetChi2()) removeThis = kTRUE;
	  break;
	}
      }
      
      cout << tr->GetChi2() << " " << tr->GetNofTrHits() << endl;
      fTracks->RemoveAt(indx);
      if (!removeThis) {
	// Add TPC hits
	TClonesArray *oldHits = trackM[iMax].GetTrHits();
	TObjArray *newHits = trackM[iMax].GetHits();
	Int_t nAdd = newHits->GetEntriesFast(), nOld = oldHits->GetEntriesFast();
	MpdEctKalmanTrack *newTrack = new ((*fTracks)[itm->first]) MpdEctKalmanTrack(trackM[iMax]);
	oldHits = newTrack->GetTrHits();
	for (Int_t ih = 0; ih < nAdd; ++ih) {
	  new ((*oldHits)[nOld++]) MpdKalmanHit (*((MpdKalmanHit*)newHits->UncheckedAt(ih)));
	}
	cout << trackM[iMax].GetChi2() << " " << trackM[iMax].GetNofHits() << " " << trackM[iMax].GetNofTrHits() << endl;
	newHits = newTrack->GetHits();
	newHits->Clear();
	for (Int_t ih = 0; ih < nOld; ++ih) newHits->Add(oldHits->UncheckedAt(ih));
	newTrack->SetNofHits(newTrack->GetNofTrHits());
	newTrack->SetParamAtHit(*trackM[iMax].GetParamAtHit());
	newTrack->SetWeightAtHit(*trackM[iMax].GetWeightAtHit());
	newTrack->SetLengAtHit(trackM[iMax].GetLengAtHit());
	newTrack->SetPosAtHit(trackM[iMax].GetPosAtHit());
      }
      //fTracks->RemoveAt(indx);
    }
  } // for (itm = matchM.begin(); ...
  trackM.clear();
  fTracks->Compress();

}

//__________________________________________________________________________
void MpdEctTrackFinderTof::PassWall(MpdEctKalmanTrack *track, Double_t thick)
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
void MpdEctTrackFinderTof::GoToBeamLine()
{
  /// Propagate tracks to the beam line (mult. scat. and energy loss corrections)

  const Int_t nR = 7;
  const Double_t rad[nR] = {27.035, 27.57, 28.105, 30.64, 33.15, 33.665, 34.178};
  const Double_t dx[nR] = {0.0031, 0.00085, 0.0031, 0.0003, 0.001, 0.00085, 0.001};
  Bool_t ok = 0;
  Int_t nReco = fTracks->GetEntriesFast();

  for (Int_t i = 0; i < nReco; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTracks->UncheckedAt(i);
    //if (track->IsFromTpc()) continue; // track from TPC
    if (track->IsFromTpc()) MergeWithTpc(track); // merge track with TPC track

    // Check track radial position
    Double_t rTr = track->GetPosNew();
    if (track->GetType() == MpdKalmanTrack::kEndcap) {
      rTr = track->GetParamNew(0) * track->GetParamNew(0) + track->GetParamNew(1) * track->GetParamNew(1);
      rTr = TMath::Sqrt (rTr);
    }

    MpdKalmanHit hit;
    hit.SetDetectorID(-999);
    hit.SetType(MpdKalmanHit::kFixedR);

    // Check if goes inward
    MpdEctKalmanTrack trTmp = *track;
    hit.SetPos(0.0);
    Double_t leng = trTmp.GetLength();
    //ok = MpdKalmanFilter::Instance()->PropagateToHit(&trTmp,&hit,kFALSE);
    ok = MpdKalmanFilter::Instance()->PropagateToHit(&trTmp,&hit,kTRUE);
    if (!ok || TMath::Abs(trTmp.GetParamNew(1)) > 149.9) {
      fTracks->RemoveAt(i);
      cout << " *** Removing track (not going inward): " << i << " "
	   << track->GetTrackID() << " " << endl;
      continue;
    } else if (trTmp.GetPosNew() >= rad[nR-1] || rTr < rad[0]) {
      // Too large DCA or inside TPC inner shell
      *track = trTmp;
      track->SetPos(trTmp.GetPosNew());
      track->SetParam(*track->GetParamNew()); // !!! track params at PCA
      track->Weight2Cov();
      if (track->IsFromTpc() || track->IsMatchTpc()) track->SetLengAtHit(trTmp.GetLength()-leng);
      else track->SetLengAtHit(trTmp.GetLength()-trTmp.GetLengAtHit());
      track->SetLength(leng);
      continue;
    }
    if (track->IsFromTpc() || track->IsMatchTpc()) track->SetLengAtHit(trTmp.GetLength()-leng);
    else track->SetLengAtHit(trTmp.GetLength()-trTmp.GetLengAtHit());

    for (Int_t j = nR-1; j >= 0; --j) {
      if (rTr < rad[j]) continue;
      hit.SetPos(rad[j]);
      //MpdKalmanFilter::Instance()->GetGeo()->SetGlobalPos(&hit,TVector3(rad[j],0.,0.),kTRUE); // put hit in GeoScheme
      //ok = MpdKalmanFilter::Instance()->PropagateToHit(track,&hit,kTRUE);
      ok = MpdKalmanFilter::Instance()->PropagateToHit(track,&hit,!track->IsFromTpc());
      //cout << ok << " " << track->GetPos() << endl;
      if (!ok) {
	fTracks->RemoveAt(i);
	cout << " *** Removing track (not passing TPC inner shell): " << i << " "
	     << track->GetTrackID() << " " << endl;
	break;
      }
      //continue;
      if (TMath::Abs(track->GetParamNew(1)) > 150.) continue; // outside TPC Z-dimensions
      // Add multiple scattering
      TMatrixDSym *cov = track->Weight2Cov();
      Double_t th = track->GetParamNew(3);
      Double_t cosTh = TMath::Cos(th);
      Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, dx[j]);
      (*cov)(2,2) += (angle2 / cosTh / cosTh);
      (*cov)(3,3) += angle2;
      //cov->Print();
      Int_t iok = 0;
      MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
      track->SetWeight(*cov);
    }
    if (!ok) continue;
    /* TMatrixDSym *cov = */ track->Weight2Cov();

    // Update track
    track->SetParam(*track->GetParamNew());
    track->SetPos(track->GetPosNew());
    // Correct for energy loss
    //*
    Double_t theta = TMath::PiOver2() - track->GetParam(3);
    Double_t pt = 1. / track->GetParam(4);
    Double_t ptCor = CorrectForLoss(TMath::Abs(pt), theta, track->GetTrackID()); // ionization loss correction
    track->SetParam(4,TMath::Sign(1./ptCor,pt));
    track->SetParamNew(*track->GetParam());
    //*/

    //*
    hit.SetPos(0.);
    hit.SetMeas(0,track->GetParam(2));
    //hit.SetRphi(track->GetParam(2)); // track Phi - check if it is necessary !!!!!!!
    //MpdKalmanFilter::Instance()->GetGeo()->SetGlobalPos(&hit,TVector3(0.,0.,0.),kTRUE); // put hit in GeoScheme
    //*/
    //cout << i << " " << track->GetTrackID() << " " << track->GetLength() << " " << ((MpdKalmanHitR*)track->GetHits()->First())->GetLength() << endl;
    Double_t pos = track->GetPosNew();
    //Double_t pos = ((MpdKalmanHitR*)track->GetHits()->Last())->GetR(); // position at last hit
    //MpdKalmanFilter::Instance()->PropagateParamR(track, &hit, kTRUE);
    //cout << track->Momentum3().Eta() << endl;
    //track->GetCovariance()->Print();
    //(track->GetWeightAtHit()->Invert()).Print();
    //MpdKalmanFilter::Instance()->PropagateToHit(track, &hit, kTRUE);
    MpdKalmanFilter::Instance()->PropagateToHit(track, &hit, !track->IsFromTpc());
    track->SetPos(pos); // save position after passing inner shell
    track->SetParam(*track->GetParamNew()); // !!! track params at PCA
  }
  fTracks->Compress();
}

//__________________________________________________________________________
Double_t MpdEctTrackFinderTof::CorrectForLoss(Double_t pt, Double_t the, Int_t id)
{
  /// Apply momentum correction due to ionization loss in the beam pipe (0.5 mmof Al)
  /// and TPC inner shell

  //return pt; ///
  const Int_t nDim = 13;
  //const Double_t mom[nDim]={0.075,0.100,0.200,0.300,0.400,0.500,0.700,0.900,
			    //1.200,1.500,2.000,2.500,3.500};
  const Double_t tkin[nDim]={0.019,0.032,0.104,0.191,0.284,0.380,0.574,0.771,
			     1.069,1.367,1.865,2.364,3.363};
  //const Double_t dedxp[nDim]={3.125,1.950,0.975,0.800,0.800,0.750,0.750,0.750,
			      //0.750,0.750,0.750,0.775,0.775}; // peak
  const Double_t dedxm[nDim]={3.195,2.085,1.100,0.952,0.918,0.898,0.898,0.890,
			      0.887,0.886,0.901,0.904,0.913}; // mean
  const Double_t piMass = 0.13957, eMass = 0.00051;
  Double_t dt, mass;

  Int_t elec = GetParticleId(id);
  if (elec) mass = eMass;
  else mass = piMass;
  Double_t mass2 = mass * mass;

  Double_t p = pt / TMath::Sin(the);
  Double_t t = TMath::Sqrt (p*p + mass2) - mass;

  if (elec) dt = 1.45; // 1.45 MeV loss for electrons
  
  else {
    //if (t < tkin[0]) dt = dedxm[0];
    if (t < tkin[0]) {
      dt = dedxm[0] + (dedxm[1]-dedxm[0])/(tkin[1]-tkin[0]) * (t-tkin[0]);
      dt = TMath::Min (dt,10.);
    }
    else if (t > tkin[nDim-1]) dt = dedxm[nDim-1];
    else {
      // Binary search
      Int_t i1 = 0, i2 = nDim-1, i = i2;
      do {
	i = i1 + (i2-i1) / 2;
	if (t > tkin[i]) i1 = i;
	else i2 = i;
	//cout << i << " " << i1 << " " << i2 << " " << t << " " << tkin[i1] << " " << tkin[i2] << endl;
      } while (i2 - i1 > 1);
      // Linear interpolation
      dt = dedxm[i1] + (dedxm[i2]-dedxm[i1])/(tkin[i2]-tkin[i1]) * (t-tkin[i1]);
    }
  }
  dt /= TMath::Sin(the);

  t += dt * 1.e-3;
  Double_t e = mass + t;
  p = TMath::Sqrt (e*e - mass2);
  pt = p * TMath::Sin(the);
  return pt;
}

//__________________________________________________________________________
Int_t MpdEctTrackFinderTof::GetParticleId(Int_t id)
{
  /// Particle ID:
  /// !!! based on MC information at the moment !!!

  MpdMCTrack* mcTr = (MpdMCTrack*) fMCTracks->UncheckedAt(id);
  Int_t pdg = mcTr->GetPdgCode();
  if (TMath::Abs(pdg) == 11) return 1;
  return 0;
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::MergeWithTpc(MpdEctKalmanTrack *track)
{
  /// Merge track with TPC track

  track->SetDirection(MpdKalmanTrack::kInward);
  cout << track->GetTrackID() << " " << track->GetNofHits() << " " << track->GetChi2() << " " << track->Pt()*track->Charge() << " " ;//<< endl;

  //Bool_t ok = MpdKalmanFilter::Instance()->Refit(track, 1); // from first to last hit
  //*
  track->GetHits()->Sort();
  track->SetChi2(0.);
  TMatrixDSym *w = track->GetWeight();
  Int_t nHits = track->GetNofHits();
  for (Int_t i = 0; i < 5; ++i) {
    for (Int_t j = i; j < 5; ++j) {
      if (j == i) (*w)(i,j) /= nHits;
      else (*w)(i,j) = (*w)(j,i) = 0.;
    }
  }

  TMatrixD param(5,1);
  TMatrixDSym weight(5), pointWeight(5);
  Int_t ibeg = 0, iend = nHits, nAcc = 0;

  MpdKalmanHit *hit = 0x0;
  for (Int_t i = ibeg; i != iend; ++i) {
    hit = (MpdKalmanHit*) track->GetHits()->UncheckedAt(i);
    if (!MpdKalmanFilter::Instance()->PropagateToHit(track, hit, kFALSE)) continue;
    ++nAcc;
    PassWall(track, 0.001);
    Double_t dChi2 = MpdKalmanFilter::Instance()->FilterHit(track, hit, pointWeight, param);
    track->SetChi2(track->GetChi2()+dChi2);
    weight = *track->GetWeight();
    weight += pointWeight;
    track->SetWeight(weight);
    track->SetParamNew(param);
    //cout << i << " " << dChi2 << " " << 1./track->GetParamNew(4) << endl;
    // Set params at hit 9-may-12
    track->SetParamAtHit(*track->GetParamNew());
    track->SetWeightAtHit(*track->GetWeight());
    track->SetLengAtHit(track->GetLength());
    track->SetPosAtHit(track->GetPosNew());
  }
  //*/
  cout << nAcc << " " << track->GetChi2() << " " << track->Pt()*track->Charge() << endl;
  
  // Propagate to TPC end-plate
  MpdKalmanHit hitTmp;
  hitTmp.SetType(MpdKalmanHit::kFixedZ);
  Double_t zEnd = 150.; // get it from geometry !!!
  //hitTmp.SetPos(TMath::Sign(zEnd,track->GetPos()));
  hitTmp.SetPos(TMath::Sign(zEnd,track->GetParam(3)));
  MpdKalmanFilter::Instance()->PropagateToHit(track, &hitTmp, kFALSE);
  Double_t thick = 0.06; // material thickness in rad. lengths
  PassWall(track,thick);

  MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) fTpcTracks->UncheckedAt(track->GetTpcIndex());
  TClonesArray *hits = tr->GetTrHits();
  /*Int_t*/ nHits = hits->GetEntriesFast();

  //TMatrixDSym pointWeight(5);
  pointWeight = 0.;
  //TMatrixD param(5,1);
  TClonesArray *oldHits = track->GetTrHits();
  Int_t nOld = oldHits->GetEntriesFast(), nAdd = 0, lay0 = -1, lay = -1; 
  nAcc = nOld;
  map<Int_t,Double_t>& stepMap = track->GetSteps();
  for (Int_t j = 0; j < nHits; ++j) {
    MpdKalmanHit *hit1 = (MpdKalmanHit*) hits->UncheckedAt(j);
    lay0 = lay;
    lay = hit1->GetLayer();
    //if (!MpdKalmanFilter::Instance()->PropagateToHit(track,hit,kFALSE)) continue;
    //cout << j << " " << hit->GetMeas(1) << endl;
    Double_t stepBack = -1.0;
    if (stepMap.find(lay) != stepMap.end() && stepMap.find(lay0) != stepMap.end()) 
      stepBack = TMath::Abs (stepMap[lay] - stepMap[lay0]);
    if (!MpdKalmanFilter::Instance()->PropagateToHit(track,hit1,kFALSE,kTRUE,stepBack)) {
      track->SetNodeNew(track->GetNode()); // restore node
      continue;
    }
    Double_t dChi2 = MpdKalmanFilter::Instance()->FilterHit(track,hit1,pointWeight,param);
    //cout << j << " " << hit->GetDist() << " " << dChi2 << endl;
    if (dChi2 > fgkChi2Cut * 4) continue;
    //track->SetParam(param);
    track->SetParamNew(param);
    (*track->GetWeight()) += pointWeight;
    track->GetHits()->Add(hit1);
    track->SetChi2(track->GetChi2()+dChi2);
    new ((*oldHits)[nOld++]) MpdKalmanHit (*hit1);
    ++nAdd;
    // Set params at hit
    track->SetParamAtHit(*track->GetParamNew());
    track->SetWeightAtHit(*track->GetWeight());
    track->SetLengAtHit(track->GetLength());
    track->SetPosAtHit(track->GetPosNew());
  }

  if (nAdd == 0) {
    cout << " Merge TPC: merging failed !!! " << endl;
    track->SetFlag(track->GetFlag() ^ 1); // clear "From_TPC" flag
    return;
  }

  /*TClonesArray *oldHits = track->GetTrHits();
  Int_t nAdd = hits->GetEntriesFast(), nOld = oldHits->GetEntriesFast();
  for (Int_t ih = 0; ih < nAdd; ++ih) {
    new ((*oldHits)[nOld++]) MpdKalmanHit (*(MpdKalmanHit*)hits->UncheckedAt(ih));
    }*/
  track->SetNofHits(track->GetNofTrHits());
  cout << " Merge TPC: " << nHits << " " << nAdd << " " << nAcc << " " << track->GetNofHits() << " " 
       << track->GetNofTrHits() << " " << tr->GetTrackID() << " " << track->GetTrackID() << " " 
       << tr->GetParam(3) << " " << track->GetChi2() << endl;
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::Smooth()
{
  /// Primary vertex constraints: smooth tracks from primary vertex (update 
  /// momentum and track length - covariance matrix is not updated !!!)

  const Double_t cutChi2 = 15.;

  MpdVertex *vtx = 0x0, *vtxOld = 0x0;
  Int_t nPart = 0;
  if (fPrimVtx) {
    vtx = (MpdVertex*) fPrimVtx->UncheckedAt(0);
    nPart = vtx->GetNTracks();
    // Save original vertex
    vtxOld = new MpdVertex();
    *vtxOld = *vtx;
  } 
  if (vtx == 0x0 || nPart < 20) {
    // Use default vertex
    TMatrixFSym cov(3);
    cov(0,0) = cov(1,1) = cov(2,2) = 1.e-6;
    MpdVertex vertDef("Def. vertex", "DefVertex", 0, 0, 0, 0, 0, 0, cov);
    if (vtx == 0x0) vtx = new MpdVertex();
    *vtx = vertDef;
    if (!fPrimVtx) {
      // No vertex
      /*
      TClonesArray tmp("MpdVertex",1);
      fPrimVtx = &tmp;
      MpdVertex *v = new (tmp[0]) MpdVertex();
      */
      TClonesArray *tmp = new TClonesArray("MpdVertex",1);
      fPrimVtx = tmp;
      MpdVertex *v = new ((*tmp)[0]) MpdVertex();
      *v = *vtx;
    }
  }

  MpdKfPrimaryVertexFinder vFinder;
  vFinder.SetVertices(fPrimVtx);
  vFinder.SetTracks(fTracks);
  vFinder.Chi2Vertex();

  // Select tracks close to the vertex
  Int_t nTr = fTracks->GetEntriesFast(), nOk = 0;
  TArrayI indsNew(nTr); //, indsOld = *(vtx->GetIndices());
  Double_t *lengs = new Double_t [nTr];
  for (Int_t itr = 0; itr < nTr; ++itr) {
    MpdKalmanTrack *track = (MpdKalmanTrack*) fTracks->UncheckedAt(itr);
    lengs[itr] = track->GetLength();
    if (track->GetChi2Vertex() > cutChi2) continue;
    indsNew[nOk++] = itr;
  }
  indsNew.Set(nOk);
  vtx->SetIndices(&indsNew);

  vFinder.Smooth();

  // Set lower limit to the track length (distance to ETOF hit)
  for (Int_t itr = 0; itr < nTr; ++itr) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTracks->UncheckedAt(itr);
    if (track->GetChi2Vertex() < cutChi2) track->SetLengAtHit(track->GetLengAtHit() + track->GetLength() - lengs[itr]);
    //if (track->GetChi2Vertex() > cutChi2) continue;
    //if (track->GetTofIndex() < 0) continue;
    Double_t dist = TMath::Abs (((MpdTofHit*) fTofHits->First())->GetZ());
    if (track->GetTofIndex() < 0) {
      if (track->GetChi2Vertex() > cutChi2) continue;
    } else {
      MpdTofHit *hit = (MpdTofHit*) fTofHits->UncheckedAt(track->GetTofIndex());
      Double_t dz = hit->GetZ() - vtx->GetZ();
      dist = hit->GetX() * hit->GetX() + hit->GetY() * hit->GetY() + dz * dz;
      dist = TMath::Sqrt (dist);
    }
    if (track->GetLength() < dist) track->SetLength(dist);
  }

  delete [] lengs;
  if (vtxOld) {
    // Restore old vertex
    *vtx = *vtxOld;
    delete vtxOld;
  } else {
    fPrimVtx->Delete();
    fPrimVtx = 0x0;
    delete vtx;
  }
}

//__________________________________________________________________________
void MpdEctTrackFinderTof::GoOutward()
{
  /// Propagate tracks outward (for matching with outer detectors)

  const Double_t cutChi2 = 15.;
  Bool_t ok = 0;
  Int_t nReco = fTracks->GetEntriesFast();
  FairMCEventHeader *mchead = (FairMCEventHeader*) FairRootManager::Instance()->GetObject("MCEventHeader.");

  for (Int_t i = 0; i < nReco; ++i) {
    MpdEctKalmanTrack *track = (MpdEctKalmanTrack*) fTracks->UncheckedAt(i);
    //if (track->IsFromTpc()) continue; // track from TPC
    MpdEctKalmanTrack tr = *track;
    tr.SetParam(*tr.GetParamAtHit());
    tr.SetParamNew(*tr.GetParamAtHit());
    tr.SetWeight(*tr.GetWeightAtHit());
    tr.SetLength(tr.GetLengAtHit());
    MpdKalmanHit *hit = (MpdKalmanHit*) tr.GetTrHits()->Last();
    //2-may-12 tr.SetPos(hit->GetPos());
    tr.SetPos(tr.GetPosAtHit());
    tr.SetPosNew(tr.GetPos());
    tr.SetDirection(MpdKalmanTrack::kOutward);
    if (hit->GetType() == MpdKalmanHit::kFixedP) {
      tr.SetNode(MpdKalmanFilter::Instance()->GetGeo()->Path(hit->GetDetectorID()));
      tr.SetNodeNew(MpdKalmanFilter::Instance()->GetGeo()->Path(hit->GetDetectorID()));
    }
    //cout << "MpdEctTrackFinderTof::GoOutward(): " << i << " " << tr.GetTrackID() << " " << tr.GetPos() << " " 
    // << tr.GetNofHits() << " " << tr.GetNofTrHits() << " " << tr.GetChi2() << endl;
    //if (tr.GetNofHits() != tr.GetNofTrHits()) exit(0);

    // Just for test - from EtofMatching
    /*
    tr.SetPos(track->GetPosNew());
    tr.SetPosNew(track->GetPosNew());
    tr.SetParam(*track->GetParam());
    tr.SetParamNew(*track->GetParam());
    tr.ReSetWeight();
    tr.SetLength(0);
    tr.SetDirection(MpdKalmanTrack::kOutward);
    MpdKalmanFilter::Instance()->Refit(&tr,-1,kTRUE);
    tr.SetPos(tr.GetPosNew());
    tr.SetParam(*tr.GetParamNew());
    */
    // Check if there are hits from TPC
    TClonesArray *trhits = tr.GetTrHits();
    trhits->Sort();
    Int_t nTpc = 0, nEct = 0, nHits = trhits->GetEntriesFast();
    for (Int_t jh = 0; jh < nHits; ++jh) {
      MpdKalmanHit *hit1 = (MpdKalmanHit*) trhits->UncheckedAt(jh);
      if (hit1->GetType() == MpdKalmanHit::kFixedR || hit1->GetType() == MpdKalmanHit::kFixedP) ++nTpc;
      else ++nEct;
    }
    //*
    TObjArray *hits = tr.GetHits();
    Bool_t ok1 = 0;
    if (nTpc > 0) {
      // TPC
      hits->Clear();
      for (Int_t jh = nEct; jh < nHits; ++jh) hits->Add(trhits->UncheckedAt(jh));
      //ok = MpdKalmanFilter::Instance()->Refit(&tr,-1,kFALSE);
      ok1 = MpdKalmanFilter::Instance()->Refit(&tr,-1,kTRUE);
      if (!ok1) cout << " Go outward: refit failed !!!" << endl; 
      // Go to TPC end-plate
      MpdKalmanHit hitTmp;
      hitTmp.SetType(MpdKalmanHit::kFixedZ);
      Double_t zEnd = 150.; // get it from geometry !!!
      hitTmp.SetPos(TMath::Sign(zEnd,track->GetParam(3)));
      //MpdKalmanFilter::Instance()->PropagateToHit(&tr, &hitTmp, kFALSE);
      MpdKalmanFilter::Instance()->PropagateToHit(&tr, &hitTmp, kTRUE);
      Double_t thick = 0.06; // material thickness in rad. lengths
      PassWall(&tr,thick);
    }
    if (nEct > 0) {
      // ECT
      hits->Clear();
      for (Int_t jh = 0; jh < nEct; ++jh) hits->Add(trhits->UncheckedAt(jh));
      if (nTpc == 0) {
	tr.SetType(MpdKalmanTrack::kEndcap);
	tr.SetChi2(0.);
	//ok = MpdKalmanFilter::Instance()->Refit(&tr,-1,kTRUE);
	TMatrixDSym *w = tr.GetWeight();
	for (Int_t k = 0; k < 5; ++k) {
	  for (Int_t l = k; l < 5; ++l) {
	    if (l == k) (*w)(k,l) /= nEct;
	    else (*w)(k,l) = (*w)(l,k) = 0.;
	  }
	}
      }
      TMatrixDSym pointWeight(5);
      TMatrixD param(5,1);
      for (Int_t jh = nEct - 1; jh >= 0; --jh) {
	MpdKalmanHit *hit1 = (MpdKalmanHit*) hits->UncheckedAt(jh);
	//if (!MpdKalmanFilter::Instance()->PropagateToHit(&tr,hit1,kFALSE)) continue;
	if (!MpdKalmanFilter::Instance()->PropagateToHit(&tr,hit1,kTRUE)) continue;
	PassWall(&tr, 0.001);
	Double_t dChi2 = MpdKalmanFilter::Instance()->FilterHit(&tr,hit1,pointWeight,param);
	//tr.SetParam(param);
	tr.SetParamNew(param);
	(*tr.GetWeight()) += pointWeight;
	tr.SetChi2(tr.GetChi2()+dChi2);
      }
    }
    printf("MpdEctTrackFinderTof::GoOutward(): %6d %2d %2d %10.3f %7.3f %1d \n", tr.GetTrackID(), 
	   nTpc, nEct, tr.GetChi2(), tr.Pt(), ok);
    // Set params at hit
    track->SetParamAtHit(*tr.GetParamNew());
    track->SetWeightAtHit(*tr.GetWeight());
    track->SetLengAtHit(tr.GetLength());
    track->SetPosAtHit(tr.GetPosNew());
    track->GetTrHits()->Sort();
    //*/

    // Below is just for debugging
    if (0) {
      tr.SetParam(*track->GetParamAtHit());
      tr.SetParamNew(*track->GetParamAtHit());
      tr.SetWeight(*track->GetWeightAtHit());
      tr.SetLength(track->GetLengAtHit());
      tr.SetPos(((MpdKalmanHit*)track->GetTrHits()->First())->GetPos());
      tr.SetPosNew(tr.GetPos());
      // Extrapolate to ETOF
      //Double_t zTof = ((MpdTofHit*)fTofHits->First())->GetZ();
      Double_t zTof = ((FairMCPoint*)fTofPoints->First())->GetZ();
      zTof = TMath::Sign (TMath::Abs(zTof),track->GetParam(3));
      MpdKalmanHit hitTmp;
      hitTmp.SetType(MpdKalmanHit::kFixedZ);
      hitTmp.SetPos(zTof);
      hitTmp.SetNofDim(2);
      //MpdKalmanFilter::Instance()->PropagateToHit(&tr, &hitTmp, kFALSE);
      MpdKalmanFilter::Instance()->PropagateToHit(&tr, &hitTmp, kTRUE);
      // Correct min. length
      if (tr.GetChi2Vertex() < cutChi2) {
	Double_t dist = zTof * zTof + tr.GetParamNew(0) * tr.GetParamNew(0);
	dist += (tr.GetParamNew(1) * tr.GetParamNew(1));
	dist = TMath::Sqrt (dist);
	//if (tr.GetLength() < dist) tr.SetLength(dist);
	tr.SetLength(dist);
      }
      // Matching
      /*
	Int_t nTof = fTofHits->GetEntriesFast();
	for (Int_t jh = 0; jh < nTof; ++jh) {
	MpdTofHit *tofHit = (MpdTofHit*) fTofHits->UncheckedAt(jh);
	Int_t nLinks = tofHit->GetNLinks();
	for (Int_t lnk = 0; lnk < nLinks; ++lnk) {
	FairLink link = tofHit->GetLink(lnk);
	if (link.GetType() != MpdTofUtils::IsTofPointIndex) continue;
	FairMCPoint *p = fTofPoints->UncheckedAt(link);
	}
	//hitTmp
	}
      */
      Int_t nTof = fTofPoints->GetEntriesFast();
      TVector3 pos, posTr(tr.GetParamNew(0),tr.GetParamNew(1),tr.GetPosNew());
      for (Int_t jh = 0; jh < nTof; ++jh) {
	FairMCPoint *p = (FairMCPoint*) fTofPoints->UncheckedAt(jh);
	if (p->GetTrackID() != tr.GetTrackID()) continue;
	Double_t dl = tr.GetLength() - p->GetLength();
	p->Position(pos);
	Double_t phi = pos.Phi();
	Double_t r = pos.Pt();
	Double_t phiTr = MpdKalmanFilter::Instance()->Proxim(phi,posTr.Phi());
	if (lun2) fprintf(lun2, "%4d %6d %2d %10.3f %10.3f % 10.3f %10.3f %10.3f %10.3f %10.3f %1d\n", mchead->GetEventID(), tr.GetTrackID(), nTpc, 
			  dl, tr.GetParamNew(0), tr.GetParamNew(1), tr.GetPosNew(), pos.Pt(), 
			  posTr.Pt()-r, (phiTr-phi)*r, tr.IsFromTpc());
	break;
      }
    } // if (0)

  } //for (Int_t i = 0; i < nReco;
}

ClassImp(MpdEctTrackFinderTof);
