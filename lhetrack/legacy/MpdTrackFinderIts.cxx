// -------------------------------------------------------------------------
// -----                  MpdTrackFinderIts source file                -----
// -----                 Created 21/07/09  by A. Zinchenko             -----
// -------------------------------------------------------------------------

/**  MpdTrackFinderIts.h
 *@author A.Zinchenko <Alexander.Zinchenko@jinr.ru>
 **
 ** Track finder in MPD Inner Tracking System (ITS) using seeds from TPC
 **/

#include "MpdStsGeoPar.h"
#include "MpdTrackFinderIts.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanTrack.h"
//#include "MpdKalmanStripHit.h"
#include "MpdStsPoint.h"
#include "MpdItsKalmanTrack.h"
#include "MpdTpcKalmanTrack.h"

#include "FairGeoNode.h"
#include "FairMCPoint.h"
#include "MpdMCTrack.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"
#include "TGeoNode.h"
#include "TGeoTube.h"
#include "TGeoManager.h"
#include "TMath.h"
//#include "TFile.h"
//#include "TLorentzVector.h"
#include "TVector2.h"
#include "TClonesArray.h"
#include <TRandom.h>

#include <iostream>
#include <map>

using std::cout;
using std::endl;
//using std::map;

const Double_t MpdTrackFinderIts::fgkChi2Cut = 20; //20; //100;
//const Double_t MpdTrackFinderIts::fgkChi2Cut = 50; // for track fitting !!!

//__________________________________________________________________________
MpdTrackFinderIts::MpdTrackFinderIts(const char *name, Int_t iVerbose )
  :FairTask(name, iVerbose),
   fGeo(0)
{
  //fKHits = new TClonesArray("MpdKalmanStripHit", 100);
  fKHits = new TClonesArray("MpdKalmanHit", 100);
  fTracks = new TClonesArray("MpdItsKalmanTrack", 100);
  fHistoDir = 0x0;
  fhLays = new TH1F("hLaysITS","ITS layers",10,0,10);
  fLayPointers = 0x0;
}


//__________________________________________________________________________
MpdTrackFinderIts::~MpdTrackFinderIts()
{
  //delete fKHits;
  //delete fTracks;
  //delete fTrackCand;
  delete [] fLayPointers;
  delete fhLays;
}

//__________________________________________________________________________
InitStatus MpdTrackFinderIts::Init()
{
  //return ReInit();
  if (ReInit() != kSUCCESS) return kERROR;

  // Read database
  FairRuntimeDb* rtdb = FairRun::Instance()->GetRuntimeDb();
  MpdStsGeoPar *geoPar = (MpdStsGeoPar*) rtdb->getContainer("MpdStsGeoPar");
  //cout << geoPar << endl;
  TString volName = "sts01 ", path = "";
  TObjArray* sensNodes = geoPar->GetGeoSensitiveNodes();
  //cout << sensNodes->GetEntriesFast() << " " << geoPar->GetGeoPassiveNodes()->GetEntriesFast() << endl;
  Int_t nLay = 5;
  Double_t size = 6.2;
  for (Int_t i = 0; i < nLay; ++i) {
    volName[5] = 97 + i; // 'a', 'b', ..
    FairGeoNode* sensVol = (FairGeoNode*) (sensNodes->FindObject(volName));
    if (sensVol == 0x0) {
      if (i == 0) fGeo = 1; // modular geometry
      break;
    }
    TArrayD* params = sensVol->getParameters();
    fRad[2*i] = params->At(0);
    fRad[2*i+1] = params->At(1);
    //dR = (rMax-rMin) / 50; 
    fDz[i] = params->At(2);
    Int_t nMods = Int_t (fDz[i] * 2. / size + 0.1);
    fZmod[i] = fDz[i] * 2. / nMods;
    cout << " *** STS sensitive volume: " << sensVol->GetName() << " " << params->At(0) 
	 << " " << fDz[i] << " " << fZmod[i] << endl;
  }
  //fStereoA[0] = -7.5;
  //fStereoA[1] = 7.5;
  fStereoA[0] = 7.5;
  fStereoA[1] = 0.0;
  //for (Int_t i = 0; i < 2; ++i) fStereoA[i] = 
  //TMath::Tan (fStereoA[i]*TMath::DegToRad()); 
  for (Int_t i = 0; i < 2; ++i) fStereoA[i] *= TMath::DegToRad();

  if (fGeo) {
    // Process modular geometry
    Double_t safety = 0.03;
    for (Int_t i = 0; i < nLay; ++i) {
      fNladders[2*i] = fNladders[2*i+1] = fNsectors[2*i] = fNsectors[2*i+1] = 0;
      volName = "sts01ladder";
      volName += (i+1);
      path = "/cave_1/sts01_0/" + volName;
      path += "_";
      TString path0 = path;
      //volName = "/cave_1/sts01_0/sts01ladder";
      //volName += (i+1);
      // Loop over all ladders to find the one with the smallest radius
      fRad[2*i+1] = fRad[2*i] = 999.9;
      for (Int_t jlad = 1; jlad < 99; ++jlad) {
	path = path0;
	path += jlad;
	if (!gGeoManager->CheckPath(path)) break;
	gGeoManager->cd(path);
	TGeoVolume *ladd = gGeoManager->GetVolume(volName);
	TGeoBBox *box = (TGeoBBox*) ladd->GetShape();
	Double_t xyzL[3] = {0}, xyzM[3];
	gGeoManager->LocalToMaster(xyzL,xyzM);
	Double_t rad = TMath::Sqrt (xyzM[0] * xyzM[0] + xyzM[1] * xyzM[1]);
	fRad[2*i] = TMath::Min (fRad[2*i],rad);
	xyzL[0] = box->GetDX();
	gGeoManager->LocalToMaster(xyzL,xyzM);
	rad = TMath::Sqrt (xyzM[0] * xyzM[0] + xyzM[1] * xyzM[1]);
	fRad[2*i] = TMath::Min (fRad[2*i],rad);
	xyzL[0] = -box->GetDX();
	gGeoManager->LocalToMaster(xyzL,xyzM);
	rad = TMath::Sqrt (xyzM[0] * xyzM[0] + xyzM[1] * xyzM[1]);
	fRad[2*i] = TMath::Min (fRad[2*i],rad);
      }
      fRad[2*i+1] = fRad[2*i];
      TGeoVolume *ladd = gGeoManager->GetVolume(volName);
      if (ladd == NULL) { nLay = i; break; }
      TGeoBBox *box = (TGeoBBox*) ladd->GetShape();
      //safety = -box->GetDY();
      //safety = box->GetDY();
      safety = 2 * box->GetDY(); // new
      fRad[2*i] -= safety;
      fRad[2*i+1] -= safety;
      //new if (i == 0) { fRad[2*i] -= safety; fRad[2*i+1] -= safety; }
    }
    FillGeoScheme();
  }

  // Get pipe radius
  fPipeR = ((TGeoTube*)gGeoManager->GetVolume("pipe1")->GetShape())->GetRmax();

  // Get cables
  TObjArray *vols = gGeoManager->GetListOfVolumes();
  Int_t nvols = vols->GetEntriesFast(), ncables = 0;
  for (Int_t i = 0; i < nvols; ++i) {
    TGeoVolume *vol = (TGeoVolume*) vols->UncheckedAt(i);
    TString cable = TString(vol->GetName());
    if (!(cable.Contains("sts") && cable.Contains("cable"))) continue;
    //cout << cable << endl;
    ++ncables;
    TGeoBBox *box = (TGeoBBox*) vol->GetShape();
    TString lad = cable;
    Int_t ip = lad.Index("cable");
    lad.Replace(ip,lad.Length()-ip+1,"");
    lad.Replace(lad.Length()-2,1,"");
    lad += "_1/";
    path = "/cave_1/sts01_0/" + lad + cable;
    path += "_1";
    gGeoManager->cd(path);
    Double_t xyzL[3] = {0}, xyzM[3];
    gGeoManager->LocalToMaster(xyzL,xyzM);
    //cout << xyzM[2] - box->GetDZ() << " " << xyzM[2] + box->GetDZ() << " " << box->GetDZ() << endl;
    if (xyzM[2] - box->GetDZ() > 0) {
      Int_t lay = TString(lad(lad.Length()-4,1)).Atoi();
      fCables[lay-1].insert(pair<Double_t,Double_t>(xyzM[2] - box->GetDZ(),xyzM[2] + box->GetDZ()));
    }
  }

  return kSUCCESS;
}

//__________________________________________________________________________
void MpdTrackFinderIts::FillGeoScheme()
{
  /// Fill Kalman filter geometry manager info

  MpdKalmanGeoScheme *geo = MpdKalmanFilter::Instance()->GetGeo();

  TGeoVolume *vSts = gGeoManager->GetVolume("sts01");

  for (Int_t layer = 1; layer < 999; ++layer) {
    // Loop over layers
    TString sladder = "sts01ladder";
    sladder += layer;
    TGeoVolume *vLay = gGeoManager->GetVolume(sladder);
    if (vLay == 0x0 && layer == 1) continue; // no first layer
    if (vLay == 0x0) break;
    sladder += "_";
    //TString sdet = "sts01detector";
    TString sdet = "sts01sensor";
    sdet += layer;
    sdet += "_";
    //TString sensor = "sts01sensor";
    TString sensor = "";
    //sensor += layer;
    //sensor += "_0";

    for (Int_t ladder = 1; ladder < 999; ++ladder) {
      // Loop over ladders
      TString sladder1 = sladder;
      sladder1 += ladder;
      TGeoNode *node = vSts->FindNode(sladder1);
      if (node == 0x0) break;
      ++fNladders[2*(layer-1)];
      ++fNladders[2*(layer-1)+1];
      TGeoVolume *vLad = node->GetVolume();
      //cout << vLad->GetNodes()->GetEntriesFast() << " " << vLad->GetNdaughters() << endl;
      Int_t nDaught = vLad->GetNdaughters(), detID = -1, detIDsts = -1;
      TObjArray *daught = vLad->GetNodes();
      //for (Int_t j = 0; j < vLad->GetNdaughters(); ++j) 
      //cout << ((TGeoNode*)(vLad->GetNodes()->UncheckedAt(j)))->GetName() << endl;

      Int_t iZ = 0;
      for (Int_t det = 0; det < nDaught; ++det) {
	// Loop over ladder daughters
	TString sdet1 = ((TGeoNode*)(daught->UncheckedAt(det)))->GetName();
	if (!sdet1.Contains("sensor") && !sdet1.Contains("sector")) continue;
	if (ladder == 1) { ++fNsectors[2*(layer-1)]; ++fNsectors[2*(layer-1)+1]; }
	Int_t det1 = TString(sdet1(sdet1.Index("_")+1,2)).Atoi();
	Int_t secType = -1;
	if (sdet1.Contains("sector")) secType = TString(sdet1(sdet1.Index("_")-2,1)).Atoi();

	++iZ;
	for (Int_t side = 0; side < 2; ++side) {
	  
	  detIDsts = fHitSts.SetDetId(secType, layer, ladder, det1, side);
	  detID = fHitSts.SetDetId(0, layer, ladder, iZ, side);
	  //detIDsts += 1000000 * ((layer-1) * 2 + side);
	  fId2Id[(layer-1) * 2 + side].insert(pair<Int_t,Int_t>(detIDsts,detID));
	  detID += 1000000 * ((layer-1) * 2 + side);

	  TString detName = sladder1 + "/" + sdet1 + "#";
	  detName += side;
	  //if (side) geo->SetDetId(detName, detID&((2<<26)-2)); // store odd side
	  //if (!side) geo->SetDetId(detName, detID); // store even side
	  geo->SetDetId(detName, detID);

	  //TString path = "/cave_1/sts01_0/" + detName + "/" + sensor;
	  TString path = "/cave_1/sts01_0/" + detName(0,detName.Length()-2);
	  //cout << detName << " " << path << endl;
	  gGeoManager->cd(path);
	  //if (side) geo->SetPath(detID&((2<<26)-2), path); // store odd side
	  //if (!side) geo->SetPath(detID, path); // store even side
	  geo->SetPath(detID, path);
          node = gGeoManager->GetCurrentNode();
	  //cout << node->GetName() << " " << detID << endl;
	  
	  TGeoVolume *vol = node->GetVolume();
	  TGeoBBox *box = (TGeoBBox*) vol->GetShape();
	  TVector2 size(box->GetDX(), box->GetDZ());
	  geo->SetSize(detID, size);
	  
	  Double_t xyzL[3] = {0}, xyzM[3], vecM[3];
	  xyzL[1] = 1;
	  gGeoManager->LocalToMasterVect(xyzL,vecM);
	  //cout << vecM[0] << " " << vecM[1] << " " << vecM[2] << endl;
	  //TVector3 norm(vecM[0], vecM[1], vecM[2]);
	  //geo->SetNormal(detID, norm);
	  
	  xyzL[1] = 0;
	  gGeoManager->LocalToMaster(xyzL,xyzM);
	  Double_t sign = TMath::Sign (1.,xyzM[0]*vecM[0]+xyzM[1]*vecM[1]);
	  if (detID % 2) xyzL[1] = 0.015 * sign;
	  else xyzL[1] = -0.015 * sign;
	  gGeoManager->LocalToMaster(xyzL,xyzM);
	  //cout << xyzM[0] << " " << xyzM[1] << " " << xyzM[2] << endl;
	  TVector3 pos(xyzM[0], xyzM[1], xyzM[2]);
	  geo->SetGlobalPos(detID, pos);

	  //xyzL[1] = sign;
	  //gGeoManager->LocalToMasterVect(xyzL,vecM);
	  //cout << vecM[0] << " " << vecM[1] << " " << vecM[2] << endl;
	  TVector3 norm(sign*vecM[0], sign*vecM[1], sign*vecM[2]);
	  geo->SetNormal(detID, norm);
	}
      }
    }
  }

}

//__________________________________________________________________________
InitStatus MpdTrackFinderIts::ReInit()
{
  fItsPoints = (TClonesArray *) FairRootManager::Instance()->GetObject("StsPoint");
  fItsHits =(TClonesArray *) FairRootManager::Instance()->GetObject("StsHit");
  if (fItsPoints == 0x0 || fItsHits == 0x0) return kERROR;
  fTpcTracks =(TClonesArray *) FairRootManager::Instance()->GetObject("TpcKalmanTrack");
  fEctTracks =(TClonesArray *) FairRootManager::Instance()->GetObject("EctTrack");
  fMCTracks =(TClonesArray *) FairRootManager::Instance()->GetObject("MCTrack");
  //fSTSTrackMatch = (TClonesArray*) FairRootManager::Instance()->GetObject("STSTrackMatch");
  //fPrimVtx =  (FairVertex *) FairRootManager::Instance() ->GetObject("PrimaryVertex");

  FairRootManager::Instance()->Register("ItsTrack", "Its", fTracks, kTRUE);

  fNPass = 1;

  return kSUCCESS;
}

//__________________________________________________________________________
void MpdTrackFinderIts::Reset() 
{
  ///
  cout << " MpdTrackFinderIts::Reset  " << endl;

  //fKHits->Clear();
  fKHits->Delete();
  fTracks->Delete();
  delete [] fLayPointers;
  fLayPointers = NULL;
}

//__________________________________________________________________________
void MpdTrackFinderIts::SetParContainers()
{
  // Get run and runtime database
  FairRunAna* run = FairRunAna::Instance();
  if ( ! run ) Fatal("SetParContainers", "No analysis run");

  FairRuntimeDb* db = run->GetRuntimeDb();
  if ( ! db ) Fatal("SetParContainers", "No runtime database");

  // Get STS geometry parameter container
  db->getContainer("MpdStsGeoPar");
}

//__________________________________________________________________________
void MpdTrackFinderIts::Finish()
{
  //Write();
}

//__________________________________________________________________________
void MpdTrackFinderIts::Exec(Option_t * option)
{

  static int eventCounter = 0;    
  cout << " - - - - \n ItsRec event " << ++eventCounter << endl;

  Reset();

  // Create Kalman hits
  if (fItsHits->GetEntriesFast() == 0) return;
  MakeKalmanHits();

  fTracks->Clear();
  for (Int_t i = 0; i < fNPass; ++i) {
    //fTracks->Clear();
    GetTrackSeeds(0); // TPC tracks
    GetTrackSeeds(1); // ECT tracks

    cout << "  Total number of hits for tracking: " << fKHits->GetEntriesFast() << endl;
    cout << "  Total number of track seeds: " << fTracks->GetEntriesFast() << endl;

    DoTracking(i);
    //StoreTracks();
    cout << "  Total number of found tracks: " << fTracks->GetEntriesFast() << endl;
    if (i != fNPass - 1) ExcludeHits(); // exclude used hits
  }
  AddHits(); // add hit objects to tracks
  cout << "  Total number of found tracks: " << fTracks->GetEntriesFast() << endl;
}

//__________________________________________________________________________
void MpdTrackFinderIts::MakeKalmanHits()
{
  /// Create Kalman hits from ITS hits.

  fhLays->Reset();
  Int_t nHits = fItsHits->GetEntriesFast(), layMax = 0, lay = 0, nKH = 0;
  Double_t r, z, xloc, errZ = 0.012, errX = 0.0023; // 120um in Z, 23um in R-Phi (local X)
  //Double_t r, z, xloc, errZ = 0.12, errX = 0.0023; // 1.2mm in Z, 23um in R-Phi (local X)
  //Double_t r, z, xloc, errZ = 50.0, errX = 0.01; // 50cm in Z, 100um in R-Phi (local X)

  for (Int_t ih = 0; ih < nHits; ++ih) {
    MpdStsHit *h = (MpdStsHit*) fItsHits->UncheckedAt(ih);
    r = TMath::Sqrt (h->GetX() * h->GetX() + h->GetY() * h->GetY());
    z = h->GetZ();
    xloc = h->GetLocalX();
    //cout << ih << " " << h->Layer()-1 << endl;
    //lay = h->Layer() * 2 + h->Side() + 1; 
    lay = (h->Layer() - 1) * 2 + h->Side(); 

    // Add error                                            
    Double_t dX = 0, dZ = 0;
    gRandom->Rannor(dX,dZ);
    //if (errZ > 2) dZ = 0.0; // 1-D case
    dZ = 0.0; // 1-D case
    Double_t meas[2] = {xloc+dX*errX, z+dZ*errZ};
    Double_t err[2] = {errX, errZ};
    Double_t cossin[2] = {TMath::Cos(fStereoA[h->Side()]), TMath::Sin(fStereoA[h->Side()])};
    //(Int_t detID, Int_t nDim, HitType hitType, Double_t *meas, Double_t *err, Double_t *cosSin, Double_t signal, Double_t dist, Int_t index)
    //MpdKalmanStripHit *hit = new ((*fKHits)[nKH++]) MpdKalmanStripHit(r, fStereoA[h->Side()],
    //			  xloc+dX*errX, z+dZ*errZ, errX, errZ, lay, ih); 
    //hit->SetType(MpdKalmanHit::kFixedR);
    MpdKalmanHit *hit = 0x0;
    //cout << h->GetDetectorID() << " " << fId2Id[lay][h->GetDetectorID()] << endl;
    //if (fGeo && h->GetUniqueID()) hit = new ((*fKHits)[nKH++]) MpdKalmanHit(lay*1000000+h->GetDetectorID(), 1, 
    if (fGeo && h->GetUniqueID()) hit = 
				    new ((*fKHits)[nKH++]) MpdKalmanHit(lay*1000000+fId2Id[lay][h->GetDetectorID()], 1, 
					MpdKalmanHit::kFixedP, meas, err, cossin, 0., r, ih);
    // Mask out sector number - sensor layout
    else if (fGeo) hit = new ((*fKHits)[nKH++]) MpdKalmanHit(lay*1000000+(h->GetDetectorID()&((2<<12)-1)), 1, 
                                                MpdKalmanHit::kFixedP, meas, err, cossin, 0., r, ih);
    else hit = new ((*fKHits)[nKH++]) MpdKalmanHit(lay*1000000+nKH-1, 1, MpdKalmanHit::kFixedR, 
						   meas, err, cossin, 0., r, ih);
    hit->SetUniqueID(0);
    // Add second measurement - just for test at the moment
    //!!!
    //hit->SetNofDim(2);
    //!!!
    layMax = TMath::Max (lay, layMax);
    fhLays->Fill(lay+0.1);
  }
  cout << " Max layer = " << layMax << ", hits: " << fKHits->GetEntriesFast() << endl;
  fKHits->Sort(); // in descending order in R
  //cout << ((MpdKalmanHit*)fKHits->UncheckedAt(0))->GetPos() << endl;
  cout << ((MpdKalmanHit*)fKHits->UncheckedAt(0))->GetDist() << endl;

  fLayPointers = new Int_t [layMax+1];
  Int_t ipos = 0;
  for (Int_t i = layMax; i > -1; --i) {
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
}

//__________________________________________________________________________
void MpdTrackFinderIts::GetTrackSeeds(Int_t iPass)
{
  /// Build ITS track seeds from TPC tracks

  TClonesArray *seeds = (iPass == 0) ? fTpcTracks : fEctTracks;
  if (seeds == NULL) return;
  Int_t nSeeds = seeds->GetEntriesFast();
  cout << " Seed tracks: " << nSeeds << endl;
 
  Int_t nCand = fTracks->GetEntriesFast();
  MpdKalmanHit hit;
  hit.SetType(MpdKalmanHit::kFixedR);
  for (Int_t itr = 0; itr < nSeeds; ++itr) {
    MpdTpcKalmanTrack *tpc = (MpdTpcKalmanTrack*) seeds->UncheckedAt(itr);
    if (tpc->GetType() != MpdKalmanTrack::kBarrel) continue;
    //tpc->GetParam()->Print();
    //MpdEctKalmanTrack *track = new ((*fTracks)[nCand++]) MpdEctKalmanTrack(itr, *tpc);
    MpdItsKalmanTrack *track = NULL;
    if (TString(tpc->ClassName()).Contains("Tpc")) track = new ((*fTracks)[nCand++]) MpdItsKalmanTrack(*tpc);
    else track = new ((*fTracks)[nCand++]) MpdItsKalmanTrack(*((MpdEctKalmanTrack*)tpc));
    track->SetUniqueID(itr+1);
    /*
    track->SetParam(*track->GetParamAtHit());
    track->SetParamNew(*track->GetParamAtHit());
    track->SetPosNew(track->GetPos());
    track->SetWeight(*track->GetWeightAtHit());
    track->SetLength(track->GetLengAtHit());
    */
    hit.SetPos(track->GetPos()); 
    track->SetParamNew(*track->GetParam());
    track->SetPos(track->GetPosNew());
    track->ReSetWeight();
    TMatrixDSym w = *track->GetWeight(); // save current weight matrix
    MpdKalmanFilter::Instance()->PropagateToHit(track,&hit,kFALSE);
    track->SetWeight(w); // restore original weight matrix (near TPC inner shell)
    //cout << nCand-1 << " " << track->GetTrackID() << endl;
    //cout << track->GetHits()->GetEntriesFast() << " " << track->GetTrHits()->GetEntriesFast() << endl;
    track->GetHits()->Clear();
    track->SetChi2Its(track->GetChi2()); // temporary storage
    track->SetChi2(0.);
  }
  cout << " Number of ITS track candidates: " << nCand << endl;
}  
  
//__________________________________________________________________________
void MpdTrackFinderIts::DoTracking(Int_t iPass)
{
  /// Run Kalman tracking
  
  Double_t vert[3] = {0.0,0.0,0.0};
  Int_t nCand = fTracks->GetEntriesFast(), iok = 0;
  Int_t lay0 = ((MpdKalmanHit*)fKHits->First())->GetLayer();
 
  for (Int_t i = 0; i < nCand; ++i) {
    MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTracks->UncheckedAt(i);
    //cout << " Track seed No. " << i << ", ID: " << track->GetTrackID() << ", Hits: " << track->GetNofTrHits() << endl;
    for (Int_t j = 0; j < track->GetNofTrHits(); ++j) {
      MpdKalmanHit *h = (MpdKalmanHit* )track->GetTrHits()->UncheckedAt(j);
      //MpdStsHit *hh = (MpdStsHit*) fItsHits->UncheckedAt(h->GetIndex());
      //Int_t id = ((FairMCPoint*) fItsPoints->UncheckedAt(hh->GetRefIndex()))->GetTrackID();
      //cout << j << " " << h->GetDist() << " " << h->GetLayer() << endl;
    }
    
    // Reset weight matrix - for debug
    /*
    Int_t nHits = track->GetNofTrHits();
    nHits *= nHits;
    for (Int_t ii = 0; ii < 5; ++ii) {
      for (Int_t j = i; j < 5; ++j) {
	if (j == ii) (*track->GetWeight())(ii,j) /= nHits;
	else (*track->GetWeight())(ii,j) = (*track->GetWeight())(j,ii) = 0.;
      }
    }
    */

    if (fGeo) iok = RunKalmanFilterMod(track, lay0); // modular geometry
    else iok = RunKalmanFilterCyl(track, lay0); // cylindrical geometry
    if (iok == -1) {
      fTracks->RemoveAt(i);
      continue;
    }

    //if (track->GetNofHits() == 0) continue; // no hits added
    // Propagate track to the beam line
    track->SetParam(*track->GetParamNew());
    track->SetPos(track->GetPosNew());
    Double_t pos = track->GetPos();
    TMatrixD par = *track->GetParam();
    TMatrixDSym cov = *track->Weight2Cov();
    Double_t leng = track->GetLength();
    TString nodeNew = track->GetNodeNew();
    //cout << " 1: " << nodeNew << ", " << track->GetNode() << endl;

    // Go to beam pipe
    MpdKalmanHit hit;
    hit.SetType(MpdKalmanHit::kFixedR);
    hit.SetPos(fPipeR);
    iok = MpdKalmanFilter::Instance()->PropagateToHit(track, &hit, kTRUE);
    if (iok != 1) {
      // Restore track
      track->SetParam(par);
      track->SetParamNew(par);
      track->SetCovariance(cov);
      track->ReSetWeight();
      track->SetPos(pos);
      track->SetPosNew(pos);
      track->SetLength(leng);
      //track->SetNode(node);
      //cout << " 2: " << nodeNew << ", " << track->GetNode() << endl;
      track->SetNodeNew(nodeNew);
    } else {
      // Add multiple scattering
      //Double_t dX = 0.05 / 8.9; // 0.5 mm of Al
      Double_t dX = 0.1 / 35.28; // 1. mm of Be
      TMatrixDSym* pcov = track->Weight2Cov();
      Double_t th = track->GetParamNew(3);
      Double_t cosTh = TMath::Cos(th);
      Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, dX);
      (*pcov)(2,2) += (angle2 / cosTh / cosTh);
      (*pcov)(3,3) += angle2;
      Int_t ok = 0;
      MpdKalmanFilter::Instance()->MnvertLocal(pcov->GetMatrixArray(), 5, 5, 5, ok);
      track->SetWeight(*pcov);
    }
    cov = *track->Weight2Cov();

    hit.SetPos(0.);
    hit.SetMeas(0,track->GetParam(2)); // track Phi
    //cout << i << " " << track->GetTrackID() << " " << track->GetLength() << " " << ((MpdKalmanHitR*)track->GetHits()->First())->GetLength() << endl;
    //Double_t pos = ((MpdKalmanHit*)track->GetHits()->Last())->GetPos();
    //MpdKalmanFilter::Instance()->PropagateParamR(track, &hit, kTRUE);
    iok = MpdKalmanFilter::Instance()->PropagateToHit(track, &hit, kTRUE);
    if (iok != 1) MpdKalmanFilter::Instance()->FindPca(track, vert);
    //track->SetPos(pos); // restore position
    track->SetParam(*track->GetParamNew()); // !!! track params at PCA
    //track->GetCovariance()->Print();
  } // for (Int_t i = 0; i < nCand;
  fTracks->Compress();
}
    
//__________________________________________________________________________
Int_t MpdTrackFinderIts::RunKalmanFilterCyl(MpdItsKalmanTrack *track, Int_t layBeg)
{
  /// Run Kalman filter in the cylindrical geometry (might not work when propagating outward!!!)

  const Int_t maxBr = 10, maxBrLay = 1000; // max number of branches

  Int_t layMax = ((MpdKalmanHit*)fKHits->First())->GetLayer();
  MpdKalmanHit *hitOK = 0x0;
  MpdKalmanHit hitTmp;
  MpdKalmanTrack::TrackDir trackDir = track->GetDirection();
  //Int_t layBeg = 0, layEnd = -1, dLay = -1, layOK = -1;
  Int_t layEnd = -1, dLay = -1, layOK = -1;
  if (trackDir == MpdKalmanTrack::kOutward) {
    layEnd = layMax + 1;
    dLay = 1;
  }
  
  TMatrixDSym pointWeight(5), pointWeightTmp(5), saveWeight(5);
  TMatrixD param(5,1), paramTmp(5,1);

  Double_t saveZ = 0.0, saveLeng = 0.0, dChi2Min = 0.0, quality[maxBrLay] = {0.0};
  //cout << " Starting hit: " << hitOK->GetLayer() << " " << hitOK->GetTrackID() << " " << hitOK->GetUsage() << endl;
  
  MpdItsKalmanTrack branch[maxBr];
  Int_t nBr = 1, nBr0 = 0, indx_lay[maxBrLay] = {0}, maxInLay = 0, ok = 0, nAdd = 0;

  branch[0] = *track;

  for (Int_t lay = layBeg; lay != layEnd; lay+=dLay) {
    Int_t nLay = GetNofHitsInLayer(lay);
    Int_t indx0 = GetHitsInLayer(lay);
    MpdKalmanHit *hitMin = 0x0;
    //cout << " lay, nLay: " << lay << " " << nLay << " " << indx0 << endl;
    Int_t indxBeg = 0, indxEnd = nLay, dIndx = 1;
    MpdItsKalmanTrack branchLay[maxBrLay];

    nBr0 = nBr;
    nBr = 0;
    for (Int_t ibr = 0; ibr < nBr0; ++ibr) {
      MpdItsKalmanTrack *curTr = &branch[ibr];
      
      // Check for missing layer (2 sides)
      Int_t lastLay = 8, frozen = 0;
      MpdKalmanHit *h = (MpdKalmanHit*) curTr->GetHits()->Last();
      if (h) lastLay = h->GetLayer();
      if (TMath::Abs(lay - lastLay) > 2) frozen = 1;

      ok = nAdd = 0;
      for (Int_t indx = indxBeg; indx != indxEnd; indx+=dIndx) {
	if (frozen) break;

	MpdKalmanHit *hit = (MpdKalmanHit*) fKHits->UncheckedAt(indx0+indx);
	Double_t posNew = 0;

	// Propagate to hit (if it is not very close to the track)
	if (TMath::Abs(hit->GetPos()-curTr->GetPosNew()) > 1.e-4) { 
	  Double_t leng = curTr->GetLength();
	  posNew = curTr->GetPosNew();
	  TMatrixD parNew = *curTr->GetParamNew();
	  TString nodeNew = curTr->GetNodeNew();
	  TString curPath = curTr->GetNode();
	  if (!MpdKalmanFilter::Instance()->PropagateToHit(curTr,hit,kTRUE,kFALSE)) { 
	    // Restore initial parameters for the failed track
	    curTr->SetPosNew(posNew);
	    curTr->SetParamNew(parNew);
	    curTr->SetLength(leng);
	    curTr->SetNodeNew(nodeNew);
	    curTr->SetNode(curPath);
	    ok = -1; 
	    break; 
	  } 

	  Double_t step = curTr->GetLength() - leng;
	  //*
	  //if (hit->GetLayer() % 2 == 10 && step > 1.e-4) {
	  if (hit->GetLayer() % 2 == 0 && step > 1.e-4) {
	    // Crossing silicon layer - add mult. scat. in the sensor
	    Double_t x0 = 9.36; // rad. length
	    TMatrixDSym *cov = curTr->Weight2Cov();
	    Double_t th = curTr->GetParamNew(3);
	    Double_t cosTh = TMath::Cos(th);
	    Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(curTr, x0, step);
	    //cout << " Scat: " << hit->GetLayer() << " " << step << " " << TMath::Sqrt(angle2) << endl;
	    (*cov)(2,2) += (angle2 / cosTh / cosTh);
	    (*cov)(3,3) += angle2;
	    Int_t iok = 0;
	    MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
	    curTr->SetWeight(*cov);
	  } else if (hit->GetLayer() % 2 != 0 && step > 1.e-4 && fCables[lay/2].size() > 0) {
	    // Crossing silicon layer - add mult. scat. in the cable
	    Double_t nCables = 0, x0 = 0.0116 * 2 / 9.36; // in rad. length - 116um cable per side
	    // Find number of cables crossed
	    TString path = gGeoManager->GetPath();
	    if (!path.Contains("sensor") && !path.Contains("sector")) {
	      cout << " !!! MpdTrackFinderIts::RunKalmanFilter - Outside detector !!! " << endl;
	      exit(0);
	    }
	    Double_t v7[3] = {curTr->GetParamNew(0), curTr->GetPosNew(), curTr->GetParamNew(1)}, v77[3];
	    gGeoManager->LocalToMaster(v7,v77);
	    Double_t zTr = TMath::Abs (v77[2]); // global Z
	    //cout << zTr << endl;
	    map<Double_t,Double_t>::iterator it;
	    for (it = fCables[lay/2].begin(); it != fCables[lay/2].end(); ++it) {
	      if (zTr < it->first || zTr > it->second) continue;
	      ++nCables;
	    }
	    //cout << " Cables: " << nCables << endl;
	    if (nCables) {
	      x0 *= nCables;
	      TMatrixDSym *cov = curTr->Weight2Cov();
	      Double_t th = curTr->GetParamNew(3);
	      Double_t cosTh = TMath::Cos(th);
	      Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(curTr, x0);
	      //cout << " Scat: " << hit->GetLayer() << " " << step << " " << TMath::Sqrt(angle2) << endl;
	      (*cov)(2,2) += (angle2 / cosTh / cosTh);
	      (*cov)(3,3) += angle2;
	      Int_t iok = 0;
	      MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
	      curTr->SetWeight(*cov);
	    }
	  }
	  //curTr->GetWeight()->Print();
	}
 
	//cout << hit->GetDetectorID() << endl;
	// Exclude used hits
	if (hit->GetFlag() != 1) continue;
	// !!! Exact ID match
	if (TrackID(hit) != track->GetTrackID()) continue;

	Double_t dChi2 = 0;
	// Simple (old) geometry
	TVector2 dist = GetDistance(curTr, hit);
	//if (track->GetTrackID() != -1335) cout << /*((TpcLheKalmanTrack*)fTpcTracks->UncheckedAt(track->GetTpcIndex()))->GetNofTrHits() <<*/ " " << dist.X() << " " << dist.Y() << " " << hit->GetLayer() << endl;
	if (dist.X() > 15.) continue; // distance in Z
	if (dist.Y() > 0.25) continue; // distance in Phi 

	hitTmp = *hit;
	hitTmp.SetMeas(0,MpdKalmanFilter::Instance()->Proxim(curTr,hit,hit->GetCosSin(0)));
	/*cout <<  hit->GetPos() << " " << hitTmp.GetMeas(0) << " " << hitTmp.GetMeas(1) << endl;
	cout << TrackID(hit) << " " << track->GetTrackID() << " " << curTr->GetPosNew() << " " 
	     << curTr->GetParamNew(0) << " " << track->GetParamNew(1) << endl; */
	dChi2 = MpdKalmanFilter::Instance()->FilterStrip(curTr,&hitTmp,pointWeight,param);
	  
	MpdItsKalmanTrack trTmp = *curTr;
	trTmp.SetParamNew(param);
	dist = GetDistance(&trTmp, hit); 
	//cout << " Chi2: " << dChi2 << " " << dist.X() << " " << dist.Y() << " " << curTr->GetParamNew(0) << " " << hit->GetLayer() << endl;
	Double_t dChi2Z = dist.X() / hit->GetErr(1); // contrib. due to Z-residual
	dChi2 += dChi2Z * dChi2Z;
	//cout << dChi2 << " " << TrackID(hit) << endl;
	if (dChi2 < 0) cout << " !!! Negative Chi2: !!! " << dChi2 << endl;

	if (TMath::Abs(dChi2) < fgkChi2Cut) {
	//if (dChi2 < fgkChi2Cut && lay % 2 == 0) {
	  if (nBr < maxBrLay) {
	    branchLay[nBr] = *curTr;
	    branchLay[nBr].GetHits()->Add(hit);
	    branchLay[nBr].SetChi2(curTr->GetChi2()+dChi2);
	    TMatrixDSym w = *curTr->GetWeight();
	    w += pointWeight;
	    branchLay[nBr].SetWeight(w);
	    branchLay[nBr].SetPosNew(curTr->GetPosNew());
	    //if (hit->GetType() != MpdKalmanHit::kFixedP) branchLay[nBr].SetPosNew(curTr->GetPosNew());
	    //else branchLay[nBr].SetPosNew(posNew); // modular geometry
	    branchLay[nBr].SetLength(curTr->GetLength());
	    branchLay[nBr].SetParamNew(param);
	    // Save track params at last hit
	    branchLay[nBr].SetLengAtHit(curTr->GetLength()); 
	    branchLay[nBr].SetParamAtHit(param);
	    branchLay[nBr].SetWeightAtHit(*branchLay[nBr].GetWeight());
	    ++nBr;
	    ++nAdd;
	  } 
	}
      } // for (Int_t indx = indxBeg; indx != indxEnd;

      // Add branch with missing layer
      if (nAdd == 0 && nBr < maxBrLay && ok != -1) {
	branchLay[nBr] = *curTr;
	++nBr;
      }

    } // for (Int_t ibr = 0; ibr < nBr0;
    maxInLay = TMath::Max (maxInLay,nBr);
    if (nBr == 0) {
      if (ok == -1) { 
	//cout << track->GetTrackID() << " " << lay << endl;
	//return ok; // stop track
	cout << branch[0].GetNode() << ", " << branch[0].GetNodeNew() << endl;
	break;
      } 
      cout << " !!! MpdTrackFinderIts::RunKalmanFilter: it can not happen |" << endl;
      exit(0);
    } else if (nBr <= maxBr) {
      for (Int_t i = 0; i < nBr; ++i) {
	branch[i].Reset();
	branch[i] = branchLay[i];
      }
    } else {
      // Too many branches - take the best ones
      for (Int_t i = 0; i < nBr; ++i) {
	Double_t c2 = TMath::Min (branchLay[i].GetChi2(),200.);
	quality[i] = branchLay[i].GetNofHits() + (0.999-c2/201);
      }
      TMath::Sort(nBr,quality,indx_lay);
      for (Int_t i = 0; i < maxBr; ++i) {
	branch[i].Reset();
        branch[i] = branchLay[indx_lay[i]];
      }
      nBr = maxBr;
    }

  } // for (Int_t lay = layOK-1; lay >= 0;

  // Select the best branch
  //cout << " Branches: " << nBr << " " << maxInLay << endl;
  Int_t ibest = 0;
  for (Int_t i = 1; i < nBr; ++i) {
    if (branch[i].GetNofHits() > branch[ibest].GetNofHits()) ibest = i;
    else if (branch[i].GetNofHits() == branch[ibest].GetNofHits() &&
	     branch[i].GetChi2() < branch[ibest].GetChi2()) ibest = i;
  }
  track->Reset();
  *track = branch[ibest];
  return 0;
}

//__________________________________________________________________________
Int_t MpdTrackFinderIts::RunKalmanFilterMod(MpdItsKalmanTrack *track, Int_t layBeg)
{
  /// Run Kalman filter in the modular geometry (might not work when propagating outward!!!)

  const Int_t maxBr = 20, maxBrLay = 1000; // max number of branches
  //const Int_t maxBr = 100, maxBrLay = 1000; // max number of branches

  //cout << fHits->GetEntriesFast() << endl;
  //Int_t layMax = ((MpdKalmanHit*)fKHits->Last())->GetLayer();
  Int_t layMax = ((MpdKalmanHit*)fKHits->First())->GetLayer();
  MpdKalmanHit *hitOK = 0x0;
  MpdKalmanHit hitTmp;
  MpdKalmanTrack::TrackDir trackDir = track->GetDirection();
  //Int_t layBeg = 0, layEnd = -1, dLay = -1, layOK = -1;
  Int_t layEnd = -1, dLay = -1, layOK = -1;
  if (trackDir == MpdKalmanTrack::kOutward) {
    layEnd = layMax + 1;
    dLay = 1;
  }
  
  TMatrixDSym pointWeight(5), pointWeightTmp(5), saveWeight(5);
  TMatrixD param(5,1), paramTmp(5,1);

  Double_t saveZ = 0.0, saveLeng = 0.0, dChi2Min = 0.0, quality[maxBrLay] = {0.0}, posNew = 0.0;
  //cout << " Starting hit: " << hitOK->GetLayer() << " " << hitOK->GetTrackID() << " " << hitOK->GetUsage() << endl;
  
  MpdItsKalmanTrack branch[maxBr];
  Int_t nBr = 1, nBr0 = 0, indx_lay[maxBrLay] = {0}, maxInLay = 0, ok = 0, nAdd = 0;

  branch[0] = *track;

  TString mass2 = "0.0194797849"; // pion mass squared
  if (fMCTracks) {
    // Get particle mass - ideal PID
    MpdMCTrack *mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(track->GetTrackID());
    TParticlePDG *pdgP = TDatabasePDG::Instance()->GetParticle(mctrack->GetPdgCode());
    if (pdgP) {
      Double_t mass = pdgP->Mass();
      if (mass < 0.1 || mass > 0.8) {
	// Electrons or protons (or heavier than protons)
	mass2 = "";
	mass2 += mass*mass;
      }
    }
  }

  for (Int_t lay = layBeg; lay != layEnd; lay+=dLay) {
    Int_t nLay = GetNofHitsInLayer(lay);
    Int_t indx0 = GetHitsInLayer(lay);
    MpdKalmanHit *hitMin = 0x0;
    //cout << " lay, nLay: " << lay << " " << nLay << " " << indx0 << endl;
    Int_t indxBeg = 0, indxEnd = nLay, dIndx = 1;
    MpdItsKalmanTrack branchLay[maxBrLay];

    nBr0 = nBr;
    nBr = 0;
    for (Int_t ibr = 0; ibr < nBr0; ++ibr) {
      MpdItsKalmanTrack *curTr = &branch[ibr], *branchTr = NULL;
      
      // Check for missing layer (2 sides)
      Int_t lastLay = 8, frozen = 0;
      //Int_t lastLay = 10, frozen = 0;
      MpdKalmanHit *h = (MpdKalmanHit*) curTr->GetHits()->Last();
      if (h) lastLay = h->GetLayer();
      if (TMath::Abs(lay - lastLay) > 2) frozen = 1;

      ok = nAdd = 0;

      Int_t firstHit = -1, skipHit = -1;
      MpdItsKalmanTrack trackBr[3] = {*curTr, *curTr, *curTr};
      map<Int_t,Int_t> trackBrM;
      Bool_t navOK = NavigateToLayer(lay, curTr, trackBr, trackBrM); // navigate to layer
      if (!navOK) { ok = -1; continue; }
      // Debug
      //cout << " Overlaps: " << trackBrM.size() << endl;
      map<Int_t,Int_t>::iterator it;
      /*
      for (it = trackBrM.begin(); it != trackBrM.end(); ++it) {
	cout << " detID " << it->first << " " << fHitSts.GetLadder(it->first%1000000) << " " 
	     << fHitSts.GetSensor(it->first%1000000) << " " << it->second << endl;
      }
      */

      for (Int_t indx = indxBeg; indx != indxEnd; indx+=dIndx) {
	//if (frozen) break;

	MpdKalmanHit *hit = (MpdKalmanHit*) fKHits->UncheckedAt(indx0+indx);

	if (hit->GetDetectorID() != firstHit) {
	  firstHit = hit->GetDetectorID();
	  if (trackBrM.find(hit->GetDetectorID()) == trackBrM.end()) {
	    // Track and hit in different detectors - can happen near sensor edges 
	    // (hits in the same layer but in different detectors)
	    skipHit = firstHit;
	    continue; // next hit
	  }

	  Double_t leng = curTr->GetLength();
	  branchTr = &trackBr[trackBrM[hit->GetDetectorID()]];
	  //cout << branchTr->GetNodeNew() << endl;
	  Double_t step = branchTr->GetLength() - leng;
	  //cout << " Step " << step << " " << hit->GetLayer() << endl;
	  //*
	  //if (hit->GetLayer() % 2 == 10 && step > 1.e-4) {
	  if (hit->GetLayer() % 2 == 0 && step > 1.e-4) {
	    // Crossing silicon layer - add mult. scat. in the sensor
	    Double_t x0 = 9.36; // rad. length
	    TMatrixDSym *cov = branchTr->Weight2Cov();
	    Double_t th = branchTr->GetParamNew(3);
	    Double_t cosTh = TMath::Cos(th);
	    //Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(branchTr, x0, step, mass2);
	    Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(branchTr, x0, step*4, mass2);
	    //cout << " Scat: " << hit->GetLayer() << " " << step << " " << TMath::Sqrt(angle2) << endl;
	    (*cov)(2,2) += (angle2 / cosTh / cosTh);
	    (*cov)(3,3) += angle2;
	    Int_t iok = 0;
	    MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
	    branchTr->SetWeight(*cov);
	  } else if (hit->GetLayer() % 2 != 0 && step > 1.e-4 && fCables[lay/2].size() > 0) {
	    /*
	    // Crossing silicon layer - add mult. scat. in the cable
	    //Double_t nCables = 0, x0 = 0.0116 * 2 / 9.36; // in rad. length - 116um cable per side
	    Double_t nCables = 0, x0 = 0.0116 * 2; // in rad. length - 116um cable per side
	    // Find number of cables crossed
	    TString path = branchTr->GetNodeNew();
	    if (!path.Contains("sensor") && !path.Contains("sector")) {
	      cout << " !!! MpdTrackFinderIts::RunKalmanFilter - Outside detector !!! " << endl;
	      exit(0);
	    }
	    gGeoManager->cd(path);
	    Double_t v7[3] = {branchTr->GetParamNew(0), branchTr->GetPosNew(), branchTr->GetParamNew(1)}, v77[3];
	    gGeoManager->LocalToMaster(v7,v77);
	    Double_t zTr = TMath::Abs (v77[2]); // global Z
	    cout <<  " zTr " << v77[0] << " " << v77[1] << " " << v77[2] << endl;
	    cout << gGeoManager->GetPath() << " " << branchTr->GetNodeNew() << endl;
	    map<Double_t,Double_t>::iterator it1;
	    for (it1 = fCables[lay/2].begin(); it1 != fCables[lay/2].end(); ++it1) {
	      if (zTr < it1->first || zTr > it1->second) continue;
	      ++nCables;
	    }
	    cout << " Cables: " << nCables << endl;
	    if (nCables) {
	      x0 *= nCables;
	      TMatrixDSym *cov = branchTr->Weight2Cov();
	      Double_t th = branchTr->GetParamNew(3);
	      Double_t cosTh = TMath::Cos(th);
	      Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(branchTr, 9.36, x0, mass2);
	      cout << " Scat: " << hit->GetLayer() << " " << cosTh << " " << TMath::Sqrt(angle2) << endl;
	      (*cov)(2,2) += (angle2 / cosTh / cosTh);
	      (*cov)(3,3) += angle2;
	      Int_t iok = 0;
	      MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
	      branchTr->SetWeight(*cov);
	    }
	    */
	  } 
	} //if (hit->GetDetectorID() != firstHit)
 
	if (hit->GetDetectorID() == skipHit) {
	  // Track and hit in different detectors
	  continue; // next hit
	}
	//cout << hit->GetDetectorID() << endl;
	// Exclude used hits
	if (hit->GetFlag() != 1) continue;
	if (frozen) continue;
	// !!! Exact ID match
	//if (TrackID(hit) != track->GetTrackID()) continue;

	Double_t dChi2 = MpdKalmanFilter::Instance()->FilterStripLocal(branchTr,hit,pointWeight,param,posNew);
	//cout << " Chi2: " << lay << " " << dChi2 << endl;
	// Add Z-contribution (if track is outside the detector)   
	Double_t sizeZ = MpdKalmanFilter::Instance()->GetGeo()->Size(hit).Y();
	//if (TMath::Abs(branchTr->GetParamNew(1)) > sizeZ) {
	if (TMath::Abs(param(1,0)) > sizeZ) {
	  // Outside detector
	  //Double_t dChi2z = (TMath::Abs(branchTr->GetParamNew(1)) - sizeZ) / hit->GetErr(1);
	  Double_t dChi2z = (TMath::Abs(param(1,0)) - sizeZ) / hit->GetErr(1);
	  //AZ!!! dChi2 += dChi2z * dChi2z;
	}
	//cout << dChi2 << " " << TrackID(hit) << endl;

	if (TMath::Abs(dChi2) < fgkChi2Cut) {
	//if (dChi2 < fgkChi2Cut && lay % 2 == 0) {
	  if (nBr < maxBrLay) {
	    branchLay[nBr] = *branchTr;
	    branchLay[nBr].GetHits()->Add(hit);
	    branchLay[nBr].SetChi2(branchTr->GetChi2()+dChi2);
	    TMatrixDSym w = *branchTr->GetWeight();
	    w += pointWeight;
	    branchLay[nBr].SetWeight(w);
	    branchLay[nBr].SetPosNew(branchTr->GetPosNew());
	    //if (hit->GetType() != MpdKalmanHit::kFixedP) branchLay[nBr].SetPosNew(curTr->GetPosNew());
	    //else branchLay[nBr].SetPosNew(posNew); // modular geometry
	    branchLay[nBr].SetLength(branchTr->GetLength());
	    branchLay[nBr].SetParamNew(param);
	    // Save track params at last hit
	    branchLay[nBr].SetLengAtHit(branchTr->GetLength()); 
	    branchLay[nBr].SetParamAtHit(param);
	    branchLay[nBr].SetWeightAtHit(*branchLay[nBr].GetWeight());
	    ++nBr;
	    ++nAdd;
	  } 
	}
      } // for (Int_t indx = indxBeg; indx != indxEnd;

      // Add branch with missing layer
      if (nAdd == 0 && nBr < maxBrLay && ok != -1) {
	it = trackBrM.begin();
	branchLay[nBr] = trackBr[it->second];
	++it;
	// Take branch with greater length
	for ( ; it != trackBrM.end(); ++it) {
	  if (trackBr[it->second].GetLength() > branchLay[nBr].GetLength()) branchLay[nBr] = trackBr[it->second];
	}
	++nBr;
      }

    } // for (Int_t ibr = 0; ibr < nBr0;
    maxInLay = TMath::Max (maxInLay,nBr);
    if (nBr == 0) {
      if (ok == -1) { 
	//cout << track->GetTrackID() << " " << lay << endl;
	//return ok; // stop track
	cout << branch[0].GetNode() << ", " << branch[0].GetNodeNew() << endl;
	break;
      } 
      cout << " !!! MpdTrackFinderIts::RunKalmanFilter: it can not happen |" << endl;
      exit(0);
    } else if (nBr <= maxBr) {
      for (Int_t i = 0; i < nBr; ++i) {
	branch[i].Reset();
	branch[i] = branchLay[i];
      }
    } else {
      // Too many branches - take the best ones
      for (Int_t i = 0; i < nBr; ++i) {
	Double_t c2 = TMath::Min (branchLay[i].GetChi2(),200.);
	quality[i] = branchLay[i].GetNofHits() + (0.999-c2/201);
      }
      TMath::Sort(nBr,quality,indx_lay);
      for (Int_t i = 0; i < maxBr; ++i) {
	branch[i].Reset();
        branch[i] = branchLay[indx_lay[i]];
      }
      nBr = maxBr;
    }

  } // for (Int_t lay = layOK-1; lay >= 0;

  // Select the best branch
  cout << " Branches: " << nBr << " " << maxInLay << endl;
  Int_t ibest = 0;
  for (Int_t i = 1; i < nBr; ++i) {
    cout << i << " " << branch[i].GetNofHits() << " " << branch[i].GetChi2() << endl;
    if (branch[i].GetNofHits() > branch[ibest].GetNofHits()) ibest = i;
    else if (branch[i].GetNofHits() == branch[ibest].GetNofHits() &&
	     branch[i].GetChi2() < branch[ibest].GetChi2()) ibest = i;
  }
  track->Reset();
  *track = branch[ibest];
  return 0;
}

//__________________________________________________________________________
Bool_t MpdTrackFinderIts::NavigateToLayer(Int_t lay, MpdItsKalmanTrack *curTr, MpdItsKalmanTrack *trackBr, 
					  std::map<Int_t,Int_t>& trackBrM)
{
  ///< Navigate track to layer in modular geometry

  const Double_t posCable[4] = {0.08, 0.12, 0.24, 0.40}; // average cable positions (1-4) w.r.t. sensor

  Int_t inode = 0, detID[3] = {-1,-1,-1}, nCables = 0; 
  TString curPath, nodeNew, mass2 = "0.0194797849"; // pion mass squared
  Double_t zTr = 0.0, x0 = 0.0, posNew = 0.0;
  MpdKalmanHit hitTmp;
  TGeoNode *node = NULL;
  TMatrixD parNew = *curTr->GetParamNew();
  //*
  if (lay % 2 == 0) {
    TString path = curTr->GetNodeNew();
    Int_t last = path.Length() - 1;
    Int_t ladd = path.Index("ladder") - 5;
    TString name = path(ladd, last-ladd+1) + "#";
    name += (lay % 2);
    detID[0] = MpdKalmanFilter::Instance()->GetGeo()->DetId(name);
    goto skip;
  }
  //*/

  posNew = curTr->GetPosNew();
  nodeNew = curTr->GetNodeNew();
  //cout << " Nodes: " << nodeNew << " " << curTr->GetNode() << " " << posNew << " " << curTr->GetPos() << endl;
  //curTr->GetParamNew()->Print();
  //curTr->GetParam()->Print();
  curPath = curTr->GetNode();
  hitTmp.SetType(MpdKalmanHit::kFixedR);
  hitTmp.SetPos(fRad[lay]);
  if (!MpdKalmanFilter::Instance()->PropagateParamR(curTr,&hitTmp,kFALSE)) {
    curTr->SetPosNew(posNew);
    curTr->SetParamNew(parNew);
    curTr->SetNodeNew(nodeNew);
    curTr->SetNode(curPath);
    //ok = -1; 
    //break;
    return kFALSE;
  }

  // Find sensor (sector)
  Double_t v7[7], v77[7];
  v77[6] = -1;
  curTr->SetNode(curTr->GetNodeNew());
  MpdKalmanFilter::Instance()->SetGeantParamB(curTr, v7, 1);
  //MpdKalmanFilter::Instance()->SetGeantParamB(curTr, v7, -1);
  node = gGeoManager->FindNode(v7[0], v7[1], v7[2]);
  //if (!TString(gGeoManager->GetPath()).Contains("sts01")) break; // outside ITS coverage
  if (!TString(gGeoManager->GetPath()).Contains("sts01")) return kFALSE; // outside ITS coverage
  zTr = v7[2]; 

  for (Int_t iover = 0; iover < 2; ++iover) {
    // Loop to find first detector
    if (iover == 0) {
      for (Int_t i2 = 0; i2 < 2; ++i2) {
	TString path = gGeoManager->GetPath();
	//cout << path << endl;
	if (!path.Contains("sensor") && !path.Contains("sector")) {
	  if (!path.Contains("ladder")) {
	    gGeoManager->SetStep(100);
	    gGeoManager->FindNextDaughterBoundary(v7,&v7[3],inode);
	    if (inode >= 0) { gGeoManager->CdDown(inode);
	      /*cout << " Node: " << inode << " " << gGeoManager->GetPath() << endl;*/ }
	    else break;
	  } else {
	    gGeoManager->MasterToLocal(v7,v77);
	    gGeoManager->MasterToLocalVect(&v7[3],&v77[3]);
	    gGeoManager->SetStep(100);
	    // Step in normal direction (in transverse plane)
	    Double_t r2 = v77[3] * v77[3] + v77[4] * v77[4];
	    v77[3] = 0;
	    v77[4] = TMath::Sqrt (r2) * TMath::Sign(1.,v77[4]);
	    //v77[4] = TMath::Sqrt (r2) * TMath::Sign(1.,-v77[4]*v77[0]); // for rotated ladders
	    gGeoManager->FindNextDaughterBoundary(v77,&v77[3],inode);
	    if (inode >= 0) {
	      gGeoManager->CdDown(inode);
	      for (Int_t i3 = 0; i3 < 3; ++i3) v77[i3] += v77[i3+3] * (gGeoManager->GetStep() + 0.04);
	      v77[6] = 1;
	    }
	    //cout << " ok " << gGeoManager->GetCurrentPoint()[0] << " " <<  gGeoManager->GetCurrentPoint()[1] << " " << gGeoManager->GetCurrentPoint()[2] << " " << gGeoManager->GetStep() << endl;
	    zTr += (gGeoManager->GetStep() * v77[5]);
	    //cout << zTr << endl;
	  }
	} else break;
      } // for (Int_t i2 = 0; i2 < 2;
    } else { // if (iover == 0)
      // Look for overlap (at first in Z)
      TString path = gGeoManager->GetPath();
      while (path.Contains("sensor") || path.Contains("sector")) {
	// find ladder
	gGeoManager->CdUp(); 
	path = gGeoManager->GetPath();
      }
      if (v77[6] < 0) {
	gGeoManager->MasterToLocal(v7,v77);
	gGeoManager->MasterToLocalVect(&v7[3],&v77[3]);
	Double_t r2 = v77[3] * v77[3] + v77[4] * v77[4];
	v77[3] = 0;
	v77[4] = TMath::Sqrt (r2) * TMath::Sign(1.,v77[4]);
	for (Int_t i3 = 0; i3 < 3; ++i3) v77[i3] += v77[i3+3] * 0.04;
      }
      gGeoManager->FindNextDaughterBoundary(v77,&v77[3],inode);
      if (inode >= 0) gGeoManager->CdDown(inode);
      else {
	// Look for overlap in adjacent ladders 
	Double_t secPos[3];
	gGeoManager->LocalToMaster(v77,secPos);
	gGeoManager->CdUp(); 
	path = gGeoManager->GetPath();
	Double_t r2 = TMath::Sqrt (v7[3] * v7[3] + v7[4] * v7[4]);
	//v7[5] = 0.0;
	//cout << secPos[0] << " " << secPos[1] << endl;
	for (Int_t i3 = 0; i3 < 3; ++i3) secPos[i3] += v7[i3+3] * 0.1 / r2;
	//cout << secPos[0] << " " << secPos[1] << " " << gGeoManager->GetStep() << endl;
	gGeoManager->SetStep(4);
	gGeoManager->FindNextDaughterBoundary(secPos,&v7[3],inode);
	if (inode >= 0) { 
	  // Found ladder - find sector (sensor)
	  for (Int_t i3 = 0; i3 < 3; ++i3) secPos[i3] += v7[i3+3] * (gGeoManager->GetStep() + 0.001);
	  gGeoManager->CdDown(inode); 
	  gGeoManager->MasterToLocal(secPos,v77);
	  gGeoManager->MasterToLocalVect(&v7[3],&v77[3]);
	  //cout << path << " " << gGeoManager->GetPath() << " " << gGeoManager->GetStep() << endl; 
	  gGeoManager->SetStep(1);
	  gGeoManager->FindNextDaughterBoundary(v77,&v77[3],inode);
	  //cout << inode << endl;
	  //if (inode >= 0) exit(0);
	  if (inode >= 0) gGeoManager->CdDown(inode); 
	}
      }
    } // if (iover == 0)
    TString path = gGeoManager->GetPath();
    if (!path.Contains("ladder")) continue;
    if (!path.Contains("sensor") && !path.Contains("sector")) {
      //if (iover) continue;
      continue;
      //Fatal(" RunKalmanFilter ", path);
    }
    //cout << path << endl;
    Int_t last = path.Length() - 1;
    //if (path.Contains("sensor") || path.Contains("frame")) last = path.Last('/') - 1;
    Int_t ladd = path.Index("ladder") - 5;
    //TString name = path(ladd, last-ladd+1);
    TString name = path(ladd, last-ladd+1) + "#";
    //name += (!track->GetDirection());
    name += (lay % 2);
    //cout << name << endl;
    detID[iover] = MpdKalmanFilter::Instance()->GetGeo()->DetId(name);
    // Debugging
    //TVector3 v3tmp;
    //MpdStsHit hitDb1(detID[iover]%1000000,v3tmp,v3tmp,0.0,-1);
    //MpdStsHit hitDb2(hit->GetDetectorID()%1000000,v3tmp,v3tmp,0.0,-1);
    //cout << iover << " " << name << " " << detID[iover] << " " << hit->GetDetectorID() << " " << hitDb1.Ladder() << " " << hitDb1.Sector() << " " << hitDb2.Ladder() << " " << hitDb2.Sector() << " " << hit->GetMeas(0) << endl;
    break;
  } // for (Int_t iover = 0; iover < 2;

  //if (detID[0] < 0 && detID[1] < 0) break;
  if (detID[0] < 0 && detID[1] < 0) return kFALSE;
  if (detID[0] < 0) { detID[0] = detID[1]; detID[1] = -1; }
	  
  // if (detID[0] >= 0 && detID[1] >= 0) exit(0); // Debugging
  curTr->SetPosNew(posNew);
  curTr->SetParamNew(parNew);
  curTr->SetNodeNew(nodeNew);
  curTr->SetNode(curPath);

  //firstDet = 0;

  if (lay % 2 != 0 && fCables[lay/2].size() > 0) {
    // Find number of cables crossed                                    
    TString path = gGeoManager->GetPath();
    if (!path.Contains("sensor") && !path.Contains("sector")) {
      cout << " !!! MpdTrackFinderIts::RunKalmanFilterMod - Outside detector !!! " << endl;
      exit(0);
    }
    zTr = TMath::Abs (zTr); // global Z
    map<Double_t,Double_t>::iterator it1;
    for (it1 = fCables[lay/2].begin(); it1 != fCables[lay/2].end(); ++it1) {
      if (zTr < it1->first || zTr > it1->second) continue;
      ++nCables;
    }
  }
  //cout <<  " zTr " << zTr << " " << nCables << endl;    
  //nCables = 0;
 
  // Crossing silicon layer - add mult. scat. in the cable            
  //Double_t x0 = 0.0116 * 2 / 9.36 * nCables; // in rad. length - 116um cable per side 
  x0 = 0.0116 * 2 * nCables; // in rad. length - 116um cable per side 
  if (fMCTracks) {
    // Get particle mass - ideal PID
    MpdMCTrack *mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(curTr->GetTrackID());
    TParticlePDG *pdgP = TDatabasePDG::Instance()->GetParticle(mctrack->GetPdgCode());
    if (pdgP) {
      Double_t mass = pdgP->Mass();
      if (mass < 0.1 || mass > 0.8) {
	// Electrons or protons (or heavier than protons)
	mass2 = "";
	mass2 += mass*mass;
      }
    }
  }

 skip:

  Int_t firstBranch = 1;

  // Propagate track to the found detector and look for the overlapping sensors and ladders
  hitTmp.SetType(MpdKalmanHit::kFixedP);

  for (Int_t iover = 0; iover < 3; ++iover) {
    if (detID[iover] >= 0) {
      // Already found detector
      hitTmp.SetDetectorID(detID[iover]);

      if (nCables) {
	// Propagate to average cable position
	MpdKalmanHit hitCable;
	hitCable.SetDetectorID(-2); // just some fake detector ID
	MpdKalmanGeoScheme *geo = MpdKalmanFilter::Instance()->GetGeo();
	TVector3 pos = geo->GlobalPos(&hitTmp);
	//pos.Print();
	pos += geo->Normal(&hitTmp) * posCable[nCables-1];
	//pos.Print();
	geo->SetGlobalPos(-2, pos, kTRUE);
	geo->SetNormal(-2, geo->Normal(&hitTmp), kTRUE);
	geo->SetPath(-2, geo->Path(detID[iover]), kTRUE);
	if (!MpdKalmanFilter::Instance()->PropagateToHit(&trackBr[iover],&hitCable,kTRUE,kTRUE)) { 
	  detID[iover] = -1;
	  continue;
	}
	TMatrixDSym *cov = trackBr[iover].Weight2Cov();
	Double_t th = trackBr[iover].GetParamNew(3);
	Double_t cosTh = TMath::Cos(th);
	x0 /= cosTh;
	Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(&trackBr[iover], 9.36, x0*2, mass2);
	//cout << TMath::Sqrt(angle2) << " " << cosTh << endl;
	(*cov)(2,2) += (angle2 / cosTh / cosTh);
	(*cov)(3,3) += angle2;
	Int_t iok = 0;
	MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
	trackBr[iover].SetWeight(*cov);
      }

      if (!MpdKalmanFilter::Instance()->PropagateToHit(&trackBr[iover],&hitTmp,kTRUE,kTRUE)) { 
	//cout << detID[iover] << " " << MpdKalmanFilter::Instance()->GetGeo()->Path(detID[iover]) << endl;
	//cout << iover << " " << detID[iover] << " " << trackBr[iover].GetNode() << " " << trackBr[iover].GetNodeNew() << " " << trackBr[iover].GetNofHits() << endl;
	detID[iover] = -1;
	continue;
      }

      if (firstBranch) {
	// Try to get the closest neighbours
	// Check if the track is close to the detector edges
	Double_t xloc = trackBr[iover].GetParamNew(0);
	Double_t zloc = trackBr[iover].GetParamNew(1);
	hitTmp.SetDetectorID(detID[iover]);
	TVector2 size = MpdKalmanFilter::Instance()->GetGeo()->Size(&hitTmp);
	Int_t idSTS = detID[iover] % 1000000;
	//cout << xloc << " " << zloc << " " << size.X() << " " << size.Y() << endl;
	
	if (TMath::Abs(TMath::Abs(xloc)-size.X()) < 0.3) {
	  Int_t ladder = fHitSts.GetLadder(idSTS);
	  if (xloc < 0) {
	    ++ladder;
	    if (ladder > fNladders[lay]) ladder = 1;
	  } else {
	    --ladder;
	    if (ladder == 0) ladder = fNladders[lay];
	  }
	  if (iover < 2 && detID[iover+1] < 0) detID[iover+1] = 1000000*lay + fHitSts.SetLadder(idSTS,ladder);
	  else if (iover < 1 && detID[iover+2] < 0) detID[iover+2] = 1000000*lay + fHitSts.SetLadder(idSTS,ladder);
	}  
	if (TMath::Abs(TMath::Abs(zloc)-size.Y()) < 0.3) {
	  Int_t sector = fHitSts.GetSensor(idSTS);
	  if (zloc < 0) --sector;
	  else ++sector;
	  if (sector > 0 && sector <= fNsectors[lay]) {
	    if (iover < 2 && detID[iover+1] < 0) detID[iover+1] = 1000000*lay + fHitSts.SetSensor(idSTS,sector);
	    else if (iover < 1 && detID[iover+2] < 0) detID[iover+2] = 1000000*lay + fHitSts.SetSensor(idSTS,sector);
	  }
	}
	firstBranch = 0;
      }	
    }
  }
  if (firstBranch) {
    // Failed propagation
    cout << " !!! Failed navigation to layer " << lay << endl;
    //exit(0);
    return kFALSE;
  }
  for (Int_t iover = 0; iover < 3; ++iover) {
    if (detID[iover] < 0) continue;
    trackBrM.insert(pair<Int_t,Int_t>(detID[iover],iover));
  }

  return kTRUE;
}

//__________________________________________________________________________
Int_t MpdTrackFinderIts::TrackID(MpdKalmanHit *hit)
{
  /// Return track ID of the hit

  FairHit *h = (FairHit*) fItsHits->UncheckedAt(hit->GetIndex());
  return ((MpdStsPoint*) fItsPoints->UncheckedAt(h->GetRefIndex()))->GetTrackID();
}

//__________________________________________________________________________
TVector2 MpdTrackFinderIts::GetDistance(MpdKalmanTrack *track, MpdKalmanHit *hit)
{
  /// Compute distance between track and hit

  Int_t lay = hit->GetLayer();
  Int_t lay2 = lay / 2;
  Int_t side = lay % 2;
  Int_t module = ((MpdStsHit*) fItsHits->UncheckedAt(hit->GetIndex()))->Module();

  Double_t zTr = track->GetParamNew(1);
  Double_t zloc = zTr + fDz[lay2];
  Int_t modTr = Int_t (zloc / fZmod[lay2]);
  Int_t dMod = modTr - module;

  Double_t dZ = 0;
  if (dMod != 0) {
    // Not in the same module - check Z-distance to module edge
    if (dMod > 0) dZ = zloc - (module+1) * fZmod[lay2];
    else dZ = zloc - module * fZmod[lay2];
    if (TMath::Abs(dMod) > 2) return TVector2(TMath::Abs(dZ),999.); // not in neighbour modules
  }

  // Translate transverse track position to local system
  Double_t xloc = track->GetParamNew(0) * hit->GetCosSin(0) + zTr * hit->GetCosSin(1);
  Double_t phTr = xloc / track->GetPosNew();

  Double_t phHit = hit->GetMeas(0) / fRad[lay]; 
  //TVector2 dist(TMath::Abs(dZ), TMath::Abs(MpdKalmanFilter::Instance()->Proxim(phTr,phHit)-phTr));
  TVector2 dist(TMath::Abs(dZ), 
		TMath::Abs(MpdKalmanFilter::Instance()->Proxim(phTr,phHit,hit->GetCosSin(0))-phTr));
  //TVector2 dist(TMath::Abs(zTr-hit->GetZ()), TMath::Abs(MpdKalmanFilter::Instance()->Proxim(track,hit)/hit->GetR()-ph));
  return dist;
}

//__________________________________________________________________________
void MpdTrackFinderIts::Write()
{
  /// Write

  TFile histoFile("ItsRec.root","RECREATE");
  Writedir2current(fHistoDir);
  histoFile.Close();
}

//__________________________________________________________________________
void MpdTrackFinderIts::Writedir2current( TObject *obj )
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
void MpdTrackFinderIts::AddHits()
{
  /// Add hit objects to tracks and compute number of wrongly assigned hits
  /// (hits with ID different from ID of starting TPC track)

  Int_t nFound = fTracks->GetEntriesFast();
  for (Int_t i = 0; i < nFound; ++i) {
    MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTracks->UncheckedAt(i);
    cout << track->GetNode() << " " << track->GetNodeNew() << endl;
    Double_t c2 = track->GetChi2();
    track->SetChi2(track->GetChi2Its());
    track->SetChi2Its(c2);
    Int_t nHits = track->GetNofHits();
    //if (nHits == 0) { fTracks->RemoveAt(i); continue; }
    TClonesArray &trHits = *track->GetTrHits();
    cout << nHits << " " << trHits.GetEntriesFast() << " " << track->GetTrackID() << endl;
    TObjArray *hits = track->GetHits();
    Int_t nWrong = 0, nMirr = 0, motherID = track->GetTrackID();
    // Get track mother ID 
    MpdMCTrack *mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(motherID);
    while (mctrack->GetMotherId() >= 0) {
      motherID = mctrack->GetMotherId();
      mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(mctrack->GetMotherId());
    }

    Int_t lastIndx = trHits.GetEntriesFast();
    for (Int_t j = 0; j < nHits; ++j) {
      MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(j);
      hit->SetUniqueID(1); // flag ITS hits
      new (trHits[lastIndx+j]) MpdKalmanHit(*hit);
      cout << " " << hit->GetLayer();
      MpdStsHit *h = (MpdStsHit*) fItsHits->UncheckedAt(hit->GetIndex());
      Int_t motherID1 = ((FairMCPoint*) fItsPoints->UncheckedAt(h->GetRefIndex()))->GetTrackID();
      cout << "-" << motherID1;
      // Get point mother ID 
      mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(motherID1);
      while (mctrack->GetMotherId() >= 0) {
        motherID1 = mctrack->GetMotherId();
        mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(mctrack->GetMotherId());
      }
      if (motherID1 != motherID) ++nWrong;
    }
    if (nHits) cout << "\n" << nWrong << endl;
    track->SetNofWrong(nWrong);
    MpdKalmanTrack *tpc = (MpdKalmanTrack*) fTpcTracks->UncheckedAt(track->GetUniqueID()-1);
    //track->SetChi2(track->GetChi2()+tpc->GetChi2());
    //track->SetLastLay();
    //track->GetParam()->Print();
    track->SetNofHits(track->GetNofTrHits()); // TPC and ITS hits
    track->SetNofIts(nHits);
    cout << nHits << " " << track->GetNofTrHits() << " " << track->GetTrackID() << " " 
	 << track->GetChi2Its() << " " << track->GetChi2() << endl;
  }
  fTracks->Compress();
}

ClassImp(MpdTrackFinderIts);
