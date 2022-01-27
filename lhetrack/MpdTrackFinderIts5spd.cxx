// -------------------------------------------------------------------------
// -----                  MpdTrackFinderIts5spd source file            -----
// -----                 Created 21/07/09  by A. Zinchenko             -----
// -------------------------------------------------------------------------
// -----   Modified by VK for ideal hit producer and 5 spd-layers geometry -------

/**
 *@author A.Zinchenko <Alexander.Zinchenko@jinr.ru>
 **
 ** Track finder in MPD Inner Tracking System (ITS) using seeds from TPC
 **/

#include "MpdStsGeoPar.h"
#include "MpdTrackFinderIts5spd.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanTrack.h"
//#include "MpdKalmanStripHit.h"
#include "MpdItsHit5spd.h"
#include "MpdStsPoint.h"
#include "MpdItsKalmanTrack.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdMCTrack.h"

#include "FairGeoNode.h"
#include "FairMCPoint.h"
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
#include "TFile.h"
//#include "TLorentzVector.h"
#include "TVector2.h"
#include "TClonesArray.h"
#include <TRandom.h>

#include <algorithm>
#include <iostream>
#include <map>

using std::cout;
using std::endl;
//using std::map;

const Double_t MpdTrackFinderIts5spd::fgkChi2Cut = 20; //20; //100;

//__________________________________________________________________________
MpdTrackFinderIts5spd::MpdTrackFinderIts5spd(Bool_t useVector, const char *name, Int_t iVerbose )
  : FairTask(name, iVerbose),
    fExact(0), //1),
    fGeo(0),
    fUseVector(useVector),
    fTpcKF(nullptr)
{
  
  if (!fUseVector) fKHits = new TClonesArray("MpdKalmanHit", 100);
  fTracks = new TClonesArray("MpdItsKalmanTrack", 100);
  fHistoDir = 0x0;
  fhLays = new TH1F("hLaysITS","ITS layers",12,0,12);
  fLayPointers = 0x0;
}


//__________________________________________________________________________
MpdTrackFinderIts5spd::~MpdTrackFinderIts5spd()
{
  //delete fKHits;
  //delete fTracks;
  //delete fTrackCand;
  delete [] fLayPointers;
  delete fhLays;
}

//__________________________________________________________________________
InitStatus MpdTrackFinderIts5spd::Init()
{

  //----- Get pipe radius
   fPipeR = ((TGeoTube*)gGeoManager->GetVolume("pipe1")->GetShape())->GetRmax();
   cout << "************** Pipe radius: " << fPipeR << " cm " << endl;

  //return ReInit();
  if (ReInit() != kSUCCESS) return kERROR;

  //----- Read database

  FairRuntimeDb* rtdb = FairRun::Instance()->GetRuntimeDb();
  MpdStsGeoPar *geoPar = (MpdStsGeoPar*) rtdb->getContainer("MpdStsGeoPar");
  // cout << geoPar << endl;
  TString volName = "sts01 ", path = "";
  TObjArray* sensNodes = geoPar->GetGeoSensitiveNodes();
  // cout << sensNodes->GetEntriesFast() << " " << geoPar->GetGeoPassiveNodes()->GetEntriesFast() << endl;

  //----- Process modular geometry

   fGeo = 1; 
   Int_t nLay = 5;
   //Double_t safety = 0.005;

   for (Int_t i = 0; i < nLay; ++i) {
     fNladders[i] = fNsectors[i] = 0;
     volName = "sts01ladder";
     volName += (i+1);
     path = "/cave_1/sts01_0/" + volName;
     path += "_";
     TString path0 = path;
     
     //----- Loop over all ladders to find the one with the smallest radius
     //fRad[i] = 999.9;
     fRad[i] = 0.;
     Int_t nladd = 0;
     
     for (Int_t jlad = 1; jlad < 99; ++jlad) {
       path = path0;
       path += jlad;
       if (!gGeoManager->CheckPath(path)) break;
       gGeoManager->cd(path);
       TGeoVolume *ladd = gGeoManager->GetVolume(volName);
       //TGeoBBox *box = (TGeoBBox*) ladd->GetShape();
       Double_t xyzL[3] = {0}, xyzM[3];
       gGeoManager->LocalToMaster(xyzL,xyzM);
       Double_t rad = TMath::Sqrt (xyzM[0] * xyzM[0] + xyzM[1] * xyzM[1]);
       fRad[i] += rad;
       ++nladd;
       /*AZ - 3.12.2017
       fRad[i] = TMath::Min (fRad[i],rad);
       xyzL[0] = box->GetDX();
       gGeoManager->LocalToMaster(xyzL,xyzM);
       rad = TMath::Sqrt (xyzM[0] * xyzM[0] + xyzM[1] * xyzM[1]);
       fRad[i] = TMath::Min (fRad[i],rad);
       xyzL[0] = -box->GetDX();
       gGeoManager->LocalToMaster(xyzL,xyzM);
       rad = TMath::Sqrt (xyzM[0] * xyzM[0] + xyzM[1] * xyzM[1]);
       fRad[i] = TMath::Min (fRad[i],rad);
       */
     }
     cout << " *** Layer # " << i << " min_R= " <<  fRad[i]  << endl;
     TGeoVolume *ladd = gGeoManager->GetVolume(volName);
     if (ladd == nullptr) { nLay = i; break; }
     fRad[i] /= nladd; // mean radius
     /*
     TGeoBBox *box = (TGeoBBox*) ladd->GetShape();
     //safety = -box->GetDY();
     //safety = box->GetDY();
     safety = 2 * box->GetDY(); // new
     fRad[i] -= safety;
     */
     
     cout << " +++ Layer # " << i << " safety min_R= " <<  fRad[i] << endl;
     
   }
   FillGeoScheme();
   
   //----- Get pipe radius
   //fPipeR = ((TGeoTube*)gGeoManager->GetVolume("pipe1")->GetShape())->GetRmax();
   //cout << "************** Pipe radius: " << fPipeR << " cm " << endl;

   //----- Get cables
   TObjArray *vols = gGeoManager->GetListOfVolumes();
   Int_t nvols = vols->GetEntriesFast(), ncables = 0;
   for (Int_t i = 0; i < nvols; ++i) {
     TGeoVolume *vol = (TGeoVolume*) vols->UncheckedAt(i);
     TString cable = TString(vol->GetName());
     if (!(cable.Contains("sts") && cable.Contains("cable"))) continue;
     //  cout << cable << endl;
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
     //  cout << xyzM[2] - box->GetDZ() << " " << xyzM[2] + box->GetDZ() << " " << box->GetDZ() << endl;
     if (xyzM[2] - box->GetDZ() > 0) {
       Int_t lay = TString(lad(lad.Length()-4,1)).Atoi();
       fCables[lay-1].insert(pair<Double_t,Double_t>(xyzM[2] - box->GetDZ(),xyzM[2] + box->GetDZ()));
     }
   }
   cout << "-I- MpdTrackFinderIts5spd: Intialization finished successfully" << endl;
   return kSUCCESS;
}

//__________________________________________________________________________
void MpdTrackFinderIts5spd::FillGeoScheme()
{
  //----- Fill Kalman filter geometry manager info

  MpdKalmanGeoScheme *geo = MpdKalmanFilter::Instance()->GetGeo();

  TGeoVolume *vSts = gGeoManager->GetVolume("sts01");

  //----- Loop over layers
  for (Int_t layer = 1; layer < 999; ++layer) {
    TString sladder = "sts01ladder";
    sladder += layer;
    TGeoVolume *vLay = gGeoManager->GetVolume(sladder);
    if (vLay == 0x0 && layer == 1) continue; // no first layer
    if (vLay == 0x0) break;
    sladder += "_";
 
    //----- Loop over ladders
    for (Int_t ladder = 1; ladder < 999; ++ladder) {
      TString sladder1 = sladder;
      sladder1 += ladder;
      TGeoNode *node = vSts->FindNode(sladder1);
      if (node == 0x0) break;
      ++fNladders[layer-1];
      TGeoVolume *vLad = node->GetVolume();
      //cout << "Node: " << layer << " " << vLad->GetNodes()->GetEntriesFast() << " " << vLad->GetNdaughters() << endl;
      Int_t nDaught = vLad->GetNdaughters(), detID = -1, detIDsts = -1;
      TObjArray *daught = vLad->GetNodes();
      
      //      for (Int_t j = 0; j < vLad->GetNdaughters(); ++j) 
      //      cout << "Name: " << ((TGeoNode*)(vLad->GetNodes()->UncheckedAt(j)))->GetName() << endl;

      //----- Loop over ladder daughters
      Int_t iZ = 0;
      for (Int_t det = 0; det < nDaught; ++det) {
	TString sdet1 = ((TGeoNode*)(daught->UncheckedAt(det)))->GetName();
	//	if (!sdet1.Contains("sensor") && !sdet1.Contains("sector")) continue;
	if (ladder == 1)  ++fNsectors[layer-1];
	Int_t det1 = TString(sdet1(sdet1.Index("_")+1,2)).Atoi();
//	Int_t secType = -1;
//	if (sdet1.Contains("sector")) secType = TString(sdet1(sdet1.Index("_")-2,1)).Atoi();

        ++iZ;

	Int_t side = 0.;

	detIDsts = fHitSts.SetDetId(layer, ladder, det1);
	detID = fHitSts.SetDetId(layer, ladder, iZ);
	//cout << "det1= " << det1 << " " << iZ << " layer= " << layer << " ladder="  << ladder << " detID= " << detID << endl;
	fId2Id[layer-1].insert(pair<Int_t,Int_t>(detIDsts,detID)); // detIDsts == detID for pixel detectors (AZ)
	detID += 1000000 * (layer-1);
	
	TString detName = sladder1 + "/" + sdet1 + "#";
	detName += side;
	
	geo->SetDetId(detName, detID);
	
	TString path = "/cave_1/sts01_0/" + detName(0,detName.Length()-2);
	//cout << "DetName: " << detName << " Path: " << path << endl;
	gGeoManager->cd(path);
	
	geo->SetPath(detID, path);
	node = gGeoManager->GetCurrentNode();
	//        cout << node->GetName() << " " << detID << endl;
	
	TGeoVolume *vol = node->GetVolume();
	TGeoBBox *box = (TGeoBBox*) vol->GetShape();
	TVector2 size(box->GetDX(), box->GetDZ());
	geo->SetSize(detID, size);
	//	  cout << "##### Node: " << node->GetName() <<  " DX= " << box->GetDX() << " DZ=" << box->GetDZ() << endl;	  
	
	Double_t xyzL[3] = {0}, xyzM[3], vecM[3];
	
	xyzL[1] = 1;
	gGeoManager->LocalToMasterVect(xyzL,vecM);
	
	xyzL[1] = 0;
	gGeoManager->LocalToMaster(xyzL,xyzM);
	Double_t sign = TMath::Sign (1.,xyzM[0]*vecM[0]+xyzM[1]*vecM[1]);
	
	TVector3 pos(xyzM[0], xyzM[1], xyzM[2]);
	geo->SetGlobalPos(detID, pos);
	
	TVector3 norm(sign*vecM[0], sign*vecM[1], sign*vecM[2]);
	geo->SetNormal(detID, norm);

      }
    }
  }

}

//__________________________________________________________________________
InitStatus MpdTrackFinderIts5spd::ReInit()
{
  fItsPoints = (TClonesArray *) FairRootManager::Instance()->GetObject("StsPoint");
  fItsHits = (TClonesArray *) FairRootManager::Instance()->GetObject("StsHit");
  if (fItsPoints == 0x0 || fItsHits == 0x0) return kERROR;
  fTpcTracks = (TClonesArray *) FairRootManager::Instance()->GetObject("TpcKalmanTrack");
  fMCTracks = (TClonesArray *) FairRootManager::Instance()->GetObject("MCTrack");
  if (fUseVector) fKHits = (TClonesArray *) FairRootManager::Instance()->GetObject("ItsKHits"); // use from VectorFinder

  FairRootManager::Instance()->Register("ItsTrackKF", "Its", fTracks, kTRUE);

  fNPass = 1;

  return kSUCCESS;
}

//__________________________________________________________________________
void MpdTrackFinderIts5spd::Reset() 
{
  cout << " MpdTrackFinderIts5spd::Reset  " << endl;

  if (!fUseVector) fKHits->Delete();
  fTracks->Delete();
  delete [] fLayPointers;
  fLayPointers = nullptr;
}

//__________________________________________________________________________
void MpdTrackFinderIts5spd::SetParContainers()
{
//----- Get run and runtime database
  FairRunAna* run = FairRunAna::Instance();
  if ( ! run ) Fatal("SetParContainers", "No analysis run");

  FairRuntimeDb* db = run->GetRuntimeDb();
  if ( ! db ) Fatal("SetParContainers", "No runtime database");

//----- Get STS geometry parameter container
  db->getContainer("MpdStsGeoPar");
}

//__________________________________________________________________________
void MpdTrackFinderIts5spd::Finish()
{
//   Write();
}

//__________________________________________________________________________
void MpdTrackFinderIts5spd::Exec(Option_t * option)
{

  static int eventCounter = 0;    
  cout << " ##### MpdTrackFinderIts5spd \n ItsRec event " << ++eventCounter << endl;

  Reset();

  //----- Create Kalman hits
  if (fItsHits->GetEntriesFast() == 0) return;
  MakeKalmanHits();

  for (Int_t i = 0; i < fNPass; ++i) {
    fTracks->Clear();
    GetTrackSeeds(i);

    cout << "  Total number of hits for tracking: " << fKHits->GetEntriesFast() << endl;
    cout << "  Total number of track seeds: " << fTracks->GetEntriesFast() << endl;

    DoTracking(i);
//  StoreTracks();
    cout << "  Total number of found tracks: " << fTracks->GetEntriesFast() << endl;
    if (i != fNPass - 1) ExcludeHits(); // exclude used hits
  }
  AddHits(); // add hit objects to tracks
  cout << "  Total number of found tracks: " << fTracks->GetEntriesFast() << endl;
}

//__________________________________________________________________________
void MpdTrackFinderIts5spd::MakeKalmanHits()
{
  //----- Create Kalman hits from ITS hits.

  fhLays->Reset();

  Int_t nHits = fItsHits->GetEntriesFast(), layMax = 0, lay = 0, nKH = 0;
  //Double_t r, z, xloc, errZ = 0.00012, errX = 0.00023; // 120um in Z, 23um in R-Phi (local X)
  //Double_t r, z, xloc, errZ = 0.012, errX = 0.0023; // 120um in Z, 23um in R-Phi (local X)
  //Double_t r, z, xloc, errZ = 0.12, errX = 0.0023; // 1.2mm in Z, 23um in R-Phi (local X)
  Double_t r, z, xloc, errZ = 0.001, errX = 0.001; // 10um in Z, 10um in R-Phi (local X)
  MpdKalmanGeoScheme *geo = MpdKalmanFilter::Instance()->GetGeo();

  if (!fUseVector) {
    for (Int_t ih = 0; ih < nHits; ++ih) {
      MpdItsHit5spd *h = (MpdItsHit5spd*) fItsHits->UncheckedAt(ih);
      r = TMath::Sqrt (h->GetX() * h->GetX() + h->GetY() * h->GetY());
      z = h->GetZ();
      xloc = h->GetLocalX();
      lay = h->Layer() - 1; 
      //    cout <<" n= " << ih << " r= " << r << " z= " << z << " layer= " << h->Layer()  << " xloc= " << xloc << endl;
      //cout << h->GetDetectorID() << " " << geo->Path(h->GetDetectorID()+1000000*lay) << endl;
      // Get local Z
      TString path = geo->Path(h->GetDetectorID()+1000000*lay);
      gGeoManager->cd(path);
      Double_t posLoc[3] = {0}, posMas[3] = {h->GetX(), h->GetY(), h->GetZ()};
      gGeoManager->MasterToLocal(posMas,posLoc);
      
      //----- Add error                                            
      Double_t dX = 0, dZ = 0;
      gRandom->Rannor(dX,dZ);
      //  if (errZ > 2) dZ = 0.0; // 1-D case
      //dZ = 0.0; // 1-D case
      //Double_t meas[2] = {xloc+dX*errX, z+dZ*errZ};
      Double_t meas[2] = {xloc+dX*errX, posLoc[2]+dZ*errZ};
      Double_t err[2] = {errX, errZ};
      Double_t cossin[2] = {TMath::Cos(fStereoA[0]), TMath::Sin(fStereoA[0])};
      
      MpdKalmanHit *hit = 0x0;
      //  cout << h->GetDetectorID() << " " << fId2Id[lay][h->GetDetectorID()] << endl;
      
      if (fGeo && h->GetUniqueID()) {
	hit = new ((*fKHits)[nKH++]) MpdKalmanHit(lay*1000000+fId2Id[lay][h->GetDetectorID()], 1, 
						  MpdKalmanHit::kFixedP, meas, err, cossin, 0., r, ih);
	//  cout << " DetID = " << h->GetDetectorID() << " " << lay*1000000+fId2Id[lay][h->GetDetectorID()] << endl;
      }
      
      hit->SetUniqueID(0);
      hit->SetNofDim(2);
    }
  }
  
  fKHits->Sort(); // in descending order in R
  layMax = ((MpdKalmanHit*)fKHits->First())->GetLayer();
  nKH = fKHits->GetEntriesFast();
  
  for (Int_t ih = 0; ih < nKH; ++ih) {
    MpdKalmanHit *hit = (MpdKalmanHit*) fKHits->UncheckedAt(ih);
    //if (hit->GetFlag() < 0) hit->SetFlag(hit->GetFlag() | MpdKalmanHit::kUsed);
    fhLays->Fill(hit->GetLayer()+0.1);
  }
  
  for (Int_t k = 0; k <= layMax; ++k) 
    cout << "Layer= " << k+1 << " Number of hits=" << fhLays->GetBinContent(k+1,0) << endl; 

  cout << " Max layer = " << layMax << " Number of hits= " << fKHits->GetEntriesFast() << endl;
  
  for (lay = 0; lay < fgkNlays; ++lay) {
    fHitMapRphi[lay].clear();
    fHitMapZ[lay].clear();
  }

  for (Int_t ihit = 0; ihit < nKH; ++ihit) {
    MpdKalmanHit *hk = (MpdKalmanHit*) fKHits->UncheckedAt(ihit);
    lay = hk->GetLayer();
    MpdItsHit5spd *h = (MpdItsHit5spd*) fItsHits->UncheckedAt(hk->GetIndex());
    TVector3 pos;
    h->Position(pos);
    z = pos.Z();
    fHitMapZ[lay].insert(pair<Double_t,Int_t>(z,ihit));
    Double_t phi = pos.Phi();
    //if (phi < 0) phi += TMath::TwoPi();
    Double_t rphi = fRad[lay] * phi;
    fHitMapRphi[lay].insert(pair<Double_t,Int_t>(rphi,ihit));
    Double_t rphiMax = fRad[lay] * TMath::Pi();
    if (TMath::Abs(rphi-rphiMax) < 3.0) {
      // Add overlap 3 cm near phi=180
      if (rphi < 0) fHitMapRphi[lay].insert(pair<Double_t,Int_t>(rphi+2*rphiMax,ihit)); // overlap 3 cm near phi=180
      else fHitMapRphi[lay].insert(pair<Double_t,Int_t>(rphi-2*rphiMax,ihit)); // overlap 3 cm near phi=180
    }
  }
	
  //  for (Int_t ihit = 0; ihit < fKHits->GetEntriesFast(); ++ihit) {
  //  cout << "ITS hit number= " << ihit << " Position= "<< ((MpdKalmanHit*)fKHits->UncheckedAt(ihit))->GetPos() << 
  //  " Distance= " << ((MpdKalmanHit*)fKHits->UncheckedAt(ihit))->GetDist() << endl;
  //  }

  fLayPointers = new Int_t [layMax+1];
  Int_t ipos = 0;
  for (Int_t i = layMax; i > -1; --i) {
    //cout << i << " " << fhLays->GetBinContent(i+1,0) << endl;
    //if (ipos) cout << ((TpcLheHit*)fHits->UncheckedAt(ipos))->GetLayer() << " "
    //     << ((TpcLheHit*)fHits->UncheckedAt(ipos-1))->GetLayer() << endl;
    Int_t cont = TMath::Nint (fhLays->GetBinContent(i+1,0));
    if (cont == 0) {
      fLayPointers[i] = -1;
      continue;
    }
    fLayPointers[i] = ipos;
    ipos += cont;
  }
}

//__________________________________________________________________________
void MpdTrackFinderIts5spd::GetTrackSeeds(Int_t iPass)
{
  /// Build ITS track seeds from TPC tracks

  Int_t nTpcTracks = fTpcTracks->GetEntriesFast();
  cout << " Seed tracks: " << nTpcTracks << endl;
 
  Int_t nCand = 0;
  MpdKalmanHit hit;
  hit.SetType(MpdKalmanHit::kFixedR);
  for (Int_t itr = 0; itr < nTpcTracks; ++itr) {
    MpdTpcKalmanTrack *tpc = (MpdTpcKalmanTrack*) fTpcTracks->UncheckedAt(itr);
    //tpc->GetParam()->Print();
    //MpdEctKalmanTrack *track = new ((*fTracks)[nCand++]) MpdEctKalmanTrack(itr, *tpc);
    MpdItsKalmanTrack *track = new ((*fTracks)[nCand++]) MpdItsKalmanTrack(*tpc);
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
    track->SetDirection(MpdKalmanTrack::kInward); //AZ - 6.12.2017
  }
  cout << " Number of ITS track candidates: " << nCand << endl;
}  
  
//__________________________________________________________________________
void MpdTrackFinderIts5spd::DoTracking(Int_t iPass)
{
//----- Run Kalman tracking
  
  Double_t vert[3] = {0.0,0.0,0.0};
  Int_t nCand = fTracks->GetEntriesFast(), iok = 0;
  Int_t lay0 = ((MpdKalmanHit*)fKHits->First())->GetLayer();
 
  for (Int_t i = 0; i < nCand; ++i) {
    MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTracks->UncheckedAt(i);
//   cout << " Track seed No. " << i << ", ID: " << track->GetTrackID() << ", Hits: " << track->GetNofTrHits() << endl;
    //for (Int_t j = 0; j < track->GetNofTrHits(); ++j) {
      //MpdKalmanHit *h = (MpdKalmanHit* )track->GetTrHits()->UncheckedAt(j);
//    cout << "Track hit number= " << j << " Distance= " << h->GetDist() << " Layer= " << h->GetLayer() << endl;
    //}
    iok = RunKalmanFilterMod(track, lay0); // modular geometry
    if (iok == -1) {
      fTracks->RemoveAt(i);
      continue;
    }

//  if (track->GetNofHits() == 0) continue; // no hits added

//------ Propagate track to the beam line
    track->SetParam(*track->GetParamNew());
    track->SetPos(track->GetPosNew());
    Double_t pos = track->GetPos();
    TMatrixD par = *track->GetParam();
    TMatrixDSym cov = *track->Weight2Cov();
    Double_t leng = track->GetLength();
    TString nodeNew = track->GetNodeNew();
    //cout << " 1: " << nodeNew << ", " << track->GetNode() << " " << track->GetHits()->GetEntriesFast() << endl;
 
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
      //cout << " 2: " << nodeNew << ", " << track->GetNode() << " " << pos << endl;
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
    if (iok != 1) {
      track->SetNodeNew(""); // for FindPca 
      MpdKalmanFilter::Instance()->FindPca(track, vert);
      track->SetCovariance(cov); //AZ-040121
    }
    //track->SetPos(pos); // restore position
    track->SetParam(*track->GetParamNew()); // !!! track params at PCA
    //track->GetCovariance()->Print();
    //cout << " 3: " << track->GetNodeNew() << ", " << track->GetNode() << endl;
  } // for (Int_t i = 0; i < nCand;
  fTracks->Compress();
}
    

//__________________________________________________________________________
Int_t MpdTrackFinderIts5spd::RunKalmanFilterMod(MpdItsKalmanTrack *track, Int_t layBeg)
{
  /// Run Kalman filter in the modular geometry (might not work when propagating outward!!!)

  const Int_t maxBr = 24, maxBrLay = 2000; // max number of branches
  //const Int_t maxBr = 100, maxBrLay = 1000; // max number of branches
  const Double_t thick[9] = {0.005, 0.005, 0.005, 0.07, 0.07}; // layer thicknesses

  MpdKalmanGeoScheme *geo = MpdKalmanFilter::Instance()->GetGeo();

  //cout << fHits->GetEntriesFast() << endl;
  //Int_t layMax = ((MpdKalmanHit*)fKHits->Last())->GetLayer();
  //Int_t layMax = ((MpdKalmanHit*)fKHits->First())->GetLayer();
  //MpdKalmanHit *hitOK = 0x0;
  MpdKalmanHit hitTmp;
  //MpdKalmanTrack::TrackDir trackDir = track->GetDirection();
  //Int_t layBeg = 0, layEnd = -1, dLay = -1, layOK = -1;
  Int_t layEnd = -1, dLay = -1;//, layOK = -1;
//  if (trackDir == MpdKalmanTrack::kOutward) {
//    layEnd = layMax + 1;
//    dLay = 1;
//  }
  
  TMatrixDSym pointWeight(5), pointWeightTmp(5), saveWeight(5);
  TMatrixD param(5,1), paramTmp(5,1);

  Double_t quality[maxBrLay] = {0.0}; //,saveZ = 0.0, saveLeng = 0.0, dChi2Min = 0.0,  posNew = 0.0;
//  cout << " Starting hit: " << hitOK->GetLayer() << " " << hitOK->GetTrackID() << " " << hitOK->GetUsage() << endl;
  
  MpdItsKalmanTrack branch[maxBr];
  Int_t nBr = 1, nBr0 = 0, indx_lay[maxBrLay] = {0}, maxInLay = 0, ok = 0, nAdd = 0;

  branch[0] = *track;

  TString mass2 = "0.0194797849"; // pion mass squared
  /*
  if (fMCTracks) {
    // Get particle mass - ideal PID
    MpdMCTrack *mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(track->GetTrackID());
    TParticlePDG *pdgP = TDatabasePDG::Instance()->GetParticle(mctrack->GetPdgCode());
    if (pdgP) {
      Double_t mass = pdgP->Mass();
      if (mass < 0.1 || mass > 0.25) {
	// Electrons or heavier than pions
	mass2 = "";
	mass2 += mass*mass;
      }
    }
  }
  */

  // cout << "##### Cycle start ##### " << endl;
  //  cout << " layBeg, LayEnd: " << layBeg << " " << layEnd << endl;
 
  for (Int_t lay = layBeg; lay != layEnd; lay+=dLay) {
    //Int_t nLay = GetNofHitsInLayer(lay);
    //Int_t indx0 = GetHitsInLayer(lay);
    //MpdKalmanHit *hitMin = 0x0;
//    cout << " lay, nLay, indx: " << lay << " " << nLay << " " << indx0 << endl;
    //Int_t indxBeg = 0, indxEnd = nLay, dIndx = 1;
    MpdItsKalmanTrack branchLay[maxBrLay];

    nBr0 = nBr;
    nBr = 0;

//    cout << " Point#1: " << endl;
    for (Int_t ibr = 0; ibr < nBr0; ++ibr) {
      MpdItsKalmanTrack *curTr = &branch[ibr], *branchTr = nullptr;
      
      // Check for missing layer (2 sides)
      Int_t lastLay = 4, frozen = 0; //VK

      MpdKalmanHit *h = (MpdKalmanHit*) curTr->GetHits()->Last();
      if (h) lastLay = h->GetLayer();
      if (TMath::Abs(lay - lastLay) > 2) frozen = 1;

      ok = nAdd = 0;

      //Int_t firstHit = -1, skipHit = -1;
      //MpdItsKalmanTrack trackBr[4] = {*curTr, *curTr, *curTr, *curTr};
      map<Int_t,MpdItsKalmanTrack> trackBr;
      multimap<Int_t,Int_t> hitsInWin;
     
      Bool_t navOK = NavigateToLayer(lay, curTr, trackBr, hitsInWin); // navigate to layer

      //AZ if (!navOK || trackBr.size() == 0) { ok = -1; continue; }
      if (!navOK) { ok = -1; continue; }
      map<Int_t,MpdItsKalmanTrack>::iterator it;
      // Debug
      /*
      cout << " Forks: " << trackBr.size() << " " << hitsInWin.size() << endl;
      for (it = trackBr.begin(); it != trackBr.end(); ++it) {
	cout << " detID " << it->first << " " << fHitSts.GetLadder(it->first%1000000) << " " 
	     << fHitSts.GetSensor(it->first%1000000) << " " 
	     << it->second.GetParamNew(0) << " " << it->second.GetParamNew(1) << endl;
      }
      */

      for (it = trackBr.begin(); it != trackBr.end(); ++it) {
	//if (frozen) break;
	//Double_t leng = curTr->GetLength(), 
	Double_t step = 0;
	branchTr = &it->second;
	/*
	step = branchTr->GetLength() - leng;
	step = 0.005 * 2; // 50um
	if (step > 1.e-4) {
	  // Crossing silicon layer - add mult. scat. in the sensor
	  Double_t x0 = 9.36; // rad. length
	  TMatrixDSym *cov = branchTr->Weight2Cov();
	  Double_t th = branchTr->GetParamNew(3);
	  Double_t cosTh = TMath::Cos(th);
	  Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(branchTr, x0, step, mass2);
	  //cout << " Scat: " << hit->GetLayer() << " " << step << " " << TMath::Sqrt(angle2) << endl;
	  (*cov)(2,2) += (angle2 / cosTh / cosTh);
	  (*cov)(3,3) += angle2;
	  Int_t iok = 0;
	  MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
	  branchTr->SetWeight(*cov);
	}
	*/
	TVector3 mom3 = branchTr->Momentum3(), norm;
	mom3.SetMag(1.0);

	typedef multimap<Int_t,Int_t>::iterator mmit;
	pair<mmit,mmit> ret = hitsInWin.equal_range(it->first);
	map<Int_t,Int_t>::iterator it1;

	// Loop over hits in one sensor
	for (it1 = ret.first; it1 != ret.second; ++it1) {
	  MpdKalmanHit *hit = (MpdKalmanHit*) fKHits->UncheckedAt(it1->second);
	  if (it1 == ret.first) {
	    norm = geo->Normal(hit);
	    //step = 0.005 / TMath::Abs(norm * mom3) * 4.0; // extra factor 4. - possible overlaps 
	    step = thick[hit->GetLayer()] / TMath::Abs(norm * mom3) * 4.0; // extra factor 4. - possible overlaps 
	  }	       
	  // Exclude used hits
	  if (hit->GetFlag() & MpdKalmanHit::kUsed) continue;
	  if (frozen) continue;
	  // !!! Exact ID match
	  if (fExact && TrackID(hit) != track->GetTrackID()) continue;

	  //cout << lay << " " << hit->GetLayer() << " " << branchTr->GetParamNew(0) << " " << hit->GetMeas(0) << " "
	  //   << branchTr->GetParamNew(1) << " " << hit->GetMeas(1) << endl;
	  Double_t dChi2 = MpdKalmanFilter::Instance()->FilterHit(branchTr,hit,pointWeight,param);

	  //cout << lay << " " << dChi2 << " " << track->GetTrackID() << " " << TrackID(hit) << endl;

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
	      /*AZ 5-jan-2018
	      branchLay[nBr].SetLengAtHit(branchTr->GetLength()); 
	      branchLay[nBr].SetParamAtHit(param);
	      branchLay[nBr].SetWeightAtHit(*branchLay[nBr].GetWeight());
	      */
	      // Add multiple scattering in the sensor
	      //*
	      Double_t x0 = 9.36; // rad. length
	      TMatrixDSym *cov = branchLay[nBr].Weight2Cov();
	      Double_t th = branchLay[nBr].GetParamNew(3);
	      Double_t cosTh = TMath::Cos(th);
	      Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(&branchLay[nBr], x0, step, mass2);
	      //cout << " Scat: " << hit->GetLayer() << " " << step << " " << TMath::Sqrt(angle2) << endl;
	      (*cov)(2,2) += (angle2 / cosTh / cosTh);
	      (*cov)(3,3) += angle2;
	      Int_t iok = 0;
	      MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
	      branchLay[nBr].SetWeight(*cov);
	      //*/
	      ++nBr;
	      ++nAdd;
	    } 
	  }
	} // for (it1 = ret.first; it1 != ret.second; - loop over hits in one sensor
      } // for (it = trackBr.begin();  - loop over forks

      //cout << " Point#3: " << endl;
      // Add branch with missing layer
      if (nAdd == 0 && nBr < maxBrLay && ok != -1) {
	hitTmp.SetType(MpdKalmanHit::kFixedR);
	hitTmp.SetPos(fRad[lay]);
	if (!MpdKalmanFilter::Instance()->PropagateToHit(curTr,&hitTmp,kTRUE,kTRUE)) { }

	// Add multiple scattering in the sensor
	//Double_t x0 = 9.36, step = 0.005 * 4.0; // rad. length
	Double_t x0 = 9.36, step = thick[lay] * 4.0; // rad. length
	TMatrixDSym *cov = curTr->Weight2Cov();
	Double_t th = curTr->GetParamNew(3);
	Double_t cosTh = TMath::Cos(th);
	Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(curTr, x0, step, mass2);
	//cout << " Scat: " << hit->GetLayer() << " " << step << " " << TMath::Sqrt(angle2) << endl;
	(*cov)(2,2) += (angle2 / cosTh / cosTh);
	(*cov)(3,3) += angle2;
	Int_t iok = 0;
	MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
	curTr->SetWeight(*cov);

	trackBr[0] = *curTr;
	it = trackBr.begin();
	branchLay[nBr] = it->second;
	++it;
	// Take branch with greater length
	for ( ; it != trackBr.end(); ++it) {
	  if (it->second.GetLength() > branchLay[nBr].GetLength()) branchLay[nBr] = it->second;
	}
	++nBr;
      }

    } // for (Int_t ibr = 0; ibr < nBr0;
    maxInLay = TMath::Max (maxInLay,nBr);
    if (nBr == 0) {
      if (ok == -1) { 
	//cout << track->GetTrackID() << " " << lay << endl;
	//return ok; // stop track
	//cout << branch[0].GetNode() << ", " << branch[0].GetNodeNew() << endl;
        //cout << " Break:  layer = " << lay << " node = " << branch[0].GetNode() <<  endl;
	break;
      } 
      cout << " !!! MpdTrackFinderIts5spd::RunKalmanFilter: it can not happen |" << endl;
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

Bool_t MpdTrackFinderIts5spd::NavigateToLayer(Int_t lay, MpdItsKalmanTrack *curTr, std::map<Int_t,MpdItsKalmanTrack> &trackBr, 
					      std::multimap<Int_t,Int_t>& hitsInWin)
{
  ///< Navigate track to layer in modular geometry

  MpdItsKalmanTrack origTr(*curTr); // original track

  Double_t posNew = curTr->GetPosNew();
  TMatrixD parNew = *curTr->GetParamNew();
  TString nodeNew = curTr->GetNodeNew();
  //  cout << " Nodes: " << nodeNew << " " << curTr->GetNode() << " Positions: " << posNew << " " << curTr->GetPos() << endl;
  //  curTr->GetParamNew()->Print();
  //  curTr->GetParam()->Print();
  TString curPath = curTr->GetNode();
  MpdKalmanHit hitTmp;
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
  
  // Find hits near the extrapolated track position
  //Double_t window = 1.0; // +- 1 cm window around extrapolated track position
  Double_t window = 2.0; // +- 1 cm window around extrapolated track position
  Double_t rphiTr = curTr->GetParamNew(0), zTr = curTr->GetParamNew(1);
  set<Int_t> indsRphi, indsZ;
  multimap<Double_t,Int_t>::iterator mitb = fHitMapRphi[lay].lower_bound(rphiTr-window), mit1;
  multimap<Double_t,Int_t>::iterator mite = fHitMapRphi[lay].upper_bound(rphiTr+window);
  for (mit1 = mitb; mit1 != mite; ++mit1) indsRphi.insert(mit1->second);
  mitb = fHitMapZ[lay].lower_bound(zTr-window);
  mite = fHitMapZ[lay].upper_bound(zTr+window);
  for (mit1 = mitb; mit1 != mite; ++mit1) indsZ.insert(mit1->second);
  vector<Int_t> indv(indsRphi.size()+indsZ.size());
  vector<Int_t>::iterator vit = set_intersection (indsRphi.begin(), indsRphi.end(), indsZ.begin(), indsZ.end(), indv.begin());
  indv.resize(vit-indv.begin());
  Int_t nWin = indv.size();
  //cout << " Hits in window: " << nWin << " " << curTr->GetTrackID() << endl;
  if (nWin == 0) {
    //cout << " !!! Navigate: No hits in window " << endl;
    //return kFALSE; // !!! FIXME - handle missing hit
    return kTRUE; // !!! FIXME - handle missing hit
    //exit(0);
  }
  
  for (Int_t j = 0; j < nWin; ++j) {
    MpdKalmanHit *h = (MpdKalmanHit*) fKHits->UncheckedAt(indv[j]);
    //cout << h->GetDetectorID() << " " << endl;
    if (trackBr.find(h->GetDetectorID()) != trackBr.end()) {
      // Sensor has already been found
      hitsInWin.insert(pair<Int_t,Int_t>(h->GetDetectorID(),indv[j]));
      continue; 
    }
    MpdItsKalmanTrack tr(origTr);
    if (!MpdKalmanFilter::Instance()->PropagateToHit(&tr,h,kTRUE,kTRUE)) {
      //cout << " !!! Navigate: failed propagation " << endl;
      continue;
    }
    // Fork
    trackBr[h->GetDetectorID()] = tr;
    hitsInWin.insert(pair<Int_t,Int_t>(h->GetDetectorID(),indv[j]));
  }
  if (trackBr.size() > 14) {
    cout << "!!! Navigate: Too many hits in window " << endl;
    //exit(0);
  }
  
  return kTRUE;
}

//__________________________________________________________________________
Int_t MpdTrackFinderIts5spd::TrackID(MpdKalmanHit *hit)
{
  /// Return track ID of the hit

  FairHit *h = (FairHit*) fItsHits->UncheckedAt(hit->GetIndex());
  return ((MpdStsPoint*) fItsPoints->UncheckedAt(h->GetRefIndex()))->GetTrackID();
}

//__________________________________________________________________________
TVector2 MpdTrackFinderIts5spd::GetDistance(MpdKalmanTrack *track, MpdKalmanHit *hit)
{
  /// Compute distance between track and hit

  Int_t lay = hit->GetLayer();
  Int_t lay2 = lay / 2;
  //Int_t side = lay % 2;
  Int_t module0 = ((MpdItsHit5spd*) fItsHits->UncheckedAt(hit->GetIndex()))->Module();

  Double_t zTr = track->GetParamNew(1);
  Double_t zloc = zTr + fDz[lay2];
  Int_t modTr = Int_t (zloc / fZmod[lay2]);
  Int_t dMod = modTr - module0;

  Double_t dZ = 0;
  if (dMod != 0) {
    // Not in the same module0 - check Z-distance to module0 edge
    if (dMod > 0) dZ = zloc - (module0+1) * fZmod[lay2];
    else dZ = zloc - module0 * fZmod[lay2];
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
void MpdTrackFinderIts5spd::Write()
{
  /// Write

  TFile histoFile("ItsRec.root","RECREATE");
  Writedir2current(fHistoDir);
  histoFile.Close();
}

//__________________________________________________________________________
void MpdTrackFinderIts5spd::Writedir2current( TObject *obj )
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
void MpdTrackFinderIts5spd::AddHits()
{
  /// Add hit objects to tracks and compute number of wrongly assigned hits
  /// (hits with ID different from ID of starting TPC track)

  Int_t nFound = fTracks->GetEntriesFast();
  for (Int_t i = 0; i < nFound; ++i) {
    MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTracks->UncheckedAt(i);
//    cout << track->GetNode() << " " << track->GetNodeNew() << endl;
    Double_t c2 = track->GetChi2();
    track->SetChi2(track->GetChi2Its());
    track->SetChi2Its(c2);
    Int_t nHits = track->GetNofHits();
    //if (nHits == 0) { fTracks->RemoveAt(i); continue; }
    TClonesArray &trHits = *track->GetTrHits();
//    cout << nHits << " " << trHits.GetEntriesFast() << " " << track->GetTrackID() << endl;
    TObjArray *hits = track->GetHits();
    Int_t nWrong = 0, motherID = track->GetTrackID(); // nMirr = 0,
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
//      cout << " " << hit->GetLayer();
      MpdItsHit5spd *h = (MpdItsHit5spd*) fItsHits->UncheckedAt(hit->GetIndex());

      Int_t motherID1 = ((FairMCPoint*) fItsPoints->UncheckedAt(h->GetRefIndex()))->GetTrackID();
//      cout << "-" << motherID1;
      // Get point mother ID 
      mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(motherID1);
      while (mctrack->GetMotherId() >= 0) {
        motherID1 = mctrack->GetMotherId();
        mctrack = (MpdMCTrack*) fMCTracks->UncheckedAt(mctrack->GetMotherId());
      }
      if (motherID1 != motherID) ++nWrong;

    }
//    if (nHits) cout << "\n" << nWrong << endl;
    track->SetNofWrong(nWrong);
    //MpdKalmanTrack *tpc = (MpdKalmanTrack*) fTpcTracks->UncheckedAt(track->GetUniqueID()-1);
    //track->SetChi2(track->GetChi2()+tpc->GetChi2());
    //track->SetLastLay();
    //track->GetParam()->Print();
    track->SetNofHits(track->GetNofTrHits()); // TPC and ITS hits
    track->SetNofIts(nHits);
//    cout << nHits << " " << track->GetNofTrHits() << " " << track->GetTrackID() << " "  << track->GetChi2Its() << " " << track->GetChi2() << endl;
//   cout << " " << endl;
  }
  fTracks->Compress();
}

//__________________________________________________________________________

Bool_t MpdTrackFinderIts5spd::Refit(MpdItsKalmanTrack *track, Double_t mass, Int_t charge)
{
  /// Refit track in ITS+TPC using track hits (toward beam line) for some
  /// particle mass and charge hypothesis 

  if (fTpcKF == nullptr) fTpcKF = (MpdTpcKalmanFilter*) FairRun::Instance()->GetTask("TPC Kalman filter");
  if (fTpcKF == nullptr) {
    cout << " !!! TPC Kalman Filter was not activated !!! " << endl;
    exit(1);
  }
  MpdKalmanGeoScheme *geoScheme = MpdKalmanFilter::Instance()->GetGeo();

  //MpdTpcKalmanTrack *tpc = (MpdTpcKalmanTrack*) fTpcTracks->UncheckedAt(track->GetUniqueID()-1);
  track->SetUniqueID(0); // reset ITS flag
  //tr.GetWeightAtHit()->Print();
  //tr.GetParamAtHit()->Print();
  //AZ fTpcKF->Refit(track, mass, charge); // this mixes up hit order (ITS and TPC) !!!
  if (track->GetNofTrHits() > track->GetNofIts()) fTpcKF->Refit(track, mass, charge); // this mixes up hit order (ITS and TPC) !!!
  //tr.GetWeight()->Print();
  //tr.GetParam()->Print();
  if (track->GetNofIts() == 0) return kTRUE;
  
  MpdKalmanHit hTmp;
  hTmp.SetType(MpdKalmanHit::kFixedR);
  hTmp.SetPos(track->GetPos()); 
  track->SetParamNew(*track->GetParam());
  track->SetPos(track->GetPosNew());
  track->ReSetWeight();
  //TMatrixDSym w = *tr.GetWeight(); // save current weight matrix
  //MpdKalmanFilter::Instance()->PropagateToHit(track,&hit,kFALSE);
  MpdKalmanFilter::Instance()->PropagateParamToHit(track,&hTmp,kFALSE);
  //tr.SetWeight(w); // restore original weight matrix (near TPC inner shell)
  track->SetDirection(MpdKalmanTrack::kInward);
  
  track->SetUniqueID(1);
  Int_t ntot = track->GetNofTrHits(), nhits = TMath::Min(15,ntot); 
  Int_t indx0 = ntot - nhits; 

  TMatrixD param(5,1);
  TMatrixDSym weight(5), pointWeight(5);
  TString mass2 = "";
  mass2 += mass * mass;
  Double_t vert[3] = {0.0,0.0,0.0};
  const Double_t thick[9] = {0.005, 0.005, 0.005, 0.07, 0.07};

  for (Int_t ih = indx0; ih < ntot; ++ih) {
    MpdKalmanHit *hit = (MpdKalmanHit*) track->GetTrHits()->UncheckedAt(ih);
    if (hit->GetUniqueID() == 0) continue; // skip TPC hits
    MpdKalmanFilter::Instance()->PropagateToHit(track, hit, kFALSE, kTRUE); // do not adjust track length for the moment
    Double_t dChi2 = MpdKalmanFilter::Instance()->FilterHit(track, hit, pointWeight, param);
    //cout << ih << " " << hit->GetPos() << " " << hit->GetUniqueID() << " " << dChi2 << "\n";

    track->SetChi2(track->GetChi2()+dChi2);
    weight = *track->GetWeight();
    weight += pointWeight;
    track->SetWeight(weight);
    track->SetParamNew(param);
    
    // Add multiple scattering in the sensor
    //*
    Double_t x0 = 9.36; //, step = 0.005 * 4.0; // rad. length
    TVector3 mom3 = track->Momentum3();
    mom3.SetMag(1.0);
    TVector3 norm = geoScheme->Normal(hit);
    Double_t cosPhi = norm.Dot(mom3);

    TMatrixDSym *cov = track->Weight2Cov();
    Double_t th = track->GetParamNew(3);
    Double_t cosTh = TMath::Cos(th);
    //Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, x0, thick[lay]/cosPhistep, mass2);
    Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, x0, 4*thick[hit->GetLayer()]/cosPhi, mass2);
    (*cov)(2,2) += (angle2 / cosTh / cosTh);
    (*cov)(3,3) += angle2;
    Int_t iok = 0;
    MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
    track->SetWeight(*cov);
    //*/
  }
  
  //------ Propagate track to the beam line
  track->SetParam(*track->GetParamNew());
  track->SetPos(track->GetPosNew());
  Double_t pos = track->GetPos();
  TMatrixD par = *track->GetParam();
  TMatrixDSym cov = *track->Weight2Cov();
  Double_t leng = track->GetLength();
  TString nodeNew = track->GetNodeNew();
  //cout << " 1: " << nodeNew << ", " << track->GetNode() << " " << track->GetHits()->GetEntriesFast() << endl;
 
  // Go to beam pipe
  hTmp.SetPos(fPipeR);
  //hTmp.SetPos(2.9); // !!! FIXME
  Int_t iok = MpdKalmanFilter::Instance()->PropagateToHit(track, &hTmp, kFALSE);
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
  track->Weight2Cov(); // rebuild covar matrix (it was overwritten by weight above)

  hTmp.SetPos(0.);
  hTmp.SetMeas(0,track->GetParam(2)); // track Phi
  iok = MpdKalmanFilter::Instance()->PropagateToHit(track, &hTmp, kFALSE);
  if (iok != 1) MpdKalmanFilter::Instance()->FindPca(track, vert);
  else track->SetNodeNew(""); //??? 3-jan-2018
  track->SetParam(*track->GetParamNew()); // !!! track params at PCA
  return kTRUE;
}

//__________________________________________________________________________

ClassImp(MpdTrackFinderIts5spd);