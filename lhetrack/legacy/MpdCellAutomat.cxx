// -------------------------------------------------------------------------
// -----                                    source file                -----
// -----                 Created 13/06/2014  by A. Zinchenko           -----
// -------------------------------------------------------------------------

/**  MpdCellAutomat.cxx
 *@author M.Strelchenko, A.Zinchenko <Alexander.Zinchenko@jinr.ru>
 **
 ** Track finder in MPD Inner Tracking System (ITS) using cellular automaton (CA)
 **/

#include "MpdStsGeoPar.h"
#include "MpdCellAutomat.h" 
#include "MpdCellTrack.h" 
#include "MpdCodeTimer.h"
#include "MpdConstField.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanTrack.h"
#include "MpdStsHit.h"
#include "MpdStsPoint.h"
#include "MpdItsKalmanTrack.h"
//#include "MpdTpcKalmanTrack.h"
//#include "MpdTpcKalmanFilter.h"

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
#include "TVector2.h"
#include "TVector3.h" 
#include "TClonesArray.h"
#include <TRandom.h>
#include <TArrayI.h>

#include <iostream>
#include <map>
#include <set>

using std::cout;
using std::endl;
//using std::map;

const Double_t MpdCellAutomat::fgkChi2Cut = 20; //20; //100; 

FILE *lunErr = 0x0; //fopen("error1.dat","w");
//FILE* lunErr= fopen ("file.txt","w");
//__________________________________________________________________________
MpdCellAutomat::MpdCellAutomat(const char *name, Int_t iVerbose ) 
 :FairTask(name, iVerbose),
  fExact(0),
  fNTracks(0),
  fNPass(2),
  fGeo(0)
{
  fKHits1 = new TClonesArray("MpdKalmanHit", 100);
  for (Int_t i = 0; i < 4; ++i) {
    fKHits[i] = new TClonesArray("MpdCellTrack", 100);
    f2DHits[i] = new TClonesArray("MpdCellTrack", 100);
  }
  
  //fKHits[0] = f2DHits[0]; //work
  fTracks = new TClonesArray("MpdItsKalmanTrack", 100);
  fTrackCand = new TClonesArray("MpdItsKalmanTrack", 100);
  fHistoDir = 0x0;
  fhLays = new TH1F("hLaysITS","ITS layers",10,0,10);
  fLayPointers = 0x0;
}


//__________________________________________________________________________
MpdCellAutomat::~MpdCellAutomat()
{
  //delete fKHits;
  //delete fTracks;
  //delete fTrackCand;
  delete [] fLayPointers;
  delete fhLays;
}

//__________________________________________________________________________
InitStatus MpdCellAutomat::Init() 
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
  Int_t nLay = 4;
  Double_t size = 6.2;
  for (Int_t i = 0; i < nLay; ++i) {
    volName[5] = 97 + i; // 'a', 'b', ..
    FairGeoNode* sensVol = (FairGeoNode*) (sensNodes->FindObject(volName));
    if (sensVol == 0x0) {
      fGeo = 1; // modular geometry
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
	Double_t xyzL[3] = {0}, xyzM[3];
	gGeoManager->LocalToMaster(xyzL,xyzM);
	Double_t rad = TMath::Sqrt (xyzM[0] * xyzM[0] + xyzM[1] * xyzM[1]);
	fRad[2*i] = TMath::Min (fRad[2*i],rad);
      }
      fRad[2*i+1] = fRad[2*i];
      TGeoVolume *ladd = gGeoManager->GetVolume(volName);
      TGeoBBox *box = (TGeoBBox*) ladd->GetShape();
      //safety = -box->GetDY();
      safety = box->GetDY();
      fRad[2*i] -= safety;
      fRad[2*i+1] -= safety;
      if (i == 0) { fRad[2*i] -= safety; fRad[2*i+1] -= safety; }
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

void MpdCellAutomat::FillGeoScheme()
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
InitStatus MpdCellAutomat::ReInit() 
{
  fItsPoints = (TClonesArray *) FairRootManager::Instance()->GetObject("StsPoint");
  fItsHits =(TClonesArray *) FairRootManager::Instance()->GetObject("StsHit");
  if (fItsPoints == 0x0 || fItsHits == 0x0) return kERROR;
  fTpcTracks =(TClonesArray *) FairRootManager::Instance()->GetObject("TpcKalmanTrack");
  fMCTracks =(TClonesArray *) FairRootManager::Instance()->GetObject("MCTrack");
  //fSTSTrackMatch = (TClonesArray*) FairRootManager::Instance()->GetObject("STSTrackMatch");
  //fPrimVtx =  (FairVertex *) FairRootManager::Instance() ->GetObject("PrimaryVertex");

  FairRootManager::Instance()->Register("ItsTrack", "Its", fTracks, kTRUE);
  FairRootManager::Instance()->Register("ItsTrackCand", "ItsCand", fTrackCand, kTRUE);
  FairRootManager::Instance()->Register("CellTrack0", "Its0", fKHits[0], kTRUE);
  FairRootManager::Instance()->Register("CellTrack1", "Its1", fKHits[1], kTRUE);
  FairRootManager::Instance()->Register("CellTrack2", "Its2", fKHits[2], kTRUE);
  FairRootManager::Instance()->Register("CellTrack3", "Its3", fKHits[3], kTRUE);

  FairRootManager::Instance()->Register("CellTrack_0", "Its0_2D", f2DHits[0], kTRUE); // 21.04
  FairRootManager::Instance()->Register("CellTrack_1", "Its1_2D", f2DHits[1], kTRUE); // 21.04
  FairRootManager::Instance()->Register("CellTrack_2", "Its2_2D", f2DHits[2], kTRUE); // 21.04
  FairRootManager::Instance()->Register("CellTrack_3", "Its3_2D", f2DHits[3], kTRUE); // 21.04
  // fNPass = 2; // 2 prochoda
  //fNPass = 2;
  fNPass = 3;
  return kSUCCESS;
}

//__________________________________________________________________________
void MpdCellAutomat::Reset() 
{
  
  cout << " MpdCellAutomat::Reset  " << endl; 
  //fKHits->Clear(); 
  for (Int_t i = 0; i < 4; ++i) {
    fKHits[i]->Delete();
    f2DHits[i]->Delete();
  }
  fKHits1->Delete();
  fTracks->Delete();
  fTrackCand->Delete();
  delete [] fLayPointers;
  fLayPointers = NULL;
  fCellMap.clear();
}

//__________________________________________________________________________
void MpdCellAutomat::SetParContainers() 
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
void MpdCellAutomat::Finish() 
{
  //Write();
}

//__________________________________________________________________________
void MpdCellAutomat::Exec(Option_t * option) 
{
  if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(),__FUNCTION__);
  static int eventCounter = 0;    
  cout << " - - - - \n ItsRec event " << ++eventCounter << endl;

  Reset();

  // Create Kalman hits
  if (fItsHits->GetEntriesFast() == 0) return;
  //MakeTrackCandidates(); // 21.02
  //ExtendCellTracks(); // 21.02
  //return;
  Build2DHits();
  MakeKalmanHits();

  for (Int_t i = 0; i < fNPass; ++i) {
    fTrackCand->Delete(); // 05.03
    fKHits[0]->Delete();// MpdCellTrack
    fKHits[1]->Delete();
    fKHits[2]->Delete();
    fKHits[3]->Delete();
   
    MakeTrackCandidates(i);// 21.02
    ExtendCellTracks(i);// 21.02
    GetTrackSeeds(i);//attention!!
    /*
      Int_t nHitsOk = 0, nHits = fKHits->GetEntriesFast();
      for (Int_t j = 0; j < nHits; ++j){
      MpdKalmanHit *hit = (MpdKalmanHit*) fKHits->Unchecked(j);
      if (hit->GetFlag() >= 0) ++nHitsOk;
      }
     */
    // cout << "  Total number of hits for tracking: " << fKHits->GetEntriesFast() << endl;
    cout << "  Total number of track seeds: " << fTracks->GetEntriesFast() << endl;
    cout << "  Total number of track seeds fTrackCand: " << fTrackCand->GetEntriesFast() << endl;
    DoTracking(i); //attention!!!
    RemoveDoubles();
    StoreTracks();
    cout << "  Total number of found tracks: " << fTrackCand->GetEntriesFast() << endl;

    //if (i != fNPass - 1) ExcludeHits(); // exclude used hits   work 
    ExcludeHits(); // exclude used hits  test
  }

  GetShortTracks();
  AddHits(); // add hit objects to tracks
  cout << "  Total number of found tracks: " << fTrackCand->GetEntriesFast() << endl;
  if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(),__FUNCTION__);
}

//__________________________________________________________________________
void MpdCellAutomat::MakeKalmanHits()
{
  /// Create Kalman hits from ITS hits.

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
    //gRandom->Rannor(dX,dZ);
    //if (errZ > 2) dZ = 0.0; // 1-D case
    dZ = 0.0; // 1-D case
    //Double_t meas[2] = {xloc+dX*errX, z+dZ*errZ};
    Double_t meas[2] = {xloc+h->GetDx()*errX, z+dZ*errZ};
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
				    new ((*fKHits1)[nKH++]) MpdKalmanHit(lay*1000000+fId2Id[lay][h->GetDetectorID()], 1, 
									 MpdKalmanHit::kFixedP, meas, err, cossin, 0., r, ih);
    // Mask out sector number - sensor layout
    else if (fGeo) hit = new ((*fKHits1)[nKH++]) MpdKalmanHit(lay*1000000+(h->GetDetectorID()&((2<<12)-1)), 1, 
                                                MpdKalmanHit::kFixedP, meas, err, cossin, 0., r, ih);
    else hit = new ((*fKHits1)[nKH++]) MpdKalmanHit(lay*1000000+nKH-1, 1, MpdKalmanHit::kFixedR, 
						    meas, err, cossin, 0., r, ih);
    hit->SetUniqueID(0);
    // Add second measurement - just for test at the moment
    //!!!
    //hit->SetNofDim(2);
    //!!!
    layMax = TMath::Max (lay, layMax);
  }
  cout << " Max layer = " << layMax << " " << fKHits1->GetEntriesFast() << endl;
  //fKHits1->Sort(); // in descending order in R
  //cout << ((MpdKalmanHit*)fKHits->UncheckedAt(0))->GetPos() << endl;
  cout << ((MpdKalmanHit*)fKHits1->UncheckedAt(0))->GetDist() << endl;
  
}

//_________________________________________________________________________
void MpdCellAutomat::Build2DHits()
{
  /// Build track candidates from ITS 2-D hits.

  fItsHits->Sort();
  Int_t nHits = fItsHits->GetEntriesFast(), layMax = 0, lay = 0, nKH = 0;
  //Double_t  errZ = 0.012, errX = 0.0023; // 120um in Z, 23um in R-Phi (local X)
  Double_t  errZ = 0.025, errX = 0.0023; // 250um in Z, 23um in R-Phi (local X)
  //Double_t r, z, xloc, errZ = 0.12, errX = 0.0023; // 1.2mm in Z, 23um in R-Phi (local X)
  //Double_t r, z, xloc, errZ = 50.0, errX = 0.01; // 50cm in Z, 100um in R-Phi (local X)
  Double_t xloc[2];
  Double_t meas[2];
  Double_t r;
  Double_t dX = 0, dZ = 0;

  //gRandom->Rannor(dX,dZ);
 
    Int_t lay0 = 1; // layer 1
    fLayBeg[0] = 0; // index for layer 1
    for (Int_t i1 = 0; i1 < nHits; ++i1) {
      MpdStsHit *h = (MpdStsHit*) fItsHits->UncheckedAt(i1);
      Int_t lay = h->Layer();
      if (lay0 != lay) { 
	fLayBeg[lay - 1] = i1;
	lay0 = lay;
      }
    } 
  
  for (Int_t i234 = 0; i234 < 4; ++i234) {  // Loop over layers 1-4
    Int_t nKH = 0, iprint = 0;
    printf("i234 =%d\n", i234);
    printf("h->lay=%d\n",lay);
    Int_t hEnd = nHits;
    if (i234 != 3) hEnd = fLayBeg[i234 + 1];

    for (Int_t ih = fLayBeg[i234]; ih < hEnd; ++ih) {
      // for (Int_t ih = fLayBeg[0]; ih < fLayBeg[1]; ++ih) { // work 28.04
      // for (Int_t ih = 0; ih < nHits; ++ih) { // old
      MpdStsHit *h = (MpdStsHit*) fItsHits->UncheckedAt(ih);
      if (h->GetFlag() < 0) continue;
      lay = (h->Layer() - 1) * 2 + h->Side();

      r = TMath::Sqrt (h->GetX() * h->GetX() + h->GetY() * h->GetY());
      // printf( "r=%d\n",r);
      //gRandom->Rannor(dX,dZ);
      xloc[0]=h->GetLocalX()+h->GetDx()*errX;
      //cout << ih << " " << h->Layer()-1 << endl;
      // Double_t dX = 0, dZ = 0;
      // gRandom->Rannor(dX,dZ);
      // dZ = 0.0;
      Double_t err[2] = {errX, errZ};
      Double_t cossin[2] = {TMath::Cos(fStereoA[h->Side()]), TMath::Sin(fStereoA[h->Side()])}; //23.08    
      MpdCellTrack *hit = 0x0;//ms
      MpdStsHit *hsts = h;
      if (h->Side()) hsts = NULL;

      for (Int_t ih1 = ih + 1; ih1 < hEnd; ++ih1) {
	// printf("Dim:xloc[0]=%f,xloc[1]=%f,nHits=%d\n",xloc[0],xloc[1],nHits);
	MpdStsHit *h1 = (MpdStsHit*) fItsHits->UncheckedAt(ih1);
	if (h1->GetFlag() < 0) continue;
	if (h->Side() == h1->Side()) continue;
	if (h->Layer() != h1->Layer()) continue;
	if (h->Ladder() != h1->Ladder()) continue;
	if (h->Detector() != h1->Detector()) continue;
	//gRandom->Rannor(dX,dZ);
	xloc[1]=h1->GetLocalX()+h1->GetDx()*errX;
	Double_t cossinn[2] = {TMath::Cos(fStereoA[h1->Side()]), TMath::Sin(fStereoA[h1->Side()])};
	if (hsts == NULL) hsts = h1;
	//Make 2-dim coordinate
	// Double_t x1=((xloc1[1]*sterSin[0]-xloc[0]*sterSin[1])/(sterCos[1]*sterSin[0]-sterCos[0]*sterSin[1]));
	Double_t x1=((xloc[1]*cossin[1]-xloc[0]*cossinn[1])/(cossinn[0]*cossin[1]-cossin[0]*cossinn[1]));
	//Double_t z1=((xloc[0]-(x1)*sterCos[0])/sterSin[0]);
	Double_t z1 = TMath::Abs(cossin[1]) > 1.e-6 ? (xloc[0] - x1 * cossin[0]) / cossin[1] : 
	  (xloc[1] - x1 * cossinn[0]) / cossinn[1];

	MpdKalmanHit hitTmp;
	//Int_t lay1 = (h1->Layer() - 1) * 2 + hsts->Side();
	Int_t lay1 = (hsts->Layer() - 1) * 2 + hsts->Side();
	hitTmp.SetDetectorID(lay1*1000000+fId2Id[lay1][hsts->GetDetectorID()]);
	Double_t sizeZ = MpdKalmanFilter::Instance()->GetGeo()->Size(&hitTmp).Y();
	if (TMath::Abs(z1) > sizeZ) continue;
	// printf("sizeZ=%f %d %d %d %d %d %d %d %d %f %f\n",sizeZ,h->Layer(),h1->Layer(),h->Ladder(),h1->Ladder(),h->Detector(),h1->Detector(),h->GetTrackID(),h1->GetTrackID(),x1,z1);
	TString path = MpdKalmanFilter::Instance()->GetGeo()->Path(hitTmp.GetDetectorID()); 
	meas[0] = x1;
	meas[1] = z1;
	gGeoManager->cd(path);
	Double_t v3[3] = {meas[0], 0.0, meas[1]}, v31[3] = {0.0}; 
	gGeoManager->LocalToMaster(v3, v31); 
	if (h->Layer() == 1 && TMath::Abs(v31[2]) > 18.0) continue; // Eta < 1.99
	else if ( h->Layer() == 2 && TMath::Abs(v31[2]) > 40.0) continue; // rad 11  
	else if ( h->Layer() == 3 && TMath::Abs(v31[2]) > 61.0) continue; //  rad 17
	else if ( h->Layer() == 4 && TMath::Abs(v31[2]) > 83.0) continue; //  rad 23

	TVector3 point1(v31);

	Float_t x1y1[3] = {v31[0],v31[1],v31[2]};
	Float_t r1 = TMath::Sqrt(x1y1[1] * x1y1[1] + x1y1[0] * x1y1[0]);
	Float_t angL = TMath::ATan(x1y1[2] / r1); //absoliyt angle
	    
	/*
      ////printf("angL=%f\n",angL);
      //(Int_t detID, Int_t nDim, HitType hitType, Double_t *meas, Double_t *err, Double_t *cosSin, Double_t signal, Double_t dist, Int_t index)
      
      //MpdKalmanStripHit *hit = new ((*fKHits)[nKH++]) MpdKalmanStripHit(r, fStereoA[h->Side()],
      //			  xloc+dX*errX, z+dZ*errZ, errX, errZ, lay, ih); 
      //hit->SetType(MpdKalmanHit::kFixedR);
      //cout << h->GetDetectorID() << " " << fId2Id[lay][h->GetDetectorID()] << endl;
      //if (fGeo && h->GetUniqueID()) hit = new ((*fKHits)[nKH++]) MpdKalmanHit(lay*1000000+h->GetDetectorID(), 1, 
      // if (fGeo && h->GetUniqueID()) hit = new ((*fKHits)[nKH++]) MpdKalmanHit(lay1*1000000+fId2Id[lay1][hsts->GetDetectorID()], 1, MpdKalmanHit::kFixedP,meas, err, cossin, 0., r, ih);//ms
      */
	
	if (fGeo && h->GetUniqueID()) hit = new ((*f2DHits[h->Layer()-1])[nKH++]) MpdCellTrack(lay1*1000000+fId2Id[lay1][hsts->GetDetectorID()], 1, MpdCellTrack::kFixedP,point1, err,
											      cossin, 0., r, ih, ih1, -1); //ms 17.06
	// Mask out sector number - sensor layout
	// else if (fGeo) hit = new ((*fKHits)[nKH++]) MpdKalmanHit(lay*1000000+(h->GetDetectorID()&((2<<12)-1)), 1, MpdKalmanHit::kFixedP,meas, err, cossin, 0., r, ih);//ms
	// else hit = new ((*fKHits)[nKH++]) MpdKalmanHit(lay*1000000+nKH-1, 1, MpdKalmanHit::kFixedR,meas, err, cossin, 0., r, ih); //ms
	else if (fGeo) hit = new ((*f2DHits[h->Layer()-1])[nKH++]) MpdCellTrack(lay*1000000+(h->GetDetectorID()&((2<<12)-1)), 1, MpdCellTrack::kFixedP,point1, err, cossin, 0., r, ih, ih1, -1); //ms 17.06
	else hit = new ((*f2DHits[h->Layer()-1])[nKH++]) MpdCellTrack(lay*1000000+nKH-1, 1, MpdCellTrack::kFixedR,point1, err, cossin, 0., r, ih, ih1, -1); //ms 17.06
	
	hit->SetUniqueID(0);
	// Add second measurement - just for test at the moment
	hit->SetNofDim(2);
	// hit->SetIndex(ih1);
	hit->SetCosSin(1,angL);
	layMax = TMath::Max (lay1, layMax);
      } // for (Int_t ih1 = 0; ih1 < nHits; ++ih1)
    } // for (Int_t ih = 0; ih < nHits; ++ih)
    printf("Schetchik:nKH=%d\n",nKH);

    
  } //   for (Int_t i234 = 0; i234 < 4; ++i234)

  //printf("fLayBeg[1]=%d\n",fLayBeg[0]);
  printf("fLayBeg[1]=%d\n",fLayBeg[1]);
  printf("fLayBeg[2]=%d\n",fLayBeg[2]);
  printf("fLayBeg[3]=%d\n",fLayBeg[3]);
  cout << "------------|                             |--------------" << endl;
  cout << "------------------ STOP Build2DHits ---------------------" << endl;
  //cout << " Max layer = " << layMax << " " << fKHits->GetEntriesFast() << endl;//ms 06.05
  /*
  for (Int_t i = 0; i < 4; ++i) {
    fKHits[i]->Sort();
  }
  */
  //fKHits->Sort(); // in descending order in R
  //cout << ((MpdKalmanHit*)fKHits->UncheckedAt(0))->GetPos() << endl;
  //cout << ((MpdCellTrack*)fKHits->UncheckedAt(0))->GetDist() << endl; //ms 06.05
}

//_________________________________________________________________________
void MpdCellAutomat::MakeTrackCandidates(Int_t iPass)
{
  // Read 2D Hits and create trackCandidates from 2D Hits at first layer

  Int_t nKH = 0, nHits = f2DHits[0]->GetEntriesFast();
  
  for (Int_t i = 0; i < nHits; ++i) {
    MpdCellTrack *cellTr = (MpdCellTrack*) f2DHits[0]->UncheckedAt(i);
    MpdStsHit *h = (MpdStsHit*) fItsHits->UncheckedAt(cellTr->GetIndex(0));
    if (h->GetFlag() < 0) continue;
    h = (MpdStsHit*) fItsHits->UncheckedAt(cellTr->GetIndex(1));
    if (h->GetFlag() < 0)continue;
    
    if (fGeo && h->GetUniqueID()) cellTr = new ((*fKHits[0])[nKH++]) MpdCellTrack(*cellTr); /// work!!
    else if (fGeo) cellTr = new ((*fKHits[0])[nKH++]) MpdCellTrack(*cellTr); //
    else cellTr = new ((*fKHits[0])[nKH++]) MpdCellTrack(*cellTr); //

    cellTr->SetCode(i);
    
  } // for (Int_t i = 0; ih < nHits; ++i)
  
  printf("Schetchik:nKH=%d\n",nKH);
  cout << "-------STOP MakeTrackCand---------" << endl;
  
}

//_________________________________________________________________________

// relizz

void MpdCellAutomat::ExtendCellTracks(Int_t iPass)
{
  /// Extend cell tracks to layers 2-4

  const Int_t nSigm = 3;
  Double_t  errZ = 0.025, errX = 0.0023; // 250um in Z, 23um in R-Phi (local X)
  
  for (Int_t i234 = 1; i234 < 4; ++i234) {  // Loop over layers 1-3
    Int_t nKH = 0, iprint = 0;
    printf("i234 =%d\n", i234);

    Int_t nHits2 = f2DHits[i234]->GetEntriesFast(); //26.05
      
    for (Int_t i = 0; i < nHits2; ++i) {
     
      MpdCellTrack *hit = (MpdCellTrack*) f2DHits[i234]->UncheckedAt(i);
      MpdStsHit *h = (MpdStsHit*) fItsHits->UncheckedAt(hit->GetIndex(0));
      if (h->GetFlag()<0 ) continue; //work
      h = (MpdStsHit*) fItsHits->UncheckedAt(hit->GetIndex(1));
      if (h->GetFlag()<0 )continue;
    
      TVector3 v311 = hit->GetMeas();

      Int_t nTracks = fKHits[i234-1]->GetEntriesFast(); //prediduchii sloi fkHits ->massiv ykazateley// work version 6.06.2014
      
      if (iprint == 0) { ++iprint; cout << "nTracks= " << nTracks << endl; }
  
      for (Int_t itr = 0; itr < nTracks; ++itr){
	Float_t angT, angL, deltaZ, dangL; //ang1->angT(transverse-"poperechnii" ); ang2->angL(longitudinal-"prodolnii");
	MpdCellTrack *track = (MpdCellTrack*) fKHits[i234-1]->UncheckedAt(itr); // work version 6.06.2014
	/// printf ("CosSin(0)1111=%f\n",track->GetCosSin(0));

	if ( track->GetPrevTrack() < 0) {
	  // Layer 2
	  TVector3 v1 = track->GetMeas();
	  Float_t x1y1[3] = {v311[0] - v1[0], v311[1] - v1[1], v311[2] - v1[2]}; //vector from layer 1 to lay 2
	  Float_t r1 = TMath::Sqrt(x1y1[1] * x1y1[1] + x1y1[0] * x1y1[0]); 
	  
	  angL = TMath::ATan (x1y1[2] / r1);
	  dangL = angL - track->GetCosSin(1); 
	  //cout << "dangL: " << dangL << endl;
	  //  if (TMath::Abs(dangL) > 0.047 * nSigm) continue;
	  if (TMath::Abs(dangL) > 0.045 * nSigm) continue; // sigma lay == 2 pt 0.05 Gev 
	  Float_t phi1 = TMath::ATan2(x1y1[1],x1y1[0]), phi0 = TMath::ATan2(v1[1],v1[0]);
	  angT = MpdKalmanFilter::Instance()->Proxim(phi0, phi1) - phi0;
	  //cout << "angT: " << r1 << " " << angT << endl;

	  if ( iPass == 0 && TMath::Abs(angT) > 0.023) continue; //Pt > 0.4 GeV/c
	  else if ( iPass == 1 && TMath::Abs(angT) > 0.043) continue; // Pt > 0.2 GeV/c
	  else if ( iPass == 2 && TMath::Abs(angT) > 0.18) continue; // Pt > 0.05 GeV/c

	  if (TMath::Abs(dangL) > nSigm * Interp(TMath::Abs(angT),0,i234-1)) continue;
	  
	} else {
	  
	  MpdCellTrack *track1 = (MpdCellTrack*) fKHits[i234-2]->UncheckedAt(track->GetPrevTrack());// trek s previos layer // work version 6.06.2014

	  TVector3 v1 = track->GetMeas();
	  TVector3 v5 = (track->GetMeas() - track1->GetMeas()); // track to previous layer
	  
	  Float_t x1y1[3] = {v311[0] - v1.X(), v311[1] - v1.Y(), v311[2]- v1.Z()}; //vector to this layer
	  Float_t r1 = TMath::Sqrt(x1y1[1] * x1y1[1] + x1y1[0] * x1y1[0]); 
	  
	  angL = TMath::ATan (x1y1[2] / r1);
	  //AZ - rough estimate
	  if (TMath::Abs(angL - track->GetCosSin(1)) > 0.045 * nSigm) continue; // sigma lay == 2 pt 0.05 Gev 

	  Double_t angLmean = angL;
	  if (h->Layer() == 3) {
	    // Take average of 2 values
	    angLmean = (angLmean + track->GetCosSin(1)) / 2;
	  } else {
	    // Take average of 3 values
	    angLmean = (angLmean + 2 * track->GetCosSin(1)) / 3;
	  }
	  //dangL = angL - track->GetCosSin(1);
	  dangL = angLmean - track->GetCosSin(1);
	  //cout << "dangL: " << dangL << endl;
	  
	  if (h->Layer() == 3 && TMath::Abs(dangL) > 0.02 * nSigm) continue; // sigma lay == 3 pt 0.05 GeV
	  else if (h->Layer() == 4 && TMath::Abs(dangL) > 0.016 * nSigm) continue; // sigma lay == 4 pt 0.05 GeV
	  Float_t phi1 = TMath::ATan2(x1y1[1],x1y1[0]), phi0 = TMath::ATan2(v5[1],v5[0]);
	  angT = MpdKalmanFilter::Instance()->Proxim(phi0, phi1) - phi0;

	  //AZ - rough estimate
	  if ( iPass == 0 && TMath::Abs(angT) > 0.023) continue; //Pt > 0.4 GeV/c
	  else if ( iPass == 1 && TMath::Abs(angT) > 0.043) continue; // Pt > 0.2 GeV/c
	  else if ( iPass == 2 && TMath::Abs(angT) > 0.18) continue; // Pt > 0.05 GeV/c
	  
	  Double_t angTmean = angT;
	    
	  if (h->Layer() == 3) {
	    // Take average of 2 values
	    angTmean = (angTmean + track->GetCosSin(0)) / 2;
	  } else {
	    // Take average of 3 values
	    angTmean = (angTmean + 2 * track->GetCosSin(0)) / 3;
	  }
	  //cout << "dangT: " << angTmean - track->GetCosSin(0) << endl;
	  //if (TMath::Abs(angT - track->GetCosSin(0)) > nSigm * Interp(TMath::Abs(angTmean),1,h->Layer()-3)) continue;///dangt 
	  if (TMath::Abs(angTmean - track->GetCosSin(0)) > nSigm * Interp(TMath::Abs(angTmean),1,i234-2)) continue;///dangt 
	  if (TMath::Abs(dangL) > nSigm * Interp(TMath::Abs(angTmean),0,i234-1)) continue;
	  
	  angL = angLmean;
	  angT = angTmean;
	}
		 
	MpdCellTrack *trNew = NULL;
	if (fGeo && h->GetUniqueID()) trNew = new ((*fKHits[h->Layer()-1])[nKH++]) MpdCellTrack(*hit); 	
	else if (fGeo) trNew = new ((*fKHits[h->Layer()-1])[nKH++]) MpdCellTrack(*hit);
	else trNew = new ((*fKHits[h->Layer()-1])[nKH++]) MpdCellTrack(*hit); 
	//printf("nomer treka=%d\n",i);
	trNew->SetUniqueID(0);
	trNew->SetNofDim(2);
	trNew->SetCosSin(0,angT);
	trNew->SetCosSin(1,angL);
	trNew->SetPrevTrack(itr);
	trNew->SetCode(track->GetCode());
	trNew->SetCode(i);
	//cout << i234 << " " << trNew->GetCode() << endl; 
	 
      } // for (Int_t itr = 0; itr < nTracks;
    } // for (Int_t i = 0; i < nHits2;

    printf("nKH=%d\n", nKH);
  } // for (Int_t i234 = 1; i234 < 4;)

}

//__________________________________________________________________________
void MpdCellAutomat::GetTrackSeeds(Int_t iPass) 
{
  /// Build ITS track seeds from Cell tracks
  if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(),__FUNCTION__);
  Int_t nCand = 0;
  TVector3 vert(0.,0.,0.),pmom;
  Int_t nTracks = fKHits[3]->GetEntriesFast();
 
  cout << "seed ITS tracks:"<< nTracks << " " << fCellMap.size() << endl;
  printf ("nTracks56789=%d\n",nTracks);
  
  MpdKalmanHit hit;
  hit.SetType(MpdKalmanHit::kFixedR);
  
  for (Int_t itr = 0; itr < nTracks; ++itr) {
    MpdCellTrack *track1 = (MpdCellTrack*) fKHits[3]->UncheckedAt(itr); // 4 layer
    // Check if this track has not been checked already during earlier passes
    if (iPass && fCellMap.find(track1->GetCode()) != fCellMap.end()) {
      cout << " Found: " << track1->GetCode() << " " << fCellMap[track1->GetCode()] << endl;
      continue;
    }
    fCellMap.insert(pair<TString,Int_t>(track1->GetCode(),1));
    //cout << itr << " " << track1->GetCode() << endl;

    MpdCellTrack *track2 = (MpdCellTrack*) fKHits[2]->UncheckedAt(track1->GetPrevTrack()); // index treka  3 layer
    MpdCellTrack *track3 = (MpdCellTrack*) fKHits[1]->UncheckedAt(track2->GetPrevTrack()); // index treka  2 layer 
    MpdCellTrack *track4 = (MpdCellTrack*) fKHits[0]->UncheckedAt(track3->GetPrevTrack()); // 1 layer
    ///tpc->GetParam()->Print();
    ///MpdEctKalmanTrack *track = new ((*fTracks)[nCand++]) MpdEctKalmanTrack(itr, *tpc);
    // MpdItsKalmanTrack *track = new ((*fTracks)[nCand++]) MpdItsKalmanTrack(*track1, vert); //sozdatia constructor iz CellTrack vvodit v ITS track constructor
    MpdItsKalmanTrack *track = new ((*fTrackCand)[nCand++]) MpdItsKalmanTrack(*track1, vert); //new 05.03.14
    // is Celltrack cozdavalsia constrycror v ITS Track!!!! ATTENTION 
    ///MpdITSKalmanTrack!!!
    Double_t pt = EvalPt(track1,track2);
    TVector3 v1 = track1->GetMeas();
    TVector3 v2 = track2->GetMeas();
    TVector3 v3 = track3->GetMeas();
    TVector3 v4 = track4->GetMeas();
    /*
    Double_t phiOut = v1.Phi(); // get azimuth angle
    Double_t phiIn = v2.Phi();
    Double_t rOut = v1.Pt();// 26.11 rastoynie ot tochki vzaimod-ya
    // Double_t phiOut = hitOut->GetMeas(0) / rOut;
    Double_t rIn = v2.Pt(); //26.11 rastoyanie
    */
    //< blok ot 25.12.2013
    Double_t phiOut = v3.Phi();
    Double_t phiIn = v4.Phi();
    Double_t rOut = v3.Pt();
    Double_t rIn = v4.Pt();
    //<end blok
    
    track->SetUniqueID(itr+1); ///
    ///<  parametrs track  in TPC:
    ///< 0: RPhi - coordinate in R-Phi direction
    ///< 1: Z - longitudinal coordinate
    ///< 2: Phi - local azimuthal angle
    ///< 3: Theta - local dip angle (angle w.r.t. the transverse plane)
    ///< 4: q/Pt - signed inverse Pt
    /*
    track->SetPos(rOut);
    track->SetParam (4, 1./pt); // q/Pt
    track->SetParam (3, track1->GetCosSin(1)); //longitudinal angL
    //track->SetParam (2, track1->GetCosSin(0)); //transverse angT
    track->SetParam (2, (v1-v2).Phi()); //transverse angT
    track->SetParam (1, v1.Z()); // Z - coordinate
    track->SetParam( 0,phiOut*rOut);
    */
    //< blok ot 25.12.2013
    track->SetPos(rIn);//<- ispravleno!
    track->SetParam (4, 1./pt); // q/Pt
    track->SetParam (3, track1->GetCosSin(1)); //longitudinal angL
    //track->SetParam (2, track1->GetCosSin(0)); //transverse angT
    track->SetParam (2, (v3-v4).Phi()); //Phi - rough estimate

    // Adjust Phi
    Double_t bz = FairRunAna::Instance()->GetField()->GetBz(0.,0.,0.);
    Double_t factor = 0.003 * bz / 10.; // 0.3 * 0.01 * 5kG / 10
    Double_t rCirc = TMath::Abs (pt / factor);
    Double_t ph = TMath::ASin ((v3-v4).Pt() / 2 / rCirc);
    track->SetParam (2, track->GetParam(2) - TMath::Sign(ph,pt));
    track->SetParam (1, v4.Z()); // Z - coordinate
    track->SetParam( 0,phiIn*rIn);
    //end blok
    //old code!
    /*
   track->SetParam(*track->GetParamAtHit());
   track->SetParamNew(*track->GetParamAtHit());
   track->SetPosNew(track->GetPos());
   track->SetWeight(*track->GetWeightAtHit());
   track->SetLength(track->GetLengAtHit());
   if (track.fHits == 0x0) return;
   //const MpdCellTrack *hitOut;// ukazatel nelzya nichego bratiat 02/12
    // const  MpdCellTrack *hitIn;// ukazatel
   // Double_t rOut = hitOut->GetPos();// GetPos() rasstoyanie ot tochki vzaimodeistviya? old
   */
    Double_t parOut[4] = {rOut,phiOut,0.,0.};
    Double_t parIn[4] = {rIn,phiIn,0.,0.};
    //EvalCovar(parOut,parIn,track,track1); //coment ot 25.12 work version
    EvalCovar(parOut,parIn,track,track4); // new ot 25.12 track4 -> 1 layer old
    track->SetPosNew(track->GetPos());
    track->SetParamNew(*track->GetParam());
    //track->ReSetWeight();
    //TMatrixDSym w = *track->GetWeight(); // save current weight matrix
    hit.SetPos(rIn-0.5);// radius tochki 4 layer 
    MpdKalmanFilter::Instance()->PropagateToHit(track,&hit,kFALSE);
    //track->SetWeight(w); // restore original weight matrix (near TPC inner shell)
    ///cout << nCand-1 << " " << track->GetTrackID() << endl;
    ///cout << track->GetHits()->GetEntriesFast() << " " << track->GetTrHits()->GetEntriesFast() << endl;
    track->GetHits()->Clear();
    track->SetChi2Its(track->GetChi2()); // temporary storage
    track->SetChi2(0.);
    track->SetDirection(MpdKalmanTrack::kOutward);
    
  }
  cout << " Number of ITS track candidates: " << nCand << endl;
  if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(),__FUNCTION__);
}  

//__________________________________________________________________________
Double_t MpdCellAutomat::EvalPt(const MpdCellTrack *track1, const MpdCellTrack *track2)
{
  /// Evaluate signed track Pt (curvature) assuming the track coming from the ///eval pt daet odin parametr a tam4
  /// primary vertex

  TVector3 v1 = track1->GetMeas();///attention!!!
  TVector3 v2 = track2->GetMeas();
  TVector2 vec1(v1.X(),v1.Y());
  TVector2 vec2(v2.X(),v2.Y());
  TVector2 vec21 = vec1 - vec2;
  Double_t cosAlpha = vec2 * vec21 / vec2.Mod() / vec21.Mod();
  Double_t rad = vec1.Mod() / 2. / TMath::Sin(TMath::ACos(cosAlpha));
  Double_t bz = FairRunAna::Instance()->GetField()->GetBz(0.,0.,0.);
  Double_t factor = 0.003 * bz / 10.; // 0.3 * 0.01 * 5kG / 10
  Double_t phi1 = vec1.Phi();
  Double_t phi2 = vec2.Phi();
  Double_t charge = phi1 - MpdKalmanFilter::Instance()->Proxim(phi1,phi2);
  if (track1->GetLayer() > track2->GetLayer()) charge = -charge;
  return factor * TMath::Abs(rad) * TMath::Sign(1., -charge);  
  }

//__________________________________________________________________________
void MpdCellAutomat::EvalCovar(Double_t *parOut, Double_t *parIn, MpdItsKalmanTrack *track,const MpdCellTrack *track1)
{ 
  /// Evaluate covariance matrix for track seed

  Double_t rIn = parIn[0], phiIn = parIn[1];
  Double_t rOut = parOut[0], phiOut = parOut[1];
  Double_t xIn = rIn * TMath::Cos(phiIn);
  Double_t yIn = rIn * TMath::Sin(phiIn);
  Double_t xOut = rOut * TMath::Cos(phiOut);
  Double_t yOut = rOut * TMath::Sin(phiOut);
  parOut[2] = xOut;
  parOut[3] = yOut;
  parIn[2] = xIn;
  parIn[3] = yIn;
  TVector3 v4 = track1->GetMeas();

  TMatrixD ww(5,5);
  ww(0,0) =  track1->GetErr(0) * track1->GetErr(0); // <RphiRphi> //26.11 error x
  ww(0,0) *= 4.; //extra factor of 4

  ww(1,1) = track1->GetErr(1) * track1->GetErr(1); // <zz> //error z

  Double_t dx = parOut[2] - parIn[2], dy = parOut[3] - parIn[3];
  Double_t dist2 = dx * dx + dy * dy;
  Double_t sinPhi = TMath::Sin (track->GetParam(2));
  Double_t cosPhi = TMath::Cos (track->GetParam(2));
  Double_t pOut = TMath::Cos(phiOut) * cosPhi + TMath::Sin(phiOut) * sinPhi;
  Double_t pIn = TMath::Cos(phiIn) * cosPhi + TMath::Sin(phiIn) * sinPhi;
  ww(2,2) = (pOut * pOut + pIn * pIn) / dist2 * ww(0,0); // <PhiPhi>
  ww(2,2) *= 2.; // extra factor of 2

  Double_t tanThe = TMath::Tan(track->GetParam(3));
  Double_t dRad = parOut[0] - parIn[0];
  Double_t denom = dRad * (1.+tanThe*tanThe);
  ww(3,3) = ww(1,1) * 2. / denom / denom; // <TheThe>
  ww(3,3) *= 1;

  ww(1,1) *= 8.; //AZ extra factor of 8

  //ww(4,4) = (track->GetParam(4)*0.5) * (track->GetParam(4)*0.5); // error 50%
  //(*fWeight)(4,4) = ((*fParam)(4,0)*0.75) * ((*fParam)(4,0)*0.75); // error 75%
  ww(4,4) = (track->GetParam(4)*1.) * (track->GetParam(4)*1.); // error 100%
  //fWeight->Print();
  //fWeight->Invert(); // weight matrix
  Int_t iok = 0;
  TMatrixD wwTmp = ww;
  MpdKalmanFilter::Instance()->MnvertLocal(ww.GetMatrixArray(), 5, 5, 5, iok);
  track->SetWeight(ww);
  //fWeight->Print();

  // Obtain errors
  //if (lunErr && idEct == tofP->GetTrackID()) fprintf(lunErr,"%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",track->GetParam(2),MpdKalmanFilter::Instance()->Proxim(track->GetParam(2),pmom.Phi()),track->GetParam(3),MpdKalmanFilter::Instance()->Proxim(track->GetParam(3),TMath::PiOver2()-pmom.Theta()),1./TMath::Sqrt((*track->GetWeight())(2,2)),1./TMath::Sqrt((*track->GetWeight())(3,3)),pmom.Mag(),rTof,pmom.Pt(),pt);
  TVector3 pos, posOut, pmom;
  // if (lunErr) fprintf(lunErr,"%-12.8s %-12.5s %-10.5s %-10.8s %-10.9s\n","|r*phi|","| Y |","| Z |","| Phi |","| Theta |");
  MpdStsHit *hit = (MpdStsHit*) fItsHits->UncheckedAt(track1->GetIndex());
  MpdStsPoint *p = (MpdStsPoint*)fItsPoints->UncheckedAt(hit->GetRefIndex());
  //hit->Position(pos);
  p->Position(pos);
  p->PositionOut(posOut);
  p->Momentum(pmom);
  Double_t phi = pmom.Phi();
  Double_t th = pmom.Theta();
  Double_t r = (pos.Pt() + posOut.Pt()) / 2;//radius
  Double_t phi1 = pos.Phi(); 
  Double_t phi2 = posOut.Phi(); 
  phi1 += MpdKalmanFilter::Instance()->Proxim(phi1, phi2);
  phi1 /= 2;
  Double_t zzz = (pos.Z() + posOut.Z()) / 2;
  /*
  if (lunErr) fprintf(lunErr,"%10.3f %10.3f %10.3f %10.3f %10.3f\n",r*phi1,TMath::PiOver2()-pmom.Theta(),
		      pos.Z(),phi,th);
  */
  /*
  if (lunErr) fprintf(lunErr,"%10.15s %10.15s %10.15s %10.15s \n","Param(0)-r*phi","Param(1)-z",
		      "Param(2)-Phi","Param(3)-Theta");
  	      
  if (lunErr&&p->GetTrackID()<1) fprintf(lunErr,"%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e\n",
		      r*(MpdKalmanFilter::Instance()->Proxim(phi1,track->GetParam(0)/r)-phi1),
		      track->GetParam(1)-zzz,MpdKalmanFilter::Instance()->Proxim(phi,track->GetParam(2))-phi,
		      TMath::PiOver2()-pmom.Theta()-track->GetParam(3),
		      TMath::Sqrt(wwTmp(0,0)),TMath::Sqrt(wwTmp(1,1)),TMath::Sqrt(wwTmp(2,2)),TMath::Sqrt(wwTmp(3,3)), phi); //9 parametr "phi" writing in "file.txt"
  */
  if (lunErr&&p->GetTrackID()<1) fprintf(lunErr,"%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e\n",
					 r*(MpdKalmanFilter::Instance()->Proxim(phi1,track->GetParam(0)/r)-phi1),
					 track->GetParam(1)-zzz,MpdKalmanFilter::Instance()->Proxim(phi,track->GetParam(2))-phi,
					 TMath::PiOver2()-pmom.Theta()-track->GetParam(3),
					 TMath::Sqrt(wwTmp(0,0)),TMath::Sqrt(wwTmp(1,1)),TMath::Sqrt(wwTmp(2,2)),TMath::Sqrt(wwTmp(3,3)));
  /*
    Int_t nHits = fItsHits->GetEntriesFast();
    for (Int_t i = 0; i < nHits; ++i){
    MpdStsHit *h = (MpdStsHit*) fItsHits->UncheckedAt(i);
    MpdStsPoint *h1 = (MpdStsPoint*) fItsPoints->UncheckedAt(h->GetRefIndex());
    ///primer coda
    }
    //if (lunErr) fclose(lunErr);
    //exit(0);
  */
  }
 
//__________________________________________________________________________
void MpdCellAutomat::DoTracking(Int_t iPass) 
{

  /// Run Kalman tracking
  if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(),__FUNCTION__);
  Double_t vert[3] = {0.0,0.0,0.0};
  //Int_t nCand = fTracks->GetEntriesFast(), iok = 0;// old
  Int_t nCand = fTrackCand->GetEntriesFast(), iok = 0;// new 05.03.14
  Int_t lay0 = ((MpdKalmanHit*)fKHits1->First())->GetLayer();//ms 06.05
  for (Int_t i = 0; i < nCand; ++i) {
    // MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTracks->UncheckedAt(i);// old
    MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTrackCand->UncheckedAt(i);// new 05.03.14
    //cout << " Track seed No. " << i << ", ID: " << track->GetTrackID() << ", Hits: " << track->GetNofTrHits() << endl;
    for (Int_t j = 0; j < track->GetNofTrHits(); ++j) {
      MpdKalmanHit *h = (MpdKalmanHit* )track->GetTrHits()->UncheckedAt(j);
      //MpdStsHit *hh = (MpdStsHit*) fItsHits->UncheckedAt(h->GetIndex());
      //Int_t id = ((FairMCPoint*) fItsPoints->UncheckedAt(hh->GetRefIndex()))->GetTrackID();
      //cout << j << " " << h->GetDist() << " " << h->GetLayer() << endl;
    }
    
    if (fGeo) iok = RunKalmanFilterCell(track); // from cell track
    //if (fGeo) iok = RunKalmanFilterCyl(track, lay0); // modular geometry
    // else iok = RunKalmanFilterCyl(track, lay0); // cylindrical geometry
    if (iok == -1) {
      // fTracks->RemoveAt(i); //old
      fTrackCand->RemoveAt(i); //new 05.03.14
      continue;
    }
    
    // Mark hits as being used 
    /*
      TObjArray *hits = track->GetHits();// track(MPdItsKalmantrack)
      for (Int_t j = 0; j < nHits; ++j) {
      hit = (MpdKalmanHit*) hits->UncheckedAt(j);
      hit->SetFlag(-1);
      }
    */
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
  //fTracks->Compress(); //old
  fTrackCand->Compress(); // new 05.03.14
  if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(),__FUNCTION__); //04.03
}
    
//__________________________________________________________________________
Int_t MpdCellAutomat::RunKalmanFilterCell(MpdItsKalmanTrack *track) 
{
  /// Run Kalman filter (fitter) for the hits from the cell track 
  /// (might not work when propagating outward!!!)

  Int_t layMax = ((MpdKalmanHit*)fKHits1->First())->GetLayer();
  
  MpdKalmanHit *hitOK = 0x0;
  MpdKalmanHit hitTmp;
  MpdKalmanTrack::TrackDir trackDir = track->GetDirection();
  //Int_t layBeg = 0, layEnd = -1, dLay = -1, layOK = -1; /// old coment
  Int_t layEnd = -1, dLay = -1, layOK = -1;
  if (trackDir == MpdKalmanTrack::kOutward) {
    layEnd = layMax + 1;
    dLay = 1;
  }
  
  TMatrixDSym pointWeight(5), pointWeightTmp(5), saveWeight(5);
  TMatrixD param(5,1), paramTmp(5,1);
  
  Double_t saveZ = 0.0, saveLeng = 0.0, dChi2Min = 0.0, posNew = 0.0;
  Int_t ok = 0;

  MpdCellTrack *cellTr = NULL;
  MpdCellTrack *track1 = (MpdCellTrack*) fKHits[3]->UncheckedAt(track->GetUniqueID()-1); //4
  MpdCellTrack *track2 = (MpdCellTrack*) fKHits[2]->UncheckedAt(track1->GetPrevTrack()); //3
  MpdCellTrack *track3 = (MpdCellTrack*) fKHits[1]->UncheckedAt(track2->GetPrevTrack()); //2
  MpdCellTrack *track4 = (MpdCellTrack*) fKHits[0]->UncheckedAt(track3->GetPrevTrack()); //1


for (Int_t lay = 0; lay < 4; ++lay) {
    // Get CellTrack from the outermost layer
  if (lay == 0) cellTr = (MpdCellTrack*) fKHits[lay]->UncheckedAt(track3->GetPrevTrack()); // 1 layer
  else if (lay == 1) cellTr = (MpdCellTrack*) fKHits[lay]->UncheckedAt(track2->GetPrevTrack()); // 2 layer
  if (lay == 2) cellTr = (MpdCellTrack*) fKHits[lay]->UncheckedAt(track1->GetPrevTrack()); // 3 layer
  else if (lay ==3) cellTr = (MpdCellTrack*) fKHits[lay]->UncheckedAt(track->GetUniqueID()-1); // 4 layer 
  /*
  for (Int_t lay = 3; lay > -1; --lay) {
    // Get CellTrack from the outermost layer
    if (lay == 3) cellTr = (MpdCellTrack*) fKHits[lay]->UncheckedAt(track->GetUniqueID()-1); // 4 layer
    else cellTr = (MpdCellTrack*) fKHits[lay]->UncheckedAt(cellTr->GetPrevTrack()); // previous layer
  */
    for (Int_t ihit = 0; ihit < 2; ++ihit) {
      // Loop over 2 sides of one layer
      
      MpdKalmanHit *hit = (MpdKalmanHit*) fKHits1->UncheckedAt(cellTr->GetIndex(ihit));
      // Exclude used hits
      // if (hit->GetFlag() < 0) continue;  //04.03
      // Propagate to hit (if it is not very close to the track)
      if (TMath::Abs(hit->GetPos()-track->GetPosNew()) > 1.e-4) { 
	Double_t leng = track->GetLength();
	Double_t posNew = track->GetPosNew();
	TMatrixD parNew = *track->GetParamNew();
	TString nodeNew = track->GetNodeNew();
	TString curPath = track->GetNode();
	if (!MpdKalmanFilter::Instance()->PropagateToHit(track,hit,kTRUE,kTRUE)) { 
	  // Restore initial parameters for the failed track ///vostanovlenie ischodnich parametrov for bad track
	  track->SetPosNew(posNew);
	  track->SetParamNew(parNew);
	  track->SetLength(leng);
	  track->SetNodeNew(nodeNew);
	  track->SetNode(curPath);
	  ok = -1; 
	  break; 
	} 
	  
	Double_t step = track->GetLength() - leng;
	Int_t lay1 = hit->GetLayer();
	
	//if (lay1 % 2 == 0 && step > 1.e-4) {
	if (lay1 % 2 != 0 && step > 1.e-4) {
	  // Crossing silicon layer - add mult. scat. in the sensor
	  Double_t x0 = 9.36; // rad. length
	  TMatrixDSym *cov = track->Weight2Cov();
	  Double_t th = track->GetParamNew(3);
	  Double_t cosTh = TMath::Cos(th);
	  Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, x0, step);
	  //cout << " Scat: " << hit->GetLayer() << " " << step << " " << TMath::Sqrt(angle2) << endl; ///old coment
	  (*cov)(2,2) += (angle2 / cosTh / cosTh);
	  (*cov)(3,3) += angle2;
	  Int_t iok = 0;
	  MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
	  track->SetWeight(*cov);
	} else if (0 && lay1 % 2 != 0 && step > 1.e-4 && fCables[lay1/2].size() > 0) {
	  // Crossing silicon layer - add mult. scat. in the cable
	  Double_t nCables = 0, x0 = 0.0116 * 2 / 9.36; // in rad. length - 116um cable per side
	  // Find number of cables crossed //poisk nomera cables skrechenich
	  TString path = gGeoManager->GetPath();
	  if (!path.Contains("sensor") && !path.Contains("sector")) {
	    cout << " !!! MpdCellAutomat::RunKalmanFilter - Outside detector !!! " << endl;
	    exit(0);
	  }
	  Double_t v7[3] = {track->GetParamNew(0), track->GetPosNew(), track->GetParamNew(1)}, v77[3];
	  gGeoManager->LocalToMaster(v7,v77);
	  Double_t zTr = TMath::Abs (v77[2]); // global Z
	  //cout << zTr << endl;
	  map<Double_t,Double_t>::iterator it;
	  for (it = fCables[lay1/2].begin(); it != fCables[lay1/2].end(); ++it) {
	    if (zTr < it->first || zTr > it->second) continue;
	    ++nCables;
	  }
	  //cout << " Cables: " << nCables << endl;
	  if (nCables) {
	    x0 *= nCables;
	    TMatrixDSym *cov = track->Weight2Cov();
	    Double_t th = track->GetParamNew(3);
	    Double_t cosTh = TMath::Cos(th);
	    Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, x0);
	    //cout << " Scat: " << hit->GetLayer() << " " << step << " " << TMath::Sqrt(angle2) << endl;
	    (*cov)(2,2) += (angle2 / cosTh / cosTh);
	    (*cov)(3,3) += angle2;
	    Int_t iok = 0;
	    MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
	    track->SetWeight(*cov);
	  }
	}
      }
 
      //cout << hit->GetDetectorID() << endl;
      // Exclude used hits!
      //if (hit->GetFlag() != 1) continue;// ? 04.03 coment work version
      if (hit->GetFlag() < 0) continue; //new 04.03 work 15.04
      // !!! Exact ID match
      if (fExact && TrackID(hit) != track->GetTrackID()) continue;
	
      Double_t dChi2 = MpdKalmanFilter::Instance()->FilterStripLocal(track,hit,pointWeight,param,posNew);
      // Add Z-contribution (if track is outside the detector)   
      Double_t sizeZ = MpdKalmanFilter::Instance()->GetGeo()->Size(hit).Y();
      //if (TMath::Abs(branchTr->GetParamNew(1)) > sizeZ) {
      if (TMath::Abs(param(1,0)) > sizeZ) {
	// Outside detector
	//Double_t dChi2z = (TMath::Abs(branchTr->GetParamNew(1)) - sizeZ) / hit->GetErr(1);
	Double_t dChi2z = (TMath::Abs(param(1,0)) - sizeZ) / hit->GetErr(1);
	dChi2 += dChi2z * dChi2z;
      }
      //if (hit->GetNofDim() == 1) cout << " lay, c2, id1: " << lay << " " << dChi2 << " " << TrackID(hit) << endl;
      //else cout << " lay, c2, id1, id2: " << lay << " " << dChi2 << " " << TrackID(hit) << " " << TrackID(hit,1) << endl;
      
      if (TMath::Abs(dChi2) < fgkChi2Cut) {
	//if (dChi2 < fgkChi2Cut && lay % 2 == 0) {
	track->GetHits()->Add(hit);
	track->SetChi2(track->GetChi2()+dChi2);
	TMatrixDSym w = *track->GetWeight();
	w += pointWeight;
	track->SetWeight(w);
	track->SetParamNew(param);
	// Save track params at last hit
	track->SetLengAtHit(track->GetLength()); 
	track->SetParamAtHit(param);
	track->SetWeightAtHit(*track->GetWeight());
      }
      Int_t lay1 = hit->GetLayer();
      if (lay1 % 2 != 0 && fCables[lay1/2].size() > 0) {
	// Crossing silicon layer - add mult. scat. in the cable
	Double_t nCables = 0, x0 = 0.0116 * 2 / 9.36; // in rad. length - 116um cable per side
	// Find number of cables crossed //poisk nomera cables skrechenich
	TString path = gGeoManager->GetPath();
	if (!path.Contains("sensor") && !path.Contains("sector")) {
	  cout << " !!! MpdCellAutomat::RunKalmanFilter - Outside detector !!! " << endl;
	  exit(0);
	}
	Double_t v7[3] = {track->GetParamNew(0), track->GetPosNew(), track->GetParamNew(1)}, v77[3];
	gGeoManager->LocalToMaster(v7,v77);
	Double_t zTr = TMath::Abs (v77[2]); // global Z
	//cout << zTr << endl;
	map<Double_t,Double_t>::iterator it;
	for (it = fCables[lay1/2].begin(); it != fCables[lay1/2].end(); ++it) {
	  if (zTr < it->first || zTr > it->second) continue;
	  ++nCables;
	}
	//cout << " Cables: " << nCables << endl;
	if (nCables) {
	  x0 *= nCables;
	  TMatrixDSym *cov = track->Weight2Cov();
	  Double_t th = track->GetParamNew(3);
	  Double_t cosTh = TMath::Cos(th);
	  Double_t angle2 = MpdKalmanFilter::Instance()->Scattering(track, x0);
	  //cout << " Scat: " << hit->GetLayer() << " " << step << " " << TMath::Sqrt(angle2) << endl;
	  (*cov)(2,2) += (angle2 / cosTh / cosTh);
	  (*cov)(3,3) += angle2;
	  Int_t iok = 0;
	  MpdKalmanFilter::Instance()->MnvertLocal(cov->GetMatrixArray(), 5, 5, 5, iok);
	  track->SetWeight(*cov);
	}
      }
    } // for (Int_t ihit = 0; ihit < 2;
  } // for (Int_t lay = 3; lay > -1;
  
  return 0;
}

//__________________________________________________________________________
void MpdCellAutomat::RemoveDoubles()
{
  /// Remove double tracks (keep the ones with better quality)

  Int_t ntracks = fTrackCand->GetEntriesFast(); //new 05.03.14
  cout << " Total tracks: " << ntracks << endl;
  MpdItsKalmanTrack *tr1, *tr2;
  MpdCellTrack *cellTr = NULL;
  TString code;
  
  for (Int_t i = 0; i < ntracks; i++) {
    tr1 = (MpdItsKalmanTrack*) fTrackCand->UncheckedAt(i);
    if (tr1 == 0x0) continue;
    for (Int_t j = i+1; j < ntracks; j++) {   			// j = 0 -> j = i+1
      //if (j == i) continue;					// add comment
      tr2 = (MpdItsKalmanTrack*) fTrackCand->UncheckedAt(j);
      if (tr2 == 0x0) continue;
      
      //Int_t nHitsCommon = GetNofCommonHits(tr1, tr2);
      //if ((float)nHitsCommon / TMath::Min(tr1->GetNofHits(),tr2->GetNofHits()) < 0.5) continue;
      if (!AreTracksDoubles(tr1, tr2)) continue;
      if (tr2->GetNofHits() < tr1->GetNofHits()) {
	fTrackCand->RemoveAt(j);
	cellTr = (MpdCellTrack*) fKHits[3]->UncheckedAt(tr2->GetUniqueID()-1);
	code = cellTr->GetCode();
	fCellMap[code] = -1; // used hit combination
	//cout << " Removed1: " << code << " " << tr2->GetNofHits() << endl;
      } else {
	if ((tr2->GetNofHits() > tr1->GetNofHits()) || (tr2->GetChi2() < tr1->GetChi2())) {
	  fTrackCand->RemoveAt(i);
	  cellTr = (MpdCellTrack*) fKHits[3]->UncheckedAt(tr1->GetUniqueID()-1);
	  code = cellTr->GetCode();
	  fCellMap[code] = -1; // used hit combination
	  //cout << " Removed2: " << code << " " << tr1->GetNofHits() << endl;
	  break;
	} else {
	  fTrackCand->RemoveAt(j);
	  cellTr = (MpdCellTrack*) fKHits[3]->UncheckedAt(tr2->GetUniqueID()-1);
	  code = cellTr->GetCode();
	  fCellMap[code] = -1; // used hit combination
	  //cout << " Removed3: " << code << " " << tr2->GetNofHits() << endl;
	}
      }
    } // for j
  } // for i

  fTrackCand->Compress();
  fNTracks = fTrackCand->GetEntriesFast();
}

//__________________________________________________________________________
Bool_t MpdCellAutomat::AreTracksDoubles(MpdItsKalmanTrack *tr1, MpdItsKalmanTrack *tr2)
{
  /// Searching common hits in 2 tracks to determine doubles
  // track1 contains fewer hits than track2

  MpdItsKalmanTrack *track1, *track2;
  if (tr1->GetNofHits() > tr2->GetNofHits()) //12.02
    track1 = tr2, track2 = tr1;
  else
    track1 = tr1, track2 = tr2;

  //Int_t limCommonPoint = (track1->GetNofHits()+1) / 2; // at least many common hits should be found //12.02
  Int_t limCommonPoint = 3; // new 25.03
  TObjArray *hits1 = track1->GetHits(), *hits2 = track2->GetHits();//12.02
  Int_t nh1 = hits1->GetEntriesFast(), nh2 = hits2->GetEntriesFast(), nHitsCommon = 0, j = 0;

  for (Int_t i = 0; i < nh1; i++){
    MpdKalmanHit *hit1 = (MpdKalmanHit*) hits1->UncheckedAt(i);
    for ( ; j < nh2; j++){
      MpdKalmanHit *hit2 = (MpdKalmanHit*) hits2->UncheckedAt(j);

      //is hit common for two tracks compared
      if (hit1 == hit2) {
        nHitsCommon++;
        break;
      }

      //if (hit2->GetLayer() < hit1->GetLayer()) break; // already closer to beam
      if (hit2->GetLayer() > hit1->GetLayer()) break; // already farther from beam
    }

    if (i+limCommonPoint-nHitsCommon > nh1) return kFALSE; // there'll be not enough common hits already
  }

  //if count of common hits is greater limit
  if (nHitsCommon < limCommonPoint) return kFALSE;

  // Test
  //MpdCellTrack *cell1 = (MpdCellTrack*) fKHits[3]->UncheckedAt(tr1->GetUniqueID()-1);
  //MpdCellTrack *cell2 = (MpdCellTrack*) fKHits[3]->UncheckedAt(tr2->GetUniqueID()-1);
  //cout << nHitsCommon << " " << cell1->GetCode() << " " << cell2->GetCode() << endl;

  return kTRUE;
}

//__________________________________________________________________________
Int_t MpdCellAutomat::TrackID(MpdKalmanHit *hit, Int_t indx) 
{
  /// Return track ID of the hit

  FairHit *h = (FairHit*) fItsHits->UncheckedAt(hit->GetIndex(indx));
  return ((MpdStsPoint*) fItsPoints->UncheckedAt(h->GetRefIndex()))->GetTrackID();
}

//__________________________________________________________________________
TVector2 MpdCellAutomat::GetDistance(MpdKalmanTrack *track, MpdKalmanHit *hit) 
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
void MpdCellAutomat::Write() 
{
  /// Write

  TFile histoFile("ItsRec.root","RECREATE");
  Writedir2current(fHistoDir);
  histoFile.Close();
}

//__________________________________________________________________________
void MpdCellAutomat::Writedir2current( TObject *obj ) 
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
void MpdCellAutomat::StoreTracks()
{
  if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(),__FUNCTION__);
  /// Transfer tracks from fTrackCand to fTracks
  
  Int_t nFound = fTracks->GetEntriesFast();
  for (Int_t i = 0; i < fNTracks; ++i) {
    MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTrackCand->UncheckedAt(i);
    //
    if (track->GetNofHits() < 7) continue; // !!! store only long tracks
    //
    track->Weight2Cov();
    new ((*fTracks)[nFound++]) MpdItsKalmanTrack(*track);
    MpdCellTrack *cellTr = (MpdCellTrack*) fKHits[3]->UncheckedAt(track->GetUniqueID()-1);
    TString code = cellTr->GetCode();
    fCellMap[code] = -1; // used hit combination
    //fTrackCand->RemoveAt(i);
  }
  fNTracks = 0;
  if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(),__FUNCTION__);
}

//__________________________________________________________________________
void MpdCellAutomat::AddHits()
{
  /// Add hit objects to tracks and compute number of wrongly assigned hits
  /// (hits with ID different from ID of starting TPC track)
  if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Start(Class()->GetName(),__FUNCTION__);
  Int_t nFound = fTracks->GetEntriesFast();//work 07
  //Int_t nFound = fTrackCand->GetEntriesFast();
  for (Int_t i = 0; i < nFound; ++i) {
    MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTracks->UncheckedAt(i);//work 07
    // MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTrackCand->UncheckedAt(i);
    cout << track->GetNode() << " " << track->GetNodeNew() << endl;
    Double_t c2 = track->GetChi2();
    track->SetChi2(track->GetChi2Its());
    track->SetChi2Its(c2);
    Int_t nHits = track->GetNofHits();
    if(nHits < 7) continue; // candidate not write track
    //if (nHits == 0) { fTracks->RemoveAt(i); continue; }
    TClonesArray &trHits = *track->GetTrHits();
    TObjArray *hits = track->GetHits();
    SetTrackID(track); // set track ID as ID of majority of its hits
    cout << nHits << " " << trHits.GetEntriesFast() << " " << track->GetTrackID() << endl;

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
      Int_t nOver = hit->Index()->GetSize();
      for (Int_t iov = 0; iov < nOver; ++iov) {
	MpdStsHit *h = (MpdStsHit*) fItsHits->UncheckedAt(hit->GetIndex(iov));
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
    }
    if (nHits) cout << "\n Wrongs: " << nWrong << endl;
    track->SetNofWrong(nWrong);
    //MpdKalmanTrack *tpc = (MpdKalmanTrack*) fTpcTracks->UncheckedAt(track->GetUniqueID()-1);
    //track->SetChi2(track->GetChi2()+tpc->GetChi2());
    //track->SetLastLay();
    //track->GetParam()->Print();
    track->SetNofHits(track->GetNofTrHits()); // TPC and ITS hits
    track->SetNofIts(nHits);
    cout << nHits << " " << track->GetNofTrHits() << " " << track->GetTrackID() << " " 
	 << track->GetChi2Its() << " " << track->GetChi2() << endl;
  }
  fTracks->Compress();//work 07
  //fTrackCand->Compress();//??
  if (MpdCodeTimer::Active()) MpdCodeTimer::Instance()->Stop(Class()->GetName(),__FUNCTION__);
}

//__________________________________________________________________________
void MpdCellAutomat::SetTrackID(MpdItsKalmanTrack* track)
{
  /// Set track ID as ID of majority of its hits

  const Int_t idMax = 9999;
  vector<Int_t> ids(idMax);
  Int_t nHits = track->GetNofHits(), locMax = 0, size = idMax, id0 = -1;
  TObjArray *hits = track->GetHits();

  for (Int_t i = 0; i < nHits; ++i) {
    MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(i);
    Int_t id = GetHitID(hit);
    if (i == 0) id0 = id; // first hit
    if (id >= size) { size = id + 1; ids.resize(size); }
    ++ids[id];
    if (ids[id] > ids[locMax]) locMax = id;
  }
  if (ids[locMax] > 1) track->SetTrackID(locMax);
  else track->SetTrackID(id0);
}

//__________________________________________________________________________
Int_t MpdCellAutomat::GetHitID(MpdKalmanHit *hit)
{
  /// Get hit ID from MCPoint ID

  Int_t nOver = hit->Index()->GetSize();
  for (Int_t iov = 0; iov < nOver; ++iov) {
    MpdStsHit *h = (MpdStsHit*) fItsHits->UncheckedAt(hit->GetIndex(iov));
    Int_t motherID1 = ((FairMCPoint*) fItsPoints->UncheckedAt(h->GetRefIndex()))->GetTrackID();
    return motherID1;
  }
}

//__________________________________________________________________________
void MpdCellAutomat::ExcludeHits()
{
  /// Exclude hits, already used for tracking, from consideration during the next passes

  Int_t nReco = fTracks->GetEntriesFast();// work 07
  //Int_t nReco = fTrackCand->GetEntriesFast();
  cout << " nReco: " << nReco << endl;
  for (Int_t i = 0; i < nReco; ++i) {
    MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTracks->UncheckedAt(i);// work 07
    //MpdItsKalmanTrack *track = (MpdItsKalmanTrack*) fTrackCand->UncheckedAt(i);
    Int_t nhitsKF = track->GetNofHits();
    TObjArray *hits = track->GetHits();
    for (Int_t j = 0; j < nhitsKF; ++j) {
      MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(j);
      MpdStsHit *h = (MpdStsHit*)fItsHits->UncheckedAt(hit->GetIndex());
      hit->SetFlag(-1);// work 15.04
      h->SetFlag(-1); // work 15.04
    }
  }
}

//__________________________________________________________________________
Double_t MpdCellAutomat::Interp(Double_t angt, Int_t choice, Int_t lay)
{
  // Perform parabolic interpolation for Angl (choice = 0) and dAngt (choice = 1)
 
  const Int_t np = 7;
  const Double_t xp[np] = {0.1/0.18, 0.1/0.069, 0.1/0.043, 0.1/0.034, 0.1/0.017, 0.1/0.0086, 0.1/0.0044};
  // const Double_t xp[np] =   {0.05,0.125,0.2, 0.25, 0.50, 1.00, 2.00}; //3.09

  Double_t dpp05[np] = {0.18, 0.069, 0.043, 0.034, 0.017, 0.0086, 0.0044}; //mean angt
  Double_t err05[np] = {0.0018, 0.00032, 0.00018, 0.00014, 0.000092, 0.00007, 0.00006}; // error mean angt

  

  const Double_t angs[2][np] = {{0.042, 0.019, 0.017, 0.016, 0.017, 0.019, 0.018}, //sigma angl
                                {0.054, 0.013, 0.0069, 0.006, 0.0036, 0.0034, 0.003}}; // sigma dangt Box(250)
  const Double_t dangs[2][np] = {{0.0015, 0.00067, 0.00069, 0.00065, 0.00073, 0.00065, 0.00068}, // error sigma angl
                                 {0.0027, 0.00045, 0.00027, 0.00018, 0.00014, 0.0001, 0.0001}}; //error sigma dangt
  /*
  const Double_t angls[3][np] = {{0.047, 0.027, 0.025, 0.023, 0.027, 0.026, 0.028},
				 {0.04, 0.022, 0.021, 0.02, 0.021, 0.021, 0.021},
				 {0.041, 0.014, 0.013, 0.012, 0.012, 0.012, 0.012}}; // sigma deltaz 2,3,4 layers
  const Double_t angts[2][np] = {{0.062, 0.013, 0.0075, 0.0058, 0.0041, 0.0039, 0.0036},
				 {0.059, 0.013, 0.0074, 0.006, 0.0035, 0.0029, 0.0026}}; //sigma dangt 3,4 layers
  */

  const Double_t angls[3][np] = {{0.045, 0.026, 0.025, 0.023, 0.027, 0.027, 0.028},
				 {0.02, 0.011, 0.01, 0.01, 0.011, 0.011, 0.01},
				 {0.016, 0.0045, 0.0037, 0.0036, 0.0033, 0.0033, 0.0031}}; // sigma deltaz 2,3,4 layers
  const Double_t angts[2][np] = {{0.031, 0.0068, 0.0037, 0.0029, 0.0021, 0.0019, 0.0018},
				 {0.019, 0.0039, 0.0021, 0.0016, 0.00091, 0.00059, 0.0005}}; //sigma dangt 3,4 layers

  Double_t x0 = 0.1 / angt;
  //Double_t x0 = angt;
  Int_t ok = 0, inds[3] = {np - 3, np - 2, np - 1};
  for (Int_t i = 0; i < np; ++i) {
    if (xp[i] < x0) continue;
    inds[1] = i;
    inds[0] = i - 1;
    inds[2] = i + 1;
    ok = 1;
    break;
  }

  if (ok) {
    if (inds[0] < 0) for (Int_t i = 0; i < 3; ++i) ++inds[i];
    if (inds[2] == np) for (Int_t i = 0; i < 3; ++i) --inds[i];
  }
  Double_t a,b,c,dy01,dy02;
  
  // Parabolic interpolation
  Double_t dx01 = xp[inds[1]] - xp[inds[0]], dx02 = xp[inds[2]] - xp[inds[0]];
  Double_t x01 = xp[inds[1]] + xp[inds[0]], x02 = xp[inds[2]] + xp[inds[0]];
  // Double_t dy01 = angs[choice][inds[1]] - angs[choice][inds[0]];
  // Double_t dy02 = angs[choice][inds[2]] - angs[choice][inds[0]];
  
  if (choice == 0) {
    dy01 = angls[lay][inds[1]] - angls[lay][inds[0]];
    dy02 = angls[lay][inds[2]] - angls[lay][inds[0]];
  }
  if (choice == 1) {
    dy01 = angts[lay][inds[1]] - angts[lay][inds[0]];
    dy02 = angts[lay][inds[2]] - angts[lay][inds[0]];
  }
  
  Double_t slope = dy01 / dx01;
  a = dy02 - slope * dx02;
  a /= (dx02 * x02 - dx02 * x01);
  
  b = slope - a * x01;
  if (choice == 0) {
    // c = angs[choice][inds[0]] - a * xp[inds[0]] * xp[inds[0]] - b * xp[inds[0]];
    c = angls[lay][inds[0]] - a * xp[inds[0]] * xp[inds[0]] - b * xp[inds[0]];
  }
  if (choice == 1) {
    // c = angs[choice][inds[0]] - a * xp[inds[0]] * xp[inds[0]] - b * xp[inds[0]];
    c = angts[lay][inds[0]] - a * xp[inds[0]] * xp[inds[0]] - b * xp[inds[0]];
  }
  //for (Int_t i = 0; i < 3; ++i) cout << xp[inds[i]] << " " << angs[choice][inds[i]] << " ";
  //cout << endl << x0 << " " << a * x0 * x0 + b * x0 + c << endl;  
  return a * x0 * x0 + b * x0 + c;
}

//__________________________________________________________________________
void MpdCellAutomat::GetShortTracks()
{
  /// Collect remaining tracks with small number of hits.

  map<TString,Int_t>::iterator it;
  const Int_t nLays = 4;
  Int_t inds[nLays];
  for (it = fCellMap.begin(); it != fCellMap.end(); ++it) {
    //cout << it->first << " " << it->second <<endl;
    if (it->second < 0) continue;
    TString code = it->first, cc;
    Int_t ibeg = 0, leng = -1;
    //cout << code << endl;
    for (Int_t i = 0; i < nLays; ++i) {
      leng = code.Index("-",leng);
      cc = code(ibeg,leng-ibeg);
      //cout << ibeg << " " << leng << endl;
      inds[i] = cc.Atoi();
      ibeg = leng + 1;
      ++leng;
      //cout << i << " " << cc << endl;
    }
    for (Int_t i = 0; i < nLays; ++i) {
      MpdCellTrack *tr = (MpdCellTrack*) f2DHits[i]->UncheckedAt(inds[i]);
      MpdStsHit *hit = (MpdStsHit*) fItsHits->UncheckedAt(tr->GetIndex());
      if (hit->GetFlag() < 0) { fCellMap[code] = -1; break; }
    }
  }

  Int_t left = 0;
  cout << " Leftovers: " << endl;
  for (it = fCellMap.begin(); it != fCellMap.end(); ++it) {
    if (it->second > 0) cout << ++left << " " << it->first << endl;
  }
}

ClassImp(MpdCellAutomat); 
