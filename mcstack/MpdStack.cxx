// -------------------------------------------------------------------------
// -----                       FairStack source file                    -----
// -----             Created 10/08/04  by D. Bertini / V. Friese       -----
// -----              adopted for NICA/MPD 29/03/10  (litvin)          -----
// -----              adopted for NICA/MPD 20/12/19  (ABychkov)        -----
// -----           added external decayer 30/04/20  (AZinchenko)       -----
// -------------------------------------------------------------------------
#include "MpdStack.h"

#include "FairDetector.h"
#include "FairMCPoint.h"
#include "FairRootManager.h"
#include "MpdDecayer.h"
#include "MpdDecayerPyt8.h"
#include "MpdMCTrack.h"

#include "TError.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TRefArray.h"
#include "TVirtualMC.h"

#include <iostream>
#include <list>

using std::cout;
using std::endl;
using std::pair;
using std::vector;

// -----   Default constructor   -------------------------------------------
MpdStack::MpdStack(Int_t size)  // {
                                //  fParticles = new TClonesArray("TParticle", size);
                                //  fTracks    = new TClonesArray("MpdMCTrack", size);
                                //  fCurrentTrack = -1;
                                //  fNPrimaries = fNParticles = fNTracks = 0;
                                //  fIndex = 0;
                                //  fStoreSecondaries = kTRUE;
                                //  fMinPoints        = 1;
                                //  fEnergyCut        = 0.;
                                //  fStoreMothers     = kTRUE;
                                //}
  :
  FairGenericStack(),
  fStack(),
  fParticles(new TClonesArray("TParticle", size)),
  fTracks(new TClonesArray("MpdMCTrack", size)),
  fStoreMap(),
  fStoreIter(),
  fIndexMap(),
  fIndexIter(),
  fPointsMap(),
  fCurrentTrack(-1),
  fNPrimaries(0),
  fNParticles(0),
  fNTracks(0),
  fIndex(0),
  fDecayFlag(0),
  fStoreSecondaries(kTRUE),
  fMinPoints(1),
  // fEnergyCut(0.),
  fEnergyCut(-0.001),  // AZ
  fStoreMothers(kTRUE),
  fMinMotherMass(0.4),
  fMaxMotherMass(6.1),
  fRadiusCut(172),  // 175 NG ecal outer XYradius, some mm smaller
  fVzCut(311),      // 296 NG barrel half length, TODO does not include EndCaps
  fNoZDC(kTRUE)     // NG
{
  fDecays = new TClonesArray("TParticle");
}
// -------------------------------------------------------------------------


// -----   Destructor   ----------------------------------------------------
MpdStack::~MpdStack() {
  if (fParticles) {
    fParticles->Delete();
    delete fParticles;
  }
  if (fTracks) {
    fTracks->Delete();
    delete fTracks;
  }
  if (fDecays) {
    fDecays->Delete();
    delete fDecays;
  }
}
// -------------------------------------------------------------------------


// -----   Virtual public method PushTrack   -------------------------------
void MpdStack::PushTrack(Int_t toBeDone,
                         Int_t parentId,
                         Int_t pdgCode,
                         Double_t px,
                         Double_t py,
                         Double_t pz,
                         Double_t e,
                         Double_t vx,
                         Double_t vy,
                         Double_t vz,
                         Double_t time,
                         Double_t polx,
                         Double_t poly,
                         Double_t polz,
                         TMCProcess proc,
                         Int_t& ntr,
                         Double_t weight,
                         Int_t is) {

  // --> Get TParticle array
  TClonesArray& partArray = *fParticles;

  // --> Create new TParticle and add it to the TParticle array
  Int_t trackId       = fNParticles;
  Int_t nPoints       = 0;
  Int_t daughter1Id   = -1;
  Int_t daughter2Id   = -1;
  TParticle* particle = new (partArray[fNParticles++])
    TParticle(pdgCode, trackId, parentId, nPoints, daughter1Id, daughter2Id, px, py, pz, e, vx, vy, vz, time);
  particle->SetPolarisation(polx, poly, polz);
  particle->SetWeight(weight);
  particle->SetUniqueID(proc);
  // AZ - global polarization transfer to decay Lambda
  if (TMath::Abs(pdgCode) == 3122 && proc == kPDecay) {
    TParticle* moth = (TParticle*) partArray.UncheckedAt(parentId);
    TVector3 polar;
    moth->GetPolarisation(polar);
    Double_t value = moth->GetWeight();
    if (TMath::Abs(moth->GetPdgCode()) == 3312)
      value *= 0.927;  // Xi+-
    else if (TMath::Abs(moth->GetPdgCode()) == 3322)
      value *= 0.900;  // Xi0
    else if (TMath::Abs(moth->GetPdgCode()) == 3212) {
      // Sigma0
      value *= (1. / 3);
      polar *= -1.0;
    }
    particle->SetPolarisation(polar);
    particle->SetWeight(value);
  }

  // --> Increment counter
  if (parentId < 0) fNPrimaries++;

  // --> Set argument variable
  ntr = trackId;

  // AZ
  if (fDecayFlag && (!particle->GetPDG()->Stable() || (particle->GetPdgCode() == 221) || (particle->GetPdgCode() == 111))
      && particle->GetPDG()->DecayList()) {
    // Decay unstable particle
    // MpdDecayer::Instance()->Decay(particle);
    // Int_t npart = MpdDecayer::Instance()->ImportParticles(fDecays);
    MpdDecayerPyt8::Instance()->Decay(particle);
    Int_t npart = MpdDecayerPyt8::Instance()->ImportParticles(fDecays);
    // Copy decay products to vector to avoid TClonesArray overwriting
    // for recursive decays
    vector<TParticle> vDecays;

    for (Int_t j = 0; j < npart; ++j) {
      TParticle tpart = *((TParticle*) fDecays->UncheckedAt(j));
      tpart.Print();
      vDecays.push_back(tpart);
    }

    for (Int_t j = 0; j < npart; ++j) {  // skip mother particle
      if (j == 0) {
        particle->SetStatusCode(11);  // decayed particle
        continue;                     // skip mother particle
      }
      // TParticle *part = (TParticle*) fDecays->UncheckedAt(j);
      TParticle* part = &vDecays[j];
      Int_t toDo      = 1;
      if (TString(gMC->GetName()) == "TGeant4") toDo = 0;
      // PushTrack(toBeDone, trackId, part->GetPdgCode(),
      PushTrack(toDo,
                trackId,
                part->GetPdgCode(),
                part->Px(),
                part->Py(),
                part->Pz(),
                part->Energy(),
                part->Vx(),
                part->Vy(),
                part->Vz(),
                time,
                polx,
                poly,
                polz,
                proc,
                ntr,
                weight,
                is);
    }
    if (parentId < 0) fNPrimaries += (npart - 1);  // treat decay products as primaries (dirty trick)
  } else {
    // --> Push particle on the stack if toBeDone is set
    if (toBeDone == 1) fStack.push(particle);
  }

  // --> Push particle on the stack if toBeDone is set
  // AZ if (toBeDone == 1) fStack.push(particle);
}

void MpdStack::PushTrack(Int_t toBeDone,
                         Int_t parentId,
                         Int_t pdgCode,
                         Double_t px,
                         Double_t py,
                         Double_t pz,
                         Double_t e,
                         Double_t vx,
                         Double_t vy,
                         Double_t vz,
                         Double_t time,
                         Double_t polx,
                         Double_t poly,
                         Double_t polz,
                         TMCProcess proc,
                         Int_t& ntr,
                         Double_t weight,
                         Int_t is,
                         Int_t secondMotherID) {

  // --> Get TParticle array
  TClonesArray& partArray = *fParticles;

  // --> Create new TParticle and add it to the TParticle array
  Int_t trackId     = fNParticles;
  Int_t nPoints     = 0;
  Int_t daughter1Id = -1;
  Int_t daughter2Id = -1;
  TParticle* particle =
    //    new(partArray[fNParticles++]) TParticle(pdgCode, trackId, parentId,
    new (partArray[fNParticles++])
      TParticle(pdgCode, trackId, secondMotherID, nPoints, daughter1Id, daughter2Id, px, py, pz, e, vx, vy, vz, time);
  particle->SetPolarisation(polx, poly, polz);
  particle->SetWeight(weight);
  particle->SetUniqueID(proc);
  // AZ - global polarization transfer to decay Lambda
  if (TMath::Abs(pdgCode) == 3122 && proc == kPDecay) {
    TParticle* moth = (TParticle*) partArray.UncheckedAt(parentId);
    TVector3 polar;
    moth->GetPolarisation(polar);
    Double_t value = moth->GetWeight();
    if (TMath::Abs(moth->GetPdgCode()) == 3312)
      value *= 0.927;  // Xi+-
    else if (TMath::Abs(moth->GetPdgCode()) == 3322)
      value *= 0.900;  // Xi0
    else if (TMath::Abs(moth->GetPdgCode()) == 3212) {
      // Sigma0
      value *= (1. / 3);
      polar *= -1.0;
    }
    particle->SetPolarisation(polar);
    particle->SetWeight(value);
  }

  // --> Increment counter
  if (parentId < 0) fNPrimaries++;

  // --> Set argument variable
  ntr = trackId;

  // AZ
  if (fDecayFlag && (!particle->GetPDG()->Stable() || (particle->GetPdgCode() == 221) || (particle->GetPdgCode() == 111))
      && particle->GetPDG()->DecayList()) {
    // Decay unstable particle
    // MpdDecayer::Instance()->Decay(particle);
    // Int_t npart = MpdDecayer::Instance()->ImportParticles(fDecays);
    MpdDecayerPyt8::Instance()->Decay(particle);
    Int_t npart = MpdDecayerPyt8::Instance()->ImportParticles(fDecays);
    // Copy decay products to vector to avoid TClonesArray overwriting
    // for recursive decays
    vector<TParticle> vDecays;

    for (Int_t j = 0; j < npart; ++j) {
      TParticle tpart = *((TParticle*) fDecays->UncheckedAt(j));
      tpart.Print();
      vDecays.push_back(tpart);
    }

    for (Int_t j = 0; j < npart; ++j) {
      if (j == 0) {
        particle->SetStatusCode(11);  // decayed particle
        continue;                     // skip mother particle
      }
      // TParticle *part = (TParticle*) fDecays->UncheckedAt(j);
      TParticle* part = &vDecays[j];
      Int_t toDo      = 1;
      if (TString(gMC->GetName()) == "TGeant4") toDo = 0;
      // PushTrack(toBeDone, trackId, part->GetPdgCode(),
      PushTrack(toDo,
                trackId,
                part->GetPdgCode(),
                part->Px(),
                part->Py(),
                part->Pz(),
                part->Energy(),
                part->Vx(),
                part->Vy(),
                part->Vz(),
                time,
                polx,
                poly,
                polz,
                proc,
                ntr,
                weight,
                is);
    }
    if (parentId < 0) fNPrimaries += (npart - 1);  // treat decay products as primaries (dirty trick)
  } else {
    // --> Push particle on the stack if toBeDone is set
    if (toBeDone == 1) fStack.push(particle);
  }

  // --> Push particle on the stack if toBeDone is set
  // AZ if (toBeDone == 1) fStack.push(particle);
}
// -------------------------------------------------------------------------


// -----   Virtual method PopNextTrack   -----------------------------------
TParticle* MpdStack::PopNextTrack(Int_t& iTrack) {

  // If end of stack: Return empty pointer
  if (fStack.empty()) {
    iTrack = -1;
    return NULL;
  }

  // If not, get next particle from stack
  TParticle* thisParticle = fStack.top();
  fStack.pop();

  if (!thisParticle) {
    iTrack = 0;
    return NULL;
  }

  fCurrentTrack = thisParticle->GetStatusCode();
  iTrack        = fCurrentTrack;

  return thisParticle;
}
// -------------------------------------------------------------------------


// -----   Virtual method PopPrimaryForTracking   --------------------------
TParticle* MpdStack::PopPrimaryForTracking(Int_t iPrim) {

  // Get the iPrimth particle from the fStack TClonesArray. This
  // should be a primary (if the index is correct).

  // Test for index
  // AZ-15.05.20 if (iPrim < 0 || iPrim >= fNPrimaries) {
  if (iPrim < 0 || iPrim >= fParticles->GetEntriesFast()) {
    cout << fParticles->GetEntriesFast() << endl;
    cout << "-E- MpdStack: Primary index out of range! " << iPrim << endl;
    Fatal("MpdStack::PopPrimaryForTracking", "Index out of range");
  }

  // Return the iPrim-th TParticle from the fParticle array. This should be
  // a primary.
  TParticle* part = (TParticle*) fParticles->At(iPrim);
  /*AZ-30.04.20
  if ( ! (part->GetMother(0) < 0) ) {
    cout << "-E- MpdStack:: Not a primary track! " << iPrim << endl;
    Fatal("MpdStack::PopPrimaryForTracking", "Not a primary track");
  }
  */
  if (part->GetStatusCode() == 11) return (TParticle*) NULL;  // AZ-30.04.20

  return part;
}
// -------------------------------------------------------------------------


// -----   Virtual public method GetCurrentTrack   -------------------------
TParticle* MpdStack::GetCurrentTrack() const {
  TParticle* currentPart = GetParticle(fCurrentTrack);
  if (!currentPart) {
    cout << "-W- MpdStack: Current track not found in stack!" << endl;
    Warning("MpdStack::GetCurrentTrack", "Track not found in stack");
  }
  return currentPart;
}
// -------------------------------------------------------------------------


// -----   Public method AddParticle   -------------------------------------
void MpdStack::AddParticle(TParticle* oldPart) {
  TClonesArray& array = *fParticles;
  TParticle* newPart  = new (array[fIndex]) TParticle(*oldPart);
  newPart->SetWeight(oldPart->GetWeight());
  newPart->SetUniqueID(oldPart->GetUniqueID());
  fIndex++;
}
// -------------------------------------------------------------------------


// -----   Public method FillTrackArray   ----------------------------------
void MpdStack::FillTrackArray() {

  cout << "-I- MpdStack: Filling MCTrack array..." << endl;

  // --> Reset index map and number of output tracks
  fIndexMap.clear();
  fNTracks = 0;

  // --> Check tracks for selection criteria
  SelectTracks();

  // --> Loop over fParticles array and copy selected tracks
  for (Int_t iPart = 0; iPart < fNParticles; iPart++) {

    fStoreIter = fStoreMap.find(iPart);
    if (fStoreIter == fStoreMap.end()) {
      cout << "-E- MpdStack: Particle " << iPart << " not found in storage map!" << endl;
      Fatal("MpdStack::FillTrackArray", "Particle not found in storage map.");
    }
    Bool_t store = (*fStoreIter).second;

    if (store) {
      MpdMCTrack* track = new ((*fTracks)[fNTracks]) MpdMCTrack(GetParticle(iPart));
      fIndexMap[iPart]  = fNTracks;
      // --> Set the number of points in the detectors for this track
      for (Int_t iDet = kSTS; iDet <= kFLOOR; iDet++) {
        pair<Int_t, Int_t> a(iPart, iDet);
        track->SetNPoints(iDet, fPointsMap[a]);
      }
      fNTracks++;
    } else
      fIndexMap[iPart] = -2;
  }

  // --> Map index for primary mothers
  fIndexMap[-1] = -1;

  // --> Screen output
  Print(0);
}
// -------------------------------------------------------------------------


// -----   Public method UpdateTrackIndex   --------------------------------
void MpdStack::UpdateTrackIndex(TRefArray* detList) {

  cout << "-I- MpdStack: Updating track indizes...";
  Int_t nColl = 0;

  // First update mother ID in MCTracks
  for (Int_t i = 0; i < fNTracks; i++) {
    MpdMCTrack* track = (MpdMCTrack*) fTracks->At(i);
    Int_t iMotherOld  = track->GetMotherId();
    if (iMotherOld < 0) continue;  // AZ
    fIndexIter = fIndexMap.find(iMotherOld);
    if (fIndexIter == fIndexMap.end()) {
      cout << "-E- MpdStack: Particle index " << iMotherOld << " not found in dex map! " << endl;
      Fatal("MpdStack::UpdateTrackIndex", "Particle index not found in map");
    }
    track->SetMotherId((*fIndexIter).second);
  }

  // Now iterate through all active detectors
  TIterator* detIter = detList->MakeIterator();
  detIter->Reset();
  FairDetector* det = NULL;
  while ((det = (FairDetector*) detIter->Next())) {


    // --> Get hit collections from detector
    Int_t iColl = 0;
    TClonesArray* hitArray;
    while ((hitArray = det->GetCollection(iColl++))) {
      nColl++;
      Int_t nPoints = hitArray->GetEntriesFast();

      // --> Update track index for all MCPoints in the collection
      for (Int_t iPoint = 0; iPoint < nPoints; iPoint++) {
        FairMCPoint* point = (FairMCPoint*) hitArray->At(iPoint);
        Int_t iTrack       = point->GetTrackID();

        fIndexIter = fIndexMap.find(iTrack);
        if (fIndexIter == fIndexMap.end()) {
          cout << "-E- MpdStack: Particle index " << iTrack << " not found in index map! " << endl;
          Fatal("MpdStack::UpdateTrackIndex", "Particle index not found in map");
        }
        point->SetTrackID((*fIndexIter).second);
        // point->SetLink(FairLink("MCTrack", (*fIndexIter).second));
      }

    }  // Collections of this detector
  }    // List of active detectors

  cout << "...stack and " << nColl << " collections updated." << endl;
}
// -------------------------------------------------------------------------


// -----   Public method Reset   -------------------------------------------
void MpdStack::Reset() {
  fIndex        = 0;
  fCurrentTrack = -1;
  fNPrimaries = fNParticles = fNTracks = 0;
  while (!fStack.empty())
    fStack.pop();
  fParticles->Clear();
  fTracks->Clear();
  fPointsMap.clear();
}
// -------------------------------------------------------------------------


// -----   Public method Register   ----------------------------------------
void MpdStack::Register() { FairRootManager::Instance()->Register("MCTrack", "Stack", fTracks, kTRUE); }
// -------------------------------------------------------------------------


// -----   Public method Print  --------------------------------------------
void MpdStack::Print(Int_t iVerbose) const {
  cout << "-I- MpdStack: Number of primaries        = " << fNPrimaries << endl;
  cout << "              Total number of particles  = " << fNParticles << endl;
  cout << "              Number of tracks in output = " << fNTracks << endl;
  if (iVerbose) {
    for (Int_t iTrack = 0; iTrack < fNTracks; iTrack++)
      ((MpdMCTrack*) fTracks->At(iTrack))->Print(iTrack);
  }
}
// -------------------------------------------------------------------------


// -----   Public method AddPoint (for current track)   --------------------
void MpdStack::AddPoint(DetectorIdMPD detId) {
  Int_t iDet = detId;
  pair<Int_t, Int_t> a(fCurrentTrack, iDet);
  if (fPointsMap.find(a) == fPointsMap.end())
    fPointsMap[a] = 1;
  else
    fPointsMap[a]++;
}
// -------------------------------------------------------------------------


// -----   Public method AddPoint (for arbitrary track)  -------------------
void MpdStack::AddPoint(DetectorIdMPD detId, Int_t iTrack) {
  if (iTrack < 0) return;
  Int_t iDet = detId;
  pair<Int_t, Int_t> a(iTrack, iDet);
  if (fPointsMap.find(a) == fPointsMap.end())
    fPointsMap[a] = 1;
  else
    fPointsMap[a]++;
}
// -------------------------------------------------------------------------


// -----   Virtual method GetCurrentParentTrackNumber   --------------------
Int_t MpdStack::GetCurrentParentTrackNumber() const {
  TParticle* currentPart = GetCurrentTrack();
  if (currentPart)
    return currentPart->GetFirstMother();
  else
    return -1;
}
// -------------------------------------------------------------------------


// -----   Public method GetParticle   -------------------------------------
TParticle* MpdStack::GetParticle(Int_t trackID) const {
  if (trackID < 0 || trackID >= fNParticles) {
    cout << "-E- MpdStack: Particle index " << trackID << " out of range." << endl;
    Fatal("MpdStack::GetParticle", "Index out of range");
  }
  return (TParticle*) fParticles->At(trackID);
}
// -------------------------------------------------------------------------


// -----   Private method SelectTracks   -----------------------------------
void MpdStack::SelectTracks() {


  //#define EDEBUG
#ifdef EDEBUG
  static Int_t lEDEBUGcounter = 0;
  if (lEDEBUGcounter <= 20) std::cout << "EDEBUG-- MpdStack::SelectTracks(): entered " << fNParticles << endl;
  lEDEBUGcounter++;
#endif

  // --> Clear storage map
  fStoreMap.clear();

  // --> Check particles in the fParticle array
  for (Int_t i = 0; i < fNParticles; i++) {

    TParticle* thisPart = GetParticle(i);
    Bool_t store        = kTRUE;

#ifdef EDEBUG
    if (lEDEBUGcounter <= 20) thisPart->Print();
    lEDEBUGcounter++;
#endif
#undef EDEBUG

    // --> Get track parameters
    Int_t iMother = thisPart->GetMother(0);
    TLorentzVector p;
    thisPart->Momentum(p);
    Double_t energy = p.E();
    Double_t mass   = thisPart->GetMass();
    Double_t eKin   = energy - mass;

    // --> Calculate number of points
    Int_t nPoints = 0, nPointsZDC = 0, nPointsNoZDC = 0;
    for (Int_t iDet = kSTS; iDet <= kFLOOR; iDet++) {
      pair<Int_t, Int_t> a(i, iDet);
      if (fPointsMap.find(a) != fPointsMap.end()) {
        nPoints += fPointsMap[a];
        if (iDet == kZDC)
          nPointsZDC += fPointsMap[a];  // NG save points in ZDC
        else
          nPointsNoZDC += fPointsMap[a];  // points not in ZDC
      }
    }

    // --> Check for cuts (store primaries in any case)
    if (iMother < 0)
      store = kTRUE;
    else {
      if (!fStoreSecondaries) store = kFALSE;
      if (nPoints < fMinPoints) store = kFALSE;
      if (eKin < fEnergyCut) store = kFALSE;
      if (fNoZDC && nPointsZDC && !nPointsNoZDC) store = kFALSE;  // NG cut many secondary zdc tracks
      // !!!!!AZ - store products of potentially interesting decays
      if (nPoints == 0) {
        store = kFALSE;
        if (fStoreMothers) {
          TParticle* mother = GetParticle(iMother);
          if (mass < 0.1 && thisPart->P() < 0.01)
            store = kFALSE;  // P < cut on e+- and gamma
          else if (mother->GetPDG()->Mass() > fMinMotherMass && mother->GetPDG()->Mass() < fMaxMotherMass &&  // 0.4 and 6.1
                   TMath::Abs(mother->GetPDG()->Mass() - 0.940) > 0.004 && mother->R() < fRadiusCut && abs(mother->Vz()) < fVzCut)
            store = kTRUE;  // except nucleons and heavy fragments
        }
      }  // no points
    }
    // --> Set storage flag
    fStoreMap[i] = store;
  }

  // AZ If flag is set, flag recursively mothers of selected tracks and their sisters
  if (fStoreMothers) {
    // Fill mother map
    std::multimap<Int_t, Int_t> moths;
    std::multimap<Int_t, Int_t>::iterator it;
    pair<std::multimap<Int_t, Int_t>::iterator, std::multimap<Int_t, Int_t>::iterator> ret;
    std::map<Int_t, Bool_t> copyMap = fStoreMap;

    for (Int_t i = 0; i < fNParticles; i++) {
      Int_t iMother = GetParticle(i)->GetMother(0);
      if (iMother >= 0) moths.insert(pair<Int_t, Int_t>(iMother, i));
    }

    for (Int_t i = 0; i < fNParticles; i++) {
      if (copyMap[i]) {  // if particle is to be stored
        Int_t iMother = GetParticle(i)->GetMother(0);

        while (iMother >= 0) {  // and it's mother is not primary
          TParticle* mother = GetParticle(iMother);
          if (mother->GetPDG()->Mass() < 0.1 && mother->P() < 0.01) break;
          // if (mother->GetPDG()->Mass() < 0.1 && mother->P() < 0.001) break;
          if (abs(mother->Vz()) > fVzCut || mother->R() > fRadiusCut) break;  // NG and doesn't originate within Rcut/Zcut
          fStoreMap[iMother] = kTRUE;                                         // store the mother
          ret                = moths.equal_range(iMother);

          for (it = ret.first; it != ret.second; ++it) {
            TParticle* daught = GetParticle(it->second);
            if (abs(daught->Vz()) > fVzCut || daught->R() > fRadiusCut) continue;  // NG and doesn't originate within Rcut/Zcut
            if (daught->GetPDG()->Mass() < 0.1 && daught->P() < 0.01) continue;
            // if (daught->GetPDG()->Mass() < 0.1 && daught->P() < 0.001) continue;
            fStoreMap[it->second] = kTRUE;  // sister
          }
          iMother = mother->GetMother(0);  // mother of the mother, loop again
        }                                  // while (iMother >= 0)

      }  // if (copyMap[i])
    }    // for (Int_t i=0; i<fNParticles;

    // AZ - If mother is not stored, find the nearest stored ancestor
    for (Int_t i = 0; i < fNParticles; ++i) {
      if (fStoreMap[i]) {
        Int_t iMother = GetParticle(i)->GetMother(0);

        while (iMother >= 0 && fStoreMap[iMother] == kFALSE) {
          TParticle* mother = GetParticle(iMother);
          iMother           = mother->GetMother(0);  // mother of the mother, loop again
        }
        if (iMother >= 0) GetParticle(i)->SetFirstMother(iMother);
      }
    }
  }  // if (fStoreMothers)

  //*/
}
// -------------------------------------------------------------------------


ClassImp(MpdStack)
