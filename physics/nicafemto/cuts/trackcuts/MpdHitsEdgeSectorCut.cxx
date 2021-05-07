/*
 * MpdHitsEdgeSectorCut.cxx
 *
 *  Created on: 14 sie 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "MpdHitsEdgeSectorCut.h"

#include <Rtypes.h>
#include <RtypesCore.h>
#include <TMath.h>
#include <TMathBase.h>
#include <TVector2.h>
#include <TVector3.h>
#include <stddef.h>

#include "NicaCut.h"
#include "NicaExpTrack.h"
#include "NicaHelix.h"
#include "NicaMpdTrack.h"
#include "NicaMpdTrackTpcPads.h"
#include "NicaPackage.h"
#include "NicaParameter.h"
#include "NicaTpcSectorGeo.h"
#include "NicaTrackTpcPads.h"

MpdHitsEdgeSectorCut::MpdHitsEdgeSectorCut() : NicaTrackCut(1) {
  SetUnitName("Hit Edge Cut [cm]", 0);
  fEdge = 1.5;
  fSec  = NULL;
}

MpdHitsEdgeSectorCut::MpdHitsEdgeSectorCut(const MpdHitsEdgeSectorCut& other) :
  NicaTrackCut(other), fEdge(other.fEdge), fSec(other.fSec) {}

Bool_t MpdHitsEdgeSectorCut::Init(Int_t format_id) {
  fSec = NicaTpcSectorGeo::Instance();
  return NicaTrackCut::Init(format_id);
}

Bool_t MpdHitsEdgeSectorCut::Pass(NicaTrack* track) {
  NicaMpdTrackTpcPads* tr = (NicaMpdTrackTpcPads*) track;
  TVector3 loc;
  Double_t dphi         = TMath::Pi() / 12.0;
  Double_t bad_points   = 0;
  Double_t total_points = 0;
  ULong64_t hit_map     = tr->GetHitMap();
  for (int iLay = 0; iLay < tr->GetTpcPadsInfo()->GetPadsNo(); iLay++) {
    if (!TESTBIT(hit_map, iLay)) continue;
    Double_t s       = tr->GetTpcPadsInfo()->GetPath(iLay);
    TVector3 pos     = tr->GetHelix().EvalPos(s);
    Int_t sector     = fSec->Sector(tr->GetTpcPadsInfo()->GetPadID(iLay));
    Double_t angle_m = fSec->SectorAngle(sector) - dphi;
    Double_t angle_p = angle_m + dphi * 2.0;
    Double_t phi_pos = pos.Phi();
    Double_t d1      = TVector2::Phi_mpi_pi(angle_m - phi_pos);
    Double_t d2      = TVector2::Phi_mpi_pi(angle_p - phi_pos);
    if (TMath::Abs(d1) < TMath::Abs(d1)) { angle_p = angle_m; }
    Double_t m = TMath::Tan(angle_p);
    Double_t d = TMath::Abs(m * pos.X() - pos.Y()) / TMath::Sqrt(1 + m * m);
    total_points++;
    if (d < fEdge) bad_points++;
  }
  //	if(bad_points!=0){
  //		std::cout<<bad_points<<" "<<total_points<<"
  //"<<std::bitset<32>(tr->GetBitMap())<<" "<<track->ClassName()<<"
  //"<<tr->GetNHits()<<std::endl;

  SetValue(bad_points / total_points);
  return Validate();
}

NicaPackage* MpdHitsEdgeSectorCut::Report() const {
  NicaPackage* report = NicaTrackCut::Report();
  report->AddObject(new NicaParameterString("Edge", fEdge));
  return report;
}

MpdHitsEdgeSectorCut::~MpdHitsEdgeSectorCut() {
  // TODO Auto-generated destructor stub
}
