//
// MpdMiniMcEvent class holds information about the Monte Carlo event
//

// MpdMiniDst information
#include "MpdMiniMessMgr.h"
#include "MpdMiniMcEvent.h"

ClassImp(MpdMiniMcEvent);

//_________________
MpdMiniMcEvent::MpdMiniMcEvent()
   : TObject(), fRunId(0), fEventId(0), fReactionPlaneAngle(0), fImpactParameter(-1), fNpart(-1), fNcoll(-1),
     fPrimaryVertexX(0), fPrimaryVertexY(0), fPrimaryVertexZ(0), fTime(0)
{
   // Default constructor
   /* empty */
}

//_________________
MpdMiniMcEvent::MpdMiniMcEvent(UInt_t runId, UInt_t eventId, Float_t rp, Float_t impactPar, Short_t npart,
                               Short_t ncoll, TVector3 vtx, Float_t time)
   : TObject(), fRunId(runId), fEventId(eventId), fReactionPlaneAngle(rp), fImpactParameter(impactPar), fNpart(npart),
     fNcoll(ncoll), fPrimaryVertexX(vtx.X()), fPrimaryVertexY(vtx.Y()), fPrimaryVertexZ(vtx.Z()), fTime(time)
{
   // Parametrized constructor
   /* empty */
}

//_________________
MpdMiniMcEvent::MpdMiniMcEvent(const MpdMiniMcEvent &copy) : TObject()
{
   // Copy constructor
   fRunId              = copy.fRunId;
   fEventId            = copy.fEventId;
   fReactionPlaneAngle = copy.fReactionPlaneAngle;
   fImpactParameter    = copy.fImpactParameter;
   fNpart              = copy.fNpart;
   fNcoll              = copy.fNcoll;
   fPrimaryVertexX     = copy.fPrimaryVertexX;
   fPrimaryVertexY     = copy.fPrimaryVertexY;
   fPrimaryVertexZ     = copy.fPrimaryVertexZ;
   fTime               = copy.fTime;
}

//_________________
MpdMiniMcEvent::~MpdMiniMcEvent()
{
   // Destructor
   /* empty */
}

//_________________
void MpdMiniMcEvent::Print(const Char_t *option __attribute__((unused))) const
{
   // Print information about MC event
   LOG_INFO << "run ID: " << fRunId << " event ID: " << fEventId << " phi (rad): " << fReactionPlaneAngle
            << " b (fm): " << fImpactParameter << " npart: " << fNpart << " ncoll: " << fNcoll << " vertex (x, y, z): ("
            << fPrimaryVertexX << "," << fPrimaryVertexY << "," << fPrimaryVertexZ << ")"
            << " time: " << fTime << endm;
}
