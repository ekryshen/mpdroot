/**
 * \class MpdMiniEvent
 * \brief Stores global MC information about the event
 *
 * The MpdMiniMcEvent class keeps variables that characterize Monte Carlo event.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \email nigmatkulov@gmail.com ; ganigmatkulov@mephi.ru
 * \date May 01, 2020
 */

#ifndef MpdMiniMcEvent_h
#define MpdMiniMcEvent_h

// ROOT headers
#include "TObject.h"
#include "TVector3.h"

//________________
class MpdMiniMcEvent : public TObject {

public:
   /// Default constructor
   MpdMiniMcEvent();
   /// Constructor with parameters
   MpdMiniMcEvent(UInt_t runId, UInt_t fEventId, Float_t fRP, Float_t fImpactPar, Short_t fNPart, Short_t fNcoll,
                  TVector3 primaryVertex, Float_t time);
   /// Copy constructor
   MpdMiniMcEvent(const MpdMiniMcEvent &event);
   /// Destructor
   virtual ~MpdMiniMcEvent();
   /// Print MC event information
   virtual void Print(const Char_t *option = "") const;

   //
   // Setters
   //

   /// Set unique run ID
   void setRunId(const UInt_t runId) { fRunId = runId; }
   /// Set unique event ID
   void setEventId(const UInt_t eventId) { fEventId = eventId; }
   /// Set reaction plane angle (rad)
   void setReactionPlaneAngle(const Float_t phi) { fReactionPlaneAngle = phi; }
   /// Set impact parameter
   void setImpactParameter(const Float_t b) { fImpactParameter = b; }
   /// Set number of participants (-1 if not set)
   void setNpart(const Short_t npart) { fNpart = npart; }
   /// Set number of binary collisions (-1 if not set)
   void setNcoll(const Short_t ncoll) { fNcoll = ncoll; }
   /// Set primary vertex position X (cm)
   void setPrimaryVertexX(const Float_t vtxX) { fPrimaryVertexX = vtxX; }
   /// Set primary vertex position Y (cm)
   void setPrimaryVertexY(const Float_t vtxY) { fPrimaryVertexY = vtxY; }
   /// Set primary vertex position Z (cm)
   void setPrimaryVertexZ(const Float_t vtxZ) { fPrimaryVertexZ = vtxZ; }
   /// Set primary vertex position as TVector3 (cm, cm, cm)
   void setPrimaryVertex(TVector3 vtx)
   {
      fPrimaryVertexX = (Float_t)vtx.X();
      fPrimaryVertexY = (Float_t)vtx.Y();
      fPrimaryVertexZ = vtx.Z();
   }
   /// Set primary vertex position (cm, cm, cm)
   void setPrimaryVertex(Double_t vtxX, Double_t vtxY, Double_t vtxZ)
   {
      fPrimaryVertexX = (Float_t)vtxX;
      fPrimaryVertexY = (Float_t)vtxY;
      fPrimaryVertexZ = (Float_t)vtxZ;
   }
   /// Set collision time (ns)
   void setTime(const Float_t time) { fTime = time; }

   //
   // Getters
   //

   /// Return unique run ID
   UInt_t runId() const { return fRunId; }
   /// Return unique event ID
   UInt_t eventId() const { return fEventId; }
   /// Return reaction plane angle (rad)
   Float_t reactionPlaneAngle() const { return fReactionPlaneAngle; }
   /// Return reaction plane angle (rad)
   Float_t phi() const { return reactionPlaneAngle(); }
   /// Return impact parameter (fm)
   Float_t impactParameter() const { return fImpactParameter; }
   /// Return impact parameter (fm)
   Float_t b() const { return impactParameter(); }
   /// Return primary vertex position (cm, cm, cm)
   TVector3 primaryVertex() const { return TVector3(fPrimaryVertexX, fPrimaryVertexY, fPrimaryVertexZ); }
   /// Return primary vertex x position (cm)
   Float_t primaryVertexX() const { return fPrimaryVertexX; }
   /// Return primary vertex y position (cm)
   Float_t primaryVertexY() const { return fPrimaryVertexY; }
   /// Return primary vertex z position (cm)
   Float_t primaryVertexZ() const { return fPrimaryVertexZ; }
   /// Return collision time (ns)
   Float_t time() const { return fTime; }

private:
   /// Unique run ID
   UInt_t fRunId;
   /// Unique event ID
   UInt_t fEventId;
   /// Reaction plane angle (rad)
   Float_t fReactionPlaneAngle;
   /// Impact parameter (fm)
   Float16_t fImpactParameter;
   /// Number of participants
   Short_t fNpart;
   /// Number of binary collisions
   Short_t fNcoll;

   /// Primary vertex x position (cm)
   Float_t fPrimaryVertexX;
   /// Primary vertex y position (cm)
   Float_t fPrimaryVertexY;
   /// Primary vertex z position (cm)
   Float_t fPrimaryVertexZ;
   /// Event time (ns)
   Float_t fTime;

   ClassDef(MpdMiniMcEvent, 1)
};

#endif // #define MpdMiniMcEvent_h
