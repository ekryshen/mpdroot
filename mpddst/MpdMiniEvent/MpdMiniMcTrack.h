/**
 * \class MpdMiniMcTrack
 * \brief Stores information about the Monte Carlo track
 *
 * The MpdMiniMcTrack class holds Monte Carlo track parameters.
 * The tracks then passed through the full reconstruction chain.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \email nigmatkulov@gmail.com ; ganigmatkulov@mephi.ru
 * \date May 01, 2020
 */

#ifndef MpdMiniMcTrack_h
#define MpdMiniMcTrack_h

// C++ headers
#include <limits>
#include <vector>

// ROOT headers
#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

//_________________
class MpdMiniMcTrack : public TObject {

 public:
  /// Default constructor
  MpdMiniMcTrack();
  /// Copy constructor
  MpdMiniMcTrack(const MpdMiniMcTrack& track);
  /// Destructor
  virtual ~MpdMiniMcTrack();
  /// Print MC track information
  virtual void Print(const Char_t* option="") const;

  //
  // Setters
  //
  /// Set particle index
  void setId(Int_t index)
  { fId = ( (index > std::numeric_limits<unsigned short>::max() ) ?
	    std::numeric_limits<unsigned short>::max() : (UShort_t)index ); }
  /// Set PdgId (pdg code)
  void setPdgId(Int_t pdg)                                 { fPdgId = pdg; }
  /// Set status
  void setStatus(Int_t status)
  { if ( status <= std::numeric_limits<char>::min() )      { fStatus = std::numeric_limits<char>::min(); }
    else if ( status >= std::numeric_limits<char>::max() ) { fStatus = std::numeric_limits<char>::max(); }
    else { fStatus = (Char_t)status; } }
  // Set parent index
  /*
  void setParentIndex(Int_t parent)
  { fParentIndex = ( ( parent > std::numeric_limits<unsigned short>::max() ) ?
		     std::numeric_limits<unsigned short>::max() : (UShort_t)parent ); }
  */
  /// Set two daughter track indices (the first and the last one)
  void setChild(Int_t child[2])                             { setFirstChild(child[0]); setLastChild(child[1]); }
  /// Set index of the first child
  void setFirstChild(Int_t child)
  { fChild[0] = ( (child > std::numeric_limits<unsigned short>::max() ) ?
		  std::numeric_limits<unsigned short>::max() : (UShort_t)child ); }
  /// Set index of the second (last) child
  void setLastChild(Int_t child)
  { fChild[1] = ( (child > std::numeric_limits<unsigned short>::max() ) ?
		  std::numeric_limits<unsigned short>::max() : (UShort_t)child ); }
  /// Set px (GeV/c)
  void setPx(Double_t px)                                   { fPx = (Float_t)px; }
  /// Set py (GeV/c)
  void setPy(Double_t py)                                   { fPy = (Float_t)py; }
  /// Set pz (GeV/c)
  void setPz(Double_t pz)                                   { fPz = (Float_t)pz; }
  /// Set energy (GeV)
  void setEnergy(Double_t e)                                { fEnergy = (Float_t)e; }
  /// Set energy (GeV)
  void setE(Double_t e)                                     { setEnergy(e); }
  /// Set four-momentum (px, py, pz, E)
  void set4momentum(TLorentzVector mom)
  { setPx(mom.Px()); setPy(mom.Py()); setPz(mom.Pz()); setE(mom.E()); }
  /// Set start x
  void setX(Double_t x)                                     { fX = (Float_t)x; }
  /// Set start y
  void setY(Double_t y)                                     { fY = (Float_t)y; }
  /// Set start z
  void setZ(Double_t z)                                     { fZ = (Float_t)z; }
  /// Set freeze-out t (fm/c)
  void setT(Double_t t)                                     { fT = (Float_t)t; }
  /// Set four-coordinate (x, y, z, t)
  void set4coordinate(TLorentzVector vec)
  { setX(vec.X()); setY(vec.Y()); setZ(vec.Z()); setT(vec.T()); }
  /// Add index of MpdMiniTrack that was reconstructed out of current MC track
  void addGlobalTrackId(UShort_t id);
  /// Add indices of MpdMiniTracks that wwere reconstructed out of current MC track
  void setGlobalTrackIds(std::vector< UShort_t > ids)      { fRecoTrackIds = ids; }


  //
  // Getters
  //

  /// Return unique ID
  UShort_t id() const                         { return fId; }
  /// Return particle status
  Int_t status() const                        { return (Int_t)fStatus; }
  /// Return particle charge
  Float_t charge() const
  { return TDatabasePDG::Instance()->GetParticle( fPdgId )->Charge(); }
  /// Return PDG code
  Int_t pdgId() const                         { return fPdgId; }
  // Return parent particle index
  // UShort_t parentIndex() const                { return fParentIndex; }
  /// Return the first child index
  UShort_t firstChildIndex() const            { return fChild[0]; }
  /// Return the second (last) child index
  UShort_t lastChildIndex() const             { return fChild[1]; }
  /// Return px (GeV/c)
  Float_t px() const                          { return fPx; }
  /// Return py (GeV/c)
  Float_t py() const                          { return fPy; }
  /// Return pz (GeV/c)
  Float_t pz() const                          { return fPz; }
  /// Return energy (GeV)
  Float_t energy() const                      { return fEnergy; }
  /// Return energy (GeV)
  Float_t e() const                           { return energy(); }
  /// Return PDG mass (GeV/c^2)
  Float_t mass() const
  { return TDatabasePDG::Instance()->GetParticle( fPdgId )->Mass(); }
  /// Return energy estimated via PDG particle mass
  Float_t energyPdg() const
  { return TMath::Sqrt( p().Mag2() + mass()*mass() ); }
  /// Return three-momentum (px, py, pz)
  TVector3 p() const                          { return TVector3(fPx, fPy, fPz); }
  /// Return three-momentum (px, py, pz)
  TVector3 momentum() const                   { return p(); }
  /// Return four-momentum
  TLorentzVector fourMomentum() const         { return TLorentzVector(fPx, fPy, fPz, fEnergy); }
  /// Return freeze-out x coordinate (fm)
  Float_t x() const                           { return fX; }
  /// Return freeze-out y coordinate (fm)
  Float_t y() const                           { return fY; }
  /// Return freeze-out z coordinate (fm)
  Float_t z() const                           { return fZ; }
  /// Return freeze-out t coordinate (fm/c)
  Float_t t() const                           { return fT; }
  /// Return freeze-out four-coordinate (x, y, z, t)
  TLorentzVector fourCoordinate() const       { return TLorentzVector(fX, fY, fZ, fT); }
  /// Return indices of the MpdMiniTracks that were reconstructed
  /// from the current MC track
  std::vector< UShort_t > recoTrackIds() const { return fRecoTrackIds; }

 private:
  /// Unique track ID
  UShort_t fId;
  /// Status of the track
  Char_t   fStatus;
  /// PDG code
  Int_t    fPdgId;
  // Parent index
  //UShort_t fParentIndex;
  /// The first and the second daughter particle indeces (-1 if not decayed)
  UShort_t fChild[2];
  /// Px (GeV/c)
  Float_t  fPx;
  /// Py (GeV/c)
  Float_t  fPy;
  /// Pz (GeV/c)
  Float_t  fPz;
  /// Energy from the generator (GeV/c)
  Float_t  fEnergy;
  /// Freeze-out x coordinate (fm)
  Float_t  fX;
  /// Freeze-out y coordinate (fm)
  Float_t  fY;
  /// Freeze-out z coordinate (fm)
  Float_t  fZ;
  /// Freeze-out t coordinate (fm/c)
  Float_t  fT;

  /// Indices of the reconstructed global tracks that were
  /// reconstructed from the current Monte Carlo track
  /// Empty when no tracks were reconstructed from the MC track.
  std::vector< UShort_t > fRecoTrackIds;

  ClassDef(MpdMiniMcTrack, 2)
};

#endif // #define MpdMiniMcTrack_h
