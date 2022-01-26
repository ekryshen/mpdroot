#ifndef MPDEMCDETPOINTKI_H
#define MPDEMCDETPOINTKI_H 1

#include "FairMCPoint.h"

#include "TObject.h"
#include "TVector3.h"
using namespace std;

class MpdEmcPointKI : public FairMCPoint
{
 public:
  /** Default constructor **/

  /** Constructor with arguments
   *@param trackID  Index of MCTrack
   *@param detID    Detector ID
   *@param pos      Ccoordinates at entrance to active volume [cm]
   *@param mom      Momentum of track at entrance [GeV]
   *@param tof      Time since event start [ns]
   *@param length   Track length since creation [cm]
   *@param ELoss    Energy deposit [GeV]
   **/
  MpdEmcPointKI(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t tof, Double_t length, Double_t ELoss);

  /** Copy constructor **/
  MpdEmcPointKI(const MpdEmcPointKI& point) { *this = point; };

  /** Destructor **/
  MpdEmcPointKI();
  virtual ~MpdEmcPointKI();

  /// \brief Check whether the points are from the same SuperParent and in the same detector volume
  /// \return True if points are the same (origin and detector), false otherwise
  Bool_t operator==(const MpdEmcPointKI& rhs) const;

  /// \brief Copy rhs to this MpdEmcPointKI
  Bool_t operator=(const MpdEmcPointKI& rhs) const;

  /// \brief Sorting points according to parent particle and detector volume
  /// \return True if this MpdEmcPointKI is smaller, false otherwise
  Bool_t operator<(const MpdEmcPointKI& rhs) const;

  /// \brief Adds energy loss from the other MpdEmcPointKI to this MpdEmcPointKI
  /// \param rhs phos::MpdEmcPointKI to add to this MpdEmcPointKI
  /// \return This MpdEmcPointKI with the summed energy loss
  MpdEmcPointKI& operator+=(const MpdEmcPointKI& rhs);

  /// \brief Creates a new MpdEmcPointKI based on this MpdEmcPointKI but adding the energy loss of the right hand side
  /// \param
  /// \return New MpdEmcPointKI based on this MpdEmcPointKI
  MpdEmcPointKI operator+(const MpdEmcPointKI& rhs) const;

  /// \brief To allow sorting
  Bool_t IsSortable() const { return kTRUE; }

  /// \brief Mathod ised for sorting Hits
  //  \param Another ModEmcPoint
  //  \return
  Int_t Compare(const TObject* obj) const;

  void AddEnergyLoss(Float_t eloss) { fELoss += eloss; }

  /** Output to screen **/
  virtual void Print(const Option_t* opt) const;

  ClassDef(MpdEmcPointKI, 1)
};

#endif
