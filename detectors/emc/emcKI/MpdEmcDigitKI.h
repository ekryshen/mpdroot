#ifndef MPDEMCDIGITKI_H
#define MPDEMCDIGITKI_H 1

#include <map>
#include "TObject.h"
class MpdEmcPointKI;

class MpdEmcDigitKI : public TObject
{
 public:
  /** Default constructor **/
  MpdEmcDigitKI();

  MpdEmcDigitKI(MpdEmcPointKI* point);

  MpdEmcDigitKI(int cellId, float energy, float time, int trackId);

  virtual ~MpdEmcDigitKI();

  // Check, if point can be added (hits from the same Tower)
  Bool_t CanAdd(MpdEmcPointKI* point) const;

  // Adds point (add energy, change time if necessary, add primary)
  void AddPoint(MpdEmcPointKI* point);

  int GetDetId() const { return fDetId; }
  Float_t GetE() const { return fE; }

  Float_t GetTime() const { return fTime; };

  void SetE(Double_t e) { fE = e; }
  void SetTime(Double_t t) { fTime = t; }

  void Print(const Option_t* opt = 0) const;

  // Methods used for event display
  double GetZcenter();
  double GetPhiCenter();

  /// To allow sorting
  Bool_t IsSortable() const { return kTRUE; }

  /// \brief Method ised for sorting Hits
  //  \param Another ModEmcDigit
  //  \return
  Int_t Compare(const TObject* obj) const;

  Int_t GetNPrimaries() const { return fNprimary; }

  void GetPrimary(Int_t i, Int_t& primId, Float_t& primEdep) const
  {
    if (i >= fNprimary) {
      primId = -1;
      primEdep = 0;
      return;
    } else {
      primId = fPrimary[i];
      primEdep = fPrimEdep[i];
    }
  }

 protected:
  UInt_t fDetId;      // Tower index
  Float_t fE;         // Full energy
  Float_t fTime;      // hit mean time
  Int_t fNprimary;    // Number of primaries
  Int_t* fPrimary;    //[fNprimary] Array of primaries
  Float_t* fPrimEdep; //[fNprimary] Array of deposited energies

  std::map<int, float> fEdep; // energy deposit for each track ID

  ClassDef(MpdEmcDigitKI, 1)
};

#endif
