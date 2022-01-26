#ifndef MPDEMCDETPOINT_H
#define MPDEMCDETPOINT_H 1


#include "FairMCPoint.h"

#include "TObject.h"
#include "TVector3.h"
using namespace std;

class MpdEmcPoint : public FairMCPoint
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
  MpdEmcPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t tof, Double_t length, Double_t ELoss);

  /** Copy constructor **/
  MpdEmcPoint(const MpdEmcPoint& point) { *this = point; };


  /** Destructor **/
  MpdEmcPoint();
  virtual ~MpdEmcPoint();

  /** Output to screen **/
  virtual void Print(const Option_t* opt) const;


  ClassDef(MpdEmcPoint,1)

};

#endif
