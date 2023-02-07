#ifndef MPD_TPC_SECTOR_H
#define MPD_TPC_SECTOR_H

#include <TObject.h>
#include <TVector2.h>
#include <TF1.h>
#include <vector>
//#include "TpcPadPlaneSignals.h"
//#include "TpcAvalanche.h"
using namespace std;

class TpcSector : public TObject {
private:
   Float_t _sectorPhi; // angle of the middle of the sector

   Float_t _minX; // half of min sector width
   Float_t _maxX; // half of max sector width
   Float_t _midX; // half of mid sector width
   Float_t _outerPadWidth;
   Float_t _outerPadHeight;
   Float_t _innerPadWidth;
   Float_t _innerPadHeight;
   Float_t _dPhi2;
   Float_t _sectInnerHeight; // height of sector inner part
   Float_t _sectOuterHeight; // height of sector outer part
   Float_t _sectHeight;      // height of whole sector
   Float_t _tpcLength;       // length from electrode to padPlane

   Int_t *_numPadsInRow;
   Int_t *_numInnerPadsInRow;
   Int_t *_numOuterPadsInRow;

   Float_t _rMin;

   UInt_t _nRows;
   UInt_t _nInnerRows;
   UInt_t _nOuterRows;
   UInt_t _nSectors;
   UInt_t _nTimeBins;

public:
   // constructors and destructors
   TpcSector();
   virtual ~TpcSector();

   //        Float_t Phi() const { return _sectorPhi; }
   Float_t GetOuterPadHeight() const { return _outerPadHeight; }
   Float_t GetOuterPadWidth() const { return _outerPadWidth; }
   Float_t GetInnerPadHeight() const { return _innerPadHeight; }
   Float_t GetInnerPadWidth() const { return _innerPadWidth; }

   UInt_t GetNumRows() const { return _nRows; }
   UInt_t GetNumInnerRows() const { return _nInnerRows; }
   UInt_t GetNumOuterRows() const { return _nOuterRows; }

   UInt_t GetNumPadsInInnerRow(UInt_t row) const { return _numInnerPadsInRow[row]; }
   UInt_t GetNumPadsInOuterRow(UInt_t row) const { return _numOuterPadsInRow[row]; }
   UInt_t GetNumPadsInRow(UInt_t row) const { return _numPadsInRow[row]; }

   Float_t GetLength() const { return _tpcLength; }
   Float_t GetSectHeight() const { return _sectHeight; }
   Float_t GetSectInnerHeight() const { return _sectInnerHeight; }

   Int_t *GetArrayPadsInRow() const { return _numPadsInRow; }
   Int_t *GetArrayPadsInInnerRow() const { return _numInnerPadsInRow; }
   Int_t *GetArrayPadsInOuterRow() const { return _numOuterPadsInRow; }

   Float_t GetRmin() const { return _rMin; }

   UInt_t  GetNSectors() const { return _nSectors; }
   UInt_t  GetNTimeBins() const { return _nTimeBins; }
   Float_t GetMaxX() const { return _maxX; }
   Float_t GetMinX() const { return _minX; }
   Float_t GetMidX() const { return _midX; }

   ClassDef(TpcSector, 9);
};

#endif // TPC_SECTOR_H_
