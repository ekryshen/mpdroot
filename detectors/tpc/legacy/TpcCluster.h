#ifndef _TPCCLUSTER_H_
#define _TPCCLUSTER_H_

#include <TObject.h>
#include <TClonesArray.h>
#include "TpcTimeBin.h"

class TpcSequence;
class TpcPadRow;

class TpcCluster : public TObject {
private:
   TpcPadRow *fPadRow;    // weak reference to padrow in that cluster appear
   TObjArray *fSequences; // array of TpcSequence objects
   TObjArray *fRegion;    // Cluster region

   // Cluster region bounds in pad row
   Int_t fLoPad; // included
   Int_t fHiPad; // excluded
   Int_t fLoBin; // included
   Int_t fHiBin; // excluded

   TpcTimeBin fBinStub;

public:
   TpcCluster(void);
   virtual ~TpcCluster(void);

   /** Accessors **/
   TpcPadRow *GetPadRow() const { return fPadRow; }
   void       SetPadRow(TpcPadRow *aPadRow) { fPadRow = aPadRow; }

   const TObjArray   *GetSequences() const { return fSequences; }
   Int_t              GetNSequences() const { return fSequences->GetEntriesFast(); }
   const TpcSequence *GetSequence(Int_t iSeq) const { return (TpcSequence *)fSequences->At(iSeq); }

   Int_t GetLoPad() const { return fLoPad; }
   Int_t GetHiPad() const { return fHiPad; }
   Int_t GetLoBin() const { return fLoBin; }
   Int_t GetHiBin() const { return fHiBin; }

   const TObjArray *GetRegion() const { return fRegion; }

   const TpcTimeBin *GetTimeBin(Int_t ipad, Int_t ibin);

   /** Modifiers **/
   void BeginCluster(TpcPadRow *aPadRow);
   void AddSequence(TpcSequence *aSequence);
   void EndCluster();

   const double ret_x() const { return 0; };
   const double ret_y() const { return 0; };
   const double ret_z() const { return 0; };
   const double ret_dx() const { return 0; };
   const double ret_dy() const { return 0; };
   const double ret_dz() const { return 0; };
   const double ret_charge() const { return 0; };

private:
   Bool_t RecalcBounds();
   void   MakeRegion();

   /** **/
   ClassDef(TpcCluster, 1);
};

#endif // _TPCCLUSTER_H_
