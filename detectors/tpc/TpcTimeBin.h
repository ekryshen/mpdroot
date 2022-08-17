#ifndef _TPCTIMEBIN_H_
#define _TPCTIMEBIN_H_

#include <TObject.h>

class TpcTimeBin : public TObject {
private:
   Int_t fQ;

public:
   TpcTimeBin(void) : fQ(0) {}
   TpcTimeBin(Long_t anQ) : fQ(anQ) {}
   virtual ~TpcTimeBin(void) {}

   // *** Accessors ***
   Int_t GetQ() const { return fQ; }
   void  SetQ(Int_t q) { fQ = q; }

   // *** Operators ***
   bool operator<(const TpcTimeBin &rhs) const { return fQ < rhs.fQ; }
   bool operator>(const TpcTimeBin &rhs) const { return rhs < *this; }
   bool operator==(const TpcTimeBin &rhs) const { return !(*this < rhs) && !(*this > rhs); }
   bool operator!=(const TpcTimeBin &rhs) const { return (*this < rhs) || (*this > rhs); }

   static Int_t CompareQ(const TObject *obj1, const TObject *obj2)
   {
      TpcTimeBin *tb1 = (TpcTimeBin *)obj1;
      TpcTimeBin *tb2 = (TpcTimeBin *)obj2;
      return tb1->GetQ() - tb2->GetQ();
   }
   ClassDef(TpcTimeBin, 1);
};
#endif // _TPCTIMEBIN_H_
