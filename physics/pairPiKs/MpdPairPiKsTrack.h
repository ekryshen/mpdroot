#ifndef MPDPAIRPIKSTRACK_H
#define MPDPAIRPIKSTRACK_H

#include "TLorentzVector.h"

class MpdPairPiKsTrack : public TLorentzVector {
public:
   MpdPairPiKsTrack() = default;
   MpdPairPiKsTrack(float px, float py, float pz, float e) : TLorentzVector(px, py, pz, e) {}
   ~MpdPairPiKsTrack() override = default;

   // unsigned int pidWord() const { return mPID; }
   // bool         pidBit(int iBit) const { return (mPID >> iBit) & 1; }

   // void setPidWord(unsigned int w) { mPID = w; }
   // void setPidBit(int iBit, int what = 1) { mPID ^= (-what ^ mPID) & (1 << iBit); }

   long int primary1() const { return mPrimary1; }
   void     setPrimary1(long int p) { mPrimary1 = p; }

   long int primary2() const { return mPrimary2; }
   void     setPrimary2(long int p) { mPrimary2 = p; }

   void     setTr1(long int tr) { mTr1 = tr; }
   long int getTr1() { return mTr1; }

   void     setTr2(long int tr) { mTr2 = tr; }
   long int getTr2() { return mTr2; }

   void setPID1(int pid) { mPID1 = pid; }
   int  getPID1() { return mPID1; }

   void setPID2(int pid) { mPID2 = pid; }
   int  getPID2() { return mPID2; }

   void setCh1(int charge) { mCh1 = charge; }
   int  getCh1() { return mCh1; }

   void setCh2(int charge) { mCh2 = charge; }
   int  getCh2() { return mCh2; }

private:
   int      mCh1      = 0;
   int      mCh2      = 0;
   int      mPID1     = -1;
   int      mPID2     = -1;
   long int mPrimary1 = -1;
   long int mPrimary2 = -1;
   long int mTr1      = -1;
   long int mTr2      = -1;

   ClassDefOverride(MpdPairPiKsTrack, 1);
};
#endif
