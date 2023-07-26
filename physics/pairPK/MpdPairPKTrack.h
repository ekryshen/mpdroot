#ifndef MPDPAIRPKTRACK_H
#define MPDPAIRPKTRACK_H

#include "TLorentzVector.h"

class MpdPairPKTrack : public TLorentzVector {
public:
   MpdPairPKTrack() = default;
   MpdPairPKTrack(float px, float py, float pz, float e) : TLorentzVector(px, py, pz, e) {}
   ~MpdPairPKTrack() override = default;

   unsigned int pidWord() const { return mPID; }
   bool         pidBit(int iBit) const { return (mPID >> iBit) & 1; }

   void setPidWord(unsigned int w) { mPID = w; }
   void setPidBit(int iBit, int what = 1) { mPID ^= (-what ^ mPID) & (1 << iBit); }

   int  primary() const { return mPrimary; }
   void setPrimary(int p) { mPrimary = p; }

   void setTr1(int tr) { mTr1 = tr; }
   int  getTr1() { return mTr1; }

   void setPID(int pid) { mPID = pid; }
   int  getPID() { return mPID; }

private:
   int mPID     = -1;
   int mPrimary = -1;
   int mTr1     = -1;

   ClassDefOverride(MpdPairPKTrack, 1);
};
#endif