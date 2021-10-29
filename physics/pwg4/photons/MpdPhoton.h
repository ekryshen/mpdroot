#ifndef MPDPHOTON_H
#define MPDPHOTON_H

#include "TLorentzVector.h"


class MpdPhoton: public TLorentzVector{
  public: 
    MpdPhoton() = default ;
    MpdPhoton(float px, float py, float pz, float e):TLorentzVector(px,py,pz,e) {}
    ~MpdPhoton() override = default ;

    unsigned int pidWord() const {return mPID ;}
    bool pidBit(int iBit)const { return (mPID>>iBit)&1 ;} 

    void setPidWord(unsigned int w) {mPID = w; }
    void setPidBit(int iBit, int what = 1) {mPID ^= (-what ^ mPID) & (1 << iBit); }

    int primary() const { return mPrimary ;}
    void setPrimary(int p) { mPrimary = p ;}

    void setTr1(int tr){mTr1 = tr ;}
    void setTr2(int tr){mTr2 = tr ;}
    int getTr1(){ return mTr1 ;}
    int getTr2(){ return mTr2 ;}

  private:
    unsigned int mPID = 0 ; // use to keep PID cuts   
    int mPrimary = 0; 
    int mTr1 = -1; 
    int mTr2 = -1; 

  ClassDefOverride(MpdPhoton, 1);

} ;
#endif 
