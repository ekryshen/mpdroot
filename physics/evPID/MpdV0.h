#ifndef MPDV0_H
#define MPDV0_H

#include "TLorentzVector.h"

class MpdV0 : public TLorentzVector {
public:
   MpdV0() = default;
   MpdV0(float px, float py, float pz, float e) : TLorentzVector(px, py, pz, e) {}
   MpdV0(MpdV0 &v0)  = default;
   ~MpdV0() override = default;

   float getBDTMom() { return mBDTmomentum; }
   void  setBDTMomentum(float p) { mBDTmomentum = p; }

   void setMatched1(int matched) { mMatched1 = matched; } // First matched track index
   void setMatched2(int matched) { mMatched2 = matched; } // Second matched track index
   int  getMatched1() { return mMatched1; }
   int  getMatched2() { return mMatched2; }

   void setCommonParent(int parent) { mCommonParent = parent; } // Index of matched particles parent (-1 for fake V0)
   int  getCommonParent() { return mCommonParent; }

   void getArmenteros(float &alpha, float &qt) const
   {
      alpha = mAlpha;
      qt    = mQt;
   } // Armenteros-Podolansky parameters
   void setArmenteros(float alpha, float qt)
   {
      mAlpha = alpha;
      mQt    = qt;
   }

   void getAsymmetry(float &asym1, float &asym2) const
   {
      asym1 = mAsym1;
      asym2 = mAsym2;
   }
   void setAsymmetry(float asym1, float asym2)
   {
      mAsym1 = asym1;
      mAsym2 = asym2;
   }

   void  setBDTValue(float v) { mBDTvalue = v; } // value estimated with TMVA BDT estimator
   float getBDTValue() { return mBDTvalue; }

   float getChi2() const { return mChi2; } // chi2 provided by Kalman fit
   void  setChi2(float chi2) { mChi2 = chi2; }

   float getDaughterDCA() const { return mDDCA; } // Minimal distance bewtween daughters
   void  setDaughterDCA(float dist) { mDDCA = dist; }

   float getMass() const { return mMass; } // calculated mass of the pair
   void  setMass(float m) { mMass = m; }

   float getCPA() const
   {
      return mCPA;
   } // cos angle between momentum and direction from primary vertex to creation point
   void setCPA(float a) { mCPA = a; }

   float getCospsi() const { return mCospsi; } // Pair orientation vrt mag field
   void  setCospsi(float psi) { mCospsi = psi; }

   float getRconv() const { return sqrt(mX * mX + mY * mY); } // Conversion radius

   void setTr1(int tr) { mTr1 = tr; } // First track index
   void setTr2(int tr) { mTr2 = tr; } // Second track index
   int  getTr1() { return mTr1; }
   int  getTr2() { return mTr2; }

   void getVertex(float &x, float &y, float &z)
   {
      x = mX;
      y = mY;
      z = mZ;
   } // Conversion position
   void setVertex(float x, float y, float z)
   {
      mX = x;
      mY = y;
      mZ = z;
   }

private:
   long int mMatched1     = -1;   // Index of MC particle matched to the first track
   long int mMatched2     = -1;   // Index of MC particle matched to the second track
   long int mCommonParent = -1;   // Index of matched particles common parent (-1 for fake V0)
   int      mTr1          = -1;   // First track index
   int      mTr2          = -1;   // Second track index
   float    mAlpha        = 0;    // Armenteros-Podolansky parameters
   float    mAsym1        = 0;    // Momentum asymmetry track 1
   float    mAsym2        = 0;    // Momentum asymmetry track 2
   float    mQt           = 0.;   // Armenteros-Podolansky parameters
   float    mCPA          = 0.;   // cos angle between momentum and direction from primary vertex to creation point
   float    mCospsi       = 0.;   // Pair orientation vrt mag field
   float    mChi2         = 999.; // chi2 provided by Kalman fit
   float    mDDCA         = 999.; // Minimal distance bewtween daughters
   float    mMass         = 0.;   // calculated mass of the pair
   float    mX            = 0;    // Conversion position
   float    mY            = 0;    // Conversion position
   float    mZ            = 0;    // Conversion position
   float    mBDTvalue     = 0.;   // Value estimated from TMVA BDT classificator
   float    mBDTmomentum  = 0.;   // Momentrum estimated from BDT regression

   ClassDefOverride(MpdV0, 1);
};
#endif
