//
// Joint Institute for Nuclear Research (JINR), Dubna, Russia, Nov-2021
// MPD TPC Clustering objects for NICA project
// author: Viktor A.Krylov JINR/LNP
// email: kryman@jinr.ru
//
#ifndef TPC_CLUSTERING_H
#define TPC_CLUSTERING_H

#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <ostream>
#include <cassert>
#include <math.h>

#include <string>
#include <vector>
#include <array>
#include <list>
#include <algorithm>
#include <iterator>
#include <climits>

#define AUXHIT_DEBUG
#define HIT_ADJUSTMENT

#define _RED_(_MSG_) "\e[1;31m" _MSG_ "\e[0m"
#define _GREEN_(_MSG_) "\e[1;32m" _MSG_ " \e[0m"
#define _YELLOW_(_MSG_) "\e[1;33m" _MSG_ "\e[0m"
//
// A Fast, Accurate, and Separable Method for Fitting a Gaussian Function
// Ibrahim Al-Nahhal, Octavia A. Dobre, Ertugrul Basar, Cecilia Moloney, Salama Ikki
// DOI:	10.1109/MSP.2019.2927685
// arXiv:1907.07241v1
//
// GUO’S algorithm without iterative procedure
// Solution of AX = B -> {A = LU} -> LUX = B -> {UX = Y} -> LY = B -> Y
// UX = Y -> X
//            | Sum(y^2)      Sum(x*y^2)     Sum(x^2*y^2) |
// MATRIX A = | Sum(x*y^2)    Sum(x^2*y^2)   Sum(x^3*y^2) |
//            | Sum(x^2*y^2)  Sum(x^3*y^2)   Sum(x^4*y^2) |
//
//            | Sum(y^2*Log(y))     |
// MATRIX B = | Sum(x*y^2*Log(y))   |
//            | Sum(x^2*y^2*Log(y)) |
//
//     | l11  0   0  |         | l11 l21 l31 |
// L = | l21 l22  0  |     U = |  0  l22 l32 |
//     | l31 l32 l33 |         |  0   0  l33 |
//
template <typename X, typename Y>
struct XY {
   X x;
   Y y;
};
template <typename X, typename Y>
double calcMatrix(const std::vector<XY<X, Y>> &rvXY, std::array<double, 5> &radSum, std::array<double, 3> &radSumLog)
{
   double dSum0 = 0, dSum1 = 0, dSum2 = 0, dSum3 = 0, dSum4 = 0;
   double dSumLog0 = 0, dSumLog1 = 0, dSumLog2 = 0;
   double yMax  = 0;
   double xMax  = 0;
   size_t nSize = rvXY.size();
   for (uint i = 0; i < nSize; i++) {
      X x  = rvXY[i].x;
      X x2 = x * x;
      X x3 = x2 * x;
      X x4 = x3 * x;

      double y = rvXY[i].y;
      if (y > yMax) {
         yMax = y;
         xMax = x;
      }
      double y2 = y * y;

      dSum0 += y2;
      dSum1 += x * y2;
      dSum2 += x2 * y2;
      dSum3 += x3 * y2;
      dSum4 += x4 * y2;

      double dYLog = y2 * log(y);
      dSumLog0 += dYLog;
      dSumLog1 += x * dYLog;
      dSumLog2 += x2 * dYLog;
   }
   radSum    = {dSum0, dSum1, dSum2, dSum3, dSum4};
   radSumLog = {dSumLog0, dSumLog1, dSumLog2};
   return xMax;
}
const double  dSQRT_2_PI = sqrt(M_PI_2);
inline double calcGaussArea(double x, double y, double d, double xMin, double xMax)
{                                // Gauss area from xMin to xMax
   double d2      = M_SQRT2 * d; // deviation
   double fErfMin = std::erf((xMin - x) / d2);
   double fErfMax = std::erf((xMax - x) / d2);
   return dSQRT_2_PI * y * d * (fErfMax - fErfMin);
}
template <typename X, typename Y>
std::array<double, 4> calcGaussRes(const std::vector<XY<X, Y>> &rvXY)
{
   assert(rvXY.size() > 2);

   std::array<double, 5> adSum;
   std::array<double, 3> adSumLog;
   double                xMax = calcMatrix<X, Y>(rvXY, adSum, adSumLog);

   double dL11 = sqrt(adSum[0]);
   double dL21 = adSum[1] / dL11;
   double dL22 = sqrt(adSum[2] - dL21 * dL21);
   double dL31 = adSum[2] / dL11;
   double dL32 = (adSum[3] - dL21 * dL31) / dL22;
   double dL33 = sqrt(adSum[4] - dL31 * dL31 - dL32 * dL32);

   double dY1 = adSumLog[0] / dL11;
   double dY2 = (adSumLog[1] - dL21 * dY1) / dL22;
   double dY3 = (adSumLog[2] - dL31 * dY1 - dL32 * dY2) / dL33;

   double dX2 = dY3 / dL33;
   double dX1 = (dY2 - dL32 * dX2) / dL22;
   double dX0 = (dY1 - dL31 * dX2 - dL21 * dX1) / dL11;

   double dX25 = -0.5 / dX2;
   double x    = dX1 * dX25;
   double y    = exp(dX0 + 0.5 * dX1 * x);
   double d    = dX25 > 0 ? sqrt(dX25) : 0; // deviation
   return {x, y, d, xMax};
}
template <typename X, typename Y>
double calcRsquared(const std::vector<XY<X, Y>> &rvXY, double x, double y, double d, double dYSum)
{ // coefficient of determination
   size_t nSize = rvXY.size();
   if (nSize > 2) { // deviation should be defined!!!
      double dYMean       = dYSum / nSize;
      double d2Deviation2 = 2 * d * d;
      double dSumTotal    = 0; // Sum total = Sum[(y - ymean)^2]
      double dSumRes      = 0; // Sum of residuals = Sum[(y - f(x))^2]
      for (uint i = 0; i < nSize; i++) {
         double dX     = rvXY[i].x - x;
         double dGauss = y * exp(-dX * dX / d2Deviation2);
         double dY     = rvXY[i].y;
         double dTot   = dY - dYMean;
         double dRes   = dY - dGauss;

         dSumTotal += dTot * dTot;
         dSumRes += dRes * dRes;
      }
      return 1 - dSumRes / dSumTotal;
   } else
      return 0;
}
namespace tpcClustering {
struct Gauss {
   std::array<double, 4>
          _array;     // array[0: mean value, 1: gaussian peak value (ADC), 2: standard deviation, 3: adjustment]
   double _dArea;     // gaussian area of the normal distribution (should be defined for reference only!)
   double _dRsquared; // r-squared value of the normal distribution (should be defined for reference only!)
public:
   Gauss() : _array{0, 0, 0, 0}, _dArea(0), _dRsquared(0) {}
   Gauss(const Gauss &r) : _dArea(r._dArea), _dRsquared(r._dRsquared) { _array = r._array; }
   Gauss &operator=(const Gauss &r)
   {
      if (this != &r) {
         _array     = r._array;
         _dArea     = r._dArea;
         _dRsquared = r._dRsquared;
      }
      return *this;
   }
   void Set(const std::array<double, 4> &r) { _array = r; }

   void   setMean(double d) { _array[0] = d; }
   double getMean() const { return _array[0]; }

   void   setMax(double d) { _array[1] = d; }
   double getMax() const { return _array[1]; }

   void   setDeviation(double d) { _array[2] = d; }
   double getDeviation() const { return _array[2]; }

   double getAdjustment() const
   {
      return _array[3];
   } // adjustment value for unreasonable fit results to max AdcSum value!!!

   void   setArea(double d) { _dArea = d; }
   double getArea() const { return _dArea; }

   void   setRsquared(double d) { _dRsquared = d; }
   double getRsquared() const { return _dRsquared; }

   std::ostringstream &Print(std::ostringstream &oss) const
   {
      oss << "hit:{ mean: " << getMean() << ", max: " << getMax() << ", deviation: " << getDeviation()
          << ", adjustment: " << getAdjustment() << ", area: " << _dArea << ", r-squared: " << _dRsquared << " }"
          << std::endl;
      return oss;
   }
};
class AdcHit {
   uint  _nTimeBin; // timebin number [0..310)
   float _fAdc;
   int   _nTrackId; //_nTrackId
public:
   AdcHit() : _fAdc(0), _nTimeBin(0), _nTrackId(-1) {} //_nTrackId
   AdcHit(int n, float f) : _nTimeBin(n), _fAdc(f), _nTrackId(-1) {}
   AdcHit(int n, float f, int nId) : _nTimeBin(n), _fAdc(f), _nTrackId(nId) {} //_nTrackId
   AdcHit(const AdcHit &r) : _nTimeBin(r._nTimeBin), _fAdc(r._fAdc), _nTrackId(r._nTrackId) {}
   ~AdcHit() {}
   const AdcHit &operator=(const AdcHit &r)
   {
      if (this != &r) {
         _nTimeBin = r._nTimeBin;
         _fAdc     = r._fAdc;
         _nTrackId = r._nTrackId; //_nTrackId
      }
      return *this;
   }
   inline float getTrackID() const { return _nTrackId; } //_nTrackId
   inline float getAdc() const { return _fAdc; }
   inline void  setAdc(float f) { _fAdc = f; }
   inline uint  getTimeBin() const { return _nTimeBin; }
   inline int   getHitWay(float fDistance, uint nBins) const
   {
      return _nTimeBin * fDistance / nBins;
   } // distance: m, nBins: 310
   inline int getHitTime(float fDistance, uint nBins, float fRate) const
   {
      return getHitWay(fDistance, nBins) / fRate;
   } // rate: m/c
   std::ostringstream &Print(std::ostringstream &oss) const
   {
      oss << _nTimeBin << ":" << _fAdc;
      return oss;
   }
   friend std::ostringstream &operator<<(std::ostringstream &, const AdcHit &);
   friend std::ostringstream &operator<<(std::ostringstream &, const AdcHit *);
};
inline std::ostringstream &operator<<(std::ostringstream &oss, const AdcHit &r)
{
   oss << r.getTimeBin() << ":" << r.getAdc();
   return oss;
}
inline std::ostringstream &operator<<(std::ostringstream &oss, const AdcHit *p)
{
   oss << p->getTimeBin() << ":" << p->getAdc();
   return oss;
}

class PadHit {
   uint  _nPad;       // Pad id number
   float _fAdcHitSum; // total ADC sum
public:
   PadHit() : _nPad(UINT_MAX), _fAdcHitSum(0) {}
   PadHit(uint nPad) : _nPad(nPad), _fAdcHitSum(0) {}
   PadHit(uint nPad, const AdcHit &r) : _nPad(nPad) { _fAdcHitSum = r.getAdc(); }
   PadHit(const PadHit &r) : _nPad(r._nPad), _fAdcHitSum(r._fAdcHitSum) {}
   ~PadHit() {}
   PadHit &operator=(const PadHit &r)
   {
      if (this != &r) {
         _nPad       = r._nPad;
         _fAdcHitSum = r._fAdcHitSum;
      }
      return *this;
   }
   PadHit &operator=(const AdcHit &r)
   {
      ((AdcHit *)this)->operator=(r);
      return *this;
   }
   inline uint         getPad() const { return _nPad; }
   inline uint         getTimeBin() const { return _nPad; }
   inline void         setAdcHitSum(float f) { _fAdcHitSum = f; }
   inline float        getAdcHitSum() const { return _fAdcHitSum; }
   inline void         AddAdc(float f) { _fAdcHitSum += f; }
   inline void         ReduceAdcSum(float f) { _fAdcHitSum -= f; }
   std::ostringstream &Print(std::ostringstream &oss) const
   {
      oss << "pad: " << _nPad << ", AdcSum: " << _fAdcHitSum << std::endl;
      return oss;
   }
   std::ostringstream &Print2(std::ostringstream &oss) const
   {
      oss << "time: " << _nPad << ", AdcSum: " << _fAdcHitSum << std::endl;
      return oss;
   }
};
inline const char *SZ_ALPHABET = "ABCDEFJHIJKLMNOPQRSTUVWXYZ";
class Cluster;
class PadCluster : public PadHit {
   Cluster          *_pCluster;          // link to a Cluster where the PadCluster belongs
   const PadCluster *_pPadClusterBefore; // link to prior PadCluster (!= null: combined cluster)
   const PadCluster *_pPadClusterAfter;  // link to next PadCluster (!= null: combined cluster)

   uint _nTimeBinMin; // timebin min range
   uint _nTimeBinMax; // timebin max range

   float _fAdcMax; // max adc value in the AdcHit vector
   float _fAdcMin; // min adc value in the AdcHit vector

   uint _nTimeBinLapMin; // overlapping timebin min range in neighbour pad
   uint _nTimeBinLapMax; // overlapping timebin max range in neighbour pad

   uint _iAdcMax; // index of the max adc value in the AdcHit vector
   uint _iAdcMin; // index of the last min adc value in the AdcHit vector

   std::vector<AdcHit> _vAdcHits; // vector of ADC values
public:
   PadCluster()
      : PadHit(), _pCluster(NULL), _pPadClusterBefore(NULL), _pPadClusterAfter(NULL), _nTimeBinMin(UINT_MAX),
        _nTimeBinMax(UINT_MAX), _fAdcMin(0), _fAdcMax(0), _nTimeBinLapMin(0), _nTimeBinLapMax(0), _iAdcMin(0),
        _iAdcMax(0)
   {
   }
   PadCluster(int nPad)
      : PadHit(nPad), _pCluster(NULL), _pPadClusterBefore(NULL), _pPadClusterAfter(NULL), _nTimeBinMin(UINT_MAX),
        _nTimeBinMax(UINT_MAX), _fAdcMin(0), _fAdcMax(0), _nTimeBinLapMin(0), _nTimeBinLapMax(0), _iAdcMin(0),
        _iAdcMax(0)
   {
   }
   PadCluster(int nPad, const AdcHit &r)
      : PadHit(nPad, r), _pCluster(NULL), _pPadClusterBefore(NULL), _pPadClusterAfter(NULL), _nTimeBinLapMin(0),
        _nTimeBinLapMax(0), _iAdcMin(0), _iAdcMax(0)
   {
      _nTimeBinMin = _nTimeBinMax = r.getTimeBin();
      _fAdcMin = _fAdcMax = r.getAdc();
      _vAdcHits.push_back(r);
   }
   PadCluster(int nPad, const AdcHit &r, Cluster *pCluster)
      : PadHit(nPad, r), _pCluster(pCluster), _pPadClusterBefore(NULL), _pPadClusterAfter(NULL), _nTimeBinLapMin(0),
        _nTimeBinLapMax(0), _iAdcMin(0), _iAdcMax(0)
   {
      _nTimeBinMin = _nTimeBinMax = r.getTimeBin();
      _fAdcMin = _fAdcMax = r.getAdc();
      _vAdcHits.push_back(r);
   }
   PadCluster(const PadCluster &r)
      : PadHit(r), _pCluster(r._pCluster), _pPadClusterBefore(r._pPadClusterBefore),
        _pPadClusterAfter(r._pPadClusterAfter), _nTimeBinMin(r._nTimeBinMin), _nTimeBinMax(r._nTimeBinMax),
        _fAdcMin(r._fAdcMin), _fAdcMax(r._fAdcMax), _nTimeBinLapMin(r._nTimeBinLapMin),
        _nTimeBinLapMax(r._nTimeBinLapMax), _iAdcMin(r._iAdcMin), _iAdcMax(r._iAdcMax)
   {
      _vAdcHits = r._vAdcHits;
   }
   virtual ~PadCluster() { _vAdcHits.clear(); }
   AdcHit operator[](size_t n) const
   {
      assert(n < _vAdcHits.size());
      return _vAdcHits[n];
   }
   const AdcHit &operator[](size_t n)
   {
      assert(n < _vAdcHits.size());
      return _vAdcHits[n];
   }
   inline uint  getClusterId(const Cluster *pCluster) const;
   inline float getClusterTime(const Cluster *pCluster) const;
   inline float getClusterPad(const Cluster *pCluster) const;

   inline const std::vector<AdcHit> &getAdcHits() const { return _vAdcHits; } //!!!!!_nTrackId

   inline size_t   getSize() const { return _vAdcHits.size(); }
   inline Cluster *getCluster() const { return _pCluster; }
   inline void     setCluster(Cluster *pCluster) { _pCluster = pCluster; }
   inline AdcHit   getAdcHitMax() const { return _vAdcHits.empty() ? AdcHit(0, 0.0) : _vAdcHits[_iAdcMax]; }
   inline AdcHit   getAdcHitMin() const { return _vAdcHits.empty() ? AdcHit(0, 0.0) : _vAdcHits[_iAdcMin]; }
   inline uint     getTimeBinMin() const { return _nTimeBinMin; }
   inline uint     getTimeBinMax() const { return _nTimeBinMax; }
   inline float    getAdcMin() const { return _fAdcMin; }
   inline float    getAdcMax() const { return _fAdcMax; }
   inline uint     getTimeBinLapMin() const { return _nTimeBinLapMin; }
   inline uint     getTimeBinLapMax() const { return _nTimeBinLapMax; }

   inline const PadCluster *getPadClusterBefore() const { return _pPadClusterBefore; }
   inline void              setPadClusterBefore(const PadCluster *pPadCluster) { _pPadClusterBefore = pPadCluster; }
   inline const PadCluster *getPadClusterAfter() const { return _pPadClusterAfter; }
   inline void              setPadClusterAfter(const PadCluster *pPadCluster) { _pPadClusterAfter = pPadCluster; }

   inline void setPadTimeOverlapping(const PadCluster *pPadCluster)
   {
      assert(!(isBefore(pPadCluster) || isAfter(pPadCluster)));

      _nTimeBinLapMin = _nTimeBinMin > pPadCluster->getTimeBinMin() ? _nTimeBinMin : pPadCluster->getTimeBinMin();
      _nTimeBinLapMax = _nTimeBinMax < pPadCluster->getTimeBinMax() ? _nTimeBinMax : pPadCluster->getTimeBinMax();
   }
   inline bool isOut(const AdcHit &r) const { return r.getTimeBin() > _nTimeBinMax + 1; }
   inline bool isCombined() const { return _pPadClusterBefore != NULL || _pPadClusterAfter != NULL; }
   inline bool isBelongTo(const PadCluster *pPadCluster) const
   {
      if (isCombined() && pPadCluster->isCombined()) { // 1:1
         uint nTimeBinMin    = _pPadClusterBefore == NULL ? _nTimeBinMin : _nTimeBinMin + 1;
         uint nTimeBinMax    = _pPadClusterAfter == NULL ? _nTimeBinMax : _nTimeBinMax - 1;
         uint nTimeBinMinExt = pPadCluster->getPadClusterBefore() == NULL ? pPadCluster->getTimeBinMin()
                                                                          : pPadCluster->getTimeBinMin() + 1;
         uint nTimeBinMaxExt =
            pPadCluster->getPadClusterAfter() == NULL ? pPadCluster->getTimeBinMax() : pPadCluster->getTimeBinMax() - 1;
         uint nAdcHitMaxTimeBin    = getAdcHitMax().getTimeBin();
         uint nAdcHitMaxTimeBinExt = pPadCluster->getAdcHitMax().getTimeBin();

         bool bMin       = nTimeBinMin < nTimeBinMinExt && pPadCluster->getPadClusterBefore() != NULL;
         bool bMax       = nTimeBinMax > nTimeBinMaxExt && pPadCluster->getPadClusterAfter() != NULL;
         bool bMinMax    = bMin || bMax;
         bool bMinExt    = nTimeBinMinExt < nTimeBinMin && _pPadClusterBefore != NULL;
         bool bMaxExt    = nTimeBinMaxExt > nTimeBinMax && _pPadClusterAfter != NULL;
         bool bMinMaxExt = bMinExt || bMaxExt;
         if (bMinMax && bMinMaxExt)
            return nAdcHitMaxTimeBin >= nTimeBinMinExt && nAdcHitMaxTimeBin <= nTimeBinMaxExt &&
                   nAdcHitMaxTimeBinExt >= nTimeBinMin && nAdcHitMaxTimeBinExt <= nTimeBinMax;
         else if (bMinMax)
            return nAdcHitMaxTimeBin >= nTimeBinMinExt && nAdcHitMaxTimeBin <= nTimeBinMaxExt;
         else if (bMinMaxExt)
            return nAdcHitMaxTimeBinExt >= nTimeBinMin && nAdcHitMaxTimeBinExt <= nTimeBinMax;
         else
            return true;
      } else if (isCombined()) { // 1:0
         uint nTimeBinMin          = _pPadClusterBefore == NULL ? _nTimeBinMin : _nTimeBinMin + 1;
         uint nTimeBinMax          = _pPadClusterAfter == NULL ? _nTimeBinMax : _nTimeBinMax - 1;
         uint nAdcHitMaxTimeBinExt = pPadCluster->getAdcHitMax().getTimeBin();
         return (_pPadClusterBefore == NULL && pPadCluster->getTimeBinMax() <= _nTimeBinMax) ||
                (_pPadClusterAfter == NULL && pPadCluster->getTimeBinMin() >= _nTimeBinMin) ||
                (nAdcHitMaxTimeBinExt >= nTimeBinMin && nAdcHitMaxTimeBinExt <= nTimeBinMax); // max hit
      } else if (pPadCluster->isCombined()) {                                                 // 0:1
         uint nTimeBinMinExt = pPadCluster->getPadClusterBefore() == NULL ? pPadCluster->getTimeBinMin()
                                                                          : pPadCluster->getTimeBinMin() + 1;
         uint nTimeBinMaxExt =
            pPadCluster->getPadClusterAfter() == NULL ? pPadCluster->getTimeBinMax() : pPadCluster->getTimeBinMax() - 1;
         uint nAdcHitMaxTimeBin = getAdcHitMax().getTimeBin();
         return (pPadCluster->getPadClusterBefore() == NULL && _nTimeBinMax <= nTimeBinMaxExt) ||
                (pPadCluster->getPadClusterAfter() == NULL && _nTimeBinMin >= nTimeBinMinExt) ||
                (nAdcHitMaxTimeBin >= nTimeBinMinExt && nAdcHitMaxTimeBin <= nTimeBinMaxExt); // max hit
      } else                                                                                  // 0:0
         return true;
   }
   inline bool isAfter(const PadCluster *pPadCluster) const { return _nTimeBinMin > pPadCluster->getTimeBinMax(); }
   inline bool isBefore(const PadCluster *pPadCluster) const { return _nTimeBinMax < pPadCluster->getTimeBinMin(); }
   AdcHit      splitAdcHit(const AdcHit &rAdcHit)
   { // proportionally split AdcHit value between the two combined clusters;
      std::size_t nSize = _vAdcHits.size();
      assert(nSize > 1);

      AdcHit adcHit   = _vAdcHits.back();
      float  adcMid   = adcHit.getAdc();
      float  adcLeft  = _vAdcHits[nSize - 2].getAdc() - adcMid;
      float  adcRight = rAdcHit.getAdc() - adcMid;
      assert(adcLeft > 0 && adcRight > 0);

      float fRatio = 1 + adcRight / adcLeft; // coefficient of proportionality (ratio: left to right)
      _vAdcHits.back().setAdc(adcMid / fRatio);
      adcHit.setAdc(adcMid - _vAdcHits.back().getAdc());
      PadHit::ReduceAdcSum(adcHit.getAdc());
      return adcHit;
   }
   inline bool isCombined(const AdcHit &r) const { return _iAdcMax < _iAdcMin && r.getAdc() > getAdcHitMin().getAdc(); }
   inline void Add(const AdcHit &r)
   {
      int nTimeBin = r.getTimeBin();
      if (_nTimeBinMin == UINT_MAX && _nTimeBinMax == UINT_MAX)
         _nTimeBinMin = _nTimeBinMax = nTimeBin;
      else if (nTimeBin < _nTimeBinMin)
         _nTimeBinMin = nTimeBin;
      else if (nTimeBin > _nTimeBinMax)
         _nTimeBinMax = nTimeBin;

      float fAdc = r.getAdc();
      if (_fAdcMin == 0 && _fAdcMax == 0)
         _fAdcMin = _fAdcMax = fAdc;
      else if (fAdc < _fAdcMin)
         _fAdcMin = fAdc;
      else if (fAdc > _fAdcMax)
         _fAdcMax = fAdc;

      size_t nSize = _vAdcHits.size();
      if (_iAdcMax > _iAdcMin) {
         if (fAdc > getAdcHitMax().getAdc())
            _iAdcMax = nSize;
         else
            _iAdcMin = nSize;
      } else {
         if (fAdc < getAdcHitMin().getAdc())
            _iAdcMin = nSize;
         else
            _iAdcMax = nSize;
      }
      PadHit::AddAdc(fAdc);
      _vAdcHits.push_back(r);
   }
   std::ostringstream &Draw(std::ostringstream &oss, uint nTimeBinMin, uint nTimeBinMax, char ch = '+') const
   {
      uint   nTimeBin = _pCluster ? lroundf(getClusterTime(_pCluster)) : UINT_MAX;
      uint   nPad     = _pCluster ? lroundf(getClusterPad(_pCluster)) : UINT_MAX;
      uint   nMax     = getAdcHitMax().getTimeBin();
      size_t nSize    = _vAdcHits.size();
      for (uint i = nTimeBinMin, j = 0; i <= nTimeBinMax; i++)
         if (j < nSize && i == _vAdcHits[j].getTimeBin()) {
            if (i == nTimeBin && nPad == PadHit::getPad())
               oss << _RED_("*"); // red cluster hit star '*'
            else if (i == nMax && nSize > 1)
               oss << '^';
            else
               oss << ch;
            j++;
         } else
            oss << ' ';
      return oss;
   }
   std::ostringstream &Draw2(std::ostringstream &oss, uint nPadMin, uint nPadMax, char ch = '+') const
   {
      uint   nTimeBin = _pCluster ? lroundf(getClusterTime(_pCluster)) : UINT_MAX;
      uint   nPad     = _pCluster ? lroundf(getClusterPad(_pCluster)) : UINT_MAX;
      uint   nMax     = getAdcHitMax().getTimeBin();
      size_t nSize    = _vAdcHits.size();
      for (uint i = nPadMin, j = 0; i <= nPadMax; i++)
         if (j < nSize && i == _vAdcHits[j].getTimeBin()) {
            if (i == nPad && nTimeBin == PadHit::getPad())
               oss << _RED_("*"); // red cluster hit star '*'
            else if (i == nMax && nSize > 1)
               oss << '^';
            else
               oss << ch;
            j++;
         } else
            oss << ' ';
      return oss;
   }
   std::ostringstream &Draw(std::ostringstream &oss, char ch = '+') const
   {
      uint   nTimeBin = _pCluster ? lroundf(getClusterTime(_pCluster)) : UINT_MAX;
      uint   nPad     = _pCluster ? lroundf(getClusterPad(_pCluster)) : UINT_MAX;
      uint   nMax     = getAdcHitMax().getTimeBin();
      size_t nSize    = _vAdcHits.size();
      for (uint i = 0; i < nSize; i++) {
         uint n = _vAdcHits[i].getTimeBin();
         if (n == nTimeBin && nPad == PadHit::getPad())
            oss << _RED_("*"); // red cluster hit star '*'
         else if (n == nMax && nSize > 1)
            oss << '^';
         else
            oss << ch;
      }
      return oss;
   }
   std::ostringstream &Print(std::ostringstream &oss) const
   {
      PadHit::Print(oss);
      uint       nClusterId  = getClusterId(_pCluster);
      char       chClusterId = _pCluster == NULL ? '.' : SZ_ALPHABET[nClusterId % strlen(SZ_ALPHABET)];
      const uint nTab        = 4;
      oss << std::string(nTab, ' ') << "PadCluster:: time: [" << _nTimeBinMin << "," << _nTimeBinMax << "]"
          << ", adcSpikes: { min[" << _iAdcMin << "], max[" << _iAdcMax << "] }"
          << ", combined: {" << (_pPadClusterBefore != NULL ? '+' : '.') << (_pPadClusterAfter != NULL ? '+' : '.')
          << "}"
          << ", overlapping: [" << _nTimeBinLapMin << "," << _nTimeBinLapMax << "]"
          << ", cluster: " << nClusterId << " '" << chClusterId << "'" << std::endl;
      oss << std::string(nTab, ' ') << "AdcHits: { ";
      for (uint i = 0; i < _vAdcHits.size(); i++) oss << _vAdcHits[i] << ' ';
      oss << '}' << std::endl;
      return oss;
   }
   std::ostringstream &Print2(std::ostringstream &oss) const
   {
      PadHit::Print2(oss);
      uint       nClusterId  = getClusterId(_pCluster);
      char       chClusterId = _pCluster == NULL ? '.' : SZ_ALPHABET[nClusterId % strlen(SZ_ALPHABET)];
      const uint nTab        = 4;
      oss << std::string(nTab, ' ') << "TimeCluster:: pads: [" << _nTimeBinMin << "," << _nTimeBinMax << "]"
          << ", adcSpikes: { min[" << _iAdcMin << "], max[" << _iAdcMax << "] }"
          << ", cluster: " << nClusterId << " '" << chClusterId << "'" << std::endl;
      oss << std::string(nTab, ' ') << "AdcHits: { ";
      for (uint i = 0; i < _vAdcHits.size(); i++) oss << _vAdcHits[i] << ' ';
      oss << '}' << std::endl;
      return oss;
   }
};
struct ClusterStat {
   float _fAdcMin;
   float _fAdcMax;
   float _fSumMin;
   float _fSumMax;
   uint  _nMin;
   uint  _nMax;
   ClusterStat() : _fAdcMin(0), _fAdcMax(0), _fSumMin(0), _fSumMax(0), _nMin(UINT_MAX), _nMax(UINT_MAX) {}
   void Reset()
   {
      _fAdcMin = _fAdcMax = _fSumMin = _fSumMax = 0;
      _nMin = _nMax = UINT_MAX;
   }
   void Calc(const std::list<const PadCluster *> &rList)
   {
      Reset();
      for (std::list<const PadCluster *>::const_iterator i = rList.begin(); i != rList.end(); i++) {
         const PadCluster *pPadCluster = *i;
         float             fAdcMin     = pPadCluster->getAdcMin();
         float             fAdcMax     = pPadCluster->getAdcMax();
         if (_fAdcMin == 0 && _fAdcMax == 0) {
            _fAdcMin = fAdcMin;
            _fAdcMax = fAdcMax;
         } else {
            if (fAdcMin < _fAdcMin) _fAdcMin = fAdcMin;
            if (fAdcMax > _fAdcMax) _fAdcMax = fAdcMax;
         }
         float fAdcSum = pPadCluster->getAdcHitSum();
         if (fAdcSum > 0) {
            if (_fSumMin == 0 && _fSumMax == 0)
               _fSumMin = _fSumMax = fAdcSum;
            else if (fAdcSum < _fSumMin)
               _fSumMin = fAdcSum;
            else if (fAdcSum > _fSumMax)
               _fSumMax = fAdcSum;
         }
         uint nTimeBinMin = pPadCluster->getTimeBinMin();
         uint nTimeBinMax = pPadCluster->getTimeBinMax();
         if (_nMin == UINT_MAX && _nMax == UINT_MAX) {
            _nMin = nTimeBinMin;
            _nMax = nTimeBinMax;
         } else {
            if (nTimeBinMin < _nMin) _nMin = nTimeBinMin;
            if (nTimeBinMax > _nMax) _nMax = nTimeBinMax;
         }
      };
   }
};
const uint nADC_DISCRECITY = 40; // should be more than 9
class ClusterHit {
   uint  _nId;     // cluster id
   float _fAdcSum; // total ADC sum of all PadClusters
public:
   Gauss _PadPlane;  // normal distribution on pad plane
   Gauss _TimePlane; // normal distribution on time plane

   ClusterHit() : _nId(0), _fAdcSum(0) {}
   ClusterHit(uint nId) : _nId(nId), _fAdcSum(0) {}
   ClusterHit(const ClusterHit &r)
      : _nId(r._nId), _PadPlane(r._PadPlane), _TimePlane(r._TimePlane), _fAdcSum(r._fAdcSum)
   {
   }
   ~ClusterHit() {}
   inline void  setId(uint nId) { _nId = nId; }
   inline uint  getId() const { return _nId; }
   inline void  setAdcSum(float f) { _fAdcSum = f; }
   inline float getAdcSum() const { return _fAdcSum; }
   inline float getPad() const { return _PadPlane.getMean(); }
   inline float getTime() const { return _TimePlane.getMean(); }

   std::ostringstream &DrawPads(std::ostringstream &oss, uint nTimeBinMin, uint nTimeBinMax,
                                const std::list<const PadCluster *> &rlstPadClusters) const
   {
      int nTable = (nADC_DISCRECITY + 3) * 2 + nTimeBinMax - nTimeBinMin + 7;
      oss << std::string(nTable, '=') << std::endl;
      oss << " pads" << std::string(nTimeBinMax - nTimeBinMin + 1, ' ')
          << "| adcMax:" << std::string(nADC_DISCRECITY - 6, ' ') << "| " << _RED_("adcSum:")
          << std::string(nADC_DISCRECITY - 6, ' ') << '|' << std::endl;
      oss << std::string(nTable, '=') << std::endl;

      ClusterStat stat;
      stat.Calc(rlstPadClusters);
      for (std::list<const PadCluster *>::const_iterator i = rlstPadClusters.begin(); i != rlstPadClusters.end(); i++) {
         const PadCluster *pPadCluster = *i;
         if (i != rlstPadClusters.begin()) oss << std::endl;
         oss << std::setw(3) << pPadCluster->getPad() << ":";

         pPadCluster->Draw(oss, nTimeBinMin, nTimeBinMax);
         oss << " |";

         float fAdc = pPadCluster->getAdcHitMax().getAdc();
         long  nAdc = 0;
         if (stat._fAdcMax > stat._fAdcMin) {
            float f = (stat._fAdcMax - stat._fAdcMin) / nADC_DISCRECITY;
            nAdc    = lroundf((fAdc - stat._fAdcMin) / f);
         } else
            nAdc = nADC_DISCRECITY / 2;
         assert(nAdc <= nADC_DISCRECITY);

         if (nAdc) oss << std::string(nAdc, '.');
         oss << (fAdc == stat._fAdcMax ? '>' : '+');
         uint nSpace = nADC_DISCRECITY - nAdc + 1;
         if (nSpace) oss << std::string(nSpace, ' ');
         oss << '|';

         float fAdcSum = pPadCluster->getAdcHitSum();
         long  nAdcSum = 0;
         if (stat._fSumMax > stat._fSumMin) {
            float f = (stat._fSumMax - stat._fSumMin) / nADC_DISCRECITY;
            nAdcSum = lroundf((fAdcSum - stat._fSumMin) / f);
         } else
            nAdcSum = nADC_DISCRECITY / 2;
         assert(nAdcSum <= nADC_DISCRECITY);

         if (nAdcSum) oss << std::string(nAdcSum, '.');
         oss << (fAdcSum == stat._fSumMax ? '>' : '+');

         nSpace = nADC_DISCRECITY - nAdcSum + 1;
         if (nSpace) oss << std::string(nSpace, ' ');
         oss << '|';
      }
      oss << std::endl << std::string(nTable, '=') << "\npad ";
      _PadPlane.Print(oss);
      return oss;
   }
   std::ostringstream &DrawTimeBins(std::ostringstream &oss, uint nPadMin, uint nPadMax,
                                    const std::list<const PadCluster *> &rlstTimeClusters) const
   {
      int nTable = (nADC_DISCRECITY + 3) * 2 + nPadMax - nPadMin + 7;
      oss << std::string(nTable, '=') << std::endl;
      oss << " time" << std::string(nPadMax - nPadMin + 1, ' ') << "| adcMax:" << std::string(nADC_DISCRECITY - 6, ' ')
          << "| " << _RED_("adcSum:") << std::string(nADC_DISCRECITY - 6, ' ') << '|' << std::endl;
      oss << std::string(nTable, '=') << std::endl;
      ClusterStat stat;
      stat.Calc(rlstTimeClusters);
      for (std::list<const PadCluster *>::const_iterator i = rlstTimeClusters.begin(); i != rlstTimeClusters.end();
           i++) {
         const PadCluster *pTimeCluster = *i;
         if (i != rlstTimeClusters.begin()) oss << std::endl;
         oss << std::setw(3) << pTimeCluster->getPad() << ":";
         pTimeCluster->Draw2(oss, nPadMin, nPadMax);
         oss << " |";

         float fAdc = pTimeCluster->getAdcHitMax().getAdc();
         long  nAdc = 0;
         if (stat._fAdcMax > stat._fAdcMin) {
            float f = (stat._fAdcMax - stat._fAdcMin) / nADC_DISCRECITY;
            nAdc    = lroundf((fAdc - stat._fAdcMin) / f);
         } else
            nAdc = nADC_DISCRECITY / 2;
         assert(nAdc <= nADC_DISCRECITY);

         if (nAdc) oss << std::string(nAdc, '.');
         oss << (fAdc == stat._fAdcMax ? '>' : '+');
         uint nSpace = nADC_DISCRECITY - nAdc + 1;
         if (nSpace) oss << std::string(nSpace, ' ');
         oss << '|';

         float fAdcSum = pTimeCluster->getAdcHitSum();
         long  nAdcSum = 0;
         if (stat._fSumMax > stat._fSumMin) {
            float f = (stat._fSumMax - stat._fSumMin) / nADC_DISCRECITY;
            nAdcSum = lroundf((fAdcSum - stat._fSumMin) / f);
         } else
            nAdcSum = nADC_DISCRECITY / 2;
         assert(nAdcSum <= nADC_DISCRECITY);

         if (nAdcSum) oss << std::string(nAdcSum, '.');
         oss << (fAdcSum == stat._fSumMax ? '>' : '+');
         nSpace = nADC_DISCRECITY - nAdcSum + 1;
         if (nSpace) oss << std::string(nSpace, ' ');
         oss << '|';
      }
      oss << std::endl << std::string(nTable, '=') << "\ntime ";
      _TimePlane.Print(oss);
      return oss;
   }
   std::ostringstream &Draw(std::ostringstream &oss, uint nTimeBinMin, uint nTimeBinMax,
                            const std::list<const PadCluster *> &rlstPadClusters, uint nPadMin, uint nPadMax,
                            const std::list<const PadCluster *> &rlstTimeClusters) const
   {
      char ch = SZ_ALPHABET[getId() % strlen(SZ_ALPHABET)];
      oss << "Cluster:: id: " << _nId << " '" << ch << "', adcSum: " << _fAdcSum << ", PadMin: " << nPadMin
          << ", PadMax: " << nPadMax << ", TimeBinMin: " << nTimeBinMin << ", TimeBinMax: " << nTimeBinMax << std::endl;
      if (!rlstPadClusters.empty()) DrawPads(oss, nTimeBinMin, nTimeBinMax, rlstPadClusters);
      if (!rlstTimeClusters.empty()) DrawTimeBins(oss, nPadMin, nPadMax, rlstTimeClusters);
      return oss;
   }
   std::ostringstream &Print(std::ostringstream &oss) const
   {
      char ch = SZ_ALPHABET[getId() % strlen(SZ_ALPHABET)];
      oss << "Cluster:: id: " << _nId << " '" << ch << "', adcSum: " << _fAdcSum << std::endl;
      oss << "    pad ";
      _PadPlane.Print(oss);
      oss << "    time ";
      _TimePlane.Print(oss);
      return oss;
   }
};
class Cluster : public ClusterHit {
   friend class PadCluster;
   bool                          _bSplitted;   // the cluster was splitted (recalc stat info)
   uint                          _nTimeBinMin; // timebin min range
   uint                          _nTimeBinMax; // timebin max range
   uint                          _nPadMin;     // pad min range
   uint                          _nPadMax;     // pad max range
   std::list<const PadCluster *> _lstPadClusters;
   std::list<const PadCluster *> _lstTimeClusters;

   void calcGaussRes(const std::vector<XY<uint, float>> &rv, Gauss &r, uint nMin, uint nMax, float fAdcSum)
   {
      size_t nSize = rv.size();
      if (nSize > 2) {
         std::array<double, 4> res = ::calcGaussRes<uint, float>(rv);
         res[0] += nMin;
         // res[3] += nMin;
         r.Set(res);
#ifdef AUXHIT_DEBUG
         if (r.getDeviation() != 0) {
            r.setArea(::calcGaussArea(r.getMean(), r.getMax(), r.getDeviation(), nMin, nMax));
            r.setRsquared(
               calcRsquared<uint, float>(rv, r.getMean() - nMin, r.getMax(), r.getDeviation(), double(fAdcSum)));
         } else { // it cannot fit the Gauss
            r.setArea(0);
            r.setRsquared(0);
         }
#endif
#ifdef HIT_ADJUSTMENT
         double d    = r.getAdjustment() + nMin; // adjust unreasonable fit results to max AdcSum value!!!
         double dMin = nMin - 0.5;               // left shift in half of the min pad
         double dMax = nMax + 0.4999;            // right shift in half of the max pad
         if (r.getMean() < dMin)
            r.setMean(d - 0.5);
         else if (r.getMean() > dMax)
            r.setMean(d + 0.4999);
         else if (r.getDeviation() == 0)
            r.setMean(d);
#endif
      } else if (nSize == 2) {
         XY<uint, float> xy0   = rv[0];
         XY<uint, float> xy1   = rv[1];
         float           fMean = (xy0.x * xy0.y + xy1.x * xy1.y) / (xy0.y + xy1.y);
         r.setMean(fMean + nMin); // shift to real value
         r.setMax(xy0.y > xy1.y ? xy0.y : xy1.y);
         r.setDeviation(0);
         r.setArea(0);
         r.setRsquared(0);
      } else if (nSize == 1) {
         r.setMean(nMin);
         r.setMax(rv[0].y);
         r.setDeviation(0);
         r.setArea(0);
         r.setRsquared(0);
      } else
         assert(false);
   }

public:
   Cluster()
      : ClusterHit(0), _bSplitted(false), _nTimeBinMin(UINT_MAX), _nTimeBinMax(UINT_MAX), _nPadMin(UINT_MAX),
        _nPadMax(UINT_MAX)
   {
   }
   Cluster(int nPad, const AdcHit &r) : ClusterHit(), _bSplitted(false)
   {
      PadCluster *pPadCluster = new PadCluster(nPad, r);
      _nTimeBinMin            = pPadCluster->getTimeBinMin();
      _nTimeBinMax            = pPadCluster->getTimeBinMax();
      _nPadMin = _nPadMax = pPadCluster->getPad();
      _lstPadClusters.push_back(pPadCluster);
   }
   Cluster(uint nId, const PadCluster *pPadCluster) : ClusterHit(nId), _bSplitted(false)
   {
      _nTimeBinMin = pPadCluster->getTimeBinMin();
      _nTimeBinMax = pPadCluster->getTimeBinMax();
      _nPadMin = _nPadMax = pPadCluster->getPad();
      _lstPadClusters.push_back(pPadCluster);
   }
   virtual ~Cluster()
   {
      while (!_lstPadClusters.empty()) {
         delete _lstPadClusters.front();
         _lstPadClusters.pop_front();
      }
      while (!_lstTimeClusters.empty()) {
         delete _lstTimeClusters.front();
         _lstTimeClusters.pop_front();
      }
   }
   inline const std::list<const PadCluster *> &getPadClusters() const { return _lstPadClusters; } //!!!!!_nTrackId
   bool                                        isSplitted() const { return _bSplitted; }
   void                                        setSplitted() { _bSplitted = true; }
   void                                        resetSplitted() { _bSplitted = false; }

   void resetStat() { _nTimeBinMin = _nTimeBinMax = _nPadMin = _nPadMax = UINT_MAX; }
   void calcStat(const PadCluster *pPadCluster)
   {
      uint nTimeBinMin = pPadCluster->getTimeBinMin();
      uint nTimeBinMax = pPadCluster->getTimeBinMax();
      if (_nTimeBinMin == UINT_MAX && _nTimeBinMax == UINT_MAX) {
         _nTimeBinMin = nTimeBinMin;
         _nTimeBinMax = nTimeBinMax;
      } else {
         if (nTimeBinMin < _nTimeBinMin) _nTimeBinMin = nTimeBinMin;
         if (nTimeBinMax > _nTimeBinMax) _nTimeBinMax = nTimeBinMax;
      }
      uint nPad = pPadCluster->getPad();
      if (_nPadMin == UINT_MAX && _nPadMax == UINT_MAX)
         _nPadMin = _nPadMax = nPad;
      else if (nPad < _nPadMin)
         _nPadMin = nPad;
      else if (nPad > _nPadMax)
         _nPadMax = nPad;
   }
   void Add(const PadCluster *pPadCluster)
   {
      calcStat(pPadCluster);
      _lstPadClusters.push_back(pPadCluster);
   }
   // void Add( const Cluster* pCluster ) { std::copy( pCluster->_lstPadClusters.begin(),
   // pCluster->_lstPadClusters.end(), _lstPadClusters.end() ); }
   void Add(const Cluster *pCluster)
   {
      std::for_each(pCluster->_lstPadClusters.begin(), pCluster->_lstPadClusters.end(),
                    [this](const PadCluster *pPadCluster) {
                       ((PadCluster *)pPadCluster)->setCluster(this);
                       this->Add(pPadCluster);
                    });
   }
   inline void Clear()
   {
      _lstPadClusters.clear();
      _lstTimeClusters.clear();
   }
   inline uint getTimeBinMin() const { return _nTimeBinMin; }
   inline uint getTimeBinMax() const { return _nTimeBinMax; }
   inline uint getPadMin() const { return _nPadMin; }
   inline uint getPadMax() const { return _nPadMax; }
   void        renewStat()
   {
      resetStat();
      for (std::list<const PadCluster *>::const_iterator i = _lstPadClusters.begin(); i != _lstPadClusters.end(); i++) {
         const PadCluster *pPadCluster = *i;
         calcStat(pPadCluster);
         ((PadCluster *)pPadCluster)->setCluster(this);
      }
   }
   Cluster *Split()
   { // split current Cluster on new ones depends on AdcHitSum values
      uint                                          nMin = 0;
      uint                                          nMax = 0;
      std::list<const PadCluster *>::const_iterator iAdcMin, iAdcMax;
      std::list<const PadCluster *>::iterator       i = _lstPadClusters.begin();
      iAdcMin = iAdcMax = i++;
      for (uint n = 1; i != _lstPadClusters.end(); i++, n++) {
         // float fAdc = (*i)->getAdc();
         float fAdc = (*i)->getAdcHitSum();
         if (nMax >= nMin) {
            if (fAdc > (*iAdcMax)->getAdcHitSum()) {
               iAdcMax = i;
               nMax    = n;
            } else {
               iAdcMin = i;
               nMin    = n;
            }
         } else if (fAdc < (*iAdcMin)->getAdcHitSum()) {
            iAdcMin = i;
            nMin    = n;
         } else {
            const PadCluster *pPadCluster = new PadCluster(*(*iAdcMin));
            Cluster          *pCluster    = new Cluster();
            pCluster->_lstPadClusters.splice(pCluster->_lstPadClusters.end(), _lstPadClusters, _lstPadClusters.begin(),
                                             i);
            pCluster->renewStat();
            _lstPadClusters.push_front(pPadCluster);
            setSplitted();
            return pCluster;
         }
      }
      return NULL;
   }
   void calcGaussRes()
   {
      float                        fClusterAdcSum = 0;
      std::vector<XY<uint, float>> vXpadYadc;
      for (std::list<const PadCluster *>::const_iterator i = _lstPadClusters.begin(); i != _lstPadClusters.end(); i++) {
         const PadCluster *pPadCluster = *i;
         // float fAdc = pPadCluster->getAdc();
         float fAdcSum = pPadCluster->getAdcHitSum();
         fClusterAdcSum += fAdcSum;

         uint                                   nPad = pPadCluster->getPad() - _nPadMin;
         std::vector<XY<uint, float>>::iterator iPad =
            std::find_if(vXpadYadc.begin(), vXpadYadc.end(), [nPad](XY<uint, float> &r) { return r.x == nPad; });
         if (iPad != vXpadYadc.end())
            (*iPad).y += fAdcSum;
         else
            vXpadYadc.push_back({.x = nPad, .y = fAdcSum});
         size_t nSize = pPadCluster->getSize();
         for (uint j = 0; j < nSize; j++) {
            AdcHit adcHit = pPadCluster->           operator[](j);
            AdcHit                                  timeHit(pPadCluster->getPad(), adcHit.getAdc());
            uint                                    nTimeBin = adcHit.getTimeBin();
            std::list<const PadCluster *>::iterator ij =
               std::find_if(_lstTimeClusters.begin(), _lstTimeClusters.end(),
                            [nTimeBin](const PadCluster *p) { return p->getPad() == nTimeBin; });
            if (ij != _lstTimeClusters.end())
               ((PadCluster *)*ij)->Add(timeHit);
            else
               _lstTimeClusters.push_back(new PadCluster(nTimeBin, timeHit, this));
         }
      }
      calcGaussRes(vXpadYadc, ClusterHit::_PadPlane, _nPadMin, _nPadMax, fClusterAdcSum);

      ClusterHit::setAdcSum(fClusterAdcSum);
      _lstTimeClusters.sort(
         [](const PadCluster *&rp1, const PadCluster *&rp2) { return rp1->getPad() < rp2->getPad(); });

      std::vector<XY<uint, float>> vXtimeYadc;
      for (std::list<const PadCluster *>::const_iterator i = _lstTimeClusters.begin(); i != _lstTimeClusters.end();
           i++) {
         const PadCluster                      *pTimeCluster = *i;
         float                                  fAdcSum      = pTimeCluster->getAdcHitSum();
         uint                                   nTime        = pTimeCluster->getPad() - _nTimeBinMin;
         std::vector<XY<uint, float>>::iterator iTime =
            std::find_if(vXtimeYadc.begin(), vXtimeYadc.end(), [nTime](XY<uint, float> &r) { return r.x == nTime; });
         if (iTime != vXtimeYadc.end())
            (*iTime).y += fAdcSum;
         else
            vXtimeYadc.push_back({.x = nTime, .y = fAdcSum});
      }
      calcGaussRes(vXtimeYadc, ClusterHit::_TimePlane, _nTimeBinMin, _nTimeBinMax, fClusterAdcSum);
   }
   std::ostringstream &Draw(std::ostringstream &oss) const
   {
      ClusterHit::Draw(oss, _nTimeBinMin, _nTimeBinMax, _lstPadClusters, _nPadMin, _nPadMax, _lstTimeClusters);
      return oss;
   }
   std::ostringstream &Print(std::ostringstream &oss, const std::list<const PadCluster *> &r) const
   {
      const uint nTAB = 4;
      if (!r.empty()) {
         oss << std::string(nTAB, ' ') << "======== ";
         for (std::list<const PadCluster *>::const_iterator i = r.begin(); i != r.end(); i++) {
            if (i != r.begin()) oss << std::string(nTAB, ' ') << "-------- ";
            if (&r == &_lstPadClusters)
               (*i)->Print(oss);
            else
               (*i)->Print2(oss);
         }
      }
      return oss;
   }
   std::ostringstream &Print(std::ostringstream &oss) const
   {
      ClusterHit::Print(oss);
      const uint nTAB = 4;
      oss << std::string(nTAB, ' ') << "PadMin: " << _nPadMin << ", PadMax: " << _nPadMax
          << ", TimeBinMin: " << _nTimeBinMin << ", TimeBinMax: " << _nTimeBinMax << std::endl;
      Print(oss, _lstPadClusters);
      Print(oss, _lstTimeClusters);
      oss << std::endl;
      return oss;
   }
};
uint PadCluster::getClusterId(const Cluster *pCluster) const
{
   return pCluster == NULL ? 0 : ((const ClusterHit *)pCluster)->getId();
}
float PadCluster::getClusterTime(const Cluster *pCluster) const
{
   return pCluster == NULL ? 0 : ((const ClusterHit *)pCluster)->getTime();
}
float PadCluster::getClusterPad(const Cluster *pCluster) const
{
   return pCluster == NULL ? 0 : ((const ClusterHit *)pCluster)->getPad();
}

inline const char *szDECADE            = "0123456789";
const uint         nTIMEBIN_SCALE_SIZE = 30;
inline void        DrawTimeBinTopScale(std::ostringstream &oss, uint nSize)
{
   oss << std::setw(4) << ":";
   for (uint i = 0; i < nSize; i++) oss << std::left << std::setw(10) << i * 10;
   oss << std::endl;
   oss << std::internal << std::setw(4) << ":";
   for (uint i = 0; i < nSize; i++) oss << szDECADE;
   oss << std::internal << std::endl;
}
inline void DrawTimeBinBottomScale(std::ostringstream &oss, uint nSize)
{
   oss << std::setw(4) << ":";
   for (uint i = 0; i < nSize; i++) oss << szDECADE;
   oss << std::endl << std::internal << std::setw(4) << ":";
   for (uint i = 0; i < nSize; i++) oss << std::left << std::setw(10) << i * 10;
   oss << std::internal << std::endl;
}
inline std::ostringstream &DrawPadClusters(std::ostringstream                             &oss,
                                           std::vector<std::list<PadCluster *>::iterator> &rviPadClusters)
{
   oss << "RowClusters::" << std::endl;
   // DrawTimeBinTopScale(oss, nTIMEBIN_SCALE_SIZE);
   assert(rviPadClusters.size() > 1);
   size_t nPadSize = rviPadClusters.size() - 1;
   for (uint i = 0; i < nPadSize; i++) {
      std::list<PadCluster *>::iterator iPadCluster     = rviPadClusters[i];
      std::list<PadCluster *>::iterator iPadClusterLast = rviPadClusters[i + 1];
      oss << std::setw(3) << (*iPadCluster)->getPad() << ":";
      uint nPos = 0;
      for (iPadCluster; iPadCluster != iPadClusterLast; iPadCluster++) {
         PadCluster *pPadCluster = *iPadCluster;
         uint        nMin        = pPadCluster->getTimeBinMin();
         if (nMin < nPos) { // shift back to the 'n' positions (for combined PadClusters!!!)
            long int nShift = nPos - nMin;
            oss.seekp(oss.tellp() - nShift);
            nPos -= nShift;
         }
         uint nGap = nMin - nPos;
         if (nGap) {
            oss << std::string(nGap, '.');
            nPos += nGap;
         }
         char ch =
            pPadCluster->getCluster() ? SZ_ALPHABET[pPadCluster->getCluster()->getId() % strlen(SZ_ALPHABET)] : '+';
         pPadCluster->Draw(oss, ch);
         nPos += pPadCluster->getSize();
      }
      oss << std::endl;
   }
   // DrawTimeBinBottomScale(oss, nTIMEBIN_SCALE_SIZE);
   return oss;
}
inline std::ostringstream &PrintPadClusters(std::ostringstream                             &oss,
                                            std::vector<std::list<PadCluster *>::iterator> &rviPadClusters)
{
   oss << "RowClusters::" << std::endl;
   assert(rviPadClusters.size() > 1);
   size_t nSize1 = rviPadClusters.size() - 1;
   for (uint i = 0; i < nSize1; i++) {
      oss << "    ---------------- id: " << i << std::endl;
      std::list<PadCluster *>::iterator iPadCluster     = rviPadClusters[i];
      std::list<PadCluster *>::iterator iPadClusterLast = rviPadClusters[i + 1];
      for (iPadCluster; iPadCluster != iPadClusterLast; iPadCluster++) {
         oss << "\t";
         (*iPadCluster)->Print(oss);
      }
   }
   oss << "    ----------------" << std::endl;
   return oss;
}
class RowClusters {
   uint                 _nRow;       // sector row number
   uint                 _nClusterId; // current cluster id for internal usage
   std::list<Cluster *> _lstClusters;
   void Joint(std::vector<std::list<PadCluster *>::iterator> &rviPadClusters, uint nPadCluster, PadCluster *pPadCluster)
   {
      Cluster                          *pCluster        = NULL;
      uint                              nTimeBinMin     = 0;
      uint                              nTimeBinMax     = 0;
      uint                              nPadMin         = 0;
      uint                              nPadMax         = 0;
      std::list<PadCluster *>::iterator iPadCluster1st  = rviPadClusters[nPadCluster - 1];
      std::list<PadCluster *>::iterator iPadClusterLast = rviPadClusters[nPadCluster];
      for (std::list<PadCluster *>::iterator iPadCluster = iPadCluster1st; iPadCluster != iPadClusterLast;
           iPadCluster++) {
         PadCluster *pPadCluster0 = (*iPadCluster);
         if (pPadCluster->isBefore(pPadCluster0))
            break;
         else if (pPadCluster->isAfter(pPadCluster0))
            continue;
         else if (pPadCluster->isBelongTo(pPadCluster0)) {
            pPadCluster->setPadTimeOverlapping(pPadCluster0);
            pCluster = pPadCluster->getCluster();
            if (pCluster == NULL) {
               pCluster = pPadCluster0->getCluster();
               assert(pCluster != NULL);

               nTimeBinMin = pCluster->getTimeBinMin();
               nTimeBinMax = pCluster->getTimeBinMax();
               nPadMin     = pCluster->getPadMin();
               nPadMax     = pCluster->getPadMax();

               pPadCluster->setCluster(pCluster);
               pCluster->Add(pPadCluster);
            } else {
               Cluster *pCluster0 = pPadCluster0->getCluster();
               assert(pCluster0);
               if (pCluster0->getTimeBinMax() >= nTimeBinMin && pCluster0->getTimeBinMin() <= nTimeBinMax &&
                   pCluster0->getPadMax() ==
                      pCluster0->getPadMin()) { // if it is truth then we should add whole PadClusters list to the
                                                // Cluster and remove Cluster0
                  pCluster->Add(pCluster0);
                  pCluster0->Clear();
                  _lstClusters.remove_if([pCluster0](Cluster *p) { return p == pCluster0; });
               }
            }
         }
      }
      if (pCluster == NULL) {
         pCluster = new Cluster(_nClusterId++, pPadCluster);
         pPadCluster->setCluster(pCluster);
         _lstClusters.push_back(pCluster);
      }
   }

public:
   RowClusters() : _nRow(0), _nClusterId(0) {}
   RowClusters(uint nRow) : _nRow(nRow), _nClusterId(0) {}
   virtual ~RowClusters()
   {
      while (!_lstClusters.empty()) {
         delete _lstClusters.front();
         _lstClusters.pop_front();
      }
   }
   inline uint                        getRow() const { return _nRow; }
   inline const std::list<Cluster *> &getClusters() const { return _lstClusters; }
   inline bool                        isEmpty() const { return _lstClusters.empty(); }
   void                               Joint(std::vector<std::list<PadCluster *>::iterator> &rviPadClusters)
   {
      assert(rviPadClusters.size() > 1);
      size_t nPadSize = rviPadClusters.size() - 1;
      for (uint i = 0; i < nPadSize; i++) {
         std::list<PadCluster *>::iterator iPadCluster1st  = rviPadClusters[i];
         std::list<PadCluster *>::iterator iPadClusterLast = rviPadClusters[i + 1];
         for (std::list<PadCluster *>::iterator iPadCluster = iPadCluster1st; iPadCluster != iPadClusterLast;
              iPadCluster++) {
            PadCluster *pPadCluster = (*iPadCluster);
            if (i)
               Joint(rviPadClusters, i, pPadCluster);
            else { // init zero pad clusters
               Cluster *pCluster = new Cluster(_nClusterId++, pPadCluster);
               pPadCluster->setCluster(pCluster);
               _lstClusters.push_back(pCluster);
            }
         }
      }
   }
   void Split()
   {
      for (std::list<Cluster *>::iterator i = _lstClusters.begin(); i != _lstClusters.end(); ++i) {
         Cluster *pCluster = *i;
         while (Cluster *pNewCluster = pCluster->Split()) {
            pNewCluster->setId(_nClusterId++);
            pNewCluster->calcGaussRes();
            _lstClusters.insert(i, pNewCluster);
         }
         if (pCluster->isSplitted()) pCluster->renewStat();
         pCluster->calcGaussRes();
      }
   }
   std::ostringstream &Draw(std::ostringstream &oss) const
   {
      oss << "RowClusters:: row: " << _nRow << std::endl;
      for (std::list<Cluster *>::const_iterator i = _lstClusters.begin(); i != _lstClusters.end(); ++i) {
         (*i)->Draw(oss);
         oss << std::endl;
      }
      return oss;
   }
   std::ostringstream &Print(std::ostringstream &oss, std::vector<std::list<PadCluster *>::iterator> &rviPadClusters)
   {
      oss << "RowClusters:: row: " << getRow() << std::endl;
      if (rviPadClusters.size() > 1) {
         size_t nSize1 = rviPadClusters.size() - 1;
         for (uint i = 0; i < nSize1; i++) {
            std::list<PadCluster *>::iterator iPadCluster     = rviPadClusters[i];
            std::list<PadCluster *>::iterator iPadClusterLast = rviPadClusters[i + 1];
            oss << "    ---------------- id: " << i << std::endl;
            for (iPadCluster; iPadCluster != iPadClusterLast; iPadCluster++) {
               oss << "\t";
               (*iPadCluster)->Print(oss);
            }
         }
      }
      oss << "    ----------------" << std::endl;
      return oss;
   }
   std::ostringstream &Print(std::ostringstream &oss) const
   {
      oss << "RowClusters:: row: " << getRow() << std::endl;
      for (std::list<Cluster *>::const_iterator i = _lstClusters.begin(); i != _lstClusters.end(); ++i)
         (*i)->Print(oss);
      return oss;
   }
};
class SectorClusters {
   int                      _nSector;
   std::list<RowClusters *> _lstRowClusters;

public:
   SectorClusters() : _nSector(0) {}
   SectorClusters(int nSector) : _nSector(nSector) {}
   virtual ~SectorClusters()
   {
      while (!_lstRowClusters.empty()) {
         delete _lstRowClusters.front();
         _lstRowClusters.pop_front();
      }
   }
   inline int                             getSector() const { return _nSector; }
   inline const std::list<RowClusters *> &getRowClusters() const { return _lstRowClusters; }
   inline void                            Add(RowClusters *pRowClusters) { _lstRowClusters.push_back(pRowClusters); }
   std::ostringstream                    &Print(std::ostringstream &oss) const
   {
      oss << "SectorClusters:: sector: " << getSector() << std::endl;
      for (std::list<RowClusters *>::const_iterator i = _lstRowClusters.begin(); i != _lstRowClusters.end(); i++)
         (*i)->Print(oss);
      return oss;
   }
};
class EventClusters {
   int                         _nEvent;
   std::list<SectorClusters *> _lstSectorClusters;

public:
   EventClusters() : _nEvent(0) {}
   EventClusters(int nEvent) : _nEvent(nEvent) {}
   virtual ~EventClusters()
   {
      while (!_lstSectorClusters.empty()) {
         delete _lstSectorClusters.front();
         _lstSectorClusters.pop_front();
      }
   }
   inline int                                getEvent() const { return _nEvent; }
   inline const std::list<SectorClusters *> &getSectorClusters() const { return _lstSectorClusters; }
   inline SectorClusters *getFront() { return _lstSectorClusters.empty() ? NULL : _lstSectorClusters.front(); }
   inline SectorClusters *getBack() { return _lstSectorClusters.empty() ? NULL : _lstSectorClusters.back(); }
   inline void            Add(SectorClusters *pSectorClusters) { _lstSectorClusters.push_back(pSectorClusters); }
   std::ostringstream    &Print(std::ostringstream &oss) const
   {
      oss << "EventClusters:: event: " << getEvent() << std::endl;
      for (std::list<SectorClusters *>::const_iterator i = _lstSectorClusters.begin(); i != _lstSectorClusters.end();
           i++)
         (*i)->Print(oss);
      return oss;
   }
};
} // namespace tpcClustering
#endif