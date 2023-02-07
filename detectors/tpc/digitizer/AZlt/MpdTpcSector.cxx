#include "MpdTpcSector.h"
#include <TMath.h>
#include <TF1.h>
#include <cassert>
#include <algorithm>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>

using namespace std;

ClassImp(TpcSector);

TpcSector::TpcSector()
{

   _tpcLength      = 170;
   _outerPadHeight = 1.8;
   _outerPadWidth  = 0.5;
   _innerPadHeight = 1.2;
   _innerPadWidth  = 0.5;
   _nSectors       = 24;
   _nTimeBins      = 512;

   _rMin = 40.3;

   _maxX            = 30.73; // 32.08 - 1.35
   _minX            = 9.295; // 10.645 - 1.35
   _sectInnerHeight = 32.4;
   _sectOuterHeight = 46.8;
   _sectHeight      = 79.2; // 80.0 - 0.4 - 0.4
   Float_t tang     = (_maxX - _minX) / _sectHeight;
   _midX            = _minX + tang * _sectInnerHeight;

   Int_t nOuterPads_min = (Int_t)(_midX / _outerPadWidth);
   Int_t nInnerPads_min = (Int_t)(_minX / _innerPadWidth);

   _nInnerRows = 27; //(UInt_t) (_sectInnerHeight / _innerPadHeight);
   _nOuterRows = 26; //(UInt_t) (_sectOuterHeight / _outerPadHeight);
   _nRows      = _nInnerRows + _nOuterRows;

   _numPadsInRow      = new Int_t[_nRows];
   _numInnerPadsInRow = new Int_t[_nInnerRows];
   _numOuterPadsInRow = new Int_t[_nOuterRows];

   for (UInt_t k = 0; k < _nInnerRows; ++k) {
      Float_t w             = k * _innerPadHeight * tang;
      Int_t   nExtraPads    = (UInt_t)(w / _innerPadWidth);
      _numPadsInRow[k]      = nInnerPads_min + nExtraPads;
      _numInnerPadsInRow[k] = nInnerPads_min + nExtraPads;
   }
   for (UInt_t k = 0; k < _nOuterRows; ++k) {
      Float_t w                      = k * _outerPadHeight * tang;
      Int_t   nExtraPads             = (UInt_t)(w / _outerPadWidth);
      _numPadsInRow[_nInnerRows + k] = nOuterPads_min + nExtraPads;
      _numOuterPadsInRow[k]          = nOuterPads_min + nExtraPads;
   }
}

TpcSector::~TpcSector()
{
   delete[] _numInnerPadsInRow;
   delete[] _numOuterPadsInRow;
}
