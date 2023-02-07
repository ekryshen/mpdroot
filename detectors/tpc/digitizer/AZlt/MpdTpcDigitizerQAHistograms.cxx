//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class TpcMWPCQAHistograms
//      see TpcMWPCQAHistograms.h for details
//
// Environment:
//      Software developed for the MPC @ NICA
//
// Author List:
//      Sergey Merts
//
//
//-----------------------------------------------------------

// Panda Headers ----------------------

// This Class' Header ------------------
#include "MpdTpcDigitizerQAHistograms.h"

MpdTpcDigitizerQAHistograms::MpdTpcDigitizerQAHistograms()
   : _hRZ_global(0), _hYZ_local(0), _hDiffuseXY(0), _hDistortXY(0), _hX_global(0), _hY_global(0), _hZ_global(0),
     _hX_local(0), _hY_local(0), _hZ_local(0), _hXY_local(0), _hXY_global(0), _hSect_dig(0), _hADC_dig(0), _hX_dig(0),
     _hY_dig(0), _hZ_dig(0), _h3D_dig(0), _h3D_el(0), _hXT_dig_1(0), _hXT_dig_5(0), _hXT_dig_10(0), _hXT_dig_20(0),
     _hXT_dig_40(0), _hXT_dig_60(0)
{
}

MpdTpcDigitizerQAHistograms::MpdTpcDigitizerQAHistograms(const std::string &suffix)
   : _hRZ_global(0), _hYZ_local(0), _hDiffuseXY(0), _hDistortXY(0), _hX_global(0), _hY_global(0), _hZ_global(0),
     _hX_local(0), _hY_local(0), _hZ_local(0), _hXY_local(0), _hXY_global(0), _hSect_dig(0), _hADC_dig(0), _hX_dig(0),
     _hY_dig(0), _hZ_dig(0), _h3D_dig(0), _h3D_el(0), _hXT_dig_1(0), _hXT_dig_5(0), _hXT_dig_10(0), _hXT_dig_20(0),
     _hXT_dig_40(0), _hXT_dig_60(0)
{
}

MpdTpcDigitizerQAHistograms::~MpdTpcDigitizerQAHistograms()
{

   delete _hRZ_global;
   delete _hYZ_local;
   delete _hX_dig;
   delete _hY_dig;
   delete _hZ_dig;
   delete _hX_local;
   delete _hY_local;
   delete _hZ_local;
   delete _hXY_local;
   delete _hY_global;
   delete _hX_global;
   delete _hZ_global;
   delete _hXY_global;
   delete _hSect_dig;
   delete _hXY_dig;
   delete _h3D_dig;
   delete _hADC_dig;
   delete _hDiffuseXY;
   delete _hDistortXY;
   delete _h3D_el;
   delete _hXT_dig_1;
   delete _hXT_dig_5;
   delete _hXT_dig_10;
   delete _hXT_dig_20;
   delete _hXT_dig_40;
   delete _hXT_dig_60;
}

void MpdTpcDigitizerQAHistograms::Initialize()
{

   UInt_t     nTimeBin = 512;
   TpcSector *sector   = new TpcSector();
   UInt_t     nSectors = sector->GetNSectors();
   UInt_t     nRows    = sector->GetNumRows();
   //    UInt_t nInRows = sector->GetNumInnerRows();
   UInt_t  nPads = sector->GetNumPadsInRow(nRows - 1) * 2;
   Float_t zMax  = sector->GetLength();
   Float_t maxX  = sector->GetMaxX();
   Float_t maxY  = sector->GetSectHeight() + 1.0;
   Float_t Rmax  = sector->GetRmin() + sector->GetSectHeight();

   _h3D_el = CreateHistogram3("3D distribution of electrons in TPC", "X (cm)", "Y (cm)", "Z (cm)", 200, -150, 150, 200,
                              -150, 150, 200, -zMax, zMax);

   _hX_local  = CreateHistogram1("local X distribution of electrons in TPC", _suffix, "X, cm", 100, -maxX, maxX);
   _hY_local  = CreateHistogram1("local Y distribution of electrons in TPC", _suffix, "Y, cm", 100, 0.0, maxY);
   _hZ_local  = CreateHistogram1("local Z distribution of electrons in TPC", _suffix, "Z, cm", 300, 0.0, zMax);
   _hXY_local = CreateHistogram2("local XY distribution of electrons in TPC", _suffix, "X, cm", "Y, cm", 500, -maxX,
                                 maxX, 500, 0.0, maxY);
   _hYZ_local = CreateHistogram2("local ZY distribution of electrons in TPC", _suffix, "Z, cm", "Y, cm", 500, 0.0, zMax,
                                 500, 0.0, maxY);

   _hX_global  = CreateHistogram1("global X distribution of electrons in TPC", _suffix, "X, cm", 200, -150, 150);
   _hY_global  = CreateHistogram1("global Y distribution of electrons in TPC", _suffix, "Y, cm", 200, -150, 150);
   _hZ_global  = CreateHistogram1("global Z distribution of electrons in TPC", _suffix, "Z, cm", 200, -zMax, zMax);
   _hXY_global = CreateHistogram2("global XY distribution of electrons in TPC", _suffix, "X, cm", "Y, cm", 800, -150,
                                  150, 800, -150, 150);
   _hRZ_global = CreateHistogram2("global ZR distribution of electrons in TPC", _suffix, "Z, cm", "R, cm", 800, -zMax,
                                  zMax, 800, -Rmax, Rmax);

   _hSect_dig =
      CreateHistogram1("Sector distribution of TPC digits", _suffix, "number of sector", nSectors, 0, nSectors);
   _hADC_dig = CreateHistogram1("Charge distribution of TPC digits", _suffix, "charge", 200, 900, 2000000);
   _hX_dig   = CreateHistogram1("X distribution of TPC digits", _suffix, "X, pads", nPads, 0, nPads);
   _hY_dig   = CreateHistogram1("Y distribution of TPC digits", _suffix, "Y, pad rows", nRows, 0, nRows);
   _hZ_dig   = CreateHistogram1("Z distribution of TPC digits", _suffix, "Z, time bins", nTimeBin, 0, nTimeBin);
   _hXY_dig  = CreateHistogram2("XY distribution of TPC digits", _suffix, "X, pads", "Y, pad rows", nPads, 0, nPads,
                                nRows, 0, nRows);
   _h3D_dig  = CreateHistogram3("XYZ distribution of TPC digits", _suffix, "X, pads", "Y, pad rows", "Z, time bins",
                                nPads, 0, nPads, nRows, 0, nRows, nTimeBin, 0, nTimeBin);
   _hDiffuseXY =
      CreateHistogram2("XY diffuse of electrons", _suffix, "X (cm)", "Y (cm)", 400, -2.0, 2.0, 400, -2.0, 2.0);
   _hDistortXY =
      CreateHistogram2("XY distortion of electrons", _suffix, "X (cm)", "Y (cm)", 400, -2.0, 2.0, 400, -2.0, 2.0);

   _hXT_dig_1  = CreateHistogram2("XT distribution of TPC digits on PadRow #1", _suffix, "X, pads", "T, time bins",
                                  nPads, 0, nPads, nTimeBin, 0, nTimeBin);
   _hXT_dig_5  = CreateHistogram2("XT distribution of TPC digits on PadRow #5", _suffix, "X, pads", "T, time bins",
                                  nPads, 0, nPads, nTimeBin, 0, nTimeBin);
   _hXT_dig_10 = CreateHistogram2("XT distribution of TPC digits on PadRow #10", _suffix, "X, pads", "T, time bins",
                                  nPads, 0, nPads, nTimeBin, 0, nTimeBin);
   _hXT_dig_20 = CreateHistogram2("XT distribution of TPC digits on PadRow #20", _suffix, "X, pads", "T, time bins",
                                  nPads, 0, nPads, nTimeBin, 0, nTimeBin);
   _hXT_dig_40 = CreateHistogram2("XT distribution of TPC digits on PadRow #40", _suffix, "X, pads", "T, time bins",
                                  nPads, 0, nPads, nTimeBin, 0, nTimeBin);
   _hXT_dig_60 = CreateHistogram2("XT distribution of TPC digits on PadRow #60", _suffix, "X, pads", "T, time bins",
                                  nPads, 0, nPads, nTimeBin, 0, nTimeBin);
}

void MpdTpcDigitizerQAHistograms::Write()
{
   _hRZ_global->Write(0, kOverwrite);
   _hYZ_local->Write(0, kOverwrite);
   _hSect_dig->Write(0, kOverwrite);
   _hX_dig->Write(0, kOverwrite);
   _hY_dig->Write(0, kOverwrite);
   _hZ_dig->Write(0, kOverwrite);
   _hX_dig->Write(0, kOverwrite);
   _hXY_dig->Write(0, kOverwrite);
   _h3D_dig->Write(0, kOverwrite);
   _hDiffuseXY->Write(0, kOverwrite);
   _hDistortXY->Write(0, kOverwrite);
   _hX_local->Write(0, kOverwrite);
   _hY_local->Write(0, kOverwrite);
   _hXY_global->Write(0, kOverwrite);
   _hZ_local->Write(0, kOverwrite);
   _hXY_local->Write(0, kOverwrite);
   _hX_global->Write(0, kOverwrite);
   _hY_global->Write(0, kOverwrite);
   _hZ_global->Write(0, kOverwrite);
   _hADC_dig->Write(0, kOverwrite);
   _h3D_el->Write(0, kOverwrite);
   _hXT_dig_1->Write(0, kOverwrite);
   _hXT_dig_5->Write(0, kOverwrite);
   _hXT_dig_10->Write(0, kOverwrite);
   _hXT_dig_20->Write(0, kOverwrite);
   _hXT_dig_40->Write(0, kOverwrite);
   _hXT_dig_60->Write(0, kOverwrite);
}

ClassImp(MpdTpcDigitizerQAHistograms);
