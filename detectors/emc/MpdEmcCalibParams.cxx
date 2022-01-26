#include "MpdEmcCalibParams.h"
#include "MpdEmcGeoUtils.h"

#include "FairLogger.h"

#include <TH2.h>

#include <iostream>


bool MpdEmcCalibParams::SetGain(TH2* h)
{
  const int MAXX = 300,
            MAXZ = 128;
  if (!h) {
    LOG(ERROR) << "no input histogam";
    return false;
  }

  if (h->GetNbinsX() != MAXX || h->GetNbinsY() != MAXZ) {
    LOG(ERROR) << "Wrong dimentions of input histogram:" << h->GetNbinsX() << "," << h->GetNbinsY() << " instead of " << MAXX << "," << MAXZ;
    return false;
  }

  auto geo = MpdEmcGeoUtils::GetInstance();
  if (!geo) {
    LOG(ERROR) << "Geometry needs to be initialized";
    return false;
  }

  for(int tower=0; tower<NCHANNELS; tower++){
    int iphi,iz  ;
    geo->DetIdToGlobalIphiIz(tower,iphi,iz) ;
    fGainCalib[tower] = h->GetBinContent(iphi+1, iz+1);
  }
  return true;
}


bool MpdEmcCalibParams::SetTimeCalib(TH2* h)
{
  const int MAXX = 300,
            MAXZ = 128;
  if (!h) {
    LOG(ERROR) << "no input histogam";
    return false;
  }

  if (h->GetNbinsX() != MAXX || h->GetNbinsY() != MAXZ) {
    LOG(ERROR) << "Wrong dimentions of input histogram:" << h->GetNbinsX() << "," << h->GetNbinsY() << " instead of " << MAXX << "," << MAXZ;
    return false;
  }

  auto geo = MpdEmcGeoUtils::GetInstance();
  if (!geo) {
    LOG(ERROR) << "Geometry needs to be initialized";
    return false;
  }

  for(int tower=0; tower<NCHANNELS; tower++){
    int iphi,iz  ;
    geo->DetIdToGlobalIphiIz(tower,iphi,iz) ;
    fTimeCalib[tower] = h->GetBinContent(iphi+1, iz+1);
  }

  return true;
}

