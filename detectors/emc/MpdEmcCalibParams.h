//--------------------------------------------------------------------
//
// Description:
//      MPD EMC calibration parameters
//
//
// Author List:
//      D.Peresunko
//
//--------------------------------------------------------------------

#ifndef MPDEMCCALIBPARAMS_H
#define MPDEMCCALIBPARAMS_H

#include <array>
#include "TObject.h"

class TH2;


class MpdEmcCalibParams: public TObject
{
 public:
  /// \brief Constructor
  MpdEmcCalibParams() = default;

  /// \brief Destructor
  ~MpdEmcCalibParams() = default;

  /// \brief Get energy calibration coefficients
  /// \param cellID Absolute ID of cell
  /// \return  energy calibration coefficient of the cell
  float GetGain(unsigned int cellID) { return fGainCalib.at(cellID); }

  /// \brief Set High Gain energy calibration coefficient
  /// \param cellID Absolute ID of cell
  /// \param c is the calibration coefficient
  void SetGain(unsigned int cellID, float c) { fGainCalib[cellID] = c; }

  /// \brief Set energy calibration coefficients in the form of 2D histogram
  /// \param 2D(phi,-z...z) histogram with calibration coefficients
  /// \return Is successful
  bool SetGain(TH2* h);

  /// \brief Get High Gain time calibration coefficients
  /// \param cellID Absolute ID of cell
  /// \return high gain time calibration coefficient of the cell
  float GetTimeCalib(unsigned int cellID) { return fTimeCalib.at(cellID); }

  /// \brief Set High Gain time calibration coefficient
  /// \param cellID Absolute ID of cell
  /// \param t is the calibration coefficient
  void SetTimeCalib(unsigned int cellID, float t) { fTimeCalib[cellID] = t; }

  /// \brief Set time calibration coefficients in the form of 2D histogram
  /// \param 2D(phi,-z..z) histogram with calibration coefficients
  /// \return Is successful
  bool SetTimeCalib(TH2* h);

 private:
  static constexpr int NCHANNELS = 38400;    ///< Number of channels 
  std::array<float, NCHANNELS> fGainCalib;   ///< Container for the gain calibration coefficients
  std::array<float, NCHANNELS> fTimeCalib;   ///< Container for the High Gain time calibration coefficients

  ClassDefNV(MpdEmcCalibParams, 1);
};
#endif
