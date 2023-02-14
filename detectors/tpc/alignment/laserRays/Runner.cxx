#include "Runner.h"
#include <exception>
#include <iostream>
#include "EigenCalculator.h"
#include "Track.h"
#include "Point.h"
#include "FileHelper.h"
#include "Enums.h"
#include "CorrectionMatrix.h"

using std::exception;
using std::runtime_error;
using std::string;

namespace TpcAlignmentLaserRays {
Runner::Runner() : vDebug(false), vRays(vDebug), vCorrectionMatrix(vDebug) {}

void Runner::SetDebugMode(bool DebugMode)
{
   vDebug.SetDebugMode(DebugMode);
}

void Runner::Calibrate(unsigned NumberOfCalibration)
{
   for (int time = 0; time < NumberOfCalibration; ++time) {
      try {
         vSLE.BuildAndSolve(vRays, vCorrectionMatrix);
      } catch (const std::exception &ex) {
         throw runtime_error(
            string("Calibrate on loop ").append(std::to_string(time)).append(" throws an exception: ").append(ex.what()));
      } catch (...) {
         throw runtime_error(string("Calibrate on loop ")
                                .append(std::to_string(time))
                                .append(" throws an exception: ")
                                .append("Unknown exception"));
      }
   }
}

void Runner::SaveAMR2Files(string outputCorectionMatrix, string outputMR, int precision)
{
   if (!outputCorectionMatrix.empty()) {
      string filepath2A = FileHelper::BuildFilePath(Solution::release, Direction::output, outputCorectionMatrix, "");
      try {
         vCorrectionMatrix.SaveMatrix2File(filepath2A, precision);
      } catch (const exception &ex) {
         throw runtime_error(string("SaveAMR2Files: ").append(ex.what()));
      }
   }

   if (!outputMR.empty()) {
      string filepath2MR = FileHelper::BuildFilePath(Solution::release, Direction::output, outputMR, "");
      try {
         FileHelper::SaveMR(filepath2MR, vSLE.GetM(), vSLE.GetR(), precision);
      } catch (const exception &ex) {
         std::cout << ex.what();
         vDebug.Print("%c", ex.what());
      }
   }
}

void Runner::LoadModelData(std::string RaysData)
{
   try {
      vRays.Read4File(RaysData.c_str());
   } catch (const exception &ex) {
      throw runtime_error(string("Rays filepath: ").append(RaysData).append(" : Error: ").append(ex.what()));
   } catch (...) {
      throw runtime_error(string("Rays filepath: ").append(RaysData).append(" : Unknowns error."));
   }
   if (vRays.Get().empty()) {
      throw runtime_error("Running FillInRays procedure get emtpy rays set.");
   }
}

void Runner::CorrectTracks(string inputFileName, string outputFileName)
{
   Rays vRays(vDebug);
   if (inputFileName.empty()) {
      throw runtime_error("inputFileName is not set");
   }
   vRays.Read4File(inputFileName.c_str());
   Rays vResult(vDebug);

   for (Track const &track : vRays.Get()) {
      vResult.AddTrack(CorrectTrack(track));
   }
   vResult.Write2File(outputFileName.c_str());
}

void Runner::LoadCorrectionMatrix(string CorrectionMatrixFile, bool FromFile)
{
   bool isSet{false};
   if (FromFile) {
      vCorrectionMatrix.LoadMatrix4File(CorrectionMatrixFile);
   } else {
      vCorrectionMatrix.SetMatrixWithDefaultValues();
   }
}

Track Runner::CorrectTrack(Track const &track)
{
   Track result(vDebug, track.Id());
   for (size_t i = 0; i < track.Get().size(); ++i) {
      result.AddPoint(vCorrectionMatrix.CorrectPoint(track.GetPoint(i)));
   }
   return result;
}
} // namespace TpcAlignmentLaserRays
