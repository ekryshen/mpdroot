#ifndef RUNNER_HH
#define RUNNER_HH
#include <vector>
#include <array>
#include <string>
#include "CorrectionMatrix.h"
#include "Rays.h"
#include "Debug.h"
#include "SystemOfLinearEquations.h"
#include "Constants.h"
#include "EigenCoefficients.h"

namespace TpcAlignmentTest {
class RunnerTest;
class RunnerTestHelper;
} // namespace TpcAlignmentTest

namespace TpcAlignmentLaserRays {
class Track;
class Point;
/// <summary>
/// Library entry point
/// </summary>
class Runner {
   friend class TpcAlignmentTest::RunnerTest;
   friend class TpcAlignmentTest::RunnerTestHelper;

private:
   Debug                   vDebug;
   Rays                    vRays;
   CorrectionMatrix        vCorrectionMatrix;
   SystemOfLinearEquations vSLE;

public:
   /// <summary>
   /// Init vDebug, vRays, vAlignHelper and IsLoadA4File variables
   /// </summary>
   /// <param name="isLoadA4File">If A-file doesn't exist then false</param>
   /// <param name="DebugMode">By default DebugMode is switched off</param>
   Runner();
   ~Runner() = default;

   void SetDebugMode(bool DebugMode);

   /// <summary>
   /// Calibrate Tpc-chamber
   /// </summary>
   /// <param name="NumberOfCalibration">Number of loop to calibrate</param>
   void Calibrate(unsigned NumberOfCalibration);

   /// <summary>
   /// Set input data source
   /// </summary>
   /// <param name="RaysData">File with model data</param>
   void LoadModelData(std::string RaysData);

   /// <summary>
   /// Correct coordinates using alignment-coefficients
   /// </summary>
   /// <param name="inputFileName">a Name of input file with Tracks</param>
   /// <param name="outputFileName">a name of output file with corrected data</param>
   void CorrectTracks(std::string inputFileName, std::string outputFileName);

   /// <summary>
   /// Read Correction coefficients(A-matrix) from a file if needed
   /// </summary>
   /// <param name="CorrectionMatrixFile">file name path</param>
   /// <param name="FromFile">From a file</param>
   void LoadCorrectionMatrix(std::string CorrectionMatrixFile, bool FromFile = true);

   /// <summary>
   /// Save alignment's coefficients to the files
   /// </summary>
   /// <param name="outputCorectionMatrix">Correction coefficients matrix</param>
   /// <param name="outputMR">M, R-coefficients</param>
   /// <param name="precision">precision to save data</param>
   void SaveAMR2Files(std::string outputCorectionMatrix, std::string outputMR, int precision);

protected:
   /// <summary>
   /// Correct Track
   /// </summary>
   /// <param name="track">Track object</param>
   /// <returns>Corrected track</returns>
   Track CorrectTrack(Track const &track);
};
} // namespace TpcAlignmentLaserRays
#endif
