#ifndef ALIGNHELPER_HH
#define ALIGNHELPER_HH
#include <cstdio>
#include <vector>
#include <array>
#include <string>
#include "Constants.h"
#include "Types.h"

namespace TpcAlignmentTest
{
	class CorrectionMatrixTest;
}

namespace TpcAlignment {
	class Debug;
	class Point;

	/// <summary>
	/// Correction procedures
	/// </summary>
	class CorrectionMatrix {
		friend class TpcAlignmentTest::CorrectionMatrixTest;
	private:
		SectorMatrix vCorrectionMatrix;
		Debug const& vDebug;

	public:

		/// <summary>
		/// Initializing object
		/// </summary>
		/// <param name="a">A matrix </param>
		/// <param name="inMode">Global debug object</param>
		CorrectionMatrix(Debug const& inMode);

		~CorrectionMatrix();

		SectorMatrix const& Get() const;

		SectorMatrix& Get();

		/// <summary>
		/// Calculate new coordinates
		/// </summary>
		/// <param name=inputPoint>input Point object</param>
		/// <returns>corrected Point</returns>
		Point CorrectPoint(Point const& inputPoint) const;

		/// <summary>
		/// Set Correction matrix to default values
		/// </summary>
		void SetMatrixWithDefaultValues();

		/// <summary>
		/// 
		/// </summary>
		/// <param name="filename"></param>
		void LoadMatrix4File(std::string filename);

		/// <summary>
		/// 
		/// </summary>
		/// <param name="filename"></param>
		/// <param name="precision"></param>
		void SaveMatrix2File(std::string filename, int precision);
	protected:
		/// <summary>
		/// Calculate new value for one of the coordinate
		/// </summary>
		/// <param name=point>Point object</param>
		/// <param name=coordinate>Coordinate's index(X=0,Y=1,Z=2)</param>
		/// <returns>New coordinate</returns>
		double  CalculateCoordinatesNewValue(Point const& point,
			unsigned coordinate) const;
	};

} // namespace TpcAlignHelper

#endif

