#ifndef POINT_HH
#define POINT_HH
#include <string>

namespace TpcAlignmentTest
{
	class PointTest;
}

namespace TpcAlignment {
	/// <summary>
	/// Helper point class
	/// </summary>
	class Point
	{
		friend class TpcAlignmentTest::PointTest;
	private:
		int vSectorId{ -1 };
	public:
		/// <summary>
		/// Coordinate X
		/// </summary>
		double X{ 0. };
		/// <summary>
		/// Coordinate Y
		/// </summary>
		double Y{ 0. };
		/// <summary>
		/// Coordinate Z
		/// </summary>
		double Z{ 0. };

	public:
		Point() = default;
		~Point() = default;
		/// <summary>
		/// Ctor with coordinates
		/// </summary>
		/// <param name="x">X</param>
		/// <param name="y">Y</param>
		/// <param name="z">Z</param>
		Point(double x, double y, double z);

		/// <summary>
		/// Set Sector Id: [0..23]
		/// </summary>
		void SetSector();

		int Sector() const;

		/// <summary>
		/// Read point object from a file
		/// </summary>
		/// <param name="fp">file descriptor</param>
		/// <param name="errorMsg">error message</param>
		void Read4File(FILE* fp, std::string errorMsg);

		/// <summary>
		/// Write point object to a file
		/// </summary>
		/// <param name="fp">file descriptor</param>
		/// <param name="errorMsg">error message</param>
		void Write2File(FILE* fp, std::string errorMsg);

		/// <summary>
		/// unary sum
		/// </summary>
		/// <param name="rhs">right value</param>
		/// <returns>Resulted point</returns>
		Point& operator+=(const Point& rhs);

		/// <summary>
		/// subtruction operator
		/// </summary>
		/// <param name="rhs">right value</param>
		/// <returns>Resulted point</returns>
		Point operator-(const Point& rhs);

		/// <summary>
		/// unary devision
		/// </summary>
		/// <param name="val"> value</param>
		/// <returns>Resulted point</returns>
		Point& operator/=(double val);

		/// <summary>
		///  devision
		/// </summary>
		/// <param name="val">value</param>
		/// <returns>Resulted point</returns>
		Point operator/(double val);

		/// <summary>
		/// Calculate Sector's Number by Point's coordinates.
		/// </summary>
		/// <param name="point">Point object</param>
		/// <returns>Segment Number</returns>
		static unsigned GetSectorNumberByPoint(Point const& point);
	protected:

		/// <summary>
		/// Calculate Sector by point XYZ-coordinates.
		/// </summary>
		/// <param name="x">X</param>
		/// <param name="y">Y</param>
		/// <param name="z">Z</param>
		/// <returns>Sector's Index</returns>
		static unsigned GetSectorNumberByCoordinate(double x, double y, double z);

		/// <summary>
		/// Write coordinate to filestream
		/// </summary>
		/// <param name="coordinate">x/y/z</param>
		/// <param name="fStream">file stream</param>
		/// <param name="errorMsg">error message</param>
		static void WriteCoordinate(double coordinate, FILE* fStream, std::string errorMsg);

		/// <summary>
		/// Read coordinate from filestream
		/// </summary>
		/// <param name="fp">filestream</param>
		/// <param name="errorMsg">error message</param>
		/// <returns>coordinate</returns>
		static double ReadCoordinate(FILE* fp, std::string errorMsg);
	};
}
#endif


