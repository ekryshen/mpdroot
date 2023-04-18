#ifndef TRACK_HH
#define TRACK_HH
#include <vector>
#include <string>
#include "Debug.h"
#include "Point.h"

namespace TpcAlignmentTest
{
	class TrackTest;
}

namespace TpcAlignment {
	/// <summary>
	/// a mark or line of marks left by a particle
	/// </summary>
	class Track
	{
		friend class TpcAlignmentTest::TrackTest;
	private:
		std::vector<Point> vTrack;
		Debug const& vDebug;
		size_t vTrackId;
	public:
		/// <summary>
		/// Init Debug mode and TrackId
		/// </summary>
		/// <param name="inDebug">Debug object</param>
		/// <param name="TrackNumber">track index in Rays collection</param>
		Track(Debug const& inDebug, size_t TrackNumber);
		~Track() = default;

		/// <summary>
		/// Add point to the track
		/// </summary>
		/// <param name="point">point object</param>
		void AddPoint(Point const & point);

		/// <summary>
		/// Track as a vector of points
		/// </summary>
		/// <returns></returns>
		std::vector<Point> const & Get() const;

		/// <summary>
		/// Track as a vector of points
		/// </summary>
		/// <returns></returns>
		std::vector<Point> & Get();

		/// <summary>
		/// Get point on the track by its index
		/// </summary>
		/// <param name="index">index</param>
		/// <returns>point</returns>
		Point const & GetPoint(size_t index) const;

		/// <summary>
		/// Get point on the track by its index
		/// </summary>
		/// <param name="index">index</param>
		/// <returns></returns>
		Point & GetPoint(size_t index);

		/// <summary>
		/// Write Track to the filestream
		/// </summary>
		/// <param name="fp">filestream</param>
		/// <param name="errorMsg">error message</param>
		void WriteTrack(FILE *fp, std::string errorMsg) const;

		/// <summary>
		/// Read track from a filestream
		/// </summary>
		/// <param name="fp">filestream</param>
		/// <param name="errorMsg">error message</param>
		void ReadTrack(FILE *fp, std::string errorMsg);

		/// <summary>
		/// Get Track Identification number
		/// This is an index in the collection of tracks
		/// </summary>
		/// <returns>ID</returns>
		size_t Id() const;

	private:
		unsigned ReadNumberOfPointsOnTrack(FILE* fStream, std::string errorMsg) const;
	};
}
#endif