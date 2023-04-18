#ifndef RAYS_HH
#define RAYS_HH
#include <vector>
#include "Debug.h"
#include <string>
#include "Types.h"
#include "Track.h"

namespace TpcAlignmentTest
{
	class RaysTest;
}
namespace TpcAlignment {
	class Track;
	/// <summary>
	/// Class Rays contains all Laser Tracks.
	/// </summary>
	class Rays
	{
		friend class TpcAlignmentTest::RaysTest;
	private:
		RaysType vRays;
		Debug const& vDebug;
	public:
		/// <summary>
		/// Init default values
		/// </summary>
		/// <param name="inDebug">Debug object</param>
		explicit Rays(Debug const& inDebug);
		~Rays() = default;

		/// <summary>
		/// Add Track to collection
		/// </summary>
		/// <param name="t">track</param>
		void AddTrack(Track const & t);

		/// <summary>
		/// Move Track to collection
		/// </summary>
		/// <param name="t">track</param>
		void AddTrack(Track && t);
		/// <summary>
		/// Return tracks collection
		/// </summary>
		/// <returns>Collection of track</returns>
		RaysType Get();

		/// <summary>
		/// Return tracks collection
		/// </summary>
		/// <returns>Collection of track</returns>
		RaysType const & Get() const;

		/// <summary>
		/// Return track object by collection's index.
		/// </summary>
		/// <param name="i">index</param>
		/// <returns>track object</returns>
		Track const & GetTrack(size_t i) const;

		/// <summary>
		/// Read laser rays data from a file
		/// </summary>
		/// <param name="RaysFilename">Filename of laser rays</param>
		void Read4File(const char* RaysFilename);

		/// <summary>
		/// Write collection of tracks into file
		/// </summary>
		/// <param name="RaysFilename">A Filename</param>
		void Write2File(const char* RaysFilename);
	protected:
		/// <summary>
		/// Write collection to filestream
		/// </summary>
		/// <param name="fp">Filestream</param>
		/// <param name="filename">Filename for error message</param>
		void SaveRays(FILE* fp, const char* filename);

		/// <summary>
		/// Read track collection from filestream
		/// </summary>
		/// <param name="RaysFilename">file name</param>
		/// <param name="fp">file descriptor</param>
		void ReadRays(const char* RaysFilename, FILE* fp);

		/// <summary>
		/// Open file 
		/// </summary>
		/// <param name="RaysFilename"></param>
		/// <param name="mode"></param>
		/// <returns>file descriptor</returns>
		static FILE* GetFileStream(const char* RaysFilename, std::string mode);
	};
}
#endif
