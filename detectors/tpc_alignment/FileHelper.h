#ifndef FILEHELPER_HH
#define FILEHELPER_HH
#include <string>
#include "Types.h"

namespace TpcAlignmentTest {
	class FileHelperTest;
}

namespace TpcAlignment {
	class Point;
	enum class Direction;
	enum class Solution;

	/// <summary>
	/// Help to deal with files.
	/// </summary>
	class FileHelper {
		friend class TpcAlignmentTest::FileHelperTest;
	private:
		const char* vTrackFilename;
		const char* vCoefficientFilename;

	public:
		FileHelper() = delete;
		/// <summary>
		/// Throw runtime_error if file doesn't exist.
		/// </summary>
		/// <param name="file">file path</param>
		/// <param name="ErrMsg">error message</param>
		static void FileExists(FILE* file, std::string ErrMsg);


		/// <summary>
		///  Build file path according rules
		/// </summary>
		/// <param name="solution">release/test folder</param>
		/// <param name="direction">input/output folder</param>
		/// <param name="filename">file name</param>
		/// <returns>filepath</returns>
		static std::string BuildFilePath(Solution solution, 
											Direction direction,
											std::string filename);

		/// <summary>
		/// Build file path according rules
		/// </summary>
		/// <param name="solution">release/test folder</param>
		/// <param name="direction">input/output</param>
		/// <param name="filename">file name</param>
		/// <param name="innerPath">prefix to file path</param>
		/// <returns>filepath</returns>
		static std::string BuildFilePath(Solution solution, Direction direction,
			std::string filename, std::string innerPath);



		/// <summary>
		/// Build file path according rules
		/// </summary>
		/// <param name="direction">input/output</param>
		/// <param name="filename">file name</param>
		/// <param name="innerPath">prefix to file path</param>
		/// <param name="track">track number</param>
		/// <returns>filepath</returns>
		static std::string BuildFilePath(Solution solution, Direction direction,
			std::string filename, std::string innerPath, unsigned track);

		/// <summary>
		/// Build file path according some rules
		/// </summary>
		/// <param name="direction"></param>
		/// <param name="track"></param>
		/// <param name="filename">file name</param>
		/// <param name="innerPath"></param>
		/// <returns>filepath</returns>
		static std::string BuildFilePath(Solution solution, Direction direction,
			unsigned track, std::string filename, std::string innerPath);

		/// <summary>
		/// Build file path according some rules
		/// </summary>
		/// <param name="direction"></param>
		/// <param name="track"></param>
		/// <param name="point"></param>
		/// <param name="filename">file name</param>
		/// <param name="innerPath"></param>
		/// <returns>filepath</returns>
		static std::string BuildFilePath(Solution solution, Direction direction,
			unsigned track, size_t point, std::string filename, std::string innerPath);

		/// <summary>
		/// File descriptor by Filename
		/// </summary>
		/// <param name="RaysFilename">file name</param>
		/// <param name="mode">add, write - file open modes</param>
		/// <returns>file descriptor</returns>
		static FILE* get_file_descriptor(const char* RaysFilename, std::string mode);

		/// <summary>
		/// Save M and R - matrices in a file
		/// </summary>
		/// <param name="filename">filepath to MR data</param>
		/// <param name="vMexpected">M-coeff</param>
		/// <param name="vRexpected">R-coeff</param>
		/// <param name="precision">double precision</param>
		static void SaveMR(std::string filename, SectorMatrix const& vM, SectorMatrix const& vR, int precision);

		/// <summary>
		/// Read MR 4 file to arguments
		/// </summary>
		/// <param name="filename">filepath to MR data</param>
		/// <param name="vMexpected">M-coeff</param>
		/// <param name="vRexpected">R-coeff</param>
		static void LoadMR(std::string filename, SectorMatrix& vM, SectorMatrix& vR);

		/// <summary>
		/// Read correction coefficient matrix from a file
		/// </summary>
		/// <param name="filename">Filepath to correction coefficient matrix data file</param>
		/// <param name="vCorrectionMatrix">correction coefficient matrix</param>
		static void LoadCorrectionMatrix(std::string filename, SectorMatrix& vCorrectionMatrix);

		/// <summary>
		/// Using fstream
		/// </summary>
		/// <param name="filename"></param>
		/// <param name="vCorrectionMatrix"></param>
		/// <param name="precision">digits saving in doubles</param>
		static void SaveCorrectionMatrix(std::string filename, SectorMatrix const& vCorrectionMatrix, int precision);
		/**********************************************************************/
		/// <summary>
		/// Load correction coefficient matrix from a file
		/// </summary>
		/// <param name="filename">Filepath to correction coefficient matrix data file</param>
		/// <param name="vCorrectionMatrix">correction coefficient matrix</param>
		static void ReadCorrectionMatrix_obsolete(std::string filename, SectorMatrix& vCorrectionMatrix);

		/// <summary>
		/// Save correction coefficient matrix to a file
		/// </summary>
		/// <param name="filename">Filepath to correction coefficient matrix data file</param>
		/// <param name="vCorrectionMatrix">correction coefficient matrix</param>
		static void SaveCorrectionMatrix_obsolete(std::string filename, SectorMatrix const& vCorrectionMatrix);


		/// <summary>
		/// Read MR 4 file to arguments
		/// </summary>
		/// <param name="filename">filepath to MR data</param>
		/// <param name="vMexpected">M-coeff</param>
		/// <param name="vRexpected">R-coeff</param>
		static void ReadMR_obsolete(std::string filename, SectorMatrix& vM, SectorMatrix& vR);
	private:
		/// <summary>
		/// Generate prefix, path, suffix by Solution and Direction
		/// </summary>
		/// <param name="solution">Release/test</param>
		/// <param name="direction">Input/output</param>
		/// <returns>prefix, path, suffix</returns>
		static std::tuple<std::string, std::string, std::string> Path(Solution solution, Direction direction);



	};
} // namespace TpcAlignment

#endif