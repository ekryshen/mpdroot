#ifndef FILEHELPER_HH
#define FILEHELPER_HH
namespace TpcAlignment {

class FileHelper {
private:
   const char *vTrackFilename;
   const char *vCoefficientFilename;

public:
   FileHelper() = delete;
   static bool FileExists(FILE *file, const char *filename, const char *msg);
   static int  ReadPointsOnTrack(FILE *fStream, const char *filename, const char *method);
   static bool WritePointsOnTrack(FILE *outputStream, int PointsOnTrack, const char *filename, const char *method);
   static bool ReadXYZ(double &x, double &y, double &z, FILE *fStream, const char *filename, const char *method,
                       int Track, int Point);
   static bool WriteXYZ(double x, double y, double z, FILE *fStream, const char *filename, const char *method,
                        int Track, int Point);
   static bool ReadCoordinate(double &coordinate, FILE *fStream, const char *filename, const char *method,
                              const char *coordinateName, int Track, int Point);
   static bool WriteCoordinate(double coordinate, FILE *fStream, const char *filename, const char *method,
                               const char *coordinateName, int Track, int Point);
};
} // namespace TpcAlignment
#endif
