#include <iostream>
#include <string>

#include "FileHelper.h"

namespace TpcAlignment {

bool FileHelper::FileExists(FILE *file, const char *filename, const char *msg)
{
   if (file == NULL) {
      std::string ErrMsg{filename};
      ErrMsg.append(" ").append(msg).append(". Failed: ");
      perror(ErrMsg.c_str());
      return false;
   }
   return true;
}

int FileHelper::ReadPointsOnTrack(FILE *fStream, const char *filename, const char *method)
{
   int PointsOnTrack{0};
   int readNumber{0};
   if (!feof(fStream)) {
      readNumber = fscanf(fStream, "%d", &PointsOnTrack);
      if (readNumber == 0) {
         std::cout << method << ": " << filename << ": can't read Points On Track's value." << '\n';
         return 0;
      }
   } else {
      std::cout << method << ": " << filename << ": is empty" << '\n';
      return 0;
   }
   return PointsOnTrack;
}

bool FileHelper::WritePointsOnTrack(FILE *outputStream, int PointsOnTrack, const char *filename, const char *method)
{
   int writtenNumber{0};
   writtenNumber = fprintf(outputStream, " %d", PointsOnTrack);
   if (writtenNumber == 0) {
      std::cout << method << ": "
                << ": can't write Points On Track's Value to the File: " << filename << '\n';
      return false;
   }
   return true;
}

bool FileHelper::ReadCoordinate(double &coordinate, FILE *fStream, const char *filename, const char *method,
                                const char *coordinateName, int Track, int Point)
{
   int readNumber{0};
   // FIXME: Does lf works on Unix?
   readNumber = fscanf(fStream, "%lf", &coordinate);
   if (readNumber == 0) {
      std::cout << method << ": " << filename << ": ReadCoordinate: Error reading coordinate " << coordinateName
                << " Track: " << Track << ", Point: " << Point << '\n';
      return false;
   }
   return true;
}

bool FileHelper::WriteCoordinate(double coordinate, FILE *fStream, const char *filename, const char *method,
                                 const char *coordinateName, int Track, int Point)
{
   int writeNumber{0};
   writeNumber = fprintf(fStream, " %f", coordinate);
   if (writeNumber == 0) {
      std::cout << method << ": " << filename << ": WriteCoordinate: Error writing coordinate " << coordinateName
                << coordinate << " Track: " << Track << ", Point: " << Point << '\n';
      return false;
   }
   return true;
}

bool FileHelper::ReadXYZ(double &x, double &y, double &z, FILE *fStream, const char *filename, const char *method,
                         int Track, int Point)
{
   if (!ReadCoordinate(x, fStream, filename, method, "X", Track, Point)) {
      return false;
   }
   if (!ReadCoordinate(y, fStream, filename, method, "Y", Track, Point)) {
      return false;
   }
   if (!ReadCoordinate(z, fStream, filename, method, "Z", Track, Point)) {
      return false;
   }
   return true;
}

bool FileHelper::WriteXYZ(double x, double y, double z, FILE *fStream, const char *filename, const char *method,
                          int Track, int Point)
{
   if (!WriteCoordinate(x, fStream, filename, method, "X", Track, Point)) {
      return false;
   }
   if (!WriteCoordinate(y, fStream, filename, method, "Y", Track, Point)) {
      return false;
   }
   if (!WriteCoordinate(z, fStream, filename, method, "Z", Track, Point)) {
      return false;
   }
   return true;
}
} // namespace TpcAlignment
