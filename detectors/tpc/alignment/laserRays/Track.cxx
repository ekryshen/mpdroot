#include "Track.h"
#include "FileHelper.h"
#include <stdexcept>
#include <string>

namespace TpcAlignmentLaserRays {
Track::Track(Debug const &inDebug, size_t TrackNumber) : vDebug(inDebug), vTrackId(TrackNumber) {}

void Track::AddPoint(Point const &point)
{
   vTrack.push_back(point);
}

std::vector<Point> const &Track::Get() const
{
   return vTrack;
}

std::vector<Point> &Track::Get()
{
   return vTrack;
}

Point const &Track::GetPoint(size_t i) const
{
   return vTrack.at(i);
}
Point &Track::GetPoint(size_t i)
{
   return vTrack[i];
}

void Track::WriteTrack(FILE *fp, std::string errorMsg) const
{
   size_t    PointsOnTrack{vTrack.size()};
   int const writeNumber = fprintf(fp, " %u", static_cast<unsigned>(PointsOnTrack));
   if (writeNumber == 0) {
      throw std::runtime_error(errorMsg.append(" Can't write Point quantaty on track."));
   }

   for (Point point : vTrack) {
      point.Write2File(fp, errorMsg);
   }
   fprintf(fp, " %c", '\n');
}

void Track::ReadTrack(FILE *fp, std::string errorMsg)
{
   int NumberOfPointsOnTrack = ReadNumberOfPointsOnTrack(fp, errorMsg);
   if (NumberOfPointsOnTrack == 0) {
      throw std::runtime_error(errorMsg.append("Empty track"));
   }
   vDebug.Print("track open\n");
   for (int i = 0; i < NumberOfPointsOnTrack; ++i) {
      Point point;
      point.Read4File(fp, "ReadTrack");
      point.SetSector();
      AddPoint(point);
   }
   char      c;
   int const readChar = fscanf(fp, "%c", &c);
   if (readChar == 0) {
      throw std::runtime_error(std::string("Can't get EOL or EOF"));
   }
}

size_t Track::Id() const
{
   return vTrackId;
}

unsigned Track::ReadNumberOfPointsOnTrack(FILE *fStream, std::string errorMsg) const
{
   unsigned PointsOnTrack{0};
   int      readNumber{0};
   if (!feof(fStream)) {
      readNumber = fscanf(fStream, "%u", &PointsOnTrack);
      if (readNumber == 0) {
         throw std::invalid_argument(errorMsg.append("can't read Points On Track's value."));
      }
   } else {
      throw std::invalid_argument(errorMsg.append("File is empty"));
   }
   return PointsOnTrack;
}
} // namespace TpcAlignmentLaserRays
