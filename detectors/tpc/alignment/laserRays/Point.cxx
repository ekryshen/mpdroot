#include "Point.h"
#include <cstdio>
#include <string>
#include <stdexcept>
#include "CorrectionMatrix.h"

namespace TpcAlignmentLaserRays {
Point::Point(double x, double y, double z) : X(x), Y(y), Z(z) {}

void Point::SetSector()
{
   vSectorId = GetSectorNumberByCoordinate(X, Y, Z);
}

int Point::Sector() const
{
   if (vSectorId < 0) {
      throw std::runtime_error("Sector is not set.");
   }
   return vSectorId;
}

void Point::Read4File(FILE *fp, std::string errorMsg)
// const char* TrackFilename,
// int readNumber,
// int vTrackNumber)
{
   X = ReadCoordinate(fp, errorMsg.append(" Coordinate: X: "));
   Y = ReadCoordinate(fp, errorMsg.append(" Coordinate: Y: "));
   Z = ReadCoordinate(fp, errorMsg.append(" Coordinate: Z: "));
}

void Point::Write2File(FILE *fp, std::string errorMsg)
{
   WriteCoordinate(X, fp, errorMsg.append(": X"));
   WriteCoordinate(Y, fp, errorMsg.append(": Y"));
   WriteCoordinate(Z, fp, errorMsg.append(": Z"));
}

Point &Point::operator+=(const Point &rhs)
{
   this->X += rhs.X;
   this->Y += rhs.Y;
   this->Z += rhs.Z;
   return *this;
}

Point Point::operator-(const Point &rhs)
{
   Point result(*this);
   result.X -= rhs.X;
   result.Y -= rhs.Y;
   result.Z -= rhs.Z;
   return result;
}

Point &Point::operator/=(double val)
{
   this->X /= val;
   this->Y /= val;
   this->Z /= val;
   return *this;
}
Point Point::operator/(double val)
{
   Point result(*this);
   result.X /= val;
   result.Y /= val;
   result.Z /= val;
   return result;
}

void Point::WriteCoordinate(double coordinate, FILE *fStream, std::string errorMsg)
{
   int writeNumber{0};
   writeNumber = fprintf(fStream, " %f", coordinate);
   if (writeNumber == 0) {
      throw std::runtime_error(errorMsg);
   }
}

double Point::ReadCoordinate(FILE *fp, std::string errorMsg)
{
   double r;
   int    readNumber = fscanf(fp, "%lf", &r);
   if (readNumber == 0) {
      throw std::invalid_argument(errorMsg.c_str());
   }
   return r;
}

unsigned Point::GetSectorNumberByPoint(Point const &point)
{
   return GetSectorNumberByCoordinate(point.X, point.Y, point.Z);
}

unsigned Point::GetSectorNumberByCoordinate(double x, double y, double z)
{
   double fi = atan2(x, y);
   if (fi > 1.E+9) {
      throw std::runtime_error(std::string("fi > 1.E+9: ").append(std::to_string(fi)));
   }
   int segment_number = static_cast<int>(floor(fi * rPI + 12.5));
   segment_number     = segment_number % 12;
   if (z < 0.) {
      segment_number = segment_number + 12;
   }
   if (segment_number < 0) {
      throw std::runtime_error(std::string("Sector couldn't be a negative number: ")
                             .append(std::to_string(segment_number)
                                        .append(" Check point: (")
                                        .append(std::to_string(x))
                                        .append(", ")
                                        .append(std::to_string(y))
                                        .append(", ")
                                        .append(std::to_string(z))
                                        .append(")")));
   }
   return static_cast<unsigned>(segment_number);
}
} // namespace TpcAlignmentLaserRays
