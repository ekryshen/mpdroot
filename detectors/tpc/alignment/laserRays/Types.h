#ifndef TYPES_HH
#define TYPES_HH
#include <array>
#include <vector>
#include "Constants.h"

namespace TpcAlignmentLaserRays {
class Track;
typedef std::array<std::vector<double>, SECTORS> SectorMatrix;
typedef std::vector<Track>                       RaysType;
} // namespace TpcAlignmentLaserRays
#endif
