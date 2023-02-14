#ifndef ENUMS_HH
#define ENUMS_HH

namespace TpcAlignmentLaserRays {
/// <summary>
/// Input/output directory
/// </summary>
enum class Direction : int { input = 0, output = 1 };
/// <summary>
/// Test or Release
/// </summary>
enum class Solution : int { test = 0, release = 1 };
} // namespace TpcAlignmentLaserRays

#endif
