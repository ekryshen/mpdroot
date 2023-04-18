#ifndef EIGENCOEFFICIENTS_HH
#define EIGENCOEFFICIENTS_HH
#include "Point.h"

namespace TpcAlignment
{
/// <summary>
/// Eigenvalues
/// </summary>
struct EigenCoefficients
{
	/// <summary>
	/// major axis of the ellipsoid of inertia
	/// </summary>
	Point P;
	/// <summary>
	/// Center of mass of track points
	/// </summary>
	Point Q;
};
}

#endif
